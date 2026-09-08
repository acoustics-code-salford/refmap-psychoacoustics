require(insight)
require(car)

# =============================================================================
# brms_vif(): collinearity diagnostics for ANY brms model
#
# Rationale: collinearity is a property of the design matrix X (the formula's
# right-hand side + data) - it does not depend on the response family, link,
# or how the model is actually fit. So a single VIF routine works for any
# brmsfit regardless of family (ordbeta, cumulative, gaussian, ...), the same
# way glmgee_vif() worked for any glmgee family - we just need the RHS
# formula and data, refit as a plain lm() (response is irrelevant to VIF, so
# a dummy response is used), and run car::vif() on that.
#
# WHY NOT performance::check_collinearity()? Confirmed TWICE in this session
# (once for glmtoolbox::glmgee, once for a brms::cumulative() ordinal model)
# that its reported "adj. VIF" for multi-df terms (factors with >2 levels,
# any interaction involving them) silently assumes Df = 1 for every term,
# regardless of the term's true Df - i.e. it reports sqrt(raw VIF) across the
# board rather than raw_VIF^(1/(2*Df)). This produces badly inflated adj. VIF
# for genuine multi-df terms (e.g. a 3-level factor's true adj. VIF getting
# reported ~5x too high). Until this is confirmed fixed upstream (this looks
# like a bug in insight's/performance's generic term-assignment path, not
# specific to one model class), this function bypasses that path entirely by
# using car::vif() directly, which gets Df right.
#
# WHAT THIS CAN'T HANDLE: brms-specific formula syntax that has no meaning in
# a plain stats::lm() call - smooth terms (s(), t2()), Gaussian processes
# (gp()), monotonic effects (mo()), measurement-error terms (me()),
# category-specific ordinal terms (cs()), or non-linear (nlf()) formulas.
# Collinearity/VIF is inherently a linear/parametric-term concept and doesn't
# generalise cleanly to smooth bases anyway - for GAM-type smooth terms, the
# analogous diagnostic is CONCURVITY, not VIF; see mgcv::concurvity() on an
# equivalent mgcv::gam() fit instead. This function will error informatively
# (via the model.matrix()/lm() failure) rather than silently mishandle these.
#
# NOT executed here - no R access in this environment. Please run and report
# back if anything errors or looks off, same as the other scripts this
# session.
# =============================================================================

# -----------------------------------------------------------------------------
# Design-matrix condition number (Belsley, Kuh & Welsch 1980)
# -----------------------------------------------------------------------------
#
# WHY THIS IS SEPARATE FROM TERM-LEVEL VIF: VIF assesses collinearity
# term-by-term (how well can THIS term be predicted from the others); the
# condition number reflects the OVERALL design matrix's numerical
# conditioning, which can be poor because of a diffuse combination across
# MANY columns jointly even when no single term's VIF is high (as confirmed
# directly in this session: mA1's terms were all Adjusted_VIF < 4, yet
# kappa(X) ~ 500). The two diagnostics are complementary, not redundant -
# neither one implies the other, so both are reported.
#
# WHY NOT A NAIVE kappa(X) CALL: computing the condition number on the raw,
# unscaled model matrix conflates two different phenomena Belsley (1991)
# calls "essential" and "non-essential" ill-conditioning - the latter being
# an artefact of predictors simply living on very different natural scales
# (e.g. Age in years next to a 0/1 dummy), which inflates the condition
# number without reflecting genuine collinearity BETWEEN predictors at all.
# The standard fix is to scale each column to unit Euclidean norm (NOT
# centre - centring the intercept column would remove the ability to detect
# non-essential ill-conditioning in the first place) before taking the SVD.
# Both the scaled (recommended) and raw kappa() figures are returned, so a
# large gap between them is itself diagnostic of how much of the raw number
# was a scaling artefact versus genuine multicollinearity.
#
# THRESHOLDS: Belsley, Kuh & Welsch (1980) suggest condition indices above
# ~30 indicate moderate collinearity and above ~100 indicate severe
# collinearity; the same thresholds are conventionally applied to the
# largest (overall) condition number of the scaled matrix.
.design_condition_diagnostics <- function(X) {
  norms <- sqrt(colSums(X^2))
  norms[norms == 0] <- 1  # guard against a degenerate constant/zero column
  X_scaled <- sweep(X, 2, norms, "/")
  
  sv <- svd(X_scaled)$d
  sv <- sv[sv > .Machine$double.eps]
  cond_indices <- max(sv) / sv
  
  list(
    condition_number_scaled      = max(cond_indices),
    condition_number_raw         = tryCatch(kappa(X, exact = TRUE), error = function(e) NA_real_),
    n_condition_indices_over_30  = sum(cond_indices > 30),
    n_condition_indices_over_100 = sum(cond_indices > 100)
  )
}


brms_vif <- function(model,
                     dpar = NULL,
                     digits = 2,
                     sort = TRUE,
                     thresholds = c(moderate = 5, high = 10),
                     seed = 1) {
  
  if (!inherits(model, "brmsfit")) {
    stop("`model` must be a fitted brms::brmsfit object.")
  }
  for (pkg in c("insight", "car")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required but not installed.", pkg))
    }
  }
  
  # --- extract the target formula (main/conditional, or a named auxiliary
  #     parameter formula such as phi/sigma via dpar) and the model's data --
  formula_list <- insight::find_formula(model)
  target_name <- if (is.null(dpar)) "conditional" else dpar
  
  if (is.null(formula_list[[target_name]])) {
    stop(sprintf(
      "Could not find a '%s' formula component on this model. ",
      target_name),
      "Available components: ", paste(names(formula_list), collapse = ", "),
      ". If checking an auxiliary parameter (e.g. phi), pass dpar to match ",
      "the name shown above (commonly matches the dpar name used in bf())."
    )
  }
  target_formula <- formula_list[[target_name]]
  
  dat <- insight::get_data(model)
  if (is.null(dat)) {
    stop("insight::get_data(model) returned NULL - could not recover the ",
         "original data used to fit `model`.")
  }
  
  # --- sanitise columns carrying non-standard S3 classes on top of numeric --
  # (e.g. datawizard::standardize()'s "dw_transformer" class) which can break
  # strict type-checking code downstream, even though arithmetic on them
  # normally still works via S3 dispatch. Factors/characters are left as-is.
  dat[] <- lapply(dat, function(col) {
    if (is.numeric(col) && length(class(col)) > 1) as.numeric(col) else col
  })
  
  # --- build an equivalent plain lm(), reusing the EXACT RHS term structure
  #     of target_formula (interactions, I(), etc.) rather than re-deparsing
  #     it, to avoid any risk of re-parsing subtly changing the formula. The
  #     response is irrelevant to VIF (a property of X only), so a fresh
  #     random dummy response is substituted in directly. -------------------
  set.seed(seed)
  dat$.brms_vif_dummy_y <- stats::rnorm(nrow(dat))
  
  lm_formula <- target_formula
  lm_formula[[2]] <- as.name(".brms_vif_dummy_y")
  environment(lm_formula) <- environment()
  
  lm_fit <- tryCatch(
    stats::lm(lm_formula, data = dat),
    error = function(e) {
      stop(
        "Could not refit the equivalent lm() for this formula component. ",
        "This usually means the formula contains brms-specific syntax with ",
        "no plain-lm() equivalent (e.g. s(), t2(), gp(), mo(), me(), cs(), ",
        "or a non-linear bf() specification). VIF/collinearity is an ",
        "inherently linear/parametric-term concept - for smooth terms, use ",
        "mgcv::concurvity() on an equivalent mgcv::gam() fit instead. ",
        "Original error: ", conditionMessage(e)
      )
    }
  )
  
  # --- run car::vif(); type="terms" is the only option that generalises to
  #     non-lm/weighted fits, and avoids car's own unconditional "consider
  #     type='predictor'" message, which is never actually usable here (see
  #     the glmgee_vif() notes earlier this session - confirmed from car's
  #     source that type='predictor' is always overridden back to 'terms'
  #     for any weighted or non-plain-lm fit anyway). ------------------------
  vif_raw <- tryCatch(
    suppressMessages(car::vif(lm_fit, type = "terms")),
    error = function(e) {
      stop("car::vif() failed - this usually means the design matrix is ",
           "rank-deficient (aliased / perfectly collinear terms), which is ",
           "plausible given how severe some VIFs have been in this session. ",
           "Try alias(lm_fit_equivalent) to find the offending term(s). ",
           "Original error: ", conditionMessage(e))
    }
  )
  
  # --- normalise car's output (plain vector vs matrix) into one data frame --
  if (is.null(dim(vif_raw))) {
    tab <- data.frame(Term = names(vif_raw), GVIF = as.numeric(vif_raw), Df = 1L)
  } else {
    tab <- as.data.frame(vif_raw)
    tab$Term <- rownames(tab)
    rownames(tab) <- NULL
  }
  if (is.null(tab$Df)) tab$Df <- 1L
  
  # --- SE-inflation-scale adjusted VIF, comparable across all terms ---------
  # (reduces to sqrt(VIF) for ordinary 1-df terms; this is the scale on which
  # the conventional ~5 / ~10 severity thresholds apply, and is computed
  # CORRECTLY here per-term Df, unlike the performance::check_collinearity()
  # bug this function exists to route around)
  tab$Adjusted_VIF <- tab$GVIF ^ (1 / (2 * tab$Df))
  tab$Tolerance     <- 1 / tab$GVIF
  
  tab$Severity <- cut(
    tab$Adjusted_VIF,
    breaks = c(-Inf, thresholds[["moderate"]], thresholds[["high"]], Inf),
    labels = c("Low", "Moderate", "High")
  )
  
  tab <- tab[, c("Term", "Df", "GVIF", "Adjusted_VIF", "Tolerance", "Severity")]
  if (sort) tab <- tab[order(-tab$Adjusted_VIF), ]
  rownames(tab) <- NULL
  
  tab$GVIF         <- round(tab$GVIF, digits)
  tab$Adjusted_VIF <- round(tab$Adjusted_VIF, digits)
  tab$Tolerance    <- round(tab$Tolerance, digits + 1)
  
  cond <- .design_condition_diagnostics(stats::model.matrix(lm_fit))
  
  structure(tab, class = c("brms_vif", "data.frame"),
            dpar_checked = target_name,
            thresholds   = thresholds,
            n_high       = sum(tab$Severity == "High"),
            n_moderate   = sum(tab$Severity == "Moderate"),
            condition    = cond)
}

# --- simple print method for nicer console output ---------------------------
print.brms_vif <- function(x, ...) {
  th <- attr(x, "thresholds")
  cat(sprintf("Collinearity check for brms model (formula component checked: '%s')\n",
              attr(x, "dpar_checked")))
  cat("VIF via equivalent lm() on the design matrix; Df-adjusted scale\n")
  cat(strrep("-", 72), "\n")
  print.data.frame(x, row.names = FALSE)
  cat(strrep("-", 72), "\n")
  n_high <- attr(x, "n_high"); n_mod <- attr(x, "n_moderate")
  cat(sprintf("Severity thresholds on Adjusted_VIF: Moderate >= %s, High >= %s\n",
              th[["moderate"]], th[["high"]]))
  if (n_high > 0) cat(sprintf("-> %d term(s) at HIGH collinearity\n", n_high))
  if (n_mod  > 0) cat(sprintf("-> %d term(s) at MODERATE collinearity\n", n_mod))
  if (n_high == 0 && n_mod == 0) cat("-> No terms flagged for collinearity concern.\n")
  
  cond <- attr(x, "condition")
  cat(strrep("-", 72), "\n")
  cat("Design matrix condition number (Belsley, Kuh & Welsch 1980; complements\n")
  cat("the term-level VIFs above - reflects OVERALL, not per-term, conditioning)\n")
  cat(sprintf("  Scaled (recommended): %.1f   |   Raw kappa(X): %.1f\n",
              cond$condition_number_scaled, cond$condition_number_raw))
  cat(sprintf("  %d of %d condition indices exceed 30 (moderate); %d exceed 100 (severe)\n",
              cond$n_condition_indices_over_30, nrow(x) + 1L,  # +1 for intercept
              cond$n_condition_indices_over_100))
  if (cond$condition_number_scaled > 100) {
    cat("  -> SEVERE overall conditioning, despite the per-term VIFs above:\n")
    cat("     this reflects a diffuse combination across MANY terms jointly,\n")
    cat("     not any single problematic pair - term-level VIF cannot detect this.\n")
  } else if (cond$condition_number_scaled > 30) {
    cat("  -> Moderate overall conditioning.\n")
  }
  invisible(x)
}

# =============================================================================
# Example usage:
#
#   # main (mu / conditional) formula - e.g. mA1 (cumulative ordinal model)
#   brms_vif(mA1)
#
#   # an auxiliary-parameter formula, e.g. phi ~ ... on an ordered-beta model
#   #   (dpar name must match how it's shown in names(insight::find_formula(m1d)))
#   names(insight::find_formula(m1d))    # check available component names first
#   brms_vif(m1d, dpar = "phi")
#
#   # any family works identically - no family-specific code path needed:
#   brms_vif(m3b)     # ordbeta
#   brms_vif(mA1)      # cumulative
# =============================================================================