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
# NOTE: not executed here (no R execution environment available) - please
# run and report back if anything errors or looks off, same as the other
# scripts this session.
# =============================================================================

# -----------------------------------------------------------------------------
# Design-matrix condition number + variance-decomposition proportions
# (Belsley, Kuh & Welsch 1980)
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
# THRESHOLDS (condition index, on the SCALED matrix): Belsley, Kuh & Welsch
# (1980) give a four-tier gradient - <=10: weak near-dependencies (no real
# concern); 10-30: moderately strong; 30-100: large/strong near-dependency;
# >100: severe. The same gradient applies to the largest (overall) condition
# number, which is simply the largest condition index.
#
# VARIANCE-DECOMPOSITION PROPORTIONS (the second half of BKW's procedure -
# a high condition index ALONE cannot identify which variables are involved,
# nor distinguish a real problem from an incidentally-large index loading on
# only one variable): for each dimension k, the proportion of variable j's
# coefficient-variance attributable to that dimension is
#   pi[j,k] = phi[j,k] / sum_k(phi[j,k]),  phi[j,k] = V[j,k]^2 / d[k]^2
# where V holds the right singular vectors and d the singular values of the
# scaled matrix. A dimension with condition index above ci_threshold (30 by
# default) on which TWO OR MORE variables have pi > vdp_threshold (0.5 by
# default, per BKW) indicates those specific variables are collinear with
# each other via that dimension. A single high-proportion variable on a
# flagged dimension is not, by itself, evidence of a collinearity problem.
.design_condition_diagnostics <- function(X, ci_threshold = 30, vdp_threshold = 0.5) {
  norms <- sqrt(colSums(X^2))
  norms[norms == 0] <- 1  # guard against a degenerate constant/zero column
  X_scaled <- sweep(X, 2, norms, "/")
  
  sv_decomp <- svd(X_scaled)
  d <- sv_decomp$d
  V <- sv_decomp$v
  
  keep <- d > .Machine$double.eps
  if (!all(keep)) {
    warning(sum(!keep), " near-zero singular value(s) dropped (the design ",
            "matrix is exactly rank-deficient) - condition indices and ",
            "variance-decomposition proportions reflect only the ",
            "non-degenerate dimensions; investigate via alias() first.")
  }
  d <- d[keep]
  V <- V[, keep, drop = FALSE]
  
  cond_indices <- max(d) / d
  
  # --- variance-decomposition proportions -----------------------------------
  phi        <- sweep(V^2, 2, d^2, "/")   # phi[j,k]: rows = variables, cols = dimensions
  row_totals <- rowSums(phi)
  row_totals[row_totals == 0] <- 1        # guard against a degenerate all-zero row
  pi_mat <- sweep(phi, 1, row_totals, "/")
  rownames(pi_mat) <- colnames(X)
  colnames(pi_mat) <- paste0("dim", seq_along(d))
  
  # --- flag dimensions with a high condition index AND >=2 implicated vars --
  flagged_dims <- which(cond_indices > ci_threshold)
  problem_groups <- lapply(flagged_dims, function(k) {
    vars <- rownames(pi_mat)[pi_mat[, k] > vdp_threshold]
    if (length(vars) >= 2) {
      list(dimension = k, condition_index = cond_indices[k], variables = vars)
    } else {
      NULL
    }
  })
  problem_groups <- Filter(Negate(is.null), problem_groups)
  
  list(
    condition_number_scaled   = max(cond_indices),
    condition_number_raw      = tryCatch(kappa(X, exact = TRUE), error = function(e) NA_real_),
    condition_indices         = cond_indices,
    n_condition_indices_total = length(cond_indices),
    n_ci_10_30                = sum(cond_indices > 10  & cond_indices <= 30),
    n_ci_30_100                = sum(cond_indices > 30  & cond_indices <= 100),
    n_ci_over_100              = sum(cond_indices > 100),
    vdp            = pi_mat,
    problem_groups = problem_groups,
    vdp_threshold  = vdp_threshold,
    ci_threshold   = ci_threshold
  )
}


brms_vif <- function(model,
                     dpar = NULL,
                     digits = 2,
                     sort = TRUE,
                     thresholds = c(moderate = sqrt(5), high = sqrt(10)),
                     ci_threshold = 30,
                     vdp_threshold = 0.5,
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
  #
  # RNG NOTE: set.seed() mutates the GLOBAL RNG state, which would otherwise
  # silently affect anything relying on random draws later in the SAME R
  # session (e.g. a subsequent brm(..., seed = ...) call's actual sampling
  # path, or any other stochastic step) - this save/restore keeps the dummy
  # response reproducible internally without leaking that side effect out.
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  on.exit({
    if (!is.null(old_seed)) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
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
  # (reduces to sqrt(VIF) for ordinary 1-df terms; per Fox (2020, "Regression
  # Diagnostics", 2nd ed.) and car::vif()'s own documentation, this quantity
  # is on the same scale as sqrt(VIF), NOT raw VIF - the conventional 5/10
  # rule-of-thumb therefore applies at sqrt(5)/sqrt(10) on THIS column, which
  # is why those are the function's default thresholds, not 5/10 directly.
  # Computed CORRECTLY here per-term Df, unlike the
  # performance::check_collinearity() bug this function exists to route
  # around.)
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
  
  cond <- .design_condition_diagnostics(stats::model.matrix(lm_fit),
                                        ci_threshold = ci_threshold,
                                        vdp_threshold = vdp_threshold)
  
  # NOTE: attributes set here (thresholds/n_high/n_moderate/condition) are
  # not guaranteed to survive generic data.frame operations like `[` or
  # rbind() on the returned object - if you subset/combine this result,
  # re-run print()/inspect the attributes on the ORIGINAL object, not a
  # derived subset of it.
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
  cat(sprintf("Severity thresholds on Adjusted_VIF: Moderate >= %.2f, High >= %.2f\n",
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
  cat(sprintf("  Condition indices (%d total): %d in (10,30] moderate | %d in (30,100] strong | %d > 100 severe\n",
              cond$n_condition_indices_total, cond$n_ci_10_30, cond$n_ci_30_100, cond$n_ci_over_100))
  
  cn <- cond$condition_number_scaled
  overall <- if (cn <= 10)       "Weak/no meaningful ill-conditioning."
  else if (cn <= 30)  "Moderate overall conditioning."
  else if (cn <= 100) "Strong overall conditioning."
  else                "SEVERE overall conditioning."
  cat("  ->", overall, "\n")
  
  if (length(cond$problem_groups) > 0) {
    cat(strrep("-", 72), "\n")
    cat(sprintf("Variance-decomposition proportions (Belsley, Kuh & Welsch 1980):\n"))
    cat(sprintf("dimensions with condition index > %g AND >= 2 variables with proportion > %.2f:\n",
                cond$ci_threshold, cond$vdp_threshold))
    for (grp in cond$problem_groups) {
      cat(sprintf("  Dimension (condition index = %.1f): %s\n",
                  grp$condition_index, paste(grp$variables, collapse = ", ")))
    }
  } else if (cond$n_ci_30_100 + cond$n_ci_over_100 > 0) {
    cat(strrep("-", 72), "\n")
    cat(sprintf("Note: %d dimension(s) have condition index > %g, but none show >= 2\n",
                cond$n_ci_30_100 + cond$n_ci_over_100, cond$ci_threshold))
    cat("variables with variance-decomposition proportion > 0.5 - i.e. no specific\n")
    cat("variable pair is clearly implicated per Belsley, Kuh & Welsch's criterion,\n")
    cat("despite the elevated condition index(es). Full proportions are in\n")
    cat("attr(x, \"condition\")$vdp if you want to inspect below-threshold loadings.\n")
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
#
#   # full variance-decomposition-proportion matrix (variables x dimensions):
#   attr(brms_vif(mA6), "condition")$vdp
#
#   # adjust BKW's default conventions if needed:
#   brms_vif(mA6, ci_threshold = 15, vdp_threshold = 0.4)
# =============================================================================