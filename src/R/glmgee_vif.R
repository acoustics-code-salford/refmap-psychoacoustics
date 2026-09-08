require(insight)
require(car)

# =============================================================================
# glmgee_vif(): collinearity diagnostics for glmtoolbox::glmgee models
#
# Rationale: collinearity is a property of the design matrix X (the formula +
# data), not of the working-correlation structure used to fit a GEE, so a
# valid and much simpler way to get VIFs for a glmgee fit is to refit its
# mean model as an ordinary glm() and run car::vif() on that - which is
# exactly what we did by hand earlier in this thread. This wraps that up as
# a single reusable function.
#
# Uses insight::find_formula() / get_data() / get_family() to pull the model
# components out of the glmgee object, since insight has explicit glmgee
# support (confirmed via its changelog). car::vif(..., type = "terms") is
# used directly (rather than the default type = "predictor") since "predictor"
# is only implemented for unweighted lm() and otherwise falls back to "terms"
# with a warning anyway - this skips that warning.
#
# NOTE: not executed here (no R access in this environment) - please run and
# report back if anything errors or looks off.
# =============================================================================

# -----------------------------------------------------------------------------
# Design-matrix condition number (Belsley, Kuh & Welsch 1980)
# -----------------------------------------------------------------------------
#
# WHY THIS IS SEPARATE FROM TERM-LEVEL VIF: VIF assesses collinearity
# term-by-term (how well can THIS term be predicted from the others); the
# condition number reflects the OVERALL design matrix's numerical
# conditioning, which can be poor because of a diffuse combination across
# MANY columns jointly even when no single term's VIF is high (confirmed
# directly this session: a related model's terms were all Adjusted_VIF < 4,
# yet kappa(X) ~ 500). The two diagnostics are complementary, not redundant -
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


glmgee_vif <- function(model,
                       digits = 2,
                       sort = TRUE,
                       thresholds = c(moderate = 5, high = 10)) {
  
  if (!inherits(model, "glmgee")) {
    stop("`model` must be a fitted glmtoolbox::glmgee object.")
  }
  for (pkg in c("insight", "car")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required but not installed.", pkg))
    }
  }
  
  # --- extract formula, data and family from the glmgee fit -----------------
  form <- insight::find_formula(model)$conditional
  dat  <- insight::get_data(model)
  fam  <- insight::get_family(model)
  
  if (is.null(form) || is.null(dat) || is.null(fam)) {
    stop("Could not extract formula/data/family from `model` via insight. ",
         "Check insight::find_formula(model), insight::get_data(model), and ",
         "insight::get_family(model) individually to see which one failed.")
  }
  
  # --- refit as an ordinary glm, purely for the design matrix / VIFs --------
  glm_fit <- tryCatch(
    stats::glm(form, data = dat, family = fam),
    error = function(e) {
      stop("Could not refit the equivalent glm(): ", conditionMessage(e))
    }
  )
  
  # --- run car::vif() ---------------------------------------------------
  # car::vif() unconditionally messages "consider setting type = 'predictor'"
  # whenever there are interaction terms and type = "terms" is used - but
  # type = "predictor" is always silently overridden back to "terms" for any
  # glm object regardless (confirmed from car's source, R/vif.R), so the
  # suggestion is never actually actionable here. Suppressed rather than left
  # to print on every call.
  vif_raw <- tryCatch(
    suppressMessages(car::vif(glm_fit, type = "terms")),
    error = function(e) {
      stop("car::vif() failed - this usually means the design matrix is ",
           "rank-deficient (aliased / perfectly collinear terms). Try ",
           "alias(glm(", deparse(form), ", data = <your data>, family = <family>)) ",
           "to find the offending term(s). Original error: ", conditionMessage(e))
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
  # the conventional ~5 / ~10 severity thresholds apply)
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
  
  cond <- .design_condition_diagnostics(stats::model.matrix(glm_fit))
  
  structure(tab, class = c("glmgee_vif", "data.frame"),
            thresholds = thresholds,
            n_high     = sum(tab$Severity == "High"),
            n_moderate = sum(tab$Severity == "Moderate"),
            condition  = cond)
}

# --- simple print method for nicer console output ---------------------------
print.glmgee_vif <- function(x, ...) {
  th <- attr(x, "thresholds")
  cat("Collinearity check for glmgee model (VIF via equivalent glm(); Df-adjusted scale)\n")
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
#   m12 <- glmtoolbox::glmgee(formula = m12formula, data = m12Data, id = ID,
#                              family = binomial(link = "logit"),
#                              corstr = "exchangeable")
#
#   glmgee_vif(m12)
#
# =============================================================================