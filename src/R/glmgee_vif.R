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
# NOTE: not executed here (no R execution environment available) - please
# run and report back if anything errors or looks off.
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


glmgee_vif <- function(model,
                       digits = 2,
                       sort = TRUE,
                       thresholds = c(moderate = sqrt(5), high = sqrt(10)),
                       ci_threshold = 30,
                       vdp_threshold = 0.5) {
  
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
  # (reduces to sqrt(VIF) for ordinary 1-df terms; per Fox (2020, "Regression
  # Diagnostics", 2nd ed.) and car::vif()'s own documentation, this quantity
  # is on the same scale as sqrt(VIF), NOT raw VIF - the conventional 5/10
  # rule-of-thumb therefore applies at sqrt(5)/sqrt(10) on THIS column, which
  # is why those are the function's default thresholds, not 5/10 directly)
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
  
  cond <- .design_condition_diagnostics(stats::model.matrix(glm_fit),
                                        ci_threshold = ci_threshold,
                                        vdp_threshold = vdp_threshold)
  
  # NOTE: attributes set here (thresholds/n_high/n_moderate/condition) are
  # not guaranteed to survive generic data.frame operations like `[` or
  # rbind() on the returned object - if you subset/combine this result,
  # re-run print()/inspect the attributes on the ORIGINAL object, not a
  # derived subset of it.
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
#   m12 <- glmtoolbox::glmgee(formula = m12formula, data = m12Data, id = ID,
#                              family = binomial(link = "logit"),
#                              corstr = "exchangeable")
#
#   glmgee_vif(m12)
#
#   # full variance-decomposition-proportion matrix (variables x dimensions),
#   # e.g. to inspect loadings below the default 0.5 threshold:
#   attr(glmgee_vif(m12), "condition")$vdp
#
#   # adjust BKW's default conventions if needed:
#   glmgee_vif(m12, ci_threshold = 15, vdp_threshold = 0.4)
# =============================================================================