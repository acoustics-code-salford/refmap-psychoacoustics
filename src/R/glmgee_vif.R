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
  
  structure(tab, class = c("glmgee_vif", "data.frame"),
            thresholds = thresholds,
            n_high     = sum(tab$Severity == "High"),
            n_moderate = sum(tab$Severity == "Moderate"))
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