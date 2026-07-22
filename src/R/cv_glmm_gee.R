# ============================================================================
# VERSION STAMP -- prints on source() so you can confirm at a glance that
# what's loaded matches what you think you just saved, rather than a stale
# copy from an earlier session. Update the date string whenever you save a
# new version of this file.
# ============================================================================
message("cv_glmmTMB.R loaded -- version 2026-07-22")

# ============================================================================
# cv_glmmTMB(): grouped K-fold cross-validation for glmmTMB models, returning
# out-of-sample, population-averaged (marginal) predictive performance.
#
# Families with built-in log predictive density support (for `elpd`): ordbeta,
# binomial, poisson, gaussian. RMSE / MAE / R2 are family-agnostic (they only
# need predicted mu vs observed y) and work for ANY glmmTMB family.
#
# Cutpoints from family_params() for ordbeta are treated as being on the
# response (0,1) probability scale and qlogis()-transformed before use --
# confirmed against glmmTMB 1.1.14 via a training-data logLik sanity check
# (see cv_glmmTMB_sanity_check() below). Re-verify on other versions.
#
# Per-row dispersion (dispformula with predictors) is handled via
# predict(fit, type = "disp"), the dispformula analogue of predicting mu.
# Fold refits now carry the ORIGINAL model's dispformula/ziformula/weights/
# control forward -- earlier versions of this function silently refit every
# fold with dispformula = ~1 regardless of what the input model used.
#
# To support a family not listed here, add an entry to `cv_glmmTMB_families`:
#   cv_glmmTMB_families$Gamma <- function(y, mu, fit, newdata) { ... }
# Signature: (y, mu, fit, newdata) -> numeric vector of per-observation log
# predictive densities. Or pass `loglik_fun` directly to cv_glmmTMB() for a
# one-off override without touching the registry.
# ============================================================================

cv_glmmTMB_families <- list(
  
  ordbeta = function(y, mu, fit, newdata) {
    if (!requireNamespace("ordbetareg", quietly = TRUE)) {
      stop("Package 'ordbetareg' is required for elpd with family = ordbeta().")
    }
    phi <- tryCatch(
      stats::predict(fit, newdata = newdata, type = "disp",
                     re.form = NA, allow.new.levels = TRUE),
      error = function(e) NULL
    )
    if (is.null(phi)) phi <- glmmTMB::sigma(fit)
    cutoffs <- qlogis(unname(glmmTMB::family_params(fit)))
    ordbetareg::dordbeta(y, mu = mu, phi = phi, cutpoints = cutoffs, log = TRUE)
  },
  
  binomial = function(y, mu, fit, newdata) {
    wvar <- tryCatch(insight::find_weights(fit), error = function(e) NULL)
    if (!is.null(wvar) && wvar %in% names(newdata)) {
      trials <- newdata[[wvar]]
    } else {
      trials <- rep(1, length(y))  # assumes binary 0/1 response
    }
    successes <- round(y * trials)
    stats::dbinom(successes, size = trials, prob = mu, log = TRUE)
  },
  
  poisson = function(y, mu, fit, newdata) {
    stats::dpois(y, lambda = mu, log = TRUE)
  },
  
  gaussian = function(y, mu, fit, newdata) {
    sd_hat <- tryCatch(
      stats::predict(fit, newdata = newdata, type = "disp",
                     re.form = NA, allow.new.levels = TRUE),
      error = function(e) NULL
    )
    if (is.null(sd_hat)) sd_hat <- glmmTMB::sigma(fit)
    stats::dnorm(y, mean = mu, sd = sd_hat, log = TRUE)
  }
)

# ----------------------------------------------------------------------------
# One-off check comparing sum(loglik_fun(...)) on the TRAINING data against
# logLik(model). Expect the same ballpark, not an exact match -- this check
# uses re.form = NA (marginal, b = 0) while logLik() integrates over the
# Laplace-approximated conditional random-effect modes. Its real purpose is
# to catch scale errors (e.g. an implausible or wrong-signed total), which
# tend to be obvious even without an exact match.
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# Default marginal-prediction engine: stats::predict() with re.form = NA
# (population-averaged) and allow.new.levels = TRUE (handles validation/new
# IDs unseen in training). Override via `predict_fn` in cv_glmmTMB() and/or
# predict_marginal_interval() if you want both to instead go through
# marginaleffects (or anything else) -- e.g.:
#
#   mae_predict_fn <- function(fit, newdata) {
#     marginaleffects::predictions(fit, newdata = newdata, re.form = NA,
#                                   type = "response",
#                                   allow.new.levels = TRUE)$estimate
#   }
#
# Whichever you choose, use the SAME predict_fn in both cv_glmmTMB() (where
# the conformal margin is calibrated) and predict_marginal_interval() (where
# it's applied) -- the margin is only valid for the prediction-generating
# process it was calibrated against.
# ----------------------------------------------------------------------------
.default_predict_fn <- function(fit, newdata) {
  stats::predict(fit, newdata = newdata, type = "response",
                 re.form = NA, allow.new.levels = TRUE)
}

cv_glmmTMB_sanity_check <- function(model, loglik_fun = NULL) {
  fam_name <- family(model)$family
  if (is.null(loglik_fun)) {
    loglik_fun <- cv_glmmTMB_families[[fam_name]]
    if (is.null(loglik_fun)) {
      stop("No built-in loglik_fun for family '", fam_name, "'; supply one via `loglik_fun`.")
    }
  }
  data <- insight::get_data(model)
  resp <- insight::find_response(model)
  mu   <- stats::predict(model, type = "response")
  y    <- data[[resp]]
  
  ll <- sum(loglik_fun(y, mu, model, data), na.rm = TRUE)
  cat("sum(loglik_fun(...)) on training data:", ll, "\n")
  cat("logLik(model):                        ", as.numeric(stats::logLik(model)), "\n")
  cat("(expect same ballpark, not exact -- marginal vs conditional-mode likelihood)\n")
  invisible(list(loglik_sum = ll, logLik = as.numeric(stats::logLik(model))))
}

# ============================================================================
cv_glmmTMB <- function(model,
                       K            = 6,
                       id_col       = NULL,
                       metrics      = c("elpd", "rmse", "mae", "r2"),
                       se           = FALSE,
                       seed         = 123,
                       workers      = NULL,
                       control      = NULL,
                       loglik_fun   = NULL,
                       backtransform_fn = NULL,
                       return_predictions = FALSE,
                       predict_fn   = NULL,
                       agg_cols     = NULL,
                       data         = NULL) {
  
  stopifnot(inherits(model, "glmmTMB"))
  
  if (!requireNamespace("insight", quietly = TRUE)) {
    stop("Package 'insight' is required (used to extract formula/data/grouping factor).")
  }
  
  fam_name <- family(model)$family
  
  # ---- extract everything needed from the model object itself -------------
  form         <- stats::formula(model)
  disp_form    <- tryCatch(stats::formula(model, component = "disp"), error = function(e) ~1)
  zi_form      <- tryCatch(stats::formula(model, component = "zi"),   error = function(e) ~0)
  # insight::get_data() only reconstructs FORMULA variables from the model's
  # internal frame -- columns bound on for convenience (e.g. StimID, not used
  # in the model) won't survive that route. Supply `data` directly (e.g. the
  # same cbind(m1Data, StimID = ...) frame you already build elsewhere) if you
  # need agg_cols that aren't actual predictors.
  if (is.null(data)) data <- insight::get_data(model)
  resp         <- insight::find_response(model)
  weight_col   <- tryCatch(insight::find_weights(model), error = function(e) NULL)
  
  if (!is.null(agg_cols) && !all(agg_cols %in% names(data))) {
    stop("agg_cols ", paste(setdiff(agg_cols, names(data)), collapse = ", "),
         " not found in `data`. If these aren't actual predictors in the model (e.g. a StimID ",
         "label), pass your own data frame via cv_glmmTMB(..., data = your_full_data_with_StimID) ",
         "instead of relying on insight::get_data(model), which only retains formula variables. ",
         "Whatever `data` is used (supplied or auto-extracted), agg_cols must be columns in it, ",
         "correctly row-aligned with everything else in that same data frame -- cv_glmmTMB()'s ",
         "internal na.omit() filtering and fold construction rely on that alignment, which is why ",
         "agg_cols isn't accepted as a separately-supplied vector.")
  }
  
  if (is.null(id_col)) {
    rand_grp <- insight::find_random(model, flatten = TRUE)
    if (length(rand_grp) == 0) {
      stop("No random-effect grouping factor found on `model`; pass id_col explicitly.")
    }
    if (length(rand_grp) > 1) {
      message("Multiple grouping factors found (", paste(rand_grp, collapse = ", "),
              "); using '", rand_grp[1], "' for fold construction. Pass id_col to override.")
    }
    id_col <- rand_grp[1]
  }
  
  # ---- density function dispatch for elpd -----------------------------------
  if ("elpd" %in% metrics && is.null(loglik_fun)) {
    builtin <- cv_glmmTMB_families[[fam_name]]
    if (is.null(builtin)) {
      warning("elpd requested but no built-in log-density for family '", fam_name,
              "'. Dropping 'elpd' from metrics. Supply `loglik_fun(y, mu, fit, newdata)`, ",
              "or add an entry to cv_glmmTMB_families, to enable it.")
      metrics <- setdiff(metrics, "elpd")
    } else {
      loglik_fun <- builtin
    }
  }
  
  # ---- clean data & build grouped folds -------------------------------------
  keep_cols  <- unique(c(all.vars(form), id_col, weight_col, agg_cols))
  clean_data <- stats::na.omit(data[, keep_cols])
  
  set.seed(seed)
  unique_ids <- unique(clean_data[[id_col]])
  id_folds   <- sample(rep(seq_len(K), length.out = length(unique_ids)))
  names(id_folds) <- unique_ids
  data_folds <- id_folds[as.character(clean_data[[id_col]])]
  
  # ---- parallel setup ---------------------------------------------------
  if (is.null(workers)) workers <- max(1, parallelly::availableCores() - 1)
  old_plan <- future::plan(future::multisession, workers = workers)
  on.exit(future::plan(old_plan), add = TRUE)
  
  # ---- per-fold fit + marginal prediction -----------------------------------
  fold_out <- future.apply::future_lapply(seq_len(K), function(k) {
    
    library(glmmTMB)
    
    train <- clean_data[data_folds != k, ]
    val   <- clean_data[data_folds == k, ]
    
    # carry the ORIGINAL model's full specification forward into each fold --
    # dispformula, ziformula, weights, and control all matter for correctness,
    # not just the conditional formula.
    fit_args <- list(
      formula     = form,
      data        = train,
      family      = family(model),
      dispformula = disp_form,
      ziformula   = zi_form,
      se          = FALSE
    )
    if (!is.null(weight_col))  fit_args$weights <- train[[weight_col]]
    if (!is.null(control))     fit_args$control <- control
    
    fit_err <- NULL
    fit <- tryCatch(do.call(glmmTMB::glmmTMB, fit_args),
                    error = function(e) { fit_err <<- conditionMessage(e); NULL })
    if (is.null(fit)) return(list(.failed = TRUE, .error = paste("fit:", fit_err)))
    
    # population-averaged (marginal) predictions -- see .default_predict_fn
    # above, or your own predict_fn, for exactly what this computes.
    pred_fn  <- if (is.null(predict_fn)) .default_predict_fn else predict_fn
    pred_err <- NULL
    mu <- tryCatch(pred_fn(fit, val),
                   error = function(e) { pred_err <<- conditionMessage(e); NULL })
    if (is.null(mu)) return(list(.failed = TRUE, .error = paste("predict:", pred_err)))
    
    y <- val[[resp]]
    
    out <- list(y = y, mu = mu, .failed = FALSE)
    if (!is.null(agg_cols)) out$agg <- val[, agg_cols, drop = FALSE]
    
    if ("elpd" %in% metrics) {
      out$ll <- tryCatch(loglik_fun(y, mu, fit, val), error = function(e) rep(NA_real_, length(y)))
    }
    
    out
  }, future.seed = TRUE)
  
  is_failure  <- vapply(fold_out, function(f) isTRUE(f$.failed), logical(1))
  fold_errors <- lapply(fold_out[is_failure], function(f) f$.error)
  fold_out    <- fold_out[!is_failure]
  n_fail      <- sum(is_failure)
  if (n_fail > 0) {
    message(n_fail, " of ", K, " folds failed to fit/predict and were dropped.")
    if (length(fold_errors) > 0) message("First error: ", fold_errors[[1]])
  }
  if (length(fold_out) == 0) stop("All folds failed.",
                                  if (length(fold_errors) > 0) paste0(" First error: ", fold_errors[[1]]) else "")
  
  # ---- pool observation-level results across folds --------------------------
  y_all  <- unlist(lapply(fold_out, `[[`, "y"),  use.names = FALSE)
  mu_all <- unlist(lapply(fold_out, `[[`, "mu"), use.names = FALSE)
  resid  <- y_all - mu_all
  n_obs  <- length(y_all)
  
  result <- list(K_requested = K, K_used = length(fold_out), n_obs = n_obs)
  
  # ELPD (sum of pointwise log predictive density, same convention as loo::loo())
  if ("elpd" %in% metrics) {
    ll_all <- unlist(lapply(fold_out, `[[`, "ll"), use.names = FALSE)
    result$elpd <- sum(ll_all, na.rm = TRUE)
    if (se) {
      # Fold-level SE, not the pointwise sqrt(n*var(ll)) convention from loo::loo():
      # that convention treats pointwise contributions as independent, which is a
      # much weaker assumption here than in single-model LOO -- residuals within a
      # K-fold CV fold share that fold's specific fitted-model estimation error.
      # Folds are only approximately independent themselves (overlapping training
      # data across folds), but that's the more defensible unit available here.
      fold_elpd <- vapply(fold_out, function(f) sum(f$ll, na.rm = TRUE), numeric(1))
      result$elpd_se <- stats::sd(fold_elpd, na.rm = TRUE) / sqrt(length(fold_elpd)) * K
      # scaled by K: fold_elpd sums are per-fold totals on unequal-sized folds: this
      # SE is on the same total-sum scale as `elpd` itself, not a per-observation scale.
    }
  }
  
  # RMSE
  fold_rmse <- NULL
  if ("rmse" %in% metrics || se) {
    fold_rmse <- vapply(fold_out, function(f) sqrt(mean((f$y - f$mu)^2, na.rm = TRUE)), numeric(1))
  }
  if ("rmse" %in% metrics) {
    result$mse  <- mean(resid^2, na.rm = TRUE)   # = CV Brier score when y_all is binary (0/1)
    result$rmse <- sqrt(result$mse)
    if (se) {
      # Fold-level SE, not sd(resid^2)/sqrt(n_obs): residuals within a fold are
      # correlated (they share that fold's fitted-model estimation error), so
      # treating all n_obs residuals as independent understates the true SE.
      result$rmse_se <- stats::sd(fold_rmse, na.rm = TRUE) / sqrt(length(fold_rmse))
    }
  }
  
  # MAE
  if ("mae" %in% metrics) {
    result$mae <- mean(abs(resid), na.rm = TRUE)
    if (se) {
      fold_mae <- vapply(fold_out, function(f) mean(abs(f$y - f$mu), na.rm = TRUE), numeric(1))
      result$mae_se <- stats::sd(fold_mae, na.rm = TRUE) / sqrt(length(fold_mae))
    }
  }
  
  # CV-R2: two variants -- Allen (1974) / Golbraikh & Tropsha (2002) style
  # PRESS-R2, and Gelman et al. (2019) variance-ratio R2 adapted from LOO to CV.
  if ("r2" %in% metrics) {
    result$r2_cv        <- 1 - sum(resid^2, na.rm = TRUE) / sum((y_all - mean(y_all))^2, na.rm = TRUE)
    result$r2_cv_gelman <- stats::var(mu_all) / (stats::var(mu_all) + stats::var(resid))
    
    if (se) {
      fold_r2 <- vapply(fold_out, function(f) {
        1 - sum((f$y - f$mu)^2) / sum((f$y - mean(f$y))^2)
      }, numeric(1))
      result$r2_cv_se <- stats::sd(fold_r2, na.rm = TRUE) / sqrt(length(fold_r2))
    }
  }
  
  # ---- optional back-transformed (original-scale) versions ------------------
  # R2 is affine-invariant so *_raw will equal the transformed-scale value for
  # a linear backtransform_fn -- included anyway for a nonlinear backtransform,
  # and so all reported metrics are consistently available on both scales.
  if (!is.null(backtransform_fn)) {
    y_raw    <- backtransform_fn(y_all)
    mu_raw   <- backtransform_fn(mu_all)
    resid_raw <- y_raw - mu_raw
    
    if ("rmse" %in% metrics) {
      mse_raw <- mean(resid_raw^2, na.rm = TRUE)
      result$rmse_raw <- sqrt(mse_raw)
      if (se) {
        mse_raw_se <- stats::sd(resid_raw^2, na.rm = TRUE) / sqrt(n_obs)
        result$rmse_raw_se <- mse_raw_se / (2 * result$rmse_raw)
      }
    }
    if ("mae" %in% metrics) {
      result$mae_raw <- mean(abs(resid_raw), na.rm = TRUE)
      if (se) result$mae_raw_se <- stats::sd(abs(resid_raw), na.rm = TRUE) / sqrt(n_obs)
    }
    if ("r2" %in% metrics) {
      result$r2_cv_raw        <- 1 - sum(resid_raw^2, na.rm = TRUE) / sum((y_raw - mean(y_raw))^2, na.rm = TRUE)
      result$r2_cv_gelman_raw <- stats::var(mu_raw) / (stats::var(mu_raw) + stats::var(resid_raw))
    }
  }
  
  if (return_predictions) {
    # native model scale (e.g. [0,1] for ordbeta) -- deliberately NOT
    # back-transformed here; see cv_conformal_margin()/predict_marginal_interval()
    # for why endpoints, not raw residuals, are what get back-transformed.
    result$y_all  <- y_all
    result$mu_all <- mu_all
    if (!is.null(agg_cols)) {
      result$agg_data <- do.call(rbind, lapply(fold_out, `[[`, "agg"))
    }
  }
  
  result
}

# ============================================================================
# Builds a counterfactual grid matching what marginaleffects::avg_predictions(
# variables = ...) constructs internally: for each combination of values in
# `variables`, a full copy of `base_data` with those columns overridden,
# leaving every OTHER column (including id_col and all other covariates) at
# each row's own OBSERVED value. This is what makes it "counterfactual"
# rather than a synthetic reference grid: real subjects, real other-covariate
# values, only the swept variable(s) changed.
#
# variables: named list, e.g.
#   list(UASLAEMaxLRScl = seq(6, 18, length.out = 25))                 -- sweep
#   list(UASLAEMaxLRScl = c(6, 12), UASType = c("H520", "T150"))       -- crossed
#
# Multiple variables are CROSSED (full Cartesian product via expand.grid),
# matching marginaleffects' default (variables = list(...) without cross =
# FALSE). If you're using cross = FALSE (matched pairs, not all combinations)
# on the marginaleffects side, this needs a different construction -- let me
# know and I'll adjust.
#
# The returned grid identifies each combination via the swept variable
# column(s) themselves (e.g. "UASLAEMaxLRScl", or "UASLAEMaxLRScl" +
# "UASType" together) -- pass names(variables) as agg_cols to
# simulate_marginal_ci() or cf_variables to cluster_boot_counterfactual().
# ============================================================================

build_counterfactual_grid <- function(base_data, variables) {
  grid <- do.call(expand.grid, c(variables, list(stringsAsFactors = FALSE)))
  do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    d <- base_data
    for (v in names(variables)) d[[v]] <- grid[i, v]
    d
  }))
}


# cv_glmmTMB(..., agg_cols = ..., return_predictions = TRUE) output.
#
# This is the cross-validated analogue of the classic "predicted vs observed
# aggregate" scatterplot (e.g. mean change in annoyance, predicted vs measured,
# per stimulus) -- but using genuinely out-of-sample predictions rather than
# a final model's in-sample fitted values. Because fold assignment is BY
# SUBJECT, not by stimulus, every subject who saw a given stimulus ends up
# contributing an out-of-sample prediction for it from whichever fold held
# THEM out (different subjects seeing the same stimulus are typically spread
# across multiple folds) -- so pooling across all K folds and then
# aggregating by stimulus reconstructs the full per-stimulus average using
# only out-of-sample individual predictions, directly comparable to the
# original in-sample version of the same plot.
# ============================================================================

cv_aggregated_metrics <- function(cv_result, backtransform_fn = NULL, bind = TRUE) {
  if (is.null(cv_result$agg_data) || is.null(cv_result$y_all) || is.null(cv_result$mu_all)) {
    stop("cv_result must come from cv_glmmTMB(..., agg_cols = ..., return_predictions = TRUE).")
  }
  
  df <- cv_result$agg_data
  df$y  <- cv_result$y_all
  df$mu <- cv_result$mu_all
  
  grp_cols <- setdiff(names(df), c("y", "mu"))
  grp_key  <- interaction(df[grp_cols], drop = TRUE)
  
  agg <- do.call(rbind, lapply(split(df, grp_key), function(sub) {
    data.frame(sub[1, grp_cols, drop = FALSE],
               y_mean  = mean(sub$y),
               mu_mean = mean(sub$mu),
               n       = nrow(sub),   # subjects contributing to this stimulus's average
               row.names = NULL)
  }))
  
  resid <- agg$y_mean - agg$mu_mean
  agg_result <- list(
    n_groups     = nrow(agg),
    mse_agg      = mean(resid^2),   # NOT "Brier score" -- y_mean is a group-averaged
    # proportion here, not a single binary draw; that
    # distinction only holds at the row level (see `mse`
    # in cv_glmgee()/cv_glmmTMB()'s own output).
    rmse_agg     = sqrt(mean(resid^2)),
    mae_agg      = mean(abs(resid)),
    r2_agg       = 1 - sum(resid^2) / sum((agg$y_mean - mean(agg$y_mean))^2),
    agg_summary  = agg   # renamed from `data` when this was a standalone return, to avoid
    # colliding/confusing with cv_result's own `agg_data` field once bound
  )
  
  if (!is.null(backtransform_fn)) {
    y_raw     <- backtransform_fn(agg$y_mean)
    mu_raw    <- backtransform_fn(agg$mu_mean)
    resid_raw <- y_raw - mu_raw
    agg_result$rmse_agg_raw <- sqrt(mean(resid_raw^2))
    agg_result$mae_agg_raw  <- mean(abs(resid_raw))
    agg_result$r2_agg_raw   <- 1 - sum(resid_raw^2) / sum((y_raw - mean(y_raw))^2)
    agg_result$agg_summary$y_mean_raw  <- y_raw
    agg_result$agg_summary$mu_mean_raw <- mu_raw
  }
  
  if (!bind) return(agg_result)   # old standalone-list behavior, for backward compatibility
  
  # default: merge into a COPY of cv_result -- one combined object per model,
  # ready for compare_cv_results() below. `data`/`agg_data` fields renamed as
  # above; everything else keeps its existing name (no collisions: cv_result's
  # row-level fields are rmse/mae/r2_cv/..., these are rmse_agg/mae_agg/r2_agg/...).
  result <- cv_result
  result[names(agg_result)] <- agg_result
  result
}

# ----------------------------------------------------------------------------
# Combines multiple cv_glmmTMB()/cv_glmgee() results (ideally already run
# through cv_aggregated_metrics(..., bind = TRUE) so both row-level and
# aggregate metrics are present) into ONE data frame, one row per model, for
# side-by-side comparison.
#
# Only known SCALAR metric fields are extracted -- non-scalar fields
# (y_all, mu_all, agg_data, agg_summary) don't belong in a comparison table
# and are skipped. A field missing from a given model's result (e.g. elpd
# for a glmgee model, which has no valid likelihood -- see cv_glmgee's own
# header comment) is filled with NA rather than erroring, so glmmTMB and
# glmgee results can sit in the same table.
# ----------------------------------------------------------------------------
compare_cv_results <- function(cv_results, model_names = NULL) {
  if (is.null(model_names)) model_names <- names(cv_results)
  if (is.null(model_names)) model_names <- paste0("model", seq_along(cv_results))
  
  scalar_fields <- c(
    "K_requested", "K_used", "n_obs",
    "elpd", "elpd_se",
    "mse", "rmse", "rmse_se", "mae", "mae_se",
    "r2_cv", "r2_cv_se", "r2_cv_gelman",
    "rmse_raw", "rmse_raw_se", "mae_raw", "mae_raw_se",
    "r2_cv_raw", "r2_cv_gelman_raw",
    "n_groups", "mse_agg", "rmse_agg", "mae_agg", "r2_agg",
    "rmse_agg_raw", "mae_agg_raw", "r2_agg_raw"
  )
  
  rows <- lapply(seq_along(cv_results), function(i) {
    r <- cv_results[[i]]
    vals <- lapply(scalar_fields, function(f) {
      v <- r[[f]]
      if (is.null(v) || length(v) != 1) NA_real_ else v
    })
    names(vals) <- scalar_fields
    data.frame(model = model_names[i], as.data.frame(vals, check.names = FALSE),
               check.names = FALSE, row.names = NULL)
  })
  
  out <- do.call(rbind, rows)
  # drop columns that are NA for every model (e.g. elpd_se if se=FALSE everywhere)
  out[, c(TRUE, !vapply(out[-1], function(col) all(is.na(col)), logical(1))), drop = FALSE]
}


#
# Distribution-free (Vovk, Gammerman & Shafer, 2005; K-fold/CV+ extension:
# Barber, Candes, Ramdas & Tibshirani, 2021, "Predictive inference with the
# jackknife+", Annals of Statistics). Uses the pooled out-of-fold residuals
# already computed by cv_glmmTMB(..., return_predictions = TRUE) -- no extra
# model fits or simulation required.
#
# This is the simpler "pooled K-fold conformal" variant (one global margin
# from all out-of-fold residuals), not full CV+/jackknife+ (which combines
# per-point leave-one-fold-out bounds individually for a tighter theoretical
# guarantee). Coverage in practice is usually close to nominal for K >= 5-10;
# swap in the full CV+ combination step later if you need the stronger
# finite-sample guarantee.
#
# Suited to continuous / bounded-continuous / count outcomes (ordbeta,
# gaussian, poisson, aggregated binomial proportions). NOT appropriate for
# strictly binary 0/1 response -- that needs conformal prediction SETS, a
# different construction.
# ============================================================================

cv_conformal_margin <- function(cv_result, level = 0.95) {
  if (is.null(cv_result$y_all) || is.null(cv_result$mu_all)) {
    stop("cv_result must come from cv_glmmTMB(..., return_predictions = TRUE).")
  }
  scores <- abs(cv_result$y_all - cv_result$mu_all)
  n      <- length(scores)
  alpha  <- 1 - level
  
  # finite-sample conformal quantile correction (Barber et al. 2021)
  idx    <- min(n, ceiling((1 - alpha) * (n + 1)))
  margin <- sort(scores)[idx]
  
  empirical_coverage <- mean(scores <= margin)
  
  list(margin = margin, level = level, n = n,
       empirical_coverage = empirical_coverage)
}

# ----------------------------------------------------------------------------
# Aggregate-level conformal margin, for prediction intervals around
# AGGREGATE predictions (e.g. a counterfactual/focal grid via
# avg_predictions(variables = ..., by = ...)), NOT individual-row
# predictions. Calibrates on the aggregate residuals (y_mean - mu_mean per
# group) already computed by cv_aggregated_metrics(), rather than pooled
# individual-row residuals -- using the individual-row margin
# (cv_conformal_margin()) directly on an aggregate prediction would
# overstate its uncertainty, since it doesn't reflect the variance reduction
# from averaging over many subjects.
#
# Two things worth being aware of, inherent to calibrating at this level:
# - `n` here is the number of GROUPS (e.g. 94 stimuli), not individual
#   observations -- a much smaller effective calibration sample, so this
#   margin is noisier than the individual-row one for the same data.
# - The conformal guarantee is tied to exchangeability of the CALIBRATION
#   groups' effective size (how many subjects went into each aggregate).
#   Applying this margin to a counterfactual aggregate built from a group of
#   a very different effective size (e.g. a small subset rather than the
#   full ~41-42 subjects per stimulus) stretches that assumption -- fine if
#   your counterfactual grid is built from the full dataset (typical group
#   sizes should roughly match), worth flagging if built from a subset.
# ----------------------------------------------------------------------------
cv_aggregated_conformal_margin <- function(agg_result, level = 0.95) {
  if (is.null(agg_result$data) || !all(c("y_mean", "mu_mean") %in% names(agg_result$data))) {
    stop("agg_result must come from cv_aggregated_metrics().")
  }
  scores <- abs(agg_result$data$y_mean - agg_result$data$mu_mean)
  n      <- length(scores)
  alpha  <- 1 - level
  
  idx    <- min(n, ceiling((1 - alpha) * (n + 1)))
  margin <- sort(scores)[idx]
  
  empirical_coverage <- mean(scores <= margin)
  
  list(margin = margin, level = level, n_groups = n,
       empirical_coverage = empirical_coverage)
}

# ----------------------------------------------------------------------------
# Auto-detects valid response bounds from the model's family, for clipping
# prediction interval endpoints to a plausible range (e.g. [0,1] for
# bounded-proportion families, [0, Inf) for count/positive-continuous
# families, unbounded for gaussian). Works for both glmmTMB and glmgee
# objects via the generic family() -- glmtoolbox implements a family()
# method for glmgee objects, so this doesn't need model-class-specific code.
#
# Family name matching is based on the string family(model)$family returns --
# worth verifying this matches what you expect for your specific model
# (e.g. run family(m1)$family and check it's in the list below) rather than
# assuming; unrecognized families fall back to no clipping with a warning
# rather than guessing wrong.
# ----------------------------------------------------------------------------
get_response_bounds <- function(model) {
  # glmtoolbox implements no family() S3 method for glmgee objects (confirmed
  # empirically -- it always errors), so branch explicitly rather than paying
  # for a guaranteed-to-fail generic call every time.
  if (inherits(model, "glmgee")) {
    fam_name <- tryCatch(model$family$family, error = function(e) NA_character_)
  } else {
    fam_name <- tryCatch(family(model)$family, error = function(e) NA_character_)
  }
  
  bounded_01     <- c("binomial", "quasibinomial", "beta_family", "ordbeta", "betabinomial", "Beta")
  bounded_below0 <- c("poisson", "quasipoisson", "Gamma", "nbinom1", "nbinom2", "genpois", "compois",
                      "inverse.gaussian", "tweedie",
                      "truncated_poisson", "truncated_nbinom1", "truncated_nbinom2")
  unbounded      <- c("gaussian")
  
  if (is.na(fam_name)) {
    warning("Could not determine model family for automatic response-bound clipping; no clipping applied.")
    return(c(-Inf, Inf))
  }
  if (fam_name %in% bounded_01)     return(c(0, 1))
  if (fam_name %in% bounded_below0) return(c(0, Inf))
  if (fam_name %in% unbounded)      return(c(-Inf, Inf))
  
  warning("Unrecognized family '", fam_name, "' for automatic response-bound clipping; no clipping applied. ",
          "Pass `clip = c(lower, upper)` explicitly if this family has known response bounds.")
  c(-Inf, Inf)
}

# ----------------------------------------------------------------------------
# Apply a margin to marginal predictions from your FINAL model (fit on all
# data) for new data.
#
# margin: from cv_conformal_margin() for individual-row predictions, or
# cv_aggregated_conformal_margin() for AGGREGATE predictions (see agg_cols
# below) -- using the wrong one of these overstates or misrepresents the
# interval, since they're calibrated at different levels.
#
# agg_cols: if supplied, predictions are aggregated (averaged) within each
# combination of these columns BEFORE the margin is applied -- e.g. for a
# counterfactual/focal grid built via build_counterfactual_grid(), pass the
# same by-columns you'd use with avg_predictions(..., by = ...). Requires
# margin to come from cv_aggregated_conformal_margin(), not
# cv_conformal_margin(), to be valid at this level.
#
# Endpoints are computed on the model's native scale first, and
# backtransform_fn -- if supplied -- is applied to fit/lower/upper
# SEPARATELY, never to (mu +/- margin) after transforming components
# individually. This is what makes the transform handling exact for any
# monotonic backtransform, not just a linear one.
#
# clip: NULL (default) auto-detects valid response bounds from the model's
# family via get_response_bounds() -- e.g. [0,1] for ordbeta/binomial,
# [0, Inf) for Gamma/poisson, unbounded for gaussian. Pass c(lower, upper)
# to override, or FALSE to disable clipping entirely.
# ----------------------------------------------------------------------------
predict_marginal_interval <- function(model, newdata, margin,
                                      backtransform_fn = NULL,
                                      clip = NULL,
                                      predict_fn = NULL,
                                      agg_cols = NULL) {
  pred_fn <- if (is.null(predict_fn)) .default_predict_fn else predict_fn
  mu <- pred_fn(model, newdata)
  
  grp_labels <- NULL
  if (!is.null(agg_cols)) {
    df <- newdata[, agg_cols, drop = FALSE]
    df$mu <- mu
    grp_key <- interaction(df[agg_cols], drop = TRUE)
    agg <- do.call(rbind, lapply(split(df, grp_key), function(sub) {
      data.frame(sub[1, agg_cols, drop = FALSE], mu = mean(sub$mu), row.names = NULL)
    }))
    mu <- agg$mu
    grp_labels <- agg[, agg_cols, drop = FALSE]
  }
  
  lower <- mu - margin
  upper <- mu + margin
  
  if (is.null(clip)) clip <- get_response_bounds(model)
  if (!isFALSE(clip)) {
    lower <- pmax(lower, clip[1])
    upper <- pmin(upper, clip[2])
  }
  
  if (!is.null(backtransform_fn)) {
    fit_out   <- backtransform_fn(mu)
    lower_out <- backtransform_fn(lower)
    upper_out <- backtransform_fn(upper)
  } else {
    fit_out <- mu; lower_out <- lower; upper_out <- upper
  }
  
  out <- data.frame(fit = fit_out, lower = lower_out, upper = upper_out)
  if (!is.null(grp_labels)) out <- cbind(grp_labels, out)
  out
}

# ----------------------------------------------------------------------------
# Monte Carlo marginal (population-averaged) prediction, for use as a
# predict_fn override. Integrates out the random effect(s) by simulation
# rather than the re.form = NA plug-in (b = 0), which avoids the
# link^-1(mean(eta)) != mean(link^-1(eta)) bias documented at
# https://github.com/vincentarelbundock/marginaleffects/issues/1231 -- can be
# severe for GLMMs with substantial between-group variance and a nonlinear
# link (logit, log, ...).
#
# Handles arbitrary random-effect structures: random slopes, correlated
# random effects, and multiple grouping factors -- draws each grouping
# factor's random effect vector from a full multivariate normal using its
# ACTUAL estimated covariance matrix (not just the intercept variance), and
# multiplies by the correct per-row random-effect design vector (e.g. for
# (1 + AmbientEnv + I(log10(UASEvents)) | ID), each draw is a 3-vector --
# intercept, AmbientEnv slope, log10(UASEvents) slope -- correlated as
# estimated, applied to each row's own covariate values for those terms).
#
# Requires MASS (ships with base R) and lme4 (for findbars(), to parse the RE
# structure the same way glmmTMB does internally).
#
# Term-name matching between the RE covariance matrix and the constructed
# design matrix is checked and errors rather than silently misaligning -- but
# given the stakes, worth spot-checking directly before trusting output on a
# new model:
#   vc <- glmmTMB::VarCorr(fit)$cond$ID   # or whatever your grouping factor is called
#   dimnames(vc)                          # should read e.g. "(Intercept)", "AmbientEnv", ...
#
# Costs ~n_draws x more compute than the plug-in default -- only worth paying
# for if the two-line mean-vs-mean check shows the plug-in bias matters for
# your data (see usage notes below).
# ----------------------------------------------------------------------------
mc_marginal_predict_fn <- function(fit, newdata, n_draws = 1000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' is required for mc_marginal_predict_fn.")
  if (!requireNamespace("lme4", quietly = TRUE)) stop("Package 'lme4' is required for mc_marginal_predict_fn (uses lme4::findbars).")
  
  eta_fixed <- stats::predict(fit, newdata = newdata, type = "link",
                              re.form = NA, allow.new.levels = TRUE)
  n <- length(eta_fixed)
  eta_draws <- matrix(eta_fixed, nrow = n, ncol = n_draws)
  
  bars <- lme4::findbars(stats::formula(fit))
  vc   <- glmmTMB::VarCorr(fit)$cond
  
  for (bar in bars) {
    grp_name <- deparse(bar[[3]])
    re_rhs   <- deparse(bar[[2]])
    re_form  <- stats::as.formula(paste("~", re_rhs))
    Z        <- stats::model.matrix(re_form, data = newdata)
    
    Sigma <- vc[[grp_name]]
    if (is.null(Sigma)) stop("Could not find grouping factor '", grp_name, "' in VarCorr(fit)$cond.")
    Sigma <- matrix(as.numeric(Sigma), nrow = nrow(Sigma), dimnames = dimnames(Sigma))
    
    ord <- match(colnames(Sigma), colnames(Z))
    if (anyNA(ord)) {
      stop("Could not match random-effect terms to design matrix columns for group '", grp_name,
           "'. VarCorr names: ", paste(colnames(Sigma), collapse = ", "),
           " | design matrix names: ", paste(colnames(Z), collapse = ", "))
    }
    Z <- Z[, ord, drop = FALSE]
    
    b_draws <- MASS::mvrnorm(n_draws, mu = rep(0, ncol(Sigma)), Sigma = Sigma)
    if (ncol(Sigma) == 1) b_draws <- matrix(b_draws, ncol = 1)  # mvrnorm drops to vector for p=1
    eta_draws <- eta_draws + Z %*% t(b_draws)
  }
  
  linkinv <- family(fit)$linkinv
  if (is.null(linkinv)) {
    linkinv <- switch(family(fit)$link,
                      logit = stats::plogis, log = exp, identity = identity,
                      stop("Unsupported link for mc_marginal_predict_fn; supply linkinv manually."))
  }
  
  rowMeans(linkinv(eta_draws))
}

# ============================================================================
# Fast, no-refit alternative to cluster_boot_avg_predictions(): simulates
# BOTH fixed-effect parameter uncertainty (beta* ~ N(beta_hat, vcov(fit)),
# same idea as marginaleffects::inferences(method="simulation") / Krinsky &
# Robb 1986) AND random-effect population uncertainty, combined, for each of
# n_draws replicates -- WITHOUT ever refitting the model. Thousands of draws
# in seconds rather than hundreds of draws in tens of minutes.
#
# UNLIKE mc_marginal_predict_fn, random-effect draws here are INDEPENDENT PER
# ACTUAL SUBJECT (via the grouping column, e.g. ID, which must be present in
# newdata), not one shared draw applied to every row. mc_marginal_predict_fn
# only ever reports row-level point estimates, for which a shared draw is
# harmless (each row's Monte Carlo average is still unbiased regardless of
# what other rows' calculations reuse). This function is different: it
# AGGREGATES across many subjects per stimulus, and if every subject shared
# one draw, the aggregate's simulated variance would never shrink with more
# contributing subjects -- it would just reflect raw between-subject
# variance instead of the variance of a mean over independent subjects. This
# was an actual bug in an earlier version of this function; a symptom to
# watch for if modifying this again: se_sim coming out roughly CONSTANT
# across groups regardless of how many subjects contribute to each is a sign
# independence has been lost somewhere.
#
# The interval is ANCHORED on the analytical point estimate you supply (e.g.
# from marginaleffects::avg_predictions() on the final model), not on the
# mean of the simulated draws -- this deliberately avoids the bootstrap-bias
# problem discussed in chat: for a nonlinear estimator (refitting a GLMM),
# E[bootstrap replicates] != the analytical estimate, and that gap does not
# reliably shrink with more replicates. Anchoring sidesteps it entirely: the
# draws are used ONLY to estimate spread (method = "normal") or a
# bias-reflected interval (method = "basic", Davison & Hinkley 1997), never
# as the center.
#
# Limitation: Sigma_hat (the random-effect covariance) is treated as
# known/fixed at its point estimate -- uncertainty in ESTIMATING the variance
# components themselves isn't propagated. Standard simplifying assumption for
# this kind of plug-in simulation approach, not unique to this implementation.
# ============================================================================

simulate_marginal_ci <- function(fit, newdata, analytical_estimate,
                                 agg_cols = NULL, n_draws = 5000, level = 0.95,
                                 method = c("normal", "basic"),
                                 backtransform_fn = NULL, seed = NULL,
                                 vcov_beta_override = NULL, clip = NULL) {
  method <- match.arg(method)
  if (!is.null(seed)) set.seed(seed)
  if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' is required.")
  if (!requireNamespace("lme4", quietly = TRUE)) stop("Package 'lme4' is required (uses lme4::findbars/nobars).")
  
  # ---- fixed-effect parameter draws (cheap: NOT full=TRUE) -------------------
  beta_hat  <- glmmTMB::fixef(fit)$cond
  
  if (!is.null(vcov_beta_override)) {
    # e.g. cov(do.call(rbind, attr(boot_result, "fixef_boot"))) -- the EMPIRICAL
    # fixed-effect covariance from cluster_boot_avg_predictions(..., return_fixef
    # = TRUE) refits, used in place of the analytical vcov(fit)$cond. Worth
    # using when the model has a high parameter-to-cluster ratio (many fixed
    # effects / a rich random-effect covariance relative to few independent
    # subjects), where Wald/delta-method asymptotics can be unreliable -- see
    # compare_fixef_to_bootstrap() for the diagnostic that would tell you this
    # is needed for a given model.
    vcov_beta <- as.matrix(vcov_beta_override)
    ord_v <- match(names(beta_hat), colnames(vcov_beta))
    if (anyNA(ord_v)) stop("vcov_beta_override column names don't match fixef(fit)$cond names.")
    vcov_beta <- vcov_beta[ord_v, ord_v, drop = FALSE]
  } else {
    vcov_beta <- tryCatch(stats::vcov(fit)$cond, error = function(e) NULL)
    if (is.null(vcov_beta)) vcov_beta <- stats::vcov(fit)  # some versions return the matrix directly
    vcov_beta <- as.matrix(vcov_beta)
  }
  
  X <- stats::model.matrix(lme4::nobars(stats::formula(fit)), data = newdata)
  ord <- match(names(beta_hat), colnames(X))
  if (anyNA(ord)) stop("Could not match fixed-effect names to design matrix columns.")
  X <- X[, ord, drop = FALSE]
  
  beta_draws <- MASS::mvrnorm(n_draws, mu = beta_hat, Sigma = vcov_beta)
  eta_draws  <- X %*% t(beta_draws)   # n x n_draws, fixed-effect uncertainty only so far
  
  # ---- random-effect population draws: INDEPENDENT per actual subject --------
  # (not one shared draw applied to every row -- that's fine for row-level point
  # estimates, as in mc_marginal_predict_fn, but WRONG here: this function
  # aggregates across many subjects per stimulus, and if every subject in a
  # draw shares the same simulated deviation, the aggregate's variance never
  # shrinks with more contributing subjects the way a real sample mean's
  # variance should -- it just reflects the raw between-subject variance.
  # Needs the actual grouping-factor column (e.g. ID) present in newdata.
  bars <- lme4::findbars(stats::formula(fit))
  vc   <- glmmTMB::VarCorr(fit)$cond
  
  for (bar in bars) {
    grp_name <- deparse(bar[[3]])
    re_form  <- stats::as.formula(paste("~", deparse(bar[[2]])))
    Z        <- stats::model.matrix(re_form, data = newdata)
    
    Sigma <- vc[[grp_name]]
    if (is.null(Sigma)) stop("Could not find grouping factor '", grp_name, "' in VarCorr(fit)$cond.")
    Sigma <- matrix(as.numeric(Sigma), nrow = nrow(Sigma), dimnames = dimnames(Sigma))
    
    ordZ <- match(colnames(Sigma), colnames(Z))
    if (anyNA(ordZ)) stop("Could not match random-effect terms to design matrix columns for group '", grp_name, "'.")
    Z <- Z[, ordZ, drop = FALSE]
    
    if (!(grp_name %in% names(newdata))) {
      stop("Grouping column '", grp_name, "' not found in newdata -- required to draw independent ",
           "random effects per subject (not just per row's covariate profile).")
    }
    grp_ids        <- newdata[[grp_name]]
    unique_grp_ids <- unique(grp_ids)
    n_unique       <- length(unique_grp_ids)
    subj_idx       <- match(grp_ids, unique_grp_ids)   # row -> which unique subject
    
    # n_draws * n_unique independent draws in one call, reshaped per RE dimension
    b_all <- MASS::mvrnorm(n_draws * n_unique, mu = rep(0, ncol(Sigma)), Sigma = Sigma)
    if (ncol(Sigma) == 1) b_all <- matrix(b_all, ncol = 1)
    
    for (k in seq_len(ncol(Sigma))) {
      bk_mat <- matrix(b_all[, k], nrow = n_unique, ncol = n_draws)  # n_unique x n_draws
      eta_draws <- eta_draws + Z[, k] * bk_mat[subj_idx, , drop = FALSE]
    }
  }
  
  linkinv <- family(fit)$linkinv
  if (is.null(linkinv)) {
    linkinv <- switch(family(fit)$link,
                      logit = stats::plogis, log = exp, identity = identity,
                      stop("Unsupported link for simulate_marginal_ci; supply linkinv manually."))
  }
  mu_draws <- linkinv(eta_draws)   # n x n_draws
  
  # ---- aggregate per draw if requested (e.g. per-stimulus mean) -------------
  if (!is.null(agg_cols)) {
    grp_key <- interaction(newdata[agg_cols], drop = TRUE)
    draws <- apply(mu_draws, 2, function(col) tapply(col, grp_key, mean))  # n_groups x n_draws
    grp_labels <- newdata[!duplicated(grp_key), agg_cols, drop = FALSE]
    grp_labels <- grp_labels[order(unique(grp_key)), , drop = FALSE]
  } else {
    draws <- mu_draws
    grp_labels <- NULL
  }
  
  if (length(analytical_estimate) != nrow(draws)) {
    stop("length(analytical_estimate) (", length(analytical_estimate), ") does not match the number of ",
         if (is.null(agg_cols)) "rows in newdata" else "groups implied by agg_cols",
         " (", nrow(draws), "). Make sure analytical_estimate is in the same order.")
  }
  
  se_sim <- apply(draws, 1, stats::sd)
  z <- stats::qnorm(1 - (1 - level) / 2)
  
  if (method == "normal") {
    lower <- analytical_estimate - z * se_sim
    upper <- analytical_estimate + z * se_sim
  } else {
    q_lo <- apply(draws, 1, stats::quantile, probs = (1 - level) / 2)
    q_hi <- apply(draws, 1, stats::quantile, probs = 1 - (1 - level) / 2)
    lower <- 2 * analytical_estimate - q_hi
    upper <- 2 * analytical_estimate - q_lo
  }
  
  if (is.null(clip)) clip <- get_response_bounds(fit)
  if (!isFALSE(clip)) {
    lower <- pmax(lower, clip[1])
    upper <- pmin(upper, clip[2])
  }
  
  fit_out   <- analytical_estimate
  lower_out <- lower
  upper_out <- upper
  if (!is.null(backtransform_fn)) {
    fit_out   <- backtransform_fn(fit_out)
    lower_out <- backtransform_fn(lower)
    upper_out <- backtransform_fn(upper)
  }
  
  out <- data.frame(fit = fit_out, lower = lower_out, upper = upper_out, se_sim = se_sim)
  if (!is.null(grp_labels)) out <- cbind(grp_labels, out)
  out
}

# group=...) proves unreliable (it's marked EXPERIMENTAL upstream).
#
# Resamples whole subjects (id_col) with replacement from `full_data` (the
# SAME augmented data frame you'd pass as `newdata` to avg_predictions() --
# i.e. your model's training data plus any extra `by`-grouping columns like
# StimID), refits the model (carrying dispformula/ziformula forward, as
# elsewhere in this file), and predicts on that SAME resampled data. Predicting
# on the original unrelabeled data here would fail: avg_predictions()'s default
# re.form = NULL needs each row's ID to exist as a fitted level in the refit
# model, and relabeled IDs (below) never match the originals.
#
# Duplicated subjects are relabeled with a unique suffix per replicate so
# glmmTMB treats them as distinct random-effect levels -- without this,
# resampling the same subject twice just reweights an existing level rather
# than simulating a new draw from the population, understating uncertainty
# the same way the row-level bootstrap did.
# ============================================================================
cluster_boot_avg_predictions <- function(model, full_data, by, id_col,
                                         R = 100, level = 0.95, seed = NULL,
                                         control = NULL, return_varcorr = FALSE,
                                         return_raw = FALSE, return_fixef = FALSE,
                                         workers = NULL) {
  
  ids       <- unique(full_data[[id_col]])
  form      <- stats::formula(model)
  disp_form <- tryCatch(stats::formula(model, component = "disp"), error = function(e) ~1)
  zi_form   <- tryCatch(stats::formula(model, component = "zi"),   error = function(e) ~0)
  fam       <- family(model)
  
  if (is.null(workers)) workers <- max(1, parallelly::availableCores() - 1)
  old_plan <- future::plan(future::multisession, workers = workers)
  on.exit(future::plan(old_plan), add = TRUE)
  
  boot_out <- future.apply::future_lapply(seq_len(R), function(r) {
    library(glmmTMB)
    
    boot_ids  <- sample(ids, length(ids), replace = TRUE)
    boot_rows <- lapply(seq_along(boot_ids), function(i) {
      d <- full_data[full_data[[id_col]] == boot_ids[i], ]
      d[[id_col]] <- paste0(d[[id_col]], "_rep", i)  # distinct RE level per draw
      d
    })
    boot_data <- do.call(rbind, boot_rows)
    
    fit_args <- list(formula = form, data = boot_data, family = fam,
                     dispformula = disp_form, ziformula = zi_form, se = FALSE)
    if (!is.null(control)) fit_args$control <- control
    
    fit_err <- NULL
    fit_r <- tryCatch(do.call(glmmTMB::glmmTMB, fit_args),
                      error = function(e) { fit_err <<- conditionMessage(e); NULL })
    if (is.null(fit_r)) return(list(.failed = TRUE, .error = paste("fit:", fit_err)))
    
    # Predict on the SAME resampled+relabeled data the model was just fit on.
    # avg_predictions()'s default re.form = NULL needs each row's ID to exist
    # as a fitted level in fit_r -- predicting on the original, unrelabeled
    # newdata here breaks that match for every row, every replicate.
    pred_err <- NULL
    pred_r <- tryCatch(
      marginaleffects::avg_predictions(fit_r, newdata = boot_data, by = by,
                                       vcov = FALSE, type = "response"),
      error = function(e) { pred_err <<- conditionMessage(e); NULL }
    )
    if (is.null(pred_r)) return(list(.failed = TRUE, .error = paste("predict:", pred_err)))
    
    out <- list(
      .failed = FALSE,
      # keep the by-column identifiers alongside estimate -- do NOT rely on
      # row position/order being stable across replicates.
      pred = as.data.frame(pred_r)[, c(by, "estimate"), drop = FALSE]
    )
    if (return_varcorr) {
      out$vc <- tryCatch(glmmTMB::VarCorr(fit_r)$cond, error = function(e) NULL)
    }
    if (return_fixef) {
      out$fx <- tryCatch(glmmTMB::fixef(fit_r)$cond, error = function(e) NULL)
    }
    out
  }, future.seed = if (is.null(seed)) TRUE else seed)
  
  is_failure   <- vapply(boot_out, function(f) isTRUE(f$.failed), logical(1))
  boot_errors  <- lapply(boot_out[is_failure], function(f) f$.error)
  boot_out     <- boot_out[!is_failure]
  
  if (length(boot_out) == 0) {
    stop("All ", R, " bootstrap replicates failed.",
         if (length(boot_errors) > 0) paste0(" First error encountered: ", boot_errors[[1]]) else "")
  }
  
  n_fail <- sum(is_failure)
  if (n_fail > 0) {
    message(n_fail, " of ", R, " bootstrap replicates failed to fit/predict.")
    # Distinct messages with counts, not just the first -- if failures share a
    # common cause (e.g. "boundary (singular) fit"), that's directly relevant
    # to whether they're excluding a particular tail of the resampling
    # distribution rather than failing at random.
    err_table <- table(unlist(boot_errors))
    err_table <- sort(err_table, decreasing = TRUE)
    message("Failure messages (count x message):")
    for (i in seq_along(err_table)) {
      message("  ", err_table[i], " x ", names(err_table)[i])
    }
  }
  
  boot_ests <- lapply(boot_out, `[[`, "pred")
  boot_vc   <- if (return_varcorr) lapply(boot_out, `[[`, "vc") else NULL
  boot_fx   <- if (return_fixef)   lapply(boot_out, `[[`, "fx") else NULL
  
  long_df <- do.call(rbind, boot_ests)
  grp_key <- interaction(long_df[by], drop = TRUE)
  alpha   <- 1 - level
  
  result <- do.call(rbind, lapply(split(long_df, grp_key), function(sub) {
    data.frame(
      sub[1, by, drop = FALSE],
      n_boot   = nrow(sub),   # how many replicates actually contributed to this row's CI
      estimate = mean(sub$estimate, na.rm = TRUE),
      sd_boot  = stats::sd(sub$estimate, na.rm = TRUE),  # actual SD of raw replicates --
      # compare directly against
      # simulate_marginal_ci()'s se_sim,
      # no percentile-range approximation
      lower    = stats::quantile(sub$estimate, probs = alpha / 2, na.rm = TRUE),
      upper    = stats::quantile(sub$estimate, probs = 1 - alpha / 2, na.rm = TRUE),
      row.names = NULL
    )
  }))
  
  rownames(result) <- NULL
  if (return_varcorr) attr(result, "varcorr_boot") <- Filter(Negate(is.null), boot_vc)
  if (return_fixef)   attr(result, "fixef_boot")    <- Filter(Negate(is.null), boot_fx)
  if (return_raw)      attr(result, "raw_boot")     <- long_df
  result
}

# ----------------------------------------------------------------------------
# Compares the ORIGINAL model's variance-component point estimate against the
# distribution of re-estimated values across bootstrap replicates (from
# cluster_boot_avg_predictions(..., return_varcorr = TRUE)). If the original
# estimate sits well outside the typical bootstrap-refit range, that's direct
# evidence it's an unstable/inflated estimate -- and a likely explanation for
# simulate_marginal_ci() (which treats it as fixed/known) coming out
# systematically wider than the bootstrap.
# ----------------------------------------------------------------------------
compare_varcorr_to_bootstrap <- function(model, boot_result, grp_name) {
  vc_boot <- attr(boot_result, "varcorr_boot")
  if (is.null(vc_boot)) stop("boot_result must come from cluster_boot_avg_predictions(..., return_varcorr = TRUE).")
  
  vc_orig <- matrix(as.numeric(glmmTMB::VarCorr(model)$cond[[grp_name]]),
                    nrow = nrow(glmmTMB::VarCorr(model)$cond[[grp_name]]),
                    dimnames = dimnames(glmmTMB::VarCorr(model)$cond[[grp_name]]))
  p <- nrow(vc_orig)
  term_names <- colnames(vc_orig)
  
  # variances (diagonal) and correlations, since raw covariances aren't
  # directly comparable in magnitude across terms on different scales
  orig_var  <- diag(vc_orig)
  orig_cor  <- stats::cov2cor(vc_orig)
  
  boot_var <- t(vapply(vc_boot, function(v) {
    m <- v[[grp_name]]
    if (is.null(m)) return(rep(NA_real_, p))
    diag(matrix(as.numeric(m), nrow = nrow(m)))
  }, numeric(p)))
  colnames(boot_var) <- term_names
  
  cat("Variance components for grouping factor '", grp_name, "':\n\n", sep = "")
  for (j in seq_len(p)) {
    b <- boot_var[, j]
    b <- b[!is.na(b)]
    cat(sprintf("  %-30s original = %.4f | bootstrap median = %.4f, range [%.4f, %.4f], n = %d\n",
                term_names[j], orig_var[j], stats::median(b), min(b), max(b), length(b)))
  }
  
  # correlations between RE terms -- matching marginal variances doesn't rule
  # out a mismatch here, and correlated slopes can matter for the aggregate
  # even when each term's own variance looks fine individually.
  if (p > 1) {
    boot_cor <- lapply(vc_boot, function(v) {
      m <- v[[grp_name]]
      if (is.null(m)) return(NULL)
      stats::cov2cor(matrix(as.numeric(m), nrow = nrow(m)))
    })
    boot_cor <- Filter(Negate(is.null), boot_cor)
    
    cat("\nCorrelations between random-effect terms:\n\n")
    for (i in 1:(p - 1)) {
      for (j in (i + 1):p) {
        b <- vapply(boot_cor, function(m) m[i, j], numeric(1))
        cat(sprintf("  %s x %s\n    original = %.3f | bootstrap median = %.3f, range [%.3f, %.3f]\n",
                    term_names[i], term_names[j], orig_cor[i, j],
                    stats::median(b), min(b), max(b)))
      }
    }
  }
  
  cat("\nIf 'original' sits outside or near the edge of the bootstrap range for any term,\n",
      "that term's variance/correlation estimate is likely unstable -- worth treating\n",
      "simulate_marginal_ci's output for that source with real caution.\n", sep = "")
  
  invisible(list(orig_var = orig_var, orig_cor = orig_cor, boot_var = boot_var))
}

# ----------------------------------------------------------------------------
# Compares vcov(model)$cond (the fixed-effect covariance simulate_marginal_ci
# draws beta* from) against the EMPIRICAL variance of fixed effects across
# bootstrap refits (from cluster_boot_avg_predictions(..., return_fixef = TRUE)).
# Random-effect variance components and correlations have already been
# checked via compare_varcorr_to_bootstrap() and found consistent -- this
# checks the other half of what simulate_marginal_ci() simulates. If
# vcov(model)$cond's diagonal is systematically larger than the bootstrap's
# empirical per-coefficient variance, that directly explains
# simulate_marginal_ci() coming out wider: it would be drawing fixed-effect
# uncertainty from a wider distribution than the data (via refitting) actually
# supports.
# ----------------------------------------------------------------------------
compare_fixef_to_bootstrap <- function(model, boot_result) {
  fx_boot <- attr(boot_result, "fixef_boot")
  if (is.null(fx_boot)) stop("boot_result must come from cluster_boot_avg_predictions(..., return_fixef = TRUE).")
  
  beta_hat  <- glmmTMB::fixef(model)$cond
  vcov_orig <- tryCatch(stats::vcov(model)$cond, error = function(e) NULL)
  if (is.null(vcov_orig)) vcov_orig <- stats::vcov(model)
  var_orig  <- diag(as.matrix(vcov_orig))
  
  fx_mat <- do.call(rbind, fx_boot)   # n_replicates x n_coef
  # align columns by name in case ordering ever differs across refits
  fx_mat <- fx_mat[, names(beta_hat), drop = FALSE]
  var_boot <- apply(fx_mat, 2, stats::var)
  
  cat("Fixed-effect coefficient variance: analytical (vcov) vs empirical (bootstrap refits)\n\n")
  for (j in seq_along(beta_hat)) {
    ratio <- var_orig[j] / var_boot[j]
    cat(sprintf("  %-30s analytical = %.5f | bootstrap empirical = %.5f | ratio = %.2f\n",
                names(beta_hat)[j], var_orig[j], var_boot[j], ratio))
  }
  cat("\nratio > 1 means vcov(model) implies MORE fixed-effect uncertainty than the bootstrap's\n",
      "own refits actually show for that coefficient -- consistently high ratios across most\n",
      "coefficients would directly explain simulate_marginal_ci() running systematically wider.\n", sep = "")
  
  invisible(list(var_orig = var_orig, var_boot = var_boot, n_boot = nrow(fx_mat)))
}

# ============================================================================
# Cluster bootstrap for COUNTERFACTUAL predictions (marginaleffects-style
# variables = list(...) sweeps/crosses), e.g. predicted response across a
# range of UASLAEMaxLRScl, or crossed with UASType.
#
# cf_variables: which column(s) get counterfactually OVERRIDDEN (swept) --
# same format as build_counterfactual_grid()'s `variables`, and as
# marginaleffects::avg_predictions(variables = ...).
#
# by: the full output GROUPING, same role as avg_predictions(by = ...).
# Defaults to names(cf_variables), but can include additional OBSERVED
# (non-swept) covariates you want broken out separately in the output rather
# than averaged over -- e.g. by = c("UASLAEMaxLRScl", "UASEvents",
# "UASOperation", "AmbientEnv") while only UASLAEMaxLRScl is swept. Getting
# this distinction right matters: leaving `by` at its default when you
# actually wanted extra grouping columns will silently average over them
# instead of keeping them separate.
#
# full_data can be the model's full training data (e.g. insight::get_data(
# model)) or any SUBSET of it -- nothing here assumes the full dataset;
# whatever rows you pass are what gets resampled from and counterfactually
# expanded, matching what you'd get from passing that same subset as
# `newdata` to avg_predictions() directly.
#
# Deliberately a SEPARATE function from cluster_boot_avg_predictions() rather
# than a modification to it, given how much that function has already been
# validated -- this reuses the same refit-per-replicate structure but builds
# the counterfactual grid from each replicate's OWN resampled+relabeled data
# (boot_data, not the original full_data), so the relabeled IDs the
# counterfactual rows inherit still match levels fit_r actually has.
#
# The refit cost is IDENTICAL to cluster_boot_avg_predictions() for the same
# R -- only the (cheap) prediction step scales with the size of the
# counterfactual grid, not the expensive refit step.
# ============================================================================

cluster_boot_counterfactual <- function(model, full_data, id_col, cf_variables,
                                        by = NULL,
                                        R = 100, level = 0.95, seed = NULL,
                                        control = NULL, workers = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(by)) by <- names(cf_variables)   # default: group by the swept variable(s) only
  
  ids       <- unique(full_data[[id_col]])
  form      <- stats::formula(model)
  disp_form <- tryCatch(stats::formula(model, component = "disp"), error = function(e) ~1)
  zi_form   <- tryCatch(stats::formula(model, component = "zi"),   error = function(e) ~0)
  fam       <- family(model)
  
  if (is.null(workers)) workers <- max(1, parallelly::availableCores() - 1)
  old_plan <- future::plan(future::multisession, workers = workers)
  on.exit(future::plan(old_plan), add = TRUE)
  
  boot_out <- future.apply::future_lapply(seq_len(R), function(r) {
    library(glmmTMB)
    
    boot_ids  <- sample(ids, length(ids), replace = TRUE)
    boot_rows <- lapply(seq_along(boot_ids), function(i) {
      d <- full_data[full_data[[id_col]] == boot_ids[i], ]
      d[[id_col]] <- paste0(d[[id_col]], "_rep", i)
      d
    })
    boot_data <- do.call(rbind, boot_rows)
    
    fit_args <- list(formula = form, data = boot_data, family = fam,
                     dispformula = disp_form, ziformula = zi_form, se = FALSE)
    if (!is.null(control)) fit_args$control <- control
    
    fit_err <- NULL
    fit_r <- tryCatch(do.call(glmmTMB::glmmTMB, fit_args),
                      error = function(e) { fit_err <<- conditionMessage(e); NULL })
    if (is.null(fit_r)) return(list(.failed = TRUE, .error = paste("fit:", fit_err)))
    
    # counterfactual grid built from THIS replicate's resampled+relabeled
    # data, so relabeled IDs match fit_r's fitted levels (same reasoning as
    # cluster_boot_avg_predictions -- predicting on the original unrelabeled
    # data here would break re.form = NULL's default level-matching).
    boot_cf <- build_counterfactual_grid(boot_data, cf_variables)
    
    pred_err <- NULL
    pred_r <- tryCatch(
      marginaleffects::avg_predictions(fit_r, newdata = boot_cf, by = by,
                                       vcov = FALSE, type = "response"),
      error = function(e) { pred_err <<- conditionMessage(e); NULL }
    )
    if (is.null(pred_r)) return(list(.failed = TRUE, .error = paste("predict:", pred_err)))
    
    list(.failed = FALSE, pred = as.data.frame(pred_r)[, c(by, "estimate"), drop = FALSE])
  }, future.seed = if (is.null(seed)) TRUE else seed)
  
  is_failure  <- vapply(boot_out, function(f) isTRUE(f$.failed), logical(1))
  boot_errors <- lapply(boot_out[is_failure], function(f) f$.error)
  boot_out    <- boot_out[!is_failure]
  
  if (length(boot_out) == 0) {
    stop("All ", R, " bootstrap replicates failed.",
         if (length(boot_errors) > 0) paste0(" First error encountered: ", boot_errors[[1]]) else "")
  }
  
  n_fail <- sum(is_failure)
  if (n_fail > 0) {
    message(n_fail, " of ", R, " bootstrap replicates failed to fit/predict.")
    err_table <- sort(table(unlist(boot_errors)), decreasing = TRUE)
    message("Failure messages (count x message):")
    for (i in seq_along(err_table)) message("  ", err_table[i], " x ", names(err_table)[i])
  }
  
  boot_ests <- lapply(boot_out, `[[`, "pred")
  long_df   <- do.call(rbind, boot_ests)
  grp_key   <- interaction(long_df[by], drop = TRUE)
  alpha     <- 1 - level
  
  result <- do.call(rbind, lapply(split(long_df, grp_key), function(sub) {
    data.frame(
      sub[1, by, drop = FALSE],
      n_boot   = nrow(sub),
      estimate = mean(sub$estimate, na.rm = TRUE),
      sd_boot  = stats::sd(sub$estimate, na.rm = TRUE),
      lower    = stats::quantile(sub$estimate, probs = alpha / 2, na.rm = TRUE),
      upper    = stats::quantile(sub$estimate, probs = 1 - alpha / 2, na.rm = TRUE),
      row.names = NULL
    )
  }))
  
  rownames(result) <- NULL
  result
}

# ============================================================================
# GEE (glmtoolbox::glmgee) support.
#
# GEE models the marginal (population-averaged) mean directly -- no random
# effect to marginalize, no re.form ambiguity, no Jensen's-gap plug-in bias.
# Trade-off: GEE is an estimating-equations method, not a full-likelihood one,
# so there is no valid elpd for it in the AIC/WAIC sense (glmtoolbox's
# logLik()/QIC() is a quasi-likelihood for model SELECTION, not a predictive
# density) and there's no subject-specific/conditional prediction available,
# ever -- only the marginal mean. See cv_glmgee_pseudo_loglik() below, though,
# for a genuine out-of-sample analogue built from the same Gaussian
# pseudo-likelihood AGPC/SGPC use internally -- evaluated on held-out
# clusters rather than as an in-sample penalized criterion.
#
# id_col is a REQUIRED argument here, not auto-detected: `id` is passed to
# glmgee() as a bare column reference, and I haven't confirmed a stable
# accessor for it on the fitted object across glmtoolbox versions -- safer to
# ask than guess, unlike the insight::find_random() route used for glmmTMB.
# ============================================================================

.glmgee_default_predict_fn <- function(fit, newdata) {
  as.numeric(stats::predict(fit, newdata = newdata, type = "response"))
}

# ----------------------------------------------------------------------------
# Out-of-sample Gaussian pseudo-log-likelihood for GEE, evaluated per CLUSTER
# on held-out data -- the elpd-style predictive score glmgee models otherwise
# lack. Built from the same per-cluster fit term AGPC/SGPC use internally
# (Carey & Wang, 2011; see glmtoolbox's geeglm.R and Vanegas, Rondon & Paula,
# 2023, The R Journal 15:105-133), WITHOUT their in-sample complexity penalty
# -- this scores genuinely held-out clusters instead, the same relationship
# elpd_loo has to WAIC/AIC.
#
# SCOPED to corstr = "exchangeable" or "independence" only. Other working
# correlation structures (AR-m, unstructured, etc.) need their own R_i
# reconstruction and are NOT implemented here -- extending this to those is
# a real undertaking, not a small addition, given how differently each
# structure is parameterized; not attempted given the added complexity/bug
# risk relative to what's actually needed here.
#
# ACCESSOR CAVEAT: phi_hat (dispersion) and rho_hat (working correlation
# parameter) are extracted via best-guess accessors on the glmgee object,
# with fallbacks, but I have NOT been able to verify these against a live
# glmtoolbox session. Run cv_glmgee_pseudo_loglik_sanity_check() (below)
# before trusting this, the same way cv_glmmTMB_sanity_check() was used to
# validate the ordbeta density scale much earlier in this file.
# ----------------------------------------------------------------------------

.glmgee_variance_fun <- function(family_name) {
  fn <- switch(family_name,
               binomial = , quasibinomial = function(mu) mu * (1 - mu),
               poisson  = , quasipoisson  = function(mu) mu,
               Gamma    = function(mu) mu^2,
               gaussian = function(mu) rep(1, length(mu)),
               inverse.gaussian = function(mu) mu^3,
               NULL
  )
  if (is.null(fn)) {
    stop("No built-in variance function for family '", family_name, "'. ",
         "Supply variance_fun directly to gee_pseudo_loglik().")
  }
  fn
}

# Per-cluster pseudo-log-likelihood (NOT the -2x, penalized IC form -- this
# is the log-likelihood itself, summable across clusters like elpd).
gee_pseudo_loglik <- function(y, mu, id, phi, rho = 0,
                              corstr = c("exchangeable", "independence"),
                              variance_fun = NULL, family_name = NULL) {
  corstr <- match.arg(corstr)
  if (is.null(variance_fun)) {
    if (is.null(family_name)) stop("Supply either variance_fun or family_name.")
    variance_fun <- .glmgee_variance_fun(family_name)
  }
  
  # id may arrive as a factor carrying ALL subject levels from the full
  # dataset, not just the ones present in this particular held-out set --
  # split() groups by level, not by value present, so an unconverted factor
  # silently creates one empty group per absent subject. Those score NA
  # (correctly, since they have no data), but at scale this is exactly what
  # produced ~168/210 NA clusters and the 210 = 42 subjects x 5 folds
  # bookkeeping artifact -- not a numerical/singularity issue as first
  # suspected. Confirmed by the 210 figure factoring exactly; the earlier
  # near-singular-Vi hypothesis is wrong for this data (every printed rho
  # sits at 0.20-0.24, nowhere near the critical value).
  clusters <- split(seq_along(y), as.character(id))
  vapply(clusters, function(idx) {
    yi <- y[idx]; mui <- mu[idx]; ni <- length(idx)
    a_sqrt <- sqrt(variance_fun(mui))
    Ri <- if (corstr == "independence" || ni == 1) diag(ni) else (1 - rho) * diag(ni) + rho * matrix(1, ni, ni)
    Vi <- diag(a_sqrt, ni) %*% Ri %*% diag(a_sqrt, ni)
    
    resid <- yi - mui
    quad    <- tryCatch(as.numeric(t(resid) %*% solve(Vi) %*% resid), error = function(e) NA_real_)
    logdetV <- tryCatch(determinant(Vi, logarithm = TRUE)$modulus[1], error = function(e) NA_real_)
    
    -0.5 * (ni * log(2 * pi) + ni * log(phi) + logdetV + quad / phi)
  }, numeric(1))
}

# Extracts phi_hat/rho_hat from a fitted glmgee object. phi confirmed correct
# (fit$dispersion matches summary()'s "Dispersion" exactly). rho: fit$corr
# returns the FULL working correlation MATRIX (confirmed empirically -- same
# values as summary()'s printed "Working correlation"), not a bare scalar --
# so it needs an off-diagonal element, not element [1] (which is the
# diagonal, trivially 1, and was the bug that produced rho = 1 and a
# singular, non-invertible reconstructed Ri).
#
# NOTE on performance: fit$dispersion returns NULL for do.call()-constructed
# glmgee objects specifically (confirmed: worked directly on m1b itself, but
# not on a do.call()-built refit) -- forcing every CV fold through the
# summary(fit)$dispersion fallback. summary.glmgee() appears to print as a
# side effect of being called at all (not just on autoprint), which is what
# produced the large per-fold output blocks seen during cv_glmgee_pseudo_loglik
# runs, and is a real, avoidable cost (large matrix formatting + console I/O)
# on top of whatever summary() computes -- wrapped in capture.output() below
# so folds run quietly. This does NOT address any other computational cost
# inside summary.glmgee() itself, only the printing.
.quiet_summary_glmgee <- function(fit) {
  smry <- NULL
  invisible(utils::capture.output(smry <- tryCatch(summary(fit), error = function(e) NULL)))
  smry
}

.glmgee_extract_phi_rho <- function(fit) {
  # fit$phi confirmed present directly on do.call()-built objects (see
  # names(fit_test) check) -- glmtoolbox's actual internal field, matching
  # standard GEE notation. Tried first now, ahead of $dispersion, which only
  # resolved for some fit paths and forced the expensive summary() fallback
  # for every CV fold otherwise.
  phi <- tryCatch(fit$phi, error = function(e) NULL)
  if (is.null(phi)) phi <- tryCatch(fit$dispersion, error = function(e) NULL)
  if (is.null(phi)) {
    smry <- .quiet_summary_glmgee(fit)
    phi  <- if (!is.null(smry)) tryCatch(smry$dispersion, error = function(e) NULL) else NULL
  }
  if (is.null(phi)) phi <- tryCatch(glmtoolbox::vcov.glmgee(fit)[1, 1], error = function(e) NULL)  # last resort, likely wrong
  if (is.null(phi)) stop("Could not extract dispersion (phi) from glmgee object -- accessor unconfirmed, see caveat.")
  
  rho_raw <- tryCatch(fit$corr, error = function(e) NULL)
  if (is.null(rho_raw)) rho_raw <- tryCatch(fit$rho, error = function(e) NULL)
  if (is.null(rho_raw)) {
    smry <- .quiet_summary_glmgee(fit)
    rho_raw <- if (!is.null(smry)) tryCatch(smry$corr, error = function(e) NULL) else NULL
  }
  
  rho <- if (is.null(rho_raw)) {
    0   # independence corstr has no rho
  } else if (is.matrix(rho_raw) && nrow(rho_raw) >= 2) {
    rho_raw[1, 2]   # off-diagonal -- valid for exchangeable, where every off-diagonal is rho
  } else {
    as.numeric(rho_raw)[1]
  }
  
  list(phi = as.numeric(phi)[1], rho = rho)
}

# One-off check: does gee_pseudo_loglik(), summed over TRAINING clusters,
# land in a sane ballpark relative to what AGPC/SGPC's own fit term implies?
# Run this BEFORE trusting cv_glmgee_pseudo_loglik() on a new model.
cv_glmgee_pseudo_loglik_sanity_check <- function(model, id_col, corstr = c("exchangeable", "independence")) {
  corstr <- match.arg(corstr)
  data <- insight::get_data(model)
  resp <- all.vars(stats::formula(model))[1]
  mu   <- stats::predict(model, type = "response")
  y    <- data[[resp]]
  id   <- data[[id_col]]
  
  pr <- .glmgee_extract_phi_rho(model)
  cat("Extracted phi =", pr$phi, " rho =", pr$rho, "-- VERIFY these against summary(model)\n")
  
  ll <- sum(gee_pseudo_loglik(y, mu, id, phi = pr$phi, rho = pr$rho,
                              corstr = corstr, family_name = model$family$family))
  cat("sum(gee_pseudo_loglik(...)) on training data:", ll, "\n")
  cat("(compare qualitatively against -0.5*(AGPC(model) - penalty) if you want a direct check;\n",
      " see the AGPC/SGPC formula in the header comment above for the penalty term)\n", sep = "")
  invisible(ll)
}

# ----------------------------------------------------------------------------
# cv_glmgee_pseudo_loglik(): out-of-sample version, K-fold, same grouped
# structure as cv_glmgee(). Standalone rather than folded into cv_glmgee()
# itself, to avoid touching that function's already-validated fold loop.
# ----------------------------------------------------------------------------
cv_glmgee_pseudo_loglik <- function(model, id_col, K = 6, seed = 123,
                                    corstr = c("exchangeable", "independence"),
                                    data = NULL, workers = NULL) {
  corstr <- match.arg(corstr)
  form <- tryCatch(model$formula, error = function(e) stats::formula(model))
  fam  <- tryCatch(model$family,  error = function(e) NULL)
  if (is.null(fam)) stop("Could not extract `family` from `model`.")
  fam_name <- fam$family
  
  # The bug this fixes: fold refits previously never received the model's
  # actual working-correlation structure, so glmgee() silently fell back to
  # its own default (independence) for every fold, regardless of what
  # `model` -- or the `corstr` argument above -- specified. Confirmed
  # directly from the printed fold summaries showing "Correlation structure:
  # Independence" for every fold despite requesting exchangeable.
  model_corstr <- tryCatch(model$corstr, error = function(e) NULL)
  if (is.null(model_corstr)) {
    stop("Could not extract `corstr` from `model` via model$corstr -- accessor unconfirmed. ",
         "Run `m1b$corstr` (or similar) directly and check what it returns before proceeding.")
  }
  if (!grepl(corstr, model_corstr, ignore.case = TRUE)) {
    warning("Requested corstr = '", corstr, "' for scoring does not obviously match model's own ",
            "corstr ('", model_corstr, "'). The fold refits will use model's actual corstr; make sure ",
            "the `corstr` argument to this function still matches, or gee_pseudo_loglik()'s closed-form ",
            "Ri reconstruction won't correspond to what was actually fit.")
  }
  
  if (is.null(data)) data <- tryCatch(insight::get_data(model), error = function(e) NULL)
  if (is.null(data)) stop("insight::get_data() could not extract data; pass `data` directly.")
  resp <- all.vars(form)[1]
  
  clean_data <- stats::na.omit(data[, unique(c(all.vars(form), id_col))])
  
  set.seed(seed)
  unique_ids <- unique(clean_data[[id_col]])
  id_folds   <- sample(rep(seq_len(K), length.out = length(unique_ids)))
  names(id_folds) <- unique_ids
  data_folds <- id_folds[as.character(clean_data[[id_col]])]
  
  if (is.null(workers)) workers <- max(1, parallelly::availableCores() - 1)
  old_plan <- future::plan(future::multisession, workers = workers)
  on.exit(future::plan(old_plan), add = TRUE)
  
  fold_ll <- future.apply::future_lapply(seq_len(K), function(k) {
    library(glmtoolbox)
    
    train <- clean_data[data_folds != k, ]
    val   <- clean_data[data_folds == k, ]
    
    fit <- tryCatch(
      do.call(glmtoolbox::glmgee,
              list(formula = form, id = as.name(id_col), family = fam, corstr = model_corstr, data = train)),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    
    mu <- tryCatch(stats::predict(fit, newdata = val, type = "response"), error = function(e) NULL)
    if (is.null(mu)) return(NULL)
    
    pr <- tryCatch(.glmgee_extract_phi_rho(fit), error = function(e) NULL)
    if (is.null(pr)) return(NULL)
    
    tryCatch(
      gee_pseudo_loglik(val[[resp]], mu, val[[id_col]], phi = pr$phi, rho = pr$rho,
                        corstr = corstr, family_name = fam_name),
      error = function(e) NULL
    )
  }, future.seed = TRUE)
  
  fold_ll <- Filter(Negate(is.null), fold_ll)
  if (length(fold_ll) == 0) stop("All folds failed.")
  n_fail <- K - length(fold_ll)
  if (n_fail > 0) message(n_fail, " of ", K, " folds failed.")
  
  all_ll <- unlist(fold_ll, use.names = FALSE)
  n_na <- sum(is.na(all_ll))
  if (n_na > 0) {
    message(n_na, " of ", length(all_ll), " scored clusters produced NA pseudo-loglik -- likely a ",
            "near-singular working covariance (Vi) for that specific fold's estimated rho and a large ",
            "cluster size, not a general failure. Excluded via na.rm when summing; if n_na is more than ",
            "a handful, that's worth investigating further rather than just accepting the total.")
  }
  list(cv_pseudo_loglik = sum(all_ll, na.rm = TRUE),
       n_clusters_scored = sum(!is.na(all_ll)),
       n_clusters_na     = n_na,
       K_used = length(fold_ll))
}

cv_glmgee <- function(model,
                      id_col,
                      K            = 6,
                      metrics      = c("rmse", "mae", "r2"),
                      se           = FALSE,
                      seed         = 123,
                      workers      = NULL,
                      backtransform_fn   = NULL,
                      return_predictions = FALSE,
                      predict_fn   = NULL,
                      agg_cols     = NULL,
                      data         = NULL) {
  
  stopifnot(inherits(model, "glmgee"))
  
  if ("elpd" %in% metrics) {
    warning("elpd is not defined for GEE models (estimating equations, not a full ",
            "joint likelihood -- see header comment above). Dropping 'elpd' from metrics.")
    metrics <- setdiff(metrics, "elpd")
  }
  
  form <- tryCatch(model$formula, error = function(e) stats::formula(model))
  fam  <- tryCatch(model$family,  error = function(e) NULL)
  cs   <- tryCatch(model$corstr,  error = function(e) NULL)
  if (is.null(fam)) stop("Could not extract `family` from `model` ($family); glmgee object structure may have changed.")
  
  # insight::get_data() only reconstructs FORMULA variables from the model's
  # internal frame -- columns bound on for convenience (e.g. StimID, not used
  # in the model) won't survive that route. Supply `data` directly (e.g. the
  # same cbind(m1Data, StimID = ...) frame used elsewhere) if you need
  # agg_cols that aren't actual predictors -- same reasoning as cv_glmmTMB.
  if (is.null(data)) data <- tryCatch(insight::get_data(model), error = function(e) NULL)
  if (is.null(data)) {
    stop("insight::get_data() could not extract data from this glmgee object; pass `data` directly.")
  }
  resp <- all.vars(form)[1]
  
  if (!is.null(agg_cols) && !all(agg_cols %in% names(data))) {
    stop("agg_cols ", paste(setdiff(agg_cols, names(data)), collapse = ", "),
         " not found in `data`. If these aren't actual predictors in the model (e.g. a StimID ",
         "label), pass your own data frame via cv_glmgee(..., data = your_full_data_with_StimID).")
  }
  
  keep_cols  <- unique(c(all.vars(form), id_col, agg_cols))
  clean_data <- stats::na.omit(data[, keep_cols])
  
  set.seed(seed)
  unique_ids <- unique(clean_data[[id_col]])
  id_folds   <- sample(rep(seq_len(K), length.out = length(unique_ids)))
  names(id_folds) <- unique_ids
  data_folds <- id_folds[as.character(clean_data[[id_col]])]
  
  if (is.null(workers)) workers <- max(1, parallelly::availableCores() - 1)
  old_plan <- future::plan(future::multisession, workers = workers)
  on.exit(future::plan(old_plan), add = TRUE)
  
  fold_out <- future.apply::future_lapply(seq_len(K), function(k) {
    library(glmtoolbox)
    
    train <- clean_data[data_folds != k, ]
    val   <- clean_data[data_folds == k, ]
    
    # id must be passed as a SYMBOL, not a pre-extracted vector: glmgee()
    # evaluates `id` via non-standard evaluation against `data` (same as your
    # own working call uses id = ID, a bare column reference, not a string or
    # vector). do.call() with as.name(id_col) reproduces that correctly,
    # generalized to whatever id_col actually is.
    fit_err <- NULL
    fit <- tryCatch(
      do.call(glmtoolbox::glmgee,
              list(formula = form, id = as.name(id_col), family = fam, corstr = cs, data = train)),
      error = function(e) { fit_err <<- conditionMessage(e); NULL }
    )
    if (is.null(fit)) return(list(.failed = TRUE, .error = paste("fit:", fit_err)))
    
    pred_fn  <- if (is.null(predict_fn)) .glmgee_default_predict_fn else predict_fn
    pred_err <- NULL
    mu <- tryCatch(pred_fn(fit, val),
                   error = function(e) { pred_err <<- conditionMessage(e); NULL })
    if (is.null(mu)) return(list(.failed = TRUE, .error = paste("predict:", pred_err)))
    
    out <- list(.failed = FALSE, y = val[[resp]], mu = mu)
    if (!is.null(agg_cols)) out$agg <- val[, agg_cols, drop = FALSE]
    out
  }, future.seed = TRUE)
  
  is_failure  <- vapply(fold_out, function(f) isTRUE(f$.failed), logical(1))
  fold_errors <- lapply(fold_out[is_failure], function(f) f$.error)
  fold_out    <- fold_out[!is_failure]
  n_fail      <- sum(is_failure)
  if (n_fail > 0) {
    message(n_fail, " of ", K, " folds failed to fit/predict and were dropped.")
    err_table <- sort(table(unlist(fold_errors)), decreasing = TRUE)
    message("Failure messages (count x message):")
    for (i in seq_along(err_table)) message("  ", err_table[i], " x ", names(err_table)[i])
  }
  if (length(fold_out) == 0) stop("All folds failed.",
                                  if (length(fold_errors) > 0) paste0(" First error: ", fold_errors[[1]]) else "")
  
  y_all  <- unlist(lapply(fold_out, `[[`, "y"),  use.names = FALSE)
  mu_all <- unlist(lapply(fold_out, `[[`, "mu"), use.names = FALSE)
  resid  <- y_all - mu_all
  n_obs  <- length(y_all)
  
  result <- list(K_requested = K, K_used = length(fold_out), n_obs = n_obs)
  
  if ("rmse" %in% metrics) {
    mse <- mean(resid^2, na.rm = TRUE)
    result$mse  <- mse   # = CV Brier score when y_all is binary (0/1)
    result$rmse <- sqrt(mse)
    if (se) {
      mse_se <- stats::sd(resid^2, na.rm = TRUE) / sqrt(n_obs)
      result$rmse_se <- mse_se / (2 * result$rmse)
    }
  }
  if ("mae" %in% metrics) {
    result$mae <- mean(abs(resid), na.rm = TRUE)
    if (se) result$mae_se <- stats::sd(abs(resid), na.rm = TRUE) / sqrt(n_obs)
  }
  if ("r2" %in% metrics) {
    result$r2_cv        <- 1 - sum(resid^2, na.rm = TRUE) / sum((y_all - mean(y_all))^2, na.rm = TRUE)
    result$r2_cv_gelman <- stats::var(mu_all) / (stats::var(mu_all) + stats::var(resid))
    if (se) {
      fold_r2 <- vapply(fold_out, function(f) {
        1 - sum((f$y - f$mu)^2) / sum((f$y - mean(f$y))^2)
      }, numeric(1))
      result$r2_cv_se <- stats::sd(fold_r2, na.rm = TRUE) / sqrt(length(fold_r2))
    }
  }
  
  if (!is.null(backtransform_fn)) {
    y_raw     <- backtransform_fn(y_all)
    mu_raw    <- backtransform_fn(mu_all)
    resid_raw <- y_raw - mu_raw
    if ("rmse" %in% metrics) {
      mse_raw <- mean(resid_raw^2, na.rm = TRUE)
      result$rmse_raw <- sqrt(mse_raw)
      if (se) result$rmse_raw_se <- (stats::sd(resid_raw^2, na.rm = TRUE) / sqrt(n_obs)) / (2 * result$rmse_raw)
    }
    if ("mae" %in% metrics) {
      result$mae_raw <- mean(abs(resid_raw), na.rm = TRUE)
      if (se) result$mae_raw_se <- stats::sd(abs(resid_raw), na.rm = TRUE) / sqrt(n_obs)
    }
    if ("r2" %in% metrics) {
      result$r2_cv_raw        <- 1 - sum(resid_raw^2, na.rm = TRUE) / sum((y_raw - mean(y_raw))^2, na.rm = TRUE)
      result$r2_cv_gelman_raw <- stats::var(mu_raw) / (stats::var(mu_raw) + stats::var(resid_raw))
    }
  }
  
  if (return_predictions) {
    result$y_all  <- y_all
    result$mu_all <- mu_all
    if (!is.null(agg_cols)) {
      result$agg_data <- do.call(rbind, lapply(fold_out, `[[`, "agg"))
    }
  }
  
  result
}

# ----------------------------------------------------------------------------
# Native GEE confidence interval for the marginal mean, via glmgee's own
# cluster-robust (sandwich) SE -- predict.glmgee(se.fit=TRUE, varest="robust").
# No bootstrap needed here, unlike the glmmTMB case.
#
# This answers a DIFFERENT question than predict_marginal_interval(): uncertainty
# in the ESTIMATED MEAN (confidence interval), not the spread of a NEW
# observation (prediction interval). GEE has no natural analogue of the
# latter -- it doesn't model between-subject variance components at all -- so
# for a genuine prediction interval on a glmgee model, use cv_glmgee(...,
# return_predictions = TRUE) + cv_conformal_margin() + predict_marginal_interval()
# with predict_fn = .glmgee_default_predict_fn, same conformal machinery as glmmTMB.
#
# Built on the link scale then back-transformed through the endpoints (not
# fit +/- margin on the response scale), for the same reason as elsewhere in
# this file: keeps bounds valid (e.g. in [0,1] for binomial) and makes
# backtransform_fn handling exact.
# ----------------------------------------------------------------------------
predict_gee_ci <- function(model, newdata, level = 0.95,
                           varest = "robust", backtransform_fn = NULL) {
  pr <- stats::predict(model, newdata = newdata, type = "link",
                       se.fit = TRUE, varest = varest)
  z  <- stats::qnorm(1 - (1 - level) / 2)
  
  eta       <- pr[, 1]
  eta_lower <- eta - z * pr[, 2]
  eta_upper <- eta + z * pr[, 2]
  
  linkinv <- tryCatch(model$family$linkinv, error = function(e) NULL)
  if (is.null(linkinv)) stop("Could not find family$linkinv on this glmgee object.")
  
  fit   <- linkinv(eta)
  lower <- linkinv(eta_lower)
  upper <- linkinv(eta_upper)
  
  if (!is.null(backtransform_fn)) {
    fit_out   <- backtransform_fn(fit)
    lower_out <- backtransform_fn(lower)
    upper_out <- backtransform_fn(upper)
  } else {
    fit_out <- fit; lower_out <- lower; upper_out <- upper
  }
  
  data.frame(fit = fit_out, lower = lower_out, upper = upper_out)
}



#
#   cv_glmmTMB_sanity_check(m1)   # run once per model before trusting elpd
#
#   res <- cv_glmmTMB(m1, K = 6, metrics = c("elpd", "rmse", "mae", "r2"),
#                      se = TRUE, control = glmmTMBctrl)
#   res
#
# If your response is a rescaled version of a raw bounded variable (e.g. a
# [-10, 10] variable transformed to [0, 1] for ordbeta), elpd and R2 need no
# adjustment, but RMSE/MAE are reported on the transformed scale by default.
# Supply the inverse transform to also get them in original units:
#
#   res <- cv_glmmTMB(m1, K = 6, se = TRUE, control = glmmTMBctrl,
#                      backtransform_fn = function(y01) y01 * 20 - 10)
#   res$rmse_raw; res$mae_raw   # in original [-10, 10] units
#
# Marginal, population-averaged prediction intervals, back-transform-safe:
#
#   # 1. Run CV once, keeping the pooled out-of-fold predictions
#   cv_res <- cv_glmmTMB(m1, K = 6, control = glmmTMBctrl,
#                         return_predictions = TRUE)
#
#   # 2. Derive the conformal margin from CV residuals (cheap, no refitting)
#   cm <- cv_conformal_margin(cv_res, level = 0.95)
#   cm$margin; cm$empirical_coverage   # sanity check: should be ~ 0.95
#
#   # 3. Fit the FINAL model on all data, then get intervals for new data
#   m1_final <- glmmTMB::glmmTMB(m1formula, dispformula = ~ 1 + AmbientEnv,
#                                 data = m1Data, family = ordbeta(link = "logit"),
#                                 control = glmmTMBctrl)
#   predict_marginal_interval(m1_final, newdata = new_data, margin = cm$margin,
#                              backtransform_fn = function(y01) y01 * 20 - 10)
#   # -> data.frame(fit, lower, upper) in original [-10, 10] units
#
# Checking whether the re.form = NA "b = 0 plug-in" default matters for your
# model (see https://github.com/vincentarelbundock/marginaleffects/issues/1231):
#
#   mean(predict(m1, type = "response", re.form = NA))
#   mean(model.frame(m1)[[1]])   # raw sample mean, for comparison
#
# If those diverge meaningfully, swap in Monte Carlo integration everywhere
# (both where the conformal margin is calibrated and where it's applied --
# same predict_fn in both, since the margin is only valid for the process it
# was calibrated against):
#
#   cv_res <- cv_glmmTMB(m1, K = 6, control = glmmTMBctrl,
#                         return_predictions = TRUE,
#                         predict_fn = mc_marginal_predict_fn)
#   cm <- cv_conformal_margin(cv_res, level = 0.95)
#   predict_marginal_interval(m1_final, newdata = new_data, margin = cm$margin,
#                              predict_fn = mc_marginal_predict_fn,
#                              backtransform_fn = function(y01) y01 * 20 - 10)
#
# Model comparison (same seed -> same folds, so differences reflect the
# models, not the split):
#
#   res1 <- cv_glmmTMB(m1, K = 6, se = TRUE, seed = 123, control = glmmTMBctrl)
#   res2 <- cv_glmmTMB(m2, K = 6, se = TRUE, seed = 123, control = glmmTMBctrl)
#   data.frame(model = c("m1", "m2"),
#              elpd  = c(res1$elpd, res2$elpd),
#              rmse  = c(res1$rmse, res2$rmse))
# ============================================================================

# ============================================================================
# VALIDATION HARNESS: which CI method recovers the KNOWN true aggregate
# sampling variance?
#
# The blockage this resolves: simulate_marginal_ci() and
# cluster_boot_avg_predictions() disagree by a small, stable margin, and we
# have no external reference to say which is right -- the bootstrap "feels"
# trustworthy (real refits) but is suspected of running narrow due to
# cloned-subject shrinkage (resampled duplicate subjects get near-identical
# BLUPs, deflating between-subject spread), while the simulation is suspected
# of running wide for reasons that survived several targeted checks.
#
# This breaks the deadlock by generating data from KNOWN parameters (m1
# itself), where the true per-group aggregate sampling distribution can be
# obtained directly from simulated responses -- no refitting, no re.form
# choice, no plug-in bias -- and compared against what each method CLAIMS.
#
# Approach:
#  (A) TRUTH: simulate n_truth fresh datasets' worth of RAW RESPONSES from
#      `model` (known parameters), each with independent new random effects.
#      The per-group aggregate mean of these raw responses, across many
#      simulated datasets, IS the true sampling distribution of the
#      aggregate, by construction -- no model-fitting involved.
#  (B) Run each candidate method ONCE on the real data and extract its implied
#      per-group SD (se_sim for simulation; sd_boot for bootstrap).
#  (C) Compare each method's implied SD to the truth from (A).
#
# CAVEAT (not fully resolved): truth here includes residual/observation-level
# noise (a PREDICTION-interval-type target), while simulate_marginal_ci() and
# cluster_boot_avg_predictions() both target uncertainty of the FITTED MEAN
# (a CONFIDENCE-interval-type target, no residual noise). These differ --
# truth_sd is likely somewhat inflated relative to either method's actual
# target as a result. This version fixes a more serious flaw (an earlier
# iteration used re.form = NA per refit, contaminating truth with the same
# Jensen's-gap bias documented earlier in this file) but doesn't claim to be
# a fully resolved, apples-to-apples validation.
#
# No refitting needed -- n_truth can be large (hundreds) cheaply; the truth
# SD itself has Monte Carlo error ~ 1/sqrt(n_truth), so more is better here.
#
# Returns per-group truth SD, each method's implied SD, and their ratios.
# Ratio ~1 = that method matches this truth reference; >1 = wider than it;
# <1 = narrower -- interpret via the caveat above, not as absolute truth.
# ============================================================================

validate_ci_methods <- function(model, full_data, by, id_col,
                                se_sim_vec, sd_boot_vec,
                                n_truth = 200, seed = 123) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(inherits(model, "glmmTMB"))
  
  resp <- insight::find_response(model)
  
  # simulate() with re.form = NULL draws FRESH random effects AND realistic
  # response values from the model's TRUE (known) parameters -- no refitting,
  # no re.form choice at prediction time, no plug-in bias. The aggregate mean
  # of these raw simulated responses, across many simulated datasets, IS the
  # true sampling distribution of the aggregate by construction.
  #
  # IMPORTANT CAVEAT, not fully resolved: this "truth" includes residual/
  # observation-level noise (a PREDICTION-interval-type target -- uncertainty
  # of a future REALIZED aggregate), whereas simulate_marginal_ci() and
  # cluster_boot_avg_predictions() both estimate uncertainty of the FITTED
  # MEAN itself (a CONFIDENCE-interval-type target, no residual noise added).
  # These are genuinely different quantities, and truth_sd here is likely
  # somewhat LARGER than either method's target as a result -- expect both
  # ratio_sim and ratio_boot to shift down relative to the previous (flawed)
  # run, not just ratio_sim. Getting a fully apples-to-apples confidence-type
  # truth reference (re-estimated fixed effects AND correctly-marginalized
  # aggregate per simulated dataset, without reusing the same mechanism as
  # the method being validated) is a harder problem than this fixes -- treat
  # this as a substantial improvement on the previous version, not a fully
  # resolved one.
  sim_responses <- stats::simulate(model, nsim = n_truth, re.form = NULL)
  
  truth_out <- lapply(seq_len(n_truth), function(s) {
    df <- full_data[, by, drop = FALSE]
    df$y <- sim_responses[[s]]
    grp_key <- interaction(df[by], drop = TRUE)
    agg <- tapply(df$y, grp_key, mean)
    agg[order(names(agg))]
  })
  
  truth_mat <- do.call(rbind, truth_out)             # n_truth x n_groups
  truth_sd  <- apply(truth_mat, 2, stats::sd)        # true aggregate SD per group
  
  ng <- length(truth_sd)
  if (length(se_sim_vec) != ng || length(sd_boot_vec) != ng) {
    warning("Length mismatch: truth has ", ng, " groups but se_sim_vec/sd_boot_vec have ",
            length(se_sim_vec), "/", length(sd_boot_vec),
            ". Ensure all three are sorted by the same `by` ordering before comparing.")
  }
  
  data.frame(
    truth_sd       = truth_sd,
    se_sim         = se_sim_vec,
    sd_boot        = sd_boot_vec,
    ratio_sim      = se_sim_vec  / truth_sd,
    ratio_boot     = sd_boot_vec / truth_sd,
    row.names = NULL
  )
}

# Usage:
#   # se_sim_vec: ci_result$se_sim (from simulate_marginal_ci, native scale)
#   # sd_boot_vec: boot_result$sd_boot (from cluster_boot_avg_predictions)
#   # All THREE must be sorted by the same `by`/StimID ordering.
#
#   val <- validate_ci_methods(m1, full_data = nd, by = "StimID", id_col = "ID",
#                              se_sim_vec  = ci_result$se_sim[order(ci_result$StimID)],
#                              sd_boot_vec = boot_result$sd_boot[order(boot_result$StimID)],
#                              n_truth = 200)   # no refitting, so this is cheap -- go big
#
#   summary(val$ratio_sim)
#   summary(val$ratio_boot)
#   # interpret RELATIVE to each other first (which is closer to 1), not the
#   # absolute values -- see the CAVEAT above about what truth_sd does and
#   # doesn't include relative to what these two methods actually target.