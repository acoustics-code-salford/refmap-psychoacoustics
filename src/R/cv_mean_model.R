require(future)         # future::plan() for parallel = TRUE (set the plan yourself)
require(future.apply)   # future.apply::future_lapply() for parallel = TRUE
require(glmmTMB)        # glmmTMB models
require(glmtoolbox)     # glmgee models (only needed if you cross-validate one)

# =============================================================================
# cv_mean_model.R  (v4)
#
# Cross-validated RMSE and R^2 for a stimulus-level model, for
#   (i)  the model: E[y] for each held-out stimulus at its own predictor values
#        (including its own experiment's contrast code), and, when the input is
#        an object from deploy_mean_model(),
#   (ii) the deployed function: E[y] through the deployed path (deployment
#        constants such as Experiment = 0, centring, frozen terms),
# from the same folds and refits. Source deploy_mean_model.R (v17+) and
# boot_mean_model.R (v14+) first: this file uses their internal functions
# (the refit worker, the design matrices and the mean map).
#
# -----------------------------------------------------------------------------
# USAGE
#   cv_mean_model(object, k = 10L, reps = 10L, strata = NULL, group = NULL,
#                 inverse = NULL, seed = NULL, parallel = FALSE, verbose = TRUE,
#                 folds = NULL, data = NULL)
#   compare_mean_models(...)   # existing cv_mean_model() results, named
#
# ARGUMENTS
#   object    Either an object returned by deploy_mean_model() (v17+): results
#             for (i) the model and (ii) the deployed function; or a fitted
#             glmmTMB or glmgee model: results for (i) only. Zero-inflated models
#             are refused.
#
#   k         Number of folds (default 10). k = number of units gives
#             leave-one-out (deterministic: reps is set to 1).
#
#   reps      Number of repeats of k-fold CV with new random folds (default 10).
#             Results are averaged over repeats; their spread is reported.
#
#   strata    NULL (default) or column name(s) in the model's stimulus-level
#             data. Folds are balanced within strata (e.g. "Experiment"), so every
#             fold's training data keep each stratum's factor levels. Each stratum
#             needs at least 2 units, and a unit (see group) must lie within one
#             stratum.
#
#   group     NULL (default: each stimulus row is a unit) or a column name in the
#             model's stimulus-level data. Rows with the same value (e.g. the same
#             sound at several levels, or the same sound in both experiments) are
#             always held out together.
#
#   inverse   NULL (identity) or a vectorised function taking values on the
#             model's response scale (the formula LHS as the model sees it) to
#             the reporting scale, e.g. function(p) 10 * p for a response stored
#             as median / 10. Applied to observed values and predictions alike
#             before RMSE and R^2 are computed.
#
#   seed      NULL (current random-number stream) or a number (reproducible; the
#             global RNG state is restored afterwards).
#
#   parallel  TRUE to refit with future.apply::future_lapply() under the current
#             future::plan(). Folds are drawn in the main session, so results do
#             not depend on the plan.
#
#   verbose   Print progress and the summary.
#
#   folds     NULL (default: drawn here) or a fold assignment to reuse: the $folds
#             matrix of an earlier result (one row per unit, named by unit; one
#             column per repeat). k, reps and seed are then ignored. An
#             alternative to giving each run the same seed, k, reps, strata and
#             group, which also gives the same folds for models fitted to the
#             same units.
#
#   data      Model input only (ignored for a deployed object): the data frame the
#             model was fitted to. NULL (default): taken from the model's call
#             (its data argument, evaluated in the formula's environment); give
#             it if that is not available. Its rows are matched to the model's
#             rows by row name, so pass it unchanged.
#
# VALUE (class "cv_mean_model")
#   $summary    for target "model" (and "deployed"): RMSE and R^2, apparent (full
#               fit, in-sample), CV mean over repeats, SD/min/max over repeats,
#               and difference (CV mean - apparent).
#   $by_stratum the same by stratum (if strata given); R^2 there uses the
#               stratum's own mean.
#   $oof        out-of-fold predictions on the reporting scale (rows = stimulus
#               rows used, cols = repeats): model, deployed (NULL for a model
#               input); and the model's
#               out-of-fold linear predictor (eta_model, for checks).
#   $observed   observed values on the reporting scale, with row, unit and
#               stratum; $apparent apparent
#               predictions; $folds fold ids
#               (units x repeats); $failures; $settings; $checks.
#
# -----------------------------------------------------------------------------
# HOW IT WORKS
#   Units (stimulus rows, or groups of rows) are assigned to k folds, within
#   strata, in each repeat. For each fold the model's original call is refitted
#   on the remaining rows (the same refit worker as boot_mean_model, with its
#   checks: convergence, coefficient names, design matrix). The held-out rows
#   are then predicted with the fold's coefficients (and ordbeta cutpoints):
#     (i)  model:    eta = X_model %*% beta, X_model = the model's own design
#                    (each row at its own factor levels; for a model input,
#                    the model's stored design matrix);
#     (ii) deployed: eta = X_dep %*% beta, X_dep = the deployed design (checked
#                    to reproduce the deployed function's linear predictor);
#   E[y] uses the full mean map (ordbeta: E = P1 + (P(y > 0) - P1) * mu), so it
#   does not depend on glmmTMB's predict(type = "response").
#   RMSE = sqrt(mean((y - yhat)^2)); R^2 = 1 - sum((y - yhat)^2) /
#   sum((y - mean(y))^2) (Kvalseth's R^2_1), computed once per repeat from the
#   pooled out-of-fold predictions (each row predicted once per repeat), not
#   averaged over folds.
#
# INTERPRETATION
#   * CV estimates are for predicting NEW STIMULI, with responses aggregated over
#     the same panel of participants (held-out medians share participants with
#     the training medians), so they are somewhat optimistic for a new panel.
#   * (ii) is scored against each stimulus's observed value, so the difference
#     between experiments that the deployment averages over counts as error:
#     what a user of the deployed equation would see.
#   * For models with random effects, (i) and (ii) are fixed-effects
#     (population-level) predictions.
#   * Model input: formula terms that depend on the data (scale(), poly(), ...)
#     are recomputed by each refit on its training rows, while the held-out rows
#     use the full-data design; the per-refit design check discards such
#     refits (many "design does not match" failures). Deploy the model (which
#     freezes such terms) to cross-validate it.
#   * Apparent values should reproduce in-sample RMSE / R^2 computed elsewhere
#     from the same predictions (a check on inverse and the designs).
#
# CHECKS (stops)
#   deploy/boot internals loaded; zero-inflated models; the model's own design
#   must reproduce its linear predictor (deployed input; glmgee model input),
#   and the deployed design the deployed function's; model input: the model's
#   rows must be found in the data by row name; inverse must be vectorised and give finite values; k, reps;
#   missing values in strata/group columns; a group spanning several strata;
#   strata with fewer than 2 units; models with offsets. Failed refits are
#   counted (warning); their rows are left out of that repeat's metrics.
#
# Not verified by execution (no R runtime was available when this was written).
#
# EXAMPLE
#   source("deploy_mean_model.R"); source("boot_mean_model.R"); source("cv_mean_model.R")
#   cvm <- cv_mean_model(dep, k = 10, reps = 10, strata = "Experiment", group = "StimNum",
#                        inverse = function(p) 10 * p, seed = 1)
#   cvm$summary
#
# -----------------------------------------------------------------------------
# compare_mean_models(...)
#   Compiles EXISTING cv_mean_model() results into one table (nothing is
#   refitted). Give the results as named arguments, compare_mean_models(m1 = cv1,
#   m2 = cv2), or as one named list; unnamed arguments are labelled by their
#   expressions. Results from deployed objects and from plain models can be
#   mixed: the deployed columns (ii) are NA for the latter. One row per model:
#   family, number of conditional-model coefficients, and for the model (i)
#   and the deployed function (ii):
#   apparent and CV RMSE and R^2, the SD of CV RMSE over repeats, and, when the
#   runs share their folds and data, the paired difference in CV RMSE from the
#   best model (mean over repeats) and the share of repeats in which the model
#   has the lowest CV RMSE.
#   Warnings (the table is still produced) when the runs differ in k, reps,
#   seed, strata, group, leave-one-out or use of inverse; when their fold
#   assignments differ (e.g. same seed but different units); or when their
#   stimuli or observed values (reporting scale) differ. The paired columns are
#   NA unless folds and observed values are identical across all runs; for
#   (ii) they compare the runs that have deployed results (at least 2 needed).
#   The SD over repeats reflects fold-to-fold randomness only, not sampling
#   uncertainty: small differences are not evidence that one model generalises
#   better.
#   Results from cv_mean_model() v1-v3 are accepted.
#   VALUE (class "mean_model_comparison"): $table, $checks.
#   EXAMPLE
#     cv1 <- cv_mean_model(dep1, k = 10, reps = 10, strata = "Experiment", group = "StimNum",
#                          inverse = function(p) 10 * p, seed = 1)
#     cv2 <- cv_mean_model(dep2, k = 10, reps = 10, strata = "Experiment", group = "StimNum",
#                          inverse = function(p) 10 * p, seed = 1)
#     cmp <- compare_mean_models(m_dA1 = cv1, m_dA2 = cv2)
#     cmp$table
#
# CHANGE LOG
#   v4.1: compare_mean_models(): paired (ii) columns computed among the deployed
#       runs when the comparison mixes deployed and plain-model results.
#   v4: cv_mean_model() also takes a fitted glmmTMB/glmgee model (results for the
#       model only; new data argument); first argument renamed object. AIC
#       removed (CV answers the out-of-sample question directly);
#       compare_mean_models() accepts mixed deployed/model results.
#   v3: compare_mean_models() now compiles existing results (no refitting),
#       warning about differing settings, folds or data; AIC stored in
#       $settings$aic.
#   v2: folds argument (reuse a fold assignment); $observed also on the model
#       scale; compare_mean_models().
# =============================================================================

cv_mean_model <- function(object, k = 10L, reps = 10L, strata = NULL, group = NULL,
                          inverse = NULL, seed = NULL, parallel = FALSE, verbose = TRUE, folds = NULL,
                          data = NULL) {
  
  # ---- deploy_mean_model.R and boot_mean_model.R must be loaded -------------
  need_fn <- c(".boot_refit_one", ".boot_design", ".boot_mean_map", ".eval_eta_codes", ".apply_scaling",
               ".ensure_raw_columns", ".aligned_rows", ".get_inverse_link")
  miss_fn <- need_fn[!vapply(need_fn, exists, logical(1), mode = "function")]
  if (length(miss_fn) > 0) {
    stop("Source deploy_mean_model.R (v17+) and boot_mean_model.R (v14+) first; missing: ",
         paste(miss_fn, collapse = ", "), call. = FALSE)
  }
  
  # ---- Model and (if deployed) the deployment ------------------------------
  has_dep <- inherits(object, "deployed_mean_model")
  if (has_dep) {
    model <- attr(object, "model")
    if (is.null(model)) stop("The deployed object has no attr(, \"model\"): re-run deploy_mean_model() (v17 or later).", call. = FALSE)
  } else {
    model <- object
  }
  cls <- if (inherits(model, "glmmTMB")) "glmmTMB" else if (inherits(model, "glmgee")) "glmgee" else class(model)[1]
  if (!cls %in% c("glmmTMB", "glmgee")) {
    stop(sprintf("object must be a deploy_mean_model() object or a glmmTMB/glmgee model, not '%s'.", cls), call. = FALSE)
  }
  if (cls == "glmmTMB" && length(tryCatch(glmmTMB::fixef(model)$zi, error = function(e) NULL)) > 0) {
    stop("Zero-inflated models are not supported.", call. = FALSE)
  }
  cf_names <- if (cls == "glmmTMB") names(glmmTMB::fixef(model)$cond) else rownames(model[["coefficients"]])
  f_full <- if (cls == "glmmTMB") stats::formula(model) else stats::formula(model[["terms"]])
  target <- all.vars(f_full[[2]])
  if (length(target) != 1L) stop("Only a single response column is supported.", call. = FALSE)
  env_eval <- tryCatch(environment(f_full), error = function(e) NULL)
  if (!is.environment(env_eval)) env_eval <- globalenv()
  fit_call <- stats::getCall(model)
  if (is.null(fit_call)) stop("Could not retrieve the model's call.", call. = FALSE)
  
  .scalar_num <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!.scalar_num(k) || k < 2 || k != round(k)) stop("k must be a single whole number >= 2.", call. = FALSE)
  if (!.scalar_num(reps) || reps < 1 || reps != round(reps)) stop("reps must be a single whole number >= 1.", call. = FALSE)
  k <- as.integer(k); reps <- as.integer(reps)
  
  dev_model <- NA_real_; dev_dep <- NA_real_; X_dep <- NULL
  if (has_dep) {
    # ---- deployed input: spec, stimulus data and designs from the deployment ----
    env <- environment(object)
    needed <- c("spec", "scaling", "numeric_codes", "constants_levels", "eta_with_constants")
    miss <- needed[!vapply(needed, exists, logical(1), envir = env, inherits = FALSE)]
    if (length(miss) > 0) stop(sprintf("The deployed function's environment lacks: %s.", paste(miss, collapse = ", ")), call. = FALSE)
    spec <- env$spec
    scaling <- env$scaling
    numeric_codes <- env$numeric_codes
    constants_levels <- env$constants_levels
    inv_link <- .get_inverse_link(spec$link)
    fam <- if (cls == "glmmTMB") spec$family_name else NA_character_
    ordb <- identical(fam, "ordbeta")
    psi0 <- if (ordb) matrix(spec$mean_info$psi, nrow = 1) else NULL
    cf0 <- spec$coefficients
    if (length(cf_names) != length(cf0)) stop("Coefficient count differs between the model and the deployed spec.", call. = FALSE)
    if (!isTRUE(spec$data_is_raw) || !is.data.frame(spec$data)) {
      stop("The model's stimulus-level data (with raw columns) could not be recovered by deploy_mean_model().", call. = FALSE)
    }
    stim <- .ensure_raw_columns(spec$data, scaling)
    rownames(stim) <- NULL
    idx0 <- .aligned_rows(spec)
    
    data_model <- .apply_scaling(stim, scaling)
    X_model <- .boot_design(spec, data_model, list())
    if (!is.null(idx0) && is.numeric(spec$fitted)) {
      dev_model <- max(abs(as.numeric(X_model[idx0, , drop = FALSE] %*% cf0) - spec$fitted))
      if (dev_model > 1e-7 * max(1, max(abs(spec$fitted)))) {
        stop(sprintf("The model's design matrix does not reproduce its own linear predictor (max|dev| = %.3g).", dev_model), call. = FALSE)
      }
    }
    data_dep <- data_model
    for (v in names(constants_levels)) {
      data_dep[[v]] <- if (!is.null(spec$xlevels[[v]])) as.character(constants_levels[[v]]) else constants_levels[[v]]
    }
    X_dep <- .boot_design(spec, data_dep, numeric_codes)
    eta_dep <- env$eta_with_constants(data_dep)
    dev_dep <- max(abs(as.numeric(X_dep %*% cf0) - eta_dep), na.rm = TRUE)
    if (dev_dep > 1e-7 * max(1, max(abs(eta_dep), na.rm = TRUE))) {
      stop(sprintf("The deployed design matrix does not reproduce the deployed linear predictor (max|dev| = %.3g).", dev_dep), call. = FALSE)
    }
  } else {
    # ---- model input: everything from the fitted model and its data ----------
    if (cls == "glmmTMB") {
      fam <- stats::family(model)$family
      inv_link <- stats::family(model)$linkinv
      cf0 <- glmmTMB::fixef(model)$cond
      frame <- model[["frame"]]
      Xf <- tryCatch(as.matrix(glmmTMB::getME(model, "X")), error = function(e) NULL)
    } else {
      fam <- NA_character_
      inv_link <- tryCatch(model[["family"]][["linkinv"]], error = function(e) NULL)
      cf0 <- stats::setNames(as.numeric(model[["coefficients"]]), cf_names)
      frame <- model[["model"]]
      Xf <- tryCatch(as.matrix(stats::model.matrix(model)), error = function(e) NULL)
    }
    if (!is.function(inv_link)) stop("Could not read the model's inverse link.", call. = FALSE)
    ordb <- identical(fam, "ordbeta")
    psi0 <- if (ordb) {
      matrix(model[["obj"]][["env"]][["parList"]](model[["fit"]][["par"]], model[["fit"]][["parfull"]])$psi, nrow = 1)
    } else NULL
    if (ordb && (length(psi0) != 2L || !all(is.finite(psi0)))) stop("Could not read the ordbeta cutpoints.", call. = FALSE)
    stim <- if (is.null(data)) tryCatch(eval(fit_call$data, env_eval), error = function(e) NULL) else data
    if (!is.data.frame(stim)) {
      stop("Could not recover the data the model was fitted to; pass it as data = (the same data frame, unchanged).", call. = FALSE)
    }
    stim <- as.data.frame(stim)
    idx0 <- match(rownames(frame), rownames(stim))
    if (is.null(frame) || anyNA(idx0)) {
      stop("The model's rows could not be matched to the data's rows by row name; pass data = the data frame the model was fitted to, unchanged.", call. = FALSE)
    }
    rownames(stim) <- NULL
    if (is.null(Xf) || nrow(Xf) != length(idx0) || ncol(Xf) != length(cf0)) {
      stop("Could not obtain the model's design matrix for its rows.", call. = FALSE)
    }
    X_model <- matrix(NA_real_, nrow(stim), ncol(Xf), dimnames = list(NULL, colnames(Xf)))
    X_model[idx0, ] <- Xf
    if (cls == "glmgee") {
      lp <- as.numeric(model[["linear.predictors"]])
      if (length(lp) == length(idx0)) {
        dev_model <- max(abs(as.numeric(Xf %*% cf0) - lp))
        if (dev_model > 1e-7 * max(1, max(abs(lp)))) {
          stop(sprintf("The model's design matrix does not reproduce its linear predictor (max|dev| = %.3g).", dev_model), call. = FALSE)
        }
      }
    }
  }
  
  # ---- Response as the model sees it, and on the reporting scale ------------
  y_model <- suppressWarnings(as.numeric(eval(f_full[[2]], stim, env_eval)))
  pool <- if (is.null(idx0)) seq_len(nrow(stim)) else idx0          # rows the model used
  pool <- pool[!is.na(y_model[pool])]
  inv <- if (is.null(inverse)) identity else inverse
  if (!is.function(inv)) stop("inverse must be a function or NULL.", call. = FALSE)
  y_obs <- tryCatch(as.numeric(inv(y_model[pool])), error = function(e) e)
  if (inherits(y_obs, "error") || length(y_obs) != length(pool) || !all(is.finite(y_obs))) {
    stop("inverse must be vectorised and give finite values for the observed responses.", call. = FALSE)
  }
  
  # ---- Units and strata -------------------------------------------------------
  if (!is.null(group)) {
    if (!(is.character(group) && length(group) == 1L && group %in% names(stim))) {
      stop("group must be the name of one column in the model's stimulus-level data.", call. = FALSE)
    }
    if (anyNA(stim[[group]][pool])) stop(sprintf("The group column '%s' has missing values.", group), call. = FALSE)
    unit <- as.character(stim[[group]][pool])
  } else {
    unit <- as.character(pool)
  }
  units <- unique(unit)
  n_u <- length(units)
  if (!is.null(strata)) {
    miss_s <- setdiff(strata, names(stim))
    if (length(miss_s) > 0) stop(sprintf("The stimulus-level data is missing strata column(s): %s", paste(miss_s, collapse = ", ")), call. = FALSE)
    if (any(vapply(strata, function(v) anyNA(stim[[v]][pool]), logical(1)))) stop("The strata columns have missing values.", call. = FALSE)
    s_row <- do.call(paste, c(lapply(strata, function(v) as.character(stim[[v]][pool])), sep = "/"))
  } else {
    s_row <- rep("all", length(pool))
  }
  s_list <- split(s_row, factor(unit, levels = units))
  if (any(lengths(lapply(s_list, unique)) > 1L)) stop("A group spans several strata; groups must lie within one stratum.", call. = FALSE)
  s_unit <- vapply(s_list, function(v) v[1], character(1), USE.NAMES = FALSE)
  if (any(table(s_unit) < 2L)) stop("Every stratum needs at least 2 units.", call. = FALSE)
  if (!is.null(folds)) {
    folds <- as.matrix(folds)
    if (is.null(rownames(folds)) || nrow(folds) != n_u || !setequal(rownames(folds), units)) {
      stop(paste0("folds must have one row per unit, named by unit, for the same units as this model ",
                  "(e.g. $folds of a cv_mean_model() result with the same data and group)."), call. = FALSE)
    }
    if (anyNA(folds) || !is.numeric(folds) || any(folds != round(folds))) stop("folds must contain whole-number fold ids.", call. = FALSE)
    folds <- folds[units, , drop = FALSE]
    storage.mode(folds) <- "integer"
    reps <- ncol(folds)
    k <- max(apply(folds, 2, function(f) length(unique(f))))
  }
  if (k > n_u) stop(sprintf("k = %d exceeds the number of units (%d).", k, n_u), call. = FALSE)
  loo <- k == n_u
  if (is.null(folds) && loo && reps > 1L) {
    if (verbose) cat("Leave-one-out (k = number of units) is deterministic: reps set to 1.\n")
    reps <- 1L
  }
  
  # ---- Folds (all random numbers drawn here) ----------------------------------
  if (is.null(folds) && !is.null(seed)) {
    seed_old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
    on.exit({
      if (!is.null(seed_old)) assign(".Random.seed", seed_old, envir = globalenv())
      else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
    }, add = TRUE)
    set.seed(seed)
  }
  # units shuffled within strata, strata concatenated, folds assigned cyclically
  # from a random start: each fold gets about k-th of every stratum
  if (is.null(folds)) folds <- matrix(vapply(seq_len(reps), function(r) {
    if (loo) return(seq_len(n_u))
    ord <- unlist(lapply(split(seq_len(n_u), s_unit), function(v) v[sample.int(length(v))]), use.names = FALSE)
    f <- integer(n_u)
    f[ord] <- (seq_along(ord) + sample.int(k, 1L) - 2L) %% k + 1L
    f
  }, integer(n_u)), nrow = n_u, dimnames = list(units, NULL))
  fold_row <- folds[match(unit, units), , drop = FALSE]                 # pool rows x reps
  
  # ---- Refit call (as in boot_mean_model) --------------------------------------
  has_offset_term <- tryCatch(!is.null(attr(stats::terms(f_full), "offset")), error = function(e) FALSE)
  if (!is.null(fit_call$offset) || has_offset_term || any(grepl("offset(", deparse(f_full), fixed = TRUE))) {
    stop("Models with an offset are not supported.", call. = FALSE)
  }
  fit_call$data <- quote(.boot_data)
  fit_call[[1]] <- if (cls == "glmmTMB") quote(glmmTMB::glmmTMB) else quote(glmtoolbox::glmgee)
  free <- setdiff(unique(unlist(lapply(as.list(fit_call)[-1], all.vars))), c(names(stim), ".boot_data"))
  ext <- list()
  for (v in free) if (exists(v, envir = env_eval)) ext[[v]] <- get(v, envir = env_eval)
  pkgs <- if (cls == "glmmTMB") "glmmTMB" else "glmtoolbox"
  
  jobs <- unlist(lapply(seq_len(reps), function(r) {
    lapply(sort(unique(fold_row[, r])), function(f) list(r = r, test = which(fold_row[, r] == f)))
  }), recursive = FALSE)
  items <- lapply(jobs, function(j) list(y = stim[[target]], rows = pool[-j$test]))
  worker_args <- list(stim = stim, target = target, fit_spec = list(call = fit_call), ext = ext, cls = cls, fam = fam,
                      cf_names = cf_names, X_train = X_model, want_disp = FALSE)
  t0 <- Sys.time()
  if (verbose) cat(sprintf("Refitting %d folds x %d repeat(s) = %d fits (%s)...\n", k, reps, length(jobs),
                           if (parallel) "future.apply" else "serial"))
  res <- if (parallel) {
    if (!requireNamespace("future.apply", quietly = TRUE)) stop("parallel = TRUE needs the future.apply package.", call. = FALSE)
    do.call(future.apply::future_lapply, c(list(X = items, FUN = .boot_refit_one), worker_args,
                                           list(future.packages = pkgs, future.seed = FALSE)))
  } else {
    do.call(lapply, c(list(X = items, FUN = .boot_refit_one), worker_args))
  }
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  # ---- Out-of-fold predictions -------------------------------------------------
  n_p <- length(pool)
  H_model <- P_model <- P_dep <- matrix(NA_real_, n_p, reps)
  reasons <- character(0)
  for (j in seq_along(jobs)) {
    out <- res[[j]]
    if (!isTRUE(out$ok)) {
      reasons <- c(reasons, if (is.null(out$reason) || is.na(out$reason)) "unknown" else out$reason)
      next
    }
    te <- jobs[[j]]$test; r <- jobs[[j]]$r
    psi <- if (ordb) matrix(out$psi, nrow = 1) else NULL
    eta_m <- X_model[pool[te], , drop = FALSE] %*% out$cf
    H_model[te, r] <- as.numeric(eta_m)
    P_model[te, r] <- as.numeric(.boot_mean_map(eta_m, inv_link, psi))
    if (has_dep) P_dep[te, r] <- as.numeric(.boot_mean_map(X_dep[pool[te], , drop = FALSE] %*% out$cf, inv_link, psi))
  }
  if (length(reasons) > 0) {
    warning(sprintf("%d of %d refits failed; their held-out rows are left out of that repeat's metrics. See $failures.",
                    length(reasons), length(jobs)), call. = FALSE)
  }
  to_rep <- function(M) {
    v <- tryCatch(as.numeric(inv(as.vector(M))), error = function(e) e)
    if (inherits(v, "error") || length(v) != length(M)) stop("inverse failed on the predictions.", call. = FALSE)
    matrix(v, nrow = nrow(M))
  }
  O_model <- to_rep(P_model)
  A_model <- as.numeric(to_rep(.boot_mean_map(X_model[pool, , drop = FALSE] %*% cf0, inv_link, psi0)))
  O_dep <- NULL; A_dep <- NULL
  if (has_dep) {
    O_dep <- to_rep(P_dep)
    A_dep <- as.numeric(to_rep(.boot_mean_map(X_dep[pool, , drop = FALSE] %*% cf0, inv_link, psi0)))
  }
  if (any(is.finite(c(P_model, if (has_dep) P_dep)) & !is.finite(c(O_model, O_dep))) || !all(is.finite(c(A_model, A_dep)))) {
    stop("inverse gave non-finite values for finite predictions.", call. = FALSE)
  }
  
  # ---- Metrics -------------------------------------------------------------------
  metrics <- function(y, p) {
    ok <- is.finite(y) & is.finite(p)
    e <- y[ok] - p[ok]
    c(RMSE = sqrt(mean(e^2)), R2 = 1 - sum(e^2) / sum((y[ok] - mean(y[ok]))^2), n = sum(ok))
  }
  targets <- if (has_dep) list(model = list(A_model, O_model), deployed = list(A_dep, O_dep)) else list(model = list(A_model, O_model))
  tg_lab <- rep(names(targets), each = 2)
  tabulate_metrics <- function(rows) {
    do.call(rbind, lapply(targets, function(tg) {
      app <- metrics(y_obs[rows], tg[[1]][rows])
      cvm <- vapply(seq_len(reps), function(r) metrics(y_obs[rows], tg[[2]][rows, r]), numeric(3))
      data.frame(metric = c("RMSE", "R2"), apparent = app[1:2],
                 cv_mean = rowMeans(cvm[1:2, , drop = FALSE]),
                 cv_sd = if (reps > 1L) apply(cvm[1:2, , drop = FALSE], 1, stats::sd) else NA_real_,
                 cv_min = apply(cvm[1:2, , drop = FALSE], 1, min), cv_max = apply(cvm[1:2, , drop = FALSE], 1, max),
                 difference = rowMeans(cvm[1:2, , drop = FALSE]) - app[1:2],
                 n_min = min(cvm[3, ]), row.names = NULL)
    }))
  }
  summ <- tabulate_metrics(seq_len(n_p))
  summ <- cbind(target = tg_lab, summ)
  by_st <- NULL
  if (!is.null(strata)) {
    by_st <- do.call(rbind, lapply(sort(unique(s_row)), function(s) {
      tb <- tabulate_metrics(which(s_row == s))
      cbind(stratum = s, target = tg_lab, tb)
    }))
  }
  
  out <- list(
    summary = summ,
    by_stratum = by_st,
    oof = list(model = O_model, deployed = O_dep, eta_model = H_model),
    observed = data.frame(row = pool, unit = unit, stratum = s_row, observed = y_obs),
    apparent = if (has_dep) data.frame(model = A_model, deployed = A_dep) else data.frame(model = A_model),
    folds = folds,
    failures = list(n_failed = length(reasons), n_fits = length(jobs), reasons = sort(table(reasons), decreasing = TRUE)),
    settings = list(k = k, reps = reps, loo = loo, strata = strata, group = group, n_units = n_u, n_rows = n_p,
                    n_coef = length(cf0), deployed = has_dep,
                    strata_units = table(s_unit), inverse = !is.null(inverse), family = if (cls == "glmmTMB") fam else "glmgee",
                    constants = if (has_dep) attr(object, "constants") else NULL, seed = seed, parallel = parallel, elapsed_sec = elapsed),
    checks = list(model_design_dev = dev_model, deployed_design_dev = dev_dep)
  )
  class(out) <- "cv_mean_model"
  if (verbose) print(out)
  invisible(out)
}

print.cv_mean_model <- function(x, digits = 4, ...) {
  s <- x$settings
  cat("<cv_mean_model>\n")
  cat(sprintf("%s: %d units (%d stimulus rows)%s%s; %d fits, %d failed (%.0f s).\n",
              if (s$loo) "Leave-one-out" else sprintf("%d-fold CV x %d repeat(s)", s$k, s$reps),
              s$n_units, s$n_rows,
              if (!is.null(s$group)) sprintf(", units = '%s'", s$group) else "",
              if (!is.null(s$strata)) sprintf(", folds balanced within %s (%s)", paste(s$strata, collapse = "/"),
                                              paste(sprintf("%s: %d", names(s$strata_units), as.integer(s$strata_units)), collapse = "; ")) else "",
              x$failures$n_fits, x$failures$n_failed, s$elapsed_sec))
  if (length(s$constants) > 0) {
    cat(sprintf("Deployed predictions use the constants (%s); the model's use each stimulus's own values.\n",
                paste(sprintf("%s = %s", names(s$constants), unlist(s$constants)), collapse = ", ")))
  }
  cat(sprintf("Scale: %s.\n\n", if (s$inverse) "reporting scale (inverse applied)" else "model response scale"))
  print(x$summary, digits = digits, row.names = FALSE)
  if (!is.null(x$by_stratum)) {
    cat("\nBy stratum (R2 relative to the stratum's own mean):\n")
    print(x$by_stratum[, c("stratum", "target", "metric", "apparent", "cv_mean", "cv_sd", "difference")], digits = digits, row.names = FALSE)
  }
  cat("\nCV estimates are for new stimuli rated by the same panel (somewhat optimistic for a new panel).\n")
  if (x$failures$n_failed > 0) cat(sprintf("NOTE: %d failed refits; see $failures.\n", x$failures$n_failed))
  invisible(x)
}

compare_mean_models <- function(...) {
  res <- list(...)
  if (length(res) == 1L && is.list(res[[1]]) && !inherits(res[[1]], "cv_mean_model")) res <- res[[1]]
  nms <- names(res)
  if (is.null(nms)) nms <- rep("", length(res))
  labels <- vapply(as.list(substitute(list(...)))[-1], function(e) paste(deparse(e), collapse = ""), character(1))
  if (length(labels) == length(res)) nms[nms == ""] <- labels[nms == ""]
  nms[nms == ""] <- paste0("model", which(nms == ""))
  if (anyDuplicated(nms)) stop("Model names must be unique.", call. = FALSE)
  names(res) <- nms
  if (length(res) < 2L) stop("Give at least 2 cv_mean_model() results.", call. = FALSE)
  if (!all(vapply(res, inherits, logical(1), "cv_mean_model"))) stop("All inputs must be cv_mean_model() results.", call. = FALSE)
  
  # ---- consistency checks (warnings only) ------------------------------------
  setting <- function(r, nm) { v <- r$settings[[nm]]; if (is.null(v)) "NULL" else paste(v, collapse = "/") }
  diff_set <- character(0)
  for (nm in c("k", "reps", "seed", "strata", "group", "loo", "inverse")) {
    v <- vapply(res, setting, character(1), nm = nm)
    if (length(unique(v)) > 1L) diff_set <- c(diff_set, sprintf("%s (%s)", nm, paste(sprintf("%s: %s", nms, v), collapse = ", ")))
  }
  if (length(diff_set) > 0) {
    warning("The CV runs differ in: ", paste(diff_set, collapse = "; "), call. = FALSE)
  }
  ref <- res[[1]]
  same_folds <- all(vapply(res, function(r) identical(r$folds, ref$folds), logical(1)))
  same_obs <- all(vapply(res, function(r) identical(r$observed$unit, ref$observed$unit) &&
                           isTRUE(all.equal(r$observed$observed, ref$observed$observed, check.attributes = FALSE)), logical(1)))
  if (!same_obs) {
    warning("The runs differ in their stimuli or observed values (reporting scale): the comparison is not like for like.", call. = FALSE)
  } else if (!same_folds) {
    warning("The runs' fold assignments differ, so their CV estimates are not paired; paired columns are NA.", call. = FALSE)
  }
  has_dep <- vapply(res, function(r) !is.null(r$oof$deployed), logical(1))
  
  # ---- table -------------------------------------------------------------------
  num_setting <- function(r, nm) { v <- r$settings[[nm]]; if (is.null(v)) NA_real_ else as.numeric(v) }
  get_val <- function(r, target, metric, col) {
    v <- r$summary[r$summary$target == target & r$summary$metric == metric, col]
    if (length(v) == 1L) as.numeric(v) else NA_real_
  }
  tab <- data.frame(model = nms,
                    family = vapply(res, function(r) as.character(r$settings$family), character(1)),
                    n_coef = vapply(res, num_setting, numeric(1), nm = "n_coef"),
                    deployed = has_dep, stringsAsFactors = FALSE)
  paired <- same_folds && same_obs
  y <- ref$observed$observed
  for (tg in c("model", "deployed")) {
    tab[[paste0(tg, "_RMSE_app")]] <- vapply(res, get_val, numeric(1), target = tg, metric = "RMSE", col = "apparent")
    tab[[paste0(tg, "_RMSE_cv")]] <- vapply(res, get_val, numeric(1), target = tg, metric = "RMSE", col = "cv_mean")
    tab[[paste0(tg, "_RMSE_cv_sd")]] <- vapply(res, get_val, numeric(1), target = tg, metric = "RMSE", col = "cv_sd")
    tab[[paste0(tg, "_R2_app")]] <- vapply(res, get_val, numeric(1), target = tg, metric = "R2", col = "apparent")
    tab[[paste0(tg, "_R2_cv")]] <- vapply(res, get_val, numeric(1), target = tg, metric = "R2", col = "cv_mean")
    tab[[paste0(tg, "_dRMSE_cv")]] <- NA_real_
    tab[[paste0(tg, "_best_share")]] <- NA_real_
    # paired columns among the runs that have this target (deployed: only runs
    # from deployed objects), when at least 2 do
    use <- if (tg == "model") rep(TRUE, length(res)) else has_dep
    if (paired && sum(use) >= 2L) {
      n_rep <- ncol(ref$folds)
      R <- matrix(vapply(res[use], function(r) apply(r$oof[[tg]], 2, function(p) {
        ok <- is.finite(y) & is.finite(p)
        sqrt(mean((y[ok] - p[ok])^2))
      }), numeric(n_rep)), nrow = n_rep)
      best <- which.min(colMeans(R))
      tab[[paste0(tg, "_dRMSE_cv")]][use] <- colMeans(R - R[, best])
      tab[[paste0(tg, "_best_share")]][use] <- colMeans(R == apply(R, 1, min))
    }
  }
  tab$failed_fits <- vapply(res, function(r) as.integer(r$failures$n_failed), integer(1))
  rownames(tab) <- NULL
  
  out <- list(table = tab,
              checks = list(settings_differ = diff_set, same_folds = same_folds, same_observed = same_obs,
                            paired = paired, all_deployed = all(has_dep)),
              settings = list(k = ref$settings$k, reps = ref$settings$reps, loo = ref$settings$loo,
                              strata = ref$settings$strata, group = ref$settings$group, n_units = ref$settings$n_units))
  class(out) <- "mean_model_comparison"
  out
}

print.mean_model_comparison <- function(x, digits = 3, ...) {
  s <- x$settings
  cat("<mean_model_comparison>\n")
  cat(sprintf("%d models, %s, %s (%d units%s%s).\n", nrow(x$table),
              if (s$loo) "leave-one-out" else sprintf("%d-fold CV x %d repeat(s)", s$k, s$reps),
              if (x$checks$paired) "same folds and data" else "NOT on the same folds and data", s$n_units,
              if (!is.null(s$group)) sprintf(", units = '%s'", s$group) else "",
              if (!is.null(s$strata)) sprintf(", balanced within %s", paste(s$strata, collapse = "/")) else ""))
  tb <- x$table
  for (tg in if (any(tb$deployed)) c("model", "deployed") else "model") {
    cols <- grep(paste0("^", tg, "_"), names(tb), value = TRUE)
    part <- tb[, c("model", if (tg == "model") c("family", "n_coef"), cols)]
    names(part) <- sub(paste0("^", tg, "_"), "", names(part))
    cat(sprintf("\n%s\n", if (tg == "model") "(i) Model: each stimulus at its own predictor values" else
      "(ii) Deployed function (deployment constants)"))
    print(part, digits = digits, row.names = FALSE)
  }
  cat("\napp = apparent (in-sample); cv = CV mean over repeats; cv_sd = SD over repeats (fold randomness only);\n")
  cat("dRMSE_cv = paired CV RMSE difference from the best model; best_share = share of repeats with the lowest CV RMSE.\n")
  if (!x$checks$paired) cat("dRMSE_cv / best_share not computed: the runs do not share folds and data.\n")
  if (any(tb$deployed) && !all(tb$deployed)) cat("Results from plain models (deployed = FALSE) have no deployed values; (ii) compares the deployed runs only.\n")
  if (any(x$table$failed_fits > 0)) cat("NOTE: some refits failed (failed_fits); see each run's $failures.\n")
  invisible(x)
}