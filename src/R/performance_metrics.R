require(insight)
require(tidyverse)
require(performance)


# GEE Rsquared (Zheng, 2000) ---------------------------------------------------
geeRsquared <- function(...){
  
  model_objects <- list(...)
  
  # if (length(list(...)) == 1){
  #   
  #   model_objects <- insight::ellipsis_info(..., ..., only_models = TRUE,
  #                                           verbose=FALSE)
  #   model_objects <- model_objects[1]
  #   
  # } else{
  #   
  #   model_objects <- insight::ellipsis_info(..., only_models = TRUE,
  #                                           verbose=FALSE)
  #   
  # }
  
  # ensure proper object names
  dot_names <- sapply(match.call(expand.dots = FALSE)[["..."]], as.character)
  check_object_names <- insight::compact_character(names(model_objects))
  if ((is.null(check_object_names) ||
       # or if length of names doesn't match number of models
       length(check_object_names) != length(model_objects) ||
       # or if names are "..1", "..2" pattern
       all(grepl("\\.\\.\\d", check_object_names))) &&
      # and length of dot-names must match length of objects
      length(model_objects) == length(dot_names)) {
    names(model_objects) <- dot_names
  }
  
  object_names <- names(model_objects)
  
  m <- mapply(function(.x, .y) {
    SSE <- sum(.x$prior.weights*(.x$y - .x$fitted.values)^2, na.rm = T)
    SST <- sum(.x$prior.weights*(.x$y - mean(.x$y, na.rm = T))^2, na.rm = T)
    r2 = 1 - SSE/SST
    
    dat <- data.frame(R2=r2)
    model_name <- gsub("\"", "", insight::safe_deparse(.y), fixed = TRUE)
    perf_df <- data.frame(Name = model_name, Model = class(.x)[1], dat, stringsAsFactors = FALSE)
    perf_df
  }, model_objects, object_names, SIMPLIFY = FALSE)
  
  dfs <- Reduce(function(x, y) merge(x, y, all = TRUE, sort = FALSE), m)
  
  return(dfs)
}


# GEE entropy reduction (Zheng, 2000) ------------------------------------------

geeEntropy <- function(..., varest = c("robust", "df-adjusted", "model", "bias-corrected")) {
  
  varest <- match.arg(varest)
  
  model_objects <- list(...)
  
  # if (length(list(...)) == 1) {
  #   model_objects <- insight::ellipsis_info(..., ..., only_models = TRUE, verbose = FALSE)
  #   model_objects <- model_objects[1]
  # } else {
  #   model_objects <- insight::ellipsis_info(..., only_models = TRUE, verbose = FALSE)
  # }
  
  dot_names <- sapply(match.call(expand.dots = FALSE)[["..."]], as.character)
  check_object_names <- insight::compact_character(names(model_objects))
  if ((is.null(check_object_names) ||
       length(check_object_names) != length(model_objects) ||
       all(grepl("\\.\\.\\d", check_object_names))) &&
      length(model_objects) == length(dot_names)) {
    names(model_objects) <- dot_names
  }
  
  object_names <- names(model_objects)
  
  m <- mapply(function(.x, .y) {
    
    y  <- .x$y
    mu <- .x$fitted.values
    
    alpha <- mean(y, na.rm = TRUE)
    
    eps <- 1e-15
    mu    <- pmin(pmax(mu, eps), 1 - eps)
    alpha <- min(max(alpha, eps), 1 - eps)
    
    numerator <- sum(mu * log(mu) + (1 - mu) * log(1 - mu), na.rm = TRUE)
    
    N <- sum(!is.na(y))   # was length(!is.na(y)), which always equals length(y)
    denominator <- N * (alpha * log(alpha) + (1 - alpha) * log(1 - alpha))
    
    h <- 1 - numerator / denominator
    
    model_name <- gsub("\"", "", insight::safe_deparse(.y), fixed = TRUE)
    data.frame(Name = model_name, Model = class(.x)[1], H = h, stringsAsFactors = FALSE)
    
  }, model_objects, object_names, SIMPLIFY = FALSE)
  
  dfs <- Reduce(function(x, y) merge(x, y, all = TRUE, sort = FALSE), m)
  
  return(dfs)
}


# epsilon-insensitive RMSE for model predictions -------------------------------
epsilonRMSE <- function(measured_data, predicted_data, measured_target=NULL,
                        measured_ci_low=NULL, measured_ci_high=NULL, fitted='fit'){
  
  # assign column names
  if (is.null(measured_target)){
    measured_target <- names(measured_data)[1]
  }
  
  if (is.null(measured_ci_low)){
    measured_ci_low <- names(measured_data)[2]
  }
  
  if (is.null(measured_ci_high)){
    measured_ci_high <- names(measured_data)[3]
  }
  
  stopifnot(nrow(measured_data) == nrow(predicted_data))
  
  lower <- measured_data[[measured_ci_low]]
  upper <- measured_data[[measured_ci_high]]
  pred  <- predicted_data[[fitted]]
  
  stopifnot(all(lower <= upper, na.rm = TRUE))
  
  rmse <- sqrt(mean((pred - measured_data[[measured_target]])^2, na.rm = TRUE))
  
  # calculate residuals and set to 0 for predictions within the confidence interval
  epsilon_error <- pmax(0, pmax(lower - pred, pred - upper))
  
  n <- sum(!is.na(epsilon_error))
  
  # calculate epsilon-insensitive RMSE
  epsilon_rmse <- sqrt(mean(epsilon_error^2, na.rm = TRUE))
  
  return(list(epsilon_rmse=epsilon_rmse, rmse=rmse))
}


# (Herron) estimated percentage correct predictions (ePCP) for GLMM models -----
glmmPCP <- function(...){
  
  model_objects <- list(...)
  
  # if (length(list(...)) == 1){
  # 
  #   model_objects <- insight::ellipsis_info(..., ..., only_models = TRUE)
  #   model_objects <- model_objects[1]
  # 
  # } else{
  # 
  #   model_objects <- insight::ellipsis_info(..., only_models = TRUE)
  # 
  # }
  
  # ensure proper object names
  dot_names <- sapply(match.call(expand.dots = FALSE)[["..."]], as.character)
  check_object_names <- insight::compact_character(names(model_objects))
  if ((is.null(check_object_names) ||
       # or if length of names doesn't match number of models
       length(check_object_names) != length(model_objects) ||
       # or if names are "..1", "..2" pattern
       all(grepl("\\.\\.\\d", check_object_names))) &&
      # and length of dot-names must match length of objects
      length(model_objects) == length(dot_names)) {
    names(model_objects) <- dot_names
  }
  
  object_names <- names(model_objects)
  
  m <- mapply(function(.x, .y) {
    y_full <- .x@resp$y
    
    n_full <- suppressWarnings(insight::n_obs(.x))
    
    pr_full <- stats::predict(.x, type = "response")
    
    pcp_full <- (sum(1 - pr_full[y_full == 0]) + sum(pr_full[y_full == 1])) / n_full
    
    dat <- data.frame(ePCP=pcp_full)
    model_name <- gsub("\"", "", insight::safe_deparse(.y), fixed = TRUE)
    perf_df <- data.frame(Name = model_name, Model = class(.x)[1], dat, stringsAsFactors = FALSE)
    
  }, model_objects, object_names, SIMPLIFY = FALSE)
  
  dfs <- Reduce(function(x, y) merge(x, y, all = TRUE, sort = FALSE), m)
  
  return(dfs)
  
}

# Cohen's f2 effect size for GEE models ---------------------------------------

geeCohensf2 <- function(geeMod){
  # loop over each variable in the model, calculate the R2 using geeRsquared
  # function for a reduced formula excluding that variable, and calculate the f2
  # effect size
  
  # first, calculate the R2 for the full model
  R2full <- geeRsquared(geeMod)$rsquare_gee
  
  # # assign variables except intercept
  # vars <- geeMod$coefficients %>%
  #   rownames() %>%
  #   .[!grepl("Intercept", .)]
  # 
  f2 <- sapply(rownames(geeMod$coefficients), function(x) {
    # create a reduced formula
    reduced_formula <- update(geeMod$formula, paste0(". ~ . - ", x))
    # fit the reduced model
    reduced_model <- update(geeMod, formula = reduced_formula)
    # calculate the R2 for the reduced model
    R2reduced <- geeRsquared(reduced_model)$rsquare_gee
    # calculate the f2 effect size
    f2 <- (R2full - R2reduced) / (1 - R2full)
    
  })
  
  return(f2)
  
}


# McKelvey & Zavoina pseudo-R^2, dispatched across Bayesian -----------------
## bayes_R2_MZ(): ------------------------------------------------------------
# model engines (not restricted to brms)
#
# SCOPE CORRECTION (v0.4): v0.3 generalized this to glmgee (a frequentist/GEE
# engine) - wrong direction. What was actually wanted was support for OTHER
# BAYESIAN model-fitting packages (e.g. rstanarm, MCMCglmm), not frequentist
# ones. The glmgee branch has been removed entirely.
#
# DESIGN: brms remains the fully-supported engine (both marginal and
# conditional R^2, via brms::posterior_linpred(), which correctly handles
# group-level/random effects for either case). For OTHER Bayesian packages,
# this uses insight::get_parameters() - a generic accessor that returns
# posterior draws of the fixed-effect coefficients across many Bayesian
# engines (rstanarm, MCMCglmm, and others insight supports) - combined with
# the model's own design matrix, to compute MARGINAL linear-predictor draws
# generically.
#
# IMPORTANT LIMITATION, stated plainly rather than glossed over: this generic
# path only covers the MARGINAL case (fixed effects only, i.e. re_formula=NA
# in brms terms). Getting the CONDITIONAL linear predictor (including
# group-level/random effects) generically across arbitrary Bayesian packages
# would require also handling each package's own random-effects draw
# structure, which is not something insight's generic accessors uniformly
# expose in a way that can be safely combined into one formula here. For any
# non-brms model, requesting the conditional scope (re_formula = NULL, the
# default) will error explicitly rather than silently returning a wrong or
# partially-marginal number - you must pass re_formula = NA for a non-brms
# model, acknowledging you're getting the marginal R^2 only.
#
# THE disc/DISPERSION-SUBMODEL LOGIC REMAINS brms-ONLY, for a genuine reason,
# not an oversight: `dpar` auxiliary-parameter formulas (phi ~ ..., disc ~ ...)
# are a brms/Stan distributional-regression concept. Other Bayesian packages
# (rstanarm, MCMCglmm) do not share this exact mechanism, so there is no
# general "disc draws" to extract for them - they fall back to the constant
# base_var_res, which is the correct behaviour for a model that has no such
# heteroskedastic-precision structure in the first place.
#
# THE FOUR brms-SPECIFIC FIXES FROM v0.2 (marginal/conditional made explicit,
# disc-aware var_res with the Jensen's-inequality correction, robust family
# detection via insight::get_family(), and the brms::fitted() export issue)
# are all preserved unchanged for the brms path - see inline comments.
#
# NOT executed here - no R access in this environment. The brms path is the
# best-tested (built directly on this session's own diagnostics); the
# generic non-brms path is new and has not been run against any actual
# rstanarm/MCMCglmm model - please verify against one before trusting it. In
# particular, confirm directly the first time this runs on a real model that
# insight::find_formula(fit)$conditional is genuinely fixed-effects-only (no
# grouping/'|' syntax reaching model.matrix()) and that the name-based column
# alignment below actually engages as intended - e.g.:
#   insight::find_formula(fit)$conditional
#   colnames(insight::get_parameters(fit, effects = "fixed"))
#   colnames(model.matrix(insight::find_formula(fit)$conditional, data = insight::get_data(fit)))


### Engine dispatch: MARGINAL linear-predictor draws, one row per posterior draw -----

.dispatch_linpred_draws <- function(fit, re_formula = NULL, ...) {
  
  if (inherits(fit, "brmsfit")) {
    # brms natively supports both marginal (re_formula = NA) and conditional
    # (re_formula = NULL, the default) via posterior_linpred() - no
    # restriction needed here.
    return(brms::posterior_linpred(fit, transform = FALSE, re_formula = re_formula, ...))
  }
  
  # --- generic path for other Bayesian packages, via insight ---------------
  if (!identical(re_formula, NA)) {
    # i.e. re_formula was left at its default (NULL) rather than explicitly
    # set to NA - the caller is (implicitly or explicitly) asking for the
    # conditional scope, which this generic path cannot provide.
    stop(
      "Conditional R^2 (including group-level/random effects) is only ",
      "supported for brmsfit models in this function. For a '",
      paste(class(fit), collapse = "/"), "' model, call with ",
      "re_formula = NA explicitly to get the MARGINAL (fixed-effects-only) ",
      "R^2 instead - see the header comment for why this limitation exists."
    )
  }
  
  if (!requireNamespace("insight", quietly = TRUE)) {
    stop("Package 'insight' is required but not installed.")
  }
  
  param_draws <- as.matrix(insight::get_parameters(fit, effects = "fixed"))
  X <- stats::model.matrix(insight::find_formula(fit)$conditional,
                           data = insight::get_data(fit))
  
  # Align by NAME, not position - a matching column count is not sufficient
  # evidence the two matrices line up (different default contrast coding or
  # term-ordering conventions across packages could give same-width,
  # different-order matrices, which would silently multiply the wrong
  # coefficient against the wrong column with no error at all). Caught via
  # independent review.
  if (!all(colnames(X) %in% colnames(param_draws)) ||
      !all(colnames(param_draws) %in% colnames(X))) {
    stop(
      "Design matrix and parameter draws have different column names - ",
      "cannot safely align them.\n",
      "colnames(X): ", paste(colnames(X), collapse = ", "), "\n",
      "colnames(param_draws): ", paste(colnames(param_draws), collapse = ", ")
    )
  }
  X <- X[, colnames(param_draws), drop = FALSE]  # explicit alignment, not assumed
  
  X %*% t(param_draws) |> t()  # ndraws x N, matching brms's draws-matrix orientation
}


### Engine dispatch: dispersion/precision (disc-style) draws, if the model has one ----

.dispatch_disc_draws <- function(fit, re_formula = NULL, ...) {
  if (inherits(fit, "brmsfit")) {
    return(tryCatch(
      brms::posterior_epred(fit, dpar = "disc", re_formula = re_formula, summary = FALSE),
      error = function(e) NULL
    ))
  }
  NULL
}


### Main function --------------------------------------------------------------

bayes_R2_MZ <- function(fit, ci = 0.95, re_formula = NULL, ...) {
  
  if (!requireNamespace("insight", quietly = TRUE)) {
    stop("Package 'insight' is required but not installed.")
  }
  
  fam_info <- insight::get_family(fit)
  fam  <- fam_info$family
  link <- fam_info$link
  
  # --- fitted values on the latent/linear scale, honouring re_formula
  y_pred  <- .dispatch_linpred_draws(fit, re_formula = re_formula, ...)
  var_fit <- apply(y_pred, 1, stats::var)  # one value per draw
  
  binary_ordinal_families <- c("cumulative", "sratio", "cratio", "acat",
                               "bernoulli", "binomial", "categorical")
  
  if (fam %in% binary_ordinal_families) {
    
    base_var_res <- if (link %in% c("probit", "probit_approx")) 1 else pi^2 / 3
    
    # --- detect a fitted disc submodel (brms-only concept), rather than
    #     assuming disc = 1 
    disc_draws <- .dispatch_disc_draws(fit, re_formula = re_formula, ...)
    
    if (!is.null(disc_draws)) {
      # CORRECTNESS NOTE: the aggregation MUST be mean(1/disc_i^2), NOT
      # 1/mean(disc_i)^2 - these are NOT interchangeable when disc genuinely
      # varies across observations. f(x) = 1/x^2 is convex, so by Jensen's
      # inequality E[1/disc^2] >= 1/E[disc]^2 (equality only if disc is
      # constant). Using 1/mean(disc)^2 systematically UNDERESTIMATES
      # var_res and therefore INFLATES R2_MZ - caught via independent
      # review, verified by direct application of Jensen's inequality
      # before being adopted here.
      #
      # CAVEAT - still not an established/citable method: McKelvey &
      # Zavoina's measure assumes one constant residual variance; this
      # per-draw mean-of-per-observation-variances is a reasoned, but not
      # peer-reviewed, extension to the heteroskedastic case.
      var_res <- base_var_res * rowMeans(1 / disc_draws^2)
    } else {
      var_res <- base_var_res
    }
    
  } else {
    if (inherits(fit, "brmsfit")) {
      sigma_draws <- as.matrix(fit, variable = "sigma")
      var_res <- as.numeric(sigma_draws)^2
    } else {
      stop(
        "Gaussian-family residual variance extraction is only implemented for ",
        "brmsfit objects currently. Extend this branch if you need it for ",
        "another engine (check whether insight::get_parameters() exposes ",
        "sigma-equivalent draws for that package first)."
      )
    }
  }
  
  R2_MZ <- var_fit / (var_fit + var_res)
  
  tail_prob <- (1 - ci) / 2
  probs <- c(tail_prob, 1 - tail_prob)
  quantiles <- stats::quantile(R2_MZ, probs)
  q_names <- paste0("Q", probs * 100)
  
  out_df <- data.frame(
    Estimate = mean(R2_MZ),
    Est.Error = stats::sd(R2_MZ),
    q1 = quantiles[1],
    q2 = quantiles[2],
    row.names = "R2",
    check.names = FALSE
  )
  colnames(out_df)[3:4] <- q_names
  
  print(out_df, digits = 3)
  invisible(out_df)
}


### Example usage -------------------------------------------------------------

if (FALSE) {
  
  # --- brms models: unchanged, full support (marginal + conditional)
  bayes_R2_MZ(mA4)                    # conditional (default)
  bayes_R2_MZ(mA4, re_formula = NA)   # marginal
  bayes_R2_MZ(mA5)                    # disc-varying model - Jensen's-corrected
  
  # --- other Bayesian packages (NEW) - marginal only
  # e.g. a model fit with rstanarm::stan_glmer() or MCMCglmm::MCMCglmm():
  # bayes_R2_MZ(some_stanreg_fit, re_formula = NA)   # required for non-brms
  # bayes_R2_MZ(some_stanreg_fit)                    # errors - conditional
  #                                                   # not supported outside brms
}


# Bayesian adaptation of Lacy (2006)'s R2O for brms ordinal models ------------
## bayes_R2o_Lacy():  ---------------
# Lacy, M.G. (2006). "An Explained Variation Measure for Ordinal Response
# Models With Comparisons to Other Ordinal R2 Measures." Sociological Methods
# & Research, 34(4), 469-520. Stata module (r2o.ado) provided the reference
# algorithm this is translated from.
#
# WHY THIS IS THE RIGHT TOOL FOR A disc-VARYING (location-scale) ORDINAL
# MODEL, WHERE McKelvey & Zavoina's extension needed contested, hand-derived
# machinery: R2O is built ENTIRELY from the model's predicted CUMULATIVE
# CATEGORY PROBABILITIES - it never invokes a latent continuous variable or
# its variance at all. Whatever a disc/precision submodel does to an
# observation's predicted category probabilities is therefore automatically
# and correctly reflected in that observation's contribution to R2O, with no
# separate disc-aggregation formula needed (no Jensen's-inequality question,
# no mean-before-or-after-squaring debate - the issue simply doesn't arise).
#
# From Lacy's own abstract: R2O "was shown to outperform various pseudo-R-
# squared measures in estimating the value of the true R-squared for a
# regression model for an underlying continuous response, even though its
# sense does not require such [a latent variable]," and is "valid regardless
# of the method used to estimate the model."
#
# MECHANICS (translated from the Stata r2o.ado algorithm):
#   SY  (total/marginal variation) = sum_j[ F_j * (1 - F_j) ], where F_j is
#       the OBSERVED cumulative relative frequency of category <= j, computed
#       ONCE from the raw data (model-independent).
#   SYX (conditional/residual variation), per posterior draw = the MEAN,
#       across observations, of sum_j[ F_ij * (1 - F_ij) ], where F_ij is
#       that observation's PREDICTED cumulative probability of category <= j
#       under that draw.
#   R2O (per draw) = 1 - SYX / SY
#
# NOT IMPLEMENTED: Lacy's bias-adjusted "ur2o" divides by a degrees-of-
# freedom correction (N-1)/(N-numcovar-1) - this doesn't map cleanly onto a
# Bayesian multilevel model (no single well-defined "number of covariates
# used"), so it's deliberately omitted here rather than guessed at. If you
# need a bias-adjusted version, that mapping would need to be worked out and
# justified separately, not assumed.
#
# NOT executed here - no R access in this environment. In particular, verify
# directly that brms::posterior_epred() for your cumulative()/ordinal family
# returns a full ndraws x N x K array of per-category probabilities (rather
# than, say, an expected category index) before trusting this - this
# function assumes that shape.

bayes_R2o_Lacy <- function(fit, re_formula = NULL, ci = 0.95, ...) {
  
  if (!requireNamespace("insight", quietly = TRUE)) {
    stop("Package 'insight' is required but not installed.")
  }
  
  # --- observed response, integer-coded 1..K, matching the data actually
  #     used by the model (not a raw/unfiltered data frame - see this
  #     session's earlier NA-alignment lessons) 
  y_raw <- insight::get_response(fit)
  y <- as.integer(factor(y_raw, levels = sort(unique(y_raw))))
  K <- max(y)
  N <- length(y)
  
  # --- SY: marginal/total variation from the OBSERVED data (model-independent,
  #     computed once) 
  p_marg <- as.numeric(table(factor(y, levels = 1:K))) / N
  F_marg <- cumsum(p_marg)
  SY <- sum(F_marg * (1 - F_marg))
  
  if (SY <= 0) {
    stop("Observed marginal variation (SY) is zero or negative - check that ",
         "the response has more than one observed category.")
  }
  
  # --- predicted category probabilities, per draw: ndraws x N x K 
  # CRITICAL: summary = FALSE is required here. Without it, posterior_epred()
  # defaults to summary = TRUE and returns a SUMMARIZED array (point estimate/
  # SE/CI per observation-category combination), NOT per-draw probabilities -
  # dim(pred_probs)[1] would then silently be N (observations), not ndraws,
  # scrambling the entire per-draw loop below without any error being thrown.
  # This exact bug was present in an earlier version of this function and
  # produced plausible-looking but meaningless R2O values - caught only via
  # a downstream diagnostic returning an obviously-impossible number.
  pred_probs <- brms::posterior_epred(fit, re_formula = re_formula,
                                      summary = FALSE, ...)
  
  if (length(dim(pred_probs)) != 3) {
    stop(
      "Expected a 3-dimensional (ndraws x N x K) array from ",
      "posterior_epred() - got dimensions: ", paste(dim(pred_probs), collapse = " x "),
      ". This function assumes a categorical-probability output shape; ",
      "confirm this matches your family before proceeding."
    )
  }
  
  ndraws <- dim(pred_probs)[1]
  r2o_draws <- numeric(ndraws)
  
  for (d in seq_len(ndraws)) {
    P    <- pred_probs[d, , ]                # N x K predicted probabilities
    Fcum <- t(apply(P, 1, cumsum))            # N x K cumulative probabilities
    iSYX <- rowSums(Fcum * (1 - Fcum))        # per-observation conditional variation
    SYX  <- mean(iSYX)
    r2o_draws[d] <- 1 - SYX / SY
  }
  
  tail_prob <- (1 - ci) / 2
  probs <- c(tail_prob, 1 - tail_prob)
  quantiles <- stats::quantile(r2o_draws, probs)
  q_names <- paste0("Q", probs * 100)
  
  out_df <- data.frame(
    Estimate = mean(r2o_draws),
    Est.Error = stats::sd(r2o_draws),
    q1 = quantiles[1],
    q2 = quantiles[2],
    Scope = if (is.null(re_formula)) "conditional" else "marginal",
    row.names = "R2O",
    check.names = FALSE
  )
  colnames(out_df)[3:4] <- q_names
  
  print(out_df, digits = 3)
  invisible(out_df)
}


# Example usage

if (FALSE) {
  bayes_R2o_Lacy(mA5)                    # conditional (default)
  bayes_R2o_Lacy(mA5, re_formula = NA)   # marginal
  
  # Direct comparison against the M&Z figures already in your table -
  # given R2O needs no disc-heteroskedasticity correction at all, this
  # comparison is itself informative about whether M&Z's extension over-
  # or under-states the disc-driven difference:
  bayes_R2o_Lacy(mA1); bayes_R2o_Lacy(mA5)
}
