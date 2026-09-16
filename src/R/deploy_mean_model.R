# =============================================================================
# deploy_mean_model.R  (v13)
#
# A model-type-agnostic "deployment" toolkit that extracts the FIXED-EFFECTS
# MEAN STRUCTURE of a fitted model (glmtoolbox::glmgee, glmmTMB, or brms) and
# turns it into a standalone R function of the predictor variables.
#
# CHANGE LOG vs v12 (this version):
#   - Fixed a real bug (not a false alarm): the internal reduced-equation-
#     vs-model.matrix cross-check compared .eval_eta(spec, cc_data) — the
#     FULL model using the REAL values of every variable in the training
#     data — against the reduced-terms evaluator, which already has any
#     requested constants folded in. Whenever deploy_mean_model() was called
#     with a constant (e.g. TrialNumberScl = 0), these were two genuinely
#     different linear predictors for any row where the real value wasn't
#     the constant, so the check spuriously "disagreed" every time —
#     nothing to do with whether the deployed function's actual predictions
#     were correct (they weren't affected; predict_fn already applied
#     constants correctly). Fixed by substituting the same constants into
#     the data used for the "robust" side of the comparison before
#     evaluating it, so both sides now describe the same scenario.
#   - Confirmed (per your test) that a constant of exactly 0 for a variable
#     with no interactions correctly leaves the intercept unchanged, since
#     its contribution is coefficient * 0 = 0 — this was working correctly
#     already; the concern was understandable but not a bug.
#
# CHANGE LOG vs v11:
#   - Addressed a real robustness concern raised about v11: explaining away
#     brms's small mu-scale deviation as an "expected Jensen's-inequality
#     gap" was correct, but risked MASKING a genuine extraction bug behind
#     that same explanation, since both look like "a small nonzero number".
#     Fixed by removing the confound entirely rather than just explaining
#     it: the self-check now compares on the LINEAR PREDICTOR (eta) scale
#     for all three model types, not the response (mu) scale.
#       - glmgee: model$linear.predictors (confirmed from glmtoolbox source
#         to be exactly eta = X %*% beta_hat + offset) instead of
#         model$fitted.values.
#       - glmmTMB: stats::predict(..., type = "link") instead of
#         type = "response".
#       - brms: stats::fitted(..., scale = "linear") instead of the default
#         scale = "response" — confirmed from brms's posterior_epred.R
#         source (which we fetched and read) to return "results on the
#         scale of the linear predictor term, that is without applying the
#         inverse link function".
#     Since eta = X %*% beta is LINEAR in beta, E[eta] over posterior draws
#     is EXACTLY X %*% E[beta] for any model — no Jensen's-inequality gap is
#     possible on this scale, for any link function. A deviation here can
#     only mean the coefficient/design-matrix extraction itself is wrong
#     (for brms, what small residual remains reflects only Monte Carlo
#     sampling noise from a finite posterior draw count, not a structural
#     approximation) — a genuinely harder-to-fool check than v11's, per the
#     concern raised. The printed report now shows both max|dev| and a
#     scale-invariant relative deviation, with thresholds tightened
#     accordingly (0.1% for point-estimate models, 1% for brms, both much
#     stricter than v11's mu-scale thresholds since there's no more
#     legitimate gap to allow room for).
#
# CHANGE LOG vs v10:
#   - m3b's fitted-value self-check now runs and reports max|dev| = 0.00251,
#     RMSE = 0.00194 (a small but nonzero deviation). This is EXPECTED and
#     not a bug: brms's fitted(re_formula=NA) averages mu = g^-1(eta) across
#     posterior draws (i.e. reports E[g^-1(eta)]), while this deployment
#     evaluates g^-1 once at the posterior MEAN of the coefficients (i.e.
#     g^-1(E[eta])). For any nonlinear link these differ by a Jensen's-
#     inequality gap that grows with posterior uncertainty/correlation and
#     link curvature — unlike glmgee/glmmTMB, which are single-point-
#     estimate models with nothing to average over, hence their exact
#     (~1e-16) match. The 1e-4 alarm threshold was calibrated for those
#     point-estimate models and was misleadingly firing for a Bayesian
#     model's expected residual. Added a point_estimate_model flag (TRUE for
#     glmgee/glmmTMB, FALSE for brms) so the printed message and threshold
#     now correctly reflect which kind of gap is expected for each model
#     type, with an explanation of Jensen's inequality shown for brms models
#     and a much larger (0.05) threshold before flagging genuine concern.
#
# CHANGE LOG vs v9:
#   - Fixed "'fitted' is not an exported object from 'namespace:brms'",
#     surfaced by v9's error-capture change. Root cause: fitted() is a base-R
#     generic owned by 'stats'; brms only registers the S3 method
#     fitted.brmsfit against it and does not re-export a standalone symbol
#     named 'fitted' from its own namespace, so brms::fitted(...) fails even
#     though the method exists and dispatches fine. Changed to
#     stats::fitted(model, ...), which is the actual generic and correctly
#     dispatches to fitted.brmsfit regardless of which package the call is
#     qualified through, as long as brms is loaded. (fixef() was unaffected
#     by this because fixef is NOT a base-R generic — brms/glmmTMB each
#     define and export it themselves.)
#
# CHANGE LOG vs v8:
#   - m3b's equation now builds and simplifies correctly (deviation 0.0074
#     vs a 0.01 tolerance) — the name-translation fix worked. The one
#     remaining gap was "Model's fitted values could not be extracted" with
#     no further detail: brms::fitted(model, re_formula = NA) was failing
#     silently, since the error was being swallowed (tryCatch(..., error =
#     function(e) NULL)) rather than surfaced. Now the actual error message
#     is captured and reported in the fitted-value self-check's "reason"
#     text, and a fallback attempt with re_formula = ~0 (lme4-style syntax,
#     in case this fork/version prefers it over NA) is tried before giving
#     up. Please re-run and share the actual error text this surfaces —
#     I can't diagnose further without seeing it.
#
# CHANGE LOG vs v7:
#   - Fixed the 'Ilog10UASEvents' error on m3b. Root cause: brms sanitizes
#     population-level coefficient names for Stan compatibility, stripping
#     parentheses (e.g. "I(log10(UASEvents))" -> "Ilog10UASEvents"). Every
#     other part of this script assumes standard stats::model.matrix()
#     column naming, so these are now translated back to that standard form
#     using a synthetic reference frame (.make_synthetic_frame_for_colnames,
#     also now used generally inside .build_column_map, decoupling column-
#     name/term-association from raw-data availability everywhere, not just
#     for brms).
#   - Per brmsfit-class.R's own documentation ("@slot data A data.frame
#     containing all variables used in the model"), model$data is now tried
#     FIRST as the raw-data source for brms, ahead of the more fragile
#     call$data-evaluation trick (which remains as a fallback). This should
#     resolve the "Could not recover the original raw fitting data" warning
#     you saw even though m3Data/m3b$data was clearly present and correct.
#
# CHANGE LOG vs v6:
#   - Fixed the 'phi_Intercept' error you hit on the ordered-beta brms model
#     (m3b, which has a phi ~ AmbientEnv precision submodel). brms::fixef()
#     returns coefficients for ALL distributional parameters in one matrix,
#     prefixed by dpar name for anything other than the primary mu parameter
#     (phi_Intercept, phi_AmbientEnvStreet, etc.). Since this tool is
#     deliberately mean(mu)-only, those rows are now filtered out before
#     column_map construction, using the brms formula object's own list of
#     additional-dpar formulas (model$formula$pforms) where available, or a
#     heuristic prefix match as a fallback. This is working as intended, not
#     a flaw in your ordbeta model — the precision/phi submodel is real and
#     estimated, it's just out of scope for a MEAN-only deployment.
#   - Widened .recover_raw_data()'s environment search to try each candidate
#     environment in turn (including globalenv() explicitly) rather than
#     picking the first non-null one and giving up — this is what caused
#     raw-data recovery to fail for m3b even though m3Data clearly existed.
#
# CHANGE LOG vs v5:
#   - Simplified coefficients are now printed as whole-number ratios, e.g.
#     "(55/9)", instead of the decimal approximation (0.03279*log10(...) is
#     no longer misleading about the fact that it's actually a ratio search
#     result). See .format_coef_display() / .format_equation(). Unsimplified
#     coefficients still print as decimals, as before.
#   - Confirmed via brms source (posterior_epred.R): fitted.brmsfit() is a
#     thin wrapper around posterior_epred(..., scale="response"), which by
#     its own documentation returns draws of E[Y|X] = mu (the mean), with
#     SMALLER variance than posterior_predict() because posterior_predict
#     additionally includes residual/observation-level noise. This confirms
#     brms::fitted(model, re_formula = NA) (already used in .extract_spec_
#     brmsfit) is the right call — posterior_predict is not what a mean-
#     prediction self-check should use, and is not referenced anywhere in
#     this script.
#
# CHANGE LOG vs v4:
#   - Fixed the large (~0.38) fitted-value deviation seen for m1 (glmmTMB),
#     confirmed resolved by your test (max|dev| = 0 for m1). Root cause:
#     stats::fitted(model) for a model with random effects returns
#     CONDITIONAL fitted values (including each group's BLUPs); our
#     deployment is deliberately fixed-effects-only, so the comparison
#     target is now the population-level/marginal prediction via
#     stats::predict(model, re.form = ~0, type = "response").
#
# CHANGE LOG vs v3:
#   - Fixed the large (~0.38) fitted-value deviation seen for m1 (glmmTMB).
#     Root cause: stats::fitted(model) for a model with random effects
#     returns CONDITIONAL fitted values (including each group's estimated
#     random effects / BLUPs), but our deployment is deliberately fixed-
#     effects-only. Comparing against fitted() therefore always shows a
#     real discrepancy reflecting between-group variation, not a bug. Now
#     compares against the POPULATION-level/marginal prediction instead,
#     via stats::predict(model, re.form = ~0, type = "response") for
#     glmmTMB and brms::fitted(model, re_formula = NA) for brms — the
#     documented convention in both ecosystems for excluding all random/
#     group-level effects. m2 (glmgee, no random effects) was already
#     essentially exact (max|dev| ~ 2e-16), which is what confirmed the
#     core reconstruction logic itself was correct and the issue was
#     specific to what glmmTMB's default fitted() represents.
#
# CHANGE LOG vs v3:
#   - Fixed the na.omit-driven row-count mismatch: the fitted-value self-
#     check now uses the model's own stored fit frame (already post-na.omit,
#     already has I(...) evaluated) directly with stats::model.matrix(),
#     instead of trying to rebuild from raw data. See .eval_eta_from_frame()
#     and the updated .validate_against_fitted().
#   - The internal reduced-equation-vs-model.matrix cross-check, and the
#     default validation data used by coef_simplify, are now filtered to
#     complete cases first, so NA rows no longer turn tolerance comparisons
#     into NA (which previously also crashed the simplification loop).
#   - .simplify_reduced() now treats an NA deviation as a stop-and-warn
#     condition instead of erroring on if(NA).
#
# CHANGE LOG vs v2:
#   - FIXED: "object 'x' not found" for I(...) terms. Root cause: the model
#     frame stored on a fitted model (model$model / stats::model.frame(model))
#     has ALREADY had transforms like I(log10(UASEvents)) evaluated, and does
#     NOT retain the raw column ("UASEvents") the transform depended on. Using
#     that stored frame as input to a fresh stats::model.frame() call (to
#     re-derive predvars) therefore fails for any I(...) term. Fixed by
#     recovering the ORIGINAL raw data frame via the model's stored call
#     (model$call$data, evaluated in the formula's environment) instead of
#     using the post-transformation model frame. See .recover_raw_data().
#   - ADDED: automatic, unconditional validation against the model's own
#     fitted values (model$fitted.values for glmgee; stats::fitted() for
#     glmmTMB; brms::fitted() for brmsfit), run every time regardless of
#     coef_simplify. Reported via attr(f, "fitted_validation") and printed.
#   - lme4::nobars() now tried via the reformulas package first (nobars has
#     moved there as of a recent lme4 release, per the deprecation warning
#     you saw), falling back to lme4::nobars(), then a regex strip.
#
# Supported classes: glmgee (glmtoolbox), glmmTMB, brmsfit (brms)
# =============================================================================

# ---- Generic ----------------------------------------------------------

deploy_mean_model <- function(model, ...,
                              coef_simplify = FALSE,
                              tolerance = 0.01,
                              tolerance_type = c("relative", "absolute"),
                              max_denominator = 500,
                              validation_data = NULL,
                              verbose = TRUE) {
  UseMethod("deploy_mean_model")
}

deploy_mean_model.default <- function(model, ...) {
  stop(sprintf(
    "deploy_mean_model() does not support objects of class '%s'. Supported: glmgee, glmmTMB, brmsfit.",
    paste(class(model), collapse = "/")
  ))
}

deploy_mean_model.glmgee <- function(model, ...,
                                     coef_simplify = FALSE,
                                     tolerance = 0.01,
                                     tolerance_type = c("relative", "absolute"),
                                     max_denominator = 500,
                                     validation_data = NULL,
                                     verbose = TRUE) {
  spec <- .extract_spec_glmgee(model)
  .deploy_from_spec(spec, ...,
                    coef_simplify = coef_simplify, tolerance = tolerance,
                    tolerance_type = tolerance_type, max_denominator = max_denominator,
                    validation_data = validation_data, verbose = verbose)
}

deploy_mean_model.glmmTMB <- function(model, ...,
                                      coef_simplify = FALSE,
                                      tolerance = 0.01,
                                      tolerance_type = c("relative", "absolute"),
                                      max_denominator = 500,
                                      validation_data = NULL,
                                      verbose = TRUE) {
  spec <- .extract_spec_glmmTMB(model)
  .deploy_from_spec(spec, ...,
                    coef_simplify = coef_simplify, tolerance = tolerance,
                    tolerance_type = tolerance_type, max_denominator = max_denominator,
                    validation_data = validation_data, verbose = verbose)
}

deploy_mean_model.brmsfit <- function(model, ...,
                                      coef_simplify = FALSE,
                                      tolerance = 0.01,
                                      tolerance_type = c("relative", "absolute"),
                                      max_denominator = 500,
                                      validation_data = NULL,
                                      verbose = TRUE,
                                      point_estimate = c("mean", "median")) {
  point_estimate <- match.arg(point_estimate)
  spec <- .extract_spec_brmsfit(model, point_estimate = point_estimate)
  .deploy_from_spec(spec, ...,
                    coef_simplify = coef_simplify, tolerance = tolerance,
                    tolerance_type = tolerance_type, max_denominator = max_denominator,
                    validation_data = validation_data, verbose = verbose)
}

# ---- Raw-data recovery (shared across all three extraction methods) -------
# The model frame stored ON the fitted object has already evaluated any
# I(...) transforms and does not retain the raw columns they depended on
# (confirmed: this was the direct cause of the "object 'UASEvents' not
# found" error). We instead recover the ORIGINAL data frame by evaluating
# the model's stored call$data in the environment the formula was created
# in — a standard technique (used by e.g. car::Anova, update()), but one I
# have not been able to execute myself; if it fails, we fall back to the
# post-transformation model frame with an explicit warning about the
# consequences (see NOTES).
.recover_raw_data <- function(model, terms_obj) {
  cl <- tryCatch(stats::getCall(model), error = function(e) NULL)
  if (is.null(cl)) cl <- tryCatch(model[["call"]], error = function(e) NULL)
  if (is.null(cl) || is.null(cl[["data"]])) return(NULL)
  
  env_candidates <- list(
    tryCatch(attr(terms_obj, ".Environment"), error = function(e) NULL),
    tryCatch(environment(model[["formula"]]), error = function(e) NULL),
    tryCatch(environment(stats::formula(model)), error = function(e) NULL),
    globalenv(),
    parent.frame(2)
  )
  for (env in env_candidates) {
    if (is.null(env) || !is.environment(env)) next
    out <- tryCatch(eval(cl[["data"]], envir = env), error = function(e) NULL)
    if (is.data.frame(out)) return(out)
  }
  NULL
}

.data_has_required_vars <- function(data, terms_obj) {
  if (is.null(data) || !is.data.frame(data)) return(FALSE)
  needed <- all.vars(stats::formula(terms_obj))
  all(needed %in% names(data))
}

.nobars <- function(formula_obj) {
  if (requireNamespace("reformulas", quietly = TRUE)) {
    return(reformulas::nobars(formula_obj))
  }
  if (requireNamespace("lme4", quietly = TRUE)) {
    return(lme4::nobars(formula_obj))
  }
  warning("Neither the 'reformulas' nor 'lme4' package is available to strip random-effect ('(...|...)') ",
          "terms from the formula; falling back to a crude regex strip. Verify the resulting equation carefully.")
  .strip_bars_regex(formula_obj)
}

.strip_bars_regex <- function(formula_obj) {
  txt <- paste(deparse(formula_obj), collapse = " ")
  txt <- gsub("\\+?\\s*\\([^()]*\\|[^()]*\\)", "", txt)
  stats::as.formula(txt, env = environment(formula_obj))
}

# ---- Extraction: glmgee (glmtoolbox) ---------------------------------
# Verified 2026-07 against glmtoolbox source (github.com/cran/glmtoolbox/
# blob/master/R/geeglm.R). See NOTES.
.extract_spec_glmgee <- function(model) {
  cf_mat <- model[["coefficients"]]
  if (is.null(cf_mat)) stop("model$coefficients not found on this glmgee object.")
  cf <- stats::setNames(as.numeric(cf_mat), rownames(cf_mat))
  
  fam <- model[["family"]]
  if (is.null(fam) || is.null(fam$link)) stop("Could not find model$family$link on this glmgee object.")
  link <- fam$link
  
  terms_full <- model[["terms"]]
  if (is.null(terms_full)) stop("model$terms not found on this glmgee object.")
  terms_obj <- stats::delete.response(terms_full)
  
  raw_data <- .recover_raw_data(model, terms_obj)
  data_is_raw <- .data_has_required_vars(raw_data, terms_obj)
  if (!data_is_raw) {
    warning("Could not recover the original raw fitting data for this glmgee model (tried evaluating ",
            "model$call$data). Falling back to the stored model frame (model$model), which will fail for ",
            "any I(...)-transformed term when used as default validation data or in the automatic fitted-",
            "value self-check. Supply 'validation_data' explicitly (with raw, untransformed columns) for ",
            "coef_simplify, or pass raw data with matching column names via that argument to enable the ",
            "self-check.")
    raw_data <- model[["model"]]
  }
  
  xlevels <- model[["levels"]]
  contrasts_arg <- model[["contrasts"]]
  # Compare on the LINEAR PREDICTOR (eta) scale, not the response (mu) scale:
  # eta = X %*% beta is linear in beta, so there is no link-function-related
  # approximation gap to worry about here (unlike comparing on the mu scale
  # for a Bayesian model — see the brms extraction below). model$linear.
  # predictors is confirmed directly from glmtoolbox source (geeglm.R) to be
  # exactly eta = X %*% beta_hat + offset.
  fitted_vals <- tryCatch(as.numeric(model[["linear.predictors"]]), error = function(e) NULL)
  
  # The frame actually used in fitting (post-na.omit, with I(...) already
  # evaluated) — exactly row-aligned with fitted_vals, and usable directly
  # with stats::model.matrix() without needing raw columns. See NOTES.
  fit_frame <- model[["model"]]
  
  list(coefficients = cf, link = link, terms = terms_obj, data = raw_data,
       xlevels = xlevels, contrasts = contrasts_arg, fitted = fitted_vals,
       data_is_raw = data_is_raw, fit_frame = fit_frame, point_estimate_model = TRUE)
}

# ---- Extraction: glmmTMB -----------------------------------------------
.extract_spec_glmmTMB <- function(model) {
  cf <- glmmTMB::fixef(model)$cond
  fam <- stats::family(model)
  if (is.null(fam$link)) stop("Could not determine link from stats::family(model)$link.")
  link <- fam$link
  
  full_formula <- tryCatch(
    stats::formula(model, component = "cond"),
    error = function(e) tryCatch(stats::formula(model), error = function(e2) NULL)
  )
  if (is.null(full_formula)) stop("Could not retrieve the conditional-model formula from this glmmTMB object.")
  fixed_formula <- .nobars(full_formula)
  terms_obj <- stats::delete.response(stats::terms(fixed_formula))
  
  raw_data <- .recover_raw_data(model, terms_obj)
  data_is_raw <- .data_has_required_vars(raw_data, terms_obj)
  if (!data_is_raw) {
    warning("Could not recover the original raw fitting data for this glmmTMB model (tried evaluating ",
            "model$call$data). Falling back to the stored model frame, which will fail for any I(...)-",
            "transformed term when used as default validation data or in the automatic fitted-value ",
            "self-check. Supply 'validation_data' explicitly (with raw, untransformed columns) for ",
            "coef_simplify, or pass raw data via that argument to enable the self-check.")
    raw_data <- tryCatch(stats::model.frame(model), error = function(e) NULL)
    if (is.null(raw_data)) raw_data <- tryCatch(model[["frame"]], error = function(e) NULL)
  }
  if (is.null(raw_data)) stop("Could not recover any training data frame from this glmmTMB object.")
  
  # IMPORTANT: stats::fitted(model) for a mixed model returns CONDITIONAL
  # fitted values (i.e. including each group's estimated random effects/
  # BLUPs). Our deployable function is deliberately fixed-effects-only, so
  # the correct comparison target is the POPULATION-level / marginal
  # prediction with random effects excluded, obtained via re.form = ~0
  # (glmmTMB's documented convention, mirroring lme4). Compare on the
  # LINEAR PREDICTOR (eta) scale via type = "link", not type = "response":
  # eta = X %*% beta is linear in beta, so there's no link-function-related
  # approximation gap possible here, giving a cleaner, harder-to-fool check
  # than comparing on the mu scale — NOT verified by execution here.
  fitted_vals <- tryCatch(
    as.numeric(stats::predict(model, re.form = ~0, type = "link")),
    error = function(e) tryCatch(
      as.numeric(stats::predict(model, re.form = NA, type = "link")),
      error = function(e2) NULL
    )
  )
  
  fit_frame <- tryCatch(stats::model.frame(model), error = function(e) NULL)
  if (is.null(fit_frame)) fit_frame <- tryCatch(model[["frame"]], error = function(e) NULL)
  
  list(coefficients = cf, link = link, terms = terms_obj, data = raw_data,
       xlevels = NULL, contrasts = NULL, fitted = fitted_vals, data_is_raw = data_is_raw,
       fit_frame = fit_frame, point_estimate_model = TRUE)
}

# ---- Extraction: brmsfit (brms) ----------------------------------------
.extract_spec_brmsfit <- function(model, point_estimate = c("mean", "median")) {
  point_estimate <- match.arg(point_estimate)
  fx <- brms::fixef(model, robust = (point_estimate == "median"))
  cf <- fx[, "Estimate"]
  names(cf) <- rownames(fx)
  names(cf)[names(cf) == "Intercept"] <- "(Intercept)"
  
  link <- tryCatch(model[["family"]][["link"]], error = function(e) NULL)
  if (is.null(link)) {
    warning("Could not read model$family$link. Defaulting to 'logit' for the mean, consistent with the ",
            "published ordered-beta-regression mean-link convention (Kubinec; verify the exact citation ",
            "yourself). Confirm this is correct for your specific model before relying on it.")
    link <- "logit"
  }
  
  bf_obj <- tryCatch(stats::formula(model), error = function(e) NULL)
  mu_formula <- tryCatch(bf_obj[["formula"]], error = function(e) NULL)
  if (is.null(mu_formula)) mu_formula <- bf_obj
  if (is.null(mu_formula)) stop("Could not retrieve the mu-formula from this brmsfit object.")
  fixed_formula <- .nobars(mu_formula)
  terms_obj <- stats::delete.response(stats::terms(fixed_formula))
  
  # brms::fixef() returns coefficients for ALL distributional parameters in
  # one matrix, prefixed by dpar name for anything other than the primary
  # mu parameter (e.g. "phi_Intercept", "phi_AmbientEnvStreet" for an
  # ordered-beta model's precision submodel). Since this tool is
  # deliberately mean(mu)-only, those must be dropped before they ever
  # reach column_map/model.matrix construction. Preferred approach: read
  # the actual additional-dpar formula names directly off the brms formula
  # object (bf_obj$pforms, a named list keyed by dpar name) — NOT verified
  # by execution, so a heuristic prefix-based fallback is used if it's
  # unavailable or structured differently than expected. See NOTES.
  other_dpars <- tryCatch(names(bf_obj[["pforms"]]), error = function(e) NULL)
  if (is.null(other_dpars)) other_dpars <- character(0)
  if (length(other_dpars) == 0) {
    known_dpars <- c("sigma", "phi", "zi", "hu", "shape", "nu", "xi", "kappa",
                     "disc", "bs", "ndt", "bias", "alpha", "beta", "quantile")
    other_dpars <- known_dpars[vapply(known_dpars, function(d) {
      any(grepl(paste0("^", d, "_"), names(cf)))
    }, logical(1))]
    if (length(other_dpars) > 0) {
      warning("Could not read additional distributional-parameter names from the brms formula object ",
              "(model$formula$pforms); inferred them heuristically from coefficient name prefixes instead (",
              paste(other_dpars, collapse = ", "), "). Please double-check that only mean(mu)-related ",
              "coefficients remain in the deployed equation.")
    }
  }
  if (length(other_dpars) > 0) {
    drop_pattern <- paste0("^(", paste(other_dpars, collapse = "|"), ")_")
    cf <- cf[!grepl(drop_pattern, names(cf))]
  }
  
  # The @data slot is documented directly in brmsfit-class.R as "A
  # data.frame containing all variables used in the model" — i.e. the raw
  # variables the formula references (including e.g. UASEvents behind
  # I(log10(UASEvents))), NOT a post-transformation model frame like
  # glmgee/glmmTMB store. So try it FIRST, before the more fragile
  # call$data-evaluation trick, which is now just a fallback.
  raw_data <- tryCatch(model[["data"]], error = function(e) NULL)
  data_is_raw <- .data_has_required_vars(raw_data, terms_obj)
  if (!data_is_raw) {
    alt <- .recover_raw_data(model, terms_obj)
    if (.data_has_required_vars(alt, terms_obj)) {
      raw_data <- alt
      data_is_raw <- TRUE
    }
  }
  if (is.null(raw_data)) stop("Could not recover any training data frame (model$data) from this brmsfit object.")
  if (!data_is_raw) {
    warning("Neither model$data nor model$call$data (evaluated) contained all the raw variables this formula ",
            "needs. coef_simplify and the fitted-value self-check may not work reliably. Supply 'validation_data' ",
            "explicitly (with raw, untransformed columns) if needed.")
  }
  
  # brms sanitizes population-level coefficient names for Stan compatibility
  # — in particular, I(...) transforms lose their parentheses (e.g.
  # "I(log10(UASEvents))" becomes "Ilog10UASEvents"). Every other part of
  # this script assumes standard stats::model.matrix() column naming, so
  # translate brms's sanitized names back to that standard form here, using
  # a synthetic reference frame (structure/levels only — doesn't need real
  # raw data to work). NOT verified by execution — see NOTES.
  standard_cols <- tryCatch({
    synth <- .make_synthetic_frame_for_colnames(terms_obj, raw_data)
    colnames(stats::model.matrix(terms_obj, synth))
  }, error = function(e) NULL)
  if (!is.null(standard_cols)) {
    sanitized_lookup <- stats::setNames(standard_cols, gsub("[()]", "", standard_cols))
    matched <- names(cf) %in% names(sanitized_lookup)
    names(cf)[matched] <- sanitized_lookup[names(cf)[matched]]
    unmatched <- setdiff(names(cf)[!matched], "(Intercept)")
    if (length(unmatched) > 0) {
      warning("Could not translate brms coefficient name(s) to standard model.matrix column names: ",
              paste(unmatched, collapse = ", "), ". The printed equation / coefficient simplification may fail ",
              "for these terms even though predictions from the deployed function remain unaffected.")
    }
  } else {
    warning("Could not build a reference model.matrix to translate brms's sanitized coefficient names; the ",
            "printed equation and coefficient simplification will likely be unavailable for this model.")
  }
  
  # NOTE: fitted is a base-R generic owned by 'stats', not something brms
  # itself exports as a standalone symbol (brms only registers the S3
  # method fitted.brmsfit against that generic) — calling brms::fitted(...)
  # therefore fails with "'fitted' is not an exported object from
  # 'namespace:brms'" even though the method exists and dispatches fine via
  # stats::fitted(...). Confirmed by your v9 test run's captured error
  # message. fixef() worked directly via brms:: earlier because fixef is
  # NOT a base-R generic — brms has to define and export that one itself.
  #
  # scale = "linear" is critical here: per brms's own posterior_epred.R
  # source (which we fetched and read directly), scale="linear" returns
  # "results on the scale of the linear predictor term, that is without
  # applying the inverse link function" — i.e. eta, not mu. Comparing on
  # the mu scale confounds two different things (our reconstruction's
  # g^-1(E[eta]) vs brms's E[g^-1(eta)] averaged over posterior draws),
  # which differ by a real but harmless Jensen's-inequality gap for any
  # nonlinear link. Comparing on the eta scale instead removes that
  # confound entirely: eta = X %*% beta is LINEAR in beta, so E[eta] over
  # posterior draws is EXACTLY X %*% E[beta] — no approximation gap, only
  # ordinary Monte Carlo sampling noise from a finite number of draws. A
  # deviation here can only mean the coefficient/design-matrix extraction
  # itself is wrong, which is a much more trustworthy signal.
  fitted_error <- NULL
  fitted_vals <- tryCatch({
    fv <- stats::fitted(model, re_formula = NA, scale = "linear", robust = (point_estimate == "median"))
    as.numeric(fv[, "Estimate"])
  }, error = function(e) {
    fitted_error <<- conditionMessage(e)
    tryCatch({
      fv2 <- stats::fitted(model, re_formula = ~0, scale = "linear", robust = (point_estimate == "median"))
      as.numeric(fv2[, "Estimate"])
    }, error = function(e2) {
      fitted_error <<- paste0(fitted_error, " | re_formula=~0 also failed: ", conditionMessage(e2))
      NULL
    })
  })
  
  fit_frame <- tryCatch(model[["data"]], error = function(e) NULL)
  
  list(coefficients = cf, link = link, terms = terms_obj, data = raw_data,
       xlevels = NULL, contrasts = NULL, fitted = fitted_vals, data_is_raw = data_is_raw,
       fit_frame = fit_frame, fitted_error = fitted_error, point_estimate_model = FALSE)
}

# ---- Robust numeric evaluator (authoritative for predictions) -----------
.eval_eta <- function(spec, data) {
  mf <- stats::model.frame(spec$terms, data, xlev = spec$xlevels, na.action = stats::na.pass)
  mm <- stats::model.matrix(spec$terms, mf, contrasts.arg = spec$contrasts)
  
  missing_cols <- setdiff(names(spec$coefficients), colnames(mm))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Model-matrix columns built from the supplied data do not include: %s.\nColumns present: %s\nThis usually means a factor is missing a level, or non-default contrasts were used.",
      paste(missing_cols, collapse = ", "), paste(colnames(mm), collapse = ", ")))
  }
  mm <- mm[, names(spec$coefficients), drop = FALSE]
  as.numeric(mm %*% spec$coefficients)
}


# Evaluates eta directly from an ALREADY-BUILT model frame (e.g. the frame
# a package used internally when fitting), skipping stats::model.frame()'s
# re-derivation of predvars. This is exactly what predict.glmgee-style code
# does when the frame is already in the right shape, and — critically — it
# does NOT need the raw columns behind any I(...) transform, since those are
# already evaluated in the stored frame. Use this only with a frame that
# genuinely came from the model (row-aligned with fitted values); for
# arbitrary new/raw data use .eval_eta() instead.
.eval_eta_from_frame <- function(spec, frame) {
  mm <- stats::model.matrix(spec$terms, frame, contrasts.arg = spec$contrasts)
  missing_cols <- setdiff(names(spec$coefficients), colnames(mm))
  if (length(missing_cols) > 0) {
    stop(sprintf("Model-matrix columns from the stored fit frame do not include: %s.",
                 paste(missing_cols, collapse = ", ")))
  }
  mm <- mm[, names(spec$coefficients), drop = FALSE]
  as.numeric(mm %*% spec$coefficients)
}

.complete_cases_for <- function(data, vars) {
  vars <- intersect(vars, names(data))
  cc <- stats::complete.cases(data[, vars, drop = FALSE])
  data[cc, , drop = FALSE]
}

.get_inverse_link <- function(link) {
  switch(link,
         logit    = stats::plogis,
         log      = exp,
         identity = identity,
         probit   = stats::pnorm,
         cloglog  = function(eta) 1 - exp(-exp(eta)),
         stop(sprintf("Inverse link for '%s' is not implemented. Add it to .get_inverse_link().", link))
  )
}

# ---- Synthetic reference frame (for column-name/assign lookups only) ------
# stats::model.matrix()'s column names and "assign" attribute depend only on
# the TERMS STRUCTURE and factor LEVELS — not on the actual row values. So
# for purposes of naming/term-association (never for real prediction), we
# can build a minimal 1-row frame: real factor columns get a value drawn
# from their actual levels (with the full level set attached, so all dummy
# columns still appear), and every other variable (including raw variables
# an I(...) transform depends on, which may be genuinely unavailable — see
# NOTES) gets an arbitrary numeric placeholder. This decouples column-name/
# term-association logic from whether raw training data was recoverable.
.make_synthetic_frame_for_colnames <- function(terms_obj, data) {
  vars <- all.vars(stats::formula(terms_obj))
  cols <- lapply(vars, function(v) {
    if (!is.null(data) && v %in% names(data) && (is.factor(data[[v]]) || is.character(data[[v]]))) {
      lv <- if (is.factor(data[[v]])) levels(data[[v]]) else sort(unique(as.character(data[[v]])))
      if (length(lv) == 0) lv <- "L1"
      factor(lv[1], levels = lv)
    } else {
      1
    }
  })
  names(cols) <- vars
  as.data.frame(cols, stringsAsFactors = FALSE)
}

# ---- Term/column parsing for the printable "reduced" equation -----------
# (unchanged from v2 — see prior version's comments; uses terms()'s own
# factors/variables/term.labels attributes and model.matrix()'s "assign"
# attribute, which are documented/official outputs, not guesswork. Assumes
# default treatment contrasts.)
.build_column_map <- function(terms_obj, coef_names, data) {
  term_labels <- attr(terms_obj, "term.labels")
  if (length(term_labels) == 0) return(stats::setNames(list(), character(0)))
  
  var_exprs <- attr(terms_obj, "variables")[-1]
  var_expr_strs <- vapply(var_exprs, function(e) paste(deparse(e), collapse = ""), character(1))
  factors_tab <- attr(terms_obj, "factors")
  
  factor_like <- stats::setNames(
    vapply(var_expr_strs, function(ve) ve %in% names(data) && (is.factor(data[[ve]]) || is.character(data[[ve]])), logical(1)),
    var_expr_strs
  )
  factor_levels <- stats::setNames(
    lapply(var_expr_strs[factor_like], function(v) {
      col <- data[[v]]
      if (is.factor(col)) levels(col) else sort(unique(as.character(col)))
    }),
    var_expr_strs[factor_like]
  )
  
  term_pieces_list <- stats::setNames(vector("list", length(term_labels)), term_labels)
  for (tl in term_labels) {
    idx <- which(factors_tab[, tl] > 0)
    piece_names <- rownames(factors_tab)[idx]
    term_pieces_list[[tl]] <- lapply(piece_names, function(pn) {
      if (isTRUE(factor_like[[pn]])) {
        list(type = "factor", raw_var = pn, expr = NA_character_)
      } else {
        rv <- all.vars(str2lang(pn))
        if (length(rv) != 1) {
          stop(sprintf(
            "Term '%s' depends on %d raw variables (%s); this function currently only supports transformations of a single raw variable (e.g. I(log10(x)), I(x^2), I(exp(x))). Multi-variable I() expressions are not supported.",
            pn, length(rv), paste(rv, collapse = ", ")))
        }
        list(type = "continuous", raw_var = rv, expr = pn)
      }
    })
  }
  
  ref_mm <- stats::model.matrix(terms_obj, .make_synthetic_frame_for_colnames(terms_obj, data))
  assign_vec <- attr(ref_mm, "assign")
  mm_colnames <- colnames(ref_mm)
  colname_to_term <- stats::setNames(vector("list", length(mm_colnames)), mm_colnames)
  for (i in seq_along(mm_colnames)) {
    if (assign_vec[i] == 0) next
    colname_to_term[[mm_colnames[i]]] <- term_labels[assign_vec[i]]
  }
  
  col_names <- setdiff(coef_names, "(Intercept)")
  column_map <- stats::setNames(vector("list", length(col_names)), col_names)
  for (cn in col_names) {
    tl <- colname_to_term[[cn]]
    if (is.null(tl)) {
      stop(sprintf(
        "Could not associate coefficient '%s' with a model term via model.matrix()'s assign attribute. Inspect colnames(stats::model.matrix(terms(model), data)) manually.",
        cn))
    }
    pieces_template <- term_pieces_list[[tl]]
    subtokens <- strsplit(cn, ":", fixed = TRUE)[[1]]
    resolved <- vector("list", length(pieces_template))
    for (i in seq_along(pieces_template)) {
      pt <- pieces_template[[i]]
      if (pt$type == "continuous") {
        resolved[[i]] <- pt
      } else {
        cand <- which(vapply(subtokens, function(s) {
          startsWith(s, pt$raw_var) && substring(s, nchar(pt$raw_var) + 1) %in% factor_levels[[pt$raw_var]]
        }, logical(1)))
        if (length(cand) != 1) {
          stop(sprintf("Could not match a factor level for '%s' within coefficient column '%s' (assumes default treatment contrasts).", pt$raw_var, cn))
        }
        lvl <- substring(subtokens[cand], nchar(pt$raw_var) + 1)
        resolved[[i]] <- list(type = "factor", raw_var = pt$raw_var, level = lvl, expr = NA_character_)
      }
    }
    column_map[[cn]] <- resolved
  }
  column_map
}

.reduce_columns <- function(coefficients, column_map, constants, data) {
  beta0 <- unname(coefficients[["(Intercept)"]])
  reduced <- list()
  
  for (cn in names(column_map)) {
    pieces <- column_map[[cn]]
    coef_val <- unname(coefficients[[cn]])
    
    fixed_pieces <- Filter(function(p) p$raw_var %in% names(constants), pieces)
    free_pieces  <- Filter(function(p) !(p$raw_var %in% names(constants)), pieces)
    
    multiplier <- 1
    for (p in fixed_pieces) {
      const_val <- constants[[p$raw_var]]
      if (p$type == "continuous") {
        env <- list2env(stats::setNames(list(const_val), p$raw_var), parent = baseenv())
        val <- eval(str2lang(p$expr), envir = env)
      } else {
        val <- as.numeric(identical(as.character(const_val), p$level))
      }
      multiplier <- multiplier * val
      if (multiplier == 0) break
    }
    
    new_coef <- coef_val * multiplier
    if (new_coef == 0) next
    
    if (length(free_pieces) == 0) {
      beta0 <- beta0 + new_coef
    } else {
      key <- paste(sort(vapply(free_pieces, function(p)
        if (p$type == "continuous") p$expr else paste0(p$raw_var, "==", p$level),
        character(1))), collapse = " : ")
      if (is.null(reduced[[key]])) {
        reduced[[key]] <- list(coef = new_coef, pieces = free_pieces)
      } else {
        reduced[[key]]$coef <- reduced[[key]]$coef + new_coef
      }
    }
  }
  list(beta0 = beta0, terms = reduced)
}

.eval_reduced_mu <- function(beta0, reduced_terms, data, link) {
  inv_link <- .get_inverse_link(link)
  n <- nrow(data)
  eta <- rep(beta0, n)
  for (term in reduced_terms) {
    contrib <- rep(1, n)
    for (p in term$pieces) {
      piece_val <- if (p$type == "continuous") {
        eval(str2lang(p$expr), envir = data)
      } else {
        as.numeric(as.character(data[[p$raw_var]]) == p$level)
      }
      contrib <- contrib * piece_val
    }
    eta <- eta + term$coef * contrib
  }
  inv_link(eta)
}

# ---- Coefficient simplification (rational approximation) -----------------
.rational_approx <- function(x, max_denominator) {
  if (x == 0) return(list(num = 0, den = 1, value = 0))
  s <- sign(x); x <- abs(x)
  h0 <- 0; h1 <- 1; k0 <- 1; k1 <- 0
  b <- x
  num <- round(x); den <- 1
  repeat {
    a <- floor(b)
    h2 <- a * h1 + h0
    k2 <- a * k1 + k0
    if (k2 > max_denominator || !is.finite(k2)) break
    num <- h2; den <- k2
    h0 <- h1; h1 <- h2
    k0 <- k1; k1 <- k2
    if (b == a) break
    b <- 1 / (b - a)
    if (!is.finite(b)) break
  }
  list(num = s * num, den = den, value = s * num / den)
}

.simplify_reduced <- function(beta0, reduced_terms, data, link, original_mu,
                              tolerance, tolerance_type, max_denominator) {
  denom_bound <- 1
  deviation <- Inf
  approx_beta0 <- NULL
  approx_terms <- NULL
  
  repeat {
    denom_bound <- denom_bound + 1
    approx_beta0 <- .rational_approx(beta0, denom_bound)
    approx_terms <- lapply(reduced_terms, function(t) {
      ap <- .rational_approx(t$coef, denom_bound)
      list(coef = ap$value, pieces = t$pieces, num = ap$num, den = ap$den)
    })
    approx_mu <- .eval_reduced_mu(approx_beta0$value, approx_terms, data, link)
    
    deviation <- if (tolerance_type == "relative") {
      max(abs((approx_mu - original_mu) / pmax(abs(original_mu), 1e-8)))
    } else {
      max(abs(approx_mu - original_mu))
    }
    if (is.na(deviation)) {
      warning("Deviation could not be computed (NA), likely due to missing values in the validation data. ",
              "Stopping simplification search at the current denominator bound; treat the result with caution.")
      break
    }
    if (deviation <= tolerance || denom_bound >= max_denominator) break
  }
  
  if (is.na(deviation) || deviation > tolerance) {
    warning(sprintf("Simplification did not reach tolerance %.4g within max_denominator = %d; achieved deviation = %.4g.",
                    tolerance, max_denominator, deviation))
  }
  list(beta0 = approx_beta0, terms = approx_terms, deviation = deviation, denom_bound = denom_bound)
}

# ---- Validation against the model's own fitted values ---------------------
# Runs unconditionally (regardless of coef_simplify): reconstructs eta (the
# LINEAR PREDICTOR, before any link function) on the raw training data using
# ONLY spec$coefficients/spec$terms (i.e. ignoring any constants requested in
# THIS deploy_mean_model() call) and compares against the model's own
# in-sample linear predictor. Comparing on the eta scale rather than the mu
# (response) scale is deliberate: eta = X %*% beta is LINEAR in beta, so
# there is no link-function-related approximation gap possible for ANY
# model type — a deviation here can only mean the coefficient/design-matrix
# extraction itself is wrong, not an artifact of a nonlinear link or of
# averaging over a posterior distribution. See the v12 changelog note.
.validate_against_fitted <- function(spec) {
  if (is.null(spec$fitted)) {
    reason <- "Model's linear predictor could not be extracted."
    if (!is.null(spec$fitted_error)) {
      reason <- paste0(reason, " Underlying error: ", spec$fitted_error)
    }
    return(list(ok = FALSE, reason = reason))
  }
  if (is.null(spec$fit_frame)) {
    return(list(ok = FALSE, reason = "The model's internal fit frame could not be recovered; self-check skipped."))
  }
  eta_hat <- tryCatch(.eval_eta_from_frame(spec, spec$fit_frame), error = function(e) e)
  if (inherits(eta_hat, "error")) {
    return(list(ok = FALSE, reason = sprintf("Reconstruction from the stored fit frame failed: %s", conditionMessage(eta_hat))))
  }
  if (length(eta_hat) != length(spec$fitted)) {
    return(list(ok = FALSE, reason = sprintf(
      "Row-count mismatch even using the model's own stored fit frame (reconstructed: %d, model's eta: %d) — this is unexpected and worth investigating; self-check skipped.",
      length(eta_hat), length(spec$fitted))))
  }
  dev <- eta_hat - spec$fitted
  max_abs_dev <- max(abs(dev))
  rel_dev <- max(abs(dev) / pmax(abs(spec$fitted), 1e-8))
  list(ok = TRUE, max_abs_dev = max_abs_dev, rel_dev = rel_dev,
       rmse = sqrt(mean(dev^2)), n = length(dev),
       point_estimate_model = isTRUE(spec$point_estimate_model))
}


# ---- Equation formatting -------------------------------------------------
.format_piece_display <- function(p) {
  if (p$type == "continuous") {
    expr <- p$expr
    if (grepl("^I\\(.*\\)$", expr)) expr <- sub("^I\\((.*)\\)$", "\\1", expr)
    expr
  } else {
    sprintf('[%s=="%s"]', p$raw_var, p$level)
  }
}

# Displays a coefficient as a whole-number ratio "(num/den)" when num/den are
# available (i.e. this coefficient came from .rational_approx()), otherwise
# falls back to a decimal. den is always a positive integer by construction
# (.rational_approx folds sign into num); den == 1 collapses to a bare
# integer rather than "(n/1)".
.format_coef_display <- function(coef_val, num = NULL, den = NULL, digits = 4) {
  if (!is.null(num) && !is.null(den)) {
    num_i <- as.integer(round(num)); den_i <- as.integer(round(den))
    if (den_i == 1) as.character(num_i) else sprintf("(%d/%d)", num_i, den_i)
  } else {
    as.character(signif(coef_val, digits))
  }
}

.format_equation <- function(beta0, reduced_terms, link, digits = 4,
                             beta0_num = NULL, beta0_den = NULL) {
  term_strs <- vapply(reduced_terms, function(t) {
    piece_strs <- vapply(t$pieces, .format_piece_display, character(1))
    coef_str <- .format_coef_display(t$coef, t[["num"]], t[["den"]], digits)
    paste0(coef_str, "*", paste(piece_strs, collapse = "*"))
  }, character(1))
  
  beta0_str <- .format_coef_display(beta0, beta0_num, beta0_den, digits)
  eta_terms <- c(beta0_str, term_strs)
  linear_predictor <- paste(eta_terms, collapse = " + ")
  linear_predictor <- gsub("\\+ -", "- ", linear_predictor)
  eta_str <- paste0("(", linear_predictor, ")")
  
  rhs <- switch(link,
                logit    = sprintf("1 / (1 + exp(-%s))", eta_str),
                log      = sprintf("exp(%s)", eta_str),
                identity = eta_str,
                probit   = sprintf("Phi(%s)", eta_str),
                cloglog  = sprintf("1 - exp(-exp(%s))", eta_str),
                sprintf("g^-1(%s)  [link = %s]", eta_str, link)
  )
  sprintf("mu = %s", rhs)
}

# ---- Orchestrator ---------------------------------------------------------
.deploy_from_spec <- function(spec, ..., coef_simplify, tolerance, tolerance_type,
                              max_denominator, validation_data, verbose) {
  constants <- list(...)
  tolerance_type <- match.arg(tolerance_type, c("relative", "absolute"))
  
  cf <- spec$coefficients
  if (!"(Intercept)" %in% names(cf)) {
    stop("No intercept term found among model coefficients; models fitted without an intercept are not currently supported.")
  }
  
  all_raw_vars <- all.vars(stats::formula(spec$terms))
  unknown_constants <- setdiff(names(constants), all_raw_vars)
  if (length(unknown_constants) > 0) {
    stop(sprintf("Variable(s) not found among model fixed-effect terms: %s", paste(unknown_constants, collapse = ", ")))
  }
  free_vars <- setdiff(all_raw_vars, names(constants))
  
  # ---- Unconditional self-check against the model's own fitted values ----
  fitted_check <- .validate_against_fitted(spec)
  
  # ---- Robust prediction function (authoritative) ----
  force(spec); force(constants); force(free_vars)
  inv_link <- .get_inverse_link(spec$link)
  predict_fn <- function(newdata) {
    stopifnot(is.data.frame(newdata))
    missing_vars <- setdiff(free_vars, names(newdata))
    if (length(missing_vars) > 0) {
      stop(sprintf("newdata is missing required column(s): %s", paste(missing_vars, collapse = ", ")))
    }
    full_data <- newdata
    for (v in names(constants)) full_data[[v]] <- constants[[v]]
    eta <- .eval_eta(spec, full_data)
    inv_link(eta)
  }
  
  # ---- Reduced-term representation (for printing / simplification) ----
  reduced_ok <- TRUE
  reduced_terms <- list()
  beta0 <- unname(cf[["(Intercept)"]])
  tryCatch({
    column_map <- .build_column_map(spec$terms, names(cf), spec$data)
    red <- .reduce_columns(cf, column_map, constants, spec$data)
    beta0 <- red$beta0
    reduced_terms <- red$terms
  }, error = function(e) {
    reduced_ok <<- FALSE
    warning(sprintf(
      "Could not build a merged/reduced printable equation (%s). Predictions from the returned function are unaffected (they use stats::model.matrix() directly), but the printed equation and coefficient simplification are unavailable.",
      conditionMessage(e)))
  })
  
  simp <- NULL
  if (reduced_ok) {
    if (!isTRUE(spec$data_is_raw) && is.null(validation_data)) {
      warning("Skipping the internal reduced-equation-vs-model.matrix cross-check: raw training data with the ",
              "original (untransformed) columns is unavailable for this model (see extraction warning above). ",
              "The returned function's predictions still use the robust model.matrix()-based path; treat the ",
              "printed equation with a little extra caution until it can be checked, e.g. by supplying ",
              "'validation_data' with raw columns.")
    } else {
      cc_data <- if (!is.null(validation_data)) validation_data else spec$data
      cc_data <- .complete_cases_for(cc_data, all_raw_vars)
      # The reduced/printable equation (beta0, reduced_terms) already has
      # any requested constants folded in. For a fair comparison, the
      # "robust" side must evaluate the SAME (constants-substituted)
      # scenario, not the raw training data's real values for those
      # variables — otherwise this check compares two different models
      # whenever constants are supplied, and will spuriously "disagree".
      cc_data_with_constants <- cc_data
      for (v in names(constants)) cc_data_with_constants[[v]] <- constants[[v]]
      check_ok <- tryCatch({
        check_mu_robust  <- inv_link(.eval_eta(spec, cc_data_with_constants))
        check_mu_reduced <- .eval_reduced_mu(beta0, reduced_terms, cc_data, spec$link)
        max(abs(check_mu_robust - check_mu_reduced)) <= 1e-6
      }, error = function(e) NA)
      if (!isTRUE(check_ok)) {
        warning("Could not confirm that the reduced/printable equation numerically matches the ",
                "stats::model.matrix()-based predictions on (complete cases of) the training data (either the ",
                "check failed to run, or the two disagreed by more than 1e-6). The returned function's ",
                "predictions always use the model.matrix()-based path and should still be correct, but treat ",
                "the printed equation and any coefficient simplification with caution until you've verified ",
                "deployed(newdata) against predict() from the original model.")
      }
    }
    
    if (coef_simplify) {
      if (!isTRUE(spec$data_is_raw) && is.null(validation_data)) {
        stop("coef_simplify = TRUE requires raw training data with the original (untransformed) columns, ",
             "which could not be recovered automatically for this model (see warning above). Supply it via ",
             "the 'validation_data' argument.")
      }
      val_data <- if (!is.null(validation_data)) validation_data else spec$data
      required_vars <- unique(unlist(lapply(reduced_terms, function(t) vapply(t$pieces, `[[`, character(1), "raw_var"))))
      missing_cols <- setdiff(required_vars, names(val_data))
      if (length(missing_cols) > 0) {
        stop(sprintf("Validation data is missing required column(s): %s", paste(missing_cols, collapse = ", ")))
      }
      val_data <- .complete_cases_for(val_data, required_vars)
      full_val <- val_data
      for (v in names(constants)) full_val[[v]] <- constants[[v]]
      original_mu <- inv_link(.eval_eta(spec, full_val))
      simp <- .simplify_reduced(beta0, reduced_terms, val_data, spec$link, original_mu,
                                tolerance, tolerance_type, max_denominator)
    }
  } else if (coef_simplify) {
    stop("coef_simplify = TRUE requires the reduced/printable equation, which could not be built for this model (see warning above).")
  }
  
  f <- predict_fn
  attr(f, "link") <- spec$link
  attr(f, "constants") <- constants
  attr(f, "free_vars") <- free_vars
  attr(f, "fitted_validation") <- fitted_check
  if (reduced_ok) {
    disp_beta0 <- if (!is.null(simp)) simp$beta0$value else beta0
    disp_terms <- if (!is.null(simp)) simp$terms else reduced_terms
    beta0_num <- if (!is.null(simp)) simp$beta0$num else NULL
    beta0_den <- if (!is.null(simp)) simp$beta0$den else NULL
    attr(f, "equation") <- .format_equation(disp_beta0, disp_terms, spec$link,
                                            beta0_num = beta0_num, beta0_den = beta0_den)
    attr(f, "coefficients") <- list(intercept = disp_beta0, terms = disp_terms)
    if (!is.null(simp)) {
      attr(f, "simplification") <- list(
        deviation = simp$deviation, denom_bound = simp$denom_bound,
        tolerance = tolerance, tolerance_type = tolerance_type,
        raw_beta0 = beta0, raw_terms = reduced_terms
      )
    }
  } else {
    attr(f, "equation") <- sprintf("mu = g^-1(eta)  [link = %s; printable equation unavailable, see warning]", spec$link)
  }
  class(f) <- c("deployed_mean_model", class(f))
  
  if (verbose) print(f)
  f
}

# ---- S3 methods -----------------------------------------------------------
print.deployed_mean_model <- function(x, ...) {
  cat("<deployed_mean_model>\n")
  cat(attr(x, "equation"), "\n")
  cat("link:", attr(x, "link"), "\n")
  cat("required input columns:", paste(attr(x, "free_vars"), collapse = ", "), "\n")
  
  co <- attr(x, "constants")
  if (length(co) > 0) {
    cat("\nConstants absorbed:\n")
    for (nm in names(co)) cat(sprintf("  %s = %s\n", nm, co[[nm]]))
  }
  
  fc <- attr(x, "fitted_validation")
  if (!is.null(fc)) {
    if (isTRUE(fc$ok)) {
      cat(sprintf("\nValidated against the model's own linear predictor (eta scale, random effects excluded, n=%d):\n",
                  fc$n))
      cat(sprintf("  max|dev| = %.4g, max relative dev = %.4g, RMSE = %.4g\n",
                  fc$max_abs_dev, fc$rel_dev, fc$rmse))
      if (isTRUE(fc$point_estimate_model)) {
        if (fc$rel_dev > 0.001) {
          cat("  ** relative deviation exceeds 0.1% — this is a point-estimate model with no averaging\n")
          cat("     involved, so this should match to floating-point precision; inspect the extraction **\n")
        }
      } else {
        cat("  Compared on the eta (linear predictor) scale specifically to avoid the Jensen's-inequality\n")
        cat("  gap that would arise from comparing on the mu scale (since eta = X %*% beta is linear in\n")
        cat("  beta, E[eta] over posterior draws equals X %*% E[beta] exactly). A small residual here\n")
        cat("  reflects only Monte Carlo sampling noise from a finite number of posterior draws, not a\n")
        cat("  structural approximation — so this check is NOT masked by that effect.\n")
        if (fc$rel_dev > 0.01) {
          cat("  ** relative deviation exceeds 1% — larger than typical MCMC noise; inspect the extraction **\n")
        }
      }
    } else {
      cat(sprintf("\nFitted-value self-check NOT performed: %s\n", fc$reason))
    }
  }
  
  s <- attr(x, "simplification")
  if (!is.null(s)) {
    cat(sprintf("\nCoefficients simplified as ratios of whole numbers (denominator <= %d)\n", s$denom_bound))
    cat(sprintf("Deviation vs original (unsimplified) model: %.4g (%s tolerance requested: %.4g)\n",
                s$deviation, s$tolerance_type, s$tolerance))
  }
  invisible(x)
}

coef.deployed_mean_model <- function(object, ...) {
  co <- attr(object, "coefficients")
  if (is.null(co)) return(NULL)
  term_disp <- vapply(co$terms, function(t) {
    piece_strs <- vapply(t$pieces, .format_piece_display, character(1))
    paste(piece_strs, collapse = ":")
  }, character(1))
  stats::setNames(c(co$intercept, vapply(co$terms, `[[`, numeric(1), "coef")),
                  c("(Intercept)", term_disp))
}

deployed_equation  <- function(x) attr(x, "equation")
deployed_deviation <- function(x) {
  s <- attr(x, "simplification")
  if (is.null(s)) NA_real_ else s$deviation
}
deployed_fit_check <- function(x) attr(x, "fitted_validation")

# =============================================================================
# NOTES / VERIFICATION STATUS (2026-07-28, v7 — fixes the phi_Intercept
# error on the ordered-beta brms model with a phi ~ AmbientEnv submodel):
#
# 0. The dpar-filtering fix relies on model$formula$pforms (a brmsformula
#    object's list of additional distributional-parameter formulas, e.g.
#    $pforms$phi for "phi ~ AmbientEnv") to know exactly which coefficient-
#    name prefixes to exclude. I have not verified this structure by
#    execution — if it's wrong or missing, the code falls back to a
#    heuristic list of common brms dpar names (sigma, phi, zi, hu, shape,
#    nu, xi, kappa, disc, bs, ndt, bias, alpha, beta, quantile) matched as
#    coefficient-name prefixes, with a warning telling you it did so. Given
#    your model has "kappa" and "xi" as further distributional parameters
#    too (per the summary() output), but those don't appear to enter fixef()
#    as prefixed *coefficients* (they're plain scalar parameters without
#    their own formula in your model), this shouldn't be an issue here —
#    but worth checking the printed equation only contains mu-relevant terms
#    once you re-run this.
#
# 1. glmgee: extraction verified against glmtoolbox source, as before.
#
# 2. Your last run showed two symptoms of the SAME cause: glmgee's internal
#    na.omit dropped 78 rows (3948 raw rows -> 3870 fitted rows). The
#    recovered raw data (3948 rows) therefore didn't row-align with
#    model$fitted.values (3870), AND contained NA rows that turned every
#    downstream numeric comparison into NA (not an error, but not a real
#    number either) — which is exactly why the cross-check warned even with
#    coef_simplify=FALSE, and why the simplification loop's if(...) choked
#    on "missing value where TRUE/FALSE needed".
#      FIX (a): the fitted-value self-check now uses the model's own stored
#        fit frame (model$model for glmgee / model.frame(model) for glmmTMB
#        / model$data for brms) directly with stats::model.matrix() — this
#        is exactly the post-na.omit, already-transformed frame the package
#        itself used, so it is guaranteed row-aligned with fitted values and
#        needs no raw columns at all.
#      FIX (b): the internal cross-check and coef_simplify's validation data
#        are now filtered to complete cases (stats::complete.cases()) over
#        the required raw variables before comparison, so NA rows no longer
#        propagate into the tolerance check.
#      FIX (c): the simplification loop now treats an NA deviation as a
#        reason to stop and warn, rather than crashing on if(NA).
#    None of this has been executed by me — please re-run and confirm the
#    fitted-value check now reports ok=TRUE with a sensible max|dev| for
#    both m1 and m2, and that coef_simplify=TRUE completes for m1.
#
# 3. stats::fitted() (glmmTMB) / brms::fitted() (brmsfit) are still not
#    verified by execution on my end.
#
# 4. Everything else from the v3 NOTES (raw-data recovery via model$call$data
#    assuming the original data object is unmodified in its environment;
#    contrasts assumed to be default treatment contrasts; multi-variable
#    I() expressions not supported; Rscript unavailable in this sandbox so
#    nothing has been executed end-to-end) still applies.
# =============================================================================