# =============================================================================
# deploy_mean_model.R  (v3)
#
# A model-type-agnostic "deployment" toolkit that extracts the FIXED-EFFECTS
# MEAN STRUCTURE of a fitted model (glmtoolbox::glmgee, glmmTMB, or brms) and
# turns it into a standalone R function of the predictor variables.
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
  
  env <- tryCatch(environment(model[["formula"]]), error = function(e) NULL)
  if (is.null(env) || !is.environment(env)) env <- tryCatch(attr(terms_obj, ".Environment"), error = function(e) NULL)
  if (is.null(env) || !is.environment(env)) env <- tryCatch(environment(stats::formula(model)), error = function(e) NULL)
  if (is.null(env) || !is.environment(env)) env <- parent.frame(2)
  
  out <- tryCatch(eval(cl[["data"]], envir = env), error = function(e) NULL)
  if (!is.data.frame(out)) return(NULL)
  out
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
  fitted_vals <- tryCatch(as.numeric(model[["fitted.values"]]), error = function(e) NULL)
  
  list(coefficients = cf, link = link, terms = terms_obj, data = raw_data,
       xlevels = xlevels, contrasts = contrasts_arg, fitted = fitted_vals,
       data_is_raw = data_is_raw)
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
  
  fitted_vals <- tryCatch(as.numeric(stats::fitted(model)), error = function(e) NULL)
  
  list(coefficients = cf, link = link, terms = terms_obj, data = raw_data,
       xlevels = NULL, contrasts = NULL, fitted = fitted_vals, data_is_raw = data_is_raw)
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
  
  raw_data <- .recover_raw_data(model, terms_obj)
  data_is_raw <- .data_has_required_vars(raw_data, terms_obj)
  if (!data_is_raw) {
    warning("Could not recover the original raw fitting data for this brmsfit model (tried evaluating ",
            "model$call$data). Falling back to model$data, which may already be in a transformed shape and ",
            "could fail for any I(...)-transformed term. Supply 'validation_data' explicitly for coef_simplify.")
    raw_data <- tryCatch(model[["data"]], error = function(e) NULL)
  }
  if (is.null(raw_data)) stop("Could not recover any training data frame (model$data) from this brmsfit object.")
  
  fitted_vals <- tryCatch({
    fv <- brms::fitted(model, robust = (point_estimate == "median"))
    as.numeric(fv[, "Estimate"])
  }, error = function(e) NULL)
  
  list(coefficients = cf, link = link, terms = terms_obj, data = raw_data,
       xlevels = NULL, contrasts = NULL, fitted = fitted_vals, data_is_raw = data_is_raw)
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
  
  ref_mm <- stats::model.matrix(terms_obj, data)
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
    if (deviation <= tolerance || denom_bound >= max_denominator) break
  }
  
  if (deviation > tolerance) {
    warning(sprintf("Simplification did not reach tolerance %.4g within max_denominator = %d; achieved deviation = %.4g.",
                    tolerance, max_denominator, deviation))
  }
  list(beta0 = approx_beta0, terms = approx_terms, deviation = deviation, denom_bound = denom_bound)
}

# ---- Validation against the model's own fitted values ---------------------
# Runs unconditionally (regardless of coef_simplify): reconstructs mu on the
# raw training data using ONLY spec$coefficients/spec$terms/spec$link (i.e.
# ignoring any constants requested in THIS deploy_mean_model() call) and
# compares against the model's own in-sample fitted values. This checks that
# coefficient/link/formula extraction is faithful, independent of whatever
# reduction the user is additionally asking for.
.validate_against_fitted <- function(spec) {
  if (is.null(spec$fitted)) {
    return(list(ok = FALSE, reason = "Model's fitted values could not be extracted."))
  }
  if (!isTRUE(spec$data_is_raw)) {
    return(list(ok = FALSE, reason = "Raw training data unavailable (see extraction warning above); self-check skipped."))
  }
  eta <- tryCatch(.eval_eta(spec, spec$data), error = function(e) e)
  if (inherits(eta, "error")) {
    return(list(ok = FALSE, reason = sprintf("Reconstruction failed: %s", conditionMessage(eta))))
  }
  inv_link <- .get_inverse_link(spec$link)
  mu_hat <- inv_link(eta)
  if (length(mu_hat) != length(spec$fitted)) {
    return(list(ok = FALSE, reason = sprintf(
      "Row-count mismatch (reconstructed: %d, model$fitted: %d) — likely due to rows dropped for missingness during fitting; self-check skipped.",
      length(mu_hat), length(spec$fitted))))
  }
  dev <- mu_hat - spec$fitted
  list(ok = TRUE, max_abs_dev = max(abs(dev)), rmse = sqrt(mean(dev^2)), n = length(dev))
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

.format_equation <- function(beta0, reduced_terms, link, digits = 4) {
  term_strs <- vapply(reduced_terms, function(t) {
    piece_strs <- vapply(t$pieces, .format_piece_display, character(1))
    paste0(signif(t$coef, digits), "*", paste(piece_strs, collapse = "*"))
  }, character(1))
  
  eta_terms <- c(as.character(signif(beta0, digits)), term_strs)
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
    check_ok <- tryCatch({
      check_mu_robust  <- inv_link(.eval_eta(spec, spec$data))
      check_mu_reduced <- .eval_reduced_mu(beta0, reduced_terms, spec$data, spec$link)
      max(abs(check_mu_robust - check_mu_reduced)) <= 1e-6
    }, error = function(e) NA)
    if (!isTRUE(check_ok)) {
      warning("Could not confirm that the reduced/printable equation numerically matches the ",
              "stats::model.matrix()-based predictions on the training data (either the check failed to run, ",
              "or the two disagreed by more than 1e-6). The returned function's predictions always use the ",
              "model.matrix()-based path and should still be correct, but treat the printed equation and any ",
              "coefficient simplification with caution until you've verified deployed(newdata) against ",
              "predict() from the original model.")
    }
    
    if (coef_simplify) {
      val_data <- if (!is.null(validation_data)) validation_data else spec$data
      if (!is.null(validation_data)) {
        # user-supplied data: nothing further required, used as-is
      } else if (!isTRUE(spec$data_is_raw)) {
        stop("coef_simplify = TRUE requires raw training data with the original (untransformed) columns, ",
             "which could not be recovered automatically for this model (see warning above). Supply it via ",
             "the 'validation_data' argument.")
      }
      required_vars <- unique(unlist(lapply(reduced_terms, function(t) vapply(t$pieces, `[[`, character(1), "raw_var"))))
      missing_cols <- setdiff(required_vars, names(val_data))
      if (length(missing_cols) > 0) {
        stop(sprintf("Validation data is missing required column(s): %s", paste(missing_cols, collapse = ", ")))
      }
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
    attr(f, "equation") <- .format_equation(disp_beta0, disp_terms, spec$link)
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
      cat(sprintf("\nValidated against model's fitted values (n=%d): max|dev| = %.4g, RMSE = %.4g\n",
                  fc$n, fc$max_abs_dev, fc$rmse))
      if (fc$max_abs_dev > 1e-4) {
        cat("  ** deviation exceeds 1e-4 on the mu scale — inspect the extraction before trusting this deployment **\n")
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
# NOTES / VERIFICATION STATUS (2026-07-27, updated after your test run):
#
# 1. glmgee: extraction verified against glmtoolbox source, as before.
#
# 2. Raw-data recovery (.recover_raw_data) is the fix for the bug you hit.
#    It works by evaluating model$call$data in the formula's environment —
#    this depends on the original data object (e.g. m2Data) still existing,
#    unmodified, in that environment. This is a common and generally
#    reliable pattern, but I have not executed it myself. If it silently
#    picks up a MODIFIED version of m2Data/m1Data (e.g. you overwrote the
#    variable after fitting), the self-check and simplification would be
#    validated against the wrong data. Given you're calling this in the same
#    session right after fitting, this should be fine — just flagging the
#    assumption.
#
# 3. The new unconditional fitted-value check compares against
#    model$fitted.values (glmgee, confirmed from source) or stats::fitted()
#    (glmmTMB) / brms::fitted() (brmsfit) — the latter two not verified by
#    execution. If either errors or returns something unexpected, the check
#    reports "NOT performed" with a reason rather than crashing the whole
#    deployment — please report what you see for m1/m2 now.
#
# 4. lme4::nobars() worked for m1 in your last run (deprecation warning
#    aside), which is a good sign the formula/terms extraction for glmmTMB
#    is sound. I've switched to preferring reformulas::nobars() to avoid
#    the noisy deprecation warning, falling back to lme4::nobars().
#
# 5. Everything else from the v2 NOTES (contrasts assumption, multi-variable
#    I() not supported, Rscript unavailable in this sandbox so nothing has
#    been executed end-to-end) still applies.
# =============================================================================