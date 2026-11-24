require(glmmTMB)

# =============================================================================
# predict_glmmTMB_ordbeta.R  (v4)
#
# v4: signature identical to predict.glmmTMB (glmmTMB source supplied by the
#     user): allow.new.levels = NULL, type includes "probs", zitype, fast,
#     debug, aggregate, do.bias.correct, bias.correct.control. Everything other
#     than ordbeta type = "response" is passed through with all arguments. For
#     ordbeta E[y]: aggregate, do.bias.correct and debug are refused (they would
#     act on eta, not E[y]); the numerical derivatives never use fast = TRUE
#     (incompatible with newparams in glmmTMB); cov.fit sets se.fit with the
#     same message as glmmTMB.
# v3: follows glmmTMB's own prediction route. newparams is the fixed-parameter
#     vector (object$fit$par), as in glmmTMB (v1/v2 wrongly took the full vector,
#     which only coincides with it in models without random effects). Standard
#     errors reproduce TMB's generalised delta method, so for models with random
#     effects they include the random-effect uncertainty, as glmmTMB's se.fit
#     does. cov.fit is supported. The v2 reproduction check is no longer needed.
# v2: tolerance fix; warning for se.fit with random effects (superseded).
#
# predict.glmmTMB.ordbeta(): predictions for glmmTMB ordered-beta (ordbeta)
# models with a correct type = "response".
#
# THE PROBLEM
#   For ordbeta, glmmTMB's predict(type = "response") returns mu = linkinv(eta),
#   the mean of the beta (continuous) component only. The expected response
#   also includes the point masses at 0 and 1 (glmmTMB issue #1359; fix in PR
#   #1373, unmerged at the time of writing; once merged, type = "response" will
#   change meaning in glmmTMB itself):
#
#     P(y = 0)     = 1 - plogis(eta - psi1)
#     P(y = 1)     = plogis(eta - psi2)
#     P(0 < y < 1) = plogis(eta - psi1) - plogis(eta - psi2)
#     E[y]         = P(y = 1) + P(0 < y < 1) * mu
#
#   (likelihood in glmmTMB.cpp, ordbeta case; the same as ordbetareg's
#   posterior_epred with component = "all"). The cutpoints psi are on the logit
#   scale, whatever the model's link; mu uses the model's link. They are read
#   from the model's parameter list, which includes cutpoints fixed via map.
#   family_params() is NOT used: it reports a transformed (and, for an
#   unidentified cutpoint, different) value.
#
# WHAT THIS FUNCTION DOES
#   * Same arguments and defaults as predict.glmmTMB. Call it directly (it is not dispatched
#     automatically, because the model's class is just "glmmTMB").
#   * type = "response": E[y] from glmmTMB's own linear predictor
#     (predict(type = "link") with the same newdata, newparams, re.form,
#     allow.new.levels and na.action), and the cutpoints at the same fixed
#     parameters (object$fit$par, or newparams).
#   * se.fit / cov.fit: TMB's generalised delta method, as used by glmmTMB
#     (sdreport documentation):
#       Var(E) = J V J' + a^2 * C
#     J: total derivative of E with respect to the fixed parameters, with the
#        random effects re-optimised at each parameter value (numerical, through
#        predict(newparams = ...), which re-optimises them as glmmTMB does);
#     V: inverse of optimHess() at the fixed parameters, as in glmmTMB's predict;
#     a: dE/d(eta) at fixed cutpoints;
#     C: the random-effect part of Var(eta), = glmmTMB's own link-scale
#        variance (from sdreport) minus J_eta V J_eta'. Only computed when the
#        model has random parameters (random effects, or fixed effects under
#        REML); otherwise C = 0 and this is the ordinary delta method.
#     E depends on the random effects only through eta, so this equals the
#     sdreport formula; the only differences are numerical derivatives.
#   * Every other type, and every non-ordbeta model, is passed straight to
#     glmmTMB's predict(). For ordbeta, type = "conditional" (and the current
#     type = "response" in glmmTMB itself) gives mu.
#   * As in glmmTMB, re.form = NULL predicts at the random effects' conditional
#     modes and re.form = NA at zero. Neither gives marginal (population-averaged)
#     predictions that integrate over the random-effects distribution.
#
# CHECKS
#   * zero-inflated ordbeta models (a ziformula on top of ordbeta's own zero
#     mass; glmmTMB permits it technically) are refused;
#   * newparams must have the length of object$fit$par;
#   * two cutpoints with psi1 < psi2;
#   * a singular Hessian gives a warning and NA standard errors;
#   * aggregate, do.bias.correct = TRUE and debug = TRUE are refused for E[y].
#   If se.fit fails in glmmTMB's fast path ("Problem with these data entries:
#   aggregate"), pass fast = FALSE; it is used for glmmTMB's own calls here.
#
# Not verified by execution (no R runtime was available when this was written).
#
# Example
#   E  <- predict.glmmTMB.ordbeta(m_dA1, newdata = nd, type = "response", re.form = NA)
#   Es <- predict.glmmTMB.ordbeta(m_dA1, newdata = nd, type = "response", re.form = NA, se.fit = TRUE)
#   Es$fit; Es$se.fit
# =============================================================================

predict.glmmTMB.ordbeta <- function(object, newdata = NULL, newparams = NULL, se.fit = FALSE, cov.fit = FALSE,
                                    re.form = NULL, allow.new.levels = NULL,
                                    type = c("link", "response", "conditional", "zprob", "zlink",
                                             "disp", "latent", "probs"),
                                    zitype = NULL, na.action = stats::na.pass, fast = NULL, debug = FALSE,
                                    aggregate = NULL, do.bias.correct = FALSE,
                                    bias.correct.control = list(sd = TRUE), ...) {
  if (!inherits(object, "glmmTMB")) stop("object must be a glmmTMB model.", call. = FALSE)
  if (!is.null(zitype)) {
    warning("zitype is deprecated: please use type instead")
    type <- zitype
  }
  type <- match.arg(type)
  
  # everything except E[y] for ordbeta: glmmTMB's own method, unchanged
  if (!identical(stats::family(object)$family, "ordbeta") || type != "response") {
    return(stats::predict(object, newdata = newdata, newparams = newparams, se.fit = se.fit, cov.fit = cov.fit,
                          re.form = re.form, allow.new.levels = allow.new.levels, type = type,
                          na.action = na.action, fast = fast, debug = debug, aggregate = aggregate,
                          do.bias.correct = do.bias.correct, bias.correct.control = bias.correct.control, ...))
  }
  
  if (length(glmmTMB::fixef(object)$zi) > 0) {
    stop("Zero-inflated ordbeta models are not supported (E[y] would need the zero-inflation component).", call. = FALSE)
  }
  if (length(aggregate) > 0 || isTRUE(do.bias.correct) || isTRUE(debug)) {
    stop("aggregate, do.bias.correct and debug are not supported for ordbeta type = \"response\".", call. = FALSE)
  }
  if (isTRUE(cov.fit)) {
    if (!isTRUE(se.fit)) message("se.fit set to TRUE because cov.fit = TRUE")
    se.fit <- TRUE
  }
  
  env <- object$obj$env
  lf <- env$lfixed()                                     # positions of the fixed parameters in the full vector
  theta0 <- if (is.null(newparams)) object$fit$par else newparams
  if (length(theta0) != sum(lf)) {
    stop("newparams must be a fixed-parameter vector like object$fit$par.", call. = FALSE)
  }
  linkinv <- stats::family(object)$linkinv
  
  eta_at <- function(theta) {
    as.numeric(stats::predict(object, newdata = newdata, newparams = theta, re.form = re.form,
                              allow.new.levels = allow.new.levels, type = "link", na.action = na.action,
                              fast = NULL, ...))                  # fast = TRUE is incompatible with newparams
  }
  psi_at <- function(theta) {
    full <- env$last.par.best
    full[lf] <- theta
    psi <- env$parList(theta, full)$psi                  # includes cutpoints fixed via map
    if (length(psi) != 2L || !all(is.finite(psi)) || !(psi[1] < psi[2])) {
      stop("Could not read two ordbeta cutpoints with psi1 < psi2 from the parameter list.", call. = FALSE)
    }
    psi
  }
  E_of <- function(eta, psi) {
    p1 <- stats::plogis(eta - psi[2])
    p1 + (stats::plogis(eta - psi[1]) - p1) * linkinv(eta)
  }
  
  # estimate: glmmTMB's own linear predictor (newparams passed through unchanged)
  eta0 <- as.numeric(stats::predict(object, newdata = newdata, newparams = newparams, re.form = re.form,
                                    allow.new.levels = allow.new.levels, type = "link", na.action = na.action,
                                    fast = fast, ...))
  psi0 <- psi_at(theta0)
  fit <- E_of(eta0, psi0)
  if (!isTRUE(se.fit)) return(fit)
  
  # V: as in glmmTMB's predict, the inverse Hessian of the objective at the fixed parameters
  H <- stats::optimHess(theta0, object$obj$fn, object$obj$gr)
  V <- tryCatch(solve(H), error = function(e) {
    warning("The Hessian of the fixed parameters is singular; standard errors are NA.", call. = FALSE)
    matrix(NA_real_, length(theta0), length(theta0))
  })
  
  # J_E, J_eta: total derivatives with respect to the fixed parameters
  J_E <- J_eta <- matrix(NA_real_, length(fit), length(theta0))
  for (j in seq_along(theta0)) {
    h <- 1e-4 * max(1, abs(theta0[j]))
    up <- theta0; up[j] <- up[j] + h
    dn <- theta0; dn[j] <- dn[j] - h
    eta_up <- eta_at(up); eta_dn <- eta_at(dn)
    J_E[, j] <- (E_of(eta_up, psi_at(up)) - E_of(eta_dn, psi_at(dn))) / (2 * h)
    J_eta[, j] <- (eta_up - eta_dn) / (2 * h)
  }
  JV_E <- J_E %*% V
  
  if (length(env$random) > 0) {
    # random-effect part of Var(eta): glmmTMB's link-scale (co)variance minus its fixed-parameter part
    lk <- stats::predict(object, newdata = newdata, newparams = newparams, se.fit = TRUE, cov.fit = cov.fit,
                         re.form = re.form, allow.new.levels = allow.new.levels, type = "link",
                         na.action = na.action, fast = fast, ...)
    hs <- 1e-5 * pmax(1, abs(eta0))
    a <- (E_of(eta0 + hs, psi0) - E_of(eta0 - hs, psi0)) / (2 * hs)
    if (isTRUE(cov.fit)) {
      if (is.null(lk$cov.fit)) stop("glmmTMB's predict did not return cov.fit for the linear predictor.", call. = FALSE)
      C <- lk$cov.fit - J_eta %*% V %*% t(J_eta)
      cov_E <- JV_E %*% t(J_E) + outer(a, a) * C
    } else {
      C <- pmax(lk$se.fit^2 - rowSums((J_eta %*% V) * J_eta), 0)
      var_E <- rowSums(JV_E * J_E) + a^2 * C
    }
  } else if (isTRUE(cov.fit)) {
    cov_E <- JV_E %*% t(J_E)
  } else {
    var_E <- rowSums(JV_E * J_E)
  }
  
  if (isTRUE(cov.fit)) {
    return(list(fit = fit, se.fit = sqrt(pmax(diag(cov_E), 0)), cov.fit = cov_E))
  }
  list(fit = fit, se.fit = sqrt(pmax(var_E, 0)))
}