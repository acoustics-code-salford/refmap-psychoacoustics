require(future)         # future::plan() for parallel = TRUE (set the plan yourself)
require(future.apply)   # future.apply::future_lapply() for parallel = TRUE
require(glmmTMB)        # glmmTMB models
require(glmtoolbox)     # glmgee models (only needed if you bootstrap one)

# =============================================================================
# boot_mean_model.R  (v15)
#
# Cluster-bootstrap uncertainty intervals for a stimulus-level model that has
# been deployed with deploy_mean_model(). Source deploy_mean_model.R (v17+)
# first: this file uses its internal functions.
#
# -----------------------------------------------------------------------------
# USAGE
#   boot_mean_model(deployed, trial_data, cluster, response,
#                   stat = c("median", "mean"), by = NULL, strata = NULL,
#                   response_transform = NULL, newdata = NULL, B = 2000L,
#                   level = 0.95, prediction_intervals = FALSE, pi_draws = 10L,
#                   parallel = FALSE, seed = NULL, back_transform = TRUE,
#                   jitter_rows = NULL, jitter_step = 1, jitter_bw = 0.25,
#                   recentre = c("mean", "median", "none"), regrid = NULL,
#                   interval = c("percentile", "bc"),
#                   pi_dispersion = c("original", "replicate"), max_fail = 0.05,
#                   allow_fixed_weights = FALSE,
#                   resample = c("clusters", "crossed"), stimulus_strata = NULL,
#                   verbose = TRUE)
#
# ARGUMENTS
#   deployed   Object returned by deploy_mean_model() (v17+, which stores the
#              model and the call's arguments). glmmTMB and glmgee models only.
#              Its constants (e.g. Experiment = 0), centring/standardisation,
#              factor contrast codes and frozen data-dependent terms are all
#              reused. The stimulus-level data are taken
#              from the model, so they are not an argument.
#
#   trial_data Long-format trial-level data, one row per trial. Must contain the
#              columns named in cluster, response, by and strata.
#
#   cluster    Name of the column identifying the resampling unit (e.g. the
#              participant ID). Each distinct value is one unit and all of its
#              trials travel together. No missing values. A person with the SAME
#              ID in several experiments is one cluster; if IDs were reused for
#              DIFFERENT people across experiments, make them unique first (the
#              printed strata sizes show how many clusters span several strata).
#
#   response   Trial-level column to aggregate, on the scale whose median/mean
#              (after response_transform) gives the model's response column.
#              Converted to numeric; NAs are ignored.
#
#   stat       "median" (default), "mean", or a function receiving the NA-free
#              values of one stimulus row and returning one number. Note: the
#              median commutes with any order-preserving transform, the mean
#              only with linear ones; the reproduction check (below) catches a
#              wrong combination.
#
#   by         Column(s) identifying a stimulus row; must exist in trial_data and
#              in the model's stimulus-level data, and identify that data's rows
#              uniquely. Default: the model's predictor variables (raw names for
#              centred/standardised ones). If those do not identify the rows
#              uniquely (rows differing only in variables not in the model), the
#              function stops: give e.g. a stimulus ID. Numeric keys are matched
#              to 12 significant digits.
#
#   strata     NULL (default) or column name(s) in trial_data. Each cluster's
#              stratum is the set of values it has across its trials (e.g.
#              "1", "2" or "1 + 2" for strata = "Experiment"), and every stratum
#              keeps its original number of clusters in each replicate. Use it
#              when group sizes were fixed by design.
#
#   response_transform
#              NULL or a function applied to the aggregated values before they
#              replace the model's response column. Needed when that column was
#              rescaled in the data (e.g. stored as median / 10), not when the
#              formula does the rescaling. If the response is a non-negative
#              quantity built from a variable that can be negative (e.g. an
#              increase defined from a [-10, 10] difference), clip it here, e.g.
#              function(x) pmax(x, 0) / 10: resampled panels can produce
#              negative medians even when the original data have none.
#
#   newdata    Data frame of raw predictor values: the deployed function's
#              required inputs (attr(deployed, "free_vars")). The deployment's
#              constants are applied. Default: the model's stimulus-level data.
#
#   B          Number of bootstrap replicates (default 2000). Percentile limits
#              with B below ~1000 are imprecise; use small B only for testing.
#
#   level      Coverage of the confidence and prediction intervals (0.95).
#
#   prediction_intervals
#              TRUE to also compute prediction intervals for a NEW stimulus-level
#              value (a panel aggregate like those fitted): glmmTMB ordbeta, beta
#              or gaussian models without random effects. The dispersion is
#              evaluated with the deployment's constants (a factor fixed at a
#              contrast code is averaged on the log-dispersion scale).
#
#   pi_draws   Draws per replicate per newdata row for the prediction intervals
#              (default 10); pooled over replicates. Memory: about
#              nrow(newdata) * B * pi_draws * 8 bytes.
#
#   parallel   TRUE to use future.apply::future_lapply() with the current
#              future::plan(), e.g. future::plan(future::multisession, workers = 6).
#              All random numbers are drawn in the main session, so results do
#              not depend on the plan or the number of workers.
#
#   seed       NULL (use the current random-number stream) or a number (the run
#              is reproducible and the global RNG state is restored afterwards).
#
#   back_transform
#              TRUE (default): predictions are returned on the scale of the
#              aggregated trial response (e.g. median dAnnoyance), by inverting
#              g(x) = formula LHS evaluated on response_transform(x). g is
#              examined over the trial response's range: flat ends (clipping)
#              are allowed and reported; an affine g is inverted exactly (expected
#              values carry over); a monotone NON-LINEAR g is inverted
#              numerically and flagged (limits are exact, but the estimate is
#              then the back-transformed mean, not the mean on that scale); a
#              non-monotone g is not inverted. FALSE: model's response scale.
#
#   jitter_rows, jitter_step, jitter_bw, recentre
#              Smoothed bootstrap for responses on a discrete grid (e.g. integer
#              ratings), whose medians stick and jump under plain resampling so
#              the bootstrap misstates their variability and bias. jitter_rows:
#              NULL (default, off), TRUE (all rows), a logical vector (one per
#              trial row) or the name of a logical column marking the
#              discrete-scale trials (e.g. Experiment 1). jitter_step: grid step
#              (1 for integer ratings). jitter_bw: jitter SD as a fraction of the
#              step (0.25). In each replicate the marked resampled values are
#              shrunk towards their stimulus's original aggregate by
#              1 / sqrt(1 + sigma^2 / var) and Gaussian jitter is added. recentre
#              ("mean", default, or TRUE; "median"; "none" or FALSE): each
#              affected stimulus's replicate aggregates are then shifted so their
#              mean (or median) equals its original aggregate. "mean" suits a
#              refit (fits respond to the average of their data, and regrid
#              preserves means, not medians); "median" matches the measurement-
#              data function, whose target is a percentile interval for a single
#              median. regrid (NULL; TRUE = jitter_step / 2, for stat =
#              "median"; or a grid step such as 0.5 for medians of integer
#              ratings; checked: the original aggregates of the smoothed stimuli
#              must lie on that grid): the
#              affected stimuli's replicate aggregates are then rounded
#              stochastically (mean-preserving) back onto that grid, so they take
#              the same values real aggregates take (including exact 0s after
#              clipping). Without it, smoothing turns exact grid values into a
#              continuum (e.g. tiny positive values instead of exact 0s), which
#              a model with a separate 0/1 component, like ordbeta, reacts to. Applied
#              before response_transform; point estimates are unchanged. Same
#              procedure as bias_corrected_jitter_median() for the measurement
#              data, applied at the stimulus level.
#
#   interval   "percentile" (default) or "bc": bias-corrected percentile
#              (Efron), which shifts the quantiles by z0 = qnorm(p_below) to
#              allow for a replicate distribution that is off-centre relative to
#              the estimate (e.g. an estimator pulled by a floor). BC uses more
#              extreme quantiles: use B >= 2000. Applies to the confidence
#              intervals (predictions and coefficients), not to the prediction
#              intervals.
#
#   pi_dispersion
#              "original" (default): prediction-interval draws use the original
#              fit's dispersion, with the replicates supplying the uncertainty in
#              eta and the cutpoints. "replicate": each replicate's own
#              dispersion. Re-aggregated replicates carry a second dose of panel
#              noise, so their dispersion is biased towards more spread and
#              "replicate" overstates the prediction intervals.
#
#   max_fail   Largest acceptable share of failed replicates (default 0.05). Above
#              it a warning is issued; above 4 x max_fail the function stops,
#              because failures can depend on the resampled data, so the
#              survivors may not represent the bootstrap distribution.
#              Replicates without a positive-definite Hessian are NOT failures:
#              they are kept and counted (with an unidentified ordbeta cutpoint,
#              e.g. psi2 when there are no 1s, almost all replicates are like that).
#
#   allow_fixed_weights
#              FALSE (default): models with weights are refused, because the
#              weights are kept from the stimulus-level data, not recomputed per
#              replicate. TRUE only for fixed design weights.
#
#   resample   "clusters" (default): resample clusters (participants) only; the
#              stimuli are fixed. "crossed": also resample the model's stimulus
#              rows with replacement, independently of the clusters, in each
#              replicate (the pigeonhole bootstrap for crossed designs; Owen 2007,
#              Ann. Appl. Stat. 1, 386-411). See INTERPRETATION.
#
#   stimulus_strata
#              resample = "crossed" only. NULL (default) or column name(s) in the
#              model's stimulus-level data; stimulus rows are resampled within
#              these strata, each keeping its original number of rows. Use it for
#              stimulus sets fixed by design, e.g. "Experiment": without it a
#              replicate can lose a whole experiment's stimuli (its factor level),
#              and such refits fail.
#
#   verbose    Print progress and the summary.
#
# VALUE (class "boot_mean_model")
#   $predictions    newdata identifiers and predictors, estimate (E[y] from the
#                   original fit), bias (mean of replicates - estimate), p_below
#                   (share of replicates below the estimate; 0.5 = centred), lower and
#                   upper confidence limits, and pi_lower / pi_upper if requested.
#                   On the aggregated trial response's scale (see back_transform).
#   $coef_model     model coefficients (plus ordbeta cutpoints psi1/psi2 and, with
#                   prediction intervals, dispersion coefficients): estimate,
#                   bootstrap SE, bootstrap bias, percentile limits, notes.
#   $coef_equation  coefficients of the deployed (printed) equation in raw units,
#                   constants absorbed, unsimplified, on the eta scale (NOT
#                   sign-reversed for logit).
#   $replicates     matrices of replicate coefficients, cutpoints, dispersion
#                   coefficients, eta and E (rows = newdata, cols = replicates).
#   $failures       discarded replicates by reason, warnings raised during refits,
#                   usable replicates without a positive-definite Hessian.
#   $settings, $checks
#
# -----------------------------------------------------------------------------
# HOW IT WORKS
#   One replicate: resample clusters with replacement (within strata); recompute
#   the stimulus-level median/mean (NAs removed) and put it, after
#   response_transform, in place of the model's response column in the model's
#   own stimulus-level data (nothing else changes: predictors, centring, factor
#   contrasts); with resample = "crossed", the stimulus rows the model used are
#   then resampled with replacement (within stimulus_strata), a row drawn twice
#   appearing twice with the same aggregate; refit with the model's original call (formula, family, start,
#   map, dispformula, control, ...); keep the conditional-model coefficients
#   (and ordbeta cutpoints; dispersion coefficients for prediction intervals).
#
#   Predictions use the DEPLOYED path (un-scaling, constants, contrast codes,
#   frozen data-dependent terms) through a design matrix built once from
#   deploy_mean_model()'s own evaluator, one unit coefficient vector per column;
#   eta for a replicate is X_new %*% beta_b, exactly what the deployed function
#   computes with beta_b. E[y] uses the full mean map with each replicate's own
#   cutpoints (ordbeta: E = P1 + (P(y > 0) - P1) * mu). Coefficients are the
#   unsimplified estimates (the printed equation's rounded ratios do not enter).
#
#   Prediction intervals pool draws from each replicate's fitted distribution
#   (a single-uniform ordbeta sampler, so P(y = 1) = P1 exactly), so they include
#   parameter uncertainty and the spread of new stimulus-level values.
#
# CHECKS (the function STOPS rather than continue silently)
#   * aggregating the ORIGINAL trial data must reproduce the model's response
#     column exactly (suggests a multiplier when the two differ by a constant);
#   * before any refitting, every replicate's response, as the model sees it
#     (formula LHS), must lie in the family's support;
#   * the training design matrix must reproduce the model's own linear predictor,
#     and the prediction design matrix the deployed function's;
#   * refused: offsets, zero-inflation, weights (unless allow_fixed_weights),
#     non-numeric responses, missing values in 'by'/'strata' columns, partial
#     jitter_rows within a stimulus, non-vectorised response_transform,
#     non-finite back-transformed values, failures above 4 x max_fail.
#   Per replicate (discarded and counted, never silently): fit errors,
#   non-convergence, changed coefficient names, non-finite coefficients, and a
#   refit design matrix (glmmTMB) / linear predictor (glmgee) that differs from
#   the precomputed one on the rows used.
#   Reported: trial rows that match no stimulus row (and clusters with none,
#   which are dropped from the resampling pool), failures above max_fail,
#   replicates without a positive-definite Hessian (kept), estimates whose
#   replicate distribution is off-centre (share of replicates below the estimate
#   outside 0.4-0.6; percentile limits are then shifted), unidentified ordbeta
#   cutpoints, small B.
#
# INTERPRETATION
#   * Confidence intervals: uncertainty from resampling clusters (participants),
#     CONDITIONAL on the stimulus set (stimuli and their predictor values are
#     fixed), for the model as specified (the cluster bootstrap does not assume
#     the model is correct; the intervals cover its best-fitting parameters).
#     Same target as a participant-clustered GEE; not generalisation to new stimuli.
#   * resample = "crossed": uncertainty from resampling both participants and
#     stimuli, i.e. a repeat of the study with a new panel AND a new set of
#     comparable stimuli (from the same design strata). Wider than "clusters";
#     still a confidence interval for the curve, so narrower than the prediction
#     intervals. Owen (2007) shows the pigeonhole bootstrap variance tends to be
#     somewhat conservative (too large). Formula terms that depend on the data
#     (e.g. scale(), poly()) would change with the resampled stimuli; the
#     per-replicate design check discards such refits, so a high failure count
#     with "design does not match" means the formula needs fixed constants.
#   * Predictions are fixed-effects predictions (random effects set to 0), as in
#     the deployed function; for a non-identity link this is not the average
#     over the random-effects distribution.
#   * Prediction intervals are for a new stimulus-level aggregate (a panel like
#     those fitted), not individual responses. They rely on the model's response
#     distribution and so inherit any misfit: report their empirical coverage.
#     With pi_dispersion = "original" they include the replicates' uncertainty in
#     eta and the cutpoints but NOT in the dispersion (deliberate; see above).
#   * Smoothing (jitter, shrinkage), recentring and regridding are approximations
#     (the shrinkage factor is a heuristic carried over from the measurement-data
#     procedure); BC adjusts the interval limits, not the point estimate, and
#     acts on the (possibly smoothed) bootstrap distribution. Calibration by
#     simulation has not been done.
#
# LIMITATIONS
#   Percentile or BC intervals (no BCa). Single-column responses only (no
#   cbind()). Zero-inflated and brms models, models with offsets, and (unless
#   allow_fixed_weights) weighted models are refused. The glmgee convergence
#   check uses fit$converged if present (not verified for all glmtoolbox
#   versions). The back-transform inspects the transform on a grid: suitable for
#   smooth transforms and clipping, not a general mathematical inverse.
#
# EXAMPLE
#   source("deploy_mean_model.R"); source("boot_mean_model.R")
#   dep <- deploy_mean_model(m_dA1, Experiment = 0,
#                            centred_vars = list(UASLAEMaxLRScl = "UASLAEMaxLR"))
#   future::plan(future::multisession, workers = 6)
#   bt <- boot_mean_model(dep, trial_data = subjDataNoAmb.med, cluster = "PersonNum",
#                         response = "dAnnoyance", stat = "median", by = "StimNum",
#                         strata = "Experiment",
#                         response_transform = function(x) pmax(x, 0) / 10,
#                         B = 2000, seed = 1, prediction_intervals = TRUE,
#                         parallel = TRUE)
#   bt; bt$predictions; bt$coef_model; bt$coef_equation
#
# -----------------------------------------------------------------------------
# PREDICTIONS FROM AN EXISTING RUN: predict(bt, newdata, model, ...)
#   predict(object, newdata = NULL, model, target = c("deployed", "original"),
#           level = NULL, interval = NULL, prediction_intervals = FALSE,
#           pi_draws = NULL, pi_dispersion = NULL, seed = NULL, ...)
#
#   Predictions with intervals for ANY newdata (stimuli, a counterfactual grid)
#   from the replicates stored in a boot_mean_model() result, with no refitting.
#   Each replicate's coefficients (and cutpoints) give E[y] at every newdata
#   row; the limits are taken across replicates, as in boot_mean_model(). The
#   resampling never depends on newdata, so the result is identical to
#   re-running boot_mean_model() with the same seed and settings and this
#   newdata (prediction-interval draws aside, which are random).
#
#   object    A boot_mean_model() result (v11+).
#   newdata   Data frame of raw predictor values. NULL (default): the stimulus
#             rows the model was fitted to.
#             target = "deployed": the deployed function's inputs
#               (attr(deployed, "free_vars")); the deployment's constants are
#               applied (e.g. Experiment = 0), as in boot_mean_model().
#             target = "original": every variable of the original model in raw
#               units, including those the deployment fixes (e.g. Experiment, as
#               one of its levels): the original model's curve at those values,
#               e.g. for one experiment.
#   model     The deploy_mean_model() object the run was made from; it carries
#             the original model (attr(, "model")) and the deployment, so it
#             serves both targets. The run stores replicate coefficients but
#             neither of these. Checked against the run: same coefficients and
#             constants.
#   target    "deployed" (default): the deployed function; "original": the
#             original model without the deployment's constants; see newdata.
#             For models with random effects, "original" is population-level
#             (random effects 0).
#   level, interval
#             NULL: as in the run. Otherwise a new level / "percentile" / "bc".
#   prediction_intervals
#             TRUE to add pi_lower / pi_upper (glmmTMB ordbeta, beta, gaussian;
#             no random effects). With pi_dispersion = "original" (the run's
#             default) this works even if the run had no prediction intervals;
#             "replicate" needs the run's stored replicate dispersion.
#   pi_draws, pi_dispersion
#             NULL: as in the run (pi_draws 10 if the run had none).
#   seed      For the prediction-interval draws only (RNG state restored).
#
#   VALUE (class "boot_mean_predictions")
#     $predictions  newdata's identifying columns, then estimate (original fit),
#                   bias, p_below, lower, upper (and pi_lower, pi_upper), on the
#                   run's back-transformed scale; outside_range names any numeric
#                   input outside the range of the model's data (extrapolation:
#                   intervals there reflect parameter uncertainty only, under
#                   the model's assumed form).
#     $replicates   E (back-transformed; rows = newdata, cols = replicates) and
#                   eta: for contrasts between rows (e.g. experiment 1 minus
#                   experiment 2 at the same predictor values), computed
#                   replicate by replicate; see the examples.
#     $settings, $checks (design check against the deployed function, or
#                   against glmmTMB's own linear predictor for target
#                   "original"; skipped (NA) if the deployment froze
#                   data-dependent terms, which glmmTMB would recompute).
#
#   EXAMPLES
#     # 1. stimuli, each at its own experiment (model curve, not the deployment)
#     p_stim <- predict(bt, newdata = stimDataNoAmb, model = dep, target = "original")
#
#     # 2. a sweep over one predictor through the deployed equation, others fixed
#     grid <- data.frame(LAELAF50diff = seq(-10, 20, by = 0.5),
#                        UASEvents = 10,
#                        UASTonalAwSHMInt05 = median(stimDataNoAmb$UASTonalAwSHMInt05))
#     # ... one column per name in attr(dep, "free_vars")
#     p_grid <- predict(bt, newdata = grid, model = dep)
#     ggplot2::ggplot(p_grid$predictions, ggplot2::aes(LAELAF50diff, estimate)) +
#       ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.3) +
#       ggplot2::geom_line()
#
#     # 3. the same sweep for each experiment, and their difference with an interval
#     grid_e <- rbind(cbind(grid, Experiment = "1"), cbind(grid, Experiment = "2"))
#     p_e <- predict(bt, newdata = grid_e, model = dep, target = "original")
#     n <- nrow(grid)
#     D <- p_e$replicates$E[1:n, ] - p_e$replicates$E[n + 1:n, ]      # replicate by replicate
#     diff_e <- data.frame(grid,
#                          estimate = p_e$predictions$estimate[1:n] - p_e$predictions$estimate[n + 1:n],
#                          lower = apply(D, 1, quantile, 0.025, na.rm = TRUE),
#                          upper = apply(D, 1, quantile, 0.975, na.rm = TRUE))
#
#     # check: identical to a re-run with newdata (same seed and settings)
#     # bt2 <- boot_mean_model(dep, ..., newdata = grid, seed = 1)
#     # all.equal(bt2$predictions$lower, p_grid$predictions$lower)
#
# CHANGE LOG
#   v15: predict() method for boot_mean_model results: predictions and intervals
#        for any newdata, deployed function or original model, from the stored
#        replicates (no refitting).
#   v14: resample = "crossed" (participants and stimuli, pigeonhole bootstrap)
#        with stimulus_strata. The worker now receives the replicate's stimulus
#        rows with its response; the per-replicate design check maps the refit's
#        rows back to the original stimulus rows.
#   v13: regrid = TRUE (= jitter_step / 2, medians only); regrid is checked
#        against the original aggregates (must lie on its grid; warning if finer
#        than jitter_step / 2; note if the whole-step grid fits).
#   v12: recentre now defaults to "mean" (TRUE = "mean").
#   v11: response to the third-party review. Stops for: zero-inflated glmmTMB
#        models; weights (unless allow_fixed_weights = TRUE); a non-numeric
#        (e.g. factor) response; missing values in 'by' or 'strata' columns; a
#        stimulus whose trials are only partly marked by jitter_rows; a
#        non-vectorised response_transform; non-finite back-transformed values;
#        zero-row newdata; failures above 4 x max_fail (warning above max_fail).
#        Clusters with no matching trials are dropped from the resampling pool.
#        Offsets also detected from the terms. Scalar argument validation.
#        glmgee: fit$converged checked if present. Interpretation notes added.
#   v10: interval = "bc" (bias-corrected percentile). pi_dispersion: prediction
#        intervals now use the original fit's dispersion by default. p_below
#        counts ties as half. Fixed (since v6): confidence limits were swapped
#        for a decreasing back-transform although they are taken from the
#        back-transformed replicates (already ordered); PI limits, which are
#        back-transformed quantiles, still swap correctly.
#   v9: regrid: mean-preserving stochastic rounding of smoothed replicate
#       aggregates back onto the real aggregates' grid (smoothing alone gave
#       fewer exact 0s and spurious tiny values, which shifted ordbeta's psi1,
#       dispersion and slope in the refits).
#   v8: back-transform: the clip kink no longer makes an affine transform
#       look non-linear (flat-end points excluded from the linearity test).
#       recentre = "median" (default; TRUE), "mean" or "none" (FALSE).
#       Replicate stimulus responses stored ($replicates$response_raw,
#       $replicates$response_original). New p_below columns (share of
#       replicates below the estimate); off-centre flags now use p_below
#       (outside 0.4-0.6) instead of mean bias / SE.
#   v7: optional smoothed bootstrap for discrete-scale responses (jitter_rows,
#       jitter_step, jitter_bw, recentre), mirroring the measurement-data
#       bias_corrected_jitter_median(). Replicate responses are now built in the
#       main session (all random numbers there); workers only refit.
#   v6: response_factor support removed (dropped from deploy_mean_model());
#       predictions are instead back-transformed by inverting response_transform
#       (and any transform in the formula LHS): affine exactly, monotone
#       non-linear numerically and flagged, non-monotone not at all. New
#       argument back_transform.
#   v5: documentation block; checks that deploy_mean_model.R is loaded, that
#       pi_draws is valid and that the model has no offset; reports unmatched
#       trial rows/clusters, model weights, bootstrap bias (new boot_bias /
#       bias columns, flagged when > 0.25 SE), unidentified ordbeta cutpoints and
#       small B; predictions now include the 'by' identifier columns.
#   v4: the refit call is passed to the worker wrapped in a list (a bare call
#       object passed through do.call() was evaluated when first used).
#   v3: the refit call carries the resampled data frame itself, not a symbol
#       (glmmTMB evaluates its model frame in environment(formula)).
#   v2: every replicate's response is checked against the family's support
#       before any refitting.
# =============================================================================

boot_mean_model <- function(deployed, trial_data, cluster, response,
                            stat = c("median", "mean"),
                            by = NULL,
                            strata = NULL,
                            response_transform = NULL,
                            newdata = NULL,
                            B = 2000L,
                            level = 0.95,
                            prediction_intervals = FALSE,
                            pi_draws = 10L,
                            parallel = FALSE,
                            seed = NULL,
                            back_transform = TRUE,
                            jitter_rows = NULL,
                            jitter_step = 1,
                            jitter_bw = 0.25,
                            recentre = c("mean", "median", "none"),
                            regrid = NULL,
                            interval = c("percentile", "bc"),
                            pi_dispersion = c("original", "replicate"),
                            max_fail = 0.05,
                            allow_fixed_weights = FALSE,
                            resample = c("clusters", "crossed"),
                            stimulus_strata = NULL,
                            verbose = TRUE) {
  
  # ---- deploy_mean_model.R must be loaded ------------------------------------
  internals <- c(".eval_eta", ".eval_eta_codes", ".apply_scaling", ".ensure_raw_columns", ".aligned_rows",
                 ".code_weights", ".contrast_matrix", ".build_column_map", ".reduce_columns",
                 ".expand_affine_pieces", ".unscale_reduced", ".format_piece_display", ".get_inverse_link")
  miss_int <- internals[!vapply(internals, exists, logical(1), mode = "function")]
  if (length(miss_int) > 0 || !exists(".DATA_SUMMARY_FNS")) {
    stop("Source deploy_mean_model.R (v17+) before using boot_mean_model(); missing: ",
         paste(c(miss_int, if (!exists(".DATA_SUMMARY_FNS")) ".DATA_SUMMARY_FNS"), collapse = ", "), call. = FALSE)
  }
  
  # ---- Inputs from the deployed object -------------------------------------
  if (!inherits(deployed, "deployed_mean_model")) stop("'deployed' must be an object returned by deploy_mean_model().", call. = FALSE)
  model <- attr(deployed, "model")
  if (is.null(model)) stop("The deployed object has no attr(, \"model\"): re-run deploy_mean_model() (v17 or later).", call. = FALSE)
  cls <- if (inherits(model, "glmmTMB")) "glmmTMB" else if (inherits(model, "glmgee")) "glmgee" else if (inherits(model, "brmsfit")) "brmsfit" else class(model)[1]
  if (cls == "brmsfit") {
    stop("brms models are not supported: refitting a Bayesian model thousands of times is impractical, and its posterior already provides uncertainty intervals.", call. = FALSE)
  }
  if (!cls %in% c("glmmTMB", "glmgee")) stop(sprintf("Models of class '%s' are not supported (glmmTMB, glmgee).", cls), call. = FALSE)
  if (cls == "glmmTMB" && length(tryCatch(glmmTMB::fixef(model)$zi, error = function(e) NULL)) > 0) {
    stop("Zero-inflated models are not supported: E[y] and the prediction draws would need the zero-inflation component.", call. = FALSE)
  }
  
  env <- environment(deployed)
  needed <- c("spec", "scaling", "numeric_codes", "constants_levels", "constants_model", "eta_with_constants")
  miss <- needed[!vapply(needed, exists, logical(1), envir = env, inherits = FALSE)]
  if (length(miss) > 0) {
    stop(sprintf("The deployed function's environment lacks: %s. Re-run deploy_mean_model() with the current version.", paste(miss, collapse = ", ")), call. = FALSE)
  }
  spec <- env$spec                       # frozen spec (data-dependent terms fixed)
  scaling <- env$scaling
  numeric_codes <- env$numeric_codes
  constants_levels <- env$constants_levels
  constants_model <- env$constants_model
  inv_link <- .get_inverse_link(spec$link)
  fam <- if (cls == "glmmTMB") spec$family_name else NA_character_
  ordb <- identical(fam, "ordbeta")
  psi0 <- if (ordb) spec$mean_info$psi else NULL
  cf0 <- spec$coefficients               # names as in the frozen spec; order = model order
  cf_names <- if (cls == "glmmTMB") names(glmmTMB::fixef(model)$cond) else rownames(model[["coefficients"]])
  if (length(cf_names) != length(cf0)) stop("Coefficient count differs between the model and the deployed spec.", call. = FALSE)
  
  stat_fn <- if (is.function(stat)) stat else switch(match.arg(stat), median = stats::median, mean = mean)
  stat_label <- if (is.function(stat)) "user function" else match.arg(stat)
  if (!is.null(response_transform) && !is.function(response_transform)) stop("response_transform must be a function or NULL.", call. = FALSE)
  .scalar_num <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!.scalar_num(B) || B < 2 || B != round(B)) stop("B must be a single whole number >= 2.", call. = FALSE)
  B <- as.integer(B)
  if (!.scalar_num(level) || !(level > 0 && level < 1)) stop("level must be a single number between 0 and 1.", call. = FALSE)
  if (!.scalar_num(max_fail) || max_fail < 0 || max_fail >= 1) stop("max_fail must be a single number in [0, 1).", call. = FALSE)
  if (!(is.numeric(jitter_step) && length(jitter_step) == 1L && jitter_step > 0)) stop("jitter_step must be a positive number.", call. = FALSE)
  if (!(is.numeric(jitter_bw) && length(jitter_bw) == 1L && jitter_bw > 0)) stop("jitter_bw must be a positive number.", call. = FALSE)
  if (isFALSE(regrid)) regrid <- NULL
  if (isTRUE(regrid)) {
    if (!identical(stat_label, "median")) {
      stop("regrid = TRUE (= jitter_step / 2, the grid of medians of step-spaced data) is only defined for stat = \"median\"; give a number instead.", call. = FALSE)
    }
    regrid <- jitter_step / 2
  }
  if (!is.null(regrid) && !(is.numeric(regrid) && length(regrid) == 1L && is.finite(regrid) && regrid > 0)) {
    stop("regrid must be NULL, TRUE or a positive number.", call. = FALSE)
  }
  interval <- match.arg(interval)
  resample <- match.arg(resample)
  pi_dispersion <- match.arg(pi_dispersion)
  recentre <- if (isTRUE(recentre)) "mean" else if (isFALSE(recentre)) "none" else match.arg(recentre)
  if (prediction_intervals && (!.scalar_num(pi_draws) || pi_draws < 1 || pi_draws != round(pi_draws))) {
    stop("pi_draws must be a single positive whole number.", call. = FALSE)
  }
  pi_draws <- as.integer(pi_draws)
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  
  # ---- Stimulus-level data and the response column ------------------------
  if (!isTRUE(spec$data_is_raw) || !is.data.frame(spec$data)) {
    stop("The model's stimulus-level data (with raw columns) could not be recovered by deploy_mean_model().", call. = FALSE)
  }
  stim_orig <- spec$data
  idx0 <- .aligned_rows(spec)             # rows of the ORIGINAL data used by the model (before renaming rows)
  stim <- .ensure_raw_columns(stim_orig, scaling)
  rownames(stim) <- NULL
  
  f_full <- if (cls == "glmmTMB") stats::formula(model) else stats::formula(model[["terms"]])
  target <- all.vars(f_full[[2]])
  if (length(target) != 1L) {
    stop(sprintf("The model's response (%s) uses %d data columns; only a single response column is supported.",
                 paste(deparse(f_full[[2]]), collapse = ""), length(target)), call. = FALSE)
  }
  if (!target %in% names(stim)) stop(sprintf("Response column '%s' is not in the model's stimulus-level data.", target), call. = FALSE)
  
  # ---- Trial-level data ----------------------------------------------------
  trial_data <- as.data.frame(trial_data)
  model_vars <- all.vars(stats::formula(spec$terms))
  if (is.null(by)) {
    by <- unique(vapply(model_vars, function(v) if (v %in% names(scaling)) scaling[[v]]$raw_name else v, character(1)))
    by_source <- "the model's predictor variables (raw names for centred/standardised ones)"
  } else {
    by_source <- "user-supplied"
  }
  need_tr <- unique(c(cluster, response, by, strata))
  miss_tr <- setdiff(need_tr, names(trial_data))
  if (length(miss_tr) > 0) stop(sprintf("trial_data is missing column(s): %s", paste(miss_tr, collapse = ", ")), call. = FALSE)
  miss_st <- setdiff(by, names(stim))
  if (length(miss_st) > 0) stop(sprintf("The stimulus-level data is missing 'by' column(s): %s", paste(miss_st, collapse = ", ")), call. = FALSE)
  
  na_st <- by[vapply(by, function(v) anyNA(stim[[v]]), logical(1))]
  if (length(na_st) > 0) stop(sprintf("The stimulus-level data have missing values in 'by' column(s): %s.", paste(na_st, collapse = ", ")), call. = FALSE)
  na_tr <- by[vapply(by, function(v) anyNA(trial_data[[v]]), logical(1))]
  if (length(na_tr) > 0) {
    stop(sprintf("trial_data has missing values in 'by' column(s): %s. Remove or fix those rows first.", paste(na_tr, collapse = ", ")), call. = FALSE)
  }
  if (!is.null(strata)) {
    na_s <- strata[vapply(strata, function(v) anyNA(trial_data[[v]]), logical(1))]
    if (length(na_s) > 0) stop(sprintf("trial_data has missing values in strata column(s): %s.", paste(na_s, collapse = ", ")), call. = FALSE)
  }
  stim_key <- .boot_key(stim, by)
  if (anyDuplicated(stim_key)) {
    stop(sprintf(paste0("The 'by' columns (%s; %s) do not identify the rows of the stimulus-level data uniquely ",
                        "(%d rows share their values with another row). Aggregating over them would merge rows the model ",
                        "treats as separate. Supply by = the column(s) that identify a stimulus row (e.g. a stimulus ID)."),
                 paste(by, collapse = ", "), by_source, sum(duplicated(stim_key) | duplicated(stim_key, fromLast = TRUE))), call. = FALSE)
  }
  tr_key <- .boot_key(trial_data, by)
  r_col <- trial_data[[response]]
  if (!(is.numeric(r_col) || is.logical(r_col))) {
    stop(sprintf("Trial-level response '%s' must be numeric or logical (it is %s); convert it explicitly first.",
                 response, paste(class(r_col), collapse = "/")), call. = FALSE)
  }
  resp <- as.numeric(r_col)
  if (all(is.na(resp))) stop(sprintf("Trial-level response '%s' has no numeric values.", response), call. = FALSE)
  
  # ---- Reproduce the model's response from the original trial data --------
  if (!is.null(response_transform)) {
    tt <- tryCatch(response_transform(c(0, 1, 2)), error = function(e) e)
    if (inherits(tt, "error") || !is.numeric(tt) || length(tt) != 3L) {
      stop("response_transform must be vectorised: given a numeric vector it must return a numeric vector of the same length.", call. = FALSE)
    }
  }
  y0 <- .boot_aggregate(resp, tr_key, seq_along(resp), stim_key, stat_fn, response_transform)
  yt <- suppressWarnings(as.numeric(stim[[target]]))
  same <- (is.na(y0) & is.na(yt)) | (!is.na(y0) & !is.na(yt) & abs(y0 - yt) <= 1e-8 * pmax(1, abs(yt)))
  if (!all(same)) {
    bad <- which(!same)
    no_trials <- sum(!stim_key %in% tr_key)
    ratio <- yt[bad] / y0[bad]
    ratio <- ratio[is.finite(ratio)]
    hint <- if (length(ratio) >= 3 && diff(range(ratio)) <= 1e-6 * max(1, abs(mean(ratio)))) {
      sprintf(" The model's column equals the aggregate times %.6g on these rows: pass response_transform = function(x) x * %.6g.", mean(ratio), mean(ratio))
    } else ""
    ex <- utils::head(bad, 5)
    stop(sprintf(paste0("Aggregating '%s' (%s per %s) over the original trial data does not reproduce the model's response column '%s' ",
                        "on %d of %d rows (e.g. rows %s: aggregate %s vs model data %s). %d stimulus rows have no matching trials.%s"),
                 response, stat_label, paste(by, collapse = " + "), target, length(bad), length(same),
                 paste(ex, collapse = ", "), paste(format(signif(y0[ex], 6)), collapse = ", "),
                 paste(format(signif(yt[ex], 6)), collapse = ", "), no_trials, hint), call. = FALSE)
  }
  
  # ---- Clusters and strata -------------------------------------------------
  cl_id <- as.character(trial_data[[cluster]])
  if (anyNA(cl_id)) stop(sprintf("Cluster column '%s' contains missing values.", cluster), call. = FALSE)
  rows_by_id <- split(seq_len(nrow(trial_data)), cl_id)
  ids <- names(rows_by_id)
  # trial rows / clusters that contribute to no stimulus row (silently ignored by
  # the aggregation, so reported)
  matched <- tr_key %in% stim_key
  n_unmatched_rows <- sum(!matched)
  has_match <- vapply(rows_by_id, function(r) any(matched[r]), logical(1))
  n_unmatched_clusters <- sum(!has_match)
  # clusters with no trial matching any stimulus row contribute nothing to any
  # aggregate; keeping them in the pool would dilute every replicate, so they are
  # dropped from the resampling population (and reported)
  ids <- ids[has_match[ids]]
  if (length(ids) < 2L) stop("Fewer than 2 clusters have trials matching the stimulus rows.", call. = FALSE)
  if (is.null(strata)) {
    st <- stats::setNames(rep("all", length(ids)), ids)
  } else {
    sv <- do.call(paste, c(lapply(strata, function(v) as.character(trial_data[[v]])), sep = "/"))
    st <- vapply(split(sv, cl_id)[ids], function(s) paste(sort(unique(s)), collapse = " + "), character(1))
  }
  strata_ids <- split(ids, st[ids])
  
  # ---- Stimulus resampling (resample = "crossed") ----------------------------
  stim_strata_rows <- NULL; sst_sizes <- NULL
  if (resample == "crossed") {
    pool <- if (is.null(idx0)) seq_len(nrow(stim)) else idx0      # stimulus rows the model used
    if (!is.null(stimulus_strata)) {
      miss_ss <- setdiff(stimulus_strata, names(stim))
      if (length(miss_ss) > 0) stop(sprintf("The stimulus-level data is missing stimulus_strata column(s): %s", paste(miss_ss, collapse = ", ")), call. = FALSE)
      na_ss <- stimulus_strata[vapply(stimulus_strata, function(v) anyNA(stim[[v]][pool]), logical(1))]
      if (length(na_ss) > 0) stop(sprintf("The stimulus-level data have missing values in stimulus_strata column(s): %s.", paste(na_ss, collapse = ", ")), call. = FALSE)
      sst <- do.call(paste, c(lapply(stimulus_strata, function(v) as.character(stim[[v]][pool])), sep = "/"))
    } else {
      sst <- rep("all", length(pool))
    }
    stim_strata_rows <- split(pool, sst)
    sst_sizes <- lengths(stim_strata_rows)
    if (any(sst_sizes < 2L)) {
      stop(sprintf("Crossed resampling needs at least 2 stimulus rows in every stimulus stratum (sizes: %s).",
                   paste(sprintf("%s: %d", names(sst_sizes), sst_sizes), collapse = "; ")), call. = FALSE)
    }
  } else if (!is.null(stimulus_strata)) {
    warning("stimulus_strata is ignored unless resample = \"crossed\".", call. = FALSE)
  }
  
  # ---- Design matrices (deployed path) -------------------------------------
  data_train <- .apply_scaling(stim, scaling)
  X_train <- .boot_design(spec, data_train, list())
  dev_train <- NA_real_
  if (!is.null(idx0) && is.numeric(spec$fitted)) {
    dev_train <- max(abs(as.numeric(X_train[idx0, , drop = FALSE] %*% cf0) - spec$fitted))
    if (dev_train > 1e-7 * max(1, max(abs(spec$fitted)))) {
      stop(sprintf("The training design matrix does not reproduce the model's own linear predictor (max|dev| = %.3g).", dev_train), call. = FALSE)
    }
  }
  
  if (is.null(newdata)) newdata <- stim
  newdata <- as.data.frame(newdata)
  if (nrow(newdata) == 0L) stop("newdata has no rows.", call. = FALSE)
  miss_nd <- setdiff(attr(deployed, "free_vars"), names(newdata))
  if (length(miss_nd) > 0) stop(sprintf("newdata is missing required column(s): %s", paste(miss_nd, collapse = ", ")), call. = FALSE)
  data_new <- .apply_scaling(newdata, scaling)
  for (v in names(constants_levels)) {
    data_new[[v]] <- if (!is.null(spec$xlevels[[v]])) as.character(constants_levels[[v]]) else constants_levels[[v]]
  }
  X_new <- .boot_design(spec, data_new, numeric_codes)
  eta0 <- as.numeric(X_new %*% cf0)
  eta_dep <- env$eta_with_constants(data_new)
  dev_new <- max(abs(eta0 - eta_dep), na.rm = TRUE)
  if (dev_new > 1e-7 * max(1, max(abs(eta_dep), na.rm = TRUE))) {
    stop(sprintf("The prediction design matrix does not reproduce the deployed linear predictor (max|dev| = %.3g).", dev_new), call. = FALSE)
  }
  
  # ---- Dispersion design (prediction intervals only) -----------------------
  X_disp <- NULL; disp0 <- NULL
  if (prediction_intervals) {
    if (cls != "glmmTMB") stop("Prediction intervals need a full response distribution: glmmTMB models only.", call. = FALSE)
    if (isTRUE(spec$has_re)) stop("Prediction intervals are not implemented for models with random effects.", call. = FALSE)
    if (!fam %in% c("ordbeta", "beta", "gaussian")) {
      stop(sprintf("Prediction intervals are implemented for the ordbeta, beta and gaussian families, not '%s'.", fam), call. = FALSE)
    }
    dform <- tryCatch(model[["modelInfo"]][["allForm"]][["dispformula"]], error = function(e) NULL)
    if (is.null(dform)) dform <- ~1
    if (any(all.names(dform) %in% .DATA_SUMMARY_FNS)) stop("The dispersion formula contains a data-dependent term (e.g. mean(x)); not supported.", call. = FALSE)
    disp0 <- tryCatch(glmmTMB::fixef(model)$disp, error = function(e) NULL)
    if (is.null(disp0) || length(disp0) == 0) stop("Could not read the dispersion coefficients (fixef(model)$disp).", call. = FALSE)
    X_disp <- .boot_disp_design(dform, data_new, stim, spec, numeric_codes)
    if (!identical(colnames(X_disp), names(disp0))) {
      stop(sprintf("Dispersion design columns (%s) do not match the model's dispersion coefficients (%s).",
                   paste(colnames(X_disp), collapse = ", "), paste(names(disp0), collapse = ", ")), call. = FALSE)
    }
  }
  
  # ---- Refit call ------------------------------------------------------------
  fit_call <- stats::getCall(model)
  if (is.null(fit_call)) stop("Could not retrieve the model's call.", call. = FALSE)
  has_offset_term <- tryCatch(!is.null(attr(stats::terms(f_full), "offset")), error = function(e) FALSE)
  if (!is.null(fit_call$offset) || has_offset_term || any(grepl("offset(", deparse(f_full), fixed = TRUE))) {
    stop("Models with an offset are not supported (the deployed linear predictor has no offset).", call. = FALSE)
  }
  has_weights <- !is.null(fit_call$weights)
  if (has_weights && !isTRUE(allow_fixed_weights)) {
    stop(paste0("The model has weights. They would be kept from the stimulus-level data, NOT recomputed for each replicate, ",
                "which is only right for fixed design weights (not, e.g., numbers of responses or inverse variances). ",
                "If they are fixed design weights, set allow_fixed_weights = TRUE."), call. = FALSE)
  }
  fit_call$data <- quote(.boot_data)
  fit_call[[1]] <- if (cls == "glmmTMB") quote(glmmTMB::glmmTMB) else quote(glmtoolbox::glmgee)
  env_eval <- tryCatch(environment(f_full), error = function(e) NULL)
  if (!is.environment(env_eval)) env_eval <- globalenv()
  free <- setdiff(unique(unlist(lapply(as.list(fit_call)[-1], all.vars))), c(names(stim), ".boot_data"))
  ext <- list()
  for (v in free) if (exists(v, envir = env_eval)) ext[[v]] <- get(v, envir = env_eval)
  pkgs <- if (cls == "glmmTMB") "glmmTMB" else "glmtoolbox"
  
  # ---- Resamples (all random numbers drawn here, in the main session) -------
  if (!is.null(seed)) {
    seed_old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
    on.exit({
      if (!is.null(seed_old)) assign(".Random.seed", seed_old, envir = globalenv())
      else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
    }, add = TRUE)
    set.seed(seed)
  }
  samples <- lapply(seq_len(B), function(b) {
    unlist(lapply(strata_ids, function(v) v[sample.int(length(v), length(v), replace = TRUE)]), use.names = FALSE)
  })
  stim_samples <- if (resample == "crossed") lapply(seq_len(B), function(b) {
    unlist(lapply(stim_strata_rows, function(v) v[sample.int(length(v), length(v), replace = TRUE)]), use.names = FALSE)
  }) else NULL
  
  # ---- Replicate responses (all built here, before any refitting) ---------
  # Optional smoothed bootstrap for responses on a discrete grid (e.g. integer
  # ratings), mirroring the measurement-data procedure: in each replicate the
  # resampled values of the marked trials are shrunk towards their stimulus's
  # ORIGINAL aggregate (factor 1 / sqrt(1 + sigma^2 / var), offsetting the
  # variance the jitter adds) and Gaussian jitter (sigma = jitter_bw * jitter_step)
  # is added, so the aggregate moves smoothly instead of sticking and jumping.
  # With recentre = "mean" (default) or "median", each stimulus containing
  # jittered trials then has its replicate aggregates shifted so that their mean
  # (or median) equals its original aggregate (removes the shift that jittering
  # introduces). Optional regrid, then response_transform (e.g. clipping), are
  # applied afterwards. Point estimates are not affected.
  stim_f <- factor(tr_key, levels = stim_key)
  raw0 <- .boot_aggregate(resp, tr_key, seq_along(resp), stim_key, stat_fn, NULL)   # original aggregates, untransformed
  jit <- .boot_jitter_rows(jitter_rows, trial_data)
  smooth_info <- NULL
  if (!is.null(jit)) {
    sigma <- jitter_bw * jitter_step
    row_stim <- match(tr_key, stim_key)
    v_stim <- vapply(split(resp[jit], stim_f[jit]), function(v) { v <- v[!is.na(v)]; if (length(v) > 1L) stats::var(v) else NA_real_ }, numeric(1))
    k_stim <- ifelse(is.finite(v_stim) & v_stim > 0, 1 / sqrt(1 + sigma^2 / v_stim), 1)
    m_row <- raw0[row_stim]
    k_row <- unname(k_stim[row_stim])
    stim_jit <- vapply(split(jit, stim_f), any, logical(1))
    stim_all <- vapply(split(jit, stim_f), function(v) length(v) == 0L || all(v), logical(1))
    if (any(stim_jit & !stim_all)) {
      stop(sprintf(paste0("jitter_rows marks only some of the trials of %d stimulus row(s); smoothing must apply to all or none of ",
                          "a stimulus's trials (the centre is the aggregate of all of them)."), sum(stim_jit & !stim_all)), call. = FALSE)
    }
    smooth_info <- list(n_rows = sum(jit), n_stimuli = sum(stim_jit), step = jitter_step, bw = jitter_bw,
                        sigma = sigma, recentre = recentre)
  }
  Yraw <- vapply(samples, function(s) {
    idx <- unlist(rows_by_id[s], use.names = FALSE)
    x <- resp[idx]
    if (!is.null(jit)) {
      j <- jit[idx] & !is.na(x) & !is.na(m_row[idx])
      x[j] <- (x[j] - m_row[idx][j]) * k_row[idx][j] + m_row[idx][j] + stats::rnorm(sum(j), 0, sigma)
    }
    .boot_aggregate(x, tr_key[idx], seq_along(x), stim_key, stat_fn, NULL)
  }, numeric(length(stim_key)))
  Yraw <- matrix(Yraw, nrow = length(stim_key))
  if (!is.null(jit) && recentre != "none") {
    centre_fn <- if (recentre == "median") stats::median else mean
    shift <- apply(Yraw, 1, centre_fn, na.rm = TRUE) - raw0
    shift[!stim_jit | !is.finite(shift)] <- 0
    Yraw <- Yraw - shift
    smooth_info$shift <- stats::setNames(shift, stim_key)
  }
  if (is.null(jit) && !is.null(regrid)) warning("regrid is ignored because no trials are smoothed (jitter_rows is NULL or marks none).", call. = FALSE)
  if (!is.null(jit) && !is.null(regrid)) {
    # the grid must be one the ORIGINAL aggregates of the smoothed stimuli lie on,
    # otherwise rounding would move replicates onto values real aggregates cannot take
    on_grid <- function(v, g) { v <- v[is.finite(v)]; all(abs(v / g - round(v / g)) < 1e-8) }
    orig_j <- raw0[stim_jit]
    if (!on_grid(orig_j, regrid)) {
      fits <- c(jitter_step, jitter_step / 2)[vapply(c(jitter_step, jitter_step / 2), function(g) on_grid(orig_j, g), logical(1))]
      off <- orig_j[is.finite(orig_j) & abs(orig_j / regrid - round(orig_j / regrid)) >= 1e-8]
      stop(sprintf("regrid = %g: the original aggregates of the smoothed stimuli are not all multiples of %g (e.g. %s). %s",
                   regrid, regrid, paste(format(utils::head(off, 3)), collapse = ", "),
                   if (length(fits)) sprintf("They do lie on multiples of %g: use that.", fits[length(fits)]) else
                     "They are not on the jitter_step grid or its half either: check jitter_step, jitter_rows and stat."), call. = FALSE)
    }
    if (regrid < jitter_step / 2 - 1e-12) {
      warning(sprintf("regrid = %g is finer than jitter_step / 2 = %g, the finest grid medians of step-spaced data can take.", regrid, jitter_step / 2), call. = FALSE)
    }
    if (abs(regrid - jitter_step / 2) < 1e-12 && on_grid(orig_j, jitter_step) && verbose) {
      cat(sprintf("Note: all original aggregates of the smoothed stimuli are multiples of jitter_step = %g; regrid = %g would match the data more exactly.\n",
                  jitter_step, jitter_step))
    }
    # Mean-preserving stochastic rounding of the smoothed stimuli's replicate
    # aggregates back onto the grid the real aggregates take (e.g. 0.5 for medians
    # of integer ratings): v -> floor(v/g)*g, or one step up with probability equal
    # to the fractional position. Restores exact grid values (including exact 0s
    # after clipping), so the refitted model sees data shaped like the original.
    rows <- which(stim_jit)
    V <- Yraw[rows, , drop = FALSE]
    lo <- floor(V / regrid) * regrid
    up <- matrix(stats::runif(length(V)), nrow(V)) < (V - lo) / regrid
    V <- lo + regrid * up
    Yraw[rows, ] <- V
    smooth_info$regrid <- regrid
  }
  Y <- if (is.null(response_transform)) Yraw else matrix(as.numeric(response_transform(as.vector(Yraw))), nrow = nrow(Yraw))
  
  # ---- Pre-check: every replicate's response must be valid for the family ----
  # The reproduction check above cannot see a transform that is a no-op on the
  # original data but matters for resampled panels (e.g. clipping at 0); such
  # replicates would fail to refit, and discarding them would bias the
  # intervals. So all replicates are checked before any refitting, on the
  # response as the model sees it (formula LHS evaluated).
  fam_check <- if (cls == "glmmTMB") fam else tryCatch(model[["family"]][["family"]], error = function(e) NA_character_)
  lhs <- f_full[[2]]
  pre <- vapply(seq_len(B), function(b) {
    sb <- stim
    sb[[target]] <- Y[, b]
    v <- suppressWarnings(as.numeric(eval(lhs, sb, env_eval)))
    inside <- .boot_in_support(fam_check, v)
    c(n_out = sum(!is.na(v) & !inside), lo = suppressWarnings(min(v, na.rm = TRUE)), hi = suppressWarnings(max(v, na.rm = TRUE)),
      n_na = sum(is.na(v)))
  }, numeric(4))
  if (any(pre["n_out", ] > 0)) {
    k <- sum(pre["n_out", ] > 0)
    stop(sprintf(paste0("In %d of %d replicates the re-aggregated response (as the model sees it, %s) falls outside what the '%s' ",
                        "family accepts (%s); across replicates it ranges from %.4g to %.4g. Those refits would fail, and ",
                        "discarding them would bias the intervals. The original data are within range, so the transform that ",
                        "defines the response must also hold for resampled panels, e.g. response_transform = function(x) ",
                        "pmax(x, 0) / 10 for a non-negative quantity that is the aggregate / 10."),
                 k, B, paste(deparse(lhs), collapse = ""), fam_check, .boot_support_label(fam_check),
                 min(pre["lo", ]), max(pre["hi", ])), call. = FALSE)
  }
  Y_list <- lapply(seq_len(B), function(b) list(y = Y[, b], rows = if (is.null(stim_samples)) NULL else stim_samples[[b]]))
  
  worker_args <- list(stim = stim, target = target, fit_spec = list(call = fit_call), ext = ext, cls = cls, fam = fam,
                      cf_names = cf_names, X_train = X_train, want_disp = prediction_intervals)
  t0 <- Sys.time()
  if (verbose) cat(sprintf("Refitting %d bootstrap replicates (%s)...\n", B, if (parallel) "future.apply" else "serial"))
  res <- if (parallel) {
    if (!requireNamespace("future.apply", quietly = TRUE)) stop("parallel = TRUE needs the future.apply package.", call. = FALSE)
    do.call(future.apply::future_lapply, c(list(X = Y_list, FUN = .boot_refit_one), worker_args,
                                           list(future.packages = pkgs, future.seed = FALSE)))
  } else {
    do.call(lapply, c(list(X = Y_list, FUN = .boot_refit_one), worker_args))
  }
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  # ---- Collect -----------------------------------------------------------------
  ok <- vapply(res, function(r) isTRUE(r$ok), logical(1))
  n_ok <- sum(ok)
  reasons <- vapply(res[!ok], function(r) if (is.null(r$reason) || is.na(r$reason)) "unknown" else r$reason, character(1))
  warn_tab <- sort(table(unlist(lapply(res, `[[`, "warnings"))), decreasing = TRUE)
  if (n_ok < 2L) stop(sprintf("Only %d of %d replicates gave a usable fit. First failure: %s", n_ok, B, reasons[1]), call. = FALSE)
  fail_rate <- (B - n_ok) / B
  if (fail_rate > 4 * max_fail) {
    stop(sprintf(paste0("%d of %d replicates (%.0f%%) failed, more than 4 x max_fail = %.0f%%. The surviving replicates may not ",
                        "represent the bootstrap distribution (failures can depend on the resampled data). Most common reason: %s"),
                 B - n_ok, B, 100 * fail_rate, 400 * max_fail, names(sort(table(reasons), decreasing = TRUE))[1]), call. = FALSE)
  }
  if (fail_rate > max_fail) {
    warning(sprintf("%d of %d replicates (%.1f%%) failed, more than max_fail = %.1f%%; the intervals rest on the surviving replicates. See $failures.",
                    B - n_ok, B, 100 * fail_rate, 100 * max_fail), call. = FALSE)
  }
  if (n_ok < B && verbose) cat(sprintf("%d of %d replicates discarded (see $failures).\n", B - n_ok, B))
  CF <- do.call(rbind, lapply(res[ok], `[[`, "cf"))
  colnames(CF) <- cf_names
  PSI <- if (ordb) do.call(rbind, lapply(res[ok], `[[`, "psi")) else NULL
  DISP <- if (prediction_intervals) do.call(rbind, lapply(res[ok], `[[`, "disp")) else NULL
  n_rows <- vapply(res, function(r) if (is.null(r$n_rows)) NA_real_ else as.numeric(r$n_rows), numeric(1))
  design_dev <- vapply(res[ok], function(r) if (is.null(r$design_dev)) NA_real_ else as.numeric(r$design_dev), numeric(1))
  n_no_pd <- if (cls == "glmmTMB") sum(!vapply(res[ok], function(r) isTRUE(r$pdHess), logical(1))) else 0L
  
  # ---- Predictions: E[y] for newdata -------------------------------------------
  ETA <- X_new %*% t(CF)                                     # rows = newdata, cols = replicates
  E <- .boot_mean_map(ETA, inv_link, PSI)
  E0 <- as.numeric(.boot_mean_map(matrix(eta0, ncol = 1), inv_link, if (ordb) matrix(psi0, nrow = 1) else NULL))
  # ---- Back-transform to the scale of the aggregated trial response --------
  # The model's response is g(aggregate) = formula LHS evaluated on
  # response_transform(aggregate). Its inverse is applied to the predictions
  # (estimate, limits, PI limits; bias recomputed from the replicates), so they
  # are on the scale of the median/mean of the trial-level 'response' column.
  bt_info <- if (isTRUE(back_transform)) .boot_back_transform(resp, target, lhs, response_transform, env_eval) else
    list(type = "none", note = "back_transform = FALSE: predictions are on the model's response scale.")
  inv <- if (is.function(bt_info$inverse)) bt_info$inverse else identity
  E_out <- matrix(inv(as.vector(E)), nrow = nrow(E))
  E0_out <- inv(E0)
  n_bad_inv <- sum(is.finite(E) & !is.finite(E_out)) + sum(is.finite(E0) & !is.finite(E0_out))
  if (n_bad_inv > 0) {
    stop(sprintf("The back-transform gave %d non-finite value(s) from finite model-scale predictions; use back_transform = FALSE.", n_bad_inv), call. = FALSE)
  }
  # limits from the back-transformed replicates themselves, so they are already
  # ordered whatever the direction of the transform (no swap needed here)
  ci_E <- .boot_limits(E_out, E0_out, probs, interval)
  pred <- data.frame(estimate = E0_out, bias = rowMeans(E_out, na.rm = TRUE) - E0_out,
                     p_below = .boot_p_below(E_out, E0_out),
                     lower = ci_E[, 1], upper = ci_E[, 2])
  pred_se <- apply(E_out, 1, stats::sd, na.rm = TRUE)
  
  if (prediction_intervals) {
    LPHI <- if (pi_dispersion == "original") {
      matrix(as.numeric(X_disp %*% disp0), nrow = nrow(X_disp), ncol = ncol(ETA))
    } else {
      X_disp %*% t(DISP)
    }
    draws <- matrix(NA_real_, nrow(ETA), ncol(ETA) * pi_draws)
    for (b in seq_len(ncol(ETA))) {
      cols <- (b - 1L) * pi_draws + seq_len(pi_draws)
      draws[, cols] <- .boot_draw(fam, ETA[, b], exp(LPHI[, b]), if (ordb) PSI[b, ] else NULL, inv_link, pi_draws)
    }
    ci_P <- matrix(inv(as.vector(.boot_q(draws, probs))), ncol = 2)
    if (isTRUE(bt_info$decreasing)) ci_P <- ci_P[, 2:1, drop = FALSE]
    pred$pi_lower <- ci_P[, 1]
    pred$pi_upper <- ci_P[, 2]
  }
  id_cols <- unique(c(intersect(by, names(newdata)), attr(deployed, "free_vars")))
  pred <- cbind(newdata[, id_cols, drop = FALSE], pred)
  
  # ---- Coefficient tables ----------------------------------------------------------
  est_model <- cf0; names(est_model) <- cf_names
  M_model <- CF
  if (ordb) {
    est_model <- c(est_model, psi1 = psi0[1], psi2 = psi0[2])
    M_model <- cbind(M_model, psi1 = PSI[, 1], psi2 = PSI[, 2])
  }
  if (prediction_intervals) {
    est_model <- c(est_model, stats::setNames(disp0, paste0("disp:", names(disp0))))
    D <- DISP
    colnames(D) <- paste0("disp:", names(disp0))
    M_model <- cbind(M_model, D)
  }
  coef_model <- .boot_coef_table(est_model, M_model, probs, interval)
  coef_model$note <- ""
  if (ordb) {
    psi_fixed <- spec$mean_info$psi_status == "fixed via map"
    coef_model$note[coef_model$term %in% c("psi1", "psi2")[psi_fixed]] <- "fixed via map"
    yl <- suppressWarnings(as.numeric(eval(f_full[[2]], stim, env_eval)))
    if (!psi_fixed[1] && sum(yl == 0, na.rm = TRUE) == 0) coef_model$note[coef_model$term == "psi1"] <- "not identified (no 0s): interval meaningless"
    if (!psi_fixed[2] && sum(yl == 1, na.rm = TRUE) == 0) coef_model$note[coef_model$term == "psi2"] <- "not identified (no 1s): interval meaningless"
  }
  off_centre <- with(coef_model, is.finite(p_below) & (p_below < 0.4 | p_below > 0.6))
  coef_model$note[off_centre & coef_model$note == ""] <- "off-centre (p_below outside 0.4-0.6)"
  
  coef_equation <- NULL; EQ <- NULL
  eq <- tryCatch({
    cmap <- .build_column_map(spec$terms, names(cf0), spec$data, spec$contrasts, spec$xlevels)
    eqf <- function(cf) .boot_equation_coefs(stats::setNames(as.numeric(cf), names(cf0)), cmap, constants_model, spec$data, scaling)
    e0 <- eqf(cf0)
    reps <- lapply(seq_len(nrow(CF)), function(i) eqf(CF[i, ]))
    keys <- unique(c(names(e0$coef), unlist(lapply(reps, function(r) names(r$coef)))))
    labs <- e0$label
    for (r in reps) labs <- c(labs, r$label[setdiff(names(r$label), names(labs))])
    M <- t(vapply(reps, function(r) { v <- r$coef[keys]; v[is.na(v)] <- 0; v }, numeric(length(keys))))
    if (length(keys) == 1L) M <- matrix(M, ncol = 1)
    est <- e0$coef[keys]; est[is.na(est)] <- 0
    tab <- .boot_coef_table(stats::setNames(est, labs[keys]), M, probs, interval)
    list(tab = tab, M = stats::setNames(as.data.frame(M), labs[keys]))
  }, error = function(e) e)
  if (inherits(eq, "error")) {
    if (verbose) cat(sprintf("Deployed-equation coefficient intervals unavailable: %s\n", conditionMessage(eq)))
  } else {
    coef_equation <- eq$tab
    EQ <- as.matrix(eq$M)
  }
  
  out <- list(
    predictions = pred,
    coef_model = coef_model,
    coef_equation = coef_equation,
    replicates = list(coef = CF, psi = PSI, disp = DISP, eta = ETA, E = E, equation = EQ,
                      response_raw = Yraw[, ok, drop = FALSE], response_original = stats::setNames(raw0, stim_key)),
    failures = list(n_failed = B - n_ok, reasons = sort(table(reasons), decreasing = TRUE),
                    warnings = warn_tab, n_without_pd_hessian = n_no_pd),
    settings = list(B = B, B_ok = n_ok, level = level, interval = interval, pi_dispersion = pi_dispersion, stat = stat_label,
                    cluster = cluster, n_clusters = length(ids), strata = strata,
                    strata_sizes = table(st), by = by, by_source = by_source, response = response,
                    target = target, family = if (cls == "glmmTMB") fam else "glmgee",
                    constants = attr(deployed, "constants"), back_transform = bt_info, smoothing = smooth_info,
                    prediction_intervals = prediction_intervals, pi_draws = pi_draws,
                    rows_per_replicate = range(n_rows, na.rm = TRUE), seed = seed,
                    parallel = parallel, elapsed_sec = elapsed,
                    unmatched_trial_rows = n_unmatched_rows, unmatched_clusters = n_unmatched_clusters,
                    has_weights = has_weights, resample = resample, stimulus_strata = stimulus_strata,
                    stimulus_strata_sizes = sst_sizes,
                    pred_off_centre = sum(is.finite(pred$p_below) & (pred$p_below < 0.4 | pred$p_below > 0.6))),
    checks = list(train_design_dev = dev_train, newdata_design_dev = dev_new,
                  replicate_design_dev_max = if (all(is.na(design_dev))) NA_real_ else max(design_dev, na.rm = TRUE),
                  replicates_design_unchecked = sum(is.na(design_dev)))
  )
  class(out) <- "boot_mean_model"
  if (verbose) print(out)
  invisible(out)
}

# ---- Worker: one replicate ----------------------------------------------------
# Self-contained (base R plus the fitting package), so it can run on future
# workers without the rest of this file.
.boot_refit_one <- function(item, stim, target, fit_spec, ext, cls, fam, cf_names, X_train, want_disp) {
  # the call arrives wrapped in a list: a bare call object passed through do.call()
  # would be EVALUATED when the argument is first used (the v3 bug)
  fit_call <- fit_spec$call
  out <- list(ok = FALSE, reason = NA_character_, warnings = character(0))
  # item: the replicate's response for every stimulus row, and (crossed
  # resampling) the stimulus rows drawn; row_map[i] = original row of row i
  stim[[target]] <- item$y
  row_map <- seq_len(nrow(stim))
  if (!is.null(item$rows)) {
    stim <- stim[item$rows, , drop = FALSE]
    row_map <- item$rows
  }
  rownames(stim) <- NULL
  out$n_rows <- sum(!is.na(stim[[target]]))
  warns <- character(0)
  # The data frame itself goes into the call: glmmTMB evaluates its model frame
  # in environment(formula), which is usually NOT where the call is evaluated,
  # so a data argument given as a symbol might not be found there.
  fit_call$data <- stim
  fenv <- list2env(ext, parent = globalenv())
  fit <- withCallingHandlers(
    tryCatch(eval(fit_call, envir = fenv), error = function(e) e),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  out$warnings <- unique(warns)
  if (inherits(fit, "error")) { out$reason <- paste("fit failed:", conditionMessage(fit)); return(out) }
  
  eta_fit <- NULL; X_fit <- NULL
  if (cls == "glmmTMB") {
    conv <- fit[["fit"]][["convergence"]]
    if (!isTRUE(conv == 0)) { out$reason <- sprintf("optimiser did not converge (code %s)", format(conv)); return(out) }
    out$pdHess <- isTRUE(fit[["sdr"]][["pdHess"]])
    cf <- glmmTMB::fixef(fit)$cond
    used <- suppressWarnings(as.integer(rownames(fit[["frame"]])))
    X_fit <- tryCatch(as.matrix(glmmTMB::getME(fit, "X")), error = function(e) NULL)
    if (is.null(X_fit)) X_fit <- tryCatch(as.matrix(fit[["obj"]][["env"]][["data"]][["X"]]), error = function(e) NULL)
    if (identical(fam, "ordbeta")) {
      psi <- tryCatch(unname(fit[["obj"]][["env"]][["parList"]](fit[["fit"]][["par"]], fit[["fit"]][["parfull"]])$psi),
                      error = function(e) NULL)
      if (length(psi) != 2L || !all(is.finite(psi))) { out$reason <- "could not read the ordbeta cutpoints"; return(out) }
      out$psi <- psi
    }
    if (want_disp) {
      d <- tryCatch(glmmTMB::fixef(fit)$disp, error = function(e) NULL)
      if (is.null(d) || !all(is.finite(d))) { out$reason <- "could not read the dispersion coefficients"; return(out) }
      out$disp <- d
    }
  } else {
    cf <- stats::setNames(as.numeric(fit[["coefficients"]]), rownames(fit[["coefficients"]]))
    used <- suppressWarnings(as.integer(rownames(fit[["model"]])))
    eta_fit <- as.numeric(fit[["linear.predictors"]])
    # convergence flag, if this glmtoolbox version stores one (not verified)
    conv_g <- fit[["converged"]]
    if (!is.null(conv_g) && !isTRUE(as.logical(conv_g))) { out$reason <- "glmgee did not converge"; return(out) }
  }
  if (!identical(names(cf), cf_names)) { out$reason <- "coefficient names differ from the original model"; return(out) }
  if (!all(is.finite(cf))) { out$reason <- "non-finite coefficients"; return(out) }
  
  # the replicate's own design (glmmTMB) or linear predictor (glmgee) must match the precomputed one
  dev <- NA_real_
  if (!anyNA(used) && length(used) > 0 && max(used) <= length(row_map)) {
    Xt <- X_train[row_map[used], , drop = FALSE]
    if (!is.null(X_fit) && identical(dim(X_fit), dim(Xt))) {
      dev <- max(abs(X_fit - Xt))
      tol <- 1e-8 * max(1, max(abs(Xt)))
    } else if (!is.null(eta_fit) && length(eta_fit) == nrow(Xt)) {
      dev <- max(abs(eta_fit - as.numeric(Xt %*% cf)))
      tol <- 1e-7 * max(1, max(abs(eta_fit)))
    }
    if (is.finite(dev) && dev > tol) {
      out$reason <- sprintf("the refit's design does not match the deployed design (max|dev| = %.3g)", dev)
      return(out)
    }
  }
  out$design_dev <- dev
  out$cf <- unname(cf)
  out$ok <- TRUE
  out
}

# ---- Helpers ------------------------------------------------------------------
.boot_key <- function(d, by) {
  do.call(paste, c(lapply(by, function(v) {
    x <- d[[v]]
    if (is.numeric(x)) formatC(x, digits = 12, format = "g") else as.character(x)
  }), sep = "\r"))
}

# stimulus-level statistic of the trial rows idx, aligned with the stimulus rows
.boot_aggregate <- function(resp, key, idx, stim_key, stat_fn, transform) {
  k <- factor(key[idx], levels = unique(stim_key))
  agg <- vapply(split(resp[idx], k), function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0L) NA_real_ else as.numeric(stat_fn(x))
  }, numeric(1))
  y <- unname(agg[stim_key])
  if (!is.null(transform)) y <- transform(y)
  as.numeric(y)
}

# Which trial rows get the smoothed (jittered) bootstrap: NULL/FALSE = none,
# TRUE = all, a logical vector (one per trial row), or the name of a logical
# column of trial_data. NA counts as FALSE.
.boot_jitter_rows <- function(jitter_rows, trial_data) {
  n <- nrow(trial_data)
  if (is.null(jitter_rows) || identical(jitter_rows, FALSE)) return(NULL)
  if (identical(jitter_rows, TRUE)) return(rep(TRUE, n))
  if (is.character(jitter_rows) && length(jitter_rows) == 1L) {
    if (!jitter_rows %in% names(trial_data)) stop(sprintf("jitter_rows: column '%s' is not in trial_data.", jitter_rows), call. = FALSE)
    jitter_rows <- trial_data[[jitter_rows]]
  }
  if (!is.logical(jitter_rows) || length(jitter_rows) != n) {
    stop("jitter_rows must be NULL, TRUE, a logical vector with one value per trial row, or the name of a logical column.", call. = FALSE)
  }
  out <- jitter_rows & !is.na(jitter_rows)
  if (!any(out)) return(NULL)
  out
}

# Inverse of g(x) = formula LHS evaluated on response_transform(x), the map from
# an aggregate of the trial response to the model's response. g is evaluated on
# a fine grid over the trial response's range. Flat ends (clipping, e.g.
# pmax(x, 0)) are allowed and reported; the strictly monotone part between them
# is tested for linearity. Affine -> exact inverse; monotone non-linear ->
# numerical inverse, flagged (the back-transformed mean is not the mean on the
# original scale); not monotone -> no back-transform.
.boot_back_transform <- function(resp, target, lhs, transform, env_eval, n_grid = 2001L) {
  none <- function(msg) list(type = "none", inverse = NULL, decreasing = FALSE,
                             note = paste0(msg, " Predictions are left on the model's response scale."))
  if (is.null(transform) && identical(lhs, as.name(target))) {
    return(list(type = "identity", inverse = identity, decreasing = FALSE,
                note = "none needed (the model's response is the aggregate itself)."))
  }
  r <- suppressWarnings(range(resp[is.finite(resp)]))
  if (!all(is.finite(r)) || r[1] == r[2]) return(none("The trial response has no usable range."))
  xs <- seq(r[1], r[2], length.out = n_grid)
  d <- list()
  d[[target]] <- if (is.null(transform)) xs else transform(xs)
  gs <- tryCatch(suppressWarnings(as.numeric(eval(lhs, d, env_eval))), error = function(e) rep(NA_real_, length(xs)))
  keep <- is.finite(gs)
  xs <- xs[keep]; gs <- gs[keep]
  if (length(gs) < 3L || length(unique(gs)) < 3L) return(none("The response transform could not be evaluated over the trial response's range, or is constant."))
  
  decreasing <- gs[length(gs)] < gs[1]
  lo_val <- if (decreasing) max(gs) else min(gs)
  hi_val <- if (decreasing) min(gs) else max(gs)
  i1 <- max(which(gs == lo_val))            # end of a flat (clipped) start, if any
  i2 <- min(which(gs == hi_val))            # start of a flat (clipped) end, if any
  # the kink of a clip lies BETWEEN grid points, so the last flat point (i1 > 1)
  # and the first flat point of a flat end (i2 < n) are left out of the
  # linearity test; the line through the remaining points still meets the
  # clipped value exactly at the clip point
  c1 <- if (i1 > 1L) i1 + 1L else i1
  c2 <- if (i2 < length(gs)) i2 - 1L else i2
  if (c2 - c1 < 2L) return(none("The response transform has no monotone range to invert."))
  core <- c1:c2
  dd <- diff(gs[core])
  if (!(all(dd > 0) || all(dd < 0))) {
    return(none("The response transform is NOT monotone over the trial response's range, so it cannot be inverted."))
  }
  clip <- character(0)
  if (i1 > 1L) clip <- c(clip, sprintf("flat (clipped) for aggregates below about %.4g", xs[i1]))
  if (i2 < length(xs)) clip <- c(clip, sprintf("flat (clipped) for aggregates above about %.4g", xs[i2]))
  clip_txt <- if (length(clip)) paste0("; the transform is ", paste(clip, collapse = " and "),
                                       ", so model values at a clipped bound map back to that bound") else ""
  
  slope <- (gs[c2] - gs[c1]) / (xs[c2] - xs[c1])
  icpt <- gs[c1] - slope * xs[c1]
  resid <- max(abs(gs[core] - (icpt + slope * xs[core])))
  if (resid <= 1e-9 * max(1, max(abs(gs)))) {
    a <- -icpt / slope; b <- 1 / slope
    return(list(type = "affine", inverse = function(y) a + b * y, intercept = a, slope = b, decreasing = decreasing,
                note = sprintf("affine inverse, aggregate = %.6g + %.6g * (model value)%s. Expected values carry over exactly.",
                               a, b, clip_txt)))
  }
  xc <- xs[i1:i2]; gc <- gs[i1:i2]
  inverse <- function(y) stats::approx(gc, xc, xout = y, rule = 2, ties = mean)$y
  list(type = "non-linear", inverse = inverse, decreasing = decreasing,
       note = paste0("NON-LINEAR monotone inverse (numerical, ", n_grid, "-point grid)", clip_txt,
                     ". Interval limits back-transform exactly (quantiles commute), but the estimate is the back-transformed ",
                     "model-scale mean, NOT the mean on the original scale; bias is computed on the original scale. ",
                     "Model values outside the transform's range are clamped to its ends."))
}

# response values each family accepts (unknown families: no restriction)
.boot_in_support <- function(fam, v) {
  if (is.na(fam)) return(rep(TRUE, length(v)))
  switch(fam,
         ordbeta = , binomial = , betabinomial = , quasibinomial = v >= 0 & v <= 1,
         beta = , beta_family = v > 0 & v < 1,
         poisson = , quasipoisson = , nbinom1 = , nbinom2 = , nbinom12 = , genpois = , compois = , bell = ,
         tweedie = v >= 0,
         Gamma = , lognormal = , inverse.gaussian = v > 0,
         rep(TRUE, length(v)))
}
.boot_support_label <- function(fam) {
  if (is.na(fam)) return("no restriction")
  switch(fam,
         ordbeta = , binomial = , betabinomial = , quasibinomial = "[0, 1]",
         beta = , beta_family = "(0, 1)",
         poisson = , quasipoisson = , nbinom1 = , nbinom2 = , nbinom12 = , genpois = , compois = , bell = , tweedie = ">= 0",
         Gamma = , lognormal = , inverse.gaussian = "> 0",
         "no restriction")
}

# design matrix of the deployed path: column j is eta with the j-th unit coefficient vector
.boot_design <- function(spec, data, codes) {
  cn <- names(spec$coefficients)
  p <- length(cn)
  cols <- lapply(seq_len(p), function(j) {
    s <- spec
    s$coefficients <- stats::setNames(as.numeric(seq_len(p) == j), cn)
    .eval_eta_codes(s, data, codes)
  })
  matrix(unlist(cols, use.names = FALSE), nrow = nrow(data), ncol = p, dimnames = list(NULL, cn))
}

# dispersion-model design for newdata; factors fixed at contrast codes are averaged
# over levels with the same weights as the linear predictor (so on the log-dispersion scale)
.boot_disp_design <- function(dform, data, stim, spec, codes) {
  dterms <- stats::delete.response(stats::terms(dform))
  dvars <- all.vars(stats::formula(dterms))
  if (length(dvars) == 0L) return(matrix(1, nrow(data), 1, dimnames = list(NULL, "(Intercept)")))
  mini <- list(xlevels = list(), contrasts = list())
  for (v in dvars) {
    col <- stim[[v]]
    if (!(is.factor(col) || is.character(col))) next
    lv <- if (!is.null(spec$xlevels[[v]])) spec$xlevels[[v]] else if (is.factor(col)) levels(col) else sort(unique(as.character(col)))
    mini$xlevels[[v]] <- lv
    if (length(lv) >= 2L) {
      mini$contrasts[[v]] <- if (!is.null(spec$contrasts[[v]])) spec$contrasts[[v]] else
        .contrast_matrix(if (is.factor(col)) attr(col, "contrasts") else NULL, lv, ordered = is.ordered(col))
    }
  }
  mm_at <- function(d) {
    for (v in names(mini$xlevels)) if (v %in% names(d) && is.factor(d[[v]])) attr(d[[v]], "contrasts") <- NULL
    mf <- stats::model.frame(dterms, d, xlev = mini$xlevels, na.action = stats::na.pass)
    stats::model.matrix(dterms, mf, contrasts.arg = if (length(mini$contrasts)) mini$contrasts else NULL)
  }
  rec <- function(d, cds) {
    if (length(cds) == 0L) return(mm_at(d))
    nm <- names(cds)[1]
    w <- .code_weights(mini, nm, cds[[nm]])
    X <- 0
    for (l in names(w)) {
      if (abs(w[[l]]) < 1e-14) next
      d2 <- d; d2[[nm]] <- l
      X <- X + w[[l]] * rec(d2, cds[-1])
    }
    X
  }
  rec(data, codes[intersect(names(codes), dvars)])
}

# E[y] from eta (matrix: rows = cases, cols = replicates); PSI has one row per replicate
.boot_mean_map <- function(ETA, inv_link, PSI = NULL) {
  MU <- array(inv_link(ETA), dim = dim(ETA))
  if (is.null(PSI)) return(MU)
  P1 <- array(stats::plogis(sweep(ETA, 2, PSI[, 2])), dim = dim(ETA))
  PP <- array(stats::plogis(sweep(ETA, 2, PSI[, 1])), dim = dim(ETA))
  P1 + (PP - P1) * MU
}

# k draws per case from the fitted response distribution (single uniform for ordbeta,
# so P(y = 1) = P1 exactly, unlike glmmTMB's simulate() before issue 1359 is fixed)
.boot_draw <- function(fam, eta, phi, psi, inv_link, k) {
  n <- length(eta)
  etaK <- rep(eta, k); muK <- rep(inv_link(eta), k); phiK <- rep(phi, k)
  y <- switch(fam,
              ordbeta = {
                p0 <- 1 - stats::plogis(etaK - psi[1])
                p1 <- stats::plogis(etaK - psi[2])
                u <- stats::runif(n * k)
                yb <- stats::rbeta(n * k, muK * phiK, (1 - muK) * phiK)
                ifelse(u < p0, 0, ifelse(u > 1 - p1, 1, yb))
              },
              beta = stats::rbeta(n * k, muK * phiK, (1 - muK) * phiK),
              gaussian = stats::rnorm(n * k, muK, phiK),
              stop(sprintf("No sampler for family '%s'.", fam)))
  matrix(y, nrow = n, ncol = k)
}

# row-wise percentile limits of a cases x replicates matrix
.boot_q <- function(M, probs) {
  q <- apply(M, 1, stats::quantile, probs = probs, na.rm = TRUE, names = FALSE)
  matrix(q, ncol = 2, byrow = TRUE)
}

# share of replicates below the estimate, ties counted half (rows = cases)
.boot_p_below <- function(M, est) {
  vapply(seq_len(nrow(M)), function(i) {
    m <- M[i, ]; m <- m[is.finite(m)]
    if (!length(m) || !is.finite(est[i])) return(NA_real_)
    (sum(m < est[i]) + 0.5 * sum(m == est[i])) / length(m)
  }, numeric(1))
}

# interval limits per case (rows = cases, cols = replicates).
# percentile: quantiles at probs. bc (bias-corrected percentile, Efron 1981):
# z0 = qnorm(share of replicates below the estimate); limits are the quantiles at
# pnorm(2 * z0 + qnorm(probs)). NA when the estimate lies outside the replicates.
.boot_limits <- function(M, est, probs, method) {
  if (method == "percentile") return(.boot_q(M, probs))
  za <- stats::qnorm(probs)
  p <- .boot_p_below(M, est)
  out <- matrix(NA_real_, nrow(M), 2)
  for (i in seq_len(nrow(M))) {
    if (!is.finite(p[i]) || p[i] <= 0 || p[i] >= 1) next
    m <- M[i, ]; m <- m[is.finite(m)]
    out[i, ] <- stats::quantile(m, stats::pnorm(2 * stats::qnorm(p[i]) + za), names = FALSE)
  }
  out
}

# estimate, bootstrap SE, bias, p_below and interval limits; M = replicates x terms
.boot_coef_table <- function(est, M, probs, method = "percentile") {
  lim <- .boot_limits(t(M), unname(est), probs, method)
  data.frame(term = names(est), estimate = unname(est), boot_se = unname(apply(M, 2, stats::sd, na.rm = TRUE)),
             boot_bias = unname(colMeans(M, na.rm = TRUE)) - unname(est),
             p_below = .boot_p_below(t(M), unname(est)),
             lower = lim[, 1], upper = lim[, 2], row.names = NULL, stringsAsFactors = FALSE)
}

# coefficients of the deployed (printed) equation for one coefficient vector:
# constants absorbed, affine shifts folded, un-scaled to raw units (not simplified)
.boot_equation_coefs <- function(cf, cmap, constants_model, data, scaling) {
  red <- .reduce_columns(cf, cmap, constants_model, data)
  ae <- .expand_affine_pieces(red$beta0, red$terms)
  if (length(scaling) > 0) ae <- .unscale_reduced(ae$beta0, ae$terms, scaling)
  keys <- names(ae$terms)
  lab <- vapply(ae$terms, function(t) paste(vapply(t$pieces, .format_piece_display, character(1)), collapse = ":"), character(1))
  list(coef = c("(Intercept)" = ae$beta0, stats::setNames(vapply(ae$terms, `[[`, numeric(1), "coef"), keys)),
       label = c("(Intercept)" = "(Intercept)", stats::setNames(lab, keys)))
}

# ---- S3 methods ------------------------------------------------------------------
print.boot_mean_model <- function(x, digits = 4, ...) {
  s <- x$settings
  cat("<boot_mean_model>\n")
  cat(sprintf("Cluster bootstrap: %d of %d replicates usable (%.0f s). Clusters: %d '%s'%s.\n",
              s$B_ok, s$B, s$elapsed_sec, s$n_clusters, s$cluster,
              if (is.null(s$strata)) ", unstratified" else sprintf(", resampled within strata of %s", paste(s$strata, collapse = "/"))))
  if (!is.null(s$strata)) {
    cat("  strata (clusters):", paste(sprintf("%s: %d", names(s$strata_sizes), as.integer(s$strata_sizes)), collapse = "; "), "\n")
  }
  if (identical(s$resample, "crossed")) {
    cat(sprintf("Crossed resampling: stimulus rows also resampled%s (%s). Target: new participants AND new comparable stimuli.\n",
                if (is.null(s$stimulus_strata)) ", unstratified" else sprintf(" within strata of %s", paste(s$stimulus_strata, collapse = "/")),
                paste(sprintf("%s: %d", names(s$stimulus_strata_sizes), as.integer(s$stimulus_strata_sizes)), collapse = "; ")))
  } else {
    cat("Stimuli fixed (resample = \"clusters\"): intervals are conditional on this stimulus set.\n")
  }
  cat(sprintf("Stimulus response: %s of trial column '%s' per %s (%s), replacing model column '%s'. Rows per replicate: %s.\n",
              s$stat, s$response, paste(s$by, collapse = " + "), s$by_source, s$target,
              paste(unique(s$rows_per_replicate), collapse = "-")))
  cat(sprintf("Intervals: %s, level %.3g.", if (s$interval == "bc") "bias-corrected percentile (BC)" else "percentile", s$level))
  if (length(s$constants) > 0) cat(sprintf(" Predictions use the deployed constants (%s).",
                                           paste(sprintf("%s = %s", names(s$constants), unlist(s$constants)), collapse = ", ")))
  cat("\n")
  cat(sprintf("Back-transform of predictions: %s\n", s$back_transform$note))
  if (!is.null(s$smoothing)) {
    sm <- s$smoothing
    cat(sprintf(paste0("Smoothed bootstrap: %d trial rows (%d stimuli) shrunk and jittered (Gaussian sigma = %.3g = %.3g x step %.3g)%s%s.\n"),
                sm$n_rows, sm$n_stimuli, sm$sigma, sm$bw, sm$step,
                if (sm$recentre != "none") sprintf("; replicate aggregates recentred (%s) on the originals (|shift| median %.3g, max %.3g)",
                                                   sm$recentre, stats::median(abs(sm$shift[sm$shift != 0])), max(abs(sm$shift))) else ", not recentred",
                if (!is.null(sm$regrid)) sprintf("; then stochastically rounded to a %.3g grid (mean-preserving)", sm$regrid) else ""))
  }
  if (s$interval == "bc" && s$B_ok < 2000) {
    cat(sprintf("  NOTE: only %d usable replicates; bias-corrected limits use more extreme quantiles and need B >= ~2000 (fine for testing).\n", s$B_ok))
  } else if (s$B_ok < 1000) {
    cat(sprintf("  NOTE: only %d usable replicates; percentile limits are imprecise below ~1000 (fine for testing).\n", s$B_ok))
  }
  if (s$interval == "bc") {
    n_na <- sum(is.na(x$predictions$lower)) + sum(is.na(x$coef_model$lower))
    if (n_na > 0) cat(sprintf("  NOTE: %d bias-corrected interval(s) are NA: the estimate lies outside all its replicates.\n", n_na))
  }
  if (s$unmatched_trial_rows > 0) {
    cat(sprintf("  NOTE: %d trial rows match no stimulus row and were ignored; %d clusters with no matching trials at all were dropped from the resampling pool. Check that this is intended (e.g. excluded conditions).\n",
                s$unmatched_trial_rows, s$unmatched_clusters))
  }
  if (x$checks$replicates_design_unchecked > 0) {
    cat(sprintf("  NOTE: the refit design could not be compared with the deployed design for %d usable replicate(s) (rows used or model matrix unavailable); they were kept unchecked.\n",
                x$checks$replicates_design_unchecked))
  }
  if (isTRUE(s$has_weights)) cat("  NOTE: the model has weights; they are kept from the stimulus-level data, NOT recomputed for each replicate.\n")
  nb <- sum(grepl("^off-centre", x$coef_model$note))
  if (nb > 0 || s$pred_off_centre > 0) {
    cat(sprintf(paste0("  NOTE: for %d coefficient(s) and %d prediction row(s) the share of replicates below the estimate (p_below) ",
                       "lies outside 0.4-0.6: the replicate distribution is off-centre, so the percentile limits are shifted ",
                       "relative to the estimate. See the p_below and bias columns.\n"), nb, s$pred_off_centre))
  }
  
  cat("\nModel coefficients (link scale):\n")
  print(format(x$coef_model, digits = digits), row.names = FALSE)
  if (!is.null(x$coef_equation)) {
    cat("\nDeployed-equation coefficients (raw units, constants absorbed, unsimplified; eta scale, i.e. NOT sign-reversed for logit):\n")
    print(format(x$coef_equation, digits = digits), row.names = FALSE)
  }
  cat(sprintf("\nPredictions (%d rows; first 6):\n", nrow(x$predictions)))
  print(format(utils::head(x$predictions), digits = digits), row.names = FALSE)
  if (s$prediction_intervals) {
    cat(sprintf(paste0("Prediction intervals: %d draws per replicate from each refit's fitted distribution (a new stimulus-level ",
                       "value, i.e. a panel aggregate like those fitted), using the %s dispersion. Dispersion evaluated with the same ",
                       "constants as the mean; a factor fixed at a contrast code is averaged on the log-dispersion scale.\n"),
                s$pi_draws, if (s$pi_dispersion == "original") "ORIGINAL fit's" else "replicates' own"))
  }
  f <- x$failures
  if (f$n_failed > 0) {
    cat(sprintf("\nDiscarded replicates: %d\n", f$n_failed))
    print(f$reasons)
  }
  if (length(f$warnings) > 0) {
    cat("\nWarnings during refits (count over replicates, top 5):\n")
    print(utils::head(f$warnings, 5))
  }
  if (f$n_without_pd_hessian > 0) {
    unident <- any(grepl("^not identified", x$coef_model$note))
    share <- f$n_without_pd_hessian / s$B_ok
    cat(sprintf("Usable replicates without a positive-definite Hessian: %d (%.0f%%); they are KEPT (excluding them would select replicates non-randomly). %s\n",
                f$n_without_pd_hessian, 100 * share,
                if (unident) "Likely cause: an unidentified ordbeta cutpoint (see notes); fixing it via start/map removes this." else
                  if (share > 0.1) "WARNING: no unidentified cutpoint explains this; check for weakly identified parameters." else ""))
  }
  invisible(x)
}

# ---- Predictions from an existing run (no refitting) -------------------------
predict.boot_mean_model <- function(object, newdata = NULL, model, target = c("deployed", "original"),
                                    level = NULL, interval = NULL, prediction_intervals = FALSE,
                                    pi_draws = NULL, pi_dispersion = NULL, seed = NULL, ...) {
  target <- match.arg(target)
  if (!inherits(object, "boot_mean_model")) stop("object must be a boot_mean_model() result.", call. = FALSE)
  if (missing(model)) stop("Give model = the deploy_mean_model() object the run was made from.", call. = FALSE)
  if (!inherits(model, "deployed_mean_model")) stop("'model' must be the deploy_mean_model() object the run was made from.", call. = FALSE)
  fit0 <- attr(model, "model")
  if (is.null(fit0)) stop("The deployed object has no attr(, \"model\"): re-run deploy_mean_model() (v17 or later).", call. = FALSE)
  env <- environment(model)
  spec <- env$spec
  scaling <- env$scaling
  numeric_codes <- env$numeric_codes
  constants_levels <- env$constants_levels
  cls <- if (inherits(fit0, "glmmTMB")) "glmmTMB" else "glmgee"
  
  # ---- the run must come from this deployment ----------------------------------
  CF <- object$replicates$coef
  if (is.null(CF)) stop("object has no stored replicate coefficients ($replicates$coef).", call. = FALSE)
  cf0 <- spec$coefficients
  est_run <- object$coef_model$estimate[match(colnames(CF), object$coef_model$term)]
  if (ncol(CF) != length(cf0) || anyNA(est_run) ||
      !isTRUE(all.equal(as.numeric(est_run), as.numeric(cf0), tolerance = 1e-8))) {
    stop("The run's coefficients do not match this deployment's model: it was not made from this deployment.", call. = FALSE)
  }
  if (!identical(object$settings$constants, attr(model, "constants"))) {
    stop("The run's deployment constants differ from this deployment's.", call. = FALSE)
  }
  inv_link <- .get_inverse_link(spec$link)
  fam <- object$settings$family
  ordb <- identical(fam, "ordbeta")
  PSI <- object$replicates$psi
  psi0 <- if (ordb) spec$mean_info$psi else NULL
  if (ordb && is.null(PSI)) stop("object has no stored replicate cutpoints ($replicates$psi).", call. = FALSE)
  
  if (is.null(level)) level <- object$settings$level
  if (!(is.numeric(level) && length(level) == 1L && level > 0 && level < 1)) stop("level must be a single number in (0, 1).", call. = FALSE)
  interval <- if (is.null(interval)) object$settings$interval else match.arg(interval, c("percentile", "bc"))
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  
  # ---- design for newdata -----------------------------------------------------------
  stim <- .ensure_raw_columns(spec$data, scaling)
  idx0 <- .aligned_rows(spec)
  if (is.null(newdata)) newdata <- if (is.null(idx0)) stim else stim[idx0, , drop = FALSE]
  newdata <- as.data.frame(newdata)
  if (nrow(newdata) == 0L) stop("newdata has no rows.", call. = FALSE)
  if (target == "deployed") {
    need <- attr(model, "free_vars")
    miss <- setdiff(need, names(newdata))
    if (length(miss) > 0) stop(sprintf("newdata is missing required column(s): %s", paste(miss, collapse = ", ")), call. = FALSE)
    data_new <- .apply_scaling(newdata, scaling)
    for (v in names(constants_levels)) {
      data_new[[v]] <- if (!is.null(spec$xlevels[[v]])) as.character(constants_levels[[v]]) else constants_levels[[v]]
    }
    codes <- numeric_codes
    X <- .boot_design(spec, data_new, codes)
    eta0 <- as.numeric(X %*% cf0)
    eta_chk <- env$eta_with_constants(data_new)
    dev <- max(abs(eta0 - eta_chk), na.rm = TRUE)
    if (dev > 1e-7 * max(1, max(abs(eta_chk), na.rm = TRUE))) {
      stop(sprintf("The prediction design does not reproduce the deployed linear predictor (max|dev| = %.3g).", dev), call. = FALSE)
    }
    id_cols <- intersect(names(newdata), unique(c(object$settings$by, need)))
  } else {
    # every model variable, in raw units: the deployed inputs plus the variables the
    # deployment fixes (factor levels AND contrast codes; scaled ones by raw name)
    fixed_vars <- unique(c(names(constants_levels), names(numeric_codes)))
    fixed_raw <- vapply(fixed_vars, function(v) if (v %in% names(scaling)) scaling[[v]]$raw_name else v, character(1))
    need <- unique(c(attr(model, "free_vars"), unname(fixed_raw)))
    miss <- setdiff(need, names(newdata))
    if (length(miss) > 0) {
      stop(sprintf("target = \"original\": newdata is missing column(s): %s (the variables the deployment fixes are needed too).",
                   paste(miss, collapse = ", ")), call. = FALSE)
    }
    for (v in intersect(names(spec$xlevels), names(newdata))) {
      bad <- setdiff(unique(as.character(newdata[[v]])), spec$xlevels[[v]])
      if (length(bad) > 0) stop(sprintf("'%s' has value(s) not among the model's levels (%s): %s", v,
                                        paste(spec$xlevels[[v]], collapse = ", "), paste(bad, collapse = ", ")), call. = FALSE)
    }
    data_new <- .apply_scaling(newdata, scaling)
    for (v in intersect(names(spec$xlevels), names(data_new))) data_new[[v]] <- as.character(data_new[[v]])
    codes <- list()
    X <- .boot_design(spec, data_new, codes)
    eta0 <- as.numeric(X %*% cf0)
    dev <- NA_real_
    frozen <- exists("frozen_info", envir = env, inherits = FALSE) && NROW(env$frozen_info) > 0
    if (cls == "glmmTMB" && !frozen) {
      # (skipped with frozen data-dependent terms: glmmTMB would recompute them on newdata)
      # (warnings suppressed: glmmTMB warns of a contrasts mismatch whenever a
      #  factor in newdata carries default contrasts; the check is unaffected)
      eta_chk <- tryCatch(suppressWarnings(as.numeric(stats::predict(fit0, newdata = data_new, type = "link", re.form = NA))),
                          error = function(e) NULL)
      if (!is.null(eta_chk) && length(eta_chk) == length(eta0)) {
        dev <- max(abs(eta0 - eta_chk), na.rm = TRUE)
        if (dev > 1e-6 * max(1, max(abs(eta_chk), na.rm = TRUE))) {
          stop(sprintf("The prediction design does not reproduce glmmTMB's own linear predictor (max|dev| = %.3g).", dev), call. = FALSE)
        }
      }
    }
    id_cols <- intersect(names(newdata), unique(c(object$settings$by, need)))
  }
  
  # ---- E[y] from the stored replicates, back-transformed as in the run ------------
  ETA <- X %*% t(CF)
  E <- .boot_mean_map(ETA, inv_link, PSI)
  E0 <- as.numeric(.boot_mean_map(matrix(eta0, ncol = 1), inv_link, if (ordb) matrix(psi0, nrow = 1) else NULL))
  bt_info <- object$settings$back_transform
  inv <- if (is.function(bt_info$inverse)) bt_info$inverse else identity
  E_out <- matrix(inv(as.vector(E)), nrow = nrow(E))
  E0_out <- inv(E0)
  if (sum(is.finite(E) & !is.finite(E_out)) + sum(is.finite(E0) & !is.finite(E0_out)) > 0) {
    stop("The run's back-transform gave non-finite values for these predictions.", call. = FALSE)
  }
  ci <- .boot_limits(E_out, E0_out, probs, interval)
  pred <- data.frame(estimate = E0_out, bias = rowMeans(E_out, na.rm = TRUE) - E0_out,
                     p_below = .boot_p_below(E_out, E0_out), lower = ci[, 1], upper = ci[, 2])
  
  # ---- prediction intervals (optional) ----------------------------------------------
  if (isTRUE(prediction_intervals)) {
    if (cls != "glmmTMB" || !fam %in% c("ordbeta", "beta", "gaussian")) {
      stop("Prediction intervals are implemented for glmmTMB ordbeta, beta and gaussian models only.", call. = FALSE)
    }
    if (isTRUE(spec$has_re)) stop("Prediction intervals are not implemented for models with random effects.", call. = FALSE)
    if (is.null(pi_dispersion)) pi_dispersion <- object$settings$pi_dispersion
    pi_dispersion <- match.arg(pi_dispersion, c("original", "replicate"))
    if (is.null(pi_draws)) pi_draws <- if (isTRUE(object$settings$prediction_intervals)) object$settings$pi_draws else 10L
    if (!(is.numeric(pi_draws) && length(pi_draws) == 1L && pi_draws >= 1 && pi_draws == round(pi_draws))) {
      stop("pi_draws must be a single positive whole number.", call. = FALSE)
    }
    dform <- tryCatch(fit0[["modelInfo"]][["allForm"]][["dispformula"]], error = function(e) NULL)
    if (is.null(dform)) dform <- ~1
    if (exists(".DATA_SUMMARY_FNS") && any(all.names(dform) %in% .DATA_SUMMARY_FNS)) {
      stop("The dispersion formula contains a data-dependent term (e.g. mean(x)); not supported.", call. = FALSE)
    }
    disp0 <- glmmTMB::fixef(fit0)$disp
    X_disp <- .boot_disp_design(dform, data_new, stim, spec, codes)
    if (!identical(colnames(X_disp), names(disp0))) stop("Dispersion design columns do not match the model's dispersion coefficients.", call. = FALSE)
    if (pi_dispersion == "replicate") {
      if (is.null(object$replicates$disp)) {
        stop("pi_dispersion = \"replicate\" needs the run's replicate dispersion: re-run with prediction_intervals = TRUE, or use \"original\".", call. = FALSE)
      }
      LPHI <- X_disp %*% t(object$replicates$disp)
    } else {
      LPHI <- matrix(as.numeric(X_disp %*% disp0), nrow = nrow(X_disp), ncol = ncol(ETA))
    }
    if (!is.null(seed)) {
      seed_old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
      on.exit({
        if (!is.null(seed_old)) assign(".Random.seed", seed_old, envir = globalenv())
        else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
      }, add = TRUE)
      set.seed(seed)
    }
    pi_draws <- as.integer(pi_draws)
    draws <- matrix(NA_real_, nrow(ETA), ncol(ETA) * pi_draws)
    for (b in seq_len(ncol(ETA))) {
      cols <- (b - 1L) * pi_draws + seq_len(pi_draws)
      draws[, cols] <- .boot_draw(fam, ETA[, b], exp(LPHI[, b]), if (ordb) PSI[b, ] else NULL, inv_link, pi_draws)
    }
    ci_P <- matrix(inv(as.vector(.boot_q(draws, probs))), ncol = 2)
    if (isTRUE(bt_info$decreasing)) ci_P <- ci_P[, 2:1, drop = FALSE]
    pred$pi_lower <- ci_P[, 1]
    pred$pi_upper <- ci_P[, 2]
  }
  
  # ---- extrapolation flag: numeric inputs outside the model data's range ------------
  rows_used <- if (is.null(idx0)) seq_len(nrow(stim)) else idx0
  num_in <- need[vapply(need, function(v) v %in% names(stim) && is.numeric(stim[[v]]) && is.numeric(newdata[[v]]), logical(1))]
  outside <- vapply(seq_len(nrow(newdata)), function(i) {
    paste(num_in[vapply(num_in, function(v) {
      r <- range(stim[[v]][rows_used], na.rm = TRUE)
      x <- newdata[[v]][i]
      is.finite(x) && (x < r[1] - 1e-9 * max(1, abs(r[1])) || x > r[2] + 1e-9 * max(1, abs(r[2])))
    }, logical(1))], collapse = ", ")
  }, character(1))
  pred$outside_range <- outside
  
  out <- list(
    predictions = cbind(newdata[, id_cols, drop = FALSE], pred),
    replicates = list(E = E_out, eta = ETA),
    settings = list(target = target, level = level, interval = interval, B_ok = ncol(ETA),
                    resample = if (is.null(object$settings$resample)) "clusters" else object$settings$resample,
                    constants = if (target == "deployed") attr(model, "constants") else NULL,
                    back_transform = bt_info$type, prediction_intervals = isTRUE(prediction_intervals),
                    n_outside_range = sum(outside != "")),
    checks = list(design_dev = dev)
  )
  class(out) <- "boot_mean_predictions"
  out
}

print.boot_mean_predictions <- function(x, digits = 4, n = 10L, ...) {
  s <- x$settings
  cat("<boot_mean_predictions>\n")
  cat(sprintf("Target: %s. %d rows; %d replicates (%s resampling); %s %.0f%% intervals; back-transform: %s.\n",
              if (s$target == "deployed") sprintf("deployed function (constants: %s)",
                                                  paste(sprintf("%s = %s", names(s$constants), unlist(s$constants)), collapse = ", "))
              else "original model at the given values (no deployment constants)",
              nrow(x$predictions), s$B_ok, s$resample, s$interval, 100 * s$level, s$back_transform))
  if (s$n_outside_range > 0) {
    cat(sprintf("%d row(s) outside the range of the model's data (see outside_range): intervals there assume the model's form holds.\n",
                s$n_outside_range))
  }
  print(utils::head(x$predictions, n), digits = digits, row.names = FALSE)
  if (nrow(x$predictions) > n) cat(sprintf("... %d more rows ($predictions).\n", nrow(x$predictions) - n))
  invisible(x)
}