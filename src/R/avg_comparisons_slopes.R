# =============================================================================
# avg_comparisons_slopes.R
#
# Manual, memory-safe replacements for marginaleffects::avg_predictions(),
# avg_comparisons(), and avg_slopes(), built for:
#   (a) Bayesian brms models requiring random-effects marginalization via
#       sample_new_levels = "gaussian" over large counterfactual grids, where
#       marginaleffects::avg_comparisons()/avg_slopes() can trigger a
#       many-to-many internal merge (see get_comparisons_data.R /
#       comparisons.R -> merge_original_data(), keyed on rowid/rowidcf/by)
#       that blows up on heavily-replicated hand-built newdata.
#   (b) Frequentist marginal (GEE) models, via simulated coefficient draws,
#       so the same downstream summarisation code works for both engines.
#
# DEPENDENCIES: posterior draws are handled as plain matrices/vectors
# throughout (not posterior::rvar) to avoid depending on Math/Ops group
# generic coverage for arbitrary user-supplied `transform` functions.
#   - brms          (Bayesian prediction)
#   - MASS          (frequentist coefficient simulation; King, Tomz &
#                     Wittenberg 2000, AJPS 44(2):341-355)
#   - bayestestR    (p_direction() - do not hand-roll this)
#   - stats         (p.adjust, quantile, model.matrix)
#   - parallel      (optional; only used if a cluster is supplied)
#
# STATISTICAL NOTES (see full discussion in source conversation):
#   - Multiple-comparisons correction: p.adjust() applied to a Bayesian
#     pd-derived p-value analog is a pragmatic hybrid, not a principled
#     Bayesian procedure. The `rope` argument implements the more
#     defensible Bayesian alternative (Kruschke 2018, AMPPS 1(2):270-280):
#     practical-equivalence testing against a domain-justified interval,
#     which does not require a repeated-testing correction at all.
#   - pd -> p-value conversion: p ~= 2(1 - pd), per Makowski, Ben-Shachar,
#     Chen & Ludecke (2019), Frontiers in Psychology, "Indices of Effect
#     Existence and Significance in the Bayesian Framework."
#   - `sample_new_levels = "gaussian"` combined with common random numbers
#     (matched seeds across compared conditions) is what makes the
#     difference/slope computations cheap and low-noise: shared simulated
#     random-effects draws cancel in the contrast rather than adding
#     unrelated variance.
#
# NOT YET IMPLEMENTED / KNOWN GAPS:
#   - lme4 / glmmTMB frequentist mixed models are NOT wired in. Do not
#     bolt on a coefficient-simulation-only path for these (unlike GEE,
#     they have random effects that must also be integrated out). Before
#     extending, check lme4::bootMer() and merTools::predictInterval(),
#     which already implement parametric-bootstrap prediction intervals
#     for lme4 objects - untested against this pipeline, but the right
#     starting point rather than re-deriving a third variant by hand.
#   - Cross-contrasts (simultaneously varying 2+ variables) are not
#     implemented; only one focal `variable` at a time.
#
# MEMORY SIZE CHECK / THRESHOLD CONFIGURATION
#   Every public function below FIRST runs a small real dry run of itself
#   (default 300 rows of newdata, same re_formula/allow_new_levels/
#   sample_new_levels/ndraws/seed as the real call), measures that dry
#   run's actual memory use via gc(), and scales the result linearly to
#   the full grid size (see .estimate_true_size_gb(), v0.6). This requires
#   no hardcoded assumptions about how a given model/family allocates
#   memory internally. If the resulting estimate exceeds the configured
#   threshold, you are shown the estimate and prompted with a numbered
#   Yes/No confirmation (utils::menu()) before the real, full-size call
#   runs. In non-interactive sessions (Rscript, knitr, etc.) this raises
#   an error instead of prompting, since there is no user available to
#   confirm.
#
#   The threshold has a script-level default but is fully configurable via
#   a single global option, WITHOUT needing to touch every function call:
#
#     options(mfx_manual.size_warn_gb = 6)   # set once, e.g. at top of script
#
#   or per-call, which overrides the global option for that call only:
#
#     avg_predictions_manual(..., size_warn_gb = 6)
#
#   Choosing a value: treat it as a fraction (roughly a third to a half) of
#   your machine's TOTAL installed RAM, not free RAM - R, the OS, and
#   everything else running need headroom too. There is no universally
#   "correct" number; this was the practical driver for making it
#   configurable rather than hardcoded. Check total RAM via Windows Task
#   Manager (Performance tab) or `memory.limit()` / `parallel::detectCores()`
#   as a rough proxy if unsure.
#
# VERSION HISTORY
#   v0.1  Initial release. Bayesian (brms) + frequentist (glmgee) engines.
#         avg_predictions_manual, avg_comparisons_manual (pairwise/
#         reference/sequential + ROPE + p.adjust), avg_slopes_manual
#         (finite-difference). Optional parallel cluster support via `cl`.
#   v0.2  Added mandatory pre-flight memory size check with interactive
#         Yes/No confirmation (utils::menu()) before any large draws
#         matrix is allocated, configurable via options(mfx_manual.size_warn_gb)
#         or the size_warn_gb argument. Non-interactive sessions now abort
#         with an informative error above threshold rather than prompting.
#   v0.3  Size estimate now corrected by an engine/family-specific overhead
#         multiplier (.memory_overhead_multiplier()) rather than the naive
#         nrow x ndraws x 8-bytes formula alone. Default ordbeta multiplier
#         (6x) is a PROVISIONAL placeholder only - see v0.4 correction below
#         before trusting it. GEE/quasibinomial defaults to 1x (single
#         linear predictor, no mixture structure), giving the same nominal
#         grid size a much looser (more honest) threshold for that engine
#         without needing a separate threshold number. Multipliers are
#         configurable via options(mfx_manual.overhead_multiplier_ordbeta =
#         <value>), _brms, and _gee.
#   v0.4  CORRECTION: the 6x ordbeta multiplier in v0.3 was calibrated from
#         a single observed case (6.0 GB estimated vs. ~35 GB process RSS)
#         that was confounded by other objects already resident in the R
#         session (multiple other brms model objects, unrelated data).
#         gc() after that call dropped usage to <6 GB, confirming most of
#         the 35 GB was reclaimable garbage rather than persistent unrelated
#         objects - but this does NOT cleanly isolate how much of that
#         garbage was actually generated by the flagged call versus
#         accumulated earlier in the session. The 6x default is retained as
#         a conservative placeholder (better to over- than under-warn) but
#         should be re-derived using profile_call_memory() below - via
#         gc(reset = TRUE) immediately before the call and gc() immediately
#         after, isolating "max used" from the pre-call baseline - in a
#         session with nothing else large resident, before being trusted
#         for real threshold decisions.
#   v0.5  A properly isolated gc(reset=TRUE)-based measurement (avg_comparisons_manual
#         on a custom/forked ordbeta family) showed a true incremental peak
#         of ~35.5 GB against a DISPLAYED estimate of 6.0 GB - i.e. still a
#         ~6x undershoot even with the v0.3/v0.4 multiplier logic active.
#         Added an explicit overhead_multiplier override argument on the
#         (WRONG - see v0.6) theory that model$family$family on a
#         custom/forked family object might not contain a plain "beta"
#         substring, silently falling through to the generic 2x branch.
#   v0.6  CORRECTION: model$family$family on the actual custom/forked
#         ordbeta family used in this project returned exactly "ordbeta"
#         (confirmed directly), which DOES match the "beta" substring check
#         - so the v0.5 diagnosis was wrong; the 6x branch almost certainly
#         fired correctly and simply wasn't a big enough multiplier. Rather
#         than continue guessing at (or asking the user to measure) a
#         constant multiplier - which has now caused two wrong diagnoses in
#         a row - the size check is REDESIGNED: .memory_overhead_multiplier(),
#         estimate_draws_size_gb(), and the overhead_multiplier argument are
#         REMOVED. avg_predictions_manual(), avg_comparisons_manual(), and
#         avg_slopes_manual() now call .estimate_true_size_gb(), which runs
#         a small REAL version of the actual call (default 300 rows, same
#         re_formula/allow_new_levels/sample_new_levels/ndraws/seed as the
#         real call) and measures its ACTUAL memory use via gc(reset=TRUE),
#         then scales that measurement linearly to the full grid size. This
#         needs no knowledge of family internals, requires no manual
#         profiling step from the user, and self-corrects if internal
#         memory behaviour changes (brms updates, different custom
#         families) without anyone needing to notice and update a constant.
#         Adds a small time cost (one small extra prediction call, typically
#         a few seconds) in exchange for not needing a trusted multiplier at
#         all. profile_call_memory() is retained as a general-purpose
#         utility for ad hoc profiling but is no longer a required step in
#         the normal workflow.
#   v0.7  BUG: the v0.6 single-sample linear extrapolation (measure at 300
#         rows, multiply by n_full/300) produced a wildly wrong estimate in
#         practice (>14,000 GB predicted for a call whose true peak was
#         ~35 GB) because it assumed peak memory is purely proportional to
#         nrow(newdata). If a substantial share of memory use is a ONE-TIME
#         fixed cost independent of grid size, scaling the whole
#         measurement by a large factor inflates that fixed cost by the
#         same factor - which is what happened. Fixed by measuring at TWO
#         sample sizes (default 100 and 2000) and fitting a line
#         (fixed_overhead + per_row_cost * n) via lm(), correctly
#         separating the two components instead of assuming one is zero.
#         dry_run_sample_size (singular) is now dry_run_sample_sizes
#         (plural, length >= 2) on all three public functions.
# =============================================================================


# -----------------------------------------------------------------------------
# Memory profiling helper (general-purpose utility; the size-check
# machinery below no longer requires calling this manually - see v0.6)
# -----------------------------------------------------------------------------

#' Measure the true incremental peak memory (MB) of one call
#'
#' Uses gc(reset = TRUE) to zero out R's peak-usage tracking immediately
#' before evaluating `expr`, then reads back the peak reached during
#' evaluation - isolating this call's footprint from whatever else is
#' already resident in the session (unlike reading Task Manager / total
#' process RSS, which bundles in everything else loaded).
#' Measure the true incremental peak memory (MB) of one call
#'
#' Uses gc(reset = TRUE) to zero out R's peak-usage tracking immediately
#' before evaluating `expr`, then reads back the peak reached during
#' evaluation - isolating this call's footprint from whatever else is
#' already resident in the session. Exposed as a standalone utility for ad
#' hoc profiling; the size-check machinery below now uses the same
#' technique internally and automatically, so calling this yourself is
#' optional, not required.
#'
#' @param expr An unevaluated expression, e.g. via quote().
#' @return A list with $result (the expression's value) and
#'   $peak_incremental_mb (the isolated peak memory increase, in MB).
profile_call_memory <- function(expr) {
  gc(reset = TRUE, full = TRUE)
  g0 <- gc()
  mb_cols0 <- which(colnames(g0) == "(Mb)")
  baseline_mb <- sum(g0[, mb_cols0[1]])
  
  result <- eval(expr, envir = parent.frame())
  
  g1 <- gc()
  mb_cols1 <- which(colnames(g1) == "(Mb)")
  peak_mb <- sum(g1[, mb_cols1[length(mb_cols1)]])
  
  list(result = result, peak_incremental_mb = peak_mb - baseline_mb)
}


# -----------------------------------------------------------------------------
# Small utility (defined locally for compatibility with R < 4.4, where
# %||% is not yet a base operator)
# -----------------------------------------------------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x


# -----------------------------------------------------------------------------
# Cluster helper
# -----------------------------------------------------------------------------

#' Build a persistent parallel cluster with the model exported once
#'
#' @param model A brmsfit or glmgee object to export to all workers.
#' @param n_workers Number of cluster workers.
#' @param extra_exports Character vector of additional object names (from the
#'   calling environment) to export to workers, e.g. custom transform helpers.
#' @return A cluster object for use as the `cl` argument below. Remember to
#'   `parallel::stopCluster(cl)` when done.
make_prediction_cluster <- function(model, n_workers = parallel::detectCores() - 1,
                                    extra_exports = character(0)) {
  cl <- parallel::makeCluster(max(1, n_workers))
  parallel::clusterEvalQ(cl, { library(brms); library(MASS) })
  parallel::clusterExport(cl, c("model", extra_exports), envir = environment())
  cl
}


# -----------------------------------------------------------------------------
# Memory size check - AUTOMATIC, empirically measured (no user-supplied
# multiplier, no family detection/guessing required)
# -----------------------------------------------------------------------------
#
# Rationale for this design (see v0.6 note in the file header): earlier
# versions tried to predict memory use analytically via a
# nrow x ndraws x 8-bytes formula corrected by a hand-set "overhead
# multiplier" per model family. This required the user (or the file's
# author) to guess or measure that multiplier separately, and a wrong
# guess (or a wrongly-diagnosed reason for a wrong guess) directly caused
# a system freeze earlier in this project's development. Rather than
# refine the guess further, the size check now RUNS A SMALL REAL VERSION
# of the actual call, on a subsample of newdata, and measures its actual
# memory use directly - then scales that measurement up linearly to the
# full grid size. This is slower by the cost of one small extra
# prediction call (typically a few seconds), but it needs no knowledge of
# how any particular family/engine allocates memory internally, and it is
# self-correcting if that internal behaviour ever changes (e.g. a brms
# update, a different custom family) without anyone needing to notice and
# update a hardcoded constant.

.gc_mb <- function(stat = c("used", "max_used")) {
  stat <- match.arg(stat)
  g <- gc()
  mb_cols <- which(colnames(g) == "(Mb)")
  col <- if (stat == "used") mb_cols[1] else mb_cols[length(mb_cols)]
  sum(g[, col])
}

#' Effective number of draws a model call will produce (used for the dry
#' run itself, not for any size formula)
.effective_draws <- function(model, ndraws = NULL, nsim = 1000) {
  if (inherits(model, "brmsfit")) {
    if (is.null(ndraws)) return(brms::ndraws(model))
    return(ndraws)
  }
  if (inherits(model, "glmgee")) return(nsim)
  stop("Cannot resolve an effective draw count for this model class.")
}

#' Empirically estimate the true peak memory (GB) of a full-size call by
#'
#' Measures at TWO proportionally-sized samples (see .choose_dry_run_sizes()),
#' not one, and fits a line through them (fixed_overhead + per_row_cost * n)
#' rather than assuming peak memory is purely proportional to nrow(newdata).
#' Pure proportionality from a single small measurement is a real bug, not
#' a hedge: if a large share of memory use is a ONE-TIME fixed cost (e.g.
#' whatever brms sets up once per posterior_epred() call, independent of
#' how many rows are being predicted), scaling that fixed cost by
#' (n_full / small_sample_size) inflates it by that same large factor,
#' producing wild overestimates for small dry-run samples relative to
#' large real grids - which is exactly what happened in practice (a
#' fixed-300-row dry run produced a >14,000 GB estimate for a call whose
#' true peak was ~35 GB). Fitting a line through two points separates the
#' fixed and variable components properly.
#'
#' @param model brmsfit or glmgee object.
#' @param newdata The FULL newdata grid the real call will use.
#' @param n_conditions Number of separate draws matrices the real call will
#'   produce in sequence (1 for a prediction, length(levels) for a
#'   comparison, 2 for a slope's hi/lo pair).
#' @param force_recalibrate If TRUE, ignores any cached calibration for
#'   this model/settings combination and re-measures from scratch.
#' @param ... Passed to .dispatch_draws() (re_formula, allow_new_levels,
#'   sample_new_levels, ndraws, nsim, seed) - MUST match what the real call
#'   will use, since these settings affect memory behaviour.
#' @return Estimated size in GB for the FULL call.
#'
#' CACHING: fitted (fixed_overhead_mb, per_row_mb) coefficients are cached
#' in .mfx_size_cache, keyed by model class/family + the settings that
#' actually affect memory behaviour (re_formula, allow_new_levels,
#' sample_new_levels, ndraws/nsim). In typical usage (many calls against
#' the same model with the same settings, only newdata/variable changing -
#' e.g. testing several different factors one after another), this means
#' the dry run only runs ONCE per session for that model/settings
#' combination; every subsequent call reuses the cached line and costs
#' essentially nothing. Call reset_size_cache() to force recalibration
#' (e.g. after updating the model or brms itself).
#'
#' SAMPLE SIZE SCALING: dry-run sample sizes are chosen as a small
#' percentage of nrow(newdata), capped, rather than fixed absolute values -
#' so the dry run never becomes a meaningful fraction of the real
#' computation, and is skipped entirely (falling back to directly
#' measuring the real call, uncached-but-cheap since the grid is already
#' small) when nrow(newdata) is too small for a meaningful two-point fit
#' to make sense in the first place.
#'
#' EXTRAPOLATION CAVEAT: the cached line is fit from small samples and
#' applied to your actual grid size. If a later call uses a MUCH larger
#' newdata than whatever calibrated the cache (e.g. calibrated on a few
#' thousand rows, later applied to hundreds of thousands), that is
#' extrapolation beyond the calibration range and carries the usual risks
#' of extrapolation - not a guarantee. Call reset_size_cache() if you
#' suspect the cached line no longer fits your current scale of grid.
.mfx_size_cache <- new.env(parent = emptyenv())

#' Clear the memory-calibration cache, forcing recalibration on next call
reset_size_cache <- function() {
  rm(list = ls(.mfx_size_cache), envir = .mfx_size_cache)
  invisible(NULL)
}

.calibration_key <- function(model, ...) {
  dots <- list(...)
  fam <- if (inherits(model, "brmsfit")) {
    tryCatch(as.character(model$family$family), error = function(e) "NA")
  } else {
    "NA"
  }
  paste(
    paste(class(model), collapse = "/"), fam,
    paste(deparse(dots$re_formula), collapse = ""),
    isTRUE(dots$allow_new_levels),
    dots$sample_new_levels %||% "uncertainty",
    dots$ndraws %||% NA, dots$nsim %||% NA,
    sep = "||"
  )
}

.choose_dry_run_sizes <- function(n_full) {
  s_lo <- max(20, min(300, ceiling(0.01 * n_full)))
  s_hi <- max(s_lo + 20, min(2000, ceiling(0.05 * n_full)))
  sort(unique(c(s_lo, s_hi)[c(s_lo, s_hi) < n_full]))
}

.estimate_true_size_gb <- function(model, newdata, n_conditions = 1,
                                   force_recalibrate = FALSE, ...) {
  n_full <- nrow(newdata)
  key <- .calibration_key(model, ...)
  
  if (!force_recalibrate && exists(key, envir = .mfx_size_cache, inherits = FALSE)) {
    cached <- get(key, envir = .mfx_size_cache, inherits = FALSE)
    predicted_mb <- cached$fixed_overhead_mb + cached$per_row_mb * n_full
    return((predicted_mb / 1024) * n_conditions)
  }
  
  .measure_at <- function(n) {
    idx <- sample(seq_len(n_full), n)
    small_data <- newdata[idx, , drop = FALSE]
    gc(reset = TRUE, full = TRUE)
    baseline_mb <- .gc_mb("used")
    invisible(.dispatch_draws(model, small_data, ...))
    peak_mb <- .gc_mb("max_used")
    gc(full = TRUE)
    peak_mb - baseline_mb
  }
  
  sizes <- .choose_dry_run_sizes(n_full)
  
  if (length(sizes) < 2) {
    # Grid too small for a meaningful two-point fit - measure the real
    # call directly (this only runs once per model/settings combination
    # per session; the result is still cached below for next time).
    incremental_mb <- .measure_at(n_full)
    assign(key, list(fixed_overhead_mb = 0, per_row_mb = incremental_mb / n_full),
           envir = .mfx_size_cache)
    return((incremental_mb / 1024) * n_conditions)
  }
  
  measured_mb <- vapply(sizes, .measure_at, numeric(1))
  fit <- stats::lm(measured_mb ~ sizes)
  fixed_overhead_mb <- max(0, unname(coef(fit)[1]))
  per_row_mb <- max(0, unname(coef(fit)[2]))
  
  assign(key, list(fixed_overhead_mb = fixed_overhead_mb, per_row_mb = per_row_mb),
         envir = .mfx_size_cache)
  
  predicted_mb <- fixed_overhead_mb + per_row_mb * n_full
  (predicted_mb / 1024) * n_conditions
}

#' Pre-flight check: warn and require confirmation above a size threshold
#'
#' @param estimated_gb Output of .estimate_true_size_gb().
#' @param size_warn_gb Threshold in GB. If NULL, falls back to
#'   getOption("mfx_manual.size_warn_gb", default = 4).
#' @param context Short string describing which call this is, shown in the
#'   prompt/error (e.g. "avg_comparisons_manual(variable = 'Sex')").
.confirm_large_computation <- function(estimated_gb, size_warn_gb = NULL, context = "") {
  threshold <- if (is.null(size_warn_gb)) {
    getOption("mfx_manual.size_warn_gb", default = 4)
  } else {
    size_warn_gb
  }
  
  if (estimated_gb <= threshold) return(invisible(TRUE))
  
  msg <- sprintf(
    paste0(
      "%s\n",
      "Estimated draws matrix size (measured via a small real dry run, ",
      "scaled to your full grid): %.1f GB (threshold: %.1f GB).\n",
      "This is the scale that has previously caused system instability on large grids.\n",
      "To change this threshold: options(mfx_manual.size_warn_gb = <value>), ",
      "or pass size_warn_gb = <value> directly to this call."
    ),
    context, estimated_gb, threshold
  )
  
  if (!interactive()) {
    stop(
      msg,
      "\nRunning in a non-interactive session - refusing to proceed automatically.",
      "\nRaise the threshold first if you are sure this is safe:",
      "\n  options(mfx_manual.size_warn_gb = <value>)  # or pass size_warn_gb= to this call",
      call. = FALSE
    )
  }
  
  choice <- utils::menu(c("Yes - proceed anyway", "No - abort"), title = msg)
  if (choice != 1) {
    stop("Aborted by user after memory size check.", call. = FALSE)
  }
  invisible(TRUE)
}


# -----------------------------------------------------------------------------
# Core prediction dispatch - one function, two engines
# -----------------------------------------------------------------------------

#' Generate a draws matrix (ndraws x nrow(newdata)) for a model + newdata
#'
#' Dispatches on model class so downstream functions are engine-agnostic.
#' Extend the if/else chain here to add new model classes; do not duplicate
#' the comparison/slope logic below per engine.
.dispatch_draws <- function(model, newdata,
                            re_formula = NULL, allow_new_levels = FALSE,
                            sample_new_levels = "uncertainty",
                            ndraws = NULL, nsim = 1000, seed = 1) {
  
  if (inherits(model, "brmsfit")) {
    set.seed(seed)
    return(brms::posterior_epred(
      model, newdata = newdata, re_formula = re_formula,
      allow_new_levels = allow_new_levels,
      sample_new_levels = sample_new_levels, ndraws = ndraws
    ))
  }
  
  if (inherits(model, "glmgee")) {
    # Marginal/GEE model: no random effects to integrate out, so simulated
    # coefficient draws alone correctly propagate parameter uncertainty
    # through the (possibly nonlinear) link. King, Tomz & Wittenberg (2000).
    set.seed(seed)
    b_sim <- MASS::mvrnorm(nsim, stats::coef(model), stats::vcov(model))
    X <- stats::model.matrix(stats::delete.response(stats::terms(model)), data = newdata)
    eta <- X %*% t(b_sim)
    pred <- model$family$linkinv(eta)
    return(t(pred))  # ndraws x nrow(newdata), matching the brms return shape
  }
  
  stop(
    "Model class '", paste(class(model), collapse = "/"),
    "' is not supported by .dispatch_draws(). ",
    "See 'NOT YET IMPLEMENTED' notes in the file header before extending."
  )
}

#' Same as .dispatch_draws() but runs multiple newdata variants in parallel
#'
#' @param newdata_list A named list of newdata data.frames, one per condition
#'   to predict (e.g. one per factor level being compared).
#' @param cl Optional cluster from make_prediction_cluster(). If NULL, runs
#'   sequentially (usually fine unless you have many levels and a slow model).
.dispatch_draws_multi <- function(model, newdata_list, cl = NULL, ...) {
  fun <- function(nd, model, ...) .dispatch_draws(model, nd, ...)
  if (is.null(cl)) {
    lapply(newdata_list, fun, model = model, ...)
  } else {
    parallel::parLapply(cl, newdata_list, fun, model = model, ...)
  }
}


# -----------------------------------------------------------------------------
# avg_predictions_manual()
# -----------------------------------------------------------------------------

#' Manual replacement for marginaleffects::avg_predictions()
#'
#' @param model brmsfit or glmgee object.
#' @param newdata Grid of predictor values (already expanded/replicated as
#'   needed - see grid-construction guidance in the source conversation).
#' @param by Character vector of column names in `newdata` to group by
#'   before averaging. NULL averages over the whole grid into one estimate.
#' @param transform Function applied to the raw draws matrix before
#'   averaging, e.g. a back-transform from a scaled/logit space.
#' @param conf_level Credible/confidence interval coverage (default 0.95).
#' @param size_warn_gb Memory size threshold (GB) for the pre-flight check.
#'   NULL (default) uses getOption("mfx_manual.size_warn_gb", 4). See file
#'   header for guidance on choosing this value.
#' @param force_recalibrate If TRUE, ignores any cached memory calibration
#'   for this model/settings combination and re-measures from scratch.
#'   Usually not needed - see reset_size_cache() to clear the whole cache.
#' @param ... Passed to .dispatch_draws() (re_formula, allow_new_levels,
#'   sample_new_levels, ndraws, nsim, seed).
avg_predictions_manual <- function(model, newdata, by = NULL,
                                   transform = identity, conf_level = 0.95,
                                   size_warn_gb = NULL,
                                   force_recalibrate = FALSE, ...) {
  est_gb <- .estimate_true_size_gb(model, newdata, n_conditions = 1,
                                   force_recalibrate = force_recalibrate, ...)
  .confirm_large_computation(est_gb, size_warn_gb, context = "avg_predictions_manual()")
  
  draws <- .dispatch_draws(model, newdata, ...)
  resp <- transform(draws)
  alpha <- 1 - conf_level
  
  summarise_cols <- function(cols) {
    draw_means <- rowMeans(resp[, cols, drop = FALSE])
    data.frame(
      estimate = mean(draw_means),
      conf.low = unname(stats::quantile(draw_means, alpha / 2)),
      conf.high = unname(stats::quantile(draw_means, 1 - alpha / 2))
    )
  }
  
  if (is.null(by)) return(summarise_cols(seq_len(ncol(resp))))
  
  grp <- interaction(newdata[by], drop = TRUE)
  col_groups <- split(seq_len(ncol(resp)), grp)
  lookup <- as.data.frame(newdata[!duplicated(grp), by, drop = FALSE])
  rownames(lookup) <- as.character(grp[!duplicated(grp)])
  
  out <- do.call(rbind, lapply(names(col_groups), function(g) {
    cbind(lookup[g, , drop = FALSE], summarise_cols(col_groups[[g]]))
  }))
  rownames(out) <- NULL
  out
}


# -----------------------------------------------------------------------------
# avg_comparisons_manual()
# -----------------------------------------------------------------------------

#' Manual replacement for marginaleffects::avg_comparisons()
#'
#' Supports 2+ factor levels via pairwise / reference / sequential contrasts,
#' optional ROPE-based practical-equivalence judgement, and optional
#' multiple-comparisons p-value adjustment (see file header caveats on the
#' latter).
#'
#' @param model brmsfit or glmgee object.
#' @param newdata Grid of predictor values, at the `variable`'s CURRENT
#'   (baseline) values - these get overwritten per level internally.
#' @param variable Name of the focal categorical column in `newdata`.
#' @param comparison_type One of "pairwise", "reference", "sequential".
#' @param levels Optional character vector restricting/ordering which
#'   levels of `variable` to compare. Defaults to all unique observed levels.
#' @param transform Function applied to each level's raw draws matrix.
#' @param conf_level Credible/confidence interval coverage (default 0.95).
#' @param p_adjust Method passed to stats::p.adjust() ("none", "holm",
#'   "BH", "bonferroni", ...). See file header caveat before using this.
#' @param rope Optional numeric length-2 vector, e.g. c(-0.2, 0.2), defining
#'   a region of practical equivalence to zero for a Kruschke (2018)-style
#'   decision column, independent of any p-value/adjustment.
#' @param cl Optional cluster from make_prediction_cluster() to parallelize
#'   prediction across levels of `variable`.
#' @param size_warn_gb Memory size threshold (GB) for the pre-flight check,
#'   applied to the TOTAL across all levels being predicted (they may be
#'   computed sequentially, but the check is conservative and sizes for the
#'   full set). NULL (default) uses getOption("mfx_manual.size_warn_gb", 4).
#' @param ... Passed to .dispatch_draws() (re_formula, allow_new_levels,
#'   sample_new_levels, ndraws, nsim, seed).
avg_comparisons_manual <- function(model, newdata, variable,
                                   comparison_type = c("pairwise", "reference", "sequential"),
                                   levels = NULL, transform = identity,
                                   conf_level = 0.95, p_adjust = "none",
                                   rope = NULL, cl = NULL,
                                   size_warn_gb = NULL,
                                   dry_run_sample_sizes = c(100, 2000), ...) {
  comparison_type <- match.arg(comparison_type)
  if (is.null(levels)) levels <- sort(unique(as.character(newdata[[variable]])))
  
  est_gb <- .estimate_true_size_gb(model, newdata, n_conditions = length(levels),
                                   sample_sizes = dry_run_sample_sizes, ...)
  .confirm_large_computation(
    est_gb, size_warn_gb,
    context = sprintf("avg_comparisons_manual(variable = '%s', %d levels)", variable, length(levels))
  )
  
  pairs <- switch(comparison_type,
                  pairwise   = t(utils::combn(levels, 2)),
                  reference  = cbind(levels[-1], levels[1]),
                  sequential = cbind(levels[-1], levels[-length(levels)])
  )
  
  # Predict once per level (not once per pair) and reuse - avoids redundant
  # prediction of shared levels across multiple pairwise comparisons.
  newdata_list <- lapply(levels, function(lv) { nd <- newdata; nd[[variable]] <- lv; nd })
  names(newdata_list) <- levels
  draws_list <- .dispatch_draws_multi(model, newdata_list, cl = cl, ...)
  mean_draws_by_level <- lapply(draws_list, function(d) rowMeans(transform(d)))
  names(mean_draws_by_level) <- levels
  
  alpha <- 1 - conf_level
  results <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
    hi <- pairs[i, 1]; lo <- pairs[i, 2]
    d <- mean_draws_by_level[[hi]] - mean_draws_by_level[[lo]]
    pd <- as.numeric(bayestestR::p_direction(d))
    row <- data.frame(
      contrast = paste(hi, "-", lo),
      estimate = mean(d),
      conf.low = unname(stats::quantile(d, alpha / 2)),
      conf.high = unname(stats::quantile(d, 1 - alpha / 2)),
      pd = pd,
      p.value = 2 * (1 - pd)
    )
    if (!is.null(rope)) {
      row$rope_decision <- if (row$conf.low > rope[2] || row$conf.high < rope[1]) {
        "outside ROPE (credible non-negligible effect)"
      } else if (row$conf.low > rope[1] && row$conf.high < rope[2]) {
        "inside ROPE (practically equivalent to zero)"
      } else {
        "overlaps ROPE (undecided)"
      }
    }
    row
  }))
  results$p.value.adj <- stats::p.adjust(results$p.value, method = p_adjust)
  results
}


# -----------------------------------------------------------------------------
# avg_slopes_manual()
# -----------------------------------------------------------------------------

#' Manual replacement for marginaleffects::avg_slopes()
#'
#' Central finite-difference slope, computed by differencing the
#' BACK-TRANSFORMED predictions (not transforming an already-computed
#' link-scale slope, which is wrong for a nonlinear transform).
#'
#' @param model brmsfit or glmgee object.
#' @param newdata Grid of predictor values.
#' @param variable Name of the focal continuous column in `newdata`.
#' @param eps Step size for the finite difference. Defaults to
#'   0.0001 * range(variable), matching marginaleffects' own default.
#' @param transform Function applied to each side's raw draws matrix.
#' @param conf_level Credible/confidence interval coverage (default 0.95).
#' @param cl Optional cluster (see avg_comparisons_manual()).
#' @param size_warn_gb Memory size threshold (GB) for the pre-flight check
#'   (sized for the hi/lo pair, i.e. n_conditions = 2). NULL (default) uses
#'   getOption("mfx_manual.size_warn_gb", 4).
#' @param ... Passed to .dispatch_draws().
avg_slopes_manual <- function(model, newdata, variable, eps = NULL,
                              transform = identity, conf_level = 0.95,
                              cl = NULL, size_warn_gb = NULL,
                              dry_run_sample_sizes = c(100, 2000), ...) {
  x <- newdata[[variable]]
  if (is.null(eps)) eps <- 0.0001 * diff(range(x))
  
  est_gb <- .estimate_true_size_gb(model, newdata, n_conditions = 2,
                                   sample_sizes = dry_run_sample_sizes, ...)
  .confirm_large_computation(
    est_gb, size_warn_gb,
    context = sprintf("avg_slopes_manual(variable = '%s')", variable)
  )
  
  nd_hi <- newdata; nd_hi[[variable]] <- x + eps / 2
  nd_lo <- newdata; nd_lo[[variable]] <- x - eps / 2
  
  draws_list <- .dispatch_draws_multi(model, list(hi = nd_hi, lo = nd_lo), cl = cl, ...)
  slope_draws <- rowMeans(transform(draws_list$hi) - transform(draws_list$lo)) / eps
  
  pd <- as.numeric(bayestestR::p_direction(slope_draws))
  alpha <- 1 - conf_level
  data.frame(
    variable = variable,
    estimate = mean(slope_draws),
    conf.low = unname(stats::quantile(slope_draws, alpha / 2)),
    conf.high = unname(stats::quantile(slope_draws, 1 - alpha / 2)),
    pd = pd,
    p.value = 2 * (1 - pd)
  )
}


# =============================================================================
# USAGE EXAMPLES
# (not run automatically - copy into your script and adapt object names)
# =============================================================================
if (FALSE) {
  
  source("avg_comparisons_slopes.R")
  
  # --- Example 0: the size check runs automatically inside every call below
  #     - there's nothing extra to run first. It performs a small real dry
  #     run of the exact call (default 300 rows), measures actual memory
  #     use, and only prompts if the scaled-up estimate exceeds threshold.
  #     If your grid is very large and even the dry run itself feels slow,
  #     reduce it: dry_run_sample_sizes = c(50, 500)
  
  # Set the threshold once for a whole session (roughly a third to a half
  # of your machine's TOTAL RAM, not free RAM). Skip this and every call
  # below uses the built-in default (4 GB).
  options(mfx_manual.size_warn_gb = 6)
  
  # A call that exceeds the threshold will show something like:
  #   avg_comparisons_manual(variable = 'NationGeo2', 2 levels)
  #   Estimated draws matrix size (measured via a small real dry run,
  #   scaled to your full grid): 8.2 GB (threshold: 6.0 GB).
  #   This is the scale that has previously caused system instability on
  #   large grids.
  #   To change this threshold: options(mfx_manual.size_warn_gb = <value>),
  #   or pass size_warn_gb = <value> directly to this call.
  #
  #   1: Yes - proceed anyway
  #   2: No - abort
  #   Selection:
  #
  # In a non-interactive session (Rscript / knitr), the same situation
  # raises an error instead of prompting, so batch jobs never hang waiting
  # for input they can't receive.
  
  # --- Example 1: average predictions (Bayesian, brms) -----------------------
  m3bpredType <- avg_predictions_manual(
    m3b, newdata = subset(grid_marg, UASEvents == 1),
    by = c("UASLAEMaxLRScl", "UASType", "AmbientEnv", "UASEvents"),
    re_formula = NULL, allow_new_levels = TRUE, sample_new_levels = "gaussian",
    ndraws = ndraws_pred, seed = 999,
    transform = \(x) betaTransform(x, -10, 10, direction = "reverse", squeeze = "none")
  )
  
  # --- Example 2: pairwise contrast with ROPE, no p.adjust (recommended
  #     as the primary/headline result - see file header notes) -------------
  sex_contrast <- avg_comparisons_manual(
    m3b, newdata = grid_test, variable = "Sex", comparison_type = "pairwise",
    re_formula = NULL, allow_new_levels = TRUE, sample_new_levels = "gaussian",
    ndraws = ndraws_pred, seed = 999,
    rope = c(-0.2, 0.2),  # <-- set this from domain knowledge, not by default
    transform = \(x) betaTransform(x, -10, 10, direction = "reverse", squeeze = "none")
  )
  
  # --- Example 3: same contrast, WITH a frequentist-style multiple-
  #     comparisons correction as a labelled supplementary sensitivity
  #     check (see file header caveat before treating this as primary) ------
  all_demo_contrasts <- do.call(rbind, lapply(
    c("Sex", "NationGeo2", "AAM_attitude2", "Home_Area2", "Area_soundscape2"),
    function(v) {
      res <- avg_comparisons_manual(
        m3b, newdata = grid_test, variable = v,
        re_formula = NULL, allow_new_levels = TRUE, sample_new_levels = "gaussian",
        ndraws = ndraws_pred, seed = 999,
        transform = \(x) betaTransform(x, -10, 10, direction = "reverse", squeeze = "none")
      )
      cbind(variable = v, res)
    }
  ))
  all_demo_contrasts$p.value.adj <- stats::p.adjust(all_demo_contrasts$p.value, method = "holm")
  
  # --- Example 4: slope of the continuous focal predictor -------------------
  lae_slope <- avg_slopes_manual(
    m3b, newdata = subset(grid_marg, UASEvents == 1), variable = "UASLAEMaxLRScl",
    re_formula = NULL, allow_new_levels = TRUE, sample_new_levels = "gaussian",
    ndraws = ndraws_pred, seed = 999,
    transform = \(x) betaTransform(x, -10, 10, direction = "reverse", squeeze = "none")
  )
  
  # --- Example 5: frequentist GEE engine - same call signature --------------
  sex_contrast_gee <- avg_comparisons_manual(
    m3Gee, newdata = grid_test, variable = "Sex",
    nsim = 1000, seed = 999,
    transform = \(x) betaTransform(x, -10, 10, direction = "reverse", squeeze = "none")
  )
  
  # --- Example 6: parallel across levels for a factor with several levels ---
  cl <- make_prediction_cluster(m3b, n_workers = 4, extra_exports = "betaTransform")
  nation_contrast <- avg_comparisons_manual(
    m3b, newdata = grid_test, variable = "NationGeo2", cl = cl,
    re_formula = NULL, allow_new_levels = TRUE, sample_new_levels = "gaussian",
    ndraws = ndraws_pred, seed = 999,
    transform = \(x) betaTransform(x, -10, 10, direction = "reverse", squeeze = "none")
  )
  parallel::stopCluster(cl)
}