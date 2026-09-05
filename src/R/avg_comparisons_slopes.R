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
#   - Cross-contrasts (difference-in-differences across two variables) are
#     supported via avg_comparisons_cross_manual() (v0.11, generalized in
#     v0.12): variable1 can have any number of levels, compared in
#     "reference" (default, matching the model's own dummy coding - one
#     row per surviving interaction coefficient) or "pairwise" style;
#     variable2 stays fixed at exactly the 2 levels defining the contrast
#     being tracked. Both variables having more than 2 levels simultaneously
#     is not supported - see that function's docstring for why this is a
#     deliberate scope decision, not just an unfinished feature.
#   - Simple effects (a 2-level contrast computed independently within each
#     level of another factor - the natural quantity for a plot of "how
#     does this contrast vary across that factor's levels", DISTINCT from
#     the interaction contrasts above, which are differences BETWEEN
#     simple effects) are supported via avg_simple_effects_manual() (v0.13).
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
#   runs.
#
#   In non-interactive sessions (Rscript, knitr, a scheduled/CI job, etc.)
#   there is no user available to answer that prompt, so the behaviour is
#   controlled separately via options(mfx_manual.noninteractive_action)
#   or the noninteractive_action argument:
#     "error" (default) - stop with an informative message. Safe default:
#       an unattended run silently proceeding past a size threshold is
#       exactly the scenario that caused system instability earlier in
#       this project.
#     "warn" - emit an immediate, visible warning and proceed anyway.
#       Appropriate once you've reviewed the sizes involved for a given
#       scheduled pipeline and are comfortable letting it run unattended.
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
#   v0.8  Two real usability problems with v0.7's fix, both raised
#         directly: (1) a FIXED dry-run size (100, 2000 rows) pays real
#         computation cost on every single call, even when calling the
#         same model with the same settings many times in a row (the
#         actual typical usage pattern - many avg_comparisons_manual()
#         calls against one model, testing different factors); (2) a fixed
#         absolute dry-run size can exceed the size of the real grid for
#         small calls, which is nonsensical. Fixed both: dry-run sample
#         sizes are now chosen as a small percentage of nrow(newdata)
#         (capped), never a meaningful fraction of the real computation;
#         and the fitted (fixed_overhead_mb, per_row_mb) calibration is now
#         CACHED per model class/family + memory-relevant settings
#         (re_formula, allow_new_levels, sample_new_levels, ndraws/nsim) in
#         .mfx_size_cache, so the dry run only actually runs once per
#         session per unique model/settings combination - every subsequent
#         call is pure arithmetic against the cached line. dry_run_sample_sizes
#         argument removed (sizing is now automatic); force_recalibrate
#         argument and reset_size_cache() added for when you want to
#         bypass or clear the cache. Extrapolating the cached line far
#         beyond the grid size it was calibrated on is still a real
#         limitation - not eliminated, just no longer the common case.
#   v0.9  BUG: avg_predictions_manual() crashed with a data.table
#         `[.data.table` error when `newdata` was a data.table rather than
#         a plain data.frame (e.g. inherited from
#         marginaleffects::datagrid(), which returns one) - data.table
#         overrides base-R single-bracket semantics, so `newdata[by]`
#         (meant as column selection) and `newdata[i, by, drop=FALSE]`
#         were being reinterpreted as join/filter operations instead.
#         Fixed by coercing newdata to a plain data.frame with
#         as.data.frame() at the top of every public function
#         (avg_predictions_manual, avg_comparisons_manual,
#         avg_slopes_manual) and inside .estimate_true_size_gb(), so
#         behaviour no longer depends on what class of object upstream
#         code happens to hand this file.
#   v0.10 Non-interactive size-threshold handling was previously fixed
#         (always stop with an error) regardless of context, which is a
#         safe default but too rigid for a scheduled/CI pipeline where
#         you've already reviewed the sizes involved. Added a
#         noninteractive_action argument on all three public functions
#         and options(mfx_manual.noninteractive_action) as its global
#         default, choosing between "error" (default, unchanged behaviour)
#         and "warn" (emit an immediate, visible warning and proceed
#         automatically instead of stopping). Interactive-session behaviour
#         (the numbered Yes/No menu prompt) is unchanged either way.
#   v0.11 Added avg_comparisons_cross_manual() for 2x2 difference-in-
#         differences (interaction) contrasts - e.g. "does the UASType gap
#         differ between AmbientEnvCore levels" - directly on the response
#         scale with a proper credible interval, rather than reading an
#         interaction coefficient off a logit-scale summary table. Reuses
#         the same .dispatch_draws_multi/.estimate_true_size_gb/
#         .confirm_large_computation machinery as avg_comparisons_manual();
#         a new, separate function rather than a modification to it, so
#         the existing tested single-variable path is untouched. Restricted
#         to exactly 2 levels per variable by design (see docstring); does
#         not attempt general NxM cross-contrasts.
#   v0.12 The strict 2x2 restriction in v0.11, on reflection, was more
#         restrictive than necessary for the common case: a k-level
#         factor's interaction with another factor is already
#         parameterized by the model itself as k-1 coefficients relative
#         to a reference level, not C(k,2) pairwise combinations - so
#         examining "the interaction across every level" of a >2-level
#         variable1 IS a well-defined request once framed that way, it
#         just wasn't what v0.11 supported. Generalized: variable1 can now
#         take any number of levels (default: all observed), compared via
#         the new comparison_type1 argument in "reference" (default - one
#         row per level vs. levels1[1], directly matching the model's own
#         dummy-coded interaction coefficients one-to-one) or "pairwise"
#         style. variable2 remains fixed at exactly 2 levels - the
#         contrast being tracked across variable1 - since letting BOTH
#         vary reintroduces the original ambiguity this function exists to
#         avoid. Also documented explicitly (not changed - already true
#         since .dispatch_draws()'s seed=1 default was first written) that
#         every leg of every multi-prediction function in this file
#         already shares one seed automatically via how ... is forwarded,
#         so correct new-level-noise cancellation does not require the
#         caller to specify a seed manually.
#   v0.13 avg_comparisons_cross_manual()'s interaction contrasts (v0.11/
#         v0.12) answer "how does the gap CHANGE relative to a reference
#         level" but were mistaken for, and cannot substitute for, "what
#         IS the gap at each level" - a plot of a 2-level contrast varying
#         across another factor's levels needs the latter (simple effects),
#         not the former (interaction contrasts). Added
#         avg_simple_effects_manual() to compute exactly that: the 2-level
#         contrast independently within each level of a second factor, as
#         its own set of rows, none expressed relative to a reference.
#         Related to the interaction contrasts by simple subtraction
#         (interaction_contrast(level) = simple_effect(level) -
#         simple_effect(reference)) but genuinely a different output shape
#         serving a different purpose - a new function, not a mode switch
#         on avg_comparisons_cross_manual(), to keep each function's output
#         unambiguous about which of the two quantities it returns.
#   v0.14 Three issues reported from real use of avg_simple_effects_manual()
#         on an 8-leg call (4 within_levels x 2 levels):
#         (1) BUG: the memory estimate (48.3 GB) was far above the actual
#         peak. Cause: .estimate_true_size_gb() multiplied the WHOLE
#         per-leg estimate (which includes large transient, family-specific
#         computation overhead - see v0.6/v0.7) by n_conditions, as if
#         every leg's worst-case transient moment occurred simultaneously.
#         In reality only one leg computes at a time; finished legs just
#         hold their much smaller final matrices. Fixed: total is now one
#         leg's full (transient-inclusive) peak PLUS (n_conditions-1) times
#         the plain, un-inflated size of one retained matrix - see
#         .total_from_leg_peak_mb() inside .estimate_true_size_gb().
#         (2) Total system memory stayed elevated after a large call
#         completed, requiring a manual gc() to reclaim it (expected R/OS
#         behaviour - see prior discussion - but avoidable). Added
#         on.exit(gc(full=TRUE), add=TRUE) to all five public functions so
#         this happens automatically before control returns to the caller.
#         (3) p.value showing exactly 0 when every posterior/draw agreed on
#         sign (pd = 1.0 exactly) - a real Monte Carlo floor artifact, not
#         a calculation error, but "p = 0" overstates precision a finite
#         draws sample can support. Added .pd_to_pvalue(), which floors the
#         reported p-value at 2/ndraws_used (the smallest non-zero value
#         that sample size could have resolved) and is now used everywhere
#         p.value is computed (avg_comparisons_manual, avg_comparisons_cross_manual,
#         avg_simple_effects_manual, avg_slopes_manual). A floored value
#         should be read/reported as "p < floor", not "p = floor".
#   v0.15 Upgraded on.exit() memory cleanup on all five public functions
#         from a single gc(full=TRUE) call to a double pass - single-pass
#         was observed to leave Task Manager showing elevated usage on at
#         least one system even with it already active (v0.14), until a
#         manually-triggered second call visibly dropped it - plausibly a
#         Windows virtual-memory-manager timing quirk rather than R itself
#         failing to mark the memory as free. Also added plain `level1`/
#         `reference1` columns to avg_comparisons_cross_manual()'s output
#         (renamed to the actual variable1 name, e.g. "AmbientEnvCore"),
#         alongside the existing formatted `contrast` string, so results
#         can be mapped directly to a plot axis without string parsing.
#   v0.16 A hand-typed level ("Streetsidesquare") that didn't exactly
#         match the model's actual factor level ("Streetside square", with
#         a space) only surfaced deep inside brms's own validate_newdata(),
#         AFTER the size-check prompt had already been shown, confirmed,
#         and real computation had started - wasting the confirmation step
#         and whatever compute ran before brms's internal check fired.
#         Added .validate_levels(), called immediately after levels/
#         levels1/levels2/within_levels are determined (whether defaulted
#         or user-supplied) in avg_comparisons_manual(),
#         avg_comparisons_cross_manual(), and avg_simple_effects_manual() -
#         before any size check or dispatch - so a typo/spacing mismatch
#         now fails instantly and for free, with a message listing the
#         actual observed levels.
#   v0.17 The hand-rolled grid-construction code (base_unique/design_cells/
#         demog_profiles/crossing pattern, refined over several turns to
#         fix a labeling inconsistency where a single "new_id_k" label was
#         being crossed against many different real demographic profiles)
#         had grown complex enough to warrant its own function rather than
#         copy-pasted, slightly-different versions per script. Added
#         build_marginal_grid(), generalizing that pattern: separates
#         within-subject design cells from between-subject demographic
#         profiles explicitly, supports independent simulate_id/
#         simulate_stim toggles (covering both the marginal-estimation use
#         case and the conditional/calibration-plot use case from the same
#         function), counterfactual variable expansion (e.g. an LAE sweep),
#         fixed_at_mean/fixed_at_mode nuisance covariates, and validates
#         that every predictor the model actually requires is accounted
#         for via insight::find_predictors() before returning anything.
#         Deliberately not a wrapper around marginaleffects::datagrid() or
#         emmeans' reference grid - see the function's docstring for the
#         three specific reasons neither tool can do this.
#   v0.18 BUG (design regression): build_marginal_grid()'s demographic
#         profile assignment used sample(replace=TRUE) to stretch the real
#         profiles up to n_new synthetic ones - copied from working code a
#         few turns earlier without questioning whether resampling was
#         actually necessary. It wasn't: the original base_unique approach
#         was fully deterministic (every real subject's own profile,
#         always included), and introducing random resampling here added
#         an unnecessary second source of Monte Carlo noise on top of the
#         one that IS required (sample_new_levels="gaussian"), with no
#         benefit - and a real downside: sample(replace=TRUE) can, by
#         chance (coupon-collector problem), omit a real profile entirely,
#         and gives an unevenly-weighted allocation rather than an exactly
#         even one. Fixed: demog_assignment now uses
#         rep_len(seq_len(nrow(demog_profiles)), n_new) - deterministic
#         cycling through all real profiles, guaranteeing every one
#         appears, with as-even-as-possible representation, and no seed
#         needed. The now-unused `seed` argument was removed (not left as
#         a silently-ignored parameter); passing it triggers an explicit
#         warning via `...` instead.
#   v0.19 The v0.18 fix's own docstring flagged, but didn't act on, a real
#         edge case: if n_new < nrow(demog_profiles), deterministic cycling
#         degrades to silently taking the first n_new real profiles in
#         whatever order they appear in the data - an arbitrary row-order
#         dependency, not a meaningful subsample. Added an explicit warning
#         for this case rather than leaving it as a documented-but-silent
#         caveat, since n_new is under direct user control and this
#         condition is trivial to detect at call time.
#   v0.20 Only one input-combination hazard (within_vars/between_vars
#         overlap) was checked; asked to audit systematically for others.
#         Added .validate_grid_inputs(), covering: id_var == stim_var;
#         id_var/stim_var also listed in within_vars/between_vars (breaks
#         the design-cell/demographic-profile separation); a
#         counterfactual_vars name also appearing anywhere else (the sweep
#         assumes the column doesn't already exist, or the cross-join
#         collides with it); an empty counterfactual_vars entry (silently
#         produces a 0-row grid, the same failure mode as the datagrid()
#         bug from earlier in this project); fixed_at_mean and
#         fixed_at_mode both listing the same variable (contradictory);
#         id_var/stim_var listed in fixed_at_mean/fixed_at_mode (mean() on
#         a character ID column returns NA); fixed_at_mean applied to a
#         non-numeric column (same NA problem); n_new < 1. All of the
#         above STOP. Separately, fixed_at_mean/fixed_at_mode also
#         appearing in within_vars/between_vars WARNS rather than stopping
#         (per explicit direction) - it's recoverable, just means the
#         fixed step overwrites the naturally-varying values with a
#         constant, which the warning states explicitly.
#   v0.21 BUG: the counterfactual_vars expansion loop passed
#         setNames(list(values), v) - a bare named list - as a positional
#         argument to tidyr::crossing(), which does not reliably give
#         crossing() an unambiguous column name+values pair. In practice
#         this produced a grid that was silently NOT expanded at all
#         (confirmed: a 52-design-cell x 150-profile grid with a 13-value
#         counterfactual sweep returned 7,800 rows = 52*150, i.e. exactly
#         the un-expanded count) while the separate verbose-message text
#         incorrectly reported "x 13 counterfactual level(s)" regardless,
#         since it computed that number independently of what actually
#         happened to the grid. Fixed by constructing a proper one-column
#         data.frame for each counterfactual variable before crossing.
#         Also added a self-check immediately after the expansion loop:
#         if nrow(grid) doesn't exactly equal n_before x (product of
#         counterfactual level counts), the function now STOPS with a
#         clear message rather than silently returning a wrong-sized grid
#         - this exact bug would have been caught immediately by that
#         check had it existed, rather than requiring the user to notice
#         the arithmetic didn't add up after the fact.
#   v0.22 Asked for a full audit rather than just confirming v0.21's fix.
#         Found two more real issues: (1) fixed_at_mode built its
#         replacement value via names(table(...)), which is ALWAYS
#         character - applying fixed_at_mode to a genuinely numeric
#         (but discrete) variable would silently turn that grid column
#         into character (e.g. "3" instead of 3), which brms::posterior_epred()
#         would not treat as the numeric predictor it needs. Fixed by
#         coercing the modal value back to the original column's type
#         (numeric/logical/factor). (2) The v0.21 self-check and the
#         verbose message each independently recomputed the expected
#         counterfactual row-count multiplier via raw length() -
#         exactly the same class of bug as v0.21 itself (a displayed/
#         checked number computed separately from what actually happens),
#         just not yet triggered: tidyr::crossing() deduplicates each
#         input's values before crossing, so a counterfactual vector with
#         any duplicate values would make the v0.21 self-check itself
#         throw a false-positive stop. Fixed by computing cf_unique_counts
#         (via length(unique(...)), matching crossing()'s real behaviour)
#         ONCE, and having both the self-check and the verbose message
#         read from that single value - making the "message disagrees
#         with reality" bug class structurally impossible here, not just
#         patched for this one instance of it.
#   v0.23 avg_comparisons_manual(comparison_type="reference") and
#         avg_simple_effects_manual() disagreed on which level their
#         shared default (levels=NULL, sort()-based) ordering treats as
#         the reference: the former follows the standard dummy-coding
#         convention (alphabetically-first level = reference, subtracted;
#         matching how the model's own coefficients are parameterized),
#         the latter computed levels[1] - levels[2] literally, giving the
#         alphabetically-first level the OPPOSITE role. Same data, same
#         default sort() ordering, opposite-signed results (confirmed:
#         "Near - Far" vs "Far - Near" for the same UASProximity contrast).
#         Fixed avg_simple_effects_manual() and avg_comparisons_cross_manual()'s
#         levels2 default to match avg_comparisons_manual()'s convention -
#         ONLY when levels/levels2 is not explicitly supplied, so existing
#         calls using explicit level order (e.g. levels = c("T150", "H520"))
#         are completely unaffected and their prior results still hold.
#         Also added an explicit console message stating the computed
#         direction whenever levels default rather than being supplied,
#         as a permanent safety net against this class of directionality
#         confusion recurring.
#   v0.24 Asked whether build_marginal_grid() has an analogous bug to
#         v0.23's. It doesn't structurally - grid construction only
#         enumerates/crosses levels, it never computes a hi-lo difference,
#         so there's no "which level is silently the reference" question
#         to get wrong. But the same FLAVOR of issue (a silent, sort()-
#         order-driven default) does exist in fixed_at_mode: which.max()
#         on a tied table() silently returns the alphabetically-first
#         level with no indication of the tie. Added an explicit warning
#         when fixed_at_mode's modal value is tied, naming the tied
#         levels and which one was picked. design_cells' "first real value
#         encountered" and demog_profiles' cycling order both also depend
#         on raw row order, but are a lower-severity category than v0.23:
#         they pick which equally-valid representative value fills a slot,
#         not which quantity gets computed/reported - not fixed, since
#         there's no direction to get backwards there.
#   v0.25 avg_comparisons_manual(), avg_comparisons_cross_manual(), and
#         avg_simple_effects_manual() default their level arguments via
#         sort(unique(as.character(newdata[[variable]]))) - always
#         character, even for a genuinely numeric variable like UASEvents
#         (which feeds I(log10(UASEvents)) in the model formula). Writing
#         a character level value into an existing numeric data.frame
#         column (nd[[variable]] <- lv) replaces the WHOLE column with a
#         coerced character vector, not just that one value - silently
#         corrupting every row of the grid, not just the ones intended to
#         differ. This was worked around manually for a UASEvents contrast
#         (direct posterior_epred() calls, bypassing these functions
#         entirely) rather than left as a standing hazard. Fixed properly:
#         added .coerce_levels_to_type(), applied to levels/levels1/
#         levels2/within_levels in all three functions immediately after
#         they're determined (whether defaulted or user-supplied),
#         coercing back to the target column's actual type (numeric/
#         logical/factor/character) before anything is written into
#         newdata. Same fix family as fixed_at_mode's type-preservation
#         fix in build_marginal_grid() (v0.22). Confirmed the internal
#         key() lookup functions (paste()-based) remain self-consistent
#         regardless of level type, since both grid construction and key
#         construction draw from the same, already-coerced level vectors.
#   v0.26 Recurring source of confusion (second occurrence, after v0.23's
#         Near/Far mixup): avg_comparisons_manual()'s default
#         comparison_type = "pairwise" computes levels[1] - levels[2]
#         LITERALLY, regardless of which level the model treats as its
#         reference - opposite of regression convention, where the
#         reference (however it's positioned in levels()) is always the
#         SUBTRACTED term. comparison_type = "reference" already does
#         match regression convention (nonreference - levels[1], treating
#         levels[1] as reference) - but achieving a regression-matching
#         sign requires listing the reference LAST under "pairwise" and
#         FIRST under "reference": two modes, opposite required level
#         ordering, same resulting sign. Not changing "pairwise"'s literal
#         behavior (would silently flip the sign of every existing
#         validated result using explicit levels, e.g. the T150-H520
#         examples throughout this project) - instead made the
#         confirmation message introduced in v0.23 UNCONDITIONAL in all
#         three comparison functions (previously only fired when levels
#         were defaulted, not when explicitly supplied - exactly the case
#         that caused this confusion). Every call now prints exactly which
#         subtraction(s) it computed before returning, regardless of how
#         comparison_type/levels were determined.
#   v0.27 Attempted to solve a 3-level sign-ordering request (all pairwise
#         UASEvents contrasts positive) by choosing a clever levels= input
#         order exploiting combn()'s enumeration pattern - got the specific
#         ordering wrong AND, on further analysis, found the whole approach
#         is structurally incapable of the general case: combn(levels, 2)
#         always makes the first input element lead 2 of 3 pairs and the
#         second lead the remaining 1, so no input reordering can give two
#         DIFFERENT levels each a turn leading multiple pairs - some
#         "make every row positive" requests are simply unreachable by
#         permuting input order, not just hard to find by hand. Added
#         auto_orient (default FALSE, non-breaking) to avg_comparisons_manual():
#         when TRUE, each row's sign/label is chosen from the ACTUAL
#         computed estimate (flip to positive, swap the label's two level
#         names together) rather than guessed from input order - correct
#         by construction regardless of how many levels or pairs are
#         involved. pd/p.value/ROPE are computed from the already-oriented
#         draws, so they remain internally consistent with the reported
#         sign. Scope note: NOT added to avg_simple_effects_manual() or
#         avg_comparisons_cross_manual() - their rows represent the SAME
#         named contrast (e.g. "T150-H520") repeated across different
#         strata, where per-row auto-orientation would let the label
#         silently flip identity row-to-row and undermine the point of
#         having one consistent column to compare across strata.
#   v0.28 BUG (serious, introduced by v0.25's own fix): v0.25 coerced the
#         WHOLE `levels`/`levels1`/`levels2`/`within_levels` vector to the
#         target column's real type immediately after determining it - for
#         a genuinely numeric column (UASEvents) this was correct and
#         necessary, but for a FACTOR column (e.g. AmbientEnvCore) it
#         turned `levels` into an actual factor object, which then got fed
#         into generic list/matrix-building code (combn(), cbind(),
#         lapply()) that does not reliably preserve factor labels through
#         those operations. Confirmed two distinct failure modes from this:
#         (1) avg_comparisons_cross_manual()'s "reference" mode
#         (cbind(levels1[-1], levels1[1])) silently dropped to underlying
#         integer codes, producing output like AmbientEnvCore="2" instead
#         of "Residential"; (2) avg_comparisons_manual()'s "pairwise" mode
#         (combn()) produced internally MISMATCHED pairs/labels/values -
#         confirmed by hand-checking m0d_env_contrasts against
#         m0d_env_means and finding estimates did not correspond to their
#         own printed labels. Fixed properly: `levels`/`levels1`/`levels2`/
#         `within_levels` are now kept as plain CHARACTER throughout all
#         pairing/labeling/naming logic in all three functions (as before
#         v0.25) - .coerce_levels_to_type() is now called only at the one
#         narrow point in each function where an individual SCALAR value
#         is written into newdata for prediction, which is the only place
#         type-safety was ever actually required. This preserves v0.25's
#         original fix (numeric columns like UASEvents still get correctly
#         typed values written into newdata) without the collateral damage
#         to factor-typed variables.
#   v0.29 A raw pd (e.g. "0.84") has no pre-existing social calibration and
#         reads as fairly high to most audiences even though it corresponds
#         to p~.33 - weak evidence conventionally. Added ndraws_used as a
#         returned column on every avg_comparisons_manual()/
#         avg_comparisons_cross_manual()/avg_simple_effects_manual()/
#         avg_slopes_manual() row, and format_pd_p_label() to build a
#         combined "pd = .., p = .." (or "p <= .." when the raw p-value is
#         at its Monte Carlo resolution floor) label from it - shows both
#         scales together so a reader can cross-check one against the
#         other, and reads ndraws_used directly from the data rather than
#         requiring a separately-remembered ndraws value (correct even if
#         a table combines rows from calls with different ndraws). Handles
#         an optional multiplicity-adjusted p-value column too; the
#         floor("<=") determination is still driven by the RAW p-value's
#         floor status even when the adjusted number is what's displayed,
#         which is a valid (via Holm's monotonicity) if not maximally
#         tight bound - documented as a deliberate simplification, not a
#         full symbolic re-derivation of adjustment applied to an
#         interval-valued input.
#   v0.30 Two additions arising from building the figures in practice:
#         (1) add_pairwise_brackets(), a plotting helper wrapping
#         ggpubr::stat_pvalue_manual() to annotate a means/estimates plot
#         with each pairwise contrast's own credible interval (not pd/p -
#         see v0.29's rationale for why). Accepts either
#         avg_comparisons_cross_manual()'s output directly, or
#         avg_comparisons_manual()'s (parsing its "contrast" string, since
#         that function has no separate group1/group2-style columns).
#         Deliberately requires a NON-coord_flip()'d base plot -
#         stat_pvalue_manual() has a confirmed, hard-to-work-around
#         interaction problem with coord_flip() (correct bracket position
#         with illegible sideways text, or legible text with misaligned
#         brackets - never both correct at once). (2) Documented, in
#         avg_comparisons_manual()'s own docstring, that "reference" and
#         "pairwise" comparison_type are MIRROR IMAGES of each other for
#         the common case of exactly 2 levels - confirmed directly:
#         levels=c("H520","T150") gives "T150-H520" under "reference" but
#         "H520-T150" under "pairwise", for the same input. Not a bug -
#         each mode is organized around a different principle (reference:
#         levels[1] is always the baseline/subtrahend; pairwise: levels[1]
#         is always the minuend) that necessarily oppose in sign once
#         there's only one comparison for both modes to produce.
#   v0.31 BUG: make_prediction_cluster()'s parallel workers failed with
#         "could not find function '.dispatch_draws'" on first real use
#         with cl= - clusterExport() only exported "model" plus whatever
#         the caller listed in extra_exports, never the internal
#         .dispatch_draws() function that every worker actually calls
#         (via .dispatch_draws_multi()'s parLapply() path). Worker
#         processes are fresh R sessions with nothing from the sourced
#         script loaded except what's explicitly exported - this was a
#         structural gap in the function itself, not something callers
#         could fix by remembering to list it in extra_exports each time.
#         Confirmed .dispatch_draws() has no further chained dependencies
#         on other internal helpers (only calls brms::posterior_epred()/
#         MASS::mvrnorm()/stats:: functions directly), so exporting just
#         this one name resolves it completely. Fixed: ".dispatch_draws"
#         is now added to the base export list unconditionally, not left
#         for the caller to supply via extra_exports.
#   v0.32 add_pairwise_brackets() failed with a "discrete value supplied
#         to continuous scale" error on a numeric x-axis (UASEvents) -
#         group1/group2 (character, from string-splitting a "contrast"
#         label) were always routed through an x_order-based integer
#         position lookup built for discrete axes, which conflicts with
#         stat_pvalue_manual() when the base plot's x-scale is actually
#         continuous. Fixed: now checks whether means_df[[x_var]] is
#         numeric - if so, group1/group2 are coerced directly to that
#         same numeric scale (no position-index translation needed, since
#         a numeric value already is its own x-position) and x_order is
#         ignored (with a message) if supplied. x_order is also now
#         optional for the discrete case: defaults to
#         levels(means_df[[x_var]]) if already a factor, otherwise
#         sort(unique(as.character(...))) - matching the default-level
#         convention used elsewhere in this file, rather than requiring
#         it as a mandatory argument every time.
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
#'   calling environment) to export to workers, e.g. custom transform helpers
#'   you reference in a `transform=` argument. Note: `.dispatch_draws` (the
#'   internal function every worker actually calls) is exported automatically
#'   below - you do NOT need to add it here. Worker processes are fresh R
#'   sessions with nothing from your script loaded except what's explicitly
#'   exported; omitting an internal helper this file depends on produces
#'   "could not find function" errors on the workers, confirmed in practice.
#' @return A cluster object for use as the `cl` argument below. Remember to
#'   `parallel::stopCluster(cl)` when done.
make_prediction_cluster <- function(model, n_workers = parallel::detectCores() - 1,
                                    extra_exports = character(0)) {
  cl <- parallel::makeCluster(max(1, n_workers))
  parallel::clusterEvalQ(cl, { library(brms); library(MASS) })
  parallel::clusterExport(cl, c("model", ".dispatch_draws", extra_exports), envir = environment())
  cl
}


# -----------------------------------------------------------------------------
# build_marginal_grid()
# -----------------------------------------------------------------------------

#' Validate build_marginal_grid()'s inputs before any grid construction runs
#'
#' Hard stops for combinations that would structurally break the grid
#' (wrong-typed columns, colliding cross-joins, or a silent empty result -
#' the same failure mode as the datagrid() 0-row bug from earlier in this
#' project, reachable a different way via an empty counterfactual_vars
#' entry). Warnings for combinations that are recoverable but almost
#' certainly not what was intended (a fixed_at_mean/fixed_at_mode variable
#' also listed as varying in within_vars/between_vars - the fixed step
#' runs AFTER crossing and silently overwrites the naturally-varying
#' values with a constant).
.validate_grid_inputs <- function(model_data, within_vars, between_vars, id_var, stim_var,
                                  counterfactual_vars, fixed_at_mean, fixed_at_mode,
                                  n_new, simulate_id) {
  cf_vars <- names(counterfactual_vars)
  
  if (identical(id_var, stim_var)) {
    stop("id_var and stim_var cannot be the same column ('", id_var, "').", call. = FALSE)
  }
  
  overlap_wb <- intersect(within_vars, between_vars)
  if (length(overlap_wb) > 0) {
    stop("Variable(s) listed in both within_vars and between_vars: ",
         paste(sQuote(overlap_wb), collapse = ", "), call. = FALSE)
  }
  
  id_stim_in_wb <- intersect(c(id_var, stim_var), c(within_vars, between_vars))
  if (length(id_stim_in_wb) > 0) {
    stop(
      "id_var/stim_var ('", id_var, "'/'", stim_var, "') must not also appear in within_vars or ",
      "between_vars - found: ", paste(sQuote(id_stim_in_wb), collapse = ", "),
      ". These are handled separately by the id_var/stim_var/simulate_id/simulate_stim ",
      "machinery; including them elsewhere breaks the design-cell/demographic-profile ",
      "separation this function relies on.",
      call. = FALSE
    )
  }
  
  cf_overlap <- intersect(cf_vars, c(within_vars, between_vars, id_var, stim_var,
                                     fixed_at_mean, fixed_at_mode))
  if (length(cf_overlap) > 0) {
    stop(
      "counterfactual_vars name(s) also appear elsewhere (within_vars/between_vars/id_var/",
      "stim_var/fixed_at_mean/fixed_at_mode): ", paste(sQuote(cf_overlap), collapse = ", "),
      ". A counterfactual sweep variable must not already exist as a column by the time the ",
      "sweep runs, or the resulting cross-join collides with (or silently duplicates) the ",
      "existing column. Remove it from wherever else it's listed.",
      call. = FALSE
    )
  }
  
  empty_cf <- cf_vars[lengths(counterfactual_vars) == 0]
  if (length(empty_cf) > 0) {
    stop(
      "counterfactual_vars entry has zero values, which would silently produce an empty ",
      "(0-row) grid - the same failure mode as the datagrid() 0-row bug from earlier in this ",
      "project: ", paste(sQuote(empty_cf), collapse = ", "), call. = FALSE
    )
  }
  
  mean_mode_overlap <- intersect(fixed_at_mean, fixed_at_mode)
  if (length(mean_mode_overlap) > 0) {
    stop(
      "Variable(s) listed in both fixed_at_mean and fixed_at_mode (contradictory - pick one): ",
      paste(sQuote(mean_mode_overlap), collapse = ", "), call. = FALSE
    )
  }
  
  id_stim_in_fixed <- intersect(c(id_var, stim_var), c(fixed_at_mean, fixed_at_mode))
  if (length(id_stim_in_fixed) > 0) {
    stop(
      "id_var/stim_var must not appear in fixed_at_mean/fixed_at_mode - found: ",
      paste(sQuote(id_stim_in_fixed), collapse = ", "),
      ". These columns are managed by simulate_id/simulate_stim; averaging or mode-ing a ",
      "character/factor ID column is not meaningful and would corrupt it (mean() on a ",
      "non-numeric column returns NA in base R without necessarily stopping).",
      call. = FALSE
    )
  }
  
  fixed_at_mean_in_data <- intersect(fixed_at_mean, names(model_data))
  non_numeric_means <- fixed_at_mean_in_data[
    !vapply(fixed_at_mean_in_data, function(v) is.numeric(model_data[[v]]), logical(1))
  ]
  if (length(non_numeric_means) > 0) {
    stop(
      "fixed_at_mean variable(s) are not numeric in the data, so mean() is not meaningful: ",
      paste(sQuote(non_numeric_means), collapse = ", "),
      ". Did you mean to list these under fixed_at_mode instead?",
      call. = FALSE
    )
  }
  
  if (simulate_id && n_new < 1) {
    stop("n_new must be at least 1 when simulate_id = TRUE (got ", n_new, ").", call. = FALSE)
  }
  
  # ---- Soft warnings: likely-unintended but recoverable overlaps ---------
  mean_wb_overlap <- intersect(fixed_at_mean, c(within_vars, between_vars))
  if (length(mean_wb_overlap) > 0) {
    warning(
      "fixed_at_mean variable(s) also appear in within_vars/between_vars: ",
      paste(sQuote(mean_wb_overlap), collapse = ", "),
      ". fixed_at_mean is applied AFTER the within/between crossing and will OVERWRITE these ",
      "columns' naturally-varying values with a single constant mean - if that's not what you ",
      "intended, remove the variable from fixed_at_mean.",
      call. = FALSE
    )
  }
  
  mode_wb_overlap <- intersect(fixed_at_mode, c(within_vars, between_vars))
  if (length(mode_wb_overlap) > 0) {
    warning(
      "fixed_at_mode variable(s) also appear in within_vars/between_vars: ",
      paste(sQuote(mode_wb_overlap), collapse = ", "),
      ". fixed_at_mode is applied AFTER the within/between crossing and will OVERWRITE these ",
      "columns' naturally-varying values with a single constant mode - if that's not what you ",
      "intended, remove the variable from fixed_at_mode.",
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}


#' Build a marginalization grid for hierarchical models with crossed random
#' effects, separating within-subject (design-cell) and between-subject
#' (demographic) variation explicitly
#'
#' Genuinely different from marginaleffects::datagrid() and emmeans'
#' reference grid, not a reimplementation of either - see conversation
#' notes for the three specific reasons: (1) datagrid()'s factor-level
#' validation rejects synthetic ID/StimID labels outright, so it cannot
#' build a grid compatible with allow_new_levels/sample_new_levels at all;
#' (2) emmeans defaults to a balanced reference grid at means/modes (the
#' "average case" approach), not the "observed value" approach (Hanmer &
#' Kalkan 2013) this project has used throughout, which requires crossing
#' against the actual empirical distribution of demographic profiles;
#' (3) neither has a notion of "within-subject design cell" vs
#' "between-subject profile" as independently-crossed entities, which is
#' what a repeated-measures factorial design with crossed random effects
#' actually needs.
#'
#' DESIGN CELLS (within_vars) are deduplicated deterministically - one row
#' per unique combination, keeping whichever real stim_var/covariate value
#' first co-occurred with it (a representative, not arbitrary, choice - see
#' the .keep_all= pattern this generalizes). DEMOGRAPHIC PROFILES
#' (between_vars) are, if simulate_id=TRUE, RESAMPLED WITH REPLACEMENT to
#' n_new - preserving the observed distribution rather than enumerating
#' every unique profile once (which would not achieve genuine Monte Carlo
#' integration over repeated new-subject draws the way resampling does).
#' The two are then crossed, so every simulated demographic profile is
#' paired with every design cell - matching the structure validated by the
#' convergence-sweep pilots earlier in this project.
#'
#' @param model brmsfit or glmgee object.
#' @param within_vars Character vector of within-subject/design-cell
#'   variable names (e.g. UASType, UASOperation, AmbientEnvCore).
#' @param between_vars Character vector of between-subject/demographic
#'   variable names (e.g. Sex, Age, NationGeo2).
#' @param id_var,stim_var Names of the subject and stimulus grouping
#'   columns. Default "ID"/"StimID".
#' @param simulate_id,simulate_stim Whether to replace id_var/stim_var
#'   with fresh synthetic labels (for use with allow_new_levels=TRUE,
#'   sample_new_levels="gaussian" downstream) or keep the real, existing
#'   values (for a conditional/existing-levels grid - e.g. the calibration
#'   plot use case from earlier, where simulate_stim = FALSE is correct).
#' @param n_new Number of synthetic profiles to generate when
#'   simulate_id = TRUE. Ignored (all real IDs kept) when FALSE. Intended
#'   to comfortably exceed the number of real profiles; if it doesn't, a
#'   warning is issued (see below).
#' @param counterfactual_vars Optional named list of variable -> vector of
#'   values to additionally cross the grid against (e.g.
#'   list(UASLAEMaxLRScl = seq(-11, 12, length.out = 11)) for an LAE
#'   sweep). Each entry multiplies the grid size by length(values) - a
#'   message is printed if the resulting expansion exceeds 5x, as an early
#'   warning before this feeds into the separate, automatic memory size
#'   check in avg_predictions_manual() etc.
#' @param fixed_at_mean Character vector of numeric covariate names to
#'   hold at their sample mean (e.g. "TrialNumberScl").
#' @param fixed_at_mode Character vector of categorical covariate names to
#'   hold at their sample mode (most frequent observed level).
#' @param newdata Optional data.frame to use instead of
#'   insight::get_data(model) - e.g. a pre-filtered subset.
#' @param verbose If TRUE (default), prints a short summary of what was
#'   built (row count, unique design cells, unique profiles, whether
#'   ID/StimID were simulated).
#' @param ... Not otherwise used. Present so that a leftover `seed=`
#'   argument from before v0.18 (when demographic profile assignment was
#'   randomly resampled, not deterministic) triggers an explicit warning
#'   rather than silently being accepted and ignored.
#' @return A data.frame ready to pass as `newdata` to
#'   avg_predictions_manual(), avg_comparisons_manual(),
#'   avg_comparisons_cross_manual(), avg_simple_effects_manual(), or
#'   avg_slopes_manual().
build_marginal_grid <- function(model, within_vars, between_vars,
                                id_var = "ID", stim_var = "StimID",
                                simulate_id = TRUE, simulate_stim = TRUE,
                                n_new = 150, counterfactual_vars = NULL,
                                fixed_at_mean = NULL, fixed_at_mode = NULL,
                                newdata = NULL, verbose = TRUE, ...) {
  dots <- list(...)
  if ("seed" %in% names(dots)) {
    warning(
      "build_marginal_grid()'s `seed` argument was removed in v0.18: demographic profile ",
      "assignment is now deterministic (cycles through all real profiles) rather than randomly ",
      "resampled, so no seed is needed or used here. This argument is being ignored.",
      call. = FALSE
    )
  }
  model_data <- if (!is.null(newdata)) as.data.frame(newdata) else as.data.frame(insight::get_data(model))
  
  # ---- Input validation --------------------------------------------------
  .validate_grid_inputs(model_data, within_vars, between_vars, id_var, stim_var,
                        counterfactual_vars, fixed_at_mean, fixed_at_mode, n_new, simulate_id)
  
  all_named <- unique(c(within_vars, between_vars, id_var, stim_var,
                        names(counterfactual_vars), fixed_at_mean, fixed_at_mode))
  missing_from_data <- setdiff(all_named, names(model_data))
  if (length(missing_from_data) > 0) {
    stop("The following variables are not present in the model's data: ",
         paste(sQuote(missing_from_data), collapse = ", "), call. = FALSE)
  }
  
  required_predictors <- tryCatch(
    insight::find_predictors(model, effects = "all", flatten = TRUE),
    error = function(e) NULL
  )
  if (!is.null(required_predictors)) {
    still_missing <- setdiff(required_predictors, all_named)
    if (length(still_missing) > 0) {
      stop(
        "The model requires the following predictor(s), not accounted for by any of ",
        "within_vars/between_vars/id_var/stim_var/counterfactual_vars/fixed_at_mean/fixed_at_mode: ",
        paste(sQuote(still_missing), collapse = ", "),
        ".\nEvery variable the model actually uses must be covered by one of these arguments, ",
        "or prediction will fail (or silently fall back to some default you didn't intend).",
        call. = FALSE
      )
    }
  } else {
    warning(
      "Could not determine the model's required predictors via insight::find_predictors() ",
      "- proceeding without this check. This can happen for custom/forked families; verify ",
      "manually that every predictor in the model formula is covered by one of the arguments above.",
      call. = FALSE
    )
  }
  
  # ---- Design cells (within-subject / stimulus-level variables) ---------
  design_cells <- model_data[!duplicated(model_data[within_vars]),
                             c(within_vars, stim_var), drop = FALSE]
  if (simulate_stim) {
    design_cells[[stim_var]] <- paste0("new_stim_", seq_len(nrow(design_cells)))
  }
  
  # ---- Demographic profiles (between-subject variables) ------------------
  demog_profiles <- model_data[!duplicated(model_data[[id_var]]),
                               c(id_var, between_vars), drop = FALSE]
  if (simulate_id) {
    # Deterministic cycling, not random resampling: with n_new typically
    # >> number of real profiles, this guarantees every real profile
    # appears (unlike sample(replace=TRUE), which can by chance omit one),
    # gives an exactly even allocation rather than one with sampling
    # noise, and needs no seed - one fewer source of randomness in a
    # pipeline that already has to carefully manage the randomness that
    # IS required (sample_new_levels="gaussian"). Only introduces an
    # arbitrary row-order dependency if n_new < nrow(demog_profiles),
    # which is not the intended use case (n_new is meant to comfortably
    # exceed the real subject count) - warn explicitly rather than let
    # that degrade silently.
    if (n_new < nrow(demog_profiles)) {
      warning(
        "n_new (", n_new, ") is smaller than the number of real profiles in '", id_var, "' (",
        nrow(demog_profiles), "). Deterministic cycling will select only the FIRST ", n_new,
        " profiles in whatever order they happen to appear in the data - an arbitrary ",
        "dependency on row order, not a meaningful or representative subsample. If you want a ",
        "representative subset instead, pre-select or pre-shuffle the demographic profiles ",
        "yourself and pass the result via `newdata`, or use n_new >= ", nrow(demog_profiles),
        " to include every real profile.",
        call. = FALSE
      )
    }
    idx <- rep_len(seq_len(nrow(demog_profiles)), n_new)
    demog_assignment <- demog_profiles[idx, , drop = FALSE]
    demog_assignment[[id_var]] <- paste0("new_id_", seq_len(n_new))
  } else {
    demog_assignment <- demog_profiles
  }
  
  # ---- Cross design cells x demographic assignment ------------------------
  grid <- as.data.frame(tidyr::crossing(design_cells, demog_assignment))
  
  # ---- Fixed-at-mean / fixed-at-mode nuisance covariates ------------------
  for (v in fixed_at_mean) grid[[v]] <- mean(model_data[[v]], na.rm = TRUE)
  for (v in fixed_at_mode) {
    tbl <- table(model_data[[v]])
    # which.max() silently returns the FIRST index on a tie, and table()
    # orders levels alphabetically for character/factor input - so a tied
    # mode picks the alphabetically-first level with no indication anything
    # was ambiguous. Same flavor of silent sort()-driven default as the
    # reference-direction inconsistency fixed in the comparison functions
    # (v0.23); warn explicitly here for the same reason.
    n_tied <- sum(tbl == max(tbl))
    if (n_tied > 1) {
      warning(
        "fixed_at_mode: '", v, "' has ", n_tied, " levels tied for most frequent (",
        paste(sQuote(names(tbl)[tbl == max(tbl)]), collapse = ", "), "). Picking '",
        names(tbl)[which.max(tbl)], "' (alphabetically first among the tied levels) - ",
        "if this matters, set this variable's value explicitly instead of relying on fixed_at_mode.",
        call. = FALSE
      )
    }
    mode_val <- names(tbl)[which.max(tbl)]
    # names(tbl) is always character - coerce back to the original column's
    # type, or a numeric fixed_at_mode variable would silently become
    # character in the returned grid (e.g. "3" instead of 3), which
    # brms::posterior_epred() would not treat the same as a numeric predictor.
    if (is.numeric(model_data[[v]])) {
      mode_val <- as.numeric(mode_val)
    } else if (is.logical(model_data[[v]])) {
      mode_val <- as.logical(mode_val)
    } else if (is.factor(model_data[[v]])) {
      mode_val <- factor(mode_val, levels = levels(model_data[[v]]))
    }
    grid[[v]] <- mode_val
  }
  
  # ---- Counterfactual expansion (e.g. an LAE sweep) ------------------------
  # cf_unique_counts is computed ONCE here and reused below (both in the
  # self-check and in the verbose message) - the earlier bug (v0.21) was
  # exactly a case of the message computing its own row-count independently
  # of what the grid actually did, and them silently disagreeing. A single
  # shared computation makes that whole class of bug structurally
  # impossible rather than just fixing the one instance of it. Uses
  # length(unique(...)) rather than raw length(...), since
  # tidyr::crossing() deduplicates each input's values before crossing -
  # using the raw length would falsely flag a mismatch if a counterfactual
  # vector ever contained duplicate values.
  cf_unique_counts <- if (!is.null(counterfactual_vars)) {
    vapply(counterfactual_vars, function(v) length(unique(v)), integer(1))
  } else {
    integer(0)
  }
  
  if (!is.null(counterfactual_vars)) {
    n_before <- nrow(grid)
    for (v in names(counterfactual_vars)) {
      # A proper one-column data.frame, not a bare list - tidyr::crossing()
      # needs an unambiguous name+values pair per argument; a raw
      # list(name = values) object passed positionally does not reliably
      # give it that, and previously produced a grid that silently wasn't
      # expanded at all (see v0.21).
      cf_df <- data.frame(x = counterfactual_vars[[v]])
      names(cf_df) <- v
      grid <- as.data.frame(tidyr::crossing(grid, cf_df))
    }
    n_after <- nrow(grid)
    expected_n_after <- n_before * prod(cf_unique_counts)
    if (n_after != expected_n_after) {
      stop(sprintf(
        "build_marginal_grid(): counterfactual expansion produced %s rows, expected %s ",
        "(%s rows before x %s counterfactual combinations). This indicates the crossing ",
        "step is not behaving as expected - stopping rather than silently returning an ",
        "incorrectly-sized grid.",
        format(n_after, big.mark = ","), format(expected_n_after, big.mark = ","),
        format(n_before, big.mark = ","), format(prod(cf_unique_counts), big.mark = ",")
      ), call. = FALSE)
    }
    if (n_after > 5 * n_before) {
      message(sprintf(
        "build_marginal_grid(): counterfactual expansion increased the grid from %s to %s rows (%.1fx).",
        format(n_before, big.mark = ","), format(n_after, big.mark = ","), n_after / n_before
      ))
    }
  }
  
  if (verbose) {
    message(sprintf(
      "build_marginal_grid(): %s rows (%s design cells x %s profile%s%s). ID %s, StimID %s.",
      format(nrow(grid), big.mark = ","),
      format(nrow(design_cells), big.mark = ","),
      format(nrow(demog_assignment), big.mark = ","),
      if (nrow(demog_assignment) == 1) "" else "s",
      if (!is.null(counterfactual_vars)) sprintf(", x %s counterfactual level(s)",
                                                 format(prod(cf_unique_counts), big.mark = ",")) else "",
      if (simulate_id) "simulated (new levels)" else "real (existing levels)",
      if (simulate_stim) "simulated (new levels)" else "real (existing levels)"
    ))
  }
  
  grid
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

#' Convert probability of direction to a two-sided p-value, floored at the
#' Monte Carlo resolution limit of the draws actually used
#'
#' bayestestR::p_direction() computes pd as a proportion of a FINITE set of
#' draws. When every single draw falls on the same side of zero, pd = 1.0
#' exactly and p = 2(1-pd) = 0 exactly - not because the true probability
#' of a sign flip is literally zero, but because the finite Monte Carlo
#' sample didn't happen to contain a single crossing. The smallest
#' non-zero p this specific sample of draws could ever have resolved is
#' 2/ndraws_used (one single draw on the other side instead of zero), so
#' reporting a floored value there - and treating it as "p < floor", not
#' "p = 0" - is the honest way to report it. Values that aren't at the
#' floor are returned unchanged.
.pd_to_pvalue <- function(pd, ndraws_used) {
  p <- 2 * (1 - pd)
  floor_p <- 2 / ndraws_used
  ifelse(p < floor_p, floor_p, p)
}

#' Format a combined "pd, p" label for plot annotation, correctly handling
#' the Monte Carlo resolution floor and, optionally, a multiplicity
#' adjustment
#'
#' Motivation: a bare pd (e.g. "0.84") has no pre-existing social
#' calibration and reads as "fairly high" to most audiences even though
#' 2*(1-0.84)=0.32 is weak evidence by conventional standards - showing
#' both numbers together lets a reader cross-check one against the other
#' rather than anchoring on whichever scale their intuition misreads.
#'
#' FLOOR HANDLING: if the RAW (pre-adjustment) p-value for a row is at its
#' Monte Carlo resolution floor (2/ndraws_used - see .pd_to_pvalue()), the
#' true p-value could be smaller than what that many draws can resolve.
#' The displayed number is switched from "p = ..." to "p <= ..." in that
#' case. When p_value_adj is supplied, that adjusted number is what's
#' displayed either way - but the "<=" flag is still driven by the RAW
#' value's floor status, not reapplied to the floor itself. This is
#' deliberately a practical, slightly conservative convention for a
#' display label, not a full symbolic re-derivation of what a
#' multiplicity adjustment does to an interval-valued input (Holm's
#' cummax step means "<= raw floor" does not simply carry through to "<=
#' floor x multiplier" in every case) - it relies only on the fact that
#' Holm's adjustment is monotonic non-decreasing in each input, so
#' "<= [the adjusted value actually computed]" remains a valid, if not
#' maximally tight, bound whenever the raw input was floored.
#'
#' @param pd Probability of direction (vectorized).
#' @param p_value Raw (unadjusted) p-value column, e.g. from
#'   avg_comparisons_manual()'s output - used only to determine floor
#'   status, not necessarily what's displayed.
#' @param ndraws_used Number of draws backing each row - pass the column
#'   of the same name now returned by every avg_*_manual() comparison
#'   function, rather than a remembered/hardcoded number, so this stays
#'   correct even if a table combines rows from calls with different
#'   ndraws.
#' @param p_value_adj Optional adjusted p-value column (e.g.
#'   p.value.adj). If supplied, this is what gets displayed instead of
#'   p_value - pass NULL (default) to display the raw p_value.
#' @param digits Decimal places for the displayed p (default 3).
#' @return Character vector of labels, e.g. "pd = 1.00, p <= .004" or
#'   "pd = 0.84, p = .328".
format_pd_p_label <- function(pd, p_value, ndraws_used, p_value_adj = NULL, digits = 3) {
  floor_p <- 2 / ndraws_used
  is_floored <- p_value <= floor_p * (1 + 1e-8)
  p_display <- if (is.null(p_value_adj)) p_value else p_value_adj
  p_fmt <- formatC(p_display, digits = digits, format = "f")
  p_str <- ifelse(is_floored, paste0("p \u2264 ", p_fmt), paste0("p = ", p_fmt))
  sprintf("pd = %.2f, %s", pd, p_str)
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

#' Fail fast, for free, on a level typo/mismatch - before the (potentially
#' expensive) size check or any prediction call runs
#'
#' Catches exactly the class of error demonstrated in practice: a
#' hand-typed level (e.g. "Streetsidesquare") that doesn't exactly match
#' what the model's factor actually contains (e.g. "Streetside square",
#' with a space) only surfaces otherwise deep inside brms's own
#' validate_newdata(), AFTER the size-check prompt has already been shown
#' and confirmed and real computation has begun - wasting both the user's
#' time answering the prompt and whatever compute ran before brms's
#' internal check fired. Checking directly against newdata's own observed
#' levels first makes this instant and free instead.
.validate_levels <- function(newdata, colname, requested_levels) {
  observed <- unique(as.character(newdata[[colname]]))
  missing <- setdiff(requested_levels, observed)
  if (length(missing) > 0) {
    stop(
      "Requested level(s) not found in newdata$", colname, ": ",
      paste(sQuote(missing), collapse = ", "), ".\n",
      "Observed levels are: ", paste(sQuote(observed), collapse = ", "), ".\n",
      "Check for typos or spacing differences - factor levels must match exactly. ",
      "Consider leaving this argument NULL to pull levels directly from the data ",
      "instead of typing them by hand.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Coerce a vector of level values back to the type of the column they'll
#' be written into
#'
#' avg_comparisons_manual()/avg_comparisons_cross_manual()/
#' avg_simple_effects_manual() default their `levels` arguments via
#' sort(unique(as.character(newdata[[variable]]))) - always producing
#' character, even for a genuinely numeric column like UASEvents (which
#' feeds I(log10(UASEvents)) in the model formula). Writing a character
#' value into an existing numeric data.frame column (nd[[variable]] <- lv)
#' replaces the WHOLE column with a coerced character vector, not just
#' that one value - silently corrupting every row, not just the ones
#' being compared. Called immediately after `levels`/`levels1`/`levels2`/
#' `within_levels` are determined (whether defaulted or user-supplied) so
#' every downstream newdata write uses a correctly-typed value.
.coerce_levels_to_type <- function(levels_val, original_column) {
  levels_chr <- as.character(levels_val)
  if (is.numeric(original_column)) {
    as.numeric(levels_chr)
  } else if (is.logical(original_column)) {
    as.logical(levels_chr)
  } else if (is.factor(original_column)) {
    factor(levels_chr, levels = levels(original_column))
  } else {
    levels_chr
  }
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
  # Defensive coercion: data.table's `[.data.table` overrides base-R bracket
  # semantics (character/logical i is treated as a join/filter key, not
  # column selection), which silently breaks the row/column indexing below
  # if newdata was ever built or passed through as a data.table (e.g. by
  # some upstream marginaleffects::datagrid() pipelines). Coercing once
  # here means every function in this file behaves identically regardless
  # of what class of object the caller happened to hand it.
  newdata <- as.data.frame(newdata)
  n_full <- nrow(newdata)
  key <- .calibration_key(model, ...)
  
  # For n_conditions > 1: only ONE leg is being actively computed at a
  # time (lapply/parLapply is sequential per worker), so only one leg's
  # transient peak (which may be dominated by large family-specific
  # scratch space - see v0.6/v0.7 notes on the ordbeta case) applies at
  # once. Previously-finished legs only contribute their much smaller,
  # un-inflated RETAINED final matrix size, since their transient scratch
  # space has already been freed by the time the next leg runs. Naively
  # multiplying the whole transient-inclusive per-leg estimate by
  # n_conditions - as done through v0.13 - overcounts substantially for
  # multi-leg calls (confirmed in practice: 48.3 GB estimated for an
  # 8-leg avg_simple_effects_manual() call whose actual peak, per gc(),
  # was nowhere near that).
  dots <- list(...)
  ndraws_eff <- .effective_draws(model, ndraws = dots$ndraws, nsim = dots$nsim %||% 1000)
  naive_leg_mb <- (as.numeric(n_full) * ndraws_eff * 8) / 1e6  # one RETAINED matrix, un-inflated
  
  .total_from_leg_peak_mb <- function(leg_peak_mb) {
    (leg_peak_mb + (n_conditions - 1) * naive_leg_mb) / 1024
  }
  
  if (!force_recalibrate && exists(key, envir = .mfx_size_cache, inherits = FALSE)) {
    cached <- get(key, envir = .mfx_size_cache, inherits = FALSE)
    leg_peak_mb <- cached$fixed_overhead_mb + cached$per_row_mb * n_full
    return(.total_from_leg_peak_mb(leg_peak_mb))
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
    return(.total_from_leg_peak_mb(incremental_mb))
  }
  
  measured_mb <- vapply(sizes, .measure_at, numeric(1))
  fit <- stats::lm(measured_mb ~ sizes)
  fixed_overhead_mb <- max(0, unname(coef(fit)[1]))
  per_row_mb <- max(0, unname(coef(fit)[2]))
  
  assign(key, list(fixed_overhead_mb = fixed_overhead_mb, per_row_mb = per_row_mb),
         envir = .mfx_size_cache)
  
  leg_peak_mb <- fixed_overhead_mb + per_row_mb * n_full
  .total_from_leg_peak_mb(leg_peak_mb)
}

#' Pre-flight check: warn and require confirmation above a size threshold
#'
#' @param estimated_gb Output of .estimate_true_size_gb().
#' @param size_warn_gb Threshold in GB. If NULL, falls back to
#'   getOption("mfx_manual.size_warn_gb", default = 4).
#' @param context Short string describing which call this is, shown in the
#'   prompt/error (e.g. "avg_comparisons_manual(variable = 'Sex')").
#' Pre-flight check: warn and require confirmation above a size threshold
#'
#' @param estimated_gb Output of .estimate_true_size_gb().
#' @param size_warn_gb Threshold in GB. If NULL, falls back to
#'   getOption("mfx_manual.size_warn_gb", default = 4).
#' @param context Short string describing which call this is, shown in the
#'   prompt/error (e.g. "avg_comparisons_manual(variable = 'Sex')").
#' @param noninteractive_action What to do when interactive() is FALSE (an
#'   Rscript/knitr/scheduled batch run - utils::menu() has no user to
#'   answer it) and the threshold is exceeded. One of "error" (default -
#'   stop, since silently running an oversized computation unattended can
#'   still cause the system instability this check exists to prevent) or
#'   "warn" (emit an immediate, visible warning to the log and proceed
#'   automatically - appropriate for a scheduled job where you've already
#'   reviewed the sizes involved, or where recovering from a crash is
#'   cheaper than a blocked pipeline). If NULL (default), falls back to
#'   getOption("mfx_manual.noninteractive_action", default = "error").
.confirm_large_computation <- function(estimated_gb, size_warn_gb = NULL, context = "",
                                       noninteractive_action = NULL) {
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
    action <- if (is.null(noninteractive_action)) {
      getOption("mfx_manual.noninteractive_action", default = "error")
    } else {
      noninteractive_action
    }
    action <- match.arg(action, c("error", "warn"))
    
    if (action == "error") {
      stop(
        msg,
        "\nRunning in a non-interactive session - refusing to proceed automatically.",
        "\nRaise the threshold first if you are sure this is safe:",
        "\n  options(mfx_manual.size_warn_gb = <value>)  # or pass size_warn_gb= to this call",
        "\nOr allow non-interactive runs to proceed with just a warning:",
        "\n  options(mfx_manual.noninteractive_action = 'warn')  # or noninteractive_action='warn' to this call",
        call. = FALSE
      )
    } else {
      warning(
        msg,
        "\nRunning in a non-interactive session with ",
        "options(mfx_manual.noninteractive_action = 'warn') set - proceeding automatically.",
        call. = FALSE, immediate. = TRUE
      )
      return(invisible(TRUE))
    }
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
  # seed defaults to 1 deliberately, not arbitrarily: every multi-leg
  # function in this file (avg_comparisons_manual, avg_comparisons_cross_manual,
  # avg_slopes_manual) forwards the same ... to every leg via
  # .dispatch_draws_multi()'s lapply/parLapply, so all legs of one call
  # always share whichever seed is in effect - explicit or this default -
  # which is REQUIRED for simulated new-ID/StimID noise to cancel correctly
  # across differences/double-differences. There is no argument anywhere in
  # this file that lets different legs of the same comparison use different
  # seeds; if you ever extend this file, preserve that property.
  
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
#' @param noninteractive_action What to do if the size threshold is
#'   exceeded while running non-interactively (see .confirm_large_computation()).
#'   NULL (default) uses getOption("mfx_manual.noninteractive_action", "error").
#' @param ... Passed to .dispatch_draws() (re_formula, allow_new_levels,
#'   sample_new_levels, ndraws, nsim, seed).
avg_predictions_manual <- function(model, newdata, by = NULL,
                                   transform = identity, conf_level = 0.95,
                                   size_warn_gb = NULL,
                                   force_recalibrate = FALSE,
                                   noninteractive_action = NULL, ...) {
  newdata <- as.data.frame(newdata)  # see .estimate_true_size_gb() note on data.table
  on.exit({ gc(full = TRUE); gc(full = TRUE) }, add = TRUE)  # double pass: a single
  # gc() call was observed to leave Task
  # Manager showing elevated usage until a
  # second, manually-triggered call actually
  # dropped it - plausibly a Windows virtual-
  # memory-manager timing quirk rather than R
  # itself failing to mark the memory free
  est_gb <- .estimate_true_size_gb(model, newdata, n_conditions = 1,
                                   force_recalibrate = force_recalibrate, ...)
  .confirm_large_computation(est_gb, size_warn_gb, context = "avg_predictions_manual()",
                             noninteractive_action = noninteractive_action)
  
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
#'   IMPORTANT for the common case of exactly 2 levels: these two modes
#'   are MIRROR IMAGES of each other for identical `levels` input, not
#'   equivalent - "pairwise" always computes levels[1]-levels[2] (literal
#'   first-minus-second); "reference" always subtracts levels[1] (treating
#'   it as the baseline), giving levels[2]-levels[1]. With >2 levels these
#'   modes produce genuinely different SETS of comparisons and the
#'   distinction feels obvious; with exactly 2 they're forced to produce
#'   the same single comparison, just organized around two different
#'   principles (first-listed-is-the-minuend vs first-listed-is-the-
#'   baseline) that oppose in sign. Confirmed directly: levels=c("H520",
#'   "T150") with comparison_type="reference" gives "T150 - H520", but
#'   the SAME levels with comparison_type="pairwise" gives "H520 - T150".
#'   Rule of thumb for 2 levels: "reference" always subtracts levels[1];
#'   "pairwise" always subtracts levels[2].
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
#' @param auto_orient If TRUE, each row's sign and 'contrast' label are
#'   chosen from the ACTUAL estimated difference (positive result reported
#'   as-is; negative result has its sign negated and its two level names
#'   swapped in the label) rather than from levels/comparison_type's input
#'   order. Default FALSE preserves the literal, order-determined behaviour
#'   documented above and used by every existing call in this project - so
#'   turning this on cannot silently change results you have already
#'   validated with the default. Useful when you want every row of a
#'   multi-pair table (e.g. all pairwise AmbientEnvCore comparisons) to
#'   read as a positive number for readability, which is NOT achievable in
#'   general by choosing a clever levels= order: with 3+ levels,
#'   comparison_type = "pairwise"'s combn(levels, 2) structurally forces
#'   the first input element to lead 2 of 3 pairs and the second to lead
#'   the remaining 1, so no reordering can give two DIFFERENT elements
#'   each a turn leading multiple pairs (e.g. wanting both "3-1" and "3-2"
#'   positive simultaneously is not reachable this way - confirmed in
#'   practice, not just in principle).
#' @param noninteractive_action What to do if the size threshold is
#'   exceeded while running non-interactively (see .confirm_large_computation()).
#'   NULL (default) uses getOption("mfx_manual.noninteractive_action", "error").
#' @param ... Passed to .dispatch_draws() (re_formula, allow_new_levels,
#'   sample_new_levels, ndraws, nsim, seed).
avg_comparisons_manual <- function(model, newdata, variable,
                                   comparison_type = c("pairwise", "reference", "sequential"),
                                   levels = NULL, transform = identity,
                                   conf_level = 0.95, p_adjust = "none",
                                   rope = NULL, cl = NULL,
                                   size_warn_gb = NULL,
                                   force_recalibrate = FALSE,
                                   auto_orient = FALSE,
                                   noninteractive_action = NULL, ...) {
  newdata <- as.data.frame(newdata)  # see .estimate_true_size_gb() note on data.table
  on.exit({ gc(full = TRUE); gc(full = TRUE) }, add = TRUE)  # double pass, see avg_predictions_manual() note
  comparison_type <- match.arg(comparison_type)
  if (is.null(levels)) levels <- sort(unique(as.character(newdata[[variable]])))
  .validate_levels(newdata, variable, levels)
  
  est_gb <- .estimate_true_size_gb(model, newdata, n_conditions = length(levels),
                                   force_recalibrate = force_recalibrate, ...)
  .confirm_large_computation(
    est_gb, size_warn_gb,
    context = sprintf("avg_comparisons_manual(variable = '%s', %d levels)", variable, length(levels)),
    noninteractive_action = noninteractive_action
  )
  
  pairs <- switch(comparison_type,
                  pairwise   = t(utils::combn(levels, 2)),
                  reference  = cbind(levels[-1], levels[1]),
                  sequential = cbind(levels[-1], levels[-length(levels)])
  )
  
  # Unconditional confirmation of exactly what's being computed - not just
  # when levels default (v0.23). "pairwise" computes levels[1]-levels[2]
  # LITERALLY regardless of which level the model treats as its reference;
  # "reference" computes (nonreference - levels[1], treating levels[1] AS
  # the reference). These give OPPOSITE required level orderings to match
  # a regression coefficient's sign convention (reference listed first
  # under "reference" mode, but listed SECOND under "pairwise" mode) -
  # confirmed as a recurring, genuine point of confusion, not something to
  # leave implicit in the output's contrast label alone.
  if (auto_orient) {
    message(sprintf(
      "avg_comparisons_manual(comparison_type = '%s', auto_orient = TRUE): initial pairing %s - final sign/label of each row is DATA-DRIVEN (flipped to positive where the estimated difference came out negative), not fixed by this initial order. Check each row's own 'contrast' label, not this message, for its actual reported direction.",
      comparison_type,
      paste(sprintf("'%s' - '%s'", pairs[, 1], pairs[, 2]), collapse = ", ")
    ))
  } else {
    message(sprintf(
      "avg_comparisons_manual(comparison_type = '%s'): computing %s",
      comparison_type,
      paste(sprintf("'%s' - '%s'", pairs[, 1], pairs[, 2]), collapse = ", ")
    ))
  }
  
  # Predict once per level (not once per pair) and reuse - avoids redundant
  # prediction of shared levels across multiple pairwise comparisons.
  # Coercion to the variable's real type happens HERE ONLY, per scalar
  # value, not on the whole `levels` vector (see v0.28) - keeps `levels`
  # plain character everywhere else (pairs, labels, names), since
  # combn()/cbind()/lapply() do not reliably preserve factor labels
  # through generic list/matrix construction.
  newdata_list <- lapply(levels, function(lv) {
    nd <- newdata
    nd[[variable]] <- .coerce_levels_to_type(lv, newdata[[variable]])
    nd
  })
  names(newdata_list) <- levels
  draws_list <- .dispatch_draws_multi(model, newdata_list, cl = cl, ...)
  mean_draws_by_level <- lapply(draws_list, function(d) rowMeans(transform(d)))
  names(mean_draws_by_level) <- levels
  
  alpha <- 1 - conf_level
  results <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
    hi <- pairs[i, 1]; lo <- pairs[i, 2]
    d <- mean_draws_by_level[[hi]] - mean_draws_by_level[[lo]]
    if (auto_orient && mean(d) < 0) {
      # Sign is DATA-DRIVEN, not guessed from input order - reordering
      # `levels` alone cannot achieve arbitrary independent signs across
      # 3+ pairwise comparisons (combn() structurally forces the first
      # input element to lead 2 of 3 pairs and the second to lead the
      # remaining 1 - no permutation makes two different elements each
      # lead multiple pairs). Flip AFTER seeing the real result instead:
      # negate the draws and swap the label together, so `contrast` and
      # `estimate` always agree, and every other quantity below (pd,
      # p.value, ROPE) is computed from the already-oriented `d`.
      d <- -d
      tmp <- hi; hi <- lo; lo <- tmp
    }
    pd <- as.numeric(bayestestR::p_direction(d))
    row <- data.frame(
      contrast = paste(hi, "-", lo),
      estimate = mean(d),
      conf.low = unname(stats::quantile(d, alpha / 2)),
      conf.high = unname(stats::quantile(d, 1 - alpha / 2)),
      pd = pd,
      p.value = .pd_to_pvalue(pd, length(d)),  # floored at this call's own MC resolution - see .pd_to_pvalue()
      ndraws_used = length(d)  # lets format_pd_p_label() recompute the exact floor per row, no need to pass ndraws separately
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
# avg_comparisons_cross_manual()
# -----------------------------------------------------------------------------

#' Difference-in-differences (interaction) contrast across two variables
#'
#' Computes, for one or more levels of `variable1`, how the `variable2`
#' contrast (hi2 - lo2) differs from that same contrast at a reference (or
#' every other) level of `variable1` - directly on the response/annoyance
#' scale, with a proper credible interval. This is the natural way to
#' formally test and report an interaction term's effect size (e.g. does
#' the UASType gap differ between AmbientEnvCore levels) in the same
#' interpretable units as avg_comparisons_manual(), rather than reading an
#' interaction coefficient off the logit-scale summary table.
#'
#' WHY variable1 CAN HAVE MANY LEVELS BUT variable2 CANNOT: a k-level
#' factor's interaction with another factor is parameterized by the model
#' as k-1 coefficients relative to a reference level (exactly what your own
#' brms coefficient table shows, e.g. AmbientEnvCorePark:UASTypeT150,
#' AmbientEnvCoreResidential:UASTypeT150, AmbientEnvCoreStreetsidesquare:UASTypeT150
#' - three terms for a 4-level AmbientEnvCore, not six pairwise
#' combinations). comparison_type1 = "reference" (the default) reproduces
#' exactly that parameterization: it holds variable2's contrast fixed at
#' the two levels you specify, and asks how that contrast differs at each
#' non-reference level of variable1 relative to variable1's reference
#' level - one row per surviving model coefficient, directly comparable to
#' it. Extending variable2 to more than 2 levels as well would reintroduce
#' genuine ambiguity (now BOTH sides would have multiple possible
#' reference/pairwise choices simultaneously), so that is intentionally
#' not supported - run this function once per variable2 contrast you want
#' to report, each time examining it across every level of variable1.
#'
#' CORRECTNESS NOTE: every (variable1 level x variable2 level) combination
#' MUST use the same seed (passed via ... exactly as for
#' avg_comparisons_manual() - and shares .dispatch_draws()'s seed=1
#' default automatically if you don't specify one, see that function's
#' comment) so that the shared simulated new-ID/StimID noise from
#' sample_new_levels="gaussian" cancels identically across every leg
#' before the differences are taken. This matters more here than for a
#' simple pairwise contrast: with only partial cancellation, uncancelled
#' new-level noise doesn't just add a little extra width, it can dominate
#' a difference-of-differences whose true magnitude is often small
#' relative to either single difference.
#'
#' @param model brmsfit or glmgee object.
#' @param newdata Grid of predictor values, at variable1/variable2's
#'   current (baseline) values - overwritten per combination internally.
#' @param variable1 Name of the focal categorical column examined across
#'   (potentially) many levels.
#' @param variable2 Name of the focal categorical column whose 2-level
#'   contrast is being tracked across variable1's levels.
#' @param levels1 Optional character vector of variable1 levels to include.
#'   Defaults to ALL unique observed levels. The first element is treated
#'   as the reference level under comparison_type1 = "reference" - order
#'   matters, so pass this explicitly if the alphabetical default isn't
#'   the reference level your model actually used.
#' @param levels2 Optional length-2 character vector selecting which two
#'   levels of variable2 define the contrast. Defaults to the first two
#'   unique observed levels - specify explicitly if variable2 has more
#'   than two levels or if the default ordering isn't what you want.
#' @param comparison_type1 "reference" (default - each non-reference
#'   variable1 level vs. levels1[1], matching the model's own dummy
#'   coding) or "pairwise" (every pair of variable1 levels - only sensible
#'   for a small number of levels, since it grows combinatorially).
#' @param transform Function applied to each combination's raw draws matrix.
#' @param conf_level Credible/confidence interval coverage (default 0.95).
#' @param rope Optional numeric length-2 vector for a Kruschke (2018)-style
#'   practical-equivalence judgement on each interaction contrast (see
#'   avg_comparisons_manual() for the same argument's full rationale).
#' @param p_adjust Method passed to stats::p.adjust() across the returned
#'   rows (see avg_comparisons_manual() for the same caveat on treating
#'   this as a hybrid rather than a purely Bayesian procedure).
#' @param cl Optional cluster from make_prediction_cluster().
#' @param size_warn_gb,force_recalibrate,noninteractive_action See
#'   avg_comparisons_manual() - identical semantics, sized for
#'   length(levels1) * 2 combinations.
#' @param ... Passed to .dispatch_draws() (re_formula, allow_new_levels,
#'   sample_new_levels, ndraws, nsim, seed).
avg_comparisons_cross_manual <- function(model, newdata, variable1, variable2,
                                         levels1 = NULL, levels2 = NULL,
                                         comparison_type1 = c("reference", "pairwise"),
                                         transform = identity, conf_level = 0.95,
                                         rope = NULL, p_adjust = "none", cl = NULL,
                                         size_warn_gb = NULL,
                                         force_recalibrate = FALSE,
                                         noninteractive_action = NULL, ...) {
  comparison_type1 <- match.arg(comparison_type1)
  newdata <- as.data.frame(newdata)
  on.exit({ gc(full = TRUE); gc(full = TRUE) }, add = TRUE)  # double pass, see avg_predictions_manual() note
  if (is.null(levels1)) levels1 <- sort(unique(as.character(newdata[[variable1]])))
  levels2_explicit <- !is.null(levels2)
  if (!levels2_explicit) {
    # Same fix and rationale as avg_simple_effects_manual(): alphabetically-
    # first level becomes the reference (subtracted) by default, matching
    # avg_comparisons_manual()'s "reference" convention. Only applies when
    # levels2 is NOT explicitly supplied - explicit levels2 (e.g.
    # c("T150", "H520")) are used in the exact order given, unaffected.
    sorted_levels2 <- sort(unique(as.character(newdata[[variable2]])))[1:2]
    levels2 <- rev(sorted_levels2)
  }
  .validate_levels(newdata, variable1, levels1)
  .validate_levels(newdata, variable2, levels2)
  if (length(levels2) != 2) {
    stop(
      "levels2 must specify exactly 2 levels of '", variable2, "' - the contrast being ",
      "tracked across levels of '", variable1, "'. Pass levels2 explicitly if '", variable2,
      "' has more than 2 observed levels."
    )
  }
  if (!levels2_explicit) {
    message(sprintf("avg_comparisons_cross_manual(): levels2 not supplied - tracking '%s' - '%s' across levels of '%s'.",
                    levels2[1], levels2[2], variable1))
  }
  
  n_conditions <- length(levels1) * 2
  est_gb <- .estimate_true_size_gb(model, newdata, n_conditions = n_conditions,
                                   force_recalibrate = force_recalibrate, ...)
  .confirm_large_computation(
    est_gb, size_warn_gb,
    context = sprintf("avg_comparisons_cross_manual('%s' [%d levels] x '%s')",
                      variable1, length(levels1), variable2),
    noninteractive_action = noninteractive_action
  )
  
  combos <- expand.grid(l1 = levels1, l2 = levels2, stringsAsFactors = FALSE)
  key <- function(l1, l2) paste(l1, l2, sep = "___")
  
  # Coercion happens HERE ONLY, per scalar value - see v0.28 rationale in
  # avg_comparisons_manual().
  newdata_list <- lapply(seq_len(nrow(combos)), function(i) {
    nd <- newdata
    nd[[variable1]] <- .coerce_levels_to_type(combos$l1[i], newdata[[variable1]])
    nd[[variable2]] <- .coerce_levels_to_type(combos$l2[i], newdata[[variable2]])
    nd
  })
  names(newdata_list) <- key(combos$l1, combos$l2)
  
  # Every leg shares one seed via ... (explicit, or .dispatch_draws()'s
  # own seed=1 default) - required for correct cancellation of shared
  # simulated new-level noise across all differences (see docstring note).
  draws_list <- .dispatch_draws_multi(model, newdata_list, cl = cl, ...)
  mean_draws <- lapply(draws_list, function(d) rowMeans(transform(d)))
  
  hi2 <- levels2[1]; lo2 <- levels2[2]
  diff_by_level1 <- lapply(levels1, function(l1) mean_draws[[key(l1, hi2)]] - mean_draws[[key(l1, lo2)]])
  names(diff_by_level1) <- levels1
  
  pairs1 <- switch(comparison_type1,
                   reference = if (length(levels1) < 2) {
                     stop("comparison_type1 = 'reference' needs at least 2 levels1.")
                   } else {
                     cbind(levels1[-1], levels1[1])
                   },
                   pairwise = t(utils::combn(levels1, 2))
  )
  
  # Same rationale as avg_comparisons_manual()'s confirmation message: state
  # exactly which (variable1 level A vs B, tracking the variable2 hi-lo gap)
  # comparisons are about to be computed, unconditionally.
  message(sprintf(
    "avg_comparisons_cross_manual(comparison_type1 = '%s'): tracking '%s' - '%s' across %s",
    comparison_type1, hi2, lo2,
    paste(sprintf("'%s' vs '%s'", pairs1[, 1], pairs1[, 2]), collapse = ", ")
  ))
  
  alpha <- 1 - conf_level
  results <- do.call(rbind, lapply(seq_len(nrow(pairs1)), function(i) {
    a <- pairs1[i, 1]; b <- pairs1[i, 2]
    dd <- diff_by_level1[[a]] - diff_by_level1[[b]]
    pd <- as.numeric(bayestestR::p_direction(dd))
    row <- data.frame(
      level1 = a, reference1 = b,  # plain columns for plotting - no string parsing needed
      contrast = sprintf("(%s=%s: %s-%s) vs (%s=%s: %s-%s)",
                         variable1, a, hi2, lo2, variable1, b, hi2, lo2),
      estimate = mean(dd),
      conf.low = unname(quantile(dd, alpha / 2)),
      conf.high = unname(quantile(dd, 1 - alpha / 2)),
      pd = pd, p.value = .pd_to_pvalue(pd, length(dd)), ndraws_used = length(dd)
    )
    if (!is.null(rope)) {
      row$rope_decision <- if (row$conf.low > rope[2] || row$conf.high < rope[1]) {
        "outside ROPE (credible non-negligible interaction)"
      } else if (row$conf.low > rope[1] && row$conf.high < rope[2]) {
        "inside ROPE (practically equivalent interaction sizes)"
      } else {
        "overlaps ROPE (undecided)"
      }
    }
    row
  }))
  names(results)[1] <- variable1  # e.g. "AmbientEnvCore" rather than the generic "level1"
  results$p.value.adj <- stats::p.adjust(results$p.value, method = p_adjust)
  results
}


# -----------------------------------------------------------------------------
# avg_simple_effects_manual()
# -----------------------------------------------------------------------------

#' Simple effects: a 2-level contrast, computed independently within each
#' level of another factor
#'
#' Distinct from avg_comparisons_cross_manual()'s INTERACTION contrasts,
#' which are the DIFFERENCES between simple effects relative to a reference
#' level. This function returns the simple effects themselves - e.g. the
#' T150-H520 gap computed separately within each of Highway/Park/
#' Residential/Streetsidesquare, four independent numbers, none expressed
#' "relative to" anything else. This is the natural quantity for a plot
#' showing how one factor's effect varies across another factor's levels
#' (e.g. "H520-T150 across all four environments"); avg_comparisons_cross_manual()
#' is the natural quantity for FORMALLY TESTING whether that variation is
#' credibly non-zero relative to a chosen reference. The two are related by
#' simple subtraction: interaction_contrast(level) = simple_effect(level) -
#' simple_effect(reference) - so running both on the same data should give
#' internally consistent numbers, worth spot-checking if you use both.
#'
#' @param model brmsfit or glmgee object.
#' @param newdata Grid of predictor values, at variable/within_variable's
#'   current (baseline) values - overwritten per combination internally.
#' @param variable Name of the focal categorical column whose 2-level
#'   contrast is computed.
#' @param within_variable Name of the categorical column whose levels the
#'   contrast is computed separately within.
#' @param levels Optional length-2 character vector selecting which two
#'   levels of `variable` define the contrast. Defaults to the first two
#'   unique observed levels - specify explicitly if `variable` has more
#'   than two levels or the default ordering isn't what you want.
#' @param within_levels Optional character vector of within_variable
#'   levels to compute the simple effect within. Defaults to ALL unique
#'   observed levels.
#' @param transform Function applied to each combination's raw draws matrix.
#' @param conf_level Credible/confidence interval coverage (default 0.95).
#' @param rope Optional numeric length-2 vector for a Kruschke (2018)-style
#'   practical-equivalence judgement on each simple effect (see
#'   avg_comparisons_manual() for the same argument's full rationale).
#' @param p_adjust Method passed to stats::p.adjust() across the returned
#'   rows (see avg_comparisons_manual() for the same caveat on treating
#'   this as a hybrid rather than a purely Bayesian procedure).
#' @param cl Optional cluster from make_prediction_cluster().
#' @param size_warn_gb,force_recalibrate,noninteractive_action See
#'   avg_comparisons_manual() - identical semantics, sized for
#'   length(within_levels) * 2 combinations.
#' @param ... Passed to .dispatch_draws() (re_formula, allow_new_levels,
#'   sample_new_levels, ndraws, nsim, seed - all legs automatically share
#'   one seed, see .dispatch_draws()'s comment; correct for this function
#'   too even though there's no double-difference to protect here, since
#'   it keeps every leg's new-level draws comparable to each other).
avg_simple_effects_manual <- function(model, newdata, variable, within_variable,
                                      levels = NULL, within_levels = NULL,
                                      transform = identity, conf_level = 0.95,
                                      rope = NULL, p_adjust = "none", cl = NULL,
                                      size_warn_gb = NULL,
                                      force_recalibrate = FALSE,
                                      noninteractive_action = NULL, ...) {
  newdata <- as.data.frame(newdata)
  on.exit({ gc(full = TRUE); gc(full = TRUE) }, add = TRUE)  # double pass, see avg_predictions_manual() note
  levels_explicit <- !is.null(levels)
  if (!levels_explicit) {
    # Match avg_comparisons_manual()'s "reference" convention: the
    # alphabetically-first level is the REFERENCE (subtracted), everything
    # else is reported as itself minus the reference - standard dummy-
    # coding direction, matching how the model's own coefficients work
    # (e.g. AmbientEnvCorePark:UASTypeT150 means Park relative to Highway).
    # Without this, the two functions' DEFAULT direction disagreed (this
    # one previously did levels[1] - levels[2] literally, giving the
    # opposite sign to avg_comparisons_manual() for the same data under
    # the same default sort() ordering). Only applies when levels is NOT
    # explicitly supplied - explicit levels are still used in the exact
    # order given, so existing calls with explicit levels (e.g.
    # levels = c("T150", "H520")) are completely unaffected by this change.
    sorted_levels <- sort(unique(as.character(newdata[[variable]])))[1:2]
    levels <- rev(sorted_levels)
  }
  if (is.null(within_levels)) within_levels <- sort(unique(as.character(newdata[[within_variable]])))
  .validate_levels(newdata, variable, levels)
  .validate_levels(newdata, within_variable, within_levels)
  if (length(levels) != 2) {
    stop(
      "levels must specify exactly 2 levels of '", variable, "' - the contrast being ",
      "computed within each level of '", within_variable, "'. Pass levels explicitly if '",
      variable, "' has more than 2 observed levels."
    )
  }
  message(sprintf(
    "avg_simple_effects_manual(): computing '%s' - '%s' (per level of '%s')%s.",
    levels[1], levels[2], within_variable,
    if (levels_explicit) "" else " [levels defaulted]"
  ))
  
  n_conditions <- length(within_levels) * 2
  est_gb <- .estimate_true_size_gb(model, newdata, n_conditions = n_conditions,
                                   force_recalibrate = force_recalibrate, ...)
  .confirm_large_computation(
    est_gb, size_warn_gb,
    context = sprintf("avg_simple_effects_manual('%s' within '%s' [%d levels])",
                      variable, within_variable, length(within_levels)),
    noninteractive_action = noninteractive_action
  )
  
  combos <- expand.grid(w = within_levels, l = levels, stringsAsFactors = FALSE)
  key <- function(w, l) paste(w, l, sep = "___")
  
  # Coercion happens HERE ONLY, per scalar value - see v0.28 rationale in
  # avg_comparisons_manual().
  newdata_list <- lapply(seq_len(nrow(combos)), function(i) {
    nd <- newdata
    nd[[within_variable]] <- .coerce_levels_to_type(combos$w[i], newdata[[within_variable]])
    nd[[variable]] <- .coerce_levels_to_type(combos$l[i], newdata[[variable]])
    nd
  })
  names(newdata_list) <- key(combos$w, combos$l)
  
  draws_list <- .dispatch_draws_multi(model, newdata_list, cl = cl, ...)
  mean_draws <- lapply(draws_list, function(d) rowMeans(transform(d)))
  
  hi <- levels[1]; lo <- levels[2]
  alpha <- 1 - conf_level
  results <- do.call(rbind, lapply(within_levels, function(w) {
    d <- mean_draws[[key(w, hi)]] - mean_draws[[key(w, lo)]]
    pd <- as.numeric(bayestestR::p_direction(d))
    row <- data.frame(
      within_level = w,
      contrast = paste(hi, "-", lo),
      estimate = mean(d),
      conf.low = unname(quantile(d, alpha / 2)),
      conf.high = unname(quantile(d, 1 - alpha / 2)),
      pd = pd, p.value = .pd_to_pvalue(pd, length(d)), ndraws_used = length(d)
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
  names(results)[1] <- within_variable
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
#' @param force_recalibrate If TRUE, ignores any cached memory calibration
#'   for this model/settings combination and re-measures from scratch.
#' @param noninteractive_action What to do if the size threshold is
#'   exceeded while running non-interactively (see .confirm_large_computation()).
#'   NULL (default) uses getOption("mfx_manual.noninteractive_action", "error").
#' @param ... Passed to .dispatch_draws().
avg_slopes_manual <- function(model, newdata, variable, eps = NULL,
                              transform = identity, conf_level = 0.95,
                              cl = NULL, size_warn_gb = NULL,
                              force_recalibrate = FALSE,
                              noninteractive_action = NULL, ...) {
  newdata <- as.data.frame(newdata)  # see .estimate_true_size_gb() note on data.table
  on.exit({ gc(full = TRUE); gc(full = TRUE) }, add = TRUE)  # double pass, see avg_predictions_manual() note
  x <- newdata[[variable]]
  if (is.null(eps)) eps <- 0.0001 * diff(range(x))
  
  est_gb <- .estimate_true_size_gb(model, newdata, n_conditions = 2,
                                   force_recalibrate = force_recalibrate, ...)
  .confirm_large_computation(
    est_gb, size_warn_gb,
    context = sprintf("avg_slopes_manual(variable = '%s')", variable),
    noninteractive_action = noninteractive_action
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
    p.value = .pd_to_pvalue(pd, length(slope_draws)),
    ndraws_used = length(slope_draws)
  )
}


# -----------------------------------------------------------------------------
# add_pairwise_brackets()
# -----------------------------------------------------------------------------

#' Annotate a means/estimates plot with pairwise-contrast credible intervals
#' as brackets, using each contrast's own CI rather than pd/p (see
#' conversation notes on why: pd and Holm-adjusted p both showed real
#' artifacts here - miscalibrated raw pd, and floored p-values distorting
#' Holm's ranking - while a contrast's own credible interval has neither
#' problem and answers "distinguishable from zero" directly).
#'
#' Accepts output from EITHER avg_comparisons_cross_manual() (which has
#' `x_var`/"reference1" columns already) or avg_comparisons_manual() (which
#' only has a formatted "contrast" string, e.g. "Park - Highway", parsed
#' apart here via strsplit on " - " - safe as long as no level name
#' contains that exact substring itself).
#'
#' Deliberately built on plain ggpubr::stat_pvalue_manual() rather than a
#' hand-rolled geom_segment()/geom_text() version: stat_pvalue_manual()
#' works fine here specifically BECAUSE the base plot is NOT coord_flip()'d
#' (coord.flip = FALSE is hardcoded below) - it has a confirmed, hard-to-
#' work-around interaction problem with coord_flip() specifically (matching
#' its own coord.flip argument to the base plot's actual flipped state
#' gives correctly-positioned brackets but illegible sideways text;
#' mismatching it gives legible text but misaligned bracket positions).
#' If you need this on a coord_flip()'d plot, use the manual geom_segment/
#' geom_text approach from the conversation instead of this function.
#'
#' Also requires stat_pvalue_manual() to see columns literally named
#' group1/group2 (an internal assertion, not something an xmin/xmax
#' argument overrides) - handled internally, no action needed from the
#' caller.
#'
#' Bracket rows are ordered by SPAN (how far apart the two x-positions
#' are, via x_order), narrowest first - stacks nested/overlapping brackets
#' outward in the standard convention and guarantees pairs sharing a
#' starting category (e.g. two comparisons both involving "Park") don't
#' land on the same visual row and overlap.
#'
#' @param base_plot A ggplot object (built from `means_df`, x-axis =
#'   `x_var`, NOT coord_flip()'d) to add brackets to.
#' @param means_df The data.frame the base plot's points/error bars were
#'   built from (e.g. output of avg_predictions_manual() or
#'   avg_simple_effects_manual()) - used only for its own conf.high range,
#'   to position brackets above the actual plotted points. Do NOT pass
#'   contrast_df here by mistake - they are on different scales (this
#'   exact mix-up produced brackets that cut through the data points in
#'   practice).
#' @param contrast_df The pairwise-contrast results to annotate with -
#'   output of avg_comparisons_manual(), avg_comparisons_cross_manual(),
#'   or avg_simple_effects_manual() applied across x_var's levels.
#' @param x_var Name of the column on the base plot's x-axis (e.g.
#'   "AmbientEnvCore" or "UASEvents") - works for either a discrete
#'   (character/factor) or continuous (numeric) x-axis; handled
#'   automatically based on the type of `means_df[[x_var]]` (see
#'   x_order below).
#' @param x_order Only used when `means_df[[x_var]]` is NOT numeric.
#'   Character vector giving the x-axis category order as displayed (e.g.
#'   env_order) - used to compute each bracket's span for the narrowest-
#'   first stacking. Optional: defaults to `levels(means_df[[x_var]])` if
#'   it's already a factor, otherwise `sort(unique(as.character(...)))`.
#'   Ignored (with a message) if supplied alongside a numeric x_var, since
#'   a numeric value already is its own x-position - no lookup needed.
#' @param label_size,step_increase,y_margin Passed through to control
#'   label text size, vertical spacing between stacked brackets, and
#'   headroom above the highest data point before the first bracket.
#' @return `base_plot` with a stat_pvalue_manual() bracket layer added.
add_pairwise_brackets <- function(base_plot, means_df, contrast_df, x_var, x_order = NULL,
                                  label_size = 4, step_increase = 0.14, y_margin = 0.3) {
  if (all(c(x_var, "reference1") %in% names(contrast_df))) {
    contrast_df$group1 <- contrast_df[[x_var]]
    contrast_df$group2 <- contrast_df$reference1
  } else if ("contrast" %in% names(contrast_df)) {
    parts <- strsplit(contrast_df$contrast, " - ", fixed = TRUE)
    contrast_df$group1 <- sapply(parts, `[`, 1)
    contrast_df$group2 <- sapply(parts, `[`, 2)
  } else {
    stop("contrast_df needs either ('", x_var, "', 'reference1') or a 'contrast' column to parse.")
  }
  
  contrast_df$ci_label <- sprintf("%.2f [%.2f, %.2f]", contrast_df$estimate, contrast_df$conf.low, contrast_df$conf.high)
  
  x_col <- means_df[[x_var]]
  if (is.numeric(x_col)) {
    # Continuous x-axis: group1/group2 are currently character (from
    # string-splitting or a factor/character column) - coerce directly to
    # the SAME numeric scale the plot already uses, rather than routing
    # through an artificial integer position index (which is what caused
    # the "discrete value supplied to continuous scale" error: character
    # xmin/xmax values on a plot whose x-scale is actually continuous).
    if (!is.null(x_order)) {
      message("add_pairwise_brackets(): x_order is ignored for a numeric x_var - ",
              "the values themselves already define their position.")
    }
    contrast_df$group1 <- as.numeric(contrast_df$group1)
    contrast_df$group2 <- as.numeric(contrast_df$group2)
    span <- abs(contrast_df$group1 - contrast_df$group2)
  } else {
    if (is.null(x_order)) {
      x_order <- if (is.factor(x_col)) levels(x_col) else sort(unique(as.character(x_col)))
    }
    x_pos <- setNames(seq_along(x_order), x_order)
    span <- abs(x_pos[contrast_df$group1] - x_pos[contrast_df$group2])
  }
  
  contrast_df <- contrast_df[order(span), ]
  
  base_plot +
    ggpubr::stat_pvalue_manual(
      contrast_df, label = "ci_label",
      y.position = max(means_df$conf.high) + y_margin,  # always the MEANS data's own range, hard-coded once here
      step.increase = step_increase, coord.flip = FALSE,
      family = "serif", size = label_size
    )
}


# =============================================================================
# USAGE EXAMPLES
# (not run automatically - copy into your script and adapt object names)
# =============================================================================
if (FALSE) {
  
  source("avg_comparisons_slopes.R")
  
  # --- Example 0: the size check runs automatically inside every call below
  #     - there's nothing extra to run first. The FIRST call against a given
  #     model/settings combination in your session pays a small dry-run
  #     cost (proportionally sized, never larger than a few percent of your
  #     grid); every SUBSEQUENT call with the same model/settings (e.g.
  #     testing several different factors one after another - your typical
  #     usage pattern) reuses the cached calibration and costs essentially
  #     nothing extra. Call reset_size_cache() if you ever want to force
  #     recalibration (e.g. after updating brms or refitting the model).
  
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
  
  # If you've refit the model or updated brms and don't trust the cached
  # calibration anymore, clear it (affects all models/settings) or force
  # a single call to recalibrate itself:
  # reset_size_cache()
  # avg_comparisons_manual(..., force_recalibrate = TRUE)
  
  # --- Example 0c: running as part of a scheduled/Rscript pipeline, where
  #     you've already reviewed the sizes involved and want oversized
  #     calls to proceed with a loud log warning rather than block the job
  options(mfx_manual.noninteractive_action = "warn")
  # or per-call, without changing the session-wide default:
  # avg_predictions_manual(..., noninteractive_action = "warn")
  # Default ("error") is unchanged for any call that doesn't set this -
  # an unattended run silently proceeding past a size threshold is exactly
  # the scenario this whole check exists to prevent, so opt in deliberately.
  
  # --- Example -1a: build the fully marginal grid (new ID and StimID) used
  #     for population-level estimates (LAE dose-response, demographic and
  #     interaction contrasts) - replaces the hand-rolled base_unique/
  #     design_cells/demog_profiles/crossing pattern from earlier in this
  #     project with one validated call. No seed needed - profile
  #     assignment is deterministic (see v0.18); the only randomness left
  #     is sample_new_levels="gaussian" itself, controlled downstream via
  #     the seed= argument to avg_predictions_manual() etc.
  grid_marg <- build_marginal_grid(
    m0d,
    within_vars = c("UASProximity", "UASType", "UASEvents", "AmbientEnvCore",
                    "AmbientVariant", "UASOperation"),
    between_vars = c("NationGeo2", "AAM_attitude2", "Home_Area2", "Age", "Sex",
                     "NoiseSensitivity", "Area_soundscape2"),
    simulate_id = TRUE, simulate_stim = TRUE, n_new = nnew_pred,
    counterfactual_vars = list(UASLAEMaxLRScl = seq(-11, 12, length.out = 11)),
    fixed_at_mean = "TrialNumberScl"
  )
  
  # --- Example -1b: build the CONDITIONAL grid (real, existing StimID) used
  #     for the calibration-against-observed-stimuli plot - same function,
  #     simulate_stim = FALSE is the only difference, making the estimand
  #     choice explicit in the call rather than implicit in a different,
  #     hand-copied block of code.
  grid_calib <- build_marginal_grid(
    m0d,
    within_vars = c("UASProximity", "UASType", "UASEvents", "AmbientEnvCore",
                    "AmbientVariant", "UASOperation"),
    between_vars = c("NationGeo2", "AAM_attitude2", "Home_Area2", "Age", "Sex",
                     "NoiseSensitivity", "Area_soundscape2"),
    simulate_id = TRUE, simulate_stim = FALSE,  # <- real StimID kept, for matching against dAnnoyanceMean
    n_new = nnew_pred,
    fixed_at_mean = "TrialNumberScl"
  )
  
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
  
  # --- Example 3b: interaction contrast - does the UASType gap (T150 vs
  #     H520) differ across ALL levels of AmbientEnvCore, relative to a
  #     reference level, in annoyance-scale units? Returns one row per
  #     non-reference AmbientEnvCore level - directly comparable one-to-one
  #     against the model's own AmbientEnvCore<level>:UASTypeT150
  #     coefficients, just expressed as an interpretable effect size with
  #     a credible interval instead of a logit-scale coefficient.
  type_by_env <- avg_comparisons_cross_manual(
    m0d, newdata = demo_grid, variable1 = "AmbientEnvCore", variable2 = "UASType",
    levels1 = c("Highway", "Park", "Residential", "Streetsidesquare"),  # Highway = reference,
    comparison_type1 = "reference",                                     # matching the model's own coding
    levels2 = c("T150", "H520"),
    re_formula = NULL, allow_new_levels = TRUE, sample_new_levels = "gaussian",
    ndraws = ndraws_pred, seed = 999,
    rope = c(-0.2, 0.2),
    transform = \(x) betaTransform(x, -10, 10, direction = "reverse", squeeze = "none")
  )
  # type_by_env now has 3 rows: Park vs Highway, Residential vs Highway,
  # Streetsidesquare vs Highway - each answering "how much bigger/smaller
  # is the T150-H520 gap here than at the Highway reference"
  
  # --- Example 3c: simple effects - the actual H520-T150 gap, computed
  #     independently within EACH of the four AmbientEnvCore levels. This
  #     - not the interaction contrasts above - is the quantity for a plot
  #     of "H520-T150 variation across all four environments".
  type_within_env <- avg_simple_effects_manual(
    m0d, newdata = demo_grid, variable = "UASType", within_variable = "AmbientEnvCore",
    levels = c("T150", "H520"),
    re_formula = NULL, allow_new_levels = TRUE, sample_new_levels = "gaussian",
    ndraws = ndraws_pred, seed = 999,
    rope = c(-0.2, 0.2),
    transform = \(x) betaTransform(x, -10, 10, direction = "reverse", squeeze = "none")
  )
  # type_within_env has 4 rows, one per AmbientEnvCore level - exactly the
  # shape for the plot you're after:
  ggplot(type_within_env, aes(x = AmbientEnvCore, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
    labs(x = NULL, y = "Predicted T150 - H520 difference in mean annoyance")
  
  # Sanity check: interaction contrasts should equal differences of simple
  # effects at the same two AmbientEnvCore levels (small numerical
  # differences expected, since each was computed from an independent set
  # of posterior/new-level draws - not a discrepancy to worry about unless
  # it's large)
  # type_within_env$estimate[type_within_env$AmbientEnvCore == "Park"] -
  #   type_within_env$estimate[type_within_env$AmbientEnvCore == "Highway"]
  # vs. type_by_env$estimate[1]  # the Park-vs-Highway interaction contrast row
  
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