# %% [markdown]
# # REFMAP laboratory listening test 1 analysis: Rough confirmatory data analysis
#
# ## Setup
#
# ### Conversion notes (dabest fork → dabest "Bingka" v2025.10.20)
#
# This notebook was migrated from an older, forked dabest API to the current
# mainline `dabest` (v2025.10.20), matching the style used in the
# `refmap_listest2_analysis_dabest` notebook. The following systematic kwarg
# renames were applied to every `.mean_diff.plot(...)` call:
#
# | Old kwarg (fork)          | New kwarg                         | Notes |
# |----------------------------|------------------------------------|-------|
# | `swarm_ylim`               | `raw_ylim`                         | direct rename |
# | `es_marker_size`           | `contrast_marker_size`             | direct rename |
# | `swarm_desat`               | `raw_desat`                        | direct rename |
# | `halfviolin_desat`          | `contrast_desat`                   | direct rename |
# | `contrast_show_deltas`     | `delta_dot`                        | controls whether individual paired-difference dots are drawn on the contrast axis |
# | `contrast_show_es`         | `delta_text`                       | controls whether the numeric effect-size value is annotated on the contrast axis |
# | `es_sf`                     | *(removed)*                        | delta-text precision is now fixed (2 significant figures); no longer configurable |
# | `bar_desat`                 | *(removed)*                        | no longer a direct top-level kwarg; the old "gapped line" group-summary desaturation doesn't have a confirmed 1:1 replacement — restyle via `group_summaries_kwargs` if needed |
# | `slopegraph_xjitter`        | folded into `slopegraph_kwargs['jitter']` | new API takes a single scalar `jitter` rather than separate x/y jitter |
# | `slopegraph_yjitter`        | *(removed)*                        | was always `0.0` in this notebook, so nothing is lost |
# | `jitter_seed`                | *(removed)*                        | no longer a distinct kwarg; reproducibility now relies on the global `np.random.seed(303)` and `dabest.load(..., random_seed=...)` already set in this notebook |
#
# To align this notebook's look with the newer `listest2` notebook's style:
# - `contrast_bars=False` is added to every `.plot()` call, matching
#   `listest2`'s convention (the new default is `True`, which would draw a
#   background bar behind each contrast/effect-size marker that wasn't part
#   of the original design).
# - `raw_bars` is left at its new default (`True`) — `listest2` never
#   overrides it, so the raw swarm points now get a background summary bar
#   too. This is a genuine visual change from the original fork's plots.
# - `show_baseline_ec=True` is added to every plot whose underlying
#   `dabest.load(...)` call used `paired='baseline'` (15 of the 46 plots in
#   this notebook), matching `listest2`'s convention for baseline-paired
#   designs. It is *not* added to `paired='sequential'` or unpaired plots,
#   since it isn't meaningful there.
#
# A matplotlib API fix was also applied throughout: `legend.legendHandles`
# (deprecated/removed in recent Matplotlib) → `legend.legend_handles`.
#
# ### Second round of updates (from hands-on testing of the first cell)
#
# - **`custom_palette` + `color_col` don't mix reliably** in the current
#   version. Every call that combined `color_col=...` with `custom_palette=`
#   has had `custom_palette` removed, so colouring is driven by `color_col`
#   alone (dabest's own default palette). Calls that use `custom_palette` on
#   its own (a `{group: colour}` dict keyed to the x-axis categories, with no
#   `color_col`) are unaffected and still use the `mycolours` palette.
# - **Delta dots are enabled everywhere**, matching the new default: every
#   `delta_dot=False` was replaced with
#   `delta_dot_kwargs={'size': 1, 'side': 'left', 'alpha': 0.25, 'zorder': 1}`.
#   This only has a visible effect on paired (slopegraph) plots; it's a
#   harmless no-op on unpaired ones. The `'size': 1` default may be worth
#   tuning down (e.g. `0.75`, as used in the first cell) on any plot where
#   dots feel crowded.
# - **The old `axs[0].findobj(mpl.text.Annotation)` hack is gone.** It was
#   specific to the very first plot in this notebook and has been replaced
#   with the same `ax.contrast_axes` + `display_round`-based results-labelling
#   pattern used everywhere else (see below), following `listest2`'s style.
# - **Every "label results" block now uses `contrast_ax = ax.contrast_axes`**
#   (rather than plotting result text on the raw-data axes) and formats the
#   95% CI with `display_round(..., digits=2, floor=False)`, plus a new
#   permutation-test *p*-value line, exactly mirroring `listest2`.
#   Correction: the *p*-value should be read from `diff[-1]`, not `diff[3]`
#   (the `statistical_tests` column layout returned by the installed dabest
#   version apparently differs slightly from what `listest2` assumed — a
#   fixed index of `3` isn't reliable, but the p-value is consistently the
#   *last* column). All blocks — including the first cell — now use
#   `display_round(diff[-1], digits=3)`, confirmed visually correct.
#   The vertical gap between the "95%" line and the new *p* line is set
#   automatically per plot as roughly 9.5% of that plot's `contrast_ylim`
#   span (matching the ratio used in the hand-tuned first cell), and the *p*
#   line's x-position is nudged `+0.09` past the 95% line's x, again mirroring
#   the first cell. These are reasonable starting points, not verified
#   renders — some plots will likely need their exact x/y nudged by hand,
#   the same way the first cell was tuned.
# - Added `from refmap_psychoacoustics.utils.format_helpers import
#   (round_trad, display_round)` to the imports, matching `listest2`.
#
# ### Third round of updates: dummy pairing-ID creation
#
# All 9 places that built the `dummyID` pairing column with
# `np.repeat(np.arange(...), ...)` have been switched to the
# `groupby(...).ngroup() + 1` approach used in `listest2`, e.g.:
#
# ```python
# grouping_cols = ['ID', 'UASType', 'UASOperation']
# data['dummyID'] = data.groupby(grouping_cols).ngroup() + 1
# ```
#
# The old `np.repeat` approach only worked because the dataframe had
# already been sorted so that every block of `len(x.unique())` consecutive
# rows belonged to one participant/condition group — fragile, since any
# change to the upstream sort order or a missing row would silently
# mis-pair data. `groupby().ngroup()` doesn't care about row order at all;
# it assigns a shared integer to every row sharing the same combination of
# `grouping_cols`, regardless of how the dataframe is sorted, so it can't
# get thrown off in the same way.
#
# `grouping_cols` was chosen per block as *every retained column except the
# paired (`x=`) variable, the outcome (`y=`), the raw stimulus-file column,
# and `dummyID` itself* — i.e. everything that should stay fixed across one
# participant's set of paired observations. A `.sort_values(by=['dummyID',
# 'ID', ...])` was added after each new assignment purely for readability
# when inspecting the dataframe; it has no bearing on correctness now.
#
# Two subtleties worth knowing about:
# - In the "Flight operation" and "UAS type and operation combined" blocks,
#   the full (`Park`+`Street`) `data` frame additionally varies over
#   `AmbientEnv`, so its `grouping_cols` includes `AmbientEnv` on top of the
#   columns used for the Park-only/Street-only subsets (which have it fixed
#   already and so don't need it).
# - In the Part B "Number of events" block, a single `grouping_cols` list
#   (`['ID', 'UASType', 'UASLAeq']`) is reused across a `for` loop over nine
#   different filtered dataframes. This is safe because grouping by a column
#   that happens to be constant within a given subset (e.g. `UASType` in
#   `dataH520`) just doesn't split the groups any further — it's a no-op,
#   not an error — so the same grouping list works whether or not that
#   column still varies in a particular subset.
#
# ### Fourth round: a real grouping bug in the first block
#
# The "Scaled vs modelled UAS LAeq level 2" block's `grouping_cols` was
# missing `UASType`. Its subset (`StimFile` containing `"_F_2"`) spans three
# UAS types (H520, M300, T150) as well as two ambient environments, so
# grouping by `['ID', 'AmbientEnv']` alone collapsed all three UAS types'
# "Modelled" rows for a given participant/environment into one group, and
# likewise for "Scaled" — silently mis-pairing data across UAS types rather
# than pairing each UAS type's own Modelled/Scaled reading. Fixed by adding
# `'UASType'` to both the `.loc[:, [...]]` column selection (it wasn't being
# retained at all) and `grouping_cols`.
#
# The other 8 `dummyID` blocks were re-checked against the same logic
# (produce/select every design factor that actually varies within that
# specific subset, besides `ID`, the paired `x=` variable, and the outcome)
# and are correct as they stand — in particular, Part B's `UASEvents` blocks
# don't need `AmbientEnv`/`UASOperation` in the grouping since those two
# factors don't vary within Part B at all (Part B testing was Park/Overflight
# only), and the "UAS type and operation combined" block deliberately
# excludes `UASType`/`UASOperation` from grouping since both are already
# fully consumed by the paired `UASTypeOp` variable.
#
# ### Fifth round: unused categories on categorical columns
#
# `AmbientEnv`, `UASLAeq`, `UASOperation` and `UASType` are all set up as
# `pd.Categorical` early in the notebook, with a fixed category list that
# includes values (e.g. `"Baseline"`, or LAeq/UASType/operation levels that
# don't apply to a given subset) which don't actually occur once the data has
# been filtered down for a particular comparison. Left as-is, this causes two
# distinct problems: `groupby()` on a categorical column doesn't limit itself
# to the categories that actually occur in the data, which distorts
# `.ngroup()`'s pairing; and separately, `dabest.load()` and `.plot()` build
# their own ordered categorical for the `x=`/`color_col=` column by calling
# `.cat.reorder_categories(all_plot_groups, ...)`, which raises a
# `ValueError` unless the column's *existing* categories are already an exact
# match for the groups being plotted.
#
# The first fix attempt used `.cat.remove_unused_categories()`, matching the
# groups-still-not-matching approach in `listest2`. In testing this turned
# out not to be reliable enough to avoid the second problem — dabest's own
# `.cat.remove_unused_categories()` call inside `_get_plot_data` doesn't
# actually take effect (its result isn't assigned back to anything), so
# whatever categories survive by the time `dabest.load()` runs still have to
# match exactly, and that match turned out to be fragile in practice.
#
# The more robust fix — now used throughout — is to convert these columns to
# plain string (`.astype(str)`) rather than trying to keep them as a trimmed
# `pd.Categorical`. This sidesteps dabest's category-reordering code path
# entirely: for a non-categorical column, dabest just builds a fresh ordered
# categorical directly from the groups being plotted, with no dependency on
# whatever categories happened to survive upstream. It also solves the
# `groupby()`/`.ngroup()` problem the same way, since grouping on a plain
# string column was never affected by unused categories in the first place.
#
# ```python
# for dataset in [dataPk, dataSt]:
#     dataset['UASLAeq'] = dataset['UASLAeq'].astype(str)
#     dataset['UASType'] = dataset['UASType'].astype(str)
#     dataset['UASOperation'] = dataset['UASOperation'].astype(str)
# ```
#
# This was applied to every column that's actually used as `x=`, `color_col=`
# or in `grouping_cols` in each block: the "Scaled vs modelled" block, both
# "segregated by ambient environment" blocks, "Flight operation" (for the
# combined `data` as well as the Park/Street subsets — including
# `UASOperation`, since it's the `x=` variable there), "UAS type and
# operation combined" (`UASLAeq` only, since `UASType`/`UASOperation` are
# fully consumed by the derived `UASTypeOp` and never used directly), and the
# Part B "Number of events" loop (`UASLAeq` and `UASType`, converted inside
# the loop itself so it applies fresh to each of the nine filtered subsets —
# this also covers every further-filtered sub-block downstream, e.g.
# `dataPk.loc[dataPk['UASLAeq'] == '60']`, since a `.loc[]` filter on an
# already-string column stays a string column).
#
# One side effect worth knowing about: a plain string column doesn't carry
# any particular display order of its own. This doesn't affect axis order
# for any `x=` variable, since `idx=(...)` always specifies that explicitly
# regardless of the column's dtype. It could in principle affect legend
# order for a `color_col=` variable that isn't also constrained elsewhere,
# but in every affected block here the values happen to sort into a sensible
# order anyway (e.g. `"42"`, `"48"`, `"54"`, `"60"`; `"H520"`, `"M300"`,
# `"T150"`) — worth a glance at the legends to confirm they still read
# top-to-bottom (or left-to-right) the way you expect.
#
# The DroneNoise blocks at the end don't need this: their `UASOperation`,
# `Location` and `UASType` columns are built with plain string assignment,
# not `pd.Categorical`, so there's no fixed category list to go stale. The
# unpaired comparisons (AAM attitude, Home residence area, Nationality,
# Area soundscape) also don't need it: they don't group on any of the
# affected columns, and their own `x=` variables (`AAM_attitude`,
# `Home_Area`, etc.) were never made categorical in the first place.
#
# ### Sixth round: Park/Street LAeq-by-type styling, from hands-on testing
#
# The "UAS LAeq segregated by ambient environment" Park and Street plots
# (`dataPkloadSQ`/`dataStloadSQ`, sequential-paired, `color_col='UASType'`)
# were restyled based on what actually rendered well: `contrast_bars=True`
# (was `False`), `delta_dot=False` instead of styling the dots via
# `delta_dot_kwargs` (they were too cluttered with 4 LAeq levels ×
# 3 UAS types), `slopegraph_kwargs` alpha raised `0.15→0.45` and jitter
# raised `0.03→0.3`, and the result-label text moved from `y=-4.5`/`y=-4.69`
# (both well outside this plot's `contrast_ylim=(-0.25, 1.75)`, so
# effectively invisible/off-plot) to `y=1.55`/`y=1.4`, which sit inside it.
#
# **This `y=-4.5`/`y=-4.69` mispositioning is not unique to this pair.** The
# same literal values were carried through unchanged from the original
# fork notebook into every "label results" block whose `ax.text()` call
# used exactly these numbers, and checking each block's own `contrast_ylim`
# shows `-4.5` falls outside the visible range in all of them:
#
# | Variable | `contrast_ylim` |
# |---|---|
# | `dataPkloadBL`, `dataPkHiloadBL` | `(-2, 0.5)` |
# | `dataloadBL`, `dataPkloadBL`, `dataStloadBL` (Flight operation) | `(-1, 1.5)` |
# | `dataloadBL`, `dataPkloadBL`, `dataStloadBL` (UAS type+operation) | `(-1, 2)` |
# | `dataloadBL` (Part B baseline) | `(0, 2.5)` |
# | `dataloadSQ` | `(-0.5, 2)` |
# | `dataLoloadSQ`, `dataHiloadSQ` | `(-1, 2)` / `(-0.5, 2)` |
# | `dataLoH520loadSQ`, `dataHiH520loadSQ`, `dataLoT150loadSQ`, `dataHiT150loadSQ`, `dataH520loadSQ`, `dataT150loadSQ` | `(-1.5, 2.5)` |
#
# This looks like a pre-existing copy/paste artifact from the original fork
# notebook (the exact same literal recurring across dozens of unrelated
# plots) rather than something introduced by the API migration, and it was
# carried through faithfully since the original y-position was preserved
# during conversion. I haven't guessed new positions for these — since I
# can't render the plots to check, a guess risks being wrong in the same way
# the very first attempt at this pair was. Worth working through them the
# same way as this pair (nudge `y` until the text sits inside `contrast_ylim`
# and reads clearly), or let me know if you'd like a best-effort automatic
# pass at repositioning them as a starting point.
#
# ### Seventh round: `ps_adjust=True` on every `dabest.load(...)` call
#
# `ps_adjust` — originally your own fork contribution, since folded into
# the mainline package as an optional parameter — was added to all 47
# `dabest.load(...)` calls in the notebook.
#
# ### Eighth round: slopegraph/delta-dot styling + annotation repositioning
# rolled out to the rest of the paired (baseline/sequential) plots
#
# Following the styling change confirmed on the "UAS type segregated by
# ambient environment" Park plot — `delta_dot=False` (rather than styling
# the dots via `delta_dot_kwargs`), `slopegraph_kwargs` alpha raised to
# `0.45` and jitter to `0.3`, and a `delta_text_kwargs={'y_coordinates':
# [...]}` override for `paired='baseline'` plots — this was applied to
# every remaining plot flagged in the previous round's `y=-4.5` table:
# both "UAS type segregated by ambient environment" sub-blocks (Park,
# Street, Park-Hi, Street-Hi — all share `contrast_ylim=(-2, 0.5)`, so the
# exact confirmed `y=-1`/`y=-1.2` position was reused across all four),
# "Flight operation" (`data`/`dataPk`/`dataSt`/Hi, `contrast_ylim=(-1,
# 1.5)`), "UAS type and operation combined" (`data`/`dataPk`/`dataSt`,
# `contrast_ylim=(-1, 2)`), and every Part B "Number of events" block
# (`dataloadBL`, `dataloadSQ`, `dataLoloadSQ`, `dataHiloadSQ`, the four
# `H520`/`T150` × `54`/`60` combinations, and the H520/T150-only blocks).
#
# **Confidence differs across these.** The Park/Street/Hi group in the
# first section shares the exact `contrast_ylim` the fix was confirmed
# against, so those four reuse the tested `y` values directly. Everywhere
# else, I don't have a confirmed position to copy — I estimated a new `y`
# for each block's own `contrast_ylim` using the same *proportions* as the
# confirmed fix (roughly 40% up from the bottom of the range for the "95%"
# line, with the *p*-line ~8–9.5% of the range's span below it, and the
# `delta_text_kwargs` y-coordinate around 60% up), rather than copying a
# number that was tuned for a different range. These are still starting
# points, not verified renders — expect some of them to need a further
# nudge once you can see them.
#
# Two blocks (`dataHiloadSQ`, and the `if ii == 0: / else:` blocks for the
# H520/T150 sub-comparisons) previously used a *different*, even-more-out-
# of-range `y` specifically for the first (`ii == 0`) comparison. I couldn't
# find a reason for that split that still applies now, so these were
# simplified to use one consistent position for every comparison, matching
# the other blocks in the same section — flagging this in case the original
# split was deliberate and you'd rather keep it.
#
# The two DroneNoise plots at the very end of the notebook still use the
# older `delta_dot_kwargs`/lower slopegraph alpha styling — they don't have
# the same `y=-4.5` annotation problem (no manual `contrast_ax.text()` in
# those cells), so they were left as-is, but let me know if you'd like the
# same `delta_dot=False`/`alpha=0.45`/`jitter=0.3` styling applied there too
# for consistency.
#
# ### Ninth round: `ps_adjust` re-verified, `contrast_bars=True` extended to
# every sequential-paired plot
#
# `ps_adjust=True` was double-checked with an AST parse of the notebook
# (rather than a text search, which can be fooled by the word appearing in
# a comment) — all 47 `dabest.load(...)` calls have it. If your local copy
# still shows some missing it, it's likely out of sync with this file —
# worth re-pulling this version.
#
# `contrast_bars=True` was applied to the *rest* of the `paired='sequential'`
# plots in Part B ("Number of events" and all its LAeq/UAS-type-segregated
# sub-blocks — `dataloadSQ`, `dataLoloadSQ`, `dataHiloadSQ`,
# `dataLoH520loadSQ`, `dataHiH520loadSQ`, `dataLoT150loadSQ`,
# `dataHiT150loadSQ`, `dataH520loadSQ`, `dataT150loadSQ`). The eighth round
# only carried this over to the very first confirmed pair
# (`dataPkloadSQ`/`dataStloadSQ`) and left every other sequential block at
# `contrast_bars=False`, which was an oversight — all of Part B's sequential
# plots are the same kind of plot and should behave the same way.
#
# `paired='baseline'` plots (the ones with `show_baseline_ec=True`) were
# **not** changed — the confirmed "UAS type segregated by ambient
# environment" fix explicitly kept `contrast_bars=False` for those, so that
# distinction (sequential → `True`, baseline → `False`) is preserved
# throughout.
#
# Everything else (the statistical/data-wrangling logic, `dabest.load(...)`
# calls, and figure-saving code) is untouched.

# %%
# import packages
import sys
import os
import numpy as np
import pandas as pd
from PyQt5.QtWidgets import QFileDialog, QApplication
from scipy import stats
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.collections import PathCollection
from matplotlib.legend_handler import HandlerPathCollection, HandlerLine2D
import seaborn as sns
import dabest
import warnings
from refmap_psychoacoustics.utils.format_helpers import (round_trad, display_round)

# Suppress FutureWarning messages to quiet pandas
warnings.simplefilter(action='ignore', category=FutureWarning)


# %%
# set plot parameters
sns.set_style('white')
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times New Roman'
mpl.rcParams.update({'font.size': 16})
mpl.rcParams['figure.autolayout'] = True
mpl.rcParams['mathtext.fontset'] = 'stix'

SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE,
       labelsize=MEDIUM_SIZE)    # fontsize of the axes title and x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

mycolours = [(0, 102, 255), (0, 204, 153), (255, 0, 102), (74, 111, 152),
             (251, 164, 49), (204, 153, 255), (90, 192, 255), (80, 245, 233),
             (255, 90, 192), (164, 201, 242), (255, 254, 139), (255, 243, 255)]
mycolours = [tuple(shade/255 for shade in colour) for colour in mycolours]

np.random.seed(303)

# enable copy-on-write mode for Pandas (will be default from Pandas 3.0)
pd.options.mode.copy_on_write = True

# check/open QApplication instance
if not QApplication.instance():
    app = QApplication(sys.argv)
else:
    app = QApplication.instance() 


# %% [markdown]
# Set the saveplots toggle to True if plot saving is desired:

# %%
saveplots = True

if saveplots:
    # select figure output save path
    outFigPath = QFileDialog.getExistingDirectory(caption=r"Select output folder to save plots in: 03 Experiment\Experiment 1\Analysis\Plots")

    # create subfolders if not already existing
    try:
        os.mkdir(os.path.join(outFigPath, "svg"))
    except FileExistsError:
        pass

    try:
        os.mkdir(os.path.join(outFigPath, "pdf"))
    except FileExistsError:
        pass


# %% [markdown]
# ## Import data and organise

# %%
# import test data
fileExts = "*.csv"

# Part A
dataBySubjAFilePath = list(QFileDialog.getOpenFileName(filter="refmap_listest1_testdataA_BySubj.csv",
                                                       caption=r"Open refmap_listest1_testdataA_BySubj.csv in: \03 Experiment\Experiment 1\Analysis\PostProcess"))[0]
dataBySubjTestA = pd.read_csv(dataBySubjAFilePath)

# Part B
dataBySubjBFilePath = list(QFileDialog.getOpenFileName(filter="refmap_listest1_testdataB_BySubj.csv",
                                                       caption=r"Open refmap_listest1_testdataB_BySubj.csv in: \03 Experiment\Experiment 1\Analysis\PostProcess"))[0]
dataBySubjTestB = pd.read_csv(dataBySubjBFilePath)

# Both parts
dataBySubjFilePath = list(QFileDialog.getOpenFileName(filter="refmap_listest1_testdata_BySubj.csv",
                                                       caption=r"Open refmap_listest1_testdata_BySubj.csv in: \03 Experiment\Experiment 1\Analysis\PostProcess"))[0]
dataBySubjTest = pd.read_csv(dataBySubjFilePath)


# %%
# categorise columns

for dataset in [dataBySubjTestA, dataBySubjTest]:
    dataset['AmbientEnv'] = pd.Categorical(dataset['AmbientEnv'], ["Park", "Street"])
    dataset['SNRlevel'] = pd.Categorical(dataset['SNRlevel'], ["Baseline", "-16", "-10", "-4", "2", "8"], ordered=True)
    dataset['UASLAeq'] = pd.Categorical(dataset['UASLAeq'], ["Baseline", "42", "48", "54", "60"], ordered=True)
    dataset['UASOperation'] = pd.Categorical(dataset['UASOperation'], ["Baseline", "Overflight", "Landing", "Takeoff"])
    dataset['UASType'] = pd.Categorical(dataset['UASType'], ["Baseline", "H520", "M300", "T150"])

for dataset in [dataBySubjTestB]:
    dataset['AmbientEnv'] = pd.Categorical(dataset['AmbientEnv'], ["Park", "Street"])
    dataset['SNRlevel'] = pd.Categorical(dataset['SNRlevel'], ["Baseline", "2", "8"], ordered=True)
    dataset['UASLAeq'] = pd.Categorical(dataset['UASLAeq'], ["Baseline", "54", "60"], ordered=True)
    dataset['UASOperation'] = pd.Categorical(dataset['UASOperation'], ["Baseline", "Overflight"])
    dataset['UASType'] = pd.Categorical(dataset['UASType'], ["Baseline", "H520", "T150"])


# %% [markdown]
# ## Setup output dataframes

# %%
out = pd.DataFrame()
outd = pd.DataFrame()

savedata = True

if savedata:
    # select data output save path
    outDataPath = QFileDialog.getExistingDirectory(caption=r"Select output folder to save output data in: \03 Experiment\Experiment 1\Analysis\Python")


# %% [markdown]
# ## Part A
# 
# ### Scaled vs modelled UAS LAeq level 2 

# %%
# select subset of data for analysis and sort
data = dataBySubjTestA[dataBySubjTestA['StimFile'].str.contains("_F_2")]
data = data.loc[:, ['ID', 'StimFile', 'UASType', 'AmbientEnv', 'Annoyance']]
data.sort_values(by=['ID', 'StimFile'], inplace=True)

# drop the "Baseline" category and convert to plain string columns. These
# are also the columns used for x=/color_col=/grouping_cols downstream, and
# dabest's own category-reordering code (inside dabest.load()/  .plot())
# turns out to be unreliable once a categorical column has had categories
# removed -- converting away from pandas' Categorical dtype entirely
# sidesteps that fragility, since dabest then just builds its own ordered
# categorical fresh from the values actually present.
data['UASType'] = data['UASType'].astype(str)
data['AmbientEnv'] = data['AmbientEnv'].astype(str)

# add column to indicate modelled or scaled, and create a dummy pairing ID
data['LvlType'] = "Modelled"
data.loc[data['StimFile'].str.contains("PwrScale"), 'LvlType'] = "Scaled"

grouping_cols = ['ID', 'UASType', 'AmbientEnv']

# ngroup() automatically assigns an identical integer to matching sets
# across the different level-type conditions
data['dummyID'] = data.groupby(grouping_cols).ngroup() + 1
data.sort_values(by=['dummyID', 'ID', 'StimFile'], inplace=True)



# %%
# assign data for processing
dataloadBL = dabest.load(ps_adjust=True, data=data, idx=("Modelled", "Scaled"), x='LvlType', y='Annoyance', paired='baseline',
                         id_col='dummyID', resamples=5000, random_seed=303)


# %%
# calculate effect sizes and plot
fig, ax = plt.subplots(figsize=(9, 4.5))
dataMD = dataloadBL.mean_diff.plot(contrast_bars=False, show_baseline_ec=True,
                                   raw_ylim=(-0.5, 10.5), color_col='AmbientEnv',
                                   slopegraph_kwargs={'linewidth': 0.25, 'alpha': 0.6, 'jitter': 0.3},
                                   delta_dot_kwargs={'size': 0.75, 'side': 'left',
                                                     'alpha': 0.25, 'zorder': 1},
                                   delta_text=True, contrast_marker_size=4, contrast_ylim=(-5.5, 4.5),
                                   legend_kwargs={'bbox_to_anchor': [-1.5, 1.3], 'fontsize': 11,
                                                  'ncol': 2, 'title': r'Ambient env'}, ax=ax)

diffs = dataloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values

contrast_ax = ax.contrast_axes

# label results
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=ii + 0.45, y=5.4, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=ii + 0.54, y=4.45, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)


if saveplots:
    
    filename = "PtALvlScaleDabest"
    
    plt.savefig(os.path.join(outFigPath, "svg", filename + ".svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", filename + ".pdf"),
                format='pdf', bbox_inches='tight')

dataMD;
dataloadBL.mean_diff.statistical_tests


# %%
out = pd.concat([out, dataloadBL.mean_diff.statistical_tests])
# calculate and display Cohen's d
outd = pd.concat([outd, dataloadBL.cohens_d.statistical_tests])
dataloadBL.cohens_d.statistical_tests


# %% [markdown]
# ### UAS LAeq segregated by ambient environment

# %%
# select subset of data for analysis and sort
data = dataBySubjTestA[~dataBySubjTestA['StimFile'].str.contains("PwrScale")]
data = data[data['UASLAeq'] != "Baseline"]
data = data.loc[:, ['ID', 'StimFile', 'AmbientEnv', 'UASLAeq', 'UASType', 'UASOperation', 'Annoyance']]
data.sort_values(by=['ID', 'AmbientEnv', 'UASOperation', 'UASType', 'UASLAeq'], inplace=True)

dataPk = data.loc[data['AmbientEnv'] == "Park", :]
dataSt = data.loc[data['AmbientEnv'] == "Street", :]

# drop the "Baseline" category by converting to plain string columns
for dataset in [dataPk, dataSt]:
    dataset['UASLAeq'] = dataset['UASLAeq'].astype(str)
    dataset['UASType'] = dataset['UASType'].astype(str)
    dataset['UASOperation'] = dataset['UASOperation'].astype(str)

# create dummy pairing ID for each participant to allow for paired analysis
grouping_cols = ['ID', 'UASType', 'UASOperation']

# ngroup() automatically assigns an identical integer to matching sets
# across the different UAS LAeq conditions
dataPk['dummyID'] = dataPk.groupby(grouping_cols).ngroup() + 1
dataSt['dummyID'] = dataSt.groupby(grouping_cols).ngroup() + 1
dataPk.sort_values(by=['dummyID', 'ID', 'StimFile'], inplace=True)
dataSt.sort_values(by=['dummyID', 'ID', 'StimFile'], inplace=True)


# %%
# assign data for processing
dataPkloadSQ = dabest.load(ps_adjust=True, data=dataPk, idx=("42", "48", "54", "60"), x='UASLAeq', y='Annoyance', paired='sequential',
                           id_col='dummyID', resamples=5000, random_seed=303)
dataStloadSQ = dabest.load(ps_adjust=True, data=dataSt, idx=("42", "48", "54", "60"), x='UASLAeq', y='Annoyance', paired='sequential',
                           id_col='dummyID', resamples=5000, random_seed=808)

# %%
# calculate effect sizes and plot
fig, ax = plt.subplots(figsize=(6, 3))
dataPkMD = dataPkloadSQ.mean_diff.plot(contrast_bars=True,
                                   raw_ylim=(-0.5, 10.5), color_col='UASType',
                                       delta_dot=False,
                                       slopegraph_kwargs={'linewidth': 0.25, 'alpha': 0.45, 'jitter': 0.3},
                                       contrast_marker_size=4,
                                       contrast_ylim=(-0.25, 1.75),
                                       float_contrast=True,
                                       delta_text=True,
                                       legend_kwargs={'loc': 'upper center', 'title': r"UAS type", 'fontsize': 12,
                                                      'frameon': True, 'title_fontsize': 12}, ax=ax)
diffs = dataPkloadSQ.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=1.55, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.75 + ii + 0.09, y=1.4, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

lgd = ax.get_legend()
for lh in lgd.legend_handles:
    lh.set_alpha(0.8)
    lh.set_linewidth(3)

if saveplots:
    filename = "PtALAeqParkDabest"
    plt.savefig(os.path.join(outFigPath, "svg", filename + ".svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", filename + ".pdf"),
                format='pdf', bbox_inches='tight')

dataPkMD;
dataPkloadSQ.mean_diff.statistical_tests


# %%
out = pd.concat([out, dataPkloadSQ.mean_diff.statistical_tests])
# calculate effect sizes - Cohen's d
outd = pd.concat([outd, dataPkloadSQ.cohens_d.statistical_tests])
dataPkloadSQ.cohens_d.statistical_tests


# %%
# calculate effect sizes and plot
fig, ax = plt.subplots(figsize=(6, 3))
dataStMD = dataStloadSQ.mean_diff.plot(contrast_bars=True,
                                   raw_ylim=(-0.5, 10.5), color_col='UASType',
                                       delta_dot=False,
                                       slopegraph_kwargs={'linewidth': 0.25, 'alpha': 0.45, 'jitter': 0.3},
                                       contrast_marker_size=4,
                                       contrast_ylim=(-0.25, 1.75),
                                       float_contrast=True,
                                       delta_text=True,
                                       legend_kwargs={'loc': 'upper center', 'title': r"UAS type", 'fontsize': 12,
                                                      'frameon': True, 'title_fontsize': 12}, ax=ax)
diffs = dataStloadSQ.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=1.55, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.75 + ii + 0.09, y=1.4, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

lgd = ax.get_legend()
for lh in lgd.legend_handles:
    lh.set_alpha(0.8)
    lh.set_linewidth(3)

if saveplots:
    filename = "PtALAeqStreetDabest"
    plt.savefig(os.path.join(outFigPath, "svg", filename + ".svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", filename + ".pdf"),
                format='pdf', bbox_inches='tight')

dataStMD;
dataStloadSQ.mean_diff.statistical_tests


# %%
out = pd.concat([out, dataStloadSQ.mean_diff.statistical_tests])
# calculate effect sizes - Cohen's d
outd = pd.concat([outd, dataStloadSQ.cohens_d.statistical_tests])
dataStloadSQ.cohens_d.statistical_tests


# %% [markdown]
# ### UAS types and operations segregated by ambient environment
# 
# #### UAS type

# %%
# select subset of data for analysis and sort
data = dataBySubjTestA[~dataBySubjTestA['StimFile'].str.contains("PwrScale")]
data = data[data['UASType'] != "Baseline"]
data = data.loc[:, ['ID', 'StimFile', 'AmbientEnv', 'UASLAeq', 'UASType', 'UASOperation', 'Annoyance']]
data.sort_values(by=['ID', 'UASLAeq', 'UASOperation', 'UASType', 'AmbientEnv'], inplace=True)

dataPk = data.loc[data['AmbientEnv'] == "Park", :]
dataSt = data.loc[data['AmbientEnv'] == "Street", :]

# drop the "Baseline" category by converting to plain string columns
for dataset in [dataPk, dataSt]:
    dataset['UASLAeq'] = dataset['UASLAeq'].astype(str)
    dataset['UASType'] = dataset['UASType'].astype(str)
    dataset['UASOperation'] = dataset['UASOperation'].astype(str)

# create dummy pairing ID for each participant to allow for paired analysis
grouping_cols = ['ID', 'UASLAeq', 'UASOperation']

# ngroup() automatically assigns an identical integer to matching sets
# across the different UAS type conditions
dataPk['dummyID'] = dataPk.groupby(grouping_cols).ngroup() + 1
dataSt['dummyID'] = dataSt.groupby(grouping_cols).ngroup() + 1
dataPk.sort_values(by=['dummyID', 'ID', 'StimFile'], inplace=True)
dataSt.sort_values(by=['dummyID', 'ID', 'StimFile'], inplace=True)


# %%
# assign data for processing
dataPkloadBL = dabest.load(ps_adjust=True, data=dataPk, idx=("H520", "M300", "T150"), x='UASType', y='Annoyance', paired='baseline',
                           id_col='dummyID', resamples=5000, random_seed=808)
dataStloadBL = dabest.load(ps_adjust=True, data=dataSt, idx=("H520", "M300", "T150"), x='UASType', y='Annoyance', paired='baseline',
                           id_col='dummyID', resamples=5000, random_seed=808)


# %%
# calculate effect sizes and plot
fig, ax = plt.subplots(figsize=(5.14, 3))
dataPkMD = dataPkloadBL.mean_diff.plot(contrast_bars=False, show_baseline_ec=True,
                                        raw_ylim=(-0.5, 10.5), color_col='UASLAeq',
                                       delta_dot=False,
                                       slopegraph_kwargs={'linewidth': 0.25, 'alpha': 0.45, 'jitter': 0.3},
                                       contrast_marker_size=4,
                                       contrast_ylim=(-2, 0.5), float_contrast=True,
                                       delta_text=True,
                                       delta_text_kwargs={'y_coordinates': [-0.5]*(len(data['UASType'].unique()) - 1)},
                                       legend_kwargs={'loc': 'upper center', 'title': r"UAS $L_\mathrm{Aeq}$", 'fontsize': 12,
                                                      'frameon': True, 'title_fontsize': 12}, ax=ax)
diffs = dataPkloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=-1, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.74 + ii + 0.09, y=-1.2, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

lgd = ax.get_legend()
for lh in lgd.legend_handles:
    lh.set_alpha(0.8)
    lh.set_linewidth(3)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtATypeParkDabest.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtATypeParkDabest.pdf"),
                format='pdf', bbox_inches='tight')

dataPkMD;
dataPkloadBL.mean_diff.statistical_tests


# %%
out = pd.concat([out, dataPkloadBL.mean_diff.statistical_tests])
# calculate effect sizes - Cohen's d
outd = pd.concat([outd, dataPkloadBL.cohens_d.statistical_tests])
dataPkloadBL.cohens_d.statistical_tests


# %%
# calculate effect sizes and plot
fig, ax = plt.subplots(figsize=(5.14, 3))
dataStMD = dataStloadBL.mean_diff.plot(contrast_bars=False, show_baseline_ec=True,
                                        raw_ylim=(-0.5, 10.5), color_col='UASLAeq',
                                       delta_dot=False,
                                       slopegraph_kwargs={'linewidth': 0.25, 'alpha': 0.45, 'jitter': 0.3},
                                       contrast_marker_size=4,
                                       contrast_ylim=(-2, 0.5), float_contrast=True,
                                       delta_text=True,
                                       delta_text_kwargs={'y_coordinates': [-0.5]*(len(data['UASType'].unique()) - 1)},
                                       legend_kwargs={'loc': 'upper center', 'title': r"UAS $L_\mathrm{Aeq}$", 'fontsize': 12,
                                                      'frameon': True, 'title_fontsize': 12}, ax=ax)

diffs = dataStloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=-1, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.74 + ii + 0.09, y=-1.2, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

lgd = ax.get_legend()
for lh in lgd.legend_handles:
    lh.set_alpha(0.8)
    lh.set_linewidth(3)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtATypeStreetDabest.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtATypeStreetDabest.pdf"),
                format='pdf', bbox_inches='tight')


dataStMD;
dataStloadBL.mean_diff.statistical_tests


# %%
out = pd.concat([out, dataStloadBL.mean_diff.statistical_tests])
# calculate effect sizes - Cohen's d
outd = pd.concat([outd, dataStloadBL.cohens_d.statistical_tests])
dataStloadBL.cohens_d.statistical_tests


# %%
# assign data for processing
dataPkHiloadBL = dabest.load(ps_adjust=True, data=dataPk.loc[dataPk['UASLAeq'] == '60'], idx=("H520", "M300", "T150"), x='UASType', y='Annoyance', paired='baseline',
                             id_col='dummyID', resamples=5000, random_seed=12345)
dataStHiloadBL = dabest.load(ps_adjust=True, data=dataSt.loc[dataSt['UASLAeq'] == '60'], idx=("H520", "M300", "T150"), x='UASType', y='Annoyance', paired='baseline',
                             id_col='dummyID', resamples=5000, random_seed=999)


# %%
# calculate effect sizes and plot
fig, ax = plt.subplots(figsize=(5.14, 3))
dataPkHiMD = dataPkHiloadBL.mean_diff.plot(contrast_bars=False, show_baseline_ec=True,
                                   raw_ylim=(-0.5, 10.5), color_col='UASLAeq',
                                           delta_dot=False,
                                           slopegraph_kwargs={'linewidth': 0.25, 'alpha': 0.45, 'jitter': 0.3},
                                           contrast_marker_size=4, contrast_ylim=(-2, 0.5), float_contrast=True,
                                           delta_text=True,
                                           delta_text_kwargs={'y_coordinates': [-0.5]*(len(data['UASType'].unique()) - 1)},
                                           legend_kwargs={'loc': 'upper center', 'title': r"UAS $L_\mathrm{Aeq}$", 'fontsize': 12,
                                                          'frameon': True, 'title_fontsize': 12}, ax=ax)
diffs = dataPkHiloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=-1, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.74 + ii + 0.09, y=-1.2, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtATypeParkHiDabest.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtATypeParkHiDabest.pdf"),
                format='pdf', bbox_inches='tight')

dataPkHiMD;
dataPkHiloadBL.mean_diff.statistical_tests


# %%
out = pd.concat([out, dataPkHiloadBL.mean_diff.statistical_tests])
# calculate effect sizes - Cohen's d
outd = pd.concat([outd, dataPkHiloadBL.cohens_d.statistical_tests])
dataPkHiloadBL.cohens_d.statistical_tests


# %%
# calculate effect sizes and plot
fig, ax = plt.subplots(figsize=(5.14, 3))
dataStHiMD = dataStHiloadBL.mean_diff.plot(contrast_bars=False, show_baseline_ec=True,
                                   raw_ylim=(-0.5, 10.5), color_col='UASLAeq',
                                           delta_dot=False,
                                           slopegraph_kwargs={'linewidth': 0.25, 'alpha': 0.45, 'jitter': 0.3},
                                           contrast_marker_size=4, contrast_ylim=(-2, 0.5), float_contrast=True,
                                           delta_text=True,
                                           delta_text_kwargs={'y_coordinates': [-0.5]*(len(data['UASType'].unique()) - 1)},
                                           legend_kwargs={'loc': 'upper center', 'title': r"UAS $L_\mathrm{Aeq}$", 'fontsize': 12,
                                                          'frameon': True, 'title_fontsize': 12}, ax=ax)
diffs = dataStHiloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=-1, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.74 + ii + 0.09, y=-1.2, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtATypeStreetHiDabest.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtATypeStreetHiDabest.pdf"),
                format='pdf', bbox_inches='tight')
out = pd.concat([out, dataStHiloadBL.mean_diff.statistical_tests])

dataStHiMD;
dataStHiloadBL.mean_diff.statistical_tests


# %%
out = pd.concat([out, dataStHiloadBL.mean_diff.statistical_tests])
# calculate effect sizes - Cohen's d
outd = pd.concat([outd, dataStHiloadBL.cohens_d.statistical_tests])
dataStHiloadBL.cohens_d.statistical_tests


# %% [markdown]
# #### Flight operation

# %%
# select subset of data for analysis and sort
data = dataBySubjTestA[~dataBySubjTestA['StimFile'].str.contains("PwrScale")]
data = data[data['UASType'] != "Baseline"]
data = data.loc[:, ['ID', 'StimFile', 'AmbientEnv', 'UASLAeq', 'UASType', 'UASOperation', 'Annoyance']]
data.sort_values(by=['ID', 'UASLAeq', 'UASType', 'AmbientEnv', 'UASOperation'], inplace=True)

dataPk = data.loc[data['AmbientEnv'] == "Park", :]
dataSt = data.loc[data['AmbientEnv'] == "Street", :]

# drop the "Baseline" category by converting to plain string columns
data['UASLAeq'] = data['UASLAeq'].astype(str)
data['UASType'] = data['UASType'].astype(str)
data['UASOperation'] = data['UASOperation'].astype(str)
for dataset in [dataPk, dataSt]:
    dataset['UASLAeq'] = dataset['UASLAeq'].astype(str)
    dataset['UASType'] = dataset['UASType'].astype(str)
    dataset['UASOperation'] = dataset['UASOperation'].astype(str)

# create dummy pairing ID for each participant to allow for paired analysis
grouping_cols_all = ['ID', 'AmbientEnv', 'UASLAeq', 'UASType']
grouping_cols_env = ['ID', 'UASLAeq', 'UASType']

# ngroup() automatically assigns an identical integer to matching sets
# across the different flight-operation conditions
data['dummyID'] = data.groupby(grouping_cols_all).ngroup() + 1
dataPk['dummyID'] = dataPk.groupby(grouping_cols_env).ngroup() + 1
dataSt['dummyID'] = dataSt.groupby(grouping_cols_env).ngroup() + 1
data.sort_values(by=['dummyID', 'ID', 'StimFile'], inplace=True)
dataPk.sort_values(by=['dummyID', 'ID', 'StimFile'], inplace=True)
dataSt.sort_values(by=['dummyID', 'ID', 'StimFile'], inplace=True)


# %%
# assign data for processing
dataloadBL = dabest.load(ps_adjust=True, data=data, idx=("Overflight", "Landing", "Takeoff"), x='UASOperation', y='Annoyance',
                         paired='baseline', id_col='dummyID', resamples=5000, random_seed=35941)
dataPkloadBL = dabest.load(ps_adjust=True, data=dataPk, idx=("Overflight", "Landing", "Takeoff"), x='UASOperation', y='Annoyance',
                           paired='baseline', id_col='dummyID', resamples=5000, random_seed=4444)
dataStloadBL = dabest.load(ps_adjust=True, data=dataSt, idx=("Overflight", "Landing", "Takeoff"), x='UASOperation', y='Annoyance',
                           paired='baseline', id_col='dummyID', resamples=5000, random_seed=6564)


# %%
# calculate effect sizes and plot
fig, ax = plt.subplots(figsize=(5.14, 3))
dataMD = dataloadBL.mean_diff.plot(contrast_bars=False, show_baseline_ec=True,
                                   raw_ylim=(-0.5, 10.5), color_col='UASLAeq',
                                   delta_dot=False,
                                   slopegraph_kwargs={'linewidth': 0.25, 'alpha': 0.45, 'jitter': 0.3},
                                   contrast_marker_size=4, contrast_ylim=(-1, 1.5), float_contrast=True,
                                   delta_text=True,
                                   delta_text_kwargs={'y_coordinates': [0.5]*(len(data['UASOperation'].unique()) - 1)},
                                   legend_kwargs={'loc': 'upper center', 'title': r"UAS $L_\mathrm{Aeq}$", 'fontsize': 12,
                                                  'frameon': True, 'title_fontsize': 12}, ax=ax)
diffs = dataloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=0.0, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.74 + ii + 0.09, y=-0.2, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

lgd = ax.get_legend()
for lh in lgd.legend_handles:
    lh.set_alpha(0.8)
    lh.set_linewidth(3)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtAOpDabest.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtAOpDabest.pdf"),
                format='pdf', bbox_inches='tight')

dataMD;
dataloadBL.mean_diff.statistical_tests


# %%
out = pd.concat([out, dataloadBL.mean_diff.statistical_tests])
# calculate effect sizes - Cohen's d
outd = pd.concat([outd, dataloadBL.cohens_d.statistical_tests])
dataloadBL.cohens_d.statistical_tests


# %%
# calculate effect sizes and plot
fig, ax = plt.subplots(figsize=(5.14, 3))
dataPkMD = dataPkloadBL.mean_diff.plot(contrast_bars=False, show_baseline_ec=True,
                                   raw_ylim=(-0.5, 10.5), color_col='UASLAeq',
                                       delta_dot=False,
                                       slopegraph_kwargs={'linewidth': 0.25, 'alpha': 0.45, 'jitter': 0.3},
                                       contrast_marker_size=4, contrast_ylim=(-1, 1.5), float_contrast=True,
                                       delta_text=True,
                                       delta_text_kwargs={'y_coordinates': [0.5]*(len(data['UASOperation'].unique()) - 1)},
                                       legend_kwargs={'loc': 'upper center', 'title': r"UAS $L_\mathrm{Aeq}$", 'fontsize': 12,
                                                      'frameon': True, 'title_fontsize': 12}, ax=ax)
diffs = dataPkloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=0.0, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.74 + ii + 0.09, y=-0.2, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

lgd = ax.get_legend()
for lh in lgd.legend_handles:
    lh.set_alpha(0.8)
    lh.set_linewidth(3)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtAOpParkDabest.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtAOpParkDabest.pdf"),
                format='pdf', bbox_inches='tight')

dataPkMD;
dataPkloadBL.mean_diff.statistical_tests


# %%
out = pd.concat([out, dataPkloadBL.mean_diff.statistical_tests])
# calculate effect sizes - Cohen's d
outd = pd.concat([outd, dataPkloadBL.cohens_d.statistical_tests])
dataPkloadBL.cohens_d.statistical_tests


# %%
# calculate effect sizes and plot
fig, ax = plt.subplots(figsize=(5.14, 3))
dataStMD = dataStloadBL.mean_diff.plot(contrast_bars=False, show_baseline_ec=True,
                                   raw_ylim=(-0.5, 10.5), color_col='UASLAeq',
                                       delta_dot=False,
                                       slopegraph_kwargs={'linewidth': 0.25, 'alpha': 0.45, 'jitter': 0.3},
                                       contrast_marker_size=4, contrast_ylim=(-1, 1.5), float_contrast=True,
                                       delta_text=True,
                                       delta_text_kwargs={'y_coordinates': [0.5]*(len(data['UASOperation'].unique()) - 1)},
                                       legend_kwargs={'loc': 'upper center', 'title': r"UAS $L_\mathrm{Aeq}$", 'fontsize': 12,
                                                      'frameon': True, 'title_fontsize': 12}, ax=ax)

diffs = dataStloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=0.0, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.74 + ii + 0.09, y=-0.2, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

lgd = ax.get_legend()
for lh in lgd.legend_handles:
    lh.set_alpha(0.8)
    lh.set_linewidth(3)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtAOpStreetDabest.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtAOpStreetDabest.pdf"),
                format='pdf', bbox_inches='tight')

dataStMD;
dataStloadBL.mean_diff.statistical_tests


# %%
out = pd.concat([out, dataStloadBL.mean_diff.statistical_tests])    
# calculate effect sizes - Cohen's d
outd = pd.concat([outd, dataStloadBL.cohens_d.statistical_tests])
dataStloadBL.cohens_d.statistical_tests


# %%
# assign data for processing
dataHiloadBL = dabest.load(ps_adjust=True, data=data.loc[data['UASLAeq'] == '60'], idx=("Overflight", "Landing", "Takeoff"),
                           x='UASOperation', y='Annoyance', paired='baseline',
                           id_col='dummyID', resamples=5000, random_seed=808)


# %%
# calculate effect sizes and plot
fig, ax = plt.subplots(figsize=(5.14, 3))
dataHiMD = dataHiloadBL.mean_diff.plot(contrast_bars=False, show_baseline_ec=True,
                                   raw_ylim=(-0.5, 10.5), color_col='UASLAeq',
                                       delta_dot=False,
                                       slopegraph_kwargs={'linewidth': 0.25, 'alpha': 0.45, 'jitter': 0.3},
                                       contrast_marker_size=4, contrast_ylim=(-1, 1.5), float_contrast=True,
                                       delta_text=True,
                                       delta_text_kwargs={'y_coordinates': [0.5]*(len(data['UASOperation'].unique()) - 1)},
                                       legend_kwargs={'loc': 'upper center', 'title': r"UAS $L_\mathrm{Aeq}$", 'fontsize': 12,
                                                      'frameon': True, 'title_fontsize': 12}, ax=ax)
diffs = dataHiloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=0.0, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.74 + ii + 0.09, y=-0.2, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

lgd = ax.get_legend()
for lh in lgd.legend_handles:
    lh.set_alpha(0.8)
    lh.set_linewidth(3)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtAOpHiDabest.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtAOpHiDabest.pdf"),
                format='pdf', bbox_inches='tight')

dataHiMD;
dataHiloadBL.mean_diff.statistical_tests


# %%
out = pd.concat([out, dataHiloadBL.mean_diff.statistical_tests])
# calculate effect sizes - Cohen's d
outd = pd.concat([outd, dataHiloadBL.cohens_d.statistical_tests])
dataHiloadBL.cohens_d.statistical_tests


# %% [markdown]
# ### UAS type and operation combined

# %%
# select subset of data for analysis and sort
data = dataBySubjTestA[~dataBySubjTestA['StimFile'].str.contains("PwrScale")]
data = data[data['UASType'] != "Baseline"]
data = data.loc[:, ['ID', 'StimFile', 'AmbientEnv', 'UASLAeq', 'UASType', 'UASOperation', 'Annoyance']]
data.sort_values(by=['AmbientEnv', 'ID', 'UASLAeq', 'UASType', 'UASOperation'], inplace=True)

# take highest LAeq in park env
dataPk = data.loc[(data['AmbientEnv'] == "Park"), :]
dataPk['UASTypeOp'] = dataPk['UASType'].astype(str) + " " + dataPk['UASOperation'].astype(str)
dataPk.sort_values(by=['ID', 'UASLAeq', 'UASTypeOp'], inplace=True)

# take highest LAeq in park env
dataSt = data.loc[(data['AmbientEnv'] == "Street"), :]
dataSt['UASTypeOp'] = dataSt['UASType'].astype(str) + " " + dataSt['UASOperation'].astype(str)
dataSt.sort_values(by=['ID', 'UASLAeq', 'UASTypeOp'], inplace=True)

# create a dummy pairing ID
data['UASTypeOp'] = data['UASType'].astype(str) + " " + data['UASOperation'].astype(str)
data.sort_values(by=['AmbientEnv', 'ID', 'UASLAeq', 'UASTypeOp'], inplace=True)

# drop the "Baseline" category by converting to a plain string column
data['UASLAeq'] = data['UASLAeq'].astype(str)
dataPk['UASLAeq'] = dataPk['UASLAeq'].astype(str)
dataSt['UASLAeq'] = dataSt['UASLAeq'].astype(str)

grouping_cols_all = ['ID', 'AmbientEnv', 'UASLAeq']
grouping_cols_env = ['ID', 'UASLAeq']

# ngroup() automatically assigns an identical integer to matching sets
# across the different UAS type/operation conditions
data['dummyID'] = data.groupby(grouping_cols_all).ngroup() + 1
dataPk['dummyID'] = dataPk.groupby(grouping_cols_env).ngroup() + 1
dataSt['dummyID'] = dataSt.groupby(grouping_cols_env).ngroup() + 1
data.sort_values(by=['dummyID', 'ID', 'StimFile'], inplace=True)
dataPk.sort_values(by=['dummyID', 'ID', 'StimFile'], inplace=True)
dataSt.sort_values(by=['dummyID', 'ID', 'StimFile'], inplace=True)


# %%
# assign data for processing

dataloadBL = dabest.load(ps_adjust=True, data=data,  idx=("T150 Overflight", "M300 Overflight", "H520 Overflight",
                                          "T150 Landing", "M300 Landing", "H520 Landing",
                                          "T150 Takeoff", "M300 Takeoff", "H520 Takeoff"), x='UASTypeOp', y='Annoyance',
                         paired='baseline', id_col='dummyID', resamples=5000, random_seed=42828412)

dataPkloadBL = dabest.load(ps_adjust=True, data=dataPk,  idx=("T150 Overflight", "M300 Overflight", "H520 Overflight",
                                              "T150 Landing", "M300 Landing", "H520 Landing",
                                              "T150 Takeoff", "M300 Takeoff", "H520 Takeoff"), x='UASTypeOp', y='Annoyance',
                           paired='baseline', id_col='dummyID', resamples=5000, random_seed=82108)
dataStloadBL = dabest.load(ps_adjust=True, data=dataSt,  idx=("T150 Overflight", "M300 Overflight", "H520 Overflight",
                                              "T150 Landing", "M300 Landing", "H520 Landing",
                                              "T150 Takeoff", "M300 Takeoff", "H520 Takeoff"), x='UASTypeOp', y='Annoyance',
                           paired='baseline', id_col='dummyID', resamples=5000, random_seed=1198232)

# %%
##### calculate effect sizes and plot
fig, ax = plt.subplots(figsize=(12, 3.5))
dataMD = dataloadBL.mean_diff.plot(contrast_bars=False, show_baseline_ec=True,
                                   raw_ylim=(-0.5, 10.5), color_col='UASLAeq',
                                   delta_dot=False,
                                   slopegraph_kwargs={'linewidth': 0.25, 'alpha': 0.45, 'jitter': 0.3},
                                   contrast_marker_size=4, contrast_ylim=(-1, 2), float_contrast=True,
                                   delta_text=True,
                                   delta_text_kwargs={'y_coordinates': [0.8]*(len(data['UASTypeOp'].unique()) - 1)},
                                   legend_kwargs={'loc': 'upper center', 'title': r"UAS $L_\mathrm{Aeq}$", 'fontsize': 12,
                                                  'frameon': True, 'title_fontsize': 12}, ax=ax)
diffs = dataloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=0.2, s="95%: [" + " ".join(vals) + "]", fontsize=9.5)
    contrast_ax.text(x=0.74 + ii + 0.09, y=-0.04, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=9.5)

lgd = ax.get_legend()
for lh in lgd.legend_handles:
    lh.set_alpha(0.8)
    lh.set_linewidth(3)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtATypeOpDabest.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtATypeOpDabest.pdf"),
                format='pdf', bbox_inches='tight')

dataMD;
dataloadBL.mean_diff.statistical_tests


# %%
out = pd.concat([out, dataloadBL.mean_diff.statistical_tests])
# calculate effect sizes - Cohen's d
outd = pd.concat([outd, dataloadBL.cohens_d.statistical_tests])
dataloadBL.cohens_d.statistical_tests


# %%
# calculate effect sizes and plot
fig, ax = plt.subplots(figsize=(12, 3.5))
dataPkMD = dataPkloadBL.mean_diff.plot(contrast_bars=False, show_baseline_ec=True,
                                   raw_ylim=(-0.5, 10.5), color_col='UASLAeq',
                                       delta_dot=False,
                                       slopegraph_kwargs={'linewidth': 0.25, 'alpha': 0.45, 'jitter': 0.3},
                                       contrast_marker_size=4, contrast_ylim=(-1, 2), float_contrast=True,
                                       delta_text=True,
                                       delta_text_kwargs={'y_coordinates': [0.8]*(len(data['UASTypeOp'].unique()) - 1)},
                                       legend_kwargs={'loc': 'upper center', 'title': r"UAS $L_\mathrm{Aeq}$", 'fontsize': 12,
                                                      'frameon': True, 'title_fontsize': 12}, ax=ax)
diffs = dataPkloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=0.2, s="95%: [" + " ".join(vals) + "]", fontsize=9.5)
    contrast_ax.text(x=0.74 + ii + 0.09, y=-0.04, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=9.5)

lgd = ax.get_legend()
for lh in lgd.legend_handles:
    lh.set_alpha(0.8)
    lh.set_linewidth(3)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtAPkTypeOpDabest.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtAPkTypeOpDabest.pdf"),
                format='pdf', bbox_inches='tight')

dataPkMD;
dataPkloadBL.mean_diff.statistical_tests


# %%
out = pd.concat([out, dataPkloadBL.mean_diff.statistical_tests])
# calculate effect sizes - Cohen's d
outd = pd.concat([outd, dataPkloadBL.cohens_d.statistical_tests])
dataPkloadBL.cohens_d.statistical_tests


# %%
# calculate effect sizes and plot
fig, ax = plt.subplots(figsize=(12, 3.5))
dataStMD = dataStloadBL.mean_diff.plot(contrast_bars=False, show_baseline_ec=True,
                                   raw_ylim=(-0.5, 10.5), color_col='UASLAeq',
                                       delta_dot=False,
                                       slopegraph_kwargs={'linewidth': 0.25, 'alpha': 0.45, 'jitter': 0.3},
                                       contrast_marker_size=4, contrast_ylim=(-1, 2), float_contrast=True,
                                       delta_text=True,
                                       delta_text_kwargs={'y_coordinates': [0.8]*(len(data['UASTypeOp'].unique()) - 1)},
                                       legend_kwargs={'loc': 'upper center', 'title': r"UAS $L_\mathrm{Aeq}$", 'fontsize': 12,
                                                      'frameon': True, 'title_fontsize': 12}, ax=ax)
diffs = dataStloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=0.2, s="95%: [" + " ".join(vals) + "]", fontsize=9.5)
    contrast_ax.text(x=0.74 + ii + 0.09, y=-0.04, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=9.5)

lgd = ax.get_legend()
for lh in lgd.legend_handles:
    lh.set_alpha(0.8)
    lh.set_linewidth(3)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtAStTypeOpDabest.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtAStTypeOpDabest.pdf"),
                format='pdf', bbox_inches='tight')

dataStMD;
dataStloadBL.mean_diff.statistical_tests


# %%
out = pd.concat([out, dataStloadBL.mean_diff.statistical_tests])
# calculate effect sizes - Cohen's d
outd = pd.concat([outd, dataStloadBL.cohens_d.statistical_tests])
dataStloadBL.cohens_d.statistical_tests


# %% [markdown]
# ### AAM attitude

# %%
# select subset of data for analysis and sort
data = dataBySubjTestA
data = data[data['UASLAeq'] != "Baseline"]
data = data.loc[:, ['ID', 'StimFile', 'UASLAeq', 'AAM_attitude', 'Annoyance']]


# %%
# assign data for processing
dataloadBL = dabest.load(ps_adjust=True, data=data, idx=['Supportive', 'Ambivalent', 'Concerned', 'Neutral'],
                         x='AAM_attitude', y='Annoyance', resamples=5000, random_seed=3616)


# %%
# calculate baseline paired effect sizes and plot
fig, ax = plt.subplots(figsize=(6, 3))
dataMD = dataloadBL.mean_diff.plot(contrast_bars=False,
                                   raw_ylim=(0, 10),
                                   raw_marker_size=0.005,
                                   swarmplot_kwargs={'alpha': 0.3},
                                   delta_dot_kwargs={'size': 1, 'side': 'left', 'alpha': 0.25, 'zorder': 1},
                                   delta_text=True,
                                   custom_palette={'Neutral': mycolours[0],
                                                   'Ambivalent': mycolours[1],
                                                   'Concerned': mycolours[2],
                                                   'Supportive': mycolours[3]},
                                   raw_desat=1,
                                   contrast_desat=1,
                                   contrast_marker_size=4,
                                   contrast_ylim=(-1, 3),
                                   ax=ax,
                                   )

diffs = dataloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=-3.9, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.7 + ii + 0.09, y=-4.28, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtAAttDabestBase.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtAAttDabestBase.pdf"),
                format='pdf', bbox_inches='tight')
out = pd.concat([out, dataloadBL.mean_diff.statistical_tests])

dataMD;
dataloadBL.mean_diff.statistical_tests


# %%
# calculate and display Cohen's d
outd = pd.concat([outd, dataloadBL.cohens_d.statistical_tests])
dataloadBL.cohens_d.statistical_tests


# %% [markdown]
# ### Home residence area

# %%
# select subset of data for analysis and sort
data = dataBySubjTestA
data = data[data['UASLAeq'] != "Baseline"]
data = data.loc[:, ['ID', 'StimFile', 'UASLAeq', 'Home_Area', 'Annoyance']]


# %%
# assign data for processing
dataloadBL = dabest.load(ps_adjust=True, data=data, idx=['Urban', 'Suburban', 'Rural'],
                         x='Home_Area', y='Annoyance', resamples=5000, random_seed=80842)


# %%
# calculate baseline paired effect sizes and plot
fig, ax = plt.subplots(figsize=(6, 3))
dataMD = dataloadBL.mean_diff.plot(contrast_bars=False,
                                   raw_ylim=(0, 10),
                                   raw_marker_size=0.005,
                                   swarmplot_kwargs={'alpha': 0.3},
                                   delta_dot_kwargs={'size': 1, 'side': 'left', 'alpha': 0.25, 'zorder': 1},
                                   delta_text=True,
                                   custom_palette={'Suburban': mycolours[0],
                                                   'Rural': mycolours[1],
                                                   'Urban': mycolours[3]},
                                   raw_desat=1,
                                   contrast_desat=1,
                                   contrast_marker_size=4,
                                   contrast_ylim=(-2, 2),
                                   ax=ax,
                                   )

diffs = dataloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=-3.9, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.7 + ii + 0.09, y=-4.28, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtAAORDabestBase.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtAAORDabestBase.pdf"),
                format='pdf', bbox_inches='tight')
out = pd.concat([out, dataloadBL.mean_diff.statistical_tests])

dataMD;
dataloadBL.mean_diff.statistical_tests


# %%
# calculate and display Cohen's d
outd = pd.concat([outd, dataloadBL.cohens_d.statistical_tests])
dataloadBL.cohens_d.statistical_tests


# %% [markdown]
# ## Part B
# 
# ### Number of events

# %%
# select subset of data for analysis and sort
data = dataBySubjTestB
data = data[data['UASLAeq'] != "Baseline"]
data['UASEvents'] = data['UASEvents'].astype(int).astype(str)
data = data.loc[:, ['ID', 'StimFile', 'UASEvents', 'UASLAeq', 'UASType', 'Annoyance']]
data.sort_values(by=['ID', 'UASType', 'UASLAeq', 'UASEvents'], inplace=True)

dataLo = data.loc[data['UASLAeq'] == "54", :]
dataHi = data.loc[data['UASLAeq'] == "60", :]
dataH520 = data.loc[data['UASType'] == "H520", :]
dataT150 = data.loc[data['UASType'] == "T150", :]
dataLoH520 = data.loc[(data['UASType'] == "H520") & (data['UASLAeq'] == "54"), :]
dataLoT150 = data.loc[(data['UASType'] == "T150") & (data['UASLAeq'] == "54"), :]
dataHiH520 = data.loc[(data['UASType'] == "H520") & (data['UASLAeq'] == "60"), :]
dataHiT150 = data.loc[(data['UASType'] == "T150") & (data['UASLAeq'] == "60"), :]

# create dummy pairing ID for each participant to allow for paired analysis.
# Grouping by ID, UASType and UASLAeq works for every subset below: any
# column that's already constant within a given subset (e.g. UASType in
# dataH520) simply doesn't split the groups any further, so one grouping
# column list can safely be reused across all of them.
grouping_cols = ['ID', 'UASType', 'UASLAeq']

# ngroup() automatically assigns an identical integer to matching sets
# across the different UAS-events conditions
for dataset in [data, dataLo, dataHi, dataH520, dataT150, dataLoH520, dataLoT150, dataHiH520, dataHiT150]:
    # drop the "Baseline"/other-level categories by converting to plain string columns
    dataset['UASLAeq'] = dataset['UASLAeq'].astype(str)
    dataset['UASType'] = dataset['UASType'].astype(str)
    dataset['dummyID'] = dataset.groupby(grouping_cols).ngroup() + 1
    dataset.sort_values(by=['dummyID', 'ID', 'StimFile'], inplace=True)


# %%
# assign data for processing
dataloadBL = dabest.load(ps_adjust=True, data=data, idx=data['UASEvents'].unique(),
                         x='UASEvents', y='Annoyance', paired='baseline',
                         id_col='dummyID', resamples=5000, random_seed=888)
dataloadSQ = dabest.load(ps_adjust=True, data=data, idx=data['UASEvents'].unique(),
                         x='UASEvents', y='Annoyance', paired='sequential',
                         id_col='dummyID', resamples=5000, random_seed=456)

dataLoloadBL = dabest.load(ps_adjust=True, data=dataLo, idx=dataLo['UASEvents'].unique(),
                           x='UASEvents', y='Annoyance', paired='baseline',
                           id_col='dummyID', resamples=5000, random_seed=4373)

dataLoloadSQ = dabest.load(ps_adjust=True, data=dataLo, idx=dataLo['UASEvents'].unique(),
                           x='UASEvents', y='Annoyance', paired='sequential',
                           id_col='dummyID', resamples=5000, random_seed=714)

dataHiloadSQ = dabest.load(ps_adjust=True, data=dataHi, idx=dataHi['UASEvents'].unique(),
                           x='UASEvents', y='Annoyance', paired='sequential',
                           id_col='dummyID', resamples=5000, random_seed=99595)

dataH520loadSQ = dabest.load(ps_adjust=True, data=dataH520, idx=dataH520['UASEvents'].unique(),
                               x='UASEvents', y='Annoyance', paired='sequential',
                               id_col='dummyID', resamples=5000, random_seed=8484984)

dataT150loadSQ = dabest.load(ps_adjust=True, data=dataT150, idx=dataT150['UASEvents'].unique(),
                               x='UASEvents', y='Annoyance', paired='sequential',
                               id_col='dummyID', resamples=5000, random_seed=2313)

dataLoH520loadSQ = dabest.load(ps_adjust=True, data=dataLoH520, idx=dataLoH520['UASEvents'].unique(),
                               x='UASEvents', y='Annoyance', paired='sequential',
                               id_col='dummyID', resamples=5000, random_seed=15619)

dataHiH520loadSQ = dabest.load(ps_adjust=True, data=dataHiH520, idx=dataHiH520['UASEvents'].unique(),
                               x='UASEvents', y='Annoyance', paired='sequential',
                               id_col='dummyID', resamples=5000, random_seed=503)

dataLoT150loadSQ = dabest.load(ps_adjust=True, data=dataLoT150, idx=dataLoT150['UASEvents'].unique(),
                               x='UASEvents', y='Annoyance', paired='sequential',
                               id_col='dummyID', resamples=5000, random_seed=4025498)

dataHiT150loadSQ = dabest.load(ps_adjust=True, data=dataHiT150, idx=dataHiT150['UASEvents'].unique(),
                               x='UASEvents', y='Annoyance', paired='sequential',
                               id_col='dummyID', resamples=5000, random_seed=117)

# %%
# calculate baseline paired effect sizes and plot
fig, ax = plt.subplots(figsize=(6, 3))
dataMD = dataloadBL.mean_diff.plot(contrast_bars=False, show_baseline_ec=True,
                                   raw_ylim=(-0.5, 10.5), color_col='UASLAeq',
                                   delta_dot=False,
                                   slopegraph_kwargs={'linewidth': 0.25, 'alpha': 0.45, 'jitter': 0.3},
                                   delta_text=True,
                                   delta_text_kwargs={'y_coordinates': [1.5]*(len(data['UASEvents'].unique()) - 1)},
                                   contrast_marker_size=4, contrast_ylim=(0, 2.5),
                                   legend_kwargs={'loc': 'upper center', 'title': r"UAS $L_\mathrm{Aeq}$", 'fontsize': 12,
                                                  'frameon': True, 'title_fontsize': 12}, ax=ax)

diffs = dataloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=1.0, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.74 + ii + 0.09, y=0.8, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

lgd = ax.get_legend()
for lh in lgd.legend_handles:
    lh.set_alpha(0.8)
    lh.set_linewidth(3)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtBDabestBase.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtBDabestBase.pdf"),
                format='pdf', bbox_inches='tight')

dataMD;
dataloadBL.mean_diff.statistical_tests


# %%
out = pd.concat([out, dataloadBL.mean_diff.statistical_tests])
# calculate and display Cohen's d
outd = pd.concat([outd, dataloadBL.cohens_d.statistical_tests])
dataloadBL.cohens_d.statistical_tests


# %%
# calculate sequential paired effect sizes and plot
fig, ax = plt.subplots(figsize=(6, 3))
dataMD = dataloadSQ.mean_diff.plot(contrast_bars=True,
                                   raw_ylim=(0, 10), color_col='UASLAeq',
                                   delta_dot=False,
                                   slopegraph_kwargs={'linewidth': 0.25, 'alpha': 0.45, 'jitter': 0.3},
                                   delta_text=True,
                                   contrast_marker_size=4, contrast_ylim=(-0.5, 2),
                                   legend_kwargs={'loc': 'upper center', 'title': r"UAS $L_\mathrm{Aeq}$", 'fontsize': 12,
                                                  'frameon': True, 'title_fontsize': 12}, ax=ax)

diffs = dataloadSQ.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=0.5, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.74 + ii + 0.09, y=0.3, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

lgd = ax.get_legend()
for lh in lgd.legend_handles:
    lh.set_alpha(0.8)
    lh.set_linewidth(3)
    
if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtBDabestSeq.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtBDabestSeq.pdf"),
                format='pdf', bbox_inches='tight')

dataMD;
dataloadSQ.mean_diff.statistical_tests


# %%
out = pd.concat([out, dataloadSQ.mean_diff.statistical_tests])
# calculate and display Cohen's d
outd = pd.concat([outd, dataloadSQ.cohens_d.statistical_tests])
dataloadSQ.cohens_d.statistical_tests


# %% [markdown]
# #### Segregated by LAeq
# 
# ##### 54 dB

# %%
# calculate sequential paired effect sizes and plot
fig, ax = plt.subplots(figsize=(6, 3))
dataLoMD = dataLoloadSQ.mean_diff.plot(contrast_bars=True,
                                   raw_ylim=(-0.5, 10.5), color_col='UASLAeq',
                                       delta_dot=False,
                                       slopegraph_kwargs={'linewidth': 0.25, 'alpha': 0.45, 'jitter': 0.3},
                                       delta_text=True,
                                       contrast_marker_size=4, contrast_ylim=(-1, 2),
                                       legend_kwargs={'loc': 'upper center', 'title': r"UAS $L_\mathrm{Aeq}$", 'fontsize': 12,
                                                      'frameon': True, 'title_fontsize': 12}, ax=ax)

diffs = dataLoloadSQ.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=0.2, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.74 + ii + 0.09, y=-0.04, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

lgd = ax.get_legend()
for lh in lgd.legend_handles:
    lh.set_alpha(0.8)
    lh.set_linewidth(3)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtB54DabestSeq.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtB54DabestSeq.pdf"),
                format='pdf', bbox_inches='tight')

dataLoMD;
dataLoloadSQ.mean_diff.statistical_tests


# %%
out = pd.concat([out, dataLoloadSQ.mean_diff.statistical_tests])
# calculate and display Cohen's d
outd = pd.concat([outd, dataLoloadSQ.cohens_d.statistical_tests])
dataLoloadSQ.cohens_d.statistical_tests


# %% [markdown]
# ##### 60 dB

# %%
# calculate sequential paired effect sizes and plot
fig, ax = plt.subplots(figsize=(6, 3))
dataHiMD = dataHiloadSQ.mean_diff.plot(contrast_bars=True,
                                   raw_ylim=(-0.5, 10.5), color_col='UASLAeq',
                                       delta_dot=False,
                                       slopegraph_kwargs={'linewidth': 0.25, 'alpha': 0.45, 'jitter': 0.3},
                                       delta_text=True,
                                       contrast_marker_size=4, contrast_ylim=(-0.5, 2),
                                       legend_kwargs={'loc': 'upper center', 'title': r"UAS $L_\mathrm{Aeq}$", 'fontsize': 12,
                                                      'frameon': True, 'title_fontsize': 12}, ax=ax)


diffs = dataHiloadSQ.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=0.5, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.74 + ii + 0.09, y=0.3, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

lgd = ax.get_legend()
for lh in lgd.legend_handles:
    lh.set_alpha(0.8)
    lh.set_linewidth(3)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtB60DabestSeq.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtB60DabestSeq.pdf"),
                format='pdf', bbox_inches='tight')

dataHiMD;
dataHiloadSQ.mean_diff.statistical_tests


# %%
out = pd.concat([out, dataHiloadSQ.mean_diff.statistical_tests])
# calculate and display Cohen's d
outd = pd.concat([outd, dataHiloadSQ.cohens_d.statistical_tests])
dataHiloadSQ.cohens_d.statistical_tests


# %% [markdown]
# #### Segregated by UAS LAeq and type
# 
# ##### 54 dB, H520

# %%
# calculate sequential paired effect sizes and plot
fig, ax = plt.subplots(figsize=(6, 3))
dataLoH520MD = dataLoH520loadSQ.mean_diff.plot(contrast_bars=True,
                                   raw_ylim=(-0.5, 10.5), color_col='UASLAeq',
                                               delta_dot=False,
                                               slopegraph_kwargs={'linewidth': 0.25, 'alpha': 0.45, 'jitter': 0.3},
                                               delta_text=True,
                                               contrast_marker_size=4, contrast_ylim=(-1.5, 2.5),
                                               legend_kwargs={'loc': 'upper center', 'title': r"H520 $L_\mathrm{Aeq}$", 'fontsize': 12,
                                                              'frameon': True, 'title_fontsize': 12}, ax=ax)

diffs = dataLoH520loadSQ.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=0.1, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.74 + ii + 0.09, y=-0.22, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

lgd = ax.get_legend()
for lh in lgd.legend_handles:
    lh.set_alpha(0.8)
    lh.set_linewidth(3)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtB54H520DabestSeq.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtB54H520DabestSeq.pdf"),
                format='pdf', bbox_inches='tight')

dataLoH520MD;
dataLoH520loadSQ.mean_diff.statistical_tests


# %%
out = pd.concat([out, dataLoH520loadSQ.mean_diff.statistical_tests])
# calculate and display Cohen's d
outd = pd.concat([outd, dataLoH520loadSQ.cohens_d.statistical_tests])
dataLoH520loadSQ.cohens_d.statistical_tests


# %% [markdown]
# ##### 60 dB, H520

# %%
# calculate sequential paired effect sizes and plot
fig, ax = plt.subplots(figsize=(6, 3))
dataHiH520MD = dataHiH520loadSQ.mean_diff.plot(contrast_bars=True,
                                   raw_ylim=(-0.5, 10.5), color_col='UASLAeq',
                                               delta_dot=False,
                                               slopegraph_kwargs={'linewidth': 0.25, 'alpha': 0.45, 'jitter': 0.3},
                                               delta_text=True,
                                               contrast_marker_size=4, contrast_ylim=(-1.5, 2.5),
                                               legend_kwargs={'loc': 'upper center', 'title': r"H520 $L_\mathrm{Aeq}$", 'fontsize': 12,
                                                              'frameon': True, 'title_fontsize': 12}, ax=ax)

diffs = dataHiH520loadSQ.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=0.1, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.74 + ii + 0.09, y=-0.22, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

lgd = ax.get_legend()
for lh in lgd.legend_handles:
    lh.set_alpha(0.8)
    lh.set_linewidth(3)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtB60H520DabestSeq.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtB60H520DabestSeq.pdf"),
                format='pdf', bbox_inches='tight')

dataHiH520MD;
dataHiH520loadSQ.mean_diff.statistical_tests


# %%
out = pd.concat([out, dataHiH520loadSQ.mean_diff.statistical_tests])
# calculate and display Cohen's d
outd = pd.concat([outd, dataHiH520loadSQ.cohens_d.statistical_tests])
dataHiH520loadSQ.cohens_d.statistical_tests


# %% [markdown]
# ##### 54 dB, T150

# %%
# calculate sequential paired effect sizes and plot
fig, ax = plt.subplots(figsize=(6, 3))
dataLoT150MD = dataLoT150loadSQ.mean_diff.plot(contrast_bars=True,
                                   raw_ylim=(-0.5, 10.5), color_col='UASLAeq',
                                               delta_dot=False,
                                               slopegraph_kwargs={'linewidth': 0.25, 'alpha': 0.45, 'jitter': 0.3},
                                               delta_text=True,
                                               contrast_marker_size=4, contrast_ylim=(-1.5, 2.5),
                                               legend_kwargs={'loc': 'upper center', 'title': r"T150 $L_\mathrm{Aeq}$", 'fontsize': 12,
                                                              'frameon': True, 'title_fontsize': 12}, ax=ax)

diffs = dataLoT150loadSQ.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=0.1, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.74 + ii + 0.09, y=-0.22, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

lgd = ax.get_legend()
for lh in lgd.legend_handles:
    lh.set_alpha(0.8)
    lh.set_linewidth(3)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtB54T150DabestSeq.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtB54T150DabestSeq.pdf"),
                format='pdf', bbox_inches='tight')

dataLoT150MD;
dataLoT150loadSQ.mean_diff.statistical_tests


# %%
out = pd.concat([out, dataLoT150loadSQ.mean_diff.statistical_tests])
# calculate and display Cohen's d
outd = pd.concat([outd, dataLoT150loadSQ.cohens_d.statistical_tests])
dataLoT150loadSQ.cohens_d.statistical_tests


# %% [markdown]
# ##### 60 dB, T150

# %%
# calculate sequential paired effect sizes and plot
fig, ax = plt.subplots(figsize=(6, 3))
dataHiT150MD = dataHiT150loadSQ.mean_diff.plot(contrast_bars=True,
                                   raw_ylim=(-0.5, 10.5), color_col='UASLAeq',
                                               delta_dot=False,
                                               slopegraph_kwargs={'linewidth': 0.25, 'alpha': 0.45, 'jitter': 0.3},
                                               delta_text=True,
                                               contrast_marker_size=4, contrast_ylim=(-1.5, 2.5),
                                               legend_kwargs={'loc': 'upper center', 'title': r"T150 $L_\mathrm{Aeq}$", 'fontsize': 12,
                                                              'frameon': True, 'title_fontsize': 12}, ax=ax)


diffs = dataHiT150loadSQ.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=0.1, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.74 + ii + 0.09, y=-0.22, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

lgd = ax.get_legend()
for lh in lgd.legend_handles:
    lh.set_alpha(0.8)
    lh.set_linewidth(3)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtB60T150DabestSeq.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtB60T150DabestSeq.pdf"),
                format='pdf', bbox_inches='tight')

dataHiT150MD;
dataHiT150loadSQ.mean_diff.statistical_tests


# %%
out = pd.concat([out, dataHiT150loadSQ.mean_diff.statistical_tests])
# calculate and display Cohen's d
outd = pd.concat([outd, dataHiT150loadSQ.cohens_d.statistical_tests])
dataHiT150loadSQ.cohens_d.statistical_tests


# %% [markdown]
# #### Segregated by UAS type
# 
# ##### H520

# %%
# calculate sequential paired effect sizes and plot
fig, ax = plt.subplots(figsize=(6, 3))
dataH520MD = dataH520loadSQ.mean_diff.plot(contrast_bars=True,
                                   raw_ylim=(-0.5, 10.5), color_col='UASLAeq',
                                               delta_dot=False,
                                               slopegraph_kwargs={'linewidth': 0.25, 'alpha': 0.45, 'jitter': 0.3},
                                               delta_text=True,
                                               contrast_marker_size=4, contrast_ylim=(-1.5, 2.5),
                                               legend_kwargs={'loc': 'upper center', 'title': r"H520 $L_\mathrm{Aeq}$", 'fontsize': 12,
                                                              'frameon': True, 'title_fontsize': 12}, ax=ax)

diffs = dataH520loadSQ.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=0.1, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.74 + ii + 0.09, y=-0.22, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

lgd = ax.get_legend()
for lh in lgd.legend_handles:
    lh.set_alpha(0.8)
    lh.set_linewidth(3)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtBH520DabestSeq.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtBH520DabestSeq.pdf"),
                format='pdf', bbox_inches='tight')

dataH520MD;
dataH520loadSQ.mean_diff.statistical_tests


# %%
out = pd.concat([out, dataH520loadSQ.mean_diff.statistical_tests])
# calculate and display Cohen's d
outd = pd.concat([outd, dataH520loadSQ.cohens_d.statistical_tests])
dataH520loadSQ.cohens_d.statistical_tests


# %% [markdown]
# ##### T150

# %%
# calculate sequential paired effect sizes and plot
fig, ax = plt.subplots(figsize=(6, 3))
dataT150MD = dataT150loadSQ.mean_diff.plot(contrast_bars=True,
                                   raw_ylim=(-0.5, 10.5), color_col='UASLAeq',
                                               delta_dot=False,
                                               slopegraph_kwargs={'linewidth': 0.25, 'alpha': 0.45, 'jitter': 0.3},
                                               delta_text=True,
                                               contrast_marker_size=4, contrast_ylim=(-1.5, 2.5),
                                               legend_kwargs={'loc': 'upper center', 'title': r"T150 $L_\mathrm{Aeq}$", 'fontsize': 12,
                                                              'frameon': True, 'title_fontsize': 12}, ax=ax)

diffs = dataT150loadSQ.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=0.1, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.74 + ii + 0.09, y=-0.22, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

lgd = ax.get_legend()
for lh in lgd.legend_handles:
    lh.set_alpha(0.8)
    lh.set_linewidth(3)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtBT150DabestSeq.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtBT150DabestSeq.pdf"),
                format='pdf', bbox_inches='tight')

dataT150MD;
dataT150loadSQ.mean_diff.statistical_tests


# %%
out = pd.concat([out, dataT150loadSQ.mean_diff.statistical_tests])
# calculate and display Cohen's d
outd = pd.concat([outd, dataT150loadSQ.cohens_d.statistical_tests])
dataT150loadSQ.cohens_d.statistical_tests


# %% [markdown]
# ### AAM attitude

# %%
# select subset of data for analysis and sort
data = dataBySubjTestB
data = data[data['UASLAeq'] != "Baseline"]
data = data.loc[:, ['ID', 'StimFile', 'UASLAeq', 'AAM_attitude', 'Annoyance']]

# separate by level
dataLo = data.loc[data['UASLAeq'] == "54", :]
dataHi = data.loc[data['UASLAeq'] == "60", :]


# %%
# assign data for processing
dataloadBL = dabest.load(ps_adjust=True, data=data, idx=['Supportive', 'Ambivalent', 'Concerned', 'Neutral'],
                         x='AAM_attitude', y='Annoyance', resamples=5000, random_seed=24624)

dataLoloadBL = dabest.load(ps_adjust=True, data=dataLo, idx=['Supportive', 'Ambivalent', 'Concerned', 'Neutral'],
                           x='AAM_attitude', y='Annoyance', resamples=5000, random_seed=8478)

dataHiloadBL = dabest.load(ps_adjust=True, data=dataHi, idx=['Supportive', 'Ambivalent', 'Concerned', 'Neutral'],
                           x='AAM_attitude', y='Annoyance', resamples=5000, random_seed=13)


# %%
# calculate baseline paired effect sizes and plot
fig, ax = plt.subplots(figsize=(6, 3))
dataMD = dataloadBL.mean_diff.plot(contrast_bars=False,
                                   raw_ylim=(0, 10),
                                   raw_marker_size=0.1,
                                   swarmplot_kwargs={'alpha': 0.3},
                                   delta_dot_kwargs={'size': 1, 'side': 'left', 'alpha': 0.25, 'zorder': 1},
                                   delta_text=True,
                                   custom_palette={'Neutral': mycolours[0],
                                                   'Ambivalent': mycolours[1],
                                                   'Concerned': mycolours[2],
                                                   'Supportive': mycolours[3]},
                                   raw_desat=1,
                                   contrast_desat=1,
                                   contrast_marker_size=4,
                                   contrast_ylim=(-1, 4),
                                   ax=ax,
                                   )

diffs = dataloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=-3.9, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.7 + ii + 0.09, y=-4.375, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtBAttDabestBase.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtBAttDabestBase.pdf"),
                format='pdf', bbox_inches='tight')
out = pd.concat([out, dataloadBL.mean_diff.statistical_tests])

dataMD;
dataloadBL.mean_diff.statistical_tests


# %%
# calculate and display Cohen's d
outd = pd.concat([outd, dataloadBL.cohens_d.statistical_tests])
dataloadBL.cohens_d.statistical_tests


# %%
# calculate baseline paired effect sizes and plot
fig, ax = plt.subplots(figsize=(6, 3))
dataLoMD = dataLoloadBL.mean_diff.plot(contrast_bars=False,
                                   raw_ylim=(0, 10),
                                       raw_marker_size=0.1,
                                       swarmplot_kwargs={'alpha': 0.3},
                                       delta_dot_kwargs={'size': 1, 'side': 'left', 'alpha': 0.25, 'zorder': 1},
                                       delta_text=True,
                                       custom_palette={'Neutral': mycolours[0],
                                                       'Ambivalent': mycolours[1],
                                                       'Concerned': mycolours[2],
                                                       'Supportive': mycolours[3]},
                                       raw_desat=1,
                                       contrast_desat=1,
                                       contrast_marker_size=4,
                                       contrast_ylim=(-1, 4),
                                       ax=ax,
                                       )

diffs = dataLoloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=-3.9, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.7 + ii + 0.09, y=-4.375, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtB54AttDabestBase.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtB54AttDabestBase.pdf"),
                format='pdf', bbox_inches='tight')
out = pd.concat([out, dataLoloadBL.mean_diff.statistical_tests])

dataLoMD;
dataLoloadBL.mean_diff.statistical_tests


# %%
# calculate and display Cohen's d
outd = pd.concat([outd, dataLoloadBL.cohens_d.statistical_tests])
dataLoloadBL.cohens_d.statistical_tests


# %%
# calculate baseline paired effect sizes and plot
fig, ax = plt.subplots(figsize=(6, 3))
dataHiMD = dataHiloadBL.mean_diff.plot(contrast_bars=False,
                                   raw_ylim=(0, 10),
                                       raw_marker_size=0.1,
                                       swarmplot_kwargs={'alpha': 0.3},
                                       delta_dot_kwargs={'size': 1, 'side': 'left', 'alpha': 0.25, 'zorder': 1},
                                       delta_text=True,
                                       custom_palette={'Neutral': mycolours[0],
                                                       'Ambivalent': mycolours[1],
                                                       'Concerned': mycolours[2],
                                                       'Supportive': mycolours[3]},
                                       raw_desat=1,
                                       contrast_desat=1,
                                       contrast_marker_size=4,
                                       contrast_ylim=(-1, 4),
                                       ax=ax,
                                       )

diffs = dataHiloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=-3.9, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.7 + ii + 0.09, y=-4.375, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtB60AttDabestBase.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtB60AttDabestBase.pdf"),
                format='pdf', bbox_inches='tight')
out = pd.concat([out, dataHiloadBL.mean_diff.statistical_tests])

dataHiMD;
dataHiloadBL.mean_diff.statistical_tests


# %%
# calculate and display Cohen's d
outd = pd.concat([outd, dataHiloadBL.cohens_d.statistical_tests])
dataHiloadBL.cohens_d.statistical_tests


# %% [markdown]
# ### Home residence area

# %%
# select subset of data for analysis and sort
data = dataBySubjTestB
data = data[data['UASLAeq'] != "Baseline"]
data = data.loc[:, ['ID', 'StimFile', 'UASLAeq', 'Home_Area', 'Annoyance']]


# %%
# assign data for processing
dataloadBL = dabest.load(ps_adjust=True, data=data, idx=['Urban', 'Suburban', 'Rural'],
                         x='Home_Area', y='Annoyance', resamples=5000, random_seed=24624)


# %%
# calculate baseline paired effect sizes and plot
fig, ax = plt.subplots(figsize=(6, 3))
dataMD = dataloadBL.mean_diff.plot(contrast_bars=False,
                                   raw_ylim=(0, 10),
                                   raw_marker_size=0.05,
                                   swarmplot_kwargs={'alpha': 0.3},
                                   delta_dot_kwargs={'size': 1, 'side': 'left', 'alpha': 0.25, 'zorder': 1},
                                   delta_text=True,
                                   custom_palette={'Suburban': mycolours[0],
                                                   'Rural': mycolours[1],
                                                   'Urban': mycolours[3]},
                                   raw_desat=1,
                                   contrast_desat=1,
                                   contrast_marker_size=4,
                                   contrast_ylim=(0, 4),
                                   ax=ax,
                                   )

diffs = dataloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=-3.9, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.7 + ii + 0.09, y=-4.28, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtBAORDabestBase.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtBAORDabestBase.pdf"),
                format='pdf', bbox_inches='tight')
out = pd.concat([out, dataloadBL.mean_diff.statistical_tests])

dataMD;
dataloadBL.mean_diff.statistical_tests


# %%
# calculate and display Cohen's d
outd = pd.concat([outd, dataloadBL.cohens_d.statistical_tests])
dataloadBL.cohens_d.statistical_tests


# %% [markdown]
# ## Parts A and B combined
# 
# ### AAM attitude

# %%
# select subset of data for analysis and sort
data = dataBySubjTest
data = data[data['UASLAeq'] != "Baseline"]
data = data.loc[:, ['ID', 'StimFile', 'AAM_attitude', 'Annoyance']]


# %%
# assign data for processing
dataloadBL = dabest.load(ps_adjust=True, data=data, idx=['Supportive', 'Ambivalent', 'Concerned', 'Neutral'],
                         x='AAM_attitude', y='Annoyance', resamples=5000, random_seed=3487)


# %%
len(data)

# %%
# calculate baseline paired effect sizes and plot
fig, ax = plt.subplots(figsize=(6, 3))
dataMD = dataloadBL.mean_diff.plot(contrast_bars=False,
                                   raw_ylim=(0, 10),
                                   raw_marker_size=0.004,
                                   swarmplot_kwargs={'alpha': 0.3},
                                   delta_dot_kwargs={'size': 1, 'side': 'left', 'alpha': 0.25, 'zorder': 1},
                                   delta_text=True,
                                   custom_palette={'Neutral': mycolours[0],
                                                   'Ambivalent': mycolours[1],
                                                   'Concerned': mycolours[2],
                                                   'Supportive': mycolours[3]},
                                   raw_desat=1,
                                   contrast_desat=1,
                                   contrast_marker_size=4,
                                   contrast_ylim=(-1, 3),
                                   ax=ax,
                                  )

diffs = dataloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=-3.9, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.7 + ii + 0.09, y=-4.28, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtsABAttDabestBase.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtsABAttDabestBase.pdf"),
                format='pdf', bbox_inches='tight')
out = pd.concat([out, dataloadBL.mean_diff.statistical_tests])

dataMD;
dataloadBL.mean_diff.statistical_tests


# %%
# calculate and display Cohen's d
outd = pd.concat([outd, dataloadBL.cohens_d.statistical_tests])
dataloadBL.cohens_d.statistical_tests


# %% [markdown]
# #### Test the change in annoyance outcome variable

# %%
# select subset of data for analysis and sort
data = dataBySubjTest
data = data[data['UASLAeq'] != "Baseline"]
data = data.loc[:, ['ID', 'StimFile', 'AAM_attitude', 'dAnnoyance']]


# %%
# assign data for processing
dataloadBL = dabest.load(ps_adjust=True, data=data, idx=['Supportive', 'Ambivalent', 'Concerned', 'Neutral'],
                         x='AAM_attitude', y='dAnnoyance', resamples=5000, random_seed=3487)


# %%
# calculate baseline paired effect sizes and plot
fig, ax = plt.subplots(figsize=(6, 3))
dataMD = dataloadBL.mean_diff.plot(contrast_bars=False,
                                   raw_ylim=(-10, 10),
                                   raw_marker_size=0.0025,
                                   swarmplot_kwargs={'alpha': 0.3},
                                   delta_dot_kwargs={'size': 1, 'side': 'left', 'alpha': 0.25, 'zorder': 1},
                                   delta_text=True,
                                   custom_palette={'Neutral': mycolours[0],
                                                   'Ambivalent': mycolours[1],
                                                   'Concerned': mycolours[2],
                                                   'Supportive': mycolours[3]},
                                   raw_desat=1,
                                   contrast_desat=1,
                                   contrast_marker_size=4,
                                   contrast_ylim=(-1, 3),
                                   ax=ax,
                                   )

diffs = dataloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=-20.9, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.7 + ii + 0.09, y=-21.279999999999998, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

#if saveplots:
#    plt.savefig(os.path.join(outFigPath, "svg", "PtsABAttDabestBase.svg"),
 #               #format='svg', bbox_inches='tight')
#    plt.savefig(os.path.join(outFigPath, "pdf", "PtsABAttDabestBase.pdf"),
 #               format='pdf', bbox_inches='tight')


dataMD;
dataloadBL.mean_diff.statistical_tests


# %%
dataloadBL.cohens_d.statistical_tests


# %% [markdown]
# ### Home residence area

# %%
# select subset of data for analysis and sort
data = dataBySubjTest
data = data[data['UASLAeq'] != "Baseline"]
data = data.loc[:, ['ID', 'StimFile', 'UASLAeq', 'Home_Area', 'Annoyance']]


# %%
# assign data for processing
dataloadBL = dabest.load(ps_adjust=True, data=data, idx=['Urban', 'Suburban', 'Rural'],
                         x='Home_Area', y='Annoyance', resamples=5000, random_seed=931)


# %%
# calculate baseline paired effect sizes and plot
fig, ax = plt.subplots(figsize=(6, 3))
dataMD = dataloadBL.mean_diff.plot(contrast_bars=False,
                                   raw_ylim=(0, 10),
                                   raw_marker_size=0.001,
                                   swarmplot_kwargs={'alpha': 0.3},
                                   delta_dot_kwargs={'size': 1, 'side': 'left', 'alpha': 0.25, 'zorder': 1},
                                   delta_text=True,
                                   custom_palette={'Suburban': mycolours[0],
                                                   'Rural': mycolours[1],
                                                   'Urban': mycolours[3]},
                                   raw_desat=1,
                                   contrast_desat=1,
                                   contrast_marker_size=4,
                                   contrast_ylim=(-2, 2),
                                   ax=ax,
                                   )

diffs = dataloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=-3.9, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.7 + ii + 0.09, y=-4.28, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtsABAORDabestBase.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtsABAORDabestBase.pdf"),
                format='pdf', bbox_inches='tight')
out = pd.concat([out, dataloadBL.mean_diff.statistical_tests])

dataMD;
dataloadBL.mean_diff.statistical_tests


# %%
# calculate and display Cohen's d
outd = pd.concat([outd, dataloadBL.cohens_d.statistical_tests])
dataloadBL.cohens_d.statistical_tests


# %% [markdown]
# #### Change in annoyance

# %%
# select subset of data for analysis and sort
data = dataBySubjTest
data = data[data['UASLAeq'] != "Baseline"]
data = data.loc[:, ['ID', 'StimFile', 'UASLAeq', 'Home_Area', 'dAnnoyance']]


# %%
# assign data for processing
dataloadBL = dabest.load(ps_adjust=True, data=data, idx=['Urban', 'Suburban', 'Rural'],
                         x='Home_Area', y='dAnnoyance', resamples=5000, random_seed=931)


# %%
# calculate baseline paired effect sizes and plot
fig, ax = plt.subplots(figsize=(6, 3))
dataMD = dataloadBL.mean_diff.plot(contrast_bars=False,
                                   raw_ylim=(-10, 10),
                                   raw_marker_size=0.001,
                                   swarmplot_kwargs={'alpha': 0.3},
                                   delta_dot_kwargs={'size': 1, 'side': 'left', 'alpha': 0.25, 'zorder': 1},
                                   delta_text=True,
                                   custom_palette={'Suburban': mycolours[0],
                                                   'Rural': mycolours[1],
                                                   'Urban': mycolours[3]},
                                   raw_desat=1,
                                   contrast_desat=1,
                                   contrast_marker_size=4,
                                   contrast_ylim=(0, 3),
                                   ax=ax,
                                   )

diffs = dataloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=-19.9, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.7 + ii + 0.09, y=-20.185, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

#if saveplots:
#    plt.savefig(os.path.join(outFigPath, "svg", "PtsABAORDabestBase.svg"),
#                format='svg', bbox_inches='tight')
#    plt.savefig(os.path.join(outFigPath, "pdf", "PtsABAORDabestBase.pdf"),
  #              format='pdf', bbox_inches='tight')
#out = pd.concat([out, dataloadBL.mean_diff.statistical_tests])

dataMD;
dataloadBL.mean_diff.statistical_tests


# %%
# calculate and display Cohen's d
#outd = pd.concat([outd, dataloadBL.cohens_d.statistical_tests])
dataloadBL.cohens_d.statistical_tests


# %% [markdown]
# ### Area soundscape character

# %%
# select subset of data for analysis and sort
data = dataBySubjTest
data = data[data['UASLAeq'] != "Baseline"]
data = data.loc[:, ['ID', 'StimFile', 'Area_soundscape', 'Annoyance']]


# %%
# assign data for processing
dataloadBL = dabest.load(ps_adjust=True, data=data, idx=['Monotonous', 'Calm', 'Chaotic', 'Vibrant'],
                         x='Area_soundscape', y='Annoyance', resamples=5000, random_seed=21064694)


# %%
# calculate baseline paired effect sizes and plot
fig, ax = plt.subplots(figsize=(6, 3))
dataMD = dataloadBL.mean_diff.plot(contrast_bars=False,
                                   raw_ylim=(0, 10),
                                   raw_marker_size=0.002,
                                   swarmplot_kwargs={'alpha': 0.3},
                                   delta_dot_kwargs={'size': 1, 'side': 'left', 'alpha': 0.25, 'zorder': 1},
                                   delta_text=True,
                                   custom_palette={'Monotonous': mycolours[0],
                                                   'Calm': mycolours[1],
                                                   'Chaotic': mycolours[2],
                                                   'Vibrant': mycolours[3]},
                                   raw_desat=1,
                                   contrast_desat=1,
                                   contrast_marker_size=4,
                                   contrast_ylim=(-1, 2),
                                   ax=ax,
                                   )

diffs = dataloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=-3.9, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.7 + ii + 0.09, y=-4.185, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtsABScapeDabestBase.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtsABScapeDabestBase.pdf"),
                format='pdf', bbox_inches='tight')
out = pd.concat([out, dataloadBL.mean_diff.statistical_tests])

dataMD;
dataloadBL.mean_diff.statistical_tests


# %%
# calculate and display Cohen's d
outd = pd.concat([outd, dataloadBL.cohens_d.statistical_tests])
dataloadBL.cohens_d.statistical_tests


# %% [markdown]
# #### Change in annoyance

# %%
# select subset of data for analysis and sort
data = dataBySubjTest
data = data[data['UASLAeq'] != "Baseline"]
data = data.loc[:, ['ID', 'StimFile', 'Area_soundscape', 'dAnnoyance']]


# %%
# assign data for processing
dataloadBL = dabest.load(ps_adjust=True, data=data, idx=['Monotonous', 'Calm', 'Chaotic', 'Vibrant'],
                         x='Area_soundscape', y='dAnnoyance', resamples=5000, random_seed=21064694)


# %%
# calculate baseline paired effect sizes and plot
fig, ax = plt.subplots(figsize=(6, 3))
dataMD = dataloadBL.mean_diff.plot(contrast_bars=False,
                                   raw_ylim=(-10, 10),
                                   raw_marker_size=0.0015,
                                   swarmplot_kwargs={'alpha': 0.3},
                                   delta_dot_kwargs={'size': 1, 'side': 'left', 'alpha': 0.25, 'zorder': 1},
                                   delta_text=True,
                                   custom_palette={'Monotonous': mycolours[0],
                                                   'Calm': mycolours[1],
                                                   'Chaotic': mycolours[2],
                                                   'Vibrant': mycolours[3]},
                                   raw_desat=1,
                                   contrast_desat=1,
                                   contrast_marker_size=4,
                                   contrast_ylim=(-1, 2),
                                   ax=ax,
                                   )

diffs = dataloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=-18.9, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.7 + ii + 0.09, y=-19.185, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

#if saveplots:
 #   plt.savefig(os.path.join(outFigPath, "svg", "PtsABScapeDabestBase.svg"),
  #              format='svg', bbox_inches='tight')
#    plt.savefig(os.path.join(outFigPath, "pdf", "PtsABScapeDabestBase.pdf"),
#                format='pdf', bbox_inches='tight')
#out = pd.concat([out, dataloadBL.mean_diff.statistical_tests])

dataMD;
dataloadBL.mean_diff.statistical_tests


# %%
# calculate and display Cohen's d
#outd = pd.concat([outd, dataloadBL.cohens_d.statistical_tests])
dataloadBL.cohens_d.statistical_tests


# %% [markdown]
# ### Nationality region

# %%
# select subset of data for analysis and sort
data = dataBySubjTest
data = data[data['UASLAeq'] != "Baseline"]
data = data.loc[:, ['ID', 'StimFile', 'UASLAeq', 'NationGeo', 'Annoyance']]


# %%
# assign data for processing
dataloadBL = dabest.load(ps_adjust=True, data=data, idx=['UK', 'Africa', 'EastAsia', 'MidEast', 'SouthAsia', 'Europe', 'SouthAmerica', 'Australasia'],
                         x='NationGeo', y='Annoyance', resamples=5000, random_seed=65494894)


# %%
# calculate baseline paired effect sizes and plot
fig, ax = plt.subplots(figsize=(13.71, 3))
dataMD = dataloadBL.mean_diff.plot(contrast_bars=False,
                                   raw_ylim=(0, 10),
                                   raw_marker_size=0.001,
                                   swarmplot_kwargs={'alpha': 0.3},
                                   delta_dot_kwargs={'size': 1, 'side': 'left', 'alpha': 0.25, 'zorder': 1},
                                   delta_text=True,
                                   custom_palette={'UK': mycolours[0],
                                                   'Africa': mycolours[1],
                                                   'EastAsia': mycolours[3],
                                                   'MidEast': mycolours[4],
                                                   'SouthAsia': mycolours[5],
                                                   'Europe': mycolours[6],
                                                   'SouthAmerica': mycolours[7],
                                                   'Australasia': mycolours[8]},
                                   raw_desat=1,
                                   contrast_desat=1,
                                   contrast_marker_size=4,
                                   contrast_ylim=(-2.5, 4.5),
                                   ax=ax,
                                   )

diffs = dataloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=-3.9, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.7 + ii + 0.09, y=-4.5649999999999995, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

dataMD;
dataloadBL.mean_diff.statistical_tests


# %%
# calculate and display Cohen's d
dataloadBL.cohens_d.statistical_tests


# %% [markdown]
# #### Change in annoyance (all)

# %%
# select subset of data for analysis and sort
data = dataBySubjTest
data = data[data['UASLAeq'] != "Baseline"]
data = data.loc[:, ['ID', 'StimFile', 'UASLAeq', 'NationGeo', 'dAnnoyance']]


# %%
# assign data for processing
dataloadBL = dabest.load(ps_adjust=True, data=data, idx=['UK', 'Africa', 'EastAsia', 'MidEast', 'SouthAsia', 'Europe', 'SouthAmerica', 'Australasia'],
                         x='NationGeo', y='dAnnoyance', resamples=5000, random_seed=65494894)


# %%
# calculate baseline paired effect sizes and plot
fig, ax = plt.subplots(figsize=(13.71, 3))
dataMD = dataloadBL.mean_diff.plot(contrast_bars=False,
                                   raw_ylim=(-10, 10),
                                   raw_marker_size=0.001,
                                   swarmplot_kwargs={'alpha': 0.3},
                                   delta_dot_kwargs={'size': 1, 'side': 'left', 'alpha': 0.25, 'zorder': 1},
                                   delta_text=True,
                                   custom_palette={'UK': mycolours[0],
                                                   'Africa': mycolours[1],
                                                   'EastAsia': mycolours[3],
                                                   'MidEast': mycolours[4],
                                                   'SouthAsia': mycolours[5],
                                                   'Europe': mycolours[6],
                                                   'SouthAmerica': mycolours[7],
                                                   'Australasia': mycolours[8]},
                                   raw_desat=1,
                                   contrast_desat=1,
                                   contrast_marker_size=4,
                                   contrast_ylim=(-4.5, 4.5),
                                   ax=ax,
                                   )

diffs = dataloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=-18.9, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.7 + ii + 0.09, y=-19.755, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

dataMD;
dataloadBL.mean_diff.statistical_tests


# %% [markdown]
# #### 2 groups: UK vs Africa

# %%
# assign data for processing
dataloadBL = dabest.load(ps_adjust=True, data=data, idx=['UK', 'Africa'],
                         x='NationGeo', y='Annoyance', resamples=5000, random_seed=65494894)

# %%
# calculate baseline paired effect sizes and plot
fig, ax = plt.subplots(figsize=(4.29, 3))
dataMD = dataloadBL.mean_diff.plot(contrast_bars=False,
                                   raw_ylim=(0, 10),
                                   raw_marker_size=0.001,
                                   swarmplot_kwargs={'alpha': 0.3},
                                   delta_dot_kwargs={'size': 1, 'side': 'left', 'alpha': 0.25, 'zorder': 1},
                                   delta_text=True,
                                   custom_palette={'UK': mycolours[0],
                                                   'Africa': mycolours[1]},
                                   raw_desat=1,
                                   contrast_desat=1,
                                   contrast_marker_size=4,
                                   contrast_ylim=(-2, 0),
                                   ax=ax,
                                   )

diffs = dataloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=1.4 + ii, y=8, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=1.4 + ii + 0.09, y=7.81, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtsABNationDabestBase.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtsABNationDabestBase.pdf"),
                format='pdf', bbox_inches='tight')
out = pd.concat([out, dataloadBL.mean_diff.statistical_tests])

dataMD;
dataloadBL.mean_diff.statistical_tests


# %%
# calculate and display Cohen's d
outd = pd.concat([outd, dataloadBL.cohens_d.statistical_tests])
dataloadBL.cohens_d.statistical_tests


# %% [markdown]
# #### Combining all non-UK regions into one 'other' group

# %%
# create new subset of data for analysis and sort
data = dataBySubjTest.copy()
data.loc[(data['NationGeo'] != "UK") & ~data['NationGeo'].isna(), 'NationGeo'] = "Other"
data = data[data['UASLAeq'] != "Baseline"]
data = data.loc[:, ['ID', 'StimFile', 'UASLAeq', 'NationGeo', 'Annoyance']]


# %%
# assign data for processing
dataloadBL = dabest.load(ps_adjust=True, data=data, idx=['UK', 'Other'],
                         x='NationGeo', y='Annoyance', resamples=5000, random_seed=5954984962)


# %%
# calculate baseline paired effect sizes and plot
fig, ax = plt.subplots(figsize=(4.29, 3))
dataMD = dataloadBL.mean_diff.plot(contrast_bars=False,
                                   raw_ylim=(0, 10),
                                   raw_marker_size=0.001,
                                   swarmplot_kwargs={'alpha': 0.3},
                                   delta_dot_kwargs={'size': 1, 'side': 'left', 'alpha': 0.25, 'zorder': 1},
                                   delta_text=True,
                                   custom_palette={'UK': mycolours[0],
                                                   'Other': mycolours[1]},
                                   raw_desat=1,
                                   contrast_desat=1,
                                   contrast_marker_size=4,
                                   contrast_ylim=(-2, 0),
                                   ax=ax,
                                   )

diffs = dataloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=1.4 + ii, y=8, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=1.4 + ii + 0.09, y=7.81, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtsABNation2DabestBase.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtsABNation2DabestBase.pdf"),
                format='pdf', bbox_inches='tight')


dataMD;
dataloadBL.mean_diff.statistical_tests


# %%
out = pd.concat([out, dataloadBL.mean_diff.statistical_tests])
# calculate and display Cohen's d
outd = pd.concat([outd, dataloadBL.cohens_d.statistical_tests])
dataloadBL.cohens_d.statistical_tests


# %% [markdown]
# ##### Change in annoyance

# %%
# create new subset of data for analysis and sort
data = dataBySubjTest.copy()
data.loc[(data['NationGeo'] != "UK") & ~data['NationGeo'].isna(), 'NationGeo'] = "Other"
data = data[data['UASLAeq'] != "Baseline"]
data = data.loc[:, ['ID', 'StimFile', 'UASLAeq', 'NationGeo', 'dAnnoyance']]


# %%
# assign data for processing
dataloadBL = dabest.load(ps_adjust=True, data=data, idx=['UK', 'Other'],
                         x='NationGeo', y='dAnnoyance', resamples=5000, random_seed=5954984962)


# %%
# calculate baseline paired effect sizes and plot
fig, ax = plt.subplots(figsize=(4.29, 3))
dataMD = dataloadBL.mean_diff.plot(contrast_bars=False,
                                   raw_ylim=(-10, 10),
                                   raw_marker_size=0.001,
                                   swarmplot_kwargs={'alpha': 0.3},
                                   delta_dot_kwargs={'size': 1, 'side': 'left', 'alpha': 0.25, 'zorder': 1},
                                   delta_text=True,
                                   custom_palette={'UK': mycolours[0],
                                                   'Other': mycolours[1]},
                                   raw_desat=1,
                                   contrast_desat=1,
                                   contrast_marker_size=4,
                                   contrast_ylim=(-2, 0),
                                   ax=ax,
                                   )

diffs = dataloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=1.4 + ii, y=8, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=1.4 + ii + 0.09, y=7.81, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

#if saveplots:
#    plt.savefig(os.path.join(outFigPath, "svg", "PtsABNation2DabestBase.svg"),
#                format='svg', bbox_inches='tight')
#    plt.savefig(os.path.join(outFigPath, "pdf", "PtsABNation2DabestBase.pdf"),
#                format='pdf', bbox_inches='tight')
#out = pd.concat([out, dataloadBL.mean_diff.statistical_tests])

dataMD;
dataloadBL.mean_diff.statistical_tests


# %%
# calculate and display Cohen's d
#outd = pd.concat([outd, dataloadBL.cohens_d.statistical_tests])
dataloadBL.cohens_d.statistical_tests


# %% [markdown]
# Combining all non-UK regions except Africa into one 'other' group

# %%
# create new subset of data for analysis and sort
data = dataBySubjTest.copy()
data.loc[(data['NationGeo'] != "UK") & (data['NationGeo'] != "Africa") & ~data['NationGeo'].isna(), 'NationGeo'] = "Other"
data = data[data['UASLAeq'] != "Baseline"]
data = data.loc[:, ['ID', 'StimFile', 'UASLAeq', 'NationGeo', 'Annoyance']]


# %%
# assign data for processing
dataloadBL = dabest.load(ps_adjust=True, data=data, idx=['UK', 'Africa', 'Other'],
                         x='NationGeo', y='Annoyance', resamples=5000, random_seed=454515665)


# %%
# calculate baseline paired effect sizes and plot
fig, ax = plt.subplots(figsize=(6, 3))
dataMD = dataloadBL.mean_diff.plot(contrast_bars=False,
                                   raw_ylim=(0, 10),
                                   raw_marker_size=0.001,
                                   swarmplot_kwargs={'alpha': 0.3},
                                   delta_dot_kwargs={'size': 1, 'side': 'left', 'alpha': 0.25, 'zorder': 1},
                                   delta_text=True,
                                   custom_palette={'UK': mycolours[0],
                                                   'Africa': mycolours[1],
                                                   'Other': mycolours[2]},
                                   raw_desat=1,
                                   contrast_desat=1,
                                   contrast_marker_size=4,
                                   contrast_ylim=(-2, 2),
                                   ax=ax,
                                   )

diffs = dataloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=-3.9, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.7 + ii + 0.09, y=-4.28, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtsABNation3DabestBase.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtsABNation3DabestBase.pdf"),
                format='pdf', bbox_inches='tight')

dataMD;
dataloadBL.mean_diff.statistical_tests


# %%
out = pd.concat([out, dataloadBL.mean_diff.statistical_tests])
# calculate and display Cohen's d
outd = pd.concat([outd, dataloadBL.cohens_d.statistical_tests])
dataloadBL.cohens_d.statistical_tests


# %% [markdown]
# #### Combine Africa and South Asia into one group

# %%
dataBySubjTest['NativeLang'].unique()


# %%
# create new subset of data for analysis and sort
data = dataBySubjTest.copy()
data.loc[(data['NationGeo'] == "Africa") | (data['NationGeo'] == "SouthAsia"), 'NationGeo'] = "Africa_SAsia"
data.loc[(data['NationGeo'] != "UK") & (data['NationGeo'] != "Africa_SAsia") & ~data['NationGeo'].isna(), 'NationGeo'] = "Other"
data = data[data['UASLAeq'] != "Baseline"]
data = data.loc[:, ['ID', 'StimFile', 'UASLAeq', 'NationGeo', 'Annoyance']]


# %%
# assign data for processing
dataloadBL = dabest.load(ps_adjust=True, data=data, idx=['UK', 'Africa_SAsia', 'Other'],
                         x='NationGeo', y='Annoyance', resamples=5000, random_seed=761278622)

# %%
# calculate baseline paired effect sizes and plot
fig, ax = plt.subplots(figsize=(6, 3))
dataMD = dataloadBL.mean_diff.plot(contrast_bars=False,
                                   raw_ylim=(0, 10),
                                   raw_marker_size=0.001,
                                   swarmplot_kwargs={'alpha': 0.3},
                                   delta_dot_kwargs={'size': 1, 'side': 'left', 'alpha': 0.25, 'zorder': 1},
                                   delta_text=True,
                                   custom_palette={'UK': mycolours[0],
                                                   'Africa_SAsia': mycolours[1],
                                                   'Other': mycolours[2]},
                                   raw_desat=1,
                                   contrast_desat=1,
                                   contrast_marker_size=4,
                                   contrast_ylim=(-2, 1.5),
                                   ax=ax,
                                   )

diffs = dataloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.25 + ii, y=-3.9, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.25 + ii + 0.09, y=-4.233, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

if saveplots:
    plt.savefig(os.path.join(outFigPath, "svg", "PtsABNation4DabestBase.svg"),
                format='svg', bbox_inches='tight')
    plt.savefig(os.path.join(outFigPath, "pdf", "PtsABNation4DabestBase.pdf"),
                format='pdf', bbox_inches='tight')


dataMD;
dataloadBL.mean_diff.statistical_tests


# %%
out = pd.concat([out, dataloadBL.mean_diff.statistical_tests])
# calculate and display Cohen's d
outd = pd.concat([outd, dataloadBL.cohens_d.statistical_tests])
dataloadBL.cohens_d.statistical_tests


# %% [markdown]
# ##### Change in annoyance

# %%
# create new subset of data for analysis and sort
data = dataBySubjTest.copy()
data.loc[(data['NationGeo'] == "Africa") | (data['NationGeo'] == "SouthAsia"), 'NationGeo'] = "Africa_SAsia"
data.loc[(data['NationGeo'] != "UK") & (data['NationGeo'] != "Africa_SAsia") & ~data['NationGeo'].isna(), 'NationGeo'] = "Other"
data = data[data['UASLAeq'] != "Baseline"]
data = data.loc[:, ['ID', 'StimFile', 'UASLAeq', 'NationGeo', 'dAnnoyance']]


# %%
# assign data for processing
dataloadBL = dabest.load(ps_adjust=True, data=data, idx=['UK', 'Africa_SAsia', 'Other'],
                         x='NationGeo', y='dAnnoyance', resamples=5000, random_seed=761278622)

# %%
# calculate baseline paired effect sizes and plot
fig, ax = plt.subplots(figsize=(6, 3))
dataMD = dataloadBL.mean_diff.plot(contrast_bars=False,
                                   raw_ylim=(-10, 10),
                                   raw_marker_size=0.001,
                                   swarmplot_kwargs={'alpha': 0.3},
                                   delta_dot_kwargs={'size': 1, 'side': 'left', 'alpha': 0.25, 'zorder': 1},
                                   delta_text=True,
                                   custom_palette={'UK': mycolours[0],
                                                   'Africa_SAsia': mycolours[1],
                                                   'Other': mycolours[2]},
                                   raw_desat=1,
                                   contrast_desat=1,
                                   contrast_marker_size=4,
                                   contrast_ylim=(-2, 1.5),
                                   ax=ax,
                                   )

diffs = dataloadBL.mean_diff.statistical_tests[['difference', 'bca_low', 'bca_high', 'pvalue_permutation']].values
# label results
contrast_ax = ax.contrast_axes
for ii, diff in enumerate(diffs):
    vals = display_round(diff[1:3], digits=2, floor=False)
    contrast_ax.text(x=0.7 + ii, y=-19.9, s="95%: [" + " ".join(vals) + "]", fontsize=10)
    contrast_ax.text(x=0.7 + ii + 0.09, y=-20.232999999999997, s=r"$p$:   " + display_round(diff[-1], digits=3), fontsize=10)

#if saveplots:
 #   plt.savefig(os.path.join(outFigPath, "svg", "PtsABNation4DabestBase.svg"),
 #               format='svg', bbox_inches='tight')
 #   plt.savefig(os.path.join(outFigPath, "pdf", "PtsABNation4DabestBase.pdf"),
 #               format='pdf', bbox_inches='tight')
#out = pd.concat([out, dataloadBL.mean_diff.statistical_tests])

dataMD;
dataloadBL.mean_diff.statistical_tests


# %%
# calculate and display Cohen's d
#outd = pd.concat([outd, dataloadBL.cohens_d.statistical_tests])
dataloadBL.cohens_d.statistical_tests


# %% [markdown]
# ## Save output data to file

# %%
if savedata:
    out.to_csv(os.path.join(outDataPath, "dabest_mean_diff.csv"))
    outd.to_csv(os.path.join(outDataPath, "dabest_cohens_d.csv"))

# %% [markdown]
# # Comparison with DroneNoise data

# %%
dataDNoiseWide = pd.read_csv(r"C:\Users\m_lot\OneDrive - University of Salford\REFMAP General\03 Experiment\Experiment 1\Analysis\Comparison_data\DroneNoise2022AnnoyWide.csv")

dataDNoiseLong = dataDNoiseWide.melt(id_vars='ID', var_name='StimFile', value_name='Annoyance')
dataDNoiseLong.head()

# %%
dataDNoiseLong['UASOperation'] = "Flyby"
dataDNoiseLong.loc[dataDNoiseLong['StimFile'].str.contains("Landing"), 'UASOperation'] = "Landing"
dataDNoiseLong.loc[dataDNoiseLong['StimFile'].str.contains("Takeoff"), 'UASOperation'] = "Takeoff"
dataDNoiseLong.loc[dataDNoiseLong['StimFile'].str.contains("Hover"), 'UASOperation'] = "Hover"
dataDNoiseLong['Location'] = "Outdoors"
dataDNoiseLong.loc[dataDNoiseLong['StimFile'].str.contains("Part_Open"), 'Location'] = "Indoors_PO"
dataDNoiseLong.loc[dataDNoiseLong['StimFile'].str.contains("Closed"), 'Location'] = "Indoors_Cl"
dataDNoiseLong['UASType'] = "Typhoon"
dataDNoiseLong.loc[dataDNoiseLong['StimFile'].str.contains("GD28X"), 'UASType'] = "GD28X"
dataDNoiseLong.loc[dataDNoiseLong['StimFile'].str.contains("M200"), 'UASType'] = "M200"

# %%
# create dummy pairing ID for each participant to allow for paired analysis
grouping_cols = ['ID', 'Location', 'UASType']

# ngroup() automatically assigns an identical integer to matching sets
# across the different flight-operation conditions
dataDNoiseLong['dummyID'] = dataDNoiseLong.groupby(grouping_cols).ngroup() + 1
dataDNoiseLong.sort_values(by=['dummyID', 'ID', 'StimFile'], inplace=True)

dataDNoiseloadBL = dabest.load(ps_adjust=True, data=dataDNoiseLong, idx=("Flyby", "Takeoff", "Landing", "Hover"),
                               x='UASOperation', y='Annoyance', paired='baseline',
                               id_col='dummyID', resamples=5000, random_seed=646)

dataDNoiseBL = dataDNoiseloadBL.mean_diff.plot(contrast_bars=False, show_baseline_ec=True,
                                   raw_ylim=(-0.5, 10.5), color_col='Location',
                                               delta_dot_kwargs={'size': 1, 'side': 'left', 'alpha': 0.25, 'zorder': 1}, slopegraph_kwargs={'linewidth': 0.25, 'alpha': 0.4},
                                               delta_text=True,
                                               contrast_marker_size=4, contrast_ylim=(0, 2.5),
                                               legend_kwargs={'loc': 'lower center'})


# %%
# create dummy pairing ID for each participant to allow for paired analysis
grouping_cols = ['ID', 'UASOperation', 'UASType']

# ngroup() automatically assigns an identical integer to matching sets
# across the different location conditions
dataDNoiseLong['dummyID'] = dataDNoiseLong.groupby(grouping_cols).ngroup() + 1
dataDNoiseLong.sort_values(by=['dummyID', 'ID', 'StimFile'], inplace=True)

dataDNoiseloadBL = dabest.load(ps_adjust=True, data=dataDNoiseLong, idx=("Outdoors", "Indoors_PO", "Indoors_Cl"),
                               x='Location', y='Annoyance', paired='baseline',
                               id_col='dummyID', resamples=5000, random_seed=4546)

dataDNoiseBL = dataDNoiseloadBL.mean_diff.plot(contrast_bars=False, show_baseline_ec=True,
                                   raw_ylim=(-0.5, 10.5), color_col='UASOperation',
                                               delta_dot_kwargs={'size': 1, 'side': 'left', 'alpha': 0.25, 'zorder': 1}, slopegraph_kwargs={'linewidth': 0.25, 'alpha': 0.4},
                                               delta_text=True,
                                               contrast_marker_size=4, contrast_ylim=(-5, -2),
                                               legend_kwargs={'loc': 'lower center'})