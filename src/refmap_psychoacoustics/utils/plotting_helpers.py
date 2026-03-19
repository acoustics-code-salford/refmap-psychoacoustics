# -*- coding: utf-8 -*-
"""
plotting_helpers.py
-------------------

Helper functions for creating various plots.

Functions
---------

spiderplot(df, *, id_column, categories=None, ax=None,
           title=None, limits=None, max_values=None,
              scalers=None, padding=1.25, annotate=True,
                offset=False, fmt='{:.3f}', labelsize=10,
                ticksize=11, legendsize=10, palette=None)
    Create a flexible spider/radar plot from a DataFrame.

violin(data, x, y, xCats=None, xjitter=0.05, yjitter=0, yjitter_type="auto",
         palette=None, figsize=(7, 4), alpha_pt=0.2, size_pt=10, seed=123, ax=None)
     Create a violin plot with optional jittered points.

violin_split(data, x, y, hue, xCats=None, xjitter=0.05, yjitter=0, yjitter_type="auto",
                med_trace=False, palette=None, alpha_pt=0.25, violin_width=0.5,
                figsize=(7, 4.65), size_pt=10, seed=123, ax=None)
        Create a split violin plot with optional jittered points and median trace.

violin_facet(data, plot_func, facet_col=None, facet_row=None, col_order=None,
             row_order=None, wrap_cols=None, sharex=True, sharey=True,
                figsize=(10, 6), panel_titles=True, title_fmt="{row} | {col}",
                legend=True, legend_loc="upper center", legend_ncol=None,
                legend_bbox_to_anchor=None, xlabel=None, ylabel=None, xticks=None,
                yticks=None, **plot_kwargs)
    Create a faceted grid of violin plots using a custom plotting function.

Requirements
------------
matplotlib
seaborn
numpy

Ownership and Quality Assurance
-------------------------------
Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
Institution: University of Salford
 
Date created: 17/03/2026
Date last modified: 17/03/2026
Python version: 3.11

Copyright statement: This file and code is part of work undertaken within
the RefMap project (www.refmap.eu), and is subject to licence as detailed
in the code repository
(https://github.com/acoustics-code-salford/refmap-psychoacoustics)

Checked by:
Date last checked:

"""

# %% imports

import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
import seaborn as sns
import pandas as pd
from scipy import stats

# %% spiderplot
def spiderplot(df, *, id_column, categories=None, ax=None,
               title=None, limits=None, max_values=None,
               scalers=None, padding=1.25, annotate=True,
               offset=False, fmt='{:.3f}', labelsize=10,
               ticksize=11, legendsize=10, titlesize=12,
               palette=None):
    """
    Flexible spider/radar plot.

    Inputs
    ------
    
    df: pandas DataFrame
        DataFrame containing the data to plot. Must include a column specified by
        `id_column` for labeling each line.
    
    id_column: str
        Column name to use for labeling each line (e.g., category or group name).
    
    categories: list of str, optional
        List of column names to include in the plot. If None, all numeric columns are used.

    ax: matplotlib Axes, optional
        Axes to plot on. If None, a new figure and axes are created.
    
    title: str, optional
        Title for the plot.

    limits: dict: {cat: (min, max)}
        Optional limits for each category. If not provided, limits are determined from the data.

    max_values: dict: {cat: max_val}
        Legacy support for specifying maximum values for categories. If provided, these will be
        used to set the upper limit for normalisation.

    scalers: dict: {cat: func}
        Optional dictionary of scaling functions for each category. If provided, these functions
        will be applied to the data before plotting.
    
    padding: float, default=1.25
        Padding for the radial axis limits.
    
    annotate: bool, default=True
        Whether to annotate each point with its actual value.

    offset: bool, default=False
        Whether to apply a small offset to series to help prevent overlap.
    
    fmt: str, default="{:.3f}"
        Format string for annotations.
    
    labelsize: int, default=10
        Size of the labels for the annotations.
    
    ticksize: int, default=11
        Size of the tick labels.
    
    legendsize: int, default=10
        Size of the legend text.

    titlesize: int, default=12
        Size of the title text.
    
    palette: list of colors, optional
        List of colors to use for the plot lines. If None, the default color cycle is used.
    
        
    Returns
    -------
    fig, ax: matplotlib Figure and Axes
        The figure and axes objects containing the plot.

    """

    def _place_annotations(ax, angles, values, actuals, color, fmt, labelsize,
                           base_offset=0.05, placed=None):
        """
        Place annotations with simple collision avoidance and angle-aware alignment.
        """
        if placed is None:
            placed = []

        for ang, val, actual in zip(angles[:-1], values[:-1], actuals):
            label = fmt.format(actual) if isinstance(actual, (float, int)) else str(actual)

            # --- Radial offset ---
            r = val + base_offset

            # --- Angle-based alignment ---
            ang_deg = np.degrees(ang)
            if 90 < ang_deg < 270:
                ha = "right"
            else:
                ha = "left"

            va = "center"

            # --- Simple collision avoidance ---
            for _ in range(5):  # limit iterations
                too_close = False
                for (prev_ang, prev_r) in placed:
                    if abs(prev_ang - ang) < 0.1 and abs(prev_r - r) < 0.05:
                        r += 0.05  # nudge outward
                        too_close = True
                        break
                if not too_close:
                    break

            placed.append((ang, r))

            ax.text(ang, r, label, size=labelsize, color=color, ha=ha, va=va)

    # --- Select categories ---
    if categories is None:
        categories = df.select_dtypes(include=np.number).columns.tolist()

    if not categories:
        raise ValueError("No numeric categories found.")

    data = df[categories]
    ids = df[id_column].tolist()

    # --- Build limits ---
    if limits is None:
        limits = {}
        for col in categories:
            col_data = data[col].values
            min_val = 0
            max_val = np.max(col_data)

            if max_values is not None:
                max_val = max_values.get(col, max_val)

            limits[col] = (min_val, padding*max_val)

    # --- Normalise function ---
    def normalise(val, vmin, vmax):
        if vmax == vmin:
            return 0.0
        return (val - vmin)/(vmax - vmin)

    # --- Apply scaling ---
    norm_data = {}
    for col in categories:
        vmin, vmax = limits[col]

        if scalers and col in scalers:
            norm_data[col] = scalers[col](data[col].values)
        else:
            norm_data[col] = np.array([normalise(v, vmin, vmax)
                                       for v in data[col].values])

    # --- Angles ---
    num_vars = len(categories)
    angles = np.linspace(0, 2*np.pi, num_vars, endpoint=False)
    angles = np.concatenate([angles, [angles[0]]])

    # --- Figure/Axes handling ---
    if ax is None:
        fig, ax = plt.subplots(subplot_kw=dict(polar=True), figsize=(6, 6))
    else:
        fig = ax.figure

    # --- Colors ---
    if palette is not None:
        ax.set_prop_cycle(cycler(color=palette))

    placed_annotations = []
    # --- Plot each row ---
    for ii, label in enumerate(ids):
        values = [norm_data[col][ii] for col in categories]
        actuals = [data[col].iloc[ii] for col in categories]

        values = np.concatenate([values, [values[0]]])

        line, = ax.plot(angles, values, label=label)
        ax.fill(angles, values, alpha=0.1)

        if offset:
            series_offset = 0.065*ii  # simple offset to help separate lines
        else:
            series_offset = 0

        if annotate:
            _place_annotations(ax, angles, values, actuals,
                               line.get_color(), fmt, labelsize, base_offset=0.05 + series_offset,
                               placed=placed_annotations)
    # --- Axis formatting ---
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories, fontsize=ticksize)
    ax.set_yticklabels([])

    # Optional radial grid (0–1 scale)
    ax.set_ylim(0, 1)

    # --- Legend & title ---
    ax.legend(fontsize=legendsize, loc="upper right", bbox_to_anchor=(1.2, 1.1))

    if title:
        ax.set_title(title, fontsize=titlesize, pad=20)

    return fig, ax


# %% violin
def violin(data, x, y, xCats=None, xjitter=0.05, yjitter=0, yjitter_type="auto",
           seed=123, palette=None, figsize=(7, 4), size_pt=10, alpha_pt=0.2,
           ax=None):
    
    # remove any rows with missing data in the relevant columns
    data = data.dropna(subset=[x, y])

    rng = np.random.default_rng(seed)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    if palette is None:
        palette = list(mpl.colors.TABLEAU_COLORS.values())

    # ordered categories
    if xCats is None:
        xCats = list(data[x].sort_values().unique())


    y_data = []
    x_positions = []

    for ii, xcat in enumerate(xCats):

        sub = data[data[x] == xcat][y]

        if sub.empty:
            y_data.append([])
            x_positions.append(ii)
            continue

        yvals = sub.to_numpy().astype(float)

        # optional y jitter
        if yjitter:
            yvals = yvals + rng.normal(0, yjitter, len(yvals))

        xvals = np.full(len(yvals), ii, dtype=float)

        # ----------------------
        # density-aware x jitter
        # ----------------------
        if xjitter:

            if yjitter_type == "discrete" or (
                yjitter_type=="auto" and len(np.unique(yvals)) < len(yvals)*0.5
            ):

                # discrete counts
                unique_vals, counts = np.unique(yvals, return_counts=True)
                max_count = counts.max()

                for val, count in zip(unique_vals, counts):

                    idx = yvals == val
                    scale = xjitter * (count / max_count)

                    jitter = stats.t(loc=0,
                                     df=30,
                                     scale=scale).rvs(len(xvals[idx]),
                                                      random_state=np.random.Generator(np.random.PCG64(seed)))

                    xvals[idx] += jitter

            else:

                # continuous KDE scaling
                kde = stats.gaussian_kde(yvals)
                dens = kde(yvals)
                dens_scaled = dens / dens.max()

                jitter = stats.t(loc=0,
                                     df=30,
                                     scale=xjitter*dens_scaled).rvs(len(xvals),
                                                                    random_state=np.random.Generator(np.random.PCG64(seed)))

                xvals += jitter

        ax.scatter(xvals, yvals, s=size_pt, facecolors='none',
                   edgecolors=palette[ii], alpha=alpha_pt,
                   linewidths=0.25)

        y_data.append(yvals)
        x_positions.append(ii)

    # ------------
    # violin layer
    # ------------

    violins = ax.violinplot(y_data, positions=x_positions,
                            widths=0.45, bw_method='scott',
                            showmeans=False, showmedians=False,
                            showextrema=False)

    for ii, pc in enumerate(violins['bodies']):

        pc.set_facecolor(palette[ii])
        pc.set_edgecolor([0.25, 0.25, 0.25])
        pc.set_linewidth(1)
        pc.set_alpha(0.2)

    # ---------------
    # boxplot overlay
    # ---------------

    medianprops = dict(linewidth=2, color=[0.2, 0.2, 0.2],
                       solid_capstyle='butt')

    boxprops = dict(linewidth=0.75, color=[0.2, 0.2, 0.2])

    ax.boxplot(y_data, positions=x_positions, showfliers=False,
               showcaps=False, medianprops=medianprops,
               whiskerprops=boxprops, boxprops=boxprops,
               widths=0.25)

    ax.set(xticks=x_positions, xticklabels=xCats)

    return fig, ax


# %% violinsplit
def violinsplit(data, x, y, hue, xCats=None, xjitter=0.05, yjitter=0,
                 yjitter_type='auto', seed=123, med_trace=False,
                 palette=None, violin_width=0.5, figsize=(7, 4.65),
                 size_pt=10, alpha_pt=0.25, ax=None):

    # remove any rows with missing data in the relevant columns
    data = data.dropna(subset=[x, y, hue])
    
    rng = np.random.default_rng(seed)
    
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    if palette is None:
        # assign default palette
        palette = list(mpl.colors.TABLEAU_COLORS.values())

    if len(palette) > 2:
        palette = palette[0:2]

    # violin layer
    sns.violinplot(data=data, x=x, y=y, hue=hue, split=True,
                   inner='quart', cut=0, width=violin_width,
                   palette=palette, ax=ax)

    plt.setp(ax.collections, alpha=0.2)

    if xCats is None:
        xCats = list(data[x].sort_values().unique())

    hueCats = list(data[hue].sort_values().unique())

    shift = 0.055

    medians_A = []
    medians_B = []
    x_positions = []

    for ii, xcat in enumerate(xCats):

        for jj, hcat in enumerate(hueCats):

            sub = data[(data[x]==xcat) & (data[hue]==hcat)][y]

            if sub.empty:
                continue

            yvals = sub.to_numpy().astype(float)

            # optional y jitter
            if yjitter:
                yvals = yvals + rng.normal(0, yjitter, len(yvals))

            # base x location at split line
            xvals = np.full(len(yvals), ii, dtype=float)

            # determine side
            side_sign = -1 if jj == 0 else 1
            color = palette[jj]

            # ----------------------
            # density-aware x jitter
            # ----------------------

            if xjitter:

                if yjitter_type == 'discrete' or (
                    yjitter_type == 'auto' and len(np.unique(yvals)) < len(yvals)*0.5
                ):

                    # discrete counts
                    unique_vals, counts = np.unique(yvals, return_counts=True)
                    max_count = counts.max()

                    for val, count in zip(unique_vals, counts):

                        idx = yvals == val
                        scale = xjitter * (count / max_count)

                        jitter = np.abs(stats.t(loc=0,
                                                df=30,
                                                scale=scale).rvs(len(xvals[idx]),
                                                                 random_state=np.random.Generator(np.random.PCG64(seed))))
                        xvals[idx] += side_sign * jitter

                else:

                    # continuous: jitter proportional to KDE density
                    kde = stats.gaussian_kde(yvals)
                    dens = kde(yvals)
                    dens_scaled = dens / dens.max()

                    jitter = np.abs(stats.t(loc=0,
                                            df=30,
                                            scale=xjitter*dens_scaled).rvs(len(xvals),
                                                                           random_state=np.random.Generator(np.random.PCG64(seed))))

                    xvals += side_sign * jitter
            
            else:
                xvals += side_sign * 0.02

            # prevent crossing split
            if jj == 0:
                xvals = np.minimum(xvals, ii)
                color = palette[0]
            else:
                xvals = np.maximum(xvals, ii)
                color = palette[1]

            ax.scatter(xvals, yvals, s=size_pt,
                       facecolors='none', edgecolors=color,
                       alpha=alpha_pt, linewidths=0.25)

            # store medians
            if jj == 0:
                medians_A.append(np.median(yvals))
            else:
                medians_B.append(np.median(yvals))

        x_positions.append(ii)

    # ------------
    # median trace
    # ------------

    if med_trace:

        if len(medians_A) == len(x_positions):
            ax.plot(
                np.array(x_positions) - 1.1*shift,
                medians_A,
                ":",
                color=palette[0],
                linewidth=2,
                alpha=0.4
            )

        if len(medians_B) == len(x_positions):
            ax.plot(
                np.array(x_positions) + 1.1*shift,
                medians_B,
                ":",
                color=palette[1],
                linewidth=2,
                alpha=0.4
            )

    ax.legend(bbox_to_anchor=(0.5, 1.2), loc='upper center', ncol=2)

    return fig, ax


# %% violin_facet
def violin_facet(data, plot_func, facet_col=None, facet_row=None,
                 col_order=None, row_order=None, wrap_cols=None,
                 sharex=True, sharey=True, figsize=(10, 6),
                 panel_titles=True, title_fmt='{row} | {col}',
                 legend=True, legend_loc='upper center',
                 legend_ncol=None, legend_bbox_to_anchor=None,
                 xlabel=None, ylabel=None, xticks=None, yticks=None,
                 **plot_kwargs):

    if facet_col is None and facet_row is None:
        raise ValueError("At least one of facet_col or facet_row must be specified")

    # consistent x category ordering
    if "x" in plot_kwargs:
        xCats_global = list(data[plot_kwargs['x']].sort_values().unique())
    else:
        xCats_global = None

    legend_handles = None
    legend_labels = None

    # =====================================================
    # WRAP MODE (single variable arranged into grid)
    # =====================================================

    if facet_row is None and wrap_cols is not None:

        facetCats = col_order if col_order else list(data[facet_col].sort_values().unique())

        ncols = wrap_cols
        nrows = int(np.ceil(len(facetCats) / wrap_cols))

        fig, axes = plt.subplots(nrows, ncols, figsize=figsize,
                                 sharex=sharex, sharey=sharey)

        axes = np.array(axes).reshape(-1)

        for ii, fcat in enumerate(facetCats):

            ax = axes[ii]

            sub = data[data[facet_col] == fcat]

            local_kwargs = plot_kwargs.copy()
            local_kwargs["ax"] = ax

            if xCats_global is not None:
                local_kwargs["xCats"] = xCats_global

            plot_func(sub, **local_kwargs)

            if legend and legend_handles is None:
                legend_handles, legend_labels = ax.get_legend_handles_labels()

            if ax.get_legend() is not None:
                ax.legend().remove()

            if panel_titles:
                ax.set_title(str(fcat))

            if xticks is not None:
                ax.set_xticks(xticks)

            if yticks is not None:
                ax.set_yticks(yticks)

            ax.set_xlabel("")
            ax.set_ylabel("")

        # hide unused axes
        for jj in range(ii + 1, len(axes)):
            axes[jj].set_visible(False)

    # =====================================================
    # GRID MODE (row × column facets)
    # =====================================================

    else:

        if facet_col:
            colCats = col_order if col_order else list(data[facet_col].sort_values().unique())
        else:
            colCats = [None]

        if facet_row:
            rowCats = row_order if row_order else list(data[facet_row].sort_values().unique())
        else:
            rowCats = [None]

        ncols = len(colCats)
        nrows = len(rowCats)

        fig, axes = plt.subplots(nrows, ncols, figsize=figsize,
                                 sharex=sharex, sharey=sharey)

        axes = np.array(axes).reshape(nrows, ncols)

        for ii, rcat in enumerate(rowCats):
            for jj, ccat in enumerate(colCats):

                ax = axes[ii, jj]

                sub = data.copy()

                if rcat is not None:
                    sub = sub[sub[facet_row] == rcat]

                if ccat is not None:
                    sub = sub[sub[facet_col] == ccat]

                local_kwargs = plot_kwargs.copy()
                local_kwargs['ax'] = ax

                if xCats_global is not None:
                    local_kwargs['xCats'] = xCats_global

                plot_func(sub, **local_kwargs)

                if legend and legend_handles is None:
                    legend_handles, legend_labels = ax.get_legend_handles_labels()

                if ax.get_legend() is not None:
                    ax.legend().remove()

                if panel_titles:

                    if facet_row and facet_col:
                        title = title_fmt.format(row=rcat, col=ccat)
                    elif facet_col:
                        title = str(ccat)
                    else:
                        title = str(rcat)

                    ax.set_title(title)

                if xticks is not None:
                    ax.set_xticks(xticks)

                if yticks is not None:
                    ax.set_yticks(yticks)

                ax.set_xlabel("")
                ax.set_ylabel("")

    # =====================================================
    # Shared legend
    # =====================================================

    if legend and legend_handles:

        if legend_ncol is None:
            legend_ncol = len(legend_labels)

        fig.legend(legend_handles, legend_labels,
                   loc=legend_loc, ncol=legend_ncol,
                   bbox_to_anchor=legend_bbox_to_anchor,
                   frameon=False)

    # =====================================================
    # Global axis labels
    # =====================================================

    if xlabel:
        fig.supxlabel(xlabel)

    if ylabel:
        fig.supylabel(ylabel)

    fig.tight_layout()

    return fig, axes
