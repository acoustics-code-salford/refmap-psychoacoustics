"""
Paired (within-subjects) cluster-bootstrap confidence intervals for the SPI.

Reference = ``target`` ratings (e.g. baseline), comparison = ``test`` ratings
(e.g. the stimulus, or pooled stimuli, with an added source). The same
participants rated both, so each bootstrap replicate resamples PARTICIPANTS with
replacement and takes, for each chosen participant, all of their target rows AND
all of their test rows. This preserves both the within-participant dependence
and the pairing between the two sets. The two sets may have different numbers
of rows (e.g. one target rating per participant, k test ratings each).

The statistic is SPI = (1 - D) * 100, with D the two-sample 2-D KS statistic
from ``soundscapy.spi.msn.ks2d(target, test)``, i.e. the same quantity as
``spi_score(target, test)`` but without integer truncation.

The interval uses Hesterberg's expanded percentile levels, with df = G - 1 or a
Satterthwaite-style effective df from leave-one-participant-out influences.

Speed
-----
Each replicate needs one D. If soundscapy's ``ks2d`` is slow (it also computes a
p-value, which the bootstrap does not need), pass ``ks2d_func=fast_ks2d``: a
vectorised D-only reimplementation of soundscapy's ``ks2d2s`` statistic (same
strict-inequality quadrant convention). Check it reproduces your installed
``ks2d`` on your own data first:

>>> from soundscapy.spi.msn import ks2d
>>> from spi_bootstrap import fast_ks2d
>>> print(ks2d(target_arr, test_arr)[0], fast_ks2d(target_arr, test_arr)[0])

Usage
-----
>>> from spi_bootstrap import spi_paired_bootstrap, fast_ks2d
>>> cols = ["ISOPleasant", "ISOEventful"]
>>> res = spi_paired_bootstrap(
...     target_df[cols], test_df[cols],
...     target_groups=target_df["participant"], test_groups=test_df["participant"],
...     n_boot=10_000, df="satterthwaite", seed=1, ks2d_func=fast_ks2d)
>>> print(res)
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from numbers import Real
from typing import Callable, Literal

import numpy as np
import pandas as pd
from scipy import stats


@dataclass
class SPIBootResult:
    """Result of :func:`spi_paired_bootstrap`."""

    spi: float  # continuous SPI = (1 - D) * 100 on the full data
    spi_int: int  # truncated SPI, as soundscapy.spi.msn.spi_score would return
    d: float  # KS2D statistic on the full data
    ci_low: float
    ci_high: float
    conf: float
    n_participants: int  # G: participants with data in both sets
    n_target: int  # target rows used
    n_test: int  # test rows used
    df: float  # df used in the expanded-percentile t quantile
    df_method: str
    percentile_levels: tuple[float, float]  # percentiles actually used (0-100)
    n_boot: int
    boot_spi: np.ndarray  # bootstrap distribution of SPI

    def __str__(self) -> str:
        lo, hi = self.percentile_levels
        n_eff = self.n_target * self.n_test / (self.n_target + self.n_test)
        return (
            f"SPI = {self.spi:.1f} (int: {self.spi_int}), "
            f"{self.conf:.0%} CI [{self.ci_low:.1f}, {self.ci_high:.1f}]\n"
            f"G = {self.n_participants} participants, rows: target = "
            f"{self.n_target}, test = {self.n_test} (n_eff = {n_eff:.1f})\n"
            f"df = {self.df:.2f} ({self.df_method}), "
            f"percentiles used = ({lo:.2f}, {hi:.2f}), B = {self.n_boot}"
        )


def fast_ks2d(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """
    Vectorised D from ``soundscapy.spi.ks2d.ks2d2s``, reproducing its convention.

    For every point of each sample taken as an origin, compare the fraction of
    each sample falling in the four quadrants, take the largest absolute
    difference, and average the maxima over the two origin sets (d1 + d2) / 2.
    As in soundscapy's ``CountQuads``, all four quadrant tests are STRICT
    inequalities, so a point sharing an x or y value with the origin (including
    the origin itself) falls in no quadrant. Returns ``(D, nan)``: the p-value is
    not computed, because the bootstrap only needs D.

    Memory is O((n + m)^2) booleans, fine up to a few thousand points.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    n, m = len(a), len(b)
    pts = np.vstack([a, b])
    gx = pts[None, :, 0] > pts[:, None, 0]  # [i, j]: x_j > x_i
    lx = pts[None, :, 0] < pts[:, None, 0]
    gy = pts[None, :, 1] > pts[:, None, 1]
    ly = pts[None, :, 1] < pts[:, None, 1]

    best = np.zeros(len(pts))  # per-origin max over the four quadrants
    for qx, qy in ((gx, gy), (lx, gy), (gx, ly), (lx, ly)):
        q = qx & qy
        diff = np.abs(q[:, :n].sum(1) / n - q[:, n:].sum(1) / m)
        np.maximum(best, diff, out=best)
    return float((best[:n].max() + best[n:].max()) / 2), float("nan")


def expanded_percentile_levels(
    n_clusters: int, df: float | None = None, conf: float = 0.95
) -> tuple[float, float]:
    """
    Hesterberg's expanded percentile levels (as percentages, 0-100).

    alpha' = Phi(-sqrt(G / (G - 1)) * t_{alpha/2, df}); the interval is the
    alpha' and 1 - alpha' percentiles of the bootstrap distribution.
    The sqrt(G/(G-1)) factor always uses the number of clusters G; only the
    t quantile uses ``df`` (default G - 1).
    """
    if n_clusters < 2:  # noqa: PLR2004
        msg = "Need at least 2 clusters."
        raise ValueError(msg)
    df = n_clusters - 1 if df is None else df
    alpha = 1 - conf
    t = stats.t.ppf(1 - alpha / 2, df)
    a = float(stats.norm.cdf(-np.sqrt(n_clusters / (n_clusters - 1)) * t))
    return 100 * a, 100 * (1 - a)


def _to_array(x: pd.DataFrame | np.ndarray, label: str) -> np.ndarray:
    arr = x.to_numpy() if isinstance(x, pd.DataFrame) else np.asarray(x)
    if arr.ndim != 2 or arr.shape[1] != 2:  # noqa: PLR2004
        msg = f"{label} must have shape (n, 2)."
        raise ValueError(msg)
    arr = arr.astype(float)
    if not np.isfinite(arr).all():
        msg = (
            f"{label} contains NaN/inf values. Drop those rows (and any matching "
            "participant ids) before calling."
        )
        raise ValueError(msg)
    return arr


def _by_group(arr: np.ndarray, groups: object, label: str) -> dict:
    groups = np.asarray(groups)
    if len(groups) != len(arr):
        msg = f"{label}_groups must have one entry per row of {label}."
        raise ValueError(msg)
    if pd.isna(groups).any():
        msg = f"{label}_groups contains missing values."
        raise ValueError(msg)
    return {k: arr[groups == k] for k in pd.unique(groups)}


def _satterthwaite_df(
    target_cl: list[np.ndarray],
    test_cl: list[np.ndarray],
    d_func: Callable,
    spi_full: float,
) -> tuple[float, bool]:
    """
    Effective df from leave-one-participant-out influence values d_g.

    df ~= 3 (sum d_g^2)^2 / sum d_g^4, capped at G - 1. Equals ~G when all
    participants contribute equally and are roughly normal. Noisy for small G
    (fourth moments from G points): treat as a diagnostic.
    Returns (df, was_capped).
    """
    g = len(target_cl)
    infl = np.empty(g)
    for i in range(g):
        a = np.vstack([c for j, c in enumerate(target_cl) if j != i])
        t = np.vstack([c for j, c in enumerate(test_cl) if j != i])
        infl[i] = spi_full - 100 * (1 - d_func(a, t))
    s2 = np.sum(infl**2)
    s4 = np.sum(infl**4)
    if s4 == 0:
        return g - 1.0, True
    raw = 3 * s2**2 / s4
    return float(min(raw, g - 1)), bool(raw >= g - 1)


def spi_paired_bootstrap(
    target: pd.DataFrame | np.ndarray,
    test: pd.DataFrame | np.ndarray,
    target_groups: pd.Series | np.ndarray | list,
    test_groups: pd.Series | np.ndarray | list,
    *,
    n_boot: int = 10_000,
    conf: float = 0.95,
    expanded: bool = True,
    df: Literal["G-1", "satterthwaite"] | float = "G-1",
    seed: int | None = None,
    ks2d_func: Callable | None = None,
) -> SPIBootResult:
    """
    Participant-level (paired) bootstrap CI for the SPI of ``test`` vs ``target``.

    Parameters
    ----------
    target
        Reference ratings, shape (n_t, 2), one row per rating. Plays the role of
        ``target`` in ``spi_score(target, test)``.
    test
        Comparison ratings, shape (n_s, 2). May contain several ratings per
        participant (e.g. pooled over stimuli); n_s need not equal n_t.
    target_groups, test_groups
        Participant id for each row of ``target`` / ``test``. Ids must match
        between the two (same type and values).
    n_boot
        Bootstrap replicates. Expanded percentiles sit in the tails for small
        G, so use >= 10,000 if you can, or >= 2,000-5,000 for G of 30-50.
    conf
        Confidence level.
    expanded
        If True use Hesterberg's expanded percentile levels; otherwise plain
        percentiles.
    df
        "G-1" (default), "satterthwaite" (leave-one-participant-out estimate),
        or a number. Only used when ``expanded`` is True.
    seed
        Seed for the bootstrap.
    ks2d_func
        Callable (a, b) -> (D, p). Defaults to ``soundscapy.spi.msn.ks2d``;
        pass :func:`fast_ks2d` for a much faster D-only version.

    Notes
    -----
    Participants present in only one of the two sets are dropped (with a
    warning), since the pairing is what the resampling relies on.
    D depends on both sample sizes (roughly through n_t * n_s / (n_t + n_s)),
    so SPI is biased low when either set is small, and comparisons of SPI between
    contrasts with different n_s are only approximately like-for-like. The
    bootstrap interval reflects participant sampling variability only, not
    variability across stimuli.
    """
    if ks2d_func is None:
        from soundscapy.spi.msn import ks2d as ks2d_func  # noqa: PLC0415

    def d_func(a: np.ndarray, b: np.ndarray) -> float:
        return float(ks2d_func(a, b)[0])

    rng = np.random.default_rng(seed)

    target_by = _by_group(_to_array(target, "target"), target_groups, "target")
    test_by = _by_group(_to_array(test, "test"), test_groups, "test")

    ids = [k for k in target_by if k in test_by]
    if not ids:
        msg = (
            "No participant ids are shared between target and test "
            "(check the id values and dtypes match)."
        )
        raise ValueError(msg)
    dropped = set(target_by) ^ set(test_by)
    if dropped:
        warnings.warn(
            f"{len(dropped)} participant id(s) appear in only one of target/test "
            f"and were dropped: {sorted(map(str, dropped))[:10]}",
            UserWarning,
            stacklevel=2,
        )

    target_cl = [target_by[k] for k in ids]
    test_cl = [test_by[k] for k in ids]
    g = len(ids)
    if g < 2:  # noqa: PLR2004
        msg = "Need at least 2 participants to bootstrap."
        raise ValueError(msg)
    if g < 10:  # noqa: PLR2004
        warnings.warn(
            f"Only {g} participants: bootstrap interval will be approximate.",
            UserWarning,
            stacklevel=2,
        )

    # --- point estimate ---
    target_all = np.vstack(target_cl)
    test_all = np.vstack(test_cl)
    d_full = d_func(target_all, test_all)
    spi_full = 100 * (1 - d_full)
    spi_int = int((1 - d_full) * 100)  # same truncation as spi_score

    # --- df ---
    if not expanded:
        df_val, df_method = float("nan"), "n/a (plain percentile)"
    elif isinstance(df, str) and df == "G-1":
        df_val, df_method = g - 1.0, "G-1"
    elif isinstance(df, str) and df == "satterthwaite":
        df_val, capped = _satterthwaite_df(target_cl, test_cl, d_func, spi_full)
        df_method = "Satterthwaite, capped at G-1" if capped else "Satterthwaite"
    elif isinstance(df, Real) and not isinstance(df, bool):
        df_val, df_method = float(df), "user-supplied"
    else:
        msg = "df must be 'G-1', 'satterthwaite' or a number."
        raise ValueError(msg)

    # --- percentile levels ---
    if expanded:
        lo, hi = expanded_percentile_levels(g, df_val, conf)
    else:
        lo, hi = 100 * (1 - conf) / 2, 100 * (1 + conf) / 2
    if n_boot * min(lo, 100 - hi) / 100 < 20:  # noqa: PLR2004
        warnings.warn(
            f"Percentile levels ({lo:.2f}, {hi:.2f}) with n_boot={n_boot} leave "
            "fewer than ~20 replicates in each tail; increase n_boot.",
            UserWarning,
            stacklevel=2,
        )

    # --- bootstrap: same participants drawn for both sets ---
    boot = np.empty(n_boot)
    for b in range(n_boot):
        pick = rng.integers(0, g, size=g)
        a = np.vstack([target_cl[i] for i in pick])
        t = np.vstack([test_cl[i] for i in pick])
        boot[b] = 100 * (1 - d_func(a, t))

    ci_low, ci_high = np.percentile(boot, [lo, hi])

    return SPIBootResult(
        spi=spi_full,
        spi_int=spi_int,
        d=d_full,
        ci_low=float(ci_low),
        ci_high=float(ci_high),
        conf=conf,
        n_participants=g,
        n_target=len(target_all),
        n_test=len(test_all),
        df=df_val,
        df_method=df_method,
        percentile_levels=(lo, hi),
        n_boot=n_boot,
        boot_spi=boot,
    )


# =============================================================================
# Chance-corrected (bias-adjusted) SPI
# =============================================================================
#
# D has a noise floor: two samples from the SAME distribution give D > 0, and
# how far above 0 depends on the sample sizes and the clustering. The adjusted
# SPI rescales D against its own no-difference expectation D0, as Cohen's kappa
# does for chance agreement:
#
#     SPI_adj = 100 * (1 - (D - D0) / (1 - D0))
#
# D0 is the mean of D over within-participant relabellings: for each
# participant, their target and test ratings are pooled and the same NUMBER of
# ratings is randomly re-assigned as "target". This keeps each participant's
# ratings together, so D0 reflects the actual design (sizes, clustering).
# No difference -> SPI_adj ~ 100 (it can exceed 100 by chance; it is not
# clipped, so it stays centred). Inside the bootstrap, D0 is recomputed on every
# resample, so the extra noise from duplicated participants is corrected too.


def _participant_pools(
    target_cl: list[np.ndarray], test_cl: list[np.ndarray]
) -> tuple[list[np.ndarray], list[int]]:
    pools = [np.vstack([a, t]) for a, t in zip(target_cl, test_cl)]
    return pools, [len(a) for a in target_cl]


def _relabel(
    pools: list[np.ndarray], n_t: list[int], idx: np.ndarray, rng: np.random.Generator
) -> tuple[np.ndarray, np.ndarray]:
    """
    Within-participant relabelling for participants ``idx`` (may repeat).

    A participant drawn more than once in a bootstrap resample gets the SAME
    relabelling for every copy, mirroring the observed resample, where the
    copies are identical (labels included).
    """
    perms: dict = {}
    a_parts, t_parts = [], []
    for i in idx:
        if i not in perms:
            perms[i] = rng.permutation(len(pools[i]))
        p, perm = pools[i], perms[i]
        a_parts.append(p[perm[: n_t[i]]])
        t_parts.append(p[perm[n_t[i] :]])
    return np.vstack(a_parts), np.vstack(t_parts)


def _d_and_null(
    target_cl: list[np.ndarray],
    test_cl: list[np.ndarray],
    pools: list[np.ndarray],
    n_t: list[int],
    idx: np.ndarray,
    n_perm: int,
    d_func: Callable,
    rng: np.random.Generator,
) -> tuple[float, np.ndarray]:
    """Observed D for participants ``idx``, plus ``n_perm`` relabelled-null Ds."""
    a = np.vstack([target_cl[i] for i in idx])
    t = np.vstack([test_cl[i] for i in idx])
    d = d_func(a, t)
    d0 = np.array([d_func(*_relabel(pools, n_t, idx, rng)) for _ in range(n_perm)])
    return d, d0


def adjust_spi(d: float | np.ndarray, d0: float) -> float | np.ndarray:
    """Linear chance-corrected SPI: 100 * (1 - (D - D0) / (1 - D0))."""
    return 100 * (1 - (d - d0) / (1 - d0))


def adjust_spi_quad(d: float | np.ndarray, d0_sq: float) -> float | np.ndarray:
    """
    Signed-quadrature corrected SPI: 100 * (1 - sign(v) * sqrt(|v|)),
    v = D^2 - E0[D^2]. Treats the noise as adding to D in quadrature
    (a heuristic). Not clipped, so it can exceed 100.
    """
    v = np.asarray(d) ** 2 - d0_sq
    out = 100 * (1 - np.sign(v) * np.sqrt(np.abs(v)))
    return float(out) if np.ndim(out) == 0 else out


@dataclass
class SPIAdjResult:
    """Result of :func:`spi_adjusted_bootstrap`."""

    spi_raw: float  # 100 (1 - D), as spi_score (continuous)
    spi_adj: float  # linear chance-corrected SPI
    spi_quad: float  # signed-quadrature corrected SPI
    d: float
    d0: float  # mean D under the within-participant null
    ci_raw: tuple[float, float]
    ci_adj: tuple[float, float]
    ci_quad: tuple[float, float]
    p_value: float  # permutation p-value for "any difference" (D >= observed)
    null_spi_range: tuple[float, float]  # 5-95% of raw SPI under the null
    null_adj_range: tuple[float, float]  # 5-95% of linear adjusted SPI, null
    null_quad_range: tuple[float, float]  # 5-95% of quadrature SPI, null
    conf: float
    n_participants: int
    n_target: int
    n_test: int
    df: float
    df_method: str
    percentile_levels: tuple[float, float]
    n_boot: int
    n_perm: int
    n_perm_boot: int
    boot_raw: np.ndarray
    boot_adj: np.ndarray
    boot_quad: np.ndarray
    boot_d0: np.ndarray  # per-replicate null mean D

    def __str__(self) -> str:
        lo, hi = self.percentile_levels
        return (
            f"G = {self.n_participants} participants; rows: target = {self.n_target}, "
            f"test = {self.n_test}\n"
            f"Raw SPI      {self.spi_raw:6.1f}   {self.conf:.0%} CI "
            f"[{self.ci_raw[0]:6.1f}, {self.ci_raw[1]:6.1f}]   "
            f"(no-difference range 5-95%: {self.null_spi_range[0]:.1f}-"
            f"{self.null_spi_range[1]:.1f})\n"
            f"Linear adj.  {self.spi_adj:6.1f}   {self.conf:.0%} CI "
            f"[{self.ci_adj[0]:6.1f}, {self.ci_adj[1]:6.1f}]   "
            f"(no-difference range 5-95%: {self.null_adj_range[0]:.1f}-"
            f"{self.null_adj_range[1]:.1f})\n"
            f"Quadrature   {self.spi_quad:6.1f}   {self.conf:.0%} CI "
            f"[{self.ci_quad[0]:6.1f}, {self.ci_quad[1]:6.1f}]   "
            f"(no-difference range 5-95%: {self.null_quad_range[0]:.1f}-"
            f"{self.null_quad_range[1]:.1f})\n"
            f"D = {self.d:.3f}, null mean D0 = {self.d0:.3f}; permutation p = "
            f"{self.p_value:.4f} ({self.n_perm} relabellings)\n"
            f"df = {self.df:.2f} ({self.df_method}), percentiles used = "
            f"({lo:.2f}, {hi:.2f}), B = {self.n_boot}, null draws per replicate = "
            f"{self.n_perm_boot}"
        )


def spi_adjusted_bootstrap(
    target: pd.DataFrame | np.ndarray,
    test: pd.DataFrame | np.ndarray,
    target_groups: pd.Series | np.ndarray | list,
    test_groups: pd.Series | np.ndarray | list,
    *,
    n_boot: int = 2000,
    n_perm: int = 2000,
    n_perm_boot: int = 20,
    conf: float = 0.95,
    expanded: bool = True,
    df: Literal["G-1", "satterthwaite"] | float = "G-1",
    seed: int | None = None,
    ks2d_func: Callable | None = None,
    progress: bool = False,
) -> SPIAdjResult:
    """
    Chance-corrected SPI with a paired participant bootstrap CI.

    Parameters
    ----------
    target, test, target_groups, test_groups
        As for :func:`spi_paired_bootstrap`.
    n_boot
        Bootstrap replicates (participants resampled with replacement).
    n_perm
        Within-participant relabellings on the full data. Gives D0 for the
        point estimate, the permutation p-value and the no-difference ranges.
    n_perm_boot
        Relabellings per bootstrap replicate, used to recompute D0 on each
        resample. Total cost is about n_boot * (1 + n_perm_boot) D evaluations.
    conf, expanded, df, seed
        As for :func:`spi_paired_bootstrap`. ``df="satterthwaite"`` uses
        leave-one-participant-out estimates, each with its own D0.
    ks2d_func
        Defaults to :func:`fast_ks2d` (same statistic as soundscapy's ks2d).

    Notes
    -----
    * The null treats each participant's baseline rating(s) as exchangeable
      with ALL their test ratings. If ``test`` pools several stimuli that
      differ from each other, the null is "baseline and every pooled stimulus
      are identical", which is strict.
    * Validation (validate_adjusted_spi.py; 15 simulated cases in this design,
      G = 43, 1 + 12 ratings per participant, 30 datasets each):
      - linear: centred at "no difference", but overstates similarity by about
        4-10 points whenever there is a real difference, and its CI
        under-covers (0.37-0.90; worst for spread changes). Not recommended as
        a measure of magnitude.
      - quadrature: within about 4 points of the population SPI when that is
        <= ~83, with conservative CIs (coverage 0.87-1.00), but very wide near
        100 (about 40 points) and biased upward by 2-11 points when the true
        value is 88-100. Differences in that top band cannot be resolved with
        this design.
    * The p-value is for "no difference at all", not a test of similarity.
    """
    ks2d_func = ks2d_func or fast_ks2d

    def d_func(a: np.ndarray, b: np.ndarray) -> float:
        return float(ks2d_func(a, b)[0])

    rng = np.random.default_rng(seed)
    target_cl, test_cl = _paired_clusters(target, test, target_groups, test_groups)
    g = len(target_cl)
    pools, n_t = _participant_pools(target_cl, test_cl)
    all_idx = np.arange(g)

    # --- point estimates and full-data null ---
    d_obs, d0_full = _d_and_null(
        target_cl, test_cl, pools, n_t, all_idx, n_perm, d_func, rng
    )
    d0m = float(d0_full.mean())
    spi_raw = 100 * (1 - d_obs)
    spi_adj = float(adjust_spi(d_obs, d0m))
    d0sq = float(np.mean(d0_full**2))
    spi_quad = adjust_spi_quad(d_obs, d0sq)
    p_value = float((1 + np.sum(d0_full >= d_obs)) / (n_perm + 1))
    null_raw = 100 * (1 - d0_full)
    null_adj = adjust_spi(d0_full, d0m)
    null_quad = adjust_spi_quad(d0_full, d0sq)

    # --- df ---
    if not expanded:
        df_val, df_method = float("nan"), "n/a (plain percentile)"
    elif isinstance(df, str) and df == "G-1":
        df_val, df_method = g - 1.0, "G-1"
    elif isinstance(df, str) and df == "satterthwaite":
        infl = np.empty(g)
        for i in range(g):
            idx = np.delete(all_idx, i)
            d_i, d0_i = _d_and_null(
                target_cl, test_cl, pools, n_t, idx, n_perm_boot, d_func, rng
            )
            infl[i] = spi_adj - adjust_spi(d_i, d0_i.mean())
        df_val = _df_from_influence(infl, g)
        df_method = (
            "Satterthwaite, capped at G-1" if df_val >= g - 1 else "Satterthwaite"
        )
    elif isinstance(df, Real) and not isinstance(df, bool):
        df_val, df_method = float(df), "user-supplied"
    else:
        msg = "df must be 'G-1', 'satterthwaite' or a number."
        raise ValueError(msg)

    if expanded:
        lo, hi = expanded_percentile_levels(g, df_val, conf)
    else:
        lo, hi = 100 * (1 - conf) / 2, 100 * (1 + conf) / 2
    if n_boot * min(lo, 100 - hi) / 100 < 20:  # noqa: PLR2004
        warnings.warn(
            f"Percentile levels ({lo:.2f}, {hi:.2f}) with n_boot={n_boot} leave "
            "fewer than ~20 replicates in each tail; increase n_boot.",
            UserWarning,
            stacklevel=2,
        )

    # --- bootstrap, with D0 recomputed on every resample ---
    iterator: object = range(n_boot)
    if progress:
        try:
            from tqdm import tqdm  # noqa: PLC0415

            iterator = tqdm(iterator, desc="bootstrap")
        except ImportError:
            pass
    boot_raw, boot_adj, boot_quad, boot_d0 = (np.empty(n_boot) for _ in range(4))
    for b in iterator:
        pick = rng.integers(0, g, size=g)
        d, d0 = _d_and_null(
            target_cl, test_cl, pools, n_t, pick, n_perm_boot, d_func, rng
        )
        boot_d0[b] = d0.mean()
        boot_raw[b] = 100 * (1 - d)
        boot_adj[b] = adjust_spi(d, boot_d0[b])
        boot_quad[b] = adjust_spi_quad(d, np.mean(d0**2))

    return SPIAdjResult(
        spi_raw=spi_raw,
        spi_adj=spi_adj,
        spi_quad=spi_quad,
        d=d_obs,
        d0=d0m,
        ci_raw=tuple(float(v) for v in np.percentile(boot_raw, [lo, hi])),
        ci_adj=tuple(float(v) for v in np.percentile(boot_adj, [lo, hi])),
        ci_quad=tuple(float(v) for v in np.percentile(boot_quad, [lo, hi])),
        p_value=p_value,
        null_spi_range=tuple(float(v) for v in np.percentile(null_raw, [5, 95])),
        null_adj_range=tuple(float(v) for v in np.percentile(null_adj, [5, 95])),
        null_quad_range=tuple(float(v) for v in np.percentile(null_quad, [5, 95])),
        conf=conf,
        n_participants=g,
        n_target=sum(n_t),
        n_test=sum(len(t) for t in test_cl),
        df=df_val,
        df_method=df_method,
        percentile_levels=(lo, hi),
        n_boot=n_boot,
        n_perm=n_perm,
        n_perm_boot=n_perm_boot,
        boot_raw=boot_raw,
        boot_adj=boot_adj,
        boot_quad=boot_quad,
        boot_d0=boot_d0,
    )


# =============================================================================
# Empirical SPI vs MSN-smoothed SPI, with a paired participant bootstrap
# =============================================================================
#
# Published SPI workflow: fit a multivariate skew-normal (MSN) to the test data,
# sample it, and compare that sample with a sample from the target MSN. The
# "empirical" SPI above skips the MSN step. ``spi_compare_bootstrap`` computes
# both on the SAME resampled participants, so the two are directly comparable
# and their difference gets its own CI.

_FIT_ERRORS = (ValueError, RuntimeError, ArithmeticError)


def msn_sample_from_params(
    xi: np.ndarray, omega: np.ndarray, alpha: np.ndarray, z: np.ndarray
) -> np.ndarray:
    """
    Transform fixed standard-normal draws ``z`` (n, k+1) into MSN samples.

    Azzalini's stochastic representation: with delta = Obar alpha /
    sqrt(1 + alpha' Obar alpha) and (U0, U) ~ N(0, [[1, delta'], [delta, Obar]]),
    Z = U if U0 > 0 else -U is standard skew-normal, and Y = xi + diag(omega)^.5 Z.
    Because ``z`` is passed in, the SAME underlying draws can be reused for every
    bootstrap replicate (common random numbers), so Monte Carlo noise does not
    add to bootstrap noise and results are reproducible.
    """
    xi = np.asarray(xi, dtype=float).ravel()
    alpha = np.asarray(alpha, dtype=float).ravel()
    omega = np.asarray(omega, dtype=float)
    k = len(xi)
    d = np.sqrt(np.diag(omega))
    obar = omega / np.outer(d, d)
    delta = obar @ alpha / np.sqrt(1.0 + alpha @ obar @ alpha)
    sigma = np.empty((k + 1, k + 1))
    sigma[0, 0] = 1.0
    sigma[0, 1:] = delta
    sigma[1:, 0] = delta
    sigma[1:, 1:] = obar
    try:
        chol = np.linalg.cholesky(sigma)
    except np.linalg.LinAlgError:
        chol = np.linalg.cholesky(sigma + 1e-10 * np.eye(k + 1))
    w = z @ chol.T
    u0, u = w[:, 0], w[:, 1:]
    zz = np.where(u0[:, None] > 0, u, -u)
    return xi + zz * d


def _fit_msn_soundscapy(arr: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fit an MSN with soundscapy (R's sn::selm) and return (xi, omega, alpha)."""
    from soundscapy.spi.msn import MultiSkewNorm  # noqa: PLC0415

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # "SPI module is experimental" etc.
        msn = MultiSkewNorm()
        try:
            msn.fit(arr)
        except Exception as e:  # rpy2's RRuntimeError need not subclass RuntimeError
            msg = f"MSN fit failed: {e}"
            raise RuntimeError(msg) from e
    return (
        np.asarray(msn.dp.xi, dtype=float).ravel(),
        np.asarray(msn.dp.omega, dtype=float),
        np.asarray(msn.dp.alpha, dtype=float).ravel(),
    )


@dataclass
class MSNSPIResult:
    """Result of :func:`msn_spi`."""

    spi: float  # MSN-route SPI, continuous
    spi_int: int  # truncated, as soundscapy's spi_score reports it
    spi_empirical: float  # 100 (1 - D) on the raw ratings, for comparison
    d: float  # KS2D between the two sets of MSN draws
    dp_target: tuple  # fitted (xi, omega, alpha)
    dp_test: tuple
    degenerate: bool  # either fit has |alpha| > alpha_flag
    n_sample: int
    seed: int | None
    shared_draws: bool

    def __float__(self) -> float:
        return self.spi

    def __str__(self) -> str:
        flag = "  (near-degenerate MSN fit: treat with caution)" if self.degenerate else ""
        return (
            f"MSN SPI = {self.spi:.1f} (int: {self.spi_int}); empirical SPI = "
            f"{self.spi_empirical:.1f}; {self.n_sample} draws per side, seed {self.seed}, "
            f"{'shared' if self.shared_draws else 'independent'} draws{flag}"
        )


def msn_spi(
    target: pd.DataFrame | np.ndarray,
    test: pd.DataFrame | np.ndarray,
    *,
    n_sample: int = 1000,
    seed: int | None = 1,
    shared_draws: bool = False,
    ks2d_func: Callable | None = None,
    fit_func: Callable | None = None,
    alpha_flag: float = 50.0,
) -> MSNSPIResult:
    """
    SPI by the published MSN route, as a single reproducible call.

    Fits an MSN to ``target`` and to ``test`` (soundscapy / R ``sn::selm`` by
    default), draws ``n_sample`` points from each fit and returns
    100 (1 - D) for the two sets of draws, plus the empirical SPI on the raw
    ratings. No bootstrap.

    With a fixed ``seed`` (default 1), the same underlying random numbers are
    used for every call, so results are reproducible and comparisons between
    environments are not affected by luck of the draw. ``seed=None`` gives
    fresh random draws each call, like soundscapy's unseeded sampling.

    Note: identical fits score about 94, not 100, at n_sample = 1000, because
    the two sets of draws differ by chance (soundscapy's spi_score has the same
    floor). Values are comparable with each other only at the same n_sample.

    ``shared_draws=True`` pushes the SAME random numbers through both fitted
    distributions. Identical fits then score exactly 100 and draw noise is much
    smaller near the top of the scale (in tests: SD 0.6-0.7 vs 1.3-1.8 points),
    so the value tracks the distance between the two fitted distributions more
    closely. It departs from the published method (independent samples), so
    values are not comparable with soundscapy's spi_score. Neither option
    removes the uncertainty from fitting the MSNs to small samples.
    """
    ks2d_func = ks2d_func or fast_ks2d
    fit_func = fit_func or _fit_msn_soundscapy
    a = _to_array(target, "target")
    t = _to_array(test, "test")
    rng = np.random.default_rng(seed)
    z_a = rng.standard_normal((n_sample, 3))
    z_t = z_a if shared_draws else rng.standard_normal((n_sample, 3))
    dp_a, dp_t = fit_func(a), fit_func(t)
    d = float(ks2d_func(msn_sample_from_params(*dp_a, z_a), msn_sample_from_params(*dp_t, z_t))[0])
    amax = max(np.abs(dp_a[2]).max(), np.abs(dp_t[2]).max())
    degenerate = bool(amax > alpha_flag)
    if degenerate:
        warnings.warn(
            f"Near-degenerate MSN fit (largest |alpha| = {amax:.3g}); the MSN SPI "
            "reflects an extreme fitted shape.",
            UserWarning,
            stacklevel=2,
        )
    return MSNSPIResult(
        spi=100 * (1 - d),
        spi_int=int((1 - d) * 100),
        spi_empirical=100 * (1 - float(ks2d_func(a, t)[0])),
        d=d,
        dp_target=dp_a,
        dp_test=dp_t,
        degenerate=degenerate,
        n_sample=n_sample,
        seed=seed,
        shared_draws=shared_draws,
    )


def _paired_clusters(
    target: object, test: object, target_groups: object, test_groups: object
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    target_by = _by_group(_to_array(target, "target"), target_groups, "target")
    test_by = _by_group(_to_array(test, "test"), test_groups, "test")
    ids = [k for k in target_by if k in test_by]
    if not ids:
        msg = (
            "No participant ids are shared between target and test "
            "(check the id values and dtypes match)."
        )
        raise ValueError(msg)
    dropped = set(target_by) ^ set(test_by)
    if dropped:
        warnings.warn(
            f"{len(dropped)} participant id(s) appear in only one of target/test "
            f"and were dropped: {sorted(map(str, dropped))[:10]}",
            UserWarning,
            stacklevel=3,
        )
    if len(ids) < 2:  # noqa: PLR2004
        msg = "Need at least 2 participants to bootstrap."
        raise ValueError(msg)
    if len(ids) < 10:  # noqa: PLR2004
        warnings.warn(
            f"Only {len(ids)} participants: bootstrap interval will be approximate.",
            UserWarning,
            stacklevel=3,
        )
    return [target_by[k] for k in ids], [test_by[k] for k in ids]


def _df_from_influence(infl: np.ndarray, g: int) -> float:
    """Satterthwaite-style df, 3 (sum d^2)^2 / sum d^4, capped at G - 1."""
    infl = infl[np.isfinite(infl)]
    if len(infl) < 3:  # noqa: PLR2004
        return g - 1.0
    s2, s4 = np.sum(infl**2), np.sum(infl**4)
    if s4 == 0:
        return g - 1.0
    return float(min(3 * s2**2 / s4, g - 1))


@dataclass
class SPICompareResult:
    """Result of :func:`spi_compare_bootstrap`. Rows: empirical, MSN, MSN - empirical."""

    names: tuple[str, ...]
    estimate: np.ndarray
    ci_low: np.ndarray
    ci_high: np.ndarray
    df: np.ndarray
    df_method: str
    percentile_levels: list[tuple[float, float]]
    bias: np.ndarray  # mean(bootstrap) - estimate
    boot: np.ndarray  # (n_boot, 3); NaN where an MSN refit failed
    n_failed: int
    n_boot: int
    n_participants: int
    n_target: int
    n_test: int
    n_sample: int
    conf: float
    dp_target: tuple  # fitted (xi, omega, alpha) on the full target data
    dp_test: tuple  # fitted (xi, omega, alpha) on the full test data
    alpha_max: np.ndarray  # per replicate: largest |alpha| of the two fits
    alpha_flag: float
    n_degenerate: int  # replicates with alpha_max > alpha_flag

    def __str__(self) -> str:
        head = (
            f"G = {self.n_participants} participants; rows: target = {self.n_target}, "
            f"test = {self.n_test}; MSN draws n = {self.n_sample} per side "
            f"(common random numbers); B = {self.n_boot}\n"
            f"{self.conf:.0%} CIs, expanded percentile, df: {self.df_method}\n"
        )
        rows = [
            f"{name:<18s}{est:7.1f}   [{lo:6.1f}, {hi:6.1f}]   df={dfv:5.1f}   "
            f"boot bias {b:+.1f}"
            for name, est, lo, hi, dfv, b in zip(
                self.names, self.estimate, self.ci_low, self.ci_high, self.df, self.bias
            )
        ]
        amax = max(np.abs(self.dp_target[2]).max(), np.abs(self.dp_test[2]).max())
        tail = (
            f"\nMSN refits failed in {self.n_failed} of {self.n_boot} replicates; "
            f"near-degenerate fits (|alpha| > {self.alpha_flag:g}) in "
            f"{self.n_degenerate} of {self.n_boot}\n"
            f"Full-data fits: largest |alpha| target = "
            f"{np.abs(self.dp_target[2]).max():.3g}, test = {np.abs(self.dp_test[2]).max():.3g}"
            + ("  <- near-degenerate" if amax > self.alpha_flag else "")
        )
        return head + "\n".join(rows) + tail


def spi_compare_bootstrap(
    target: pd.DataFrame | np.ndarray,
    test: pd.DataFrame | np.ndarray,
    target_groups: pd.Series | np.ndarray | list,
    test_groups: pd.Series | np.ndarray | list,
    *,
    n_boot: int = 2000,
    conf: float = 0.95,
    df: Literal["G-1", "satterthwaite"] | float = "G-1",
    n_sample: int = 1000,
    seed: int | None = None,
    shared_draws: bool = False,
    ks2d_func: Callable | None = None,
    fit_func: Callable | None = None,
    alpha_flag: float = 50.0,
    progress: bool = False,
) -> SPICompareResult:
    """
    Paired participant bootstrap of the empirical SPI AND the MSN-smoothed SPI.

    For every replicate (and for the full data), participants are resampled once;
    then, on that same resample:

    * empirical SPI  = 100 (1 - D(target rows, test rows))
    * MSN SPI        = fit an MSN to the target rows and to the test rows, draw
                       ``n_sample`` points from each fitted MSN, and compute
                       100 (1 - D(target draw, test draw)) (the published workflow)
    * difference     = MSN SPI - empirical SPI

    Draws come from fixed standard-normal arrays transformed with each
    replicate's fitted parameters (common random numbers), so the intervals
    reflect participant sampling and fit uncertainty, not Monte Carlo noise.

    Notes
    -----
    * The MSN fits treat ratings as independent (as ``selm`` does); the cluster
      bootstrap is what accounts for the dependence in the intervals.
    * The two SPIs sit on different scales even for identical distributions:
      the empirical SPI's ceiling falls with small n (roughly n_target * n_test /
      (n_target + n_test)), while the MSN SPI's ceiling is set by ``n_sample``.
      The difference is therefore partly a noise-floor effect, not only misfit.
    * Each replicate needs two MSN fits (R via rpy2 with the default fitter).
      Time one fit before choosing ``n_boot``. A replicate where a fit fails
      keeps its empirical value and gets NaN for the MSN columns.
    * ``fit_func(array (n, 2)) -> (xi, omega, alpha)`` can be replaced. Any
      exception raised by the default fitter (including R errors) is converted
      to RuntimeError and counted as a failed replicate.
    * Near-degenerate fits: skew-normal maximum likelihood can drive the shape
      parameter alpha towards infinity on small samples (a half-normal-like
      fit). Replicates where either fit has |alpha| > ``alpha_flag`` are
      counted (``n_degenerate``) and a warning is issued if the full-data fit
      is one of them; the MSN SPI then reflects the model more than the data.
    * Draws use one fixed set of standard-normal numbers (common random
      numbers), so the reported MSN point estimate carries one Monte Carlo
      realisation of draw noise (about +/-1-2 SPI points at n_sample = 1000;
      change ``seed`` to see it). The CIs exclude this noise.
    * ``shared_draws=True`` uses the same random numbers for both sides (see
      :func:`msn_spi`): no 1,000-draw floor, but not comparable with
      soundscapy's spi_score.
    * Draws are untruncated (as MultiSkewNorm.sample); soundscapy's
      sample_mtsn truncates to [-1, 1] and would give different values.
    * The default ``df="G-1"`` is applied to all three rows; "satterthwaite"
      estimates df separately per row via leave-one-participant-out refits
      (2G extra fits).
    """
    if ks2d_func is None:
        ks2d_func = fast_ks2d
    if fit_func is None:
        fit_func = _fit_msn_soundscapy

    def d_func(a: np.ndarray, b: np.ndarray) -> float:
        return float(ks2d_func(a, b)[0])

    rng = np.random.default_rng(seed)
    target_cl, test_cl = _paired_clusters(target, test, target_groups, test_groups)
    g = len(target_cl)
    z_a = rng.standard_normal((n_sample, 3))
    z_t = z_a if shared_draws else rng.standard_normal((n_sample, 3))

    def stats_for(a: np.ndarray, t: np.ndarray) -> np.ndarray:
        """[empirical, MSN, MSN - empirical, max |alpha|]; NaN if a fit fails."""
        emp = 100 * (1 - d_func(a, t))
        try:
            fa, ft = fit_func(a), fit_func(t)
            sa = msn_sample_from_params(*fa, z_a)
            st = msn_sample_from_params(*ft, z_t)
            msn = 100 * (1 - d_func(sa, st))
        except _FIT_ERRORS:
            return np.array([emp, np.nan, np.nan, np.nan])
        amax = max(np.abs(fa[2]).max(), np.abs(ft[2]).max())
        return np.array([emp, msn, msn - emp, amax])

    # --- point estimates (errors propagate here so a broken R setup is obvious) ---
    target_all = np.vstack(target_cl)
    test_all = np.vstack(test_cl)
    dp_a = fit_func(target_all)
    dp_t = fit_func(test_all)
    emp = 100 * (1 - d_func(target_all, test_all))
    msn = 100 * (
        1
        - d_func(
            msn_sample_from_params(*dp_a, z_a), msn_sample_from_params(*dp_t, z_t)
        )
    )
    est = np.array([emp, msn, msn - emp])
    amax_full = max(np.abs(dp_a[2]).max(), np.abs(dp_t[2]).max())
    if amax_full > alpha_flag:
        warnings.warn(
            f"Full-data MSN fit is near-degenerate (largest |alpha| = {amax_full:.3g}); "
            "the MSN SPI reflects an extreme fitted shape. Inspect dp_target / dp_test.",
            UserWarning,
            stacklevel=2,
        )

    # --- df per row ---
    if isinstance(df, str) and df == "G-1":
        dfs, df_method = np.full(3, g - 1.0), "G-1"
    elif isinstance(df, str) and df == "satterthwaite":
        infl = np.empty((g, 3))
        for i in range(g):
            a = np.vstack([c for j, c in enumerate(target_cl) if j != i])
            t = np.vstack([c for j, c in enumerate(test_cl) if j != i])
            infl[i] = est - stats_for(a, t)[:3]
        dfs = np.array([_df_from_influence(infl[:, c], g) for c in range(3)])
        df_method = "Satterthwaite per row (capped at G-1)"
    elif isinstance(df, Real) and not isinstance(df, bool):
        dfs, df_method = np.full(3, float(df)), "user-supplied"
    else:
        msg = "df must be 'G-1', 'satterthwaite' or a number."
        raise ValueError(msg)

    # --- bootstrap: the same participants drawn for both sets, both statistics ---
    iterator: object = range(n_boot)
    if progress:
        try:
            from tqdm import tqdm  # noqa: PLC0415

            iterator = tqdm(iterator, desc="bootstrap")
        except ImportError:
            pass
    boot4 = np.empty((n_boot, 4))
    for b in iterator:
        pick = rng.integers(0, g, size=g)
        a = np.vstack([target_cl[i] for i in pick])
        t = np.vstack([test_cl[i] for i in pick])
        boot4[b] = stats_for(a, t)
    boot, alpha_max = boot4[:, :3], boot4[:, 3]

    n_failed = int(np.isnan(boot[:, 1]).sum())
    if n_failed > 0.05 * n_boot:
        warnings.warn(
            f"{n_failed} of {n_boot} MSN refits failed; the MSN and difference "
            "intervals use the remaining replicates and may be unreliable.",
            UserWarning,
            stacklevel=2,
        )

    levels, lo_ci, hi_ci = [], np.empty(3), np.empty(3)
    for c in range(3):
        lo, hi = expanded_percentile_levels(g, dfs[c], conf)
        levels.append((lo, hi))
        lo_ci[c], hi_ci[c] = np.nanpercentile(boot[:, c], [lo, hi])

    return SPICompareResult(
        names=("empirical", "MSN-smoothed", "MSN - empirical"),
        estimate=est,
        ci_low=lo_ci,
        ci_high=hi_ci,
        df=dfs,
        df_method=df_method,
        percentile_levels=levels,
        bias=np.nanmean(boot, axis=0) - est,
        boot=boot,
        n_failed=n_failed,
        n_boot=n_boot,
        n_participants=g,
        n_target=len(target_all),
        n_test=len(test_all),
        n_sample=n_sample,
        conf=conf,
        dp_target=dp_a,
        dp_test=dp_t,
        alpha_max=alpha_max,
        alpha_flag=alpha_flag,
        n_degenerate=int(np.nansum(alpha_max > alpha_flag)),
    )