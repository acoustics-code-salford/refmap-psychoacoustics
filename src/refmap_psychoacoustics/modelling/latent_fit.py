"""
latent_fit.py: fit a model of your choice, optionally followed by a sigmoid link,
to observed responses, in ONE call.

    fit = fit_latent(df, predictors, obs="Y", lower=0, upper=10, se="Y_SE",
                     latent_model=my_function, sigmoid="gompertz")
    print(fit)          # parameters, fit statistics
    fit.table           # observed, predicted, latent score, latent score from observations

The chain is   predictors -> latent model -> latent score -> sigmoid -> predicted response.
The latent score is never observed; it only exists inside the model. With
sigmoid=None there is no link: the latent model IS the predicted response, and
this is simply a weighted, multi-seed, global-then-local curve fit of any function.

What you choose, and where
--------------------------
response      df, predictors, obs, se, lower, upper are arguments of fit_latent().
              Keep one dict per response type and unpack it: fit_latent(**RESPONSES["x"], ...)
sigmoid       "gompertz", "richards", "logistic", or None for no link (see SIGMOIDS).
latent model  a function you pass in, or a name in LATENT_MODELS. Its parameters are
              read from its signature, so changing the model means editing that one
              function and nothing else:
                * arguments named in `predictors` are data (columns of df)
                * every other argument is a parameter to fit
                * a parameter with a default value is held fixed at that value,
                  unless you give it bounds
bounds        explicit bounds={"name": (low, high)} win. Otherwise the function's own
              default bounds are used, attached with the @with_bounds decorator:

                  @with_bounds({"_a": (0, 50), "_b": (-4, 4), "k": (0, 1)})
                  def my_model(x, y, size_a, size_b, k): ...

              A key matches a parameter name exactly, or (if it starts with "_") the
              end of a name. A parameter with no bounds from either place is an error.
              Hold a parameter with fixed={"name": value}.

The sigmoid's growth rate
-------------------------
The sigmoid sees u = growth * latent score. If the latent model has free coefficients
that can absorb an overall scale, the growth rate is redundant with them, and fits
run off to infinity (coefficients up, growth down). So growth is fixed at 1 by
default. Use free_growth=True only for a latent model with no free scale.

How the fit works
-----------------
Stage 1, differential_evolution: a global search that lands in the right region
without needing a starting guess. Stage 2, least_squares: polishes that solution,
and is where `loss` (squared error, or a robust loss) is applied. This is repeated
for n_seeds seeds and the lowest cost is kept, because a single seed can land in a
poor local optimum.

Adding your own sigmoid: SIGMOIDS["name"] = (forward, inverse, [(param, low, high)]),
where forward(u, *params) returns a fraction of the span (0 to 1) and inverse(f, *params)
undoes it.
"""
from __future__ import annotations

import inspect
import re
import textwrap
import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution, least_squares

# --------------------------------------------------------------------------
# Sigmoids. Each is written on u = growth * latent score and returns a fraction
# of the span (0 to 1), with its inverse and its own parameters (name, low, high).
#   gompertz  f = exp(-c exp(-u))          inflection at u = ln(c), f = 1/e
#   richards  f = (1 + q exp(-u))^(-1/v)   v = 1 is logistic, v -> 0 is Gompertz
#   logistic  f = 1 / (1 + q exp(-u))      symmetric, midpoint at u = ln(q)
# --------------------------------------------------------------------------
def _gompertz(u, c):
    return np.exp(-np.exp(np.clip(np.log(c) - u, -700, 700)))

def _gompertz_inverse(f, c):
    return -np.log(-np.log(f) / c)

def _richards(u, q, v):
    return (1 + q * np.exp(np.clip(-u, -700, 700))) ** (-1 / v)

def _richards_inverse(f, q, v):
    return -np.log((f ** (-v) - 1) / q)

def _logistic(u, q):
    return 1 / (1 + q * np.exp(np.clip(-u, -700, 700)))

def _logistic_inverse(f, q):
    return -np.log((1 / f - 1) / q)

SIGMOIDS = {
    "gompertz": (_gompertz, _gompertz_inverse, [("sig_c", 1e-3, 100.0)]),
    "richards": (_richards, _richards_inverse, [("sig_q", 1e-3, 100.0), ("sig_v", 0.05, 10.0)]),
    "logistic": (_logistic, _logistic_inverse, [("sig_q", 1e-3, 100.0)]),
}
GROWTH_BOUNDS = (1e-4, 5.0)      # used only with free_growth=True

# Models selectable by name. Families register themselves here (see pa_models.py).
LATENT_MODELS = {}


def with_bounds(rules):
    """Decorator: attach default bounds to a latent-model function (see module docstring)."""
    def attach(func):
        func.default_bounds = dict(rules)
        return func
    return attach


def with_equation(text):
    """Decorator: attach a readable formula for the latent score z to a model function.
    Parameter names in the text are replaced by fitted values in LatentFit.equation()."""
    def attach(func):
        func.equation = text
        return func
    return attach


SIGMOID_TEXT = {"gompertz": "exp(-sig_c * exp(-u))",
                "richards": "(1 + sig_q * exp(-u)) ** (-1 / sig_v)",
                "logistic": "1 / (1 + sig_q * exp(-u))"}


def _default_bounds(func, name):
    for key, bnd in getattr(func, "default_bounds", {}).items():
        if name == key or (key.startswith("_") and name.endswith(key)):
            return bnd
    raise ValueError(f"'{name}' is treated as a parameter to fit, but has no bounds. If it "
                     f"is data, add it to `predictors`. Otherwise pass "
                     f"bounds={{'{name}': (low, high)}}, hold it with fixed={{'{name}': value}}, "
                     f"or attach default bounds to the function with @with_bounds.")


# --------------------------------------------------------------------------
# Internals
# --------------------------------------------------------------------------
class _Problem:
    """Arrays, model and parameter bookkeeping shared by every optimiser call."""

    def __init__(self, X, y, lower, upper, func, names, link, link_names, free, fixed):
        self.X, self.y, self.lower, self.upper = X, y, lower, upper
        self.func, self.names = func, names
        self.link, self.link_names = link, link_names     # link is None or (forward, inverse)
        self.free, self.fixed = free, fixed

    def values(self, theta):
        v = dict(self.fixed)
        v.update(zip(self.free, theta))
        return v

    def latent(self, v, X=None):
        return self.func(**(self.X if X is None else X), **{n: v[n] for n in self.names})

    def response(self, v, X=None):
        z = self.latent(v, X)
        if self.link is None:
            return z
        f = self.link[0](v["sig_b"] * z, *[v[n] for n in self.link_names])
        return self.lower + (self.upper - self.lower) * f

    def residuals(self, theta, weights):
        with np.errstate(all="ignore"):
            r = self.response(self.values(theta)) - self.y
        if weights is not None:
            r = r * np.sqrt(weights)
        # extreme trial parameters can give inf/NaN: a huge residual steers the
        # optimiser away instead of crashing it
        return np.where(np.isfinite(r), r, 1e6)


def _optimise(problem, bounds, weights, seeds, loss, f_scale, de_kwargs):
    """DE then least_squares for each seed; returns the lowest-cost result."""
    lb, ub = np.array(bounds, float).T
    settings = dict(maxiter=100, popsize=15, tol=1e-8, polish=False)
    settings.update(de_kwargs or {})
    runs = []
    for s in seeds:
        stage1 = differential_evolution(
            lambda t: np.sum(problem.residuals(t, weights) ** 2),
            list(zip(lb, ub)), seed=s, **settings)
        runs.append(least_squares(problem.residuals, stage1.x, bounds=(lb, ub),
                                  args=(weights,), loss=loss, f_scale=f_scale,
                                  x_scale="jac"))
    costs = [float(r.cost) for r in runs]
    return runs[int(np.argmin(costs))], costs


def _invert(problem, v, y, zero_as_eps=True):
    """Response scale -> latent scale. With a sigmoid this is defined only strictly
    between the asymptotes: values equal to `lower` become a tiny epsilon (if
    zero_as_eps), anything else on or beyond an asymptote becomes NaN. With no
    sigmoid the latent score is the response, so y is returned unchanged."""
    y = np.asarray(y, dtype=float)
    if problem.link is None:
        return y
    f = (y - problem.lower) / (problem.upper - problem.lower)
    if zero_as_eps:
        f = np.where(f == 0, np.finfo(float).eps, f)
    bad = ~((f > 0) & (f < 1)) & ~np.isnan(y)
    if bad.any():
        warnings.warn(f"{bad.sum()} value(s) on or beyond the asymptotes: returning "
                      "NaN for them (the latent score is undefined there).")
    f = np.where((f > 0) & (f < 1), f, np.nan)
    with np.errstate(all="ignore"):
        u = problem.link[1](f, *[v[n] for n in problem.link_names])
    return u / v["sig_b"]


@dataclass
class LatentFit:
    """Result of fit_latent()."""
    params: pd.Series              # every parameter, fitted and fixed
    fitted: list                   # names that were optimised
    table: pd.DataFrame            # observed, predicted, residual, latent, latent_from_obs, ...
    stats: dict                    # rmse, r2, reduced_chi2, n_obs, ...
    tau: float | None              # model-error term added to the SEs (weights only)
    weights: np.ndarray | None
    seed_costs: list               # cost of each seed's fit
    result: object                 # scipy OptimizeResult of the chosen fit
    latent_model: str
    sigmoid: str | None
    at_bound: list = field(default_factory=list)
    _problem: object = field(default=None, repr=False)
    _predictors: dict = field(default=None, repr=False)

    def _arrays(self, df):
        return {k: df[c].to_numpy(float) for k, c in self._predictors.items()}

    def predict(self, df, scale="response"):
        """Predictions for new data with the same predictor columns.
        scale="response" gives the fitted response, scale="latent" the latent score."""
        v, X = self.params.to_dict(), self._arrays(df)
        return self._problem.response(v, X) if scale == "response" else self._problem.latent(v, X)

    def latent_from_response(self, y, zero_as_eps=True):
        """Map observed responses back to the latent scale with the fitted sigmoid."""
        return _invert(self._problem, self.params.to_dict(), y, zero_as_eps)

    def equation(self, values=True):
        """The fitted model as readable equations. values=False keeps parameter names.
        The formula for z is the model function's `equation` (see with_equation),
        or its source code if it has none."""
        p, pr = self.params.to_dict(), self._problem

        def sub(text):
            if not values:
                return text
            for n in sorted(p, key=len, reverse=True):
                v = f"({p[n]:.4g})" if p[n] < 0 else f"{p[n]:.4g}"
                text = re.sub(rf"\b{n}\b", v, text)
            return text

        if pr.link is None:
            lines = ["response = z"]
        else:
            lines = [f"response = {pr.lower:g} + ({pr.upper:g} - {pr.lower:g}) * f(u)",
                     sub("u = sig_b * z"),
                     sub("f(u) = " + SIGMOID_TEXT.get(self.sigmoid, f"<{self.sigmoid}>"))]
        text = getattr(pr.func, "equation", None)
        if text:
            lines.append(sub("z = " + text))
        else:
            try:
                src = textwrap.indent(inspect.getsource(pr.func), "    ")
            except (OSError, TypeError):
                src = f"    {pr.func.__name__}(...)  (source not available)\n"
            fitted = ", ".join(f"{n}={p[n]:.4g}" for n in pr.names)
            lines.append(f"z = {pr.func.__name__}(...), where\n{src}  with {fitted}")
        return "\n".join(lines)

    def __str__(self):
        s = self.stats
        head = (f"latent_model={self.latent_model}  sigmoid={self.sigmoid}  n={s['n_obs']}"
                f"  free parameters={s['n_free']}\n"
                f"RMSE {s['rmse']:.3f}  R2 {s['r2']:.4f}")
        if "reduced_chi2" in s:
            head += f"  reduced chi2 {s['reduced_chi2']:.2f}  tau {self.tau:.4f}"
        head += f"\nseed costs: {[round(float(c), 4) for c in self.seed_costs]}\n"
        status = ["at bound" if n in self.at_bound else "fitted" if n in self.fitted
                  else "fixed" for n in self.params.index]
        table = pd.DataFrame({"value": self.params.round(6), "status": status})
        return head + table.to_string()


# --------------------------------------------------------------------------
# The one function you call
# --------------------------------------------------------------------------
def fit_latent(df, predictors, obs, lower=None, upper=None, se=None, *, latent_model,
               sigmoid="gompertz", weighting=None, floor_frac=0.1, tau=None,
               free_growth=False, bounds=None, fixed=None, loss="linear",
               f_scale=1.0, seed=0, n_seeds=3, de_kwargs=None):
    """Fit a latent model, plus an optional sigmoid link, to one response. Returns a LatentFit.

    df, predictors  the data, and {name the model function uses: column in df}
    obs, se         columns of the observed values and their standard errors
    lower, upper    the sigmoid's fixed asymptotes (the limits of the response);
                    not needed when sigmoid=None
    latent_model    a function, or a name in LATENT_MODELS (see the module docstring)
    sigmoid         "gompertz", "richards", "logistic", or None for no link
    weighting       None (default: "se_tau" if se is given, else "none"),
                    "se_tau", "se" or "none". Weights are 1 / (SE_floored^2 + tau^2).
                    The SE floor is floor_frac * median positive SE, so a zero SE
                    cannot give an infinite weight. tau (model error beyond rater
                    noise) is estimated from an unweighted pilot fit unless you
                    pass a number. Treat that estimate as a heuristic.
    free_growth     fit the sigmoid's growth rate too (see the module docstring)
    bounds, fixed   {name: (low, high)} and {name: value}, by parameter name
    loss, f_scale   passed to scipy's least_squares (e.g. "soft_l1"); with weights
                    f_scale is in standard errors, otherwise in response units
    seed, n_seeds   n_seeds fits from seeds seed, seed+1, ...; the best is kept
    de_kwargs       overrides for differential_evolution (maxiter, popsize, ...)
    """
    # ---- resolve the pieces -------------------------------------------------
    if isinstance(latent_model, str):
        if latent_model not in LATENT_MODELS:
            raise ValueError(f"latent_model '{latent_model}' is not registered. Registered: "
                             f"{list(LATENT_MODELS)}. Pass a function, or import the module "
                             "that registers the family (e.g. pa_models).")
        func, label = LATENT_MODELS[latent_model], latent_model
    else:
        func, label = latent_model, getattr(latent_model, "__name__", "custom")
    if sigmoid is None:
        link, link_spec = None, []
        if free_growth:
            raise ValueError("free_growth needs a sigmoid.")
    else:
        if sigmoid not in SIGMOIDS:
            raise ValueError(f"sigmoid must be one of {list(SIGMOIDS)} or None.")
        link, link_spec = SIGMOIDS[sigmoid][:2], SIGMOIDS[sigmoid][2]
        if lower is None or upper is None:
            raise ValueError("A sigmoid needs lower= and upper= (the response's limits).")
        if upper <= lower:
            raise ValueError("upper must be greater than lower.")

    # ---- parameters come from the model function's signature -----------------
    needed, names, defaults = [], [], {}
    for name, p in inspect.signature(func).parameters.items():
        if name in predictors:
            needed.append(name)
        else:
            names.append(name)
            if p.default is not inspect.Parameter.empty:
                defaults[name] = p.default
    if not needed:
        raise ValueError(f"None of the model function's arguments are in predictors "
                         f"{list(predictors)}.")

    user_bounds, fixed = dict(bounds or {}), dict(fixed or {})
    link_names = [n for n, _, _ in link_spec]
    valid = names + link_names + (["sig_b"] if link else [])
    for key in list(user_bounds) + list(fixed):
        if key not in valid:
            raise ValueError(f"Unknown parameter '{key}'. Valid names: {valid}")
    both = set(user_bounds) & set(fixed)
    if both:
        raise ValueError(f"{sorted(both)} given in both bounds and fixed.")
    if "sig_b" in user_bounds and not free_growth:
        raise ValueError("sig_b is fixed at 1 unless free_growth=True.")

    free, bnd = [], {}
    for name in names:
        if name in fixed:
            continue
        if name in user_bounds:
            bnd[name] = user_bounds[name]
        elif name in defaults:
            fixed[name] = defaults[name]
            continue
        else:
            bnd[name] = _default_bounds(func, name)
        free.append(name)
    for name, lo, hi in link_spec:
        if name not in fixed:
            bnd[name] = user_bounds.get(name, (lo, hi))
            free.append(name)
    if link:
        if free_growth and "sig_b" not in fixed:
            bnd["sig_b"] = user_bounds.get("sig_b", GROWTH_BOUNDS)
            free.append("sig_b")
        else:
            fixed.setdefault("sig_b", 1.0)

    # ---- data: only the columns the model uses; drop incomplete rows ---------
    if weighting is None:
        weighting = "se_tau" if se is not None else "none"
    if weighting not in ("se_tau", "se", "none"):
        raise ValueError("weighting must be None, 'se_tau', 'se' or 'none'.")
    if weighting != "none" and se is None:
        raise ValueError(f"weighting='{weighting}' needs se=<column of standard errors>.")
    for col in [predictors[k] for k in needed] + [obs] + ([se] if weighting != "none" else []):
        if col not in df.columns:
            raise ValueError(f"Column '{col}' is not in df.")
    X = {k: df[predictors[k]].to_numpy(float) for k in needed}
    y = df[obs].to_numpy(float)
    s = df[se].to_numpy(float) if weighting != "none" else None
    keep = np.isfinite(y) & np.all([np.isfinite(a) for a in X.values()], axis=0)
    if s is not None:
        keep &= np.isfinite(s)
    if (~keep).any():
        warnings.warn(f"Dropped {(~keep).sum()} row(s) with missing values.")
    if keep.sum() <= len(free):
        raise ValueError(f"Only {keep.sum()} usable rows for {len(free)} free parameters.")
    X, y = {k: a[keep] for k, a in X.items()}, y[keep]
    s = s[keep] if s is not None else None
    at_limit = (y <= lower) | (y >= upper) if link else np.zeros(len(y), bool)
    if at_limit.any():
        warnings.warn(f"{at_limit.sum()} observation(s) on or beyond the asymptotes: they are "
                      "fitted, but the latent score cannot be recovered from them (see "
                      "table['obs_at_limit']).")

    problem = _Problem(X, y, lower, upper, func, names, link, link_names, free, fixed)
    bounds_list = [bnd[n] for n in free]
    try:   # fail early, with the real error, if the model does not run
        problem.response(problem.values([np.mean(b) for b in bounds_list]))
    except Exception as e:
        raise RuntimeError(f"The latent model failed to run: {e}") from e

    # ---- weights --------------------------------------------------------------
    weights, tau_used = None, None
    if weighting != "none":
        positive = s[s > 0]
        if positive.size == 0:
            raise ValueError("All standard errors are zero; cannot build weights.")
        se_eff = np.maximum(s, floor_frac * np.median(positive))
        if weighting == "se":
            tau_used = 0.0
        elif tau is not None:
            tau_used = float(tau)
        else:
            pilot, _ = _optimise(problem, bounds_list, None, [seed], loss, f_scale, de_kwargs)
            pred0 = problem.response(problem.values(pilot.x))
            msr = np.sum((pred0 - y) ** 2) / max(len(y) - len(free), 1)
            tau_used = float(np.sqrt(max(0.0, msr - np.mean(se_eff ** 2))))
        weights = 1.0 / (se_eff ** 2 + tau_used ** 2)

    # ---- the fit ---------------------------------------------------------------
    best, costs = _optimise(problem, bounds_list, weights,
                            [seed + i for i in range(n_seeds)], loss, f_scale, de_kwargs)
    if max(costs) - min(costs) > 0.05 * max(min(costs), 1e-12):
        warnings.warn(f"Seeds gave costs differing by more than 5% {[round(c, 4) for c in costs]}; "
                      "the lowest was kept. Large differences mean local optima: check "
                      "bounds, add seeds, or simplify the model.")
    lb, ub = np.array(bounds_list, float).T
    at_bound = [n for n, hit in zip(free, np.isclose(best.x, lb) | np.isclose(best.x, ub)) if hit]
    if at_bound:
        warnings.warn(f"Parameters at a bound (consider widening it): {', '.join(at_bound)}")

    # ---- outputs ---------------------------------------------------------------
    v = problem.values(best.x)
    pred = problem.response(v)
    latent = problem.latent(v)
    if not np.all(np.isfinite(pred)):
        warnings.warn("Some fitted predictions are not finite (check the model's inputs).")
    table = pd.DataFrame({"observed": y, "predicted": pred, "residual": pred - y,
                          "latent": latent}, index=df.index[keep])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table["latent_from_obs"] = _invert(problem, v, y)
    table["obs_at_limit"] = at_limit
    if weights is not None:
        table["weight"] = weights
        table["std_residual"] = table["residual"] * np.sqrt(weights)
    sse = np.sum((pred - y) ** 2)
    stats = dict(n_obs=len(y), n_free=len(free), n_dropped=int((~keep).sum()),
                 rmse=float(np.sqrt(sse / len(y))),
                 r2=float(1 - sse / np.sum((y - y.mean()) ** 2)))
    if weights is not None:
        stats["reduced_chi2"] = float(np.sum(best.fun ** 2) / max(len(y) - len(free), 1))
    order = names + link_names + (["sig_b"] if link else [])
    return LatentFit(params=pd.Series({n: v[n] for n in order}), fitted=free, table=table,
                     stats=stats, tau=tau_used, weights=weights, seed_costs=costs,
                     result=best, latent_model=label, sigmoid=sigmoid, at_bound=at_bound,
                     _problem=problem, _predictors={k: predictors[k] for k in needed})