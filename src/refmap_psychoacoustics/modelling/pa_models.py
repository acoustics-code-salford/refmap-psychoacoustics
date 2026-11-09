"""
pa_models.py: the psychoacoustic annoyance (PA) family of latent models for latent_fit.

Importing this module registers the models by name, so either of these works:

    import pa_models
    fit_latent(..., latent_model="pa_4term")
    fit_latent(..., latent_model=pa_models.pa_4term)

Each model is a plain function. Arguments named in `predictors` are data and every
other argument is a parameter to fit. The @with_bounds rules give the default bounds
by name (a key starting with "_" matches the end of a name).

"""
import numpy as np

from .latent_fit import LATENT_MODELS, with_bounds, with_equation

PA_BOUNDS = {"_a": (0.0, 50.0), "_b": (-4.0, 4.0)}


@with_equation("""loud * (1 + sqrt(sharp_w**2 + tonal_w**2 + roughness_w**2 + fluctuation_w**2 + impulsive_w**2))
    where sharp_w = sharp_a * sharp / loud**sharp_b
          tonal_w = tonal_a * tonal / loud**tonal_b
          roughness_w = roughness_a * roughness / roughness**tonal_b
          fluctuation_w = fluctuation_a * fluctuation / fluctuation**tonal_b
          impulsive_w = impulsive_a * impulsive / loud**impulsive_b""")
@with_bounds(PA_BOUNDS)
def pa_5term(loud, sharp, tonal, rough, fluctuation, impulsive,
               sharp_a, sharp_b, tonal_a, tonal_b, rough_a, rough_b,
               fluctuation_a, fluctuation_b, impulsive_a, impulsive_b):
    """loud * (1 + sqrt(sum of squared weighted terms))."""
    sharp_w = sharp_a * sharp / loud**sharp_b
    tonal_w = tonal_a * tonal / loud**tonal_b
    rough_w = rough_a * rough / loud**rough_b
    fluctuation_w = fluctuation_a * fluctuation / loud**fluctuation_b
    impulsive_w = impulsive_a * impulsive / loud**impulsive_b
    return loud * (1 + np.sqrt(
        sharp_w**2 + tonal_w**2 + rough_w**2 + fluctuation_w**2 + impulsive_w**2))


@with_equation("""loud * (1 + sqrt(sharp_w**2 + tonal_w**2 + modulation_w**2 + impulsive_w**2))
    where sharp_w = sharp_a * sharp / loud**sharp_b
          tonal_w = tonal_a * tonal / loud**tonal_b
          modulation_w = modulation_a * (0.4*fluctuation + 0.6*rough) / loud**modulation_b
          impulsive_w = impulsive_a * impulsive / loud**impulsive_b""")
@with_bounds(PA_BOUNDS)
def pa_4term(loud, sharp, tonal, rough, fluctuation, impulsive,
               sharp_a, sharp_b, tonal_a, tonal_b, modulation_a, modulation_b,
               impulsive_a, impulsive_b):
    """As pa_5term, with roughness and fluctuation merged into one
    modulation term (0.4 * fluctuation + 0.6 * roughness)."""
    sharp_w = sharp_a * sharp / loud**sharp_b
    tonal_w = tonal_a * tonal / loud**tonal_b
    modulation_w = modulation_a * (0.4 * fluctuation + 0.6 * rough) / loud**modulation_b
    impulsive_w = impulsive_a * impulsive / loud**impulsive_b
    return loud * (1 + np.sqrt(
        sharp_w**2 + tonal_w**2 + modulation_w**2 + impulsive_w**2))


@with_equation("""loud * (1 + sqrt(sharp_w**2 + tonal_w**2 + modulation_w**2 + impulsive_w**2))
    where sharp_w = sharp_a * max(sharp - sharp_c, 0) *sharp_b * log10(loud + 10)
          tonal_w = tonal_a * tonal / loud**tonal_b
          modulation_w = modulation_a * (0.4*fluctuation + 0.6*rough) / loud**modulation_b
          impulsive_w = impulsive_a * impulsive / loud**impulsive_b""")
@with_bounds(PA_BOUNDS)
def pa_original(loud, sharp, tonal, rough, fluctuation, impulsive,
                   sharp_a, sharp_b, sharp_c, tonal_a, tonal_b, modulation_a, modulation_b,
                   impulsive_a, impulsive_b):
    """loud * (1 + sqrt(...)). Has no free overall scale, so use free_growth=True."""

    sharp_w = sharp_a * np.maximum(sharp - sharp_c, 0) * sharp_b * np.log10(loud + 10)
    tonal_w = tonal_a * tonal / loud**tonal_b
    modulation_w = modulation_a * (0.4 * fluctuation + 0.6 * rough) / loud**modulation_b
    impulsive_w = impulsive_a * impulsive / loud**impulsive_b
    return loud * (1 + np.sqrt(
        sharp_w**2 + tonal_w**2 + modulation_w**2 + impulsive_w**2))


@with_equation("""loud * (1 + sqrt(sharp_w**2 + tonal_w**2 + modulation_w**2 + impulsive_w**2))
    where sharp_w = sharp_a * sharp / loud**sharp_b
          tonal_w = tonal_a * tonal / loud**tonal_b
          modulation_w = modulation_a * (0.4*fluctuation + 0.6*rough) / loud**modulation_b
          impulsive_w = impulsive_a * impulsive / loud**impulsive_b""")
@with_bounds(PA_BOUNDS)
def pa_4term_log(loud, sharp, tonal, rough, fluctuation, impulsive,
                 sharp_a, sharp_b, tonal_a, tonal_b, modulation_a, modulation_b,
                 impulsive_a, impulsive_b):
    """log of loud * (1 + sqrt(...)). Use with free_growth=True: the growth rate
    then acts as a power-law exponent on PA."""
    pa = pa_4term(loud, sharp, tonal, rough, fluctuation, impulsive,
                  sharp_a, sharp_b, tonal_a, tonal_b, modulation_a, modulation_b,
                  impulsive_a, impulsive_b)
    return np.log(pa)

 
LATENT_MODELS.update({"pa_4term": pa_4term, "pa_5term": pa_5term,
                      "pa_original": pa_original, "pa_4term_log": pa_4term_log})