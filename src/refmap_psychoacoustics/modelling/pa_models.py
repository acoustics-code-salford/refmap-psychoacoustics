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

PA_BOUNDS = {"_a": (0.0, 20.0), "_b": (-3.0, 3.0)}


@with_equation("""loudness * (1 + sqrt(sharpness_w**2 + tonality_w**2 + roughness_w**2 + fluctuation_w**2 + impulsive_w**2))
    where sharpness_w = sharp_a * sharpness / loudness**sharp_b
          tonality_w = tonal_a * tonality / loudness**tonal_b
          roughness_w = rough_a * roughness / loudness**rough_b
          fluctuation_w = fluct_a * fluctuation / loudness**fluct_b
          impulsive_w = impuls_a * impulsive / loudness**impuls_b""")
@with_bounds(PA_BOUNDS)
def pa_5term(loudness, sharpness, tonality, roughness, fluctuation, impulsive,
               sharp_a, sharp_b, tonal_a, tonal_b, rough_a, rough_b,
               fluct_a, fluct_b, impuls_a, impuls_b):
    """loud * (1 + sqrt(sum of squared weighted terms))."""
    sharpness_w = sharp_a * sharpness / loudness**sharp_b
    tonality_w = tonal_a * tonality / loudness**tonal_b
    roughness_w = rough_a * roughness / loudness**rough_b
    fluctuation_w = fluct_a * fluctuation / loudness**fluct_b
    impulsive_w = impuls_a * impulsive / loudness**impuls_b
    return loudness * (1 + np.sqrt(
        sharpness_w**2 + tonality_w**2 + roughness_w**2 + fluctuation_w**2 + impulsive_w**2))


@with_equation("""loudness * (1 + sqrt(sharpness_w**2 + tonality_w**2 + modulation_w**2 + impulsive_w**2))
    where sharpness_w = sharp_a * sharpness / loudness**sharp_b
          tonality_w = tonal_a * tonality / loudness**tonal_b
          modulation_w = modulat_a * (0.4*fluctuation + 0.6*roughness) / loudness**modulat_b
          impulsive_w = impuls_a * impulsive / loudness**impuls_b""")
@with_bounds(PA_BOUNDS)
def pa_4term(loudness, sharpness, tonality, roughness, fluctuation, impulsive,
               sharp_a, sharp_b, tonal_a, tonal_b, modulat_a, modulat_b,
               impuls_a, impuls_b):
    """As pa_5term, with roughness and fluctuation merged into one
    modulation term (0.4 * fluctuation + 0.6 * roughness)."""
    sharpness_w = sharp_a * sharpness / loudness**sharp_b
    tonality_w = tonal_a * tonality / loudness**tonal_b
    modulation_w = modulat_a * (0.4 * fluctuation + 0.6 * roughness) / loudness**modulat_b
    impulsive_w = impuls_a * impulsive / loudness**impuls_b
    return loudness * (1 + np.sqrt(
        sharpness_w**2 + tonality_w**2 + modulation_w**2 + impulsive_w**2))


@with_equation("""loudness * (1 + sqrt(sharpness_w**2 + tonality_w**2 + modulation_w**2 + impulsive_w**2))
    where sharpness_w = sharp_a * max(sharpness - sharp_c, 0) *sharp_b * log10(loudness + 10)
          tonality_w = tonal_a * tonality / loudness**tonal_b
          modulation_w = modulat_a * (0.4*fluctuation + 0.6*roughness) / loudness**modulat_b
          impulsive_w = impuls_a * impulsive / loudness**impuls_b""")
@with_bounds(PA_BOUNDS)
def pa_original(loudness, sharpness, tonality, roughness, fluctuation, impulsive,
                   sharp_a, sharp_b, sharp_c, tonal_a, tonal_b, modulat_a, modulat_b,
                   impuls_a, impuls_b):
    """loudness * (1 + sqrt(...)). Has no free overall scale, so use free_growth=True."""

    sharpness_w = sharp_a * np.maximum(sharpness - sharp_c, 0) * sharp_b * np.log10(loudness + 10)
    tonality_w = tonal_a * tonality / loudness**tonal_b
    modulation_w = modulat_a * (0.4 * fluctuation + 0.6 * roughness) / loudness**modulat_b
    impulsive_w = impuls_a * impulsive / loudness**impuls_b
    return loudness * (1 + np.sqrt(
        sharpness_w**2 + tonality_w**2 + modulation_w**2 + impulsive_w**2))


@with_equation("""loudness * (1 + sqrt(sharpness_w**2 + tonality_w**2 + modulation_w**2 + impulsive_w**2))
    where sharpness_w = sharp_a * sharpness / loudness**sharp_b
          tonality_w = tonal_a * tonality / loudness**tonal_b
          modulation_w = modulat_a * (0.4*fluctuation + 0.6*roughness) / loudness**modulat_b
          impulsive_w = impuls_a * impulsive / loudness**impuls_b""")
@with_bounds(PA_BOUNDS)
def pa_4term_log(loudness, sharpness, tonality, roughness, fluctuation, impulsive,
                 sharp_a, sharp_b, tonal_a, tonal_b, modulat_a, modulat_b,
                 impuls_a, impuls_b):
    """log of loudness * (1 + sqrt(...)). Use with free_growth=True: the growth rate
    then acts as a power-law exponent on PA."""
    pa = pa_4term(loudness, sharpness, tonality, roughness, fluctuation, impulsive,
                  sharp_a, sharp_b, tonal_a, tonal_b, modulat_a, modulat_b,
                  impuls_a, impuls_b)
    return np.log(pa)

 
LATENT_MODELS.update({"pa_4term": pa_4term, "pa_5term": pa_5term,
                      "pa_original": pa_original, "pa_4term_log": pa_4term_log})