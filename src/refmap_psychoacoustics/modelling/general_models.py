"""
general_models.py: general models for latent_fit.

Importing this module registers the models by name, so either of these works:

    import general_models
    fit_latent(..., latent_model="model_name")
    fit_latent(..., latent_model=general_models.model_name)

Each model is a plain function. Arguments named in `predictors` are data and every
other argument is a parameter to fit. The @with_bounds rules give the default bounds
by name (a key starting with "_" matches the end of a name).

"""
import numpy as np

from .latent_fit import LATENT_MODELS, with_bounds, with_equation

BOUNDS = {"_a": (-50.0, 50.0)}


@with_equation("""beta_a * predictor""")
@with_bounds(BOUNDS)
def lin_1term(predictor,
               beta_a):
    """Linear model with one term."""
    return beta_a * predictor


@with_equation("""beta1_a * predictor1 + beta2_a * predictor2""")
@with_bounds(BOUNDS)
def lin_2term(predictor1, predictor2,
               beta1_a, beta2_a):
    """Linear model with one term."""
    return beta1_a * predictor1 + beta2_a * predictor2
 
LATENT_MODELS.update({"lin_1term": lin_1term, "lin_2term": lin_2term})