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

BOUNDS = {"_a": (-75.0, 75.0)}


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


@with_equation("""beta1_a * predictor1 + beta2_a * predictor2 + beta3_a * predictor3""")
@with_bounds(BOUNDS)
def lin_3term(predictor1, predictor2, predictor3,
               beta1_a, beta2_a, beta3_a):
    """Linear model with three terms."""
    return beta1_a * predictor1 + beta2_a * predictor2 + beta3_a * predictor3


@with_equation("""beta1_a * predictor1 + beta2_a * predictor2 + beta3_a * predictor3 + beta4_a * predictor4""")
@with_bounds(BOUNDS)
def lin_4term(predictor1, predictor2, predictor3, predictor4,
               beta1_a, beta2_a, beta3_a, beta4_a):
    """Linear model with four terms."""
    return beta1_a * predictor1 + beta2_a * predictor2 + beta3_a * predictor3 + beta4_a * predictor4
 

@with_equation("""beta1_a * predictor1 + beta2_a * predictor2 + beta3_a * predictor3 + beta4_a * predictor4 + beta5_a * predictor5""")
@with_bounds(BOUNDS)
def lin_5term(predictor1, predictor2, predictor3, predictor4, predictor5,
               beta1_a, beta2_a, beta3_a, beta4_a, beta5_a):
    """Linear model with five terms."""
    return beta1_a * predictor1 + beta2_a * predictor2 + beta3_a * predictor3 + beta4_a * predictor4 + beta5_a * predictor5


@with_equation("""beta1_a * predictor1 + beta2_a * predictor2 + beta3_a * predictor3
+ beta4_a * predictor4 + beta5_a * predictor5 + beta6_a * predictor6""")
@with_bounds(BOUNDS)
def lin_6term(predictor1, predictor2, predictor3, predictor4, predictor5, predictor6,
               beta1_a, beta2_a, beta3_a, beta4_a, beta5_a, beta6_a):
    """Linear model with six terms."""
    return beta1_a * predictor1 + beta2_a * predictor2 + beta3_a * predictor3 + beta4_a * predictor4 + beta5_a * predictor5 + beta6_a * predictor6


@with_equation("""beta1_a * predictor1 + beta2_a * predictor2 + beta3_a * predictor3
+ beta4_a * predictor4 + beta5_a * predictor5 + beta6_a * predictor6 + beta7_a * predictor7""")
@with_bounds(BOUNDS)
def lin_7term(predictor1, predictor2, predictor3, predictor4, predictor5, predictor6, predictor7,
              beta1_a, beta2_a, beta3_a, beta4_a, beta5_a, beta6_a, beta7_a):
    """Linear model with seven terms."""
    return beta1_a * predictor1 + beta2_a * predictor2 + beta3_a * predictor3 + beta4_a * predictor4 + beta5_a * predictor5 + beta6_a * predictor6 + beta7_a * predictor7


LATENT_MODELS.update({"lin_1term": lin_1term, "lin_2term": lin_2term,
                      "lin_3term": lin_3term, "lin_4term": lin_4term,
                      "lin_5term": lin_5term, "lin_6term": lin_6term,
                      "lin_7term": lin_7term})