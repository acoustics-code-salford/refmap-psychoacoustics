# brms_utilities.R
require(brms)
require(ggridges)
require(tidyverse)


# Prior predictive distribution ppd_ridgeplot ----------------------------------
## Adapted from: https://bruno.nicenboim.me/posts/posts/2026-01-09-ordinal-models/index.html

ppd_ridgeplot <- function(fit, title = "Prior Predictive Distribution", 
                          subtitle = NULL, ndraws = 500) {
  
  yrep <- brms::posterior_predict(fit, ndraws = ndraws)
  
  # extract observed response categories from model
  response_name <- all.vars(fit$formula$formula)[1]
  
  observed_y <- fit$data[[response_name]]
  
  # preserve ordinal ordering
  response_levels <- sort(unique(observed_y))
  
  proportions_per_draw <- lapply(seq_len(nrow(yrep)), function(i) {
    props <- table(factor(yrep[i, ], levels = response_levels)) / ncol(yrep)
    data.frame(
      draw = i,
      response = factor(response_levels,
                        levels = response_levels,
                        ordered = TRUE),
      proportion = as.numeric(props)
    )
  })
  
  ppd_data <- do.call(rbind, proportions_per_draw)
  
  # reference line for equal occupancy
  equal_prob <- 1 / length(response_levels)
  
  ggplot(ppd_data, aes(x = proportion, y = response)) +
    ggridges::geom_density_ridges(fill = "steelblue", alpha = 0.7, scale = 0.9, stat = "binline") +
    scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
    geom_vline(xintercept = equal_prob, linetype = "dashed") +
    labs(title = title,
         subtitle = subtitle,
         x = "Predicted category proportion",
         y = "Response category") +
    theme_minimal(base_size = 12) +
    coord_flip()
  
}

# Bayesian R2 table extraction function ----------------------------------------

get_Br2_table <- function(...) {
  model_names <- match.call(expand.dots = FALSE)$... |> as.character()
  models <- list(...)
  
  # Helper to safely extract metrics or return NAs matching the column shape
  get_metric <- function(model, type) {
    metric <- model$criteria[[type]]
    
    if (!is.null(metric)) {
      return(tibble::as_tibble(metric))
    }
    
    # Dynamically find column names from a model that actually has data
    valid_mod <- purrr::keep(models, ~ !is.null(.x$criteria[[type]])) |> first()
    col_names <- if (!is.null(valid_mod)) colnames(valid_mod$criteria[[type]]) else c("Estimate", "Est.Error", "Q2.5", "Q97.5")
    
    # Return a single row of NAs with correct columns
    tibble::as_tibble(matrix(NA_real_, nrow = 1, ncol = length(col_names), dimnames = list(NULL, col_names)))
  }
  
  # Loop through models, extract, flatten to one row, and sort
  map2_dfr(models, model_names, function(mod, name) {
    dplyr::bind_rows(
      get_metric(mod, "bayes_R2_marginal")    |> tibble::add_column(Type = "Marginal", .before = 1),
      get_metric(mod, "bayes_R2_conditional") |> tibble::add_column(Type = "Conditional", .before = 1),
      get_metric(mod, "bayes_R2_MZ_marginal") |> tibble::add_column(Type = "McKelvey&Zavoina Marginal", .before = 1),
      get_metric(mod, "bayes_R2_MZ_conditional") |> tibble::add_column(Type = "McKelvey&Zavoina Conditional", .before = 1)
    ) |> 
      tibble::add_column(Model = name, .before = 1)
  }) |> 
    # Pivot columns wider so each model has exactly one row
    tidyr::pivot_wider(
      names_from = Type, 
      values_from = -c(Model, Type),
      names_glue = "{Type}_{.value}"
    ) |> 
    
    # Reorder columns so Marginal statistics block together, then Conditional block together
    dplyr::relocate(starts_with("McKelvey&Zavoina Conditional"), .after = Model) |>
    dplyr::relocate(starts_with("McKelvey&Zavoina Marginal"), .after = Model) |>
    dplyr::relocate(starts_with("Conditional"), .after = Model) |>
    dplyr::relocate(starts_with("Marginal"), .after = Model) |> 
    
    # Reorder rows by the marginal R2 Estimate descending
    dplyr::arrange(desc(Marginal_Estimate))
}



# build_bf(): generalized brms formula builder ----------------------------
#
# FIX vs. the original: the auxiliary-parameter block was hardcoded to
# "sigma" as the LHS, so it could only ever build sigma ~ ... formulas -
# unusable for phi (beta/ordbeta), disc (cumulative), zi, etc. This version
# takes a named list `dpars`, one entry per auxiliary parameter you want a
# formula for, keyed by whatever name that family actually uses (confirm via
# names(insight::find_formula(your_model)) if unsure). Also supports building
# MORE THAN ONE auxiliary formula at once (e.g. disc AND zi together), which
# the single-block design could not do regardless of naming.
#
# OTHER FIXES (see accompanying chat message for full explanation of each):
#   - re_list MUST now be a named list; an unnamed list previously produced
#     a silently broken formula referencing a literal "NA" grouping factor.
#   - added `tag` support per group, for the "(term | tag | group)" syntax
#     used to correlate group-level effects ACROSS formulas (e.g. correlating
#     mu's and disc's by-ID intercepts) - the original had no way to express
#     this at all, only the |/|| toggle for within-formula correlation.
#   - fixed the guard that decides whether to build an auxiliary formula at
#     all: previously only checked fixed/re, silently ignoring an
#     intercept-only request (e.g. wanting `disc ~ 0` with no other terms).
#   - added validation for `y` (response variable name).
#
# NOT executed here - no R access in this environment. Please run and report
# back if anything errors or looks off.

build_bf <- function(y,
                     fixed = NULL,
                     intercept = TRUE,
                     re = NULL,
                     dpars = list(),   # named list: dpars$disc, dpars$phi, dpars$sigma, ...
                     ...) {
  
  if (is.null(y) || !nzchar(y)) {
    stop("`y` (the response variable name) must be a non-empty string.")
  }
  

  # Internal engine to assemble one formula string (used for mu AND for every
  # entry in `dpars`) - unchanged in logic from the original, plus the fixes
  # listed above.

  construct_str <- function(lhs, fixed_vec, inc_intercept, re_list) {
    
    # 1. Fixed effects handling (unchanged from original)
    if (length(fixed_vec) > 0) {
      fixed_body <- paste(fixed_vec, collapse = " + ")
      fixed_str <- if (inc_intercept) fixed_body else paste("0 +", fixed_body)
    } else {
      fixed_str <- if (inc_intercept) "1" else "0"
    }
    
    # 2. Random effects handling
    re_str <- ""
    if (is.list(re_list) && length(re_list) > 0) {
      
      if (is.null(names(re_list)) || any(!nzchar(names(re_list)))) {
        stop(
          "`re` (and each `dpars$<name>$re`) must be a NAMED list, one entry ",
          "per grouping factor, e.g. list(ID = list(slopes = c('x1','x2'))). ",
          "An unnamed list silently produced a broken formula (grouping by a ",
          "literal 'NA') in the previous version of this function - this is ",
          "now a hard error instead."
        )
      }
      
      re_parts <- character(length(re_list))
      
      for (i in seq_along(re_list)) {
        grp_name <- names(re_list)[i]
        item <- re_list[[i]]
        
        # Parse list or atomic vector options
        if (is.list(item)) {
          slopes <- item$slopes
          cor <- if (!is.null(item$cor)) item$cor else TRUE
          re_inc_intercept <- if (!is.null(item$intercept)) item$intercept else TRUE
          tag <- item$tag  # NEW: arbitrary string to correlate this group's
          # effects with the SAME tag used elsewhere (e.g.
          # in the mu formula's re, or another dpar's re) -
          # see brms's "(term | tag | group)" syntax.
        } else {
          slopes <- item
          cor <- TRUE
          re_inc_intercept <- TRUE
          tag <- NULL
        }
        
        # Select the group-identifier portion: plain "group", "|| group"
        # (handled via pipe below), or "| tag | group" if a tag is given.
        # A tag always implies correlated effects, so `cor` is ignored (with
        # a warning) if both are supplied inconsistently.
        if (!is.null(tag)) {
          if (isFALSE(cor)) {
            warning(sprintf(
              "Group '%s': `tag` was supplied, which implies correlated ",
              "effects; ignoring cor = FALSE.", grp_name
            ))
          }
          pipe_group <- sprintf("%s | %s", tag, grp_name)
        } else {
          pipe_group <- grp_name
        }
        pipe <- if (is.null(tag) && !cor) "||" else "|"
        
        int_prefix <- if (re_inc_intercept) "1" else "0"
        
        if (length(slopes) > 0) {
          re_parts[i] <- sprintf("(%s + %s %s %s)", int_prefix,
                                 paste(slopes, collapse = " + "),
                                 pipe, pipe_group)
        } else {
          if (!re_inc_intercept) {
            warning(sprintf(
              "Group '%s' specified no intercept and no slopes; defaulting to (1 | %s)",
              grp_name, grp_name
            ))
            int_prefix <- "1"
          }
          re_parts[i] <- sprintf("(%s %s %s)", int_prefix, pipe, pipe_group)
        }
      }
      re_str <- paste0(" + ", paste(re_parts, collapse = " + "))
    }
    
    sprintf("%s ~ %s%s", lhs, fixed_str, re_str)
  }
  

  # Build the mu (location) formula - always required

  f_mu <- as.formula(construct_str(y, fixed, intercept, re))
  

  # Build zero or more auxiliary-parameter formulas from `dpars`

  aux_formulas <- list()
  
  for (dpar_name in names(dpars)) {
    spec <- dpars[[dpar_name]]
    
    spec_fixed     <- spec$fixed
    spec_intercept <- if (!is.null(spec$intercept)) spec$intercept else TRUE
    spec_re        <- spec$re
    
    # FIX: also build the formula if intercept was explicitly turned off,
    # even with no fixed terms/re specified (e.g. a bare `disc ~ 0` request) -
    # the original only checked fixed/re and silently dropped this case.
    has_content <- length(spec_fixed) > 0 || length(spec_re) > 0 || isFALSE(spec_intercept)
    
    if (has_content) {
      aux_formulas[[dpar_name]] <- as.formula(
        construct_str(dpar_name, spec_fixed, spec_intercept, spec_re)
      )
    }
  }
  
  if (length(aux_formulas) == 0) {
    return(brms::bf(f_mu, ...))
  }
  
  do.call(brms::bf, c(list(f_mu), unname(aux_formulas), list(...)))
}


## Example usage ---------------------------------

if (FALSE) {
  
  ## --- ordered-beta phi model (this session's m1d) ---------------------------
  build_bf(
    y = "dAnnoyanceOrdBetaScl",
    fixed = c("TrialNumberScl", "UASProximity", "UASLAEMaxLRScl*UASEvents"),
    re = list(ID = list(slopes = c("UASProximity", "UASLAEMaxLRScl"))),
    dpars = list(
      phi = list(fixed = c("AmbientEnvCore", "UASOperation", "UASType"),
                 re = list(ID = list(slopes = NULL)))
    )
  )
  
  ## cumulative disc model (this session's identifiability requirement) ----
  # NB: intercept = FALSE is required here per Bürkner & Vuorre (2019) - one
  # level of a categorical disc predictor must anchor disc = 1, or the model
  # is unidentified. This function will happily build "disc ~ 0 + ..." but
  # will NOT catch a missing intercept=FALSE for you - that's a modelling
  # decision, not something a formula-string builder can validate.
  build_bf(
    y = "AnnoyanceOrd",
    fixed = c("UASLAEMaxLRScl*I(log10(UASEvents))", "AmbientEnv*UASOperation*UASLAEMaxLRScl"),
    re = list(ID = list(slopes = c("UASLAEMaxLRScl", "AmbientEnv"), tag = "q")),
    dpars = list(
      disc = list(fixed = c("AmbientEnv", "UASOperation", "UASType"),
                  intercept = FALSE,
                  re = list(ID = list(slopes = NULL, tag = "q")))  # same tag "q"
      # as mu's ID group above -> correlates mu's and disc's
      # by-ID intercepts, per brms's shared-tag syntax.
    )
  )
  
  # --- multiple auxiliary formulas at once (not possible in the original) ---
  build_bf(
    y = "y", fixed = "x1",
    dpars = list(
      phi = list(fixed = "x2"),
      zi  = list(fixed = "x3", intercept = FALSE)
    )
  )
}