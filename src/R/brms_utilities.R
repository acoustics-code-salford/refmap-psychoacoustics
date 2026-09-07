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
      get_metric(mod, "bayes_R2_conditional") |> tibble::add_column(Type = "Conditional", .before = 1)
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
    dplyr::relocate(starts_with("Marginal"), .after = Model) |> 
    
    # Reorder rows by the marginal R2 Estimate descending
    dplyr::arrange(desc(Marginal_Estimate))
}


# brms formula builder -------------------------------------------------
# Helper function to build the formula programmatically
build_bf <- function(y, 
                     fixed = NULL, 
                     intercept = TRUE,
                     re = NULL, 
                     sigma_fixed = NULL, 
                     sigma_intercept = TRUE,
                     sigma_re = NULL) {
  
  # Internal engine to assemble formula strings
  construct_str <- function(lhs, fixed_vec, inc_intercept, re_list) {
    # 1. Fixed effects handling
    if (length(fixed_vec) > 0) {
      fixed_body <- paste(fixed_vec, collapse = " + ")
      fixed_str <- if (inc_intercept) fixed_body else paste("0 +", fixed_body)
    } else {
      fixed_str <- if (inc_intercept) "1" else "0"
    }
    
    # 2. Random effects handling
    re_str <- ""
    if (is.list(re_list) && length(re_list) > 0) {
      re_parts <- character(length(re_list))
      
      for (i in seq_along(re_list)) {
        grp_name <- names(re_list)[i]
        item <- re_list[[i]]
        
        # Parse list or atomic vector options
        if (is.list(item)) {
          slopes <- item$slopes
          cor <- if (!is.null(item$cor)) item$cor else TRUE
          re_inc_intercept <- if (!is.null(item$intercept)) item$intercept else TRUE
        } else {
          slopes <- item
          cor <- TRUE
          re_inc_intercept <- TRUE
        }
        
        # Select pipe operator
        pipe <- if (cor) "|" else "||"
        
        # Build random term LHS
        int_prefix <- if (re_inc_intercept) "1" else "0"
        
        if (length(slopes) > 0) {
          re_parts[i] <- sprintf("(%s + %s %s %s)", int_prefix, paste(slopes, collapse = " + "), pipe, grp_name)
        } else {
          # If no slopes and no intercept, default to intercept (0 | group is invalid syntax)
          if (!re_inc_intercept) {
            warning(sprintf("Group '%s' specified no intercept and no slopes; defaulting to (1 | %s)", grp_name, grp_name))
            int_prefix <- "1"
          }
          re_parts[i] <- sprintf("(%s %s %s)", int_prefix, pipe, grp_name)
        }
      }
      re_str <- paste0(" + ", paste(re_parts, collapse = " + "))
    }
    
    if (is.null(lhs)) {
      return(sprintf("%s ~ %s%s", y, fixed_str, re_str))
    } else {
      return(sprintf("%s ~ %s%s", lhs, fixed_str, re_str))
    }
  }
  
  # Build location formula
  f_mu <- as.formula(construct_str(lhs = NULL, fixed, intercept, re))
  
  if (is.null(sigma_fixed) && is.null(sigma_re)) {
    return(brms::bf(f_mu))
  }
  
  # Build scale formula
  f_sigma <- as.formula(construct_str(lhs = "sigma", sigma_fixed, sigma_intercept, sigma_re))
  
  return(brms::bf(f_mu, f_sigma))
}
