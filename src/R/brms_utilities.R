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