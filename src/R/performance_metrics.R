require(insight)
require(dplyr)
require(performance)


## GEE Rsquared (Zheng, 2000)
geeRsquared <- function(...){
  
  model_objects <- list(...)
  
  # if (length(list(...)) == 1){
  #   
  #   model_objects <- insight::ellipsis_info(..., ..., only_models = TRUE,
  #                                           verbose=FALSE)
  #   model_objects <- model_objects[1]
  #   
  # } else{
  #   
  #   model_objects <- insight::ellipsis_info(..., only_models = TRUE,
  #                                           verbose=FALSE)
  #   
  # }
  
  # ensure proper object names
  dot_names <- sapply(match.call(expand.dots = FALSE)[["..."]], as.character)
  check_object_names <- insight::compact_character(names(model_objects))
  if ((is.null(check_object_names) ||
       # or if length of names doesn't match number of models
       length(check_object_names) != length(model_objects) ||
       # or if names are "..1", "..2" pattern
       all(grepl("\\.\\.\\d", check_object_names))) &&
      # and length of dot-names must match length of objects
      length(model_objects) == length(dot_names)) {
    names(model_objects) <- dot_names
  }
  
  object_names <- names(model_objects)
  
  m <- mapply(function(.x, .y) {
    SSE <- sum(.x$prior.weights*(.x$y - .x$fitted.values)^2, na.rm = T)
    SST <- sum(.x$prior.weights*(.x$y - mean(.x$y, na.rm = T))^2, na.rm = T)
    r2 = 1 - SSE/SST
    
    dat <- data.frame(R2=r2)
    model_name <- gsub("\"", "", insight::safe_deparse(.y), fixed = TRUE)
    perf_df <- data.frame(Name = model_name, Model = class(.x)[1], dat, stringsAsFactors = FALSE)
    perf_df
  }, model_objects, object_names, SIMPLIFY = FALSE)
  
  dfs <- Reduce(function(x, y) merge(x, y, all = TRUE, sort = FALSE), m)
  
  return(dfs)
}


## GEE entropy reduction (Zheng, 2000)

geeEntropy <- function(..., varest = c("robust", "df-adjusted", "model", "bias-corrected")) {
  
  varest <- match.arg(varest)
  
  model_objects <- list(...)
  
  # if (length(list(...)) == 1) {
  #   model_objects <- insight::ellipsis_info(..., ..., only_models = TRUE, verbose = FALSE)
  #   model_objects <- model_objects[1]
  # } else {
  #   model_objects <- insight::ellipsis_info(..., only_models = TRUE, verbose = FALSE)
  # }
  
  dot_names <- sapply(match.call(expand.dots = FALSE)[["..."]], as.character)
  check_object_names <- insight::compact_character(names(model_objects))
  if ((is.null(check_object_names) ||
       length(check_object_names) != length(model_objects) ||
       all(grepl("\\.\\.\\d", check_object_names))) &&
      length(model_objects) == length(dot_names)) {
    names(model_objects) <- dot_names
  }
  
  object_names <- names(model_objects)
  
  m <- mapply(function(.x, .y) {
    
    y  <- .x$y
    mu <- .x$fitted.values
    
    alpha <- mean(y, na.rm = TRUE)
    
    eps <- 1e-15
    mu    <- pmin(pmax(mu, eps), 1 - eps)
    alpha <- min(max(alpha, eps), 1 - eps)
    
    numerator <- sum(mu * log(mu) + (1 - mu) * log(1 - mu), na.rm = TRUE)
    
    N <- sum(!is.na(y))   # was length(!is.na(y)), which always equals length(y)
    denominator <- N * (alpha * log(alpha) + (1 - alpha) * log(1 - alpha))
    
    h <- 1 - numerator / denominator
    
    model_name <- gsub("\"", "", insight::safe_deparse(.y), fixed = TRUE)
    data.frame(Name = model_name, Model = class(.x)[1], H = h, stringsAsFactors = FALSE)
    
  }, model_objects, object_names, SIMPLIFY = FALSE)
  
  dfs <- Reduce(function(x, y) merge(x, y, all = TRUE, sort = FALSE), m)
  
  return(dfs)
}


## epsilon-insensitive RMSE for model predictions
epsilonRMSE <- function(measured_data, predicted_data, measured_target=NULL,
                        measured_ci_low=NULL, measured_ci_high=NULL, fitted='fit'){
  
  # assign column names
  if (is.null(measured_target)){
    measured_target <- names(measured_data)[1]
  }
  
  if (is.null(measured_ci_low)){
    measured_ci_low <- names(measured_data)[2]
  }
  
  if (is.null(measured_ci_high)){
    measured_ci_high <- names(measured_data)[3]
  }
  
  stopifnot(nrow(measured_data) == nrow(predicted_data))
  
  lower <- measured_data[[measured_ci_low]]
  upper <- measured_data[[measured_ci_high]]
  pred  <- predicted_data[[fitted]]
  
  stopifnot(all(lower <= upper, na.rm = TRUE))
  
  rmse <- sqrt(mean((pred - measured_data[[measured_target]])^2, na.rm = TRUE))
  
  # calculate residuals and set to 0 for predictions within the confidence interval
  epsilon_error <- pmax(0, pmax(lower - pred, pred - upper))
  
  n <- sum(!is.na(epsilon_error))
  
  # calculate epsilon-insensitive RMSE
  epsilon_rmse <- sqrt(mean(epsilon_error^2, na.rm = TRUE))
  
  return(list(epsilon_rmse=epsilon_rmse, rmse=rmse))
}


## (Herron) estimated percentage correct predictions (ePCP) for GLMM models
glmmPCP <- function(...){
  
  model_objects <- list(...)
  
  # if (length(list(...)) == 1){
  #   
  #   model_objects <- insight::ellipsis_info(..., ..., only_models = TRUE)
  #   model_objects <- model_objects[1]
  #   
  # } else{
  #   
  #   model_objects <- insight::ellipsis_info(..., only_models = TRUE)
  #   
  # }
  
  # ensure proper object names
  dot_names <- sapply(match.call(expand.dots = FALSE)[["..."]], as.character)
  check_object_names <- insight::compact_character(names(model_objects))
  if ((is.null(check_object_names) ||
       # or if length of names doesn't match number of models
       length(check_object_names) != length(model_objects) ||
       # or if names are "..1", "..2" pattern
       all(grepl("\\.\\.\\d", check_object_names))) &&
      # and length of dot-names must match length of objects
      length(model_objects) == length(dot_names)) {
    names(model_objects) <- dot_names
  }
  
  object_names <- names(model_objects)
  
  m <- mapply(function(.x, .y) {
    y_full <- .x@resp$y
    
    n_full <- suppressWarnings(insight::n_obs(.x))
    
    pr_full <- stats::predict(.x, type = "response")
    
    pcp_full <- (sum(1 - pr_full[y_full == 0]) + sum(pr_full[y_full == 1])) / n_full
    
    dat <- data.frame(ePCP=pcp_full)
    model_name <- gsub("\"", "", insight::safe_deparse(.y), fixed = TRUE)
    perf_df <- data.frame(Name = model_name, Model = class(.x)[1], dat, stringsAsFactors = FALSE)
    
  }, model_objects, object_names, SIMPLIFY = FALSE)
  
  dfs <- Reduce(function(x, y) merge(x, y, all = TRUE, sort = FALSE), m)
  
  return(dfs)
  
}

## Cohen's f2 effect size for GEE models
geeCohensf2 <- function(geeMod){
  # loop over each variable in the model, calculate the R2 using geeRsquared
  # function for a reduced formula excluding that variable, and calculate the f2
  # effect size
  
  # first, calculate the R2 for the full model
  R2full <- geeRsquared(geeMod)$rsquare_gee
  
  # # assign variables except intercept
  # vars <- geeMod$coefficients %>%
  #   rownames() %>%
  #   .[!grepl("Intercept", .)]
  # 
  f2 <- sapply(rownames(geeMod$coefficients), function(x) {
    # create a reduced formula
    reduced_formula <- update(geeMod$formula, paste0(". ~ . - ", x))
    # fit the reduced model
    reduced_model <- update(geeMod, formula = reduced_formula)
    # calculate the R2 for the reduced model
    R2reduced <- geeRsquared(reduced_model)$rsquare_gee
    # calculate the f2 effect size
    f2 <- (R2full - R2reduced) / (1 - R2full)
    
  })
  
  return(f2)
  
}