require(SelectBoost)
require(SelectBoost.beta)
require(glmnet)
require(reformulas)
require(rsample)
require(future.apply)

## Variable selection stability analysis using elastic net in SelectBoost with clustered cross-validation
varStabilityClust <- function(formula, data, cluster_id, family, resamples=50, cv_fold=5,
                              alpha=0.5, seed=808, select_threshold=0.6,
                              use.parallel=FALSE){
  
  # ensure no random effects in formula
  fixed_formula <- reformulas::nobars(formula)
  
  set.seed(seed)
  folds <- rsample::group_vfold_cv(data, group={{cluster_id}}, v=cv_fold)
  
  selection_list <- list()
  sb_matrices <- list()
  
  worker_fun <- function(split) {
    
    train_data <- rsample::analysis(split)
    
    X_train <- model.matrix(fixed_formula, data = train_data)[, -1]
    y_train <- train_data[[all.vars(fixed_formula)[1]]]
    
    # --- correlation / c0 ---
    corr_mat <- abs(cor(X_train, use = "pairwise.complete.obs"))
    corr_mat[is.na(corr_mat)] <- 0
    corr_vals <- corr_mat[upper.tri(corr_mat)]
    q_vals <- quantile(corr_vals, probs = c(0.9, 1))
    
    steps.seq <- seq(mean(q_vals), max(q_vals), by = 0.05)
    if (length(steps.seq) <= 3) {
      steps.seq <- seq(mean(q_vals), max(q_vals), length.out = 5)
    }
    
    # --- SelectBoost ---
    if (grepl("beta", family$family, ignore.case = TRUE)) {
      
      squeeze <- identical(family$family, "ordbeta")
      
      sb_matrix.raw <- SelectBoost.beta::sb_beta(
        X = X_train,
        Y = y_train,
        selector = betareg_glmnet,
        B = resamples,
        steps.seq = steps.seq,
        squeeze = squeeze,
        use.parallel = FALSE,
        alpha = alpha,
        choose = "cv",
        nfolds = cv_fold
      )
      
      sb_matrix <- SelectBoost::force.non.inc(sb_matrix.raw)
      
    } else {
      
      sb_matrix.raw <- SelectBoost::fastboost(
        X = X_train,
        Y = y_train,
        func = function(X, Y) glmnet_selector(
          X, Y,
          family = family,
          alpha = alpha,
          nfolds = cv_fold
        ),
        B = resamples,
        steps.seq = steps.seq,
        use.parallel = FALSE
      )
      
      sb_matrix <- SelectBoost::force.non.inc(sb_matrix.raw)
      
    }
    
    c0_vals <- as.numeric(sub("c0 = ", "", rownames(sb_matrix)))
    
    sb_valid <- sb_matrix[!is.na(c0_vals), , drop = FALSE]
    c0_valid <- as.numeric(sub("c0 = ", "", rownames(sb_valid)))
    
    print("c0_valid:")
    print(c0_valid)
    
    confidence <- sapply(seq_len(ncol(sb_valid)), function(j) {
      x <- sb_valid[, j]
      
      selected_idx <- which(x > 0)
      
      print("selected_idx:")
      print(selected_idx)
      
      if (length(selected_idx) == 0) return(0)
      
      1 - min(c0_valid[selected_idx])
    })
    
    names(confidence) <- colnames(sb_valid)
    
    print("confidence:")
    print(confidence)
    
    sel_freq <- colMeans(sb_matrix)
    
    print("sel_freq:")
    print(sel_freq)
    
    list(sel_freq = sel_freq, sb_matrix = sb_matrix, confidence = confidence)
  }
  
  if (use.parallel) {
    
    results <- future.apply::future_lapply(
      folds$splits,
      worker_fun,
      future.packages = c("SelectBoost", "SelectBoost.beta", "glmnet", "reformulas"),
      future.seed = TRUE
    )
    
  } else {
    
    results <- lapply(
      folds$splits,
      worker_fun
    )
  }
  
  selection_list <- lapply(results, `[[`, "sel_freq")
  confidence_list <- lapply(results, `[[`, "confidence")
  
  print("selection_list:")
  print(selection_list)
  print("confidence_list:")
  print(confidence_list)
  
  sb_matrices   <- lapply(results, `[[`, "sb_matrix")
  
  # Align variables across folds
  all_vars <- unique(unlist(lapply(selection_list, names)))
  
  sel_mat <- sapply(selection_list, function(x) {
    out <- setNames(rep(NA, length(all_vars)), all_vars)
    out[names(x)] <- x
    out
  })
  
  conf_mat <- sapply(confidence_list, function(x) {
    out <- setNames(rep(NA, length(all_vars)), all_vars)
    out[names(x)] <- x
    out
  })
  
  print("sel_mat:")
  print(sel_mat)
  print("conf_mat:")
  print(conf_mat)
  
  
  # Aggregate
  sel_freq_mean <- rowMeans(sel_mat, na.rm = TRUE)
  sel_freq_mean <- sort(sel_freq_mean, decreasing = TRUE)
  confidence_mean <- rowMeans(conf_mat, na.rm = TRUE)
  confidence_mean <- sort(confidence_mean, decreasing = TRUE)
  
  selected_names <- names(confidence_mean)[
    (confidence_mean >= 0.3)
  ]
  
  main_effects <- unique(gsub(":(.*)", "", selected_names))
  
  final_vars <- unique(c(selected_names, main_effects))
  
  return(list(
    sel_freq_mean = sel_freq_mean,
    confidence_mean = confidence_mean,
    sel_mat = sel_mat,
    sb_matrices = sb_matrices,
    selected_names = selected_names,
    final_vars = final_vars
  ))
}


### Elastic net selector function

glmnet_selector <- function(x, y, family = "binomial", alpha = 0.5, nfolds=5, ...) {
  
  fit <- glmnet::cv.glmnet(
    x, y,
    family = family,
    alpha = alpha,
    nfolds = nfolds
  )
  
  coef_vec <- as.vector(coef(fit, s = "lambda.min"))
  names(coef_vec) <- rownames(coef(fit))
  
  return(coef_vec)
}