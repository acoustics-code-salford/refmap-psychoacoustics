require(glmtoolbox)
require(broom)
require(dplyr)
require(ggplot2)

## Compare GEE correlation structures by model performance
geeCorStr <- function(formula, data, family, id,
                      corstr = c("Exchangeable", "AR-M-dependent", "Stationary-M-dependent",
                                 "Independence", "Unstructured", "Non-Stationary-M-dependent"),
                      m_lag = 1, ..., verbose = FALSE) {
  
  # Resolve id and inject it into data under a fixed name, so glmgee's
  # internal eval(id, envir = data, ...) finds it regardless of environment
  id_vec <- eval(substitute(id), data, parent.frame())
  data <- data
  data[[".geeCorStr_id"]] <- id_vec
  
  corstr_choices <- match.arg(corstr, several.ok = TRUE)
  corstr_formatted <- ifelse(grepl("M-dependent", corstr_choices),
                             paste0(corstr_choices, "(", m_lag, ")"),
                             corstr_choices)
  names(corstr_formatted) <- corstr_choices
  
  model_list <- setNames(vector("list", length(corstr_formatted)), names(corstr_formatted))
  
  for (i in seq_along(corstr_formatted)) {
    cs <- corstr_formatted[[i]]
    
    model_list[[i]] <- tryCatch({
      warn_msgs <- character(0)
      
      fit <- withCallingHandlers({
        glmtoolbox::glmgee(formula = formula, data = data, family = family,
                           id = .geeCorStr_id, corstr = cs, ...)
      }, warning = function(w) {
        warn_msgs <<- c(warn_msgs, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
      
      bad_patterns <- c("Iteration limit exceeded", "not positive definite")
      if (any(sapply(bad_patterns, grepl, x = paste(warn_msgs, collapse = " ")))) {
        stop(paste("Non-convergence detected:", paste(unique(warn_msgs), collapse = "; ")))
      }
      
      # Tag with the true family name directly, since we already know it here.
      # Avoids both (a) inspecting .x$call$family post-hoc, which is fragile
      # once family is passed as a variable, and (b) do.call(), which was
      # found to change how glmgee internally stores $y (vector vs matrix)
      # and broke performance:: metrics downstream for true binomial fits.
      attr(fit, "gee_family") <- family$family
      fit
    }, error = function(e) {
      warning(sprintf("Model failed for corstr = '%s': %s", cs, conditionMessage(e)))
      NULL
    })
  }
  
  model_list <- Filter(Negate(is.null), model_list)
  
  if (length(model_list) == 0) {
    stop("All correlation structures failed to fit — no models available to compare. ",
         "Check the warnings above for the specific failure reasons.")
  }
  
  # Splice list into separate ... arguments, preserving names as
  # geeModelPerformance expects (mirrors calling it with each model named)
  model_perf <- do.call(geeModelPerformance, c(model_list, list(verbose = verbose)))
  
  # replace the generic model names in the Name column with the correlation structure names
  model_perf$Name <- names(model_list)
  
  return(model_perf)
  
}


## Compiling GEE GLM model performance parameters
geeModelPerformance <- function(..., verbose=FALSE){
  
  model_objects <- list(...)
  
  # if (length(list(...)) == 1){
  #   
  #   model_objects <- insight::ellipsis_info(..., ..., only_models = TRUE,
  #                                           verbose=verbose)
  #   model_objects <- model_objects[1]
  #   
  # } else{
  #   
  #   model_objects <- insight::ellipsis_info(..., only_models = TRUE,
  #                                           verbose=verbose)
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
  
  # function to get the correct family name
  get_family_name <- function(x) {
    # Preferred: an explicit tag set by geeCorStr at fit time. Sidesteps
    # call/family-object inspection entirely, which is fragile once family
    # is passed as a variable inside a wrapper function.
    tagged <- attr(x, "gee_family")
    if (!is.null(tagged)) return(tagged)
    
    fam <- x$call$family
    if (is.list(fam) && !is.null(fam$family)) {
      return(fam$family)
    }
    val <- tryCatch(eval(fam, envir = environment(x$formula)), error = function(e) NULL)
    if (is.list(val) && !is.null(val$family)) {
      return(val$family)
    }
    x$family$family
  }
  
  m <- mapply(function(.x, .y) {
    fam_name <- get_family_name(.x)
    if (fam_name %in% c("gaussian", "quasibinomial")) {
      r2 <- geeRsquared(.x)
      R2 <- r2$R2
      H <- NA
      D <- NA
      AUC <- NA
      Lloss <- NA
      PCP <- data.frame(pcp_model=NA, model_ci_low=NA, model_ci_high=NA, lrt_chisq=NA)
      
    } else if (fam_name == "binomial") {
      
      if (!is.null(.x$varest)) {
        varest <- .x$varest
      } else {
        varest <- "robust"
      }
      
      n_pos <- sum(.x$y == 1)
      n_neg <- sum(.x$y == 0)
      r <- rank(.x$fitted.values)
      auc <- (sum(r[.x$y == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
      
      r2 <- geeRsquared(.x)
      R2 <- r2$R2
      h <- geeEntropy(.x, varest=varest)
      H <- h$H
      D <- performance::r2_tjur(.x)
      AUC <- auc
      Lloss <- performance::performance_logloss(.x)
      PCP <- performance::performance_pcp(.x)
      
    } else {
      R2 <- NA
      H <- NA
      D <- NA
      AUC <- NA
      Lloss <- NA
      PCP <- data.frame(pcp_model=NA, model_ci_low=NA, model_ci_high=NA, lrt_chisq=NA)
      
    }
    dat <- data.frame(R2=R2,
                      RMSE = sqrt(mean((.x$y - .x$fitted.values)^2, na.rm = TRUE)),
                      MSE  = mean((.x$y - .x$fitted.values)^2, na.rm = TRUE),
                      MAE  = mean(abs(.x$y - .x$fitted.values), na.rm = TRUE),
                      H=H,
                      D=D,
                      GPAIC=glmtoolbox::AGPC(.x, verbose=verbose),
                      GPSBIC=glmtoolbox::SGPC(.x, verbose=verbose),
                      AUC=AUC,
                      Lloss=Lloss,
                      ePCP=PCP$pcp_model,
                      ePCP95Lo=PCP$model_ci_low,
                      ePCP95Hi=PCP$model_ci_high,
                      ePCPChiSq=PCP$lrt_chisq
    )
    # omit NA columns
    dat <- Filter(function(x)!all(is.na(x)), dat)
    
    model_name <- gsub("\"", "", insight::safe_deparse(.y), fixed = TRUE)
    perf_df <- data.frame(Name = model_name, Model = class(.x)[1], dat, stringsAsFactors = FALSE)
    
  }, model_objects, object_names, SIMPLIFY = FALSE)
  
  # merge dataframes
  dfs <- dplyr::tibble(Reduce(function(x, y) merge(x, y, all = TRUE, sort = FALSE), m))
  
  # add AIC and SBIC weights
  wGPAIC <- exp(-0.5*(dfs$GPAIC - min(dfs$GPAIC)))/sum(exp(-0.5*(dfs$GPAIC - min(dfs$GPAIC))))
  wGPSBIC <- exp(-0.5*(dfs$GPSBIC - min(dfs$GPSBIC)))/sum(exp(-0.5*(dfs$GPSBIC - min(dfs$GPSBIC))))
  
  dfs$wGPAIC <- format(round(wGPAIC, 6), scientific=FALSE, format=f, digits=6)
  dfs$wGPSBIC <- format(round(wGPSBIC, 6), scientific=FALSE, format=f, digits=6)
  
  # calculate mean of AIC and SBIC weights for each model in the frame (row-wise)
  wICMean <- rowMeans(cbind(wGPAIC, wGPSBIC), na.rm = TRUE)
  dfs$wICMean <- format(round(wICMean, 6), scientific=FALSE, format=f, digits=6)
  
  # move columns
  dfs <- dfs |> dplyr::relocate(wGPAIC, .after=GPAIC)
  dfs <- dfs |> dplyr::relocate(wGPSBIC, .after=GPSBIC)
  dfs <- dfs |> dplyr::relocate(wICMean, .after=wGPSBIC)
  
  # print output if verbose input argument is true
  if (verbose){
    print(dfs)
  }
  
  return(dfs)
}


## Formatting GEE parameters
geeParams <- function(geeMod, ci=0.95, exponentiate=FALSE, varest='bias-corrected',
                      round=NULL, verbose=TRUE){
  
  geeModParams <- broom::tidy(geeMod, exponentiate=exponentiate, conf.level=ci, conf.int=TRUE, varest=varest)
  
  # Shannon information or 'surprisal' S-value with base 2 (Greenland, 2019)
  geeModParams$S2.value <- -log2(geeModParams$p.value)
  
  # Bayes factor upper bound (Benjamin & Berger, 2016) - evidence in favour of the alternative hypothesis over the null hypothesis, based on the p-value
  geeModParams$Bfb <- 1/(-exp(1)*geeModParams$p.value*log(geeModParams$p.value))
  
  # Posterior probability from Bayes factor upper bound - probability that the alternative hypothesis is true given the data, assuming a prior probability of 0.5 for the alternative hypothesis
  geeModParams$Pp <- geeModParams$Bfb/(1 + geeModParams$Bfb)
  
  if (!is.null(round)){
    
    parOut <- dplyr::tibble(cbind(geeModParams[, names(geeModParams)=='term'],
                            round(geeModParams[, names(geeModParams)!='term'], round)))
    
  } else{
    
    parOut <- geeModParams
    
  }
  
  if (!is.null(geeMod$param_map)) {
    parOut$term <- geeMod$param_map[parOut$term]
  }
  
  if (verbose){
    print(parOut)
  }
  
  return(parOut)
  
}


## glmgee ggplot2 cluster-specific (positive) variable plot
geeCluster_stemplot <- function(glmgeeMod, clusterVar, colour='blue', alpha=1){
  # check the lengths of clusterVar and ids match
  if (length(clusterVar) != length(unique(glmgeeMod$ids))){
    stop("Cluster variable length does not match number of clusters")
  }
  
  # ensure cluster indexes are factors
  ids_ordered <- factor(glmgeeMod$ids, levels=sort(unique(glmgeeMod$model$`(id)`)))
  
  df <- data.frame(ID=ids_ordered, values=as.vector(clusterVar))
  
  df <- df |> dplyr::arrange(ID)
  
  ggplot(data=df) +
    
    geom_segment(aes(x=ID, xend=ID, yend=values), y=0, 
                 colour=colour, alpha=alpha) +
    
    geom_point(aes(x=ID, y=values), colour=colour, size=2)
  
}



## glmgee ggplot2 residuals plot
gee_residplot <- function(glmgeeMod, residual_vals, colour='blue', shape=1, size=2, alpha=1){
  
  # check the lengths of residuals and fitted values match
  if (length(residual_vals) != length(glmgeeMod$fitted.values)){
    stop("Residuals length does not match fitted values length")
  }
  
  ylimits <- c(min(-3.5, min(residual_vals)), max(+3.5, max(residual_vals)))
  
  df <- data.frame(fitted_vals=as.vector(glmgeeMod$fitted.values), resid_vals=as.vector(residual_vals))
  
  ggplot(data=df) +
    
    geom_point(aes(x=fitted_vals, y=resid_vals), colour=colour, shape=shape, size=size, alpha=alpha) +
    
    scale_y_continuous(limits=ylimits, breaks=seq(-20, 20, by=1),
                       labels=seq(-20, 20, by=1)) +
    
    geom_hline(yintercept=3, linetype="dotted", color = "black") +
    geom_hline(yintercept=-3, linetype="dotted", color = "black")
  
  
}
