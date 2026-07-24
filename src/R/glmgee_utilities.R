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
      
      # Known signatures of a fit that "ran" but didn't actually converge
      bad_patterns <- c("Iteration limit exceeded",
                        "not positive definite")
      if (any(sapply(bad_patterns, grepl, x = paste(warn_msgs, collapse = " ")))) {
        stop(paste("Non-convergence detected:", paste(unique(warn_msgs), collapse = "; ")))
      }
      
      fit
    }, error = function(e) {
      warning(sprintf("Model failed for corstr = '%s': %s", cs, conditionMessage(e)))
      NULL
    })
  }
  
  model_list <- Filter(Negate(is.null), model_list)
  
  # Splice list into separate ... arguments, preserving names as
  # geeModelPerformance expects (mirrors calling it with each model named)
  model_perf <- do.call(geeModelPerformance, c(model_list, list(verbose = verbose)))
  
  # replace the generic model names in the Name column with the correlation structure names
  model_perf$Name <- names(model_list)
  
  return(model_perf)
  
}



## Formatting GEE parameters
geeParams <- function(geeMod, ci=0.95, exponentiate=FALSE, varest='bias-corrected',
                      round=NULL, verbose=TRUE){
  
  geeModParams <- broom::tidy(geeMod, exponentiate=exponentiate, conf.level=ci, conf.int=TRUE, varest=varest)
  
  # Shannon information or 'surprisal' S-value with base 2 (Greenland, 2019)
  geeModParams$S2.value <- -log2(geeModParams$p.value)
  
  # Bayes factor upper bound (Benjamin & Berger, 2016) - evidence in favour of the alternative hypothesis over the null hypothesis, based on the     p-value
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
