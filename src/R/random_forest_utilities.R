# random_forest_tools.R

require(party)
require(permimp)
require(caret)
require(cowplot)
require(tidyverse)

# Multiple iteration aggregated conditional permutation importance random forest regression function -----
multi_crfReg <- function(dataIn, iVars, dVar, seeds, ntree, mtry, permImpCondThres=0.95, minsplit=20, minbucket=7, nperm=1){
  # initialise variables
  crfOOBErrAll <- 0
  crfOOBRMSE <- 0
  crfOOBMAE <- 0
  crfOOBErrR2 <- 0
  crfMarPermImpVals <- 0
  crfConPermImpVals <- 0
  crfMarPermImpValsPerTree <- data.frame()
  crfConPermImpValsPerTree <- data.frame()
  
  for (iters in 1:length(seeds)){
    
    # formula for regression
    formVars <- reformulate(iVars, dVar)
    
    # set random seed
    set.seed(seeds[iters])
    # train crf model
    crfModel <- party::cforest(formVars, data=dataIn,
                               controls=party::cforest_unbiased(ntree=ntree,
                                                                mtry=mtry,
                                                                minsplit=minsplit,
                                                                minbucket=minbucket))
    
    # get OOB predictions
    crfModelOOB <- predict(crfModel, OOB=TRUE, type='response')
    
    # get OOB error
    crfModelOOBErr <- as.numeric(as.matrix(as.numeric(as.matrix(crfModelOOB))
                                           - as.numeric(as.matrix(crfModel@data@env$response[[names(crfModel@data@env$response)]]))))
    
    # OOB RMSE, MAE and Rsquared
    crfOOBRMSE <- crfOOBRMSE + sqrt(mean(crfModelOOBErr^2))
    crfOOBMAE <- crfOOBMAE + mean(abs(crfModelOOBErr))
    crfOOBErrR2 <- crfOOBErrR2 + cor(as.numeric(as.matrix(crfModelOOB)),
                                     as.numeric(as.matrix(crfModel@data@env$response[[names(crfModel@data@env$response)]])))^2
    
    # set random seed
    set.seed(seeds[iters])
    
    # set random seed
    set.seed(seeds[iters])
    # conditional variable permutation importance
    crfConPermImp <- permimp::permimp(crfModel, nperm=nperm, conditional=TRUE, threshold=permImpCondThres, progressBar=FALSE)
    
    crfConPermImpVals <- crfConPermImpVals + crfConPermImp$values
    crfConPermImpValsPerTree <- rbind(crfConPermImpValsPerTree, crfConPermImp$perTree)
  }
  
  # average metrics
  crfOOBErrAll <- crfOOBErrAll/length(seeds)
  crfOOBRMSE <- crfOOBRMSE/length(seeds)
  crfOOBMAE <- crfOOBMAE/length(seeds)
  crfOOBErrR2 <- crfOOBErrR2/length(seeds)
  crfConPermImpVals <- data.frame(CondPermImp=sort(crfConPermImpVals/length(seeds), decreasing=TRUE))
  crfConPermImpValsQtl <- data.frame(apply(crfConPermImpValsPerTree, 2, quantile, probs=c(0.25, 0.50, 0.75)))
  
  resultsOut <- list('OOB_RMSE'=crfOOBRMSE, 'OOB_MAE'=crfOOBMAE, 'Rsquared'=crfOOBErrR2, 'conditional_permimp'=crfConPermImpVals,                      'conditional_permimp_perTree'=crfConPermImpValsPerTree, 'conditional_permimp_qtl'=crfConPermImpValsQtl)
  return(resultsOut)
}


# Comparing random forest rankings from two seeds  -----
crfReg <- function(dataIn, iVars, dVar, seeds, ntree, mtry, permImpCondThres=0.95, minsplit=20, minbucket=7, nperm=1){
  # initialise variables
  crfOOBErrAll <- 0
  crfOOBRMSE <- 0
  crfOOBMAE <- 0
  crfOOBErrR2 <- 0
  crfMarPermImpVals <- 0
  crfConPermImpVals <- 0
  crfMarPermImpValsPerTree <- data.frame()
  crfConPermImpValsPerTree <- data.frame()
  
  # formula for regression
  formVars <- reformulate(iVars, dVar)
  
  for (iters in 1:length(seeds)){
    
    # set random seed
    set.seed(seeds[iters])
    # train crf model
    crfModel <- party::cforest(formVars, data=dataIn,
                               controls=party::cforest_unbiased(ntree=ntree,
                                                                mtry=mtry,
                                                                minsplit=minsplit,
                                                                minbucket=minbucket))
    
    # conditional variable permutation importance
    crfConPermImp <- permimp::permimp(crfModel, nperm=nperm, conditional=TRUE, threshold=permImpCondThres, progressBar=FALSE)
    
    crfConPermImpVals <- crfConPermImp$values
    
    if (iters == 1){
      crfConPermImpVals1 <- data.frame(CondPermImp=sort(crfConPermImpVals, decreasing=TRUE))
      crfConPermImpValsPerTree1 <- crfConPermImp$perTree
      crfConPermImpValsQtl1 <- data.frame(apply(crfConPermImpValsPerTree1, 2, quantile, probs=c(0.25, 0.50, 0.75)))
      
      # get OOB predictions
      crfModelOOB <- predict(crfModel, OOB=TRUE, type='response')
      
      # get OOB error
      crfModelOOBErr <- as.numeric(as.matrix(as.numeric(as.matrix(crfModelOOB))
                                             - as.numeric(as.matrix(crfModel@data@env$response[[names(crfModel@data@env$response)]]))))
      
      # OOB RMSE, error quartiles and Rsquared
      crfOOBRMSE <- sqrt(mean(crfModelOOBErr^2))
      crfOOBMAE <- crfOOBMAE + mean(abs(crfModelOOBErr))
      crfOOBErrR2 <- cor(as.numeric(as.matrix(crfModelOOB)),
                         as.numeric(as.matrix(crfModel@data@env$response[[names(crfModel@data@env$response)]])))^2
      
    }
    
    else{
      crfConPermImpValsN <- data.frame(CondPermImp=sort(crfConPermImpVals, decreasing=TRUE))
      
      nVarImpChecks <- min(max(sum(crfConPermImpVals1 >= crfConPermImpVals1$CondPermImp[1]*0.1),
                               sum(crfConPermImpValsN >= crfConPermImpValsN$CondPermImp[1]*0.1)), 4)  # number of variable importance values with a value at least 10% of the highest importance
      if (sum(rownames(crfConPermImpVals1)[1:nVarImpChecks] != rownames(crfConPermImpValsN)[1:nVarImpChecks]) > 0){
        warning("Permutation importance rank order within 10% of max differs between seeds: increase number of trees ('ntree') or number of permutations ('nperm'), or subsample of features ('mtry')")
      }
      else{
        resultsOut <- list('OOB_errors'=crfModelOOBErr, 'OOB_RMSE'=crfOOBRMSE, 'OOB_MAE'=crfOOBMAE, 'Rsquared'=crfOOBErrR2, 'conditional_permimp'=crfConPermImpVals1, 'conditional_permimp_perTree'=crfConPermImpValsPerTree1, 'conditional_permpimp_qtl'=crfConPermImpValsQtl1)
        return(resultsOut)
      }
      
    }
    
  }
  
}

# Hyperparameter tuning and plotting -------

mtryTune <- function(dataIn, iVars, dVar, seeds, ntrees, minsplit=20, minbucket=7){
  
  formVars <- reformulate(iVars, dVar)
  
  # set mtry values and corresponding iVars/mtry ratios
  if (length(iVars) > 9){
    iVars_mtrys <- c(10.5, 5.25, 3.5, 2.75, 2.25, 1.75, 1.5, 1.25)
    mtrys <- as.integer(length(iVars)/iVars_mtrys)
  }
  else{
    mtrys <- seq(2, length(iVars) - 3, by=1)
    iVars_mtrys <- length(iVars)/mtrys
  }
  iVars_mtrys <- iVars_mtrys[mtrys >= 2]  # remove 0 or 1 values
  mtrys <- mtrys[mtrys >= 2]  # remove 0 or 1 values
  
  # remove any duplicated values
  iVars_mtrys <- iVars_mtrys[!(duplicated(mtrys) | duplicated(mtrys, fromLast = TRUE))]
  mtrys <- mtrys[!(duplicated(mtrys) | duplicated(mtrys, fromLast = TRUE))]
  
  # ensure mtrys is less than length(iVars)
  iVars_mtrys <- iVars_mtrys[mtrys < length(iVars)]
  mtrys <- mtrys[mtrys < length(iVars)]
  
  resRMSEMap = matrix(data=0, nrow=length(mtrys), ncol=length(ntrees))
  resRsquaredMap = matrix(data=0, nrow=length(mtrys), ncol=length(ntrees))
  resMAEMap = matrix(data=0, nrow=length(mtrys), ncol=length(ntrees))
  
  
  for (ii in seq(1, length(ntrees))){
    
    tuneMod.results <- data.frame(RMSE=numeric(length(mtrys)),
                                  Rsquared=numeric(length(mtrys)),
                                  MAE=numeric(length(mtrys)))
    
    for (seed in seeds){
      set.seed(seed)
      ntree = ntrees[ii]
      tuneMod <- caret::train(formVars,
                              data=dataIn,
                              method='cforest',
                              controls=party::cforest_unbiased(ntree=ntree,
                                                               minsplit=minsplit,
                                                               minbucket=minbucket),
                              tuneGrid=data.frame(.mtry=mtrys),
                              trControl = trainControl(method = "oob"))
      
      
      
      # accumulate results
      resRMSEMap[, ii] <- resRMSEMap[, ii] + tuneMod$results$RMSE
      resRsquaredMap[, ii] <- resRsquaredMap[, ii] + tuneMod$results$Rsquared
      resMAEMap[, ii] <- resMAEMap[, ii] + tuneMod$results$MAE
      
      tuneMod.results <- tuneMod.results + tuneMod$results[, which(names(tuneMod$results) != 'mtry')]
    }
    
    # average results
    tuneMod.results <- tuneMod.results/length(seeds)
    tuneMod.results <- cbind(tuneMod.results, data.frame(mtry=mtrys), data.frame(iVars_mtry=iVars_mtrys))
    
    print(tuneMod.results)
    
  }
  
  # average results
  resRMSEMap <- resRMSEMap/length(seeds)
  resRsquaredMap <- resRsquaredMap/length(seeds)
  resMAEMap <- resMAEMap/length(seeds)
  
  
  # convert to data frames with mtry as row names and ntree as column names, and convert to long format using tidyverse
  resdfRMSEMap <- as.data.frame(resRMSEMap)
  rownames(resdfRMSEMap) <- mtrys
  colnames(resdfRMSEMap) <- ntrees
  resdfRsquaredMap <- as.data.frame(resRsquaredMap)
  rownames(resdfRsquaredMap) <- mtrys
  colnames(resdfRsquaredMap) <- ntrees
  resdfMAEMap <- as.data.frame(resMAEMap)
  rownames(resdfMAEMap) <- mtrys
  colnames(resdfMAEMap) <- ntrees
  
  
  # convert dataframes to long format using tidyverse
  resdfRMSEMap <- resdfRMSEMap |>
    rownames_to_column('mtry') |>
    gather(key='ntree', value='RMSE', -mtry)
  
  resdfRsquaredMap <- resdfRsquaredMap |>
    rownames_to_column('mtry') |>
    gather(key='ntree', value='Rsquared', -mtry)
  
  resdfMAEMap <- resdfMAEMap |>
    rownames_to_column('mtry') |>
    gather(key='ntree', value='MAE', -mtry)
  
  # ensure ntree and mtry columns are ordered factors
  resdfRMSEMap$ntree <- factor(resdfRMSEMap$ntree, levels=as.character(ntrees))
  resdfRMSEMap$mtry <- factor(resdfRMSEMap$mtry, levels=as.character(mtrys))
  
  resdfRsquaredMap$ntree <- factor(resdfRsquaredMap$ntree, levels=as.character(ntrees))
  resdfRsquaredMap$mtry <- factor(resdfRsquaredMap$mtry, levels=as.character(mtrys))
  
  resdfMAEMap$ntree <- factor(resdfMAEMap$ntree, levels=as.character(ntrees))
  resdfMAEMap$mtry <- factor(resdfMAEMap$mtry, levels=as.character(mtrys))
  
  # plot heatmaps using ggplot, with extreme (min or max) value plotted as overlaid point using annotate and colourbar scale reversed
  pHeatmapRMSE <- ggplot(resdfRMSEMap) +
    geom_tile(aes(x=ntree, y=mtry, fill=RMSE)) +
    scale_fill_viridis(option="viridis", direction=-1) +
    geom_point(data=resdfRMSEMap[which(resdfRMSEMap$RMSE
                                       == min(resdfRMSEMap$RMSE),
                                       arr.ind = TRUE),],
               aes(x=ntree, y=mtry), colour="red", size=2) +
    guides(colour = guide_colourbar(reverse=TRUE)) +
    labs(x="ntree", y="mtry", fill="RMSE") +
    theme(text = element_text(family = "serif"))
  
  pHeatmapRsquared <- ggplot(resdfRsquaredMap) +
    geom_tile(aes(x=ntree, y=mtry, fill=Rsquared)) +
    scale_fill_viridis(option="viridis", direction=1) +
    geom_point(data=resdfRsquaredMap[which(resdfRsquaredMap$Rsquared
                                           == max(resdfRsquaredMap$Rsquared),
                                           arr.ind = TRUE),],
               aes(x=ntree, y=mtry), colour="red", size=2) +
    guides(colour = guide_colourbar(reverse=TRUE)) +
    labs(x="ntree", y="mtry", fill="Rsquared") +
    theme(text = element_text(family = "serif"))
  
  pHeatmapMAE <- ggplot(resdfMAEMap) +
    geom_tile(aes(x=ntree, y=mtry, fill=MAE)) +
    scale_fill_viridis(option="viridis", direction=-1) +
    geom_point(data=resdfMAEMap[which(resdfMAEMap$MAE
                                      == min(resdfMAEMap$MAE),
                                      arr.ind = TRUE),],
               aes(x=ntree, y=mtry), colour="red", size=2) +
    guides(colour = guide_colourbar(reverse=TRUE)) +
    labs(x="ntree", y="mtry", fill="MAE") +
    theme(text = element_text(family = "serif"))
  
  p <-  cowplot::plot_grid(pHeatmapRMSE,
                           pHeatmapRsquared,
                           pHeatmapMAE,
                           ncol=3, nrow=1)
  
  return(p)
  
}  # end of function