
# import libraries
library(naniar)
library(corrplot)
library(ggplot2)
library(reshape2)


windowsFonts()
par(family="serif")

# directory paths
dirpath = "C:\\Users\\m_lot\\OneDrive - University of Salford\\REFMAP General\\03 Experiment\\Experiment 1\\Analysis\\"
infilepath = paste(dirpath, "PostProcess", sep="")
outplotpath = paste(dirpath, "Plots", sep="")

# load Part A data
dfA <- read.csv("C:\\Users\\m_lot\\OneDrive - University of Salford\\REFMAP General\\03 Experiment\\Experiment 1\\Analysis\\PostProcess\\refmap_listest1_testdataA_BySubj.csv", header=TRUE)
dfB <- read.csv("C:\\Users\\m_lot\\OneDrive - University of Salford\\REFMAP General\\03 Experiment\\Experiment 1\\Analysis\\PostProcess\\refmap_listest1_testdataB_BySubj.csv", header=TRUE)

# form subset of data for correlation analysis
dropcols <- c("SessionPart", "ID.", "TrialNumber", "StimFile",
              "UAS_noticed", "Recording", "StimDuration",
              "UASOperation", "UASType", "AmbientEnv",
              "Hearing_impaired", "Sex", "NationGeo",
              "NativeLang", "Home.Area", "Area.soundscape",
              "ResponseFactorLoud",	"ResponseFactorChar",
              "ResponseFactorDuration",	"ResponseFactorQuantRep",
              "ResponseFactorProximity",	"ResponseFactorAmb",
              "ResponseFactorExper", "ResponseFactorSens",
              "ResponseFactorComfort",	"ResponseFactorPrivSafe",
              "AAMExperience",	"AAM.attitude")

subdfA <- dfA[, !(names(dfA) %in% dropcols)]
subdfB <- dfB[, !(names(dfB) %in% dropcols)]

# replace "No UAS" and "No answer entries with NA
subdfA <- subdfA |> replace_with_na(replace=list(UASLAeq="No UAS", SNRlevel="No UAS", Age="No answer"))
subdfB <- subdfB |> replace_with_na(replace=list(UASLAeq="No UAS", SNRlevel="No UAS", Age="No answer"))


# convert character variables to numeric
cols2num <- c("UASLAeq", "SNRlevel", "Age")
subdfA[cols2num] <- sapply(subdfA[cols2num], as.numeric)
subdfB[cols2num] <- sapply(subdfB[cols2num], as.numeric)

# correlation matrices
pearsonA <- cor(subdfA, method = "pearson", use="complete.obs")
spearmanA <- cor(subdfA, method = "spearman", use="complete.obs")
kendallA <- cor(subdfA, method = "kendall", use="complete.obs")
pearsonB <- cor(subdfB, method = "pearson", use="complete.obs")
spearmanB <- cor(subdfB, method = "spearman", use="complete.obs")
kendallB <- cor(subdfB, method = "kendall", use="complete.obs")

# plotting
# reverse colour palette
redblueRev <- colorRampPalette(c('#053061', '#2166AC', '#4393C3',
                                 '#92C5DE', '#D1E5F0', '#FFFFFF',
                                 '#FDDBC7', '#F4A582', '#D6604D',
                                 '#B2182B', '#67001F'))

pdf(file=paste(outplotpath, "PtACorrPearson.pdf", sep="\\"), width=48, height=30)
corrplot(pearsonA, type='upper', col=redblueRev(100), addCoef.col='black', number.cex = 0.5,
         method='ellipse', tl.cex=1, na.label=" ", diag=FALSE, family="serif")
dev.off()

pdf(file=paste(outplotpath, "PtACorrSpearman.pdf", sep="\\"), width=48, height=30)
corrplot(spearmanA, type='upper', col=redblueRev(100), addCoef.col='black', number.cex = 0.5,
         method='ellipse', tl.cex=1, na.label=" ", diag=FALSE, family="serif")
dev.off()

pdf(file=paste(outplotpath, "PtACorrKendall.pdf", sep="\\"), width=48, height=30)
corrplot(kendallA, type='upper', col=redblueRev(100), addCoef.col='black', number.cex = 0.5,
         method='ellipse', tl.cex=1, na.label=" ", diag=FALSE, family="serif")
dev.off()

pdf(file=paste(outplotpath, "PtBCorrPearson.pdf", sep="\\"), width=48, height=30)
corrplot(pearsonB, type='upper', col=redblueRev(100), addCoef.col='black', number.cex = 0.5,
         method='ellipse', tl.cex=1, na.label=" ", diag=FALSE, family="serif")
dev.off()

pdf(file=paste(outplotpath, "PtBCorrSpearman.pdf", sep="\\"), width=48, height=30)
corrplot(spearmanB, type='upper', col=redblueRev(100), addCoef.col='black', number.cex = 0.5,
         method='ellipse', tl.cex=1, na.label=" ", diag=FALSE, family="serif")
dev.off()

pdf(file=paste(outplotpath, "PtBCorrKendall.pdf", sep="\\"), width=48, height=30)
corrplot(kendallB, type='upper', col=redblueRev(100), addCoef.col='black', number.cex = 0.5,
         method='ellipse', tl.cex=1, na.label=" ", diag=FALSE, family="serif")
dev.off()

# plot bar charts of response correlations
row.names(data.frame(pearsonA[!(names(data.frame(pearsonA))
                                %in% c("Annoyance")),"Annoyance"]))


library(party)
#library(partykit)
#require(parallel)
require(flexplot)
library(randomForest)
library(tidyr)
set.seed(1)
nCores <- detectCores()
clust <- makeCluster(nCores)
data(avengers)
flexplot(ptsd~1, data=avengers)
rf_model <- cforest(ptsd~., data=avengers)#, applyfun=parLapply(clust))
estimates(rf_model)

stimDataA <- read.csv("C:\\Users\\m_lot\\OneDrive - University of Salford\\REFMAP General\\03 Experiment\\Experiment 1\\Analysis\\PostProcess\\refmap_listest1_testdataA_ByStim.csv", header=TRUE)

encode_ordinal <- function(x, order = unique(x)) {
  x <- as.numeric(factor(x, levels = order, exclude = NULL))
  x
}

SNRCats <- c("No UAS", "-16", "-10", "-4", "2", "8")
UASLAeqCats <- c("No UAS", "42", "48", "54", "60")

stimDataANum <- stimDataA 

stimDataANum[['SNRlevel']] <- encode_ordinal(stimDataA[['SNRlevel']], order=SNRCats)
stimDataANum[['UASLAeq']] <- encode_ordinal(stimDataA[['UASLAeq']], order=UASLAeqCats)

stimDataANum <- cbind(stimDataANum[, which(colnames(stimDataANum)=="UASLAeq"):which(colnames(stimDataANum)=="SNRlevel")],
                      stimDataANum[, which(colnames(stimDataANum)=="UASLAeqMaxLR"):which(colnames(stimDataANum)=="UASImpulsHMSMaxMaxLR")],
                      stimDataANum[, which(colnames(stimDataANum)=="UASLAeqdiffAmbLAF90"):which(colnames(stimDataANum)=="UASPartLoudGMSTPowAvg")])
colSelect = c('ValenceMedian', 'ArousalMedian', 'AnnoyMedian', 'HighAnnoyProportion')
stimDataANum <- cbind(stimDataANum, stimDataA[, colSelect])
stimDataANum <- subset(stimDataANum, select = -c(UASLAEMaxLR))
includeVars <- names(stimDataANum)[1:39]
formVar <- reformulate(includeVars, 'AnnoyMedian')
crf_model <- cforest(formVar, data=stimDataANum,
                     controls=cforest_unbiased(ntree=40000, mtry=as.integer(length(includeVars)/3)))
results <- estimates(crf_model)
results
impCrf <- varimp(crf_model, nperm=10)
impCRFOrdered <- round(sort(impCrf, decreasing=F), digits=4)
par(mai=c(1,3,1,1))
barplot(height=impCRFOrdered, names.arg=rownames(impCRFOrdered),
        horiz=TRUE, las=1)

X = drop_na(stimDataANum[, includeVars])
y = stimDataANum[rowSums(is.na(stimDataANum), 2) == 0, 'AnnoyMedian']

rf_model <- randomForest(x=X, y=y, ntree=20000,
                         mtry=as.integer(length(includeVars)/3),
                         importance=TRUE, nPerm=5, localImp=T, proximity = F)
rf_model$localImportance
imp <- as.data.frame(rf_model$importance)
impMSEOrdered <- imp[order(imp$`%IncMSE`, decreasing = F),]
impPurOrdered <- imp[order(imp$IncNodePurity, decreasing = F),]

par(mai=c(1,2,1,1))
barplot(height=impMSEOrdered$`%IncMSE`, names.arg=rownames(impMSEOrdered),
        horiz=TRUE, las=1)

par(mai=c(1,2,1,1))
barplot(height=impPurOrdered$IncNodePurity, names.arg=rownames(impPurOrdered),
        horiz=TRUE, las=1)
