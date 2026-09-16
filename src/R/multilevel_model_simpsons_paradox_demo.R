library(lme4)
library(ggplot2)
library(ggnewscale)
library(ciTools)
library(dplyr)

# open a file save dialog to save the plot to a selected folder
# toggle to save plots
saveplots = TRUE

if (saveplots){
  # set output plot directory
  choose.files(caption="Just cancel this", filters=matrix(data=c(" ", " "), ncol=2))  # workaround for bug in RTerm choose.dir
  outFigPath <- utils::choose.dir(caption="Select output folder to save plots")
  
  if (!dir.exists(file.path(outFigPath, "svg"))){dir.create(file.path(outFigPath, "svg"))}
  if (!dir.exists(file.path(outFigPath, "pdf"))){dir.create(file.path(outFigPath, "pdf"))}
  
}


# Adapted from https://rpubs.com/cschwarz/mlm2w1d1 
# Schwarz, C (2020). Advanced Multilevel Modeling (ICPSR MLM2): Fixed and Random Effects

# plot colour scheme

mycolourlist = list(c(0, 102, 255), c(0, 204, 153), c(255, 0, 102), c(74, 111, 152), c(251, 164, 49), c(204, 153, 255), c(90, 192, 255), c(80, 245, 233), c(255, 90, 192), c(164, 201, 242), c(255, 254, 139), c(255, 243, 255))
mycolours = matrix()

for (ii in 1:length(mycolourlist)){
  mycolours[ii] = rgb(mycolourlist[[ii]][1]/255,
                      mycolourlist[[ii]][2]/255,
                      mycolourlist[[ii]][3]/255)
}

set.seed(1234)
n_groups <- 6
nobs <- n_groups*167
shift_val <- 10

grp <- rep_len(1:n_groups, nobs)
Clusters <- factor(seq(1, nobs) %% n_groups + 1)
grp_means <- -seq(-n_groups, n_groups, length = n_groups) + rnorm(n_groups, sd = 1.2) + shift_val
Predictor <- rnorm(nobs, grp_means[grp], sd = 1)
err <- rnorm(nobs)
Outcome <- Predictor + 3.5 * grp + err

data <- data.frame(Outcome, Predictor, grp, Clusters)
true_means <- aggregate(cbind(Predictor, Outcome) ~ Clusters, data = data, FUN = mean)

# Models
re <- lmer(Outcome ~ Predictor + (1 | Clusters), data = data)
fe <- lm(Outcome ~ Predictor + Clusters, data = data)

# 1. Base grid for marginal/population level lines
marginal_grid <- data.frame(
  Predictor = seq(mean(data$Predictor) - 6, mean(data$Predictor) + 6, length.out = 50)
)

# 2. LMM Marginal Predictions
pred_grid_lmer <- ciTools::add_ci(df = marginal_grid, fit = re, alpha = 0.05, type = "parametric", includeRanef = FALSE)
pred_grid_lmer <- pred_grid_lmer |> dplyr::rename(Outcome = pred, conf.low = LCB0.025, conf.high = UCB0.975)

# 3. LM Fixed-Effect Marginal Predictions via Averaging
# Create a grid crossing all predictors with ALL clusters
fe_grid <- expand.grid(
  Predictor = marginal_grid$Predictor,
  Clusters = levels(data$Clusters)
)

# Extract standard linear predictions with confidence intervals
fe_preds <- cbind(fe_grid, as.data.frame(predict(fe, newdata = fe_grid, interval = 'confidence')))

# Mathematically average across the clusters to find the true marginal mean
pred_grid_lm <- fe_preds |>
  dplyr::group_by(Predictor) |>
  dplyr::summarise(
    Outcome = mean(fit),
    conf.low = mean(lwr),
    conf.high = mean(upr)
  )

# 4. Cluster-Specific Lines (Random Effects Tracking)
pred_list <- by(data, data$Clusters, function(sub_df) {
  data.frame(
    Predictor = seq(min(sub_df$Predictor), max(sub_df$Predictor), length.out = 50),
    Clusters = unique(sub_df$Clusters)
  )
})
pred_grid_re <- do.call(rbind, pred_list)
pred_grid_re$Outcome <- predict(re, newdata = pred_grid_re, re.form = ~(1 | Clusters))

pred_grid_fe <- do.call(rbind, pred_list)
pred_grid_fe$Outcome <- predict(fe, newdata = pred_grid_fe)

# Plotting block with fixed mapping elements
ggplot(data, aes(x = Predictor, y = Outcome, colour = Clusters)) + 
  geom_point(alpha = 0.2) + 
  theme(panel.grid = element_blank(), text = element_text(family = "serif")) +
  scale_colour_manual(values = mycolours) +
  
  # Mean Points
  geom_point(data=true_means, aes(x=Predictor, y=Outcome, shape="Cluster mean"),
             size=3, stroke=3) +
  scale_shape_manual(name = "", values = c("Cluster mean" = 3)) +
  
  # Standard OLS Line
  geom_smooth(aes(linetype = "Uncorrelated OLS mean"), method = "lm", se = FALSE, colour = "black", linewidth = 1.25) +
  
  # LMM Population Line with Confidence Ribbon
  geom_line(data = pred_grid_lmer, aes(x = Predictor, y = Outcome, linetype = "Correlated LMM mean"), 
            colour = "grey50", linewidth = 1.25, inherit.aes = FALSE) +
  
  # Unified Scales
  scale_linetype_manual(name = "Models", values = c("Uncorrelated OLS mean" = "solid",
                                                    "Correlated LMM mean" = "dashed")) +
  
  # Hierarchical Layer Ordering
  guides(colour = guide_legend(order = 1, ncol=2),
         shape = guide_legend(order = 2), 
         linetype = guide_legend(order = 3)) +
  
  # Add the random effect predictions back over the top
  geom_line(data = pred_grid_re, aes(x = Predictor, y = Outcome, colour = Clusters), linewidth = 1, alpha = 0.85) + 
  scale_colour_manual(values = mycolours) +
  
  scale_x_continuous(limits = c(0, max(data$Predictor) + 1)) +
  scale_y_continuous(limits = c(min(data$Outcome) - 1, max(data$Outcome) + 1))

filename <- "SimpsonsMultilevelDemo"

if (saveplots){
  ggsave(filename=paste(filename, ".svg", sep=""), width=8, height=4.5,
         path=file.path(outFigPath, "svg"), system_fonts = list('serif'="Times New Roman"))
  unlink(paste(filename, ".svg", sep=""))
  
  ggsave(filename=paste(filename, ".pdf", sep=""), width=8, height=4.5, path=file.path(outFigPath, "pdf"))
  unlink(paste(filename, ".pdf", sep=""))
}


# Plotting block with fixed mapping elements
ggplot(data, aes(x = Predictor, y = Outcome, colour = Clusters)) + 
  geom_point(alpha = 0.2) + 
  theme(panel.grid = element_blank(), text = element_text(family = "serif")) +
  scale_colour_manual(values = mycolours) +
  
  # Reference Points
  geom_point(data=true_means, aes(x=Predictor, y=Outcome, shape="Cluster mean"),
             size=3, stroke=3) +
  scale_shape_manual(name = "", values = c("Cluster mean" = 3)) +

  # LMM Population Line with Confidence Ribbon
  geom_ribbon(data = pred_grid_lmer, aes(x = Predictor, ymin = conf.low, ymax = conf.high, linetype = "Correlated LMM mean & CIs"), 
              inherit.aes = FALSE, alpha = 0.15, fill=NA, linewidth=0.5, colour='grey50') +
  geom_line(data = pred_grid_lmer, aes(x = Predictor, y = Outcome, linetype = "Correlated LMM mean & CIs"), 
            colour = "grey50", linewidth = 1.25, inherit.aes = FALSE) +
  
  # Fixed Effects Marginal Line with Confidence Ribbon
  geom_ribbon(data = pred_grid_lm, aes(x = Predictor, ymin = conf.low, ymax = conf.high, linetype = "Fixed cluster effect mean & CIs"), 
              inherit.aes = FALSE, alpha = 0.15, fill=NA, linewidth=0.5, colour='blue') +
  geom_line(data = pred_grid_lm, aes(x = Predictor, y = Outcome, linetype = "Fixed cluster effect mean & CIs"), 
            colour = "blue", linewidth = 1.25, inherit.aes = FALSE) +
  
  # Unified Scales
  scale_linetype_manual(name = "Models", values = c("Correlated LMM mean & CIs" = 'dashed',
                                                    "Fixed cluster effect mean & CIs" = 'dotted')) +
  # Hierarchical Layer Ordering
  guides(colour = guide_legend(order = 1, ncol=2),
         shape = guide_legend(order = 2), 
         linetype = guide_legend(order = 3)) +
  
  # Add the random effect predictions back over the top
  geom_line(data = pred_grid_re, aes(x = Predictor, y = Outcome, colour = Clusters), linewidth = 1, alpha = 0.05) + 
  geom_line(data = pred_grid_re, aes(x = Predictor, y = Outcome, colour = Clusters), linetype='dotted', linewidth = 1.5, alpha = 0.95) + 
  scale_colour_manual(values = mycolours) +
  
  scale_x_continuous(limits = c(0, max(data$Predictor) + 1)) +
  scale_y_continuous(limits = c(min(pred_grid_lmer$conf.low), max(pred_grid_lmer$conf.high)))

filename <- "SimpsonsMultilevelDemo_FEconf"

if (saveplots){
  ggsave(filename=paste(filename, ".svg", sep=""), width=8, height=4.5,
         path=file.path(outFigPath, "svg"), system_fonts = list('serif'="Times New Roman"))
  unlink(paste(filename, ".svg", sep=""))
  
  ggsave(filename=paste(filename, ".pdf", sep=""), width=8, height=4.5, path=file.path(outFigPath, "pdf"))
  unlink(paste(filename, ".pdf", sep=""))
}
