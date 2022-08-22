#!/usr/bin/env Rscript 
# -*- coding: utf-8 -*-
# By Elena


library("optparse")
library("survival")
library("survminer") 


# Get arguments
######################################################

option_list = list(
    make_option(c("--input"), type="character", default=NULL, 
              help="Prepped input file with GrimAge + Episcore projections, and sex and age, as well as time to death/death status, for testing LBC data", metavar="character"),
    make_option(c("--cpg_bage"), type="character", default=NULL, 
              help="Option to obtain predictor directly from methylation data (T/F)", metavar="character"),
    make_option(c("--cpg_file"), type="character", default=NULL, 
              help="File location of prepped methylation (beta) data", metavar="character"),
    make_option(c("--cpg_subset"), type="character", default=NULL, 
              help="Option to obtain only use a subset of CpGs in training of predictor (T/F)", metavar="character"),
    make_option(c("--combo_bage"), type="character", default=NULL, 
              help="Option to obtain predictor from both episcores and methylation data (T/F)", metavar="character"),
    make_option(c("--pheno"), type="character", default=NULL, 
              help="Prepped input file with GrimAge + Episcore projections, and sex and age, as well as time to death/death status, for testing LBC data", metavar="character"),
    make_option(c("--training"), type="character", default=NULL, 
              help="Training data set basic info, to obtain age distribution & transform predictor to age scale", metavar="character"),
    make_option(c("--weights"), type="character", default=NULL, 
              help="Result file with weights and selected features/episcores from elnet model", metavar="character"),
    make_option(c("--name"), type="character", default=NULL, 
              help="Run name", metavar="character"),
    make_option(c("--out"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# Note that to make predictor be on age scale we need to force variance and mean to be equal to those from training data!
# Make sure that these are not sig dif between waves in GS!

output_dir <- opt$out
#output_dir <- "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/biological_age/elasticnet_models/elnet_test/w1w3/random/"

# CpGs instead of Episcores/other covars as predictors?
######################################################

cpg_bage <- opt$cpg_bage
cpg_subset <- opt$cpg_subset
if (cpg_bage == "T") {
    cpg_file <- opt$cpg_file
    # cpg_file <- "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/chronological_age/data_prep/w1w3w4/methylation_training_noscale.rds"
}
combo_bage <- opt$combo_bage

# Load training data info + model weights
######################################################

training_info <- read.delim(opt$training, row.names = 2)
weights <- read.delim(opt$weights, row.names = 1)
# training_info <- read.delim("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/chronological_age/elasticnet_models/cv_folds/w1w3/random/gs_basic_folds_random.tsv", row.names = 2) # Terminal
# weights <- read.delim("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/biological_age/elasticnet_models/elnet_train/w1w3w4/random_cpgpredictor_cpgsubset100/elnet_coefficients_random_cpgpredictor_cpgsubset100.tsv", row.names = 1) # Terminal

# Import + prep testing data
######################################################

testing_info <- read.delim(opt$pheno, row.names = 20)
if (cpg_bage == "F") {
  testing <- readRDS(opt$input)
  testing <- testing[rownames(testing_info),]
}
if ((cpg_bage == "T")) {
  testing <- readRDS(cpg_file)
  death <- readRDS(opt$input)
  #death <- death[,c("tte", "dead")]
  samples <- Reduce(intersect, list(rownames(testing), rownames(death))) # 1331/1342
  testing <- testing[samples,]
  death <- death[samples,]
  testing_info <- testing_info[samples,]
  testing <- cbind(death, testing)
}

# testing <- readRDS("/Cluster_Filespace/Marioni_Group/Elena/data/lbc_data/bage_variables.rds")
# testing_info <- read.delim("/Cluster_Filespace/Marioni_Group/Elena/data/lbc_data/lbc_testing_info.tsv", row.names = 20)

# Remove people that have missing or strange time-to-event values (negative)
testing <- testing[!is.na(testing$tte) & testing$tte>0, ]

# Filter test data
x_test <- testing[,rownames(weights)]
x_test <- t(x_test)
print(rownames(x_test))
print(dim(x_test))


# Prediction
######################################################

pred <- x_test * weights[,"Coefficient"]
pred_pp <- colSums(pred)

# Scale to same scale as age in training
scale_pred <- function(x, mean_pred, sd_pred, mean_train, sd_train) { 
  scaled <- mean_train + (x - mean_pred)*(sd_train/sd_pred)
  return(scaled)
}

mean_pred <- mean(pred_pp)
mean_train <- mean(training_info$age)
sd_pred <- sd(pred_pp)
sd_train <- sd(training_info$age)
pred_pp_scaled <- scale_pred(pred_pp, mean_pred, sd_pred, mean_train, sd_train)


# Fuse everything into df
######################################################

pred_df <- data.frame(pred_pp_scaled, testing_info[colnames(x_test),"DNAmGrimAge",drop=FALSE], testing_info[colnames(x_test),"age",drop=FALSE], testing_info[colnames(x_test),"sex",drop=FALSE], testing[colnames(x_test),"tte",drop=FALSE], testing[colnames(x_test),"dead",drop=FALSE], testing_info[colnames(x_test),"cohort"])
names(pred_df) <- c("bAge", "GrimAge", "Age", "Sex", "TTE", "Dead", "Cohort")
cat("Have prepped testing data.\n")


# Age acceleration estimate (regress age out of score)
# For both our bAge + GrimAge
######################################################

# bAge residuals
pred_df$GrimAgeAccel <- resid(lm(GrimAge ~ Age, data=pred_df, na.action=na.exclude))
pred_df$bAgeAccel <- resid(lm(bAge ~ Age, data=pred_df, na.action=na.exclude))


# Get bottom and top quantiles for GrimAge + bAge accel
# Bottom quantile 0, top quantile 1, everything else NA
######################################################

pred_df$GrimAgeAccel_lohi <- NA 
pred_df$bAgeAccel_lohi <- NA

pred_df$GrimAgeAccel_lohi[pred_df$GrimAgeAccel <= quantile(pred_df$GrimAgeAccel, 0.25)] <- 0
pred_df$GrimAgeAccel_lohi[pred_df$GrimAgeAccel >= quantile(pred_df$GrimAgeAccel, 0.75)] <- 1

pred_df$bAgeAccel_lohi[pred_df$bAgeAccel <= quantile(pred_df$bAgeAccel, 0.25)] <- 0
pred_df$bAgeAccel_lohi[pred_df$bAgeAccel >= quantile(pred_df$bAgeAccel, 0.75)] <- 1

# Per cohort
for (cohort in unique(pred_df$Cohort)) {
  col_g <- paste0("GrimAgeAccel_lohi_", cohort)
  col_b <- paste0("bAgeAccel_lohi_", cohort)
  subs <- rownames(pred_df[pred_df["Cohort"] == cohort,])
  print(length(subs))
  pred_df[col_g] <- NA
  
  pred_df[subs, col_g][pred_df[subs, "GrimAgeAccel"] <= quantile(pred_df[subs, "GrimAgeAccel"], 0.25)] <- 0
  pred_df[subs, col_g][pred_df[subs, "GrimAgeAccel"] >= quantile(pred_df[subs, "GrimAgeAccel"], 0.75)] <- 1

  pred_df[subs, col_b][pred_df[subs, "bAgeAccel"] <= quantile(pred_df[subs, "bAgeAccel"], 0.25)] <- 0
  pred_df[subs, col_b][pred_df[subs, "bAgeAccel"] >= quantile(pred_df[subs, "bAgeAccel"], 0.75)] <- 1
}


# Output pred_df table to make plots in python
######################################################

# Output table
write.table(data.frame("Sample" = rownames(pred_df), pred_df), file = paste0(output_dir, "predictions.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)


# Cox models + export results to make plots in python
######################################################

grimage_fit <- coxph(Surv(TTE, Dead) ~ Age + factor(Sex) + scale(GrimAgeAccel), data=pred_df)
coxzph_grim <- cox.zph(grimage_fit)
print(coxzph_grim)
pdf(paste0(output_dir, "coxph_assumptions_grimage.pdf"))
print(ggcoxzph(coxzph_grim))
dev.off() 
grim_stuff <- cbind(summary(grimage_fit)$coefficients, exp(confint(grimage_fit)))
rownames(grim_stuff) <- c("GrimAge_Age", "GrimAge_Sex", "GrimAge_GrimAgeAccel")
for (cohort in unique(pred_df$Cohort)) {
  print(cohort)
  grimage_fit_cohort <- coxph(Surv(TTE, Dead) ~ Age + factor(Sex) + scale(GrimAgeAccel), data=pred_df[pred_df$Cohort == cohort,])
  grim_stuff_cohort <- cbind(summary(grimage_fit_cohort)$coefficients, exp(confint(grimage_fit_cohort)))
  rownames(grim_stuff_cohort) <- c(paste0(cohort, "_GrimAge_Age"), paste0(cohort, "_GrimAge_Sex"), paste0(cohort, "_GrimAge_GrimAgeAccel"))
  # Append
  grim_stuff <- rbind(grim_stuff, grim_stuff_cohort)
  # Cox PH assumptions
  coxzph_cohort_grim <- cox.zph(grimage_fit_cohort)
  print(coxzph_cohort_grim)
  pdf(paste0(output_dir, "coxph_assumptions_grimage_", cohort, ".pdf"))
  print(ggcoxzph(coxzph_cohort_grim))
  dev.off() 
}
colnames(grim_stuff) <- c("logHR", "HR", "SE", "Z", "p", "HR_CI95_Low", "HR_CI95_High")

bage_fit <- coxph(Surv(TTE, Dead) ~ Age + factor(Sex) + scale(bAgeAccel), data=pred_df)
coxzph_bage <- cox.zph(bage_fit)
print(coxzph_bage)
pdf(paste0(output_dir, "coxph_assumptions_bage.pdf"))
print(ggcoxzph(coxzph_bage))
dev.off() 
bage_stuff <- cbind(summary(bage_fit)$coefficients, exp(confint(bage_fit)))
rownames(bage_stuff) <- c("bAge_Age", "bAge_Sex", "bAge_bAgeAccel")
for (cohort in unique(pred_df$Cohort)) {
  print(cohort)
  bage_fit_cohort <- coxph(Surv(TTE, Dead) ~ Age + factor(Sex) + scale(bAgeAccel), data=pred_df[pred_df$Cohort == cohort,])
  bage_stuff_cohort <- cbind(summary(bage_fit_cohort)$coefficients, exp(confint(bage_fit_cohort)))
  rownames(bage_stuff_cohort) <- c(paste0(cohort, "_bAge_Age"), paste0(cohort, "_bAge_Sex"), paste0(cohort, "_bAge_bAgeAccel"))
  # Append
  bage_stuff <- rbind(bage_stuff, bage_stuff_cohort)
  # Cox PH assumptions
  coxzph_cohort_bage <- cox.zph(bage_fit_cohort)
  print(coxzph_cohort_bage)
  pdf(paste0(output_dir, "coxph_assumptions_bage_", cohort, ".pdf"))
  print(ggcoxzph(coxzph_cohort_bage))
  dev.off() 
}
colnames(bage_stuff) <- c("logHR", "HR", "SE", "Z", "p", "HR_CI95_Low", "HR_CI95_High")

cox <- rbind(bage_stuff, grim_stuff)

# Export!
write.table(data.frame("Variable" = rownames(cox), cox), file = paste0(output_dir, "coxph_testing.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)

