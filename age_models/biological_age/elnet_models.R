#!/usr/bin/env Rscript 
# -*- coding: utf-8 -*-
# By Elena


library("optparse")
library("glmnet")
library("survival")
library("coxme")


# Get arguments
######################################################

option_list = list(
    make_option(c("--input"), type="character", default=NULL, 
              help="Prepped input file with GrimAge + Episcore projections, and sex and age, as well as time to death/death", metavar="character"),
    make_option(c("--cpg_bage"), type="character", default=NULL, 
              help="Option to obtain predictor directly from methylation data (T/F)", metavar="character"),
    make_option(c("--cpg_file"), type="character", default=NULL, 
              help="File location of prepped methylation (beta) data", metavar="character"),
    make_option(c("--cpg_subset"), type="character", default=NULL, 
              help="Option to obtain only use a subset of CpGs in training of predictor (T/F)", metavar="character"),
    make_option(c("--cpg_subset_file"), type="character", default=NULL, 
              help="File location of CpG subset to be used in training of predictor", metavar="character"),
    make_option(c("--combo_bage"), type="character", default=NULL, 
              help="Option to obtain predictor from both episcores and methylation data (T/F)", metavar="character"),
    make_option(c("--agecause"), type="character", default=NULL, 
              help="File with age at death and cause of death", metavar="character"),
    make_option(c("--folds"), type="character", default=NULL, 
              help="File with folds assigned", metavar="character"),
    make_option(c("--sex"), type="character", default=NULL, 
              help="Sex-stratified option (T/F)", metavar="character"),
    make_option(c("--squared"), type="character", default=NULL, 
              help="Non-linear exploration option (T/F)", metavar="character"),
    make_option(c("--name"), type="character", default=NULL, 
              help="Run name", metavar="character"),
    make_option(c("--kinship"), type="character", default=NULL, 
              help="Run LMM with kinship matrix option", metavar="character"),
    make_option(c("--kinship_file"), type="character", default=NULL, 
              help="Location of kinship matrix", metavar="character"),
    make_option(c("--deathage"), type="character", default=NULL, 
              help="If not null, option to limit training to deaths over given age", metavar="character"),
    make_option(c("--deathcause"), type="character", default=NULL, 
              help="If not null, option to limit training to certain causes for death (listed ICD10 codes separated by commas)", metavar="character"),
    make_option(c("--out"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


set.seed(1234) # Set seed to ensure fold variation minimised 
seed <- 1234


# Sex stratified models?
######################################################

sex_s <- opt$sex


# Include non-linear exploration?
######################################################

epi_s <- opt$squared


# CpGs instead of Episcores/other covars as predictors?
######################################################

cpg_bage <- opt$cpg_bage
if (cpg_bage == "T") {
    cpg_file <- opt$cpg_file
    # cpg_file <- "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/chronological_age/data_prep/w1w3w4/methylation_training_noscale.rds"
}
combo_bage <- opt$combo_bage

# CpGs subsets?
######################################################

cpg_subset <- opt$cpg_subset
if (cpg_subset == "T") {
    cpg_subset_file <- opt$cpg_subset_file
    # cpg_subset_file <- "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/survival_ewas/coxph/subsets/coxph_survival_ewas_3000.tsv"
}


# Train just on given age subset or cause of deaths?
######################################################

deathage <- opt$deathage
if (opt$deathage != "F") {
    deathage <- as.numeric(opt$deathage)
}

deathcause <- opt$deathcause
if (opt$deathcause != "F") {
    deathcause <- unlist(strsplit(opt$deathcause, ","))
}


# Kinship matrix in model (LMM)?
######################################################

kin <- opt$kinship
kin_f <- opt$kinship_file
#kin <- "T"
#kin_f <- "/Cluster_Filespace/Marioni_Group/Elena/data/gs_data/gs_kinship_20k.RDS"


# Import data and prep to run
######################################################

folds <- read.delim(opt$folds, header = TRUE, row.names = 2)
# folds <- read.delim("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/chronological_age/elasticnet_models/cv_folds/w1w3w4/random/gs_basic_folds_random.tsv", header = TRUE, row.names = 2)
if (cpg_bage == "F") {
    df <- read.table(opt$input, header = TRUE, check.names = FALSE, row.names = 1)
    df <- df[match(rownames(folds), rownames(df)),] # Match rownames
} else {
    death <- read.table(opt$input, header = TRUE, check.names = FALSE, row.names = 1)
    # death <- read.table("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/biological_age/data_prep/gs_w1w3w4_bage_variables_scaledgrim.tsv", header = TRUE, check.names = FALSE, row.names = 1)
    if (combo_bage == "F") {
        death <- death[,c("tte", "dead", "Sample_Name")]
    } 
    df <- readRDS(cpg_file)
    death <- death[match(rownames(folds), rownames(death)),] # Match rownames
    df <- df[match(rownames(folds), rownames(df)),] # Match rownames
    df <- cbind(death, df)
}

# Import causes and ages of death
agecause <- read.delim(opt$agecause, row.names = 2)

# Limit to over given age if selected
if (deathage != "F") {
    samples_dead <- rownames(agecause[(agecause$age_event >= deathage) & (agecause$dead == 1),])
    samples_alive <- rownames(agecause[agecause$dead == 0,])
    samples <- c(samples_dead, samples_alive)
    df <- df[rownames(df) %in% samples,]
    cat(paste("\nTotal number of deaths remaining in training: ", length(rownames(df[df$dead == 1,]))))
    cat(paste("\nTotal number of individuals remaining in training: ", length(rownames(df))))
} else {
    cat(paste("\nTotal number of deaths in training: ", length(rownames(df[df$dead == 1,]))))
    cat(paste("\nTotal number of individuals in training: ", length(rownames(df))))
}


# Limit to certain causes of death if selected
#if (causeage != "F") {
#    samples_dead <- rownames(agecause[(agecause$cause %in% deathcause) & (agecause$dead == 1),])
#    samples_alive <- rownames(agecause[agecause$dead == 0,])
#    samples <- c(samples_dead, samples_alive)
#    df <- df[rownames(df) %in% samples,]
#    cat(paste("\nTotal number of deaths remaining in training: ", length(rownames(df[df$dead == 1,]))))
#    cat(paste("\nTotal number of individuals remaining in training: ", length(rownames(df))))
#}

# Limit CpGs if option selected
if ((cpg_bage == "T") & (cpg_subset == "T")) {
    sub_file <- read.table(cpg_subset_file, header = TRUE, row.names = 1)
    # sub_file <- read.table("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/survival_ewas/coxph/subsets/coxph_survival_ewas_3000.tsv", header = TRUE, row.names = 1)
    df <- df[,(colnames(df) %in% rownames(sub_file)) | (colnames(df) %in% colnames(death))]
}



# Remove people that have missing or strange time-to-event values (negative)
df <- df[!is.na(df$tte) & df$tte>0, ]
cat(paste0("\nFinal df dimensions: ", dim(df)))

# Number of deaths per fold
folds <- folds[rownames(df),]
folds["dead"] <- df["dead"]
for (fold in 1:length(unique(folds$Fold))) {
    samples <- folds[folds["Fold"] == fold,]
    samples <- rownames(samples[samples["dead"] == 1,])
    cat(paste0("\nNumber of deaths in fold ", fold, ": ", length(samples)))
}


folds <- folds$Fold
cat("\nImported and prepped data!\n")


# Fit Cox models
######################################################

if (kin == "F") {
    cat("Fitting elastic net model (Cox PH) for time to death (not accounting for relatedness)...\n")
    # Variables to keep: age, sex, grimage components, and 109 episcores (everything but dead status and tte, first two columns)
    x <- as.matrix(df[4:ncol(df)])
    y <- Surv(df$tte,df$dead) # Time to event, and wether the event has happened or not
    # Obtain lambda
    cv <- cv.glmnet(x, y, family = "cox", type.measure = "C", seed = seed, foldid = folds, nfolds = length(unique(folds)), trace.it = 1, thresh = 1e-5) # Harrell's C index to obtain best parameters (kind of like residuals to minimize in OLS)
    lambda <- cv$lambda.min
    # Obtain coefs
    fit <- glmnet(x, y, family = "cox", trace.it = 1, thresh = 1e-5)
    #fit <- glmnet(x, y, family = "cox", trace.it = 1, thresh = 1e-5)
    results <- as.data.frame(cbind("lambda" = fit$lambda, t(as.matrix(fit$beta))))

    cat("Now extracting info of interest.\n")
    #coefs <- coef(fit) # Get betas
    coefs <- as.data.frame(t(results[results$lambda == lambda,][1,]))
    coefs <- coefs[2:nrow(coefs),,drop=FALSE]
    colnames(coefs) <- "Coefficient" 
    coefs <- coefs[coefs$Coefficient !=0,,drop=FALSE] # Remove coeficients that are 0 
    coefs["Variable"] <- rownames(coefs)
    coefs <- coefs[,c("Variable", "Coefficient")]
} else { # YET TO FINISH
    cat("Fitting elastic net model (Cox PH) for time to death (with kinship matrix)...\n")
    kinship_matrix <- readRDS(kin_f)
    # Cox with Kinship matrix
    cox <- coxme(Surv(df$tte, df$dead) ~ (1|df$id), varlist = kinship_matrix*2)
    ran_effecs <- ranef(cox)
}


# Export 
######################################################

cat("Exporting!\n")
filename <- paste0(opt$out, "elnet_coefficients_", opt$name, ".tsv")
#filename <- paste0("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/biological_age/elasticnet_models/elnet_train/w1w3/random/", "elnet_coefficients_random", ".tsv")
write.table(coefs, file = filename, row.names = FALSE, sep = "\t", quote = FALSE)
