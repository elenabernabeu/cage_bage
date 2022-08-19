#!/usr/bin/env Rscript 
# -*- coding: utf-8 -*-

# NOTE: Working directory should contain model weights files, if not alter script accordingly
# NOTE 2: All outputs will be stored in working directory
# Required packages
library(tidyverse)
library(survival)

###### (0) User input
#########################################################################################################
#########################################################################################################

methylationTable <- "/Cluster_Filespace/Marioni_Group/Elena/data/lbc_data/lbc_bvals_cagetesting.rds" # Path to methylation table. Each column corresponds to an individual. Each row corresponds to CpG. First column is the CpG name. If RDS, assumed first column is rownames. 
methylationTable_format <- "rds" # rds, tsv, or csv

phenotypeTable <- "/Cluster_Filespace/Marioni_Group/Elena/data/lbc_data/lbc_testing_info.tsv" # Path to phenotype table. Each column corresponds to a variable. Each row corresponds to an individual/sample. First column is sample name. If RDS, assumed first column is rownames. Should include DNAmAge.
phenotypeTable_format <- "rds" # rds, tsv, or csv

grimageTable <- "/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/GrimAge_output.csv" # Path to GrimAge results table. Each column represents a GrimAge component. Each row corresponds to an individual/sample. First column is sample name.
grimageTable_format <- "csv" # rds, tsv, or csv

# Insert column names for each of the variables corresponding to the names in phenotypeTable
ageColname <- "Age" # Age in years
sexColname <- "Female" # Sex variable is binary, 1 = females, 0 = males
tteColname <- "tte" # Time to event (death) column
deathColname <- "dead" # Death variable is binary, 1 = dead, 0 = alive

###### (1) Data loading
#########################################################################################################
#########################################################################################################

message("1. Loading data") 
message("1.1 Loading methylation data - rows to be CpGs and columns to be individuals") 

## Loading in model coefficients. Make sure these files are present in the current working directory or change the paths to the correct directory.
coefficients <- read.delim("bage_coefficients.tsv")

## Loading in CpG coefficients for episcore projection
cpgs <- read.delim("cpg_episcore_weights.tsv")

## Loading methylation data
if (tolower(methylationTable_format) == "rds") {
  data <- readRDS(methylationTable)
} else if (tolower(methylationTable_format) == "tsv") {
  data <- read.delim(methylationTable, sep = "\t", row.names = 1)
} else if (tolower(methylationTable_format) == "csv") {
  data <- read.csv(methylationTable, sep = ",", row.names = 1)
} else {
  message("Unrecognized methylation data format. Accepted formats: rds, tsv, and csv.")
}

## Loading phenotype data
if (tolower(phenotypeTable_format) == "rds") {
  pheno <- readRDS(phenotypeTable)
} else if (tolower(phenotypeTable_format) == "tsv") {
  pheno <- read.delim(phenotypeTable, sep = "\t", row.names = 1)
} else if (tolower(phenotypeTable_format) == "csv") {
  pheno <- read.csv(phenotypeTable, sep = ",", row.names = 1)
} else {
  message("Unrecognized phenotype data format. Accepted formats: rds, tsv, and csv.")
}

## Loading GrimAge data
if (tolower(grimageTable_format) == "rds") {
  grim <- readRDS(grimageTable)
} else if (tolower(grimageTable_format) == "tsv") {
  grim <- read.delim(grimageTable, sep = "\t", row.names = 2)
} else if (tolower(grimageTable_format) == "csv") {
  grim <- read.csv(grimageTable, sep = ",", row.names = 2)
} else {
  message("Unrecognized GrimAge data format. Accepted formats: rds, tsv, and csv.")
}


###### (2) QC and data prep for Episcore calculation
#########################################################################################################
#########################################################################################################

## Check if Data needs to be Transposed
message("2. Quality Control and data preparation") 
message("2.1 Checking if row names are CpG sites") 

if(ncol(data) > nrow(data)){
  message("It seems that individuals are rows - data will be transposed!")
  data <- t(data) 
}

## Keep just individuals with non-weird TTEs
message("2.2 Removing individuals with missing TTEs") 
pheno <- pheno[!is.na(pheno[tteColname]) & pheno[tteColname]>0,]

## Keep just individuals in pheno table
message("2.3 Subsetting samples to those in phenotype table") 
data <- data[,which(colnames(data) %in% rownames(pheno))]
grim <- grim[which(rownames(grim) %in% rownames(pheno)),]

## Subset CpG sites to those present on list for predictors 
message("2.4 Subsetting CpG sites to those required for predictor calculation") 
coef <- data[intersect(rownames(data), cpgs$CpG_Site),]

## Check if Beta or M Values
message("2.5 Checking of Beta or M-values are present") 
m_to_beta <- function(val) {
  beta <- 2^val/(2^val + 1)
  return(beta)
}

coef <- if((range(coef, na.rm = T) > 1)[[2]] == "TRUE") {
    message("Suspect that M Values are present. Converting to Beta Values")
    m_to_beta(coef)
  } else {
    message("Suspect that Beta Values are present");
    coef
  }

## Scale data if needed
message("2.6 Scaling data (if needed)") 
ids <- colnames(coef)
scaled <- apply(coef, 1, function(x) sd(x, na.rm = T)) 

coef <-  if(range(scaled)[1] == 1 & range(scaled)[2] == 1) { 
    coef
  } else { 
    coef_scale <- apply(coef, 1, scale)
    coef_scale <- t(coef_scale)
    coef_scale <- as.data.frame(coef_scale)
    colnames(coef_scale) <- ids
    coef_scale
  }

## Identify CpGs missing from input dataframe, include them and provide values as mean methylation value at that site
message("2.7 Find CpGs not present in uploaded file, add these with mean Beta Value for CpG site from training sample") 
coef <- if (nrow(coef)==length(unique(cpgs$CpG_Site))) { 
  message("All sites present")
  coef
  } else if (nrow(coef)==0) { 
  message("There Are No Necessary CpGs in The dataset - All Individuals Would Have Same Values For Predictors. Analysis Is Not Informative!")
  } else { 
    missing_cpgs = cpgs[-which(cpgs$CpG_Site %in% rownames(coef)), c("CpG_Site", "Mean_Beta_Value")]
    message(paste(length(unique(missing_cpgs$CpG_Site)), "unique sites are missing - adding to dataset with mean Beta Value from training sample", sep = " "))
    mat = matrix(nrow=length(unique(missing_cpgs$CpG_Site)), ncol = ncol(coef))
    row.names(mat) <- unique(missing_cpgs$CpG_Site)
    colnames(mat) <- colnames(coef) 
    mat[is.na(mat)] <- 1
    missing_cpgs1 <- if (length(which(duplicated(missing_cpgs$CpG_Site))) > 1) { 
      missing_cpgs[-which(duplicated(missing_cpgs$CpG_Site)),]
    } else {
      missing_cpgs
    }  
  ids = unique(row.names(mat))
  missing_cpgs1 = missing_cpgs1[match(ids,missing_cpgs1$CpG_Site),]
  mat = mat*missing_cpgs1$Mean_Beta_Value
  coef = rbind(coef,mat)
} 

message("2.8 Convert NA Values to mean for each probe") 
## Convert NAs to Mean Value for all individuals across each probe 
na_to_mean <-function(methyl) {
  methyl[is.na(methyl)] <- mean(methyl, na.rm=T)
  return(methyl)
}

coef <- t(apply(coef,1,function(x) na_to_mean(x)))


###### (3) Calculate Episcores 
#########################################################################################################
#########################################################################################################

message("3. Calculating Episcores") 
loop <- unique(cpgs$Predictor)
out <- data.frame()
for(i in loop){ 
  tmp=coef[intersect(row.names(coef),cpgs[cpgs$Predictor %in% i,"CpG_Site"]),]
  tmp_coef = cpgs[cpgs$Predictor %in% i, ]
  if(nrow(tmp_coef) > 1) { 
    tmp_coef = tmp_coef[match(row.names(tmp),tmp_coef$CpG_Site),]
    out[colnames(coef),i]=colSums(tmp_coef$Coefficient*tmp)
  } else {
    tmp2 = as.matrix(tmp)*tmp_coef$Coefficient 
    out[colnames(coef),i] = tmp2[,1]
  }
} 

## Save file
message("3.1. Exporting Episcores")  
write.table(out, "episcore_projections.tsv", sep = "\t", quote = FALSE)


###### (4) Scale GrimAge components
#########################################################################################################
#########################################################################################################

message("4. Prepping GrimAge components") 

samples <- rownames(out)
grim_pred <- grim[samples, c("DNAmGrimAge"), drop = FALSE]
grim <- grim[samples, c("DNAmADM", "DNAmB2M", "DNAmCystatinC", "DNAmGDF15", "DNAmLeptin", "DNAmPACKYRS", "DNAmPAI1", "DNAmTIMP1")]

message("4.1. Scale GrimAge components (if needed)") 
scaled_grim <- apply(grim, 2, function(x) sd(x, na.rm = T)) 
ids <- colnames(grim)

grim <-  if(range(scaled_grim)[1] == 1 & range(scaled_grim)[2] == 1) { 
    grim
  } else { 
    grim_scale <- scale(grim)
    grim_scale <- as.data.frame(grim_scale)
    grim_scale
  }


###### (5) bAge prediction
#########################################################################################################
#########################################################################################################

message("5. Obtaining bAge predictions") 
## Change pheno table columns
names(pheno)[names(pheno) == ageColname] <- "Age"
names(pheno)[names(pheno) == tteColname] <- "TTE"
names(pheno)[names(pheno) == deathColname] <- "Dead"
names(pheno)[names(pheno) == sexColname] <- "Sex"

## Fuse all Episcores and GrimAge components, along with other stuff, for each sample
scores <- cbind(pheno[samples, c("TTE", "Dead", "Age", "Sex")], grim, out)

## Filter to elements in predictor
scores <- scores[, coefficients$Variable]

## Calculate bAge
message("5.1. Calculating bAge") 
scores <- t(scores)
pred <- scores * coefficients[,"Coefficient"]
pred_pp <- colSums(pred)

## Scale to same scale as age in training
message("5.2. Scaling bAge") 
scale_pred <- function(x, mean_pred, sd_pred, mean_train, sd_train) { 
  scaled <- mean_train + (x - mean_pred)*(sd_train/sd_pred)
  return(scaled)
}

mean_pred <- mean(pred_pp)
mean_train <- 47.5 # Mean age in training data
sd_pred <- sd(pred_pp)
sd_train <- 14.9 # SD age in training data
pred_pp_scaled <- scale_pred(pred_pp, mean_pred, sd_pred, mean_train, sd_train)

## Make df with everything
pred_df <- data.frame(pred_pp_scaled, grim_pred, pheno[samples, c("Age", "Sex", "TTE", "Dead")])
names(pred_df) <- c("bAge", "GrimAge", "Age", "Sex", "TTE", "Dead")

## Obtain bAgeAccel and GrimAgeAccel
message("5.3. Obtaining bAgeAccel") 
pred_df$GrimAgeAccel <- resid(lm(GrimAge ~ Age, data=pred_df, na.action=na.exclude))
pred_df$bAgeAccel <- resid(lm(bAge ~ Age, data=pred_df, na.action=na.exclude))

## Export
message("5.4. Exporting predictions") 
write.table(data.frame("Sample" = rownames(pred_df), pred_df), file = paste0("predictions.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)


###### (6) bAge vs GrimAge prediction
#########################################################################################################
#########################################################################################################

message("6. Assessing predictive ability of bAge") 
## CoxPH models, compared with GrimAge
grimage_fit <- coxph(Surv(TTE, Dead) ~ Age + factor(Sex) + scale(GrimAgeAccel), data=pred_df)
grim_stuff <- cbind(summary(grimage_fit)$coefficients, exp(confint(grimage_fit)))
rownames(grim_stuff) <- c("GrimAge_Age", "GrimAge_Sex", "GrimAge_GrimAgeAccel")
for (cohort in unique(pred_df$Cohort)) {
  grimage_fit_cohort <- coxph(Surv(TTE, Dead) ~ Age + factor(Sex) + scale(GrimAgeAccel), data=pred_df[pred_df$Cohort == cohort,])
  grim_stuff_cohort <- cbind(summary(grimage_fit_cohort)$coefficients, exp(confint(grimage_fit_cohort)))
  rownames(grim_stuff_cohort) <- c(paste0(cohort, "_GrimAge_Age"), paste0(cohort, "_GrimAge_Sex"), paste0(cohort, "_GrimAge_GrimAgeAccel"))
  # Append
  grim_stuff <- rbind(grim_stuff, grim_stuff_cohort)
}
colnames(grim_stuff) <- c("logHR", "HR", "SE", "Z", "p", "HR_CI95_Low", "HR_CI95_High")

bage_fit <- coxph(Surv(TTE, Dead) ~ Age + factor(Sex) + scale(bAgeAccel), data=pred_df)
bage_stuff <- cbind(summary(bage_fit)$coefficients, exp(confint(bage_fit)))
rownames(bage_stuff) <- c("bAge_Age", "bAge_Sex", "bAge_bAgeAccel")
for (cohort in unique(pred_df$Cohort)) {
  bage_fit_cohort <- coxph(Surv(TTE, Dead) ~ Age + factor(Sex) + scale(bAgeAccel), data=pred_df[pred_df$Cohort == cohort,])
  bage_stuff_cohort <- cbind(summary(bage_fit_cohort)$coefficients, exp(confint(bage_fit_cohort)))
  rownames(bage_stuff_cohort) <- c(paste0(cohort, "_bAge_Age"), paste0(cohort, "_bAge_Sex"), paste0(cohort, "_bAge_bAgeAccel"))
  # Append
  bage_stuff <- rbind(bage_stuff, bage_stuff_cohort)
}
colnames(bage_stuff) <- c("logHR", "HR", "SE", "Z", "p", "HR_CI95_Low", "HR_CI95_High")
cox <- rbind(bage_stuff, grim_stuff)

message("6.1. Exporting Cox PH results") 
write.table(data.frame("Variable" = rownames(cox), cox), file = paste0("coxph_testing.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)

