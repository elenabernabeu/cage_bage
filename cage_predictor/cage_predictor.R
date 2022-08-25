#!/usr/bin/env Rscript 
# -*- coding: utf-8 -*-

# NOTE: Working directory should contain model weights files in data/ subdir, if not alter script accordingly
# NOTE 2: All outputs will be stored in working directory
# Required packages
library(tidyverse)


###### (0) User input
#########################################################################################################
#########################################################################################################

#methylationTable <- "" # Path to methylation table. Each column corresponds to an individual. Each row corresponds to CpG. First column is the CpG name. If RDS, assumed first column is rownames. 
#methylationTable_format <- "" # rds, tsv, or csv
methylationTable <- "/Cluster_Filespace/Marioni_Group/Elena/data/geo_data/test_methylation.RDS"
methylationTable_format <- "rds"


###### (1) Data loading
#########################################################################################################
#########################################################################################################

message("1. Loading data") 
message("1.1 Loading methylation data - rows to be CpGs and columns to be individuals") 

## Loading in model coefficients. Make sure these files are present in the current working directory or change the paths to the correct directory.
coefficients <- read.delim("./data/elnet_coefficients_linear.tsv", sep = "\t")
coefficients_log <- read.delim("./data/elnet_coefficients_log.tsv", sep = "\t")
row.names(coefficients) <- coefficients[,"CpG_Site"]
row.names(coefficients_log) <- coefficients_log[,"CpG_Site"]

# Beta means
means <- read.delim("./data/cpg_meanbeta_gs20k.tsv")
row.names(means) <- means[,"cpg"]

# Intercepts
intercept <- read.table("./data/intercept_linear.txt", sep = "")[1, "V1"]
intercept_log <- read.table("./data/intercept_log.txt", sep = "")[1, "V1"]

# Coefficients for log models
coef_log_2 <- coefficients_log[rownames(coefficients_log)[grep('_2', rownames(coefficients_log))],,drop=FALSE]
coef_log_2_simp <- gsub('_2', '', rownames(coef_log_2))
coef_log <- coefficients_log[which(!(rownames(coefficients_log) %in% rownames(coef_log_2))),,drop=FALSE]

# Coefficients for non-log models
coef_2 <- coefficients[rownames(coefficients)[grep('_2', rownames(coefficients))],,drop=FALSE]
coef_2_simp <- gsub('_2', '', rownames(coef_2))
coef <- coefficients[rownames(coefficients)[which(!(rownames(coefficients) %in% rownames(coef_2)))],,drop=FALSE]

# Total CpGs
cpgs_linear <- union(rownames(coef), rownames(coef_log))
cpgs_squared <- union(coef_2_simp, coef_log_2_simp)
all_cpgs <- union(cpgs_linear, cpgs_squared)

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


###### (2) QC and data prep
#########################################################################################################
#########################################################################################################

## Check if Data needs to be Transposed
message("2. Quality Control and data preparation") 
message("2.1 Checking if row names are CpG sites") 

if(ncol(data) > nrow(data)){
  message("It seems that individuals are rows - data will be transposed!")
  data <- t(data) 
}

## Subset CpG sites to those present on list for predictors 
message("2.2 Subsetting CpG sites to those required for predictor calculation") 
coef_data <- data[intersect(rownames(data), all_cpgs),]

## Check if Beta or M Values
message("2.3 Checking of Beta or M-values are present") 
m_to_beta <- function(val) {
  beta <- 2^val/(2^val + 1)
  return(beta)
}

coef_data <- if((range(coef_data, na.rm = T) > 1)[[2]] == "TRUE") {
    message("Suspect that M Values are present. Converting to Beta Values")
    m_to_beta(coef_data)
  } else {
    message("Suspect that Beta Values are present");
    coef_data
  }

## Identify CpGs missing from input dataframe, include them and provide values as mean methylation value at that site
message("2.4 Find CpGs not present in uploaded file, add these with mean Beta Value for CpG site from training sample") 
coef_data <- if (nrow(coef_data)==length(all_cpgs)) { 
  message("All sites present")
  coef_data
  } else if (nrow(coef_data)==0) { 
  message("There Are No Necessary CpGs in The dataset - All Individuals Would Have Same Values For Predictors. Analysis Is Not Informative!")
  } else {
    missing_cpgs <- all_cpgs[which(!(all_cpgs %in% rownames(coef_data)))]
    message(paste(length(missing_cpgs), "unique sites are missing - adding to dataset with mean Beta Value from GS training sample (N = 18,413)", sep = " "))
    mat = matrix(nrow=length(missing_cpgs), ncol = ncol(coef_data))
    row.names(mat) <- missing_cpgs
    colnames(mat) <- colnames(coef_data)
    mat[is.na(mat)] <- 1
  ids <- unique(row.names(mat))
  missing_cpg_means <- means[ids,]
  mat <- mat*missing_cpg_means[,"mean"]
  coef_data <- rbind(coef_data,mat)
} 

message("2.5 Convert NA Values to mean for each probe") 
## Convert NAs to Mean Value for all individuals across each probe 
na_to_mean <-function(methyl) {
  methyl[is.na(methyl)] <- mean(methyl, na.rm=T)
  return(methyl)
}

coef_data <- t(apply(coef_data,1,function(x) na_to_mean(x)))


###### (3) cAge prediction
#########################################################################################################
#########################################################################################################

message("4. Obtaining cAge predictions") 
message("4.1. Preparing data for prediction") 

## Prep for linear predictor
scores <- coef_data[rownames(coef),]
scores_quadratic <- coef_data[coef_2_simp,]**2
rownames(scores_quadratic) <- paste0(rownames(scores_quadratic), "_2")
scores_linear <- rbind(scores, scores_quadratic)

## Prep for log predictor
scores_log <- coef_data[rownames(coef_log),]
scores_quadratic_log <- coef_data[coef_log_2_simp,]**2
rownames(scores_quadratic_log) <- paste0(rownames(scores_quadratic_log), "_2")
scores_log <- rbind(scores_log, scores_quadratic_log)

## Calculate cAge with linear model
message("4.2. Calculating cAge using model trained on linear age") 
coefficients <- coefficients[rownames(scores_linear),]
pred_linear <- scores_linear * coefficients[,"Coefficient"]
pred_linear_pp <- colSums(pred_linear)
pred_linear_pp <- pred_linear_pp + intercept

## Identify any individuals predicted as under 20s, and re-run with model trained on log(age)
over20s <- names(pred_linear_pp[pred_linear_pp > 20])
pred_linear_pp <- pred_linear_pp[over20s]
under20s <- names(pred_linear_pp[pred_linear_pp < 20])

## Now re-run model for those
coefficients_log <- coefficients_log[rownames(scores_log),]
pred_log <- scores_log * coefficients_log[,"Coefficient"]
pred_log_pp <- colSums(pred_log)
pred_log_pp <- pred_log_pp + intercept_log
pred_log_pp <- exp(pred_log_pp[under20s])


###### (5) Export results
#########################################################################################################
#########################################################################################################

message("Exporting predictions to working directory!") 
predictions <- c(pred_log_pp, pred_linear_pp)
results <- data.frame("Sample" = names(predictions), "Predicted_Age" = predictions)
write.table(results, file = "cage_predictions.tsv", quote = FALSE, sep = "\t", row.names = FALSE)