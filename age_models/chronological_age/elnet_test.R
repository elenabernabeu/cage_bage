#!/usr/bin/env Rscript 
# -*- coding: utf-8 -*-
# By Elena


library("optparse")
library("stringr")
library("imputeTS")

# Get arguments
######################################################

option_list = list(
    make_option(c("--meth"), type="character", default=NULL, 
              help="Methylation file for test data", metavar="character"),
    make_option(c("--pheno"), type="character", default=NULL, 
              help="Pheno file with age and sex for test data", metavar="character"),
    make_option(c("--weights_dir"), type="character", default=NULL, 
              help="Directory where file with weights and selected CpGs from elnet model is located", metavar="character"),
    make_option(c("--name"), type="character", default=NULL, 
              help="Run name", metavar="character"),
    make_option(c("--out"), type="character", default=NULL, 
              help="Output directory", metavar="character"),
    make_option(c("--cpg"), type="character", default=NULL, 
              help="Scale option per CpG site (T/F)", metavar="character"),
    make_option(c("--sample"), type="character", default=NULL, 
              help="Scale option per sample (T/F)", metavar="character"),
    make_option(c("--test_cohort"), type="character", default=NULL, 
              help="Cohorts to include in testing, comma delimited ('all' to include all)", metavar="character"),
    make_option(c("--logage"), type="character", default=NULL, 
              help="Option to run models with log(age) instead of age (T/F)", metavar="character"),
    make_option(c("--logage20"), type="character", default=NULL, 
              help="Option to run models with log(age) instead of age for under 20s (Y), over 20s (O), or none/all (F) based on --logage flag.", metavar="character"),
    make_option(c("--sex"), type="character", default=NULL, 
              help="Sex stratified option (T/F)", metavar="character"),
    make_option(c("--squared"), type="character", default=NULL, 
              help="CpG squared option (T/F)", metavar="character"),
    make_option(c("--squaredsubset"), type="character", default=NULL, 
              help="CpG squared subset option (T/F)", metavar="character"),
    make_option(c("--squaredsubset_n"), type="character", default=NULL, 
              help="Number of squared CpGs that were subsetted on in training (from age^2 EWAS OR CpG^2 lm models)", metavar="character"),
    make_option(c("--squaredcpg2"), type="character", default=NULL, 
              help="Obtain squared CpGs from CpG2 linear model", metavar="character"),
    make_option(c("--squaredsubset_n_cpg2"), type="character", default=NULL, 
              help="Number of squared CpGs that were subsetted on in training (from CpG^2 lm models)", metavar="character"),
    make_option(c("--age2cpg2"), type="character", default=NULL, 
              help="Use subsets from from age^2 EWAS and from age ~ CpG^2 linear models", metavar="character"),
    make_option(c("--external"), type="character", default=NULL, 
              help="External data used in training option (T/F)", metavar="character"),
    make_option(c("--lasso"), type="character", default=NULL, 
              help="Lasso option in training (T/F)", metavar="character"),
    make_option(c("--loo"), type="character", default=NULL, 
              help="Whether the leave one out option was chosen in training (T/F)", metavar="character"),
    make_option(c("--lbc"), type="character", default=NULL, 
              help="Whether to include LBC in testing (T/F)", metavar="character"),
    make_option(c("--lbc_df"), type="character", default=NULL, 
              help="LBC target file", metavar="character"),
    make_option(c("--lbc_meth"), type="character", default=NULL, 
              help="LBC methylation file", metavar="character"),
    make_option(c("--random"), type="character", default=NULL, 
              help="Totally random option (T/F) - ignored assigned folds and randomize data when training", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# For testing
######################################################

#weights_dir <- "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/chronological_age/elasticnet_models/elnet_train/w1w3w4/random_noscalesample_noscalecpg_subset20K/"
#output_dir <- "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/chronological_age/elasticnet_models/elnet_test/w1w3w4/random_noscalesample_noscalecpg_subset20K/"
#meth_data <- readRDS("/Cluster_Filespace/Marioni_Group/Elena/data/geo_data/test_methylation.RDS")
#sample_data <- readRDS("/Cluster_Filespace/Marioni_Group/Elena/data/geo_data/test_samples.RDS")
#scale_cpg <- "F"
#scale_sample <- "F"
#sex_strat <- "F"
#ext <- "F"
#loo <- "F"
#lasso <- "F"
#rand <- "F"
#cpg_s <- "F"
#cpg_s_subset <- "F"
#name <- "random_noscalesample_noscalecpg_subset20K"
#lbc <- "T"
#lbc_df <- read.delim("/Cluster_Filespace/Marioni_Group/Elena/data/lbc_data/lbc_targets_cagetesting.tsv", row.names = 1, sep = " ")
#lbc_meth <- readRDS("/Cluster_Filespace/Marioni_Group/Elena/data/lbc_data/lbc_bvals_cagetesting.rds")
#w <- paste0(weights_dir, "elnet_coefficients_", name, squared_f, squaredsubset_f, random_f, lasso_f, ".tsv") 

# Log age option
######################################################

logage <- opt$logage
logage20 <- opt$logage20

# Scale options
######################################################

scale_cpg <- opt$cpg
scale_sample <- opt$sample


# Squared CpG option
######################################################

cpg_s <- opt$squared
cpg_s_subset <- opt$squaredsubset
cpg_s_cpg2 <- opt$squaredcpg2
cpg_s_age2cpg2 <- opt$age2cpg2
if ((cpg_s_subset == "T") & (cpg_s_age2cpg2 == "F")) {
    cpg_subset_n <- opt$squaredsubset_n
}
if ((cpg_s_subset == "T") & (cpg_s_age2cpg2 == "T")) {
    cpg_subset_n <- opt$squaredsubset_n
    cpg_subset_n_cpg2 <- opt$squaredsubset_n_cpg2

}


# Sex stratified option
######################################################

sex_strat <- opt$sex


# Was external data used in training?
######################################################

ext <- opt$external


# Cohort selected for testing?
######################################################

test_cohort <- opt$test_cohort
if (test_cohort != "all") {
    test_cohort <- unlist(str_split(test_cohort, ","))
}
print(test_cohort)

# Were folds randomized in training?
######################################################

rand <- opt$random


# Was lasso used in training?
######################################################

lasso <- opt$lasso


# Was the LOO framework used in training?
######################################################

loo <- opt$loo


# Use LBC in testing?
######################################################

lbc <- opt$lbc


# Stuff for output/input
######################################################

lasso_f <- ""
if (lasso == "T") {
    lasso_f <- "_lasso"
}

random_f <- ""
if (rand == "T") {
    random_f <- "_randomizedfolds"
}

squared_f <- ""
if (cpg_s == "T") {
    squared_f <- "_squared"
}

sex_f <- ""
if (sex_strat == "T") {
    sex_f <- "_sexstratified"
}

loo_f <- ""
if (loo == "T") {
    loo_f <- "_loo"
}

lbc_f <- ""
if (lbc == "T") {
    lbc_f <- "_lbc"
}

squaredsubset_f <- ""
if (cpg_s_subset == "T") {
    squaredsubset_f <- paste0("_squaredsubset",cpg_subset_n)
}

squared_cpg_f <- ""
if ((cpg_s_cpg2 == "T") & (cpg_s_age2cpg2 == "F")) {
    squaredsubset_f <- paste0("_squaredsubset",cpg_subset_n,"K")
    squared_cpg_f <- "_cpg2lm"
}

if ((cpg_s_cpg2 == "T") & (cpg_s_age2cpg2 == "T")) {
    squaredsubset_f <- paste0("_squaredsubset", cpg_subset_n, "_squaredsubsetcpg2lm", cpg_subset_n_cpg2, "K")
}

logage20_f <- ""
if (logage20 %in% c("Y", "O")) {
    logage20_f <- paste0("_logage20_", logage20)
}

test_cohort_f <- ""
if (test_cohort != "all") {
    test_cohort_f <- paste0("_test", opt$test_cohort)
}

name <- opt$name
weights_dir <- opt$weights_dir


# Import + prep data
######################################################

meth_data <- readRDS(opt$meth)
sample_data <- readRDS(opt$pheno)
output_dir <- opt$out

if (lbc == "T") {
    lbc_df <- read.delim(opt$lbc_df, row.names = 1, sep = " ")
    lbc_meth <- readRDS(opt$lbc_meth)
    common <- Reduce(intersect, list(colnames(meth_data), colnames(lbc_meth)))
    # Keep common
    lbc_meth <- lbc_meth[,common]
    meth_data <- meth_data[,common]
    # Bind data
    meth_data <- rbind(meth_data, lbc_meth)
    sample_data <- rbind(sample_data, lbc_df)
}

geos <- unique(sample_data$GEO)
cat("\nHave imported testing methylation + sample data.\n")

# Make sure pheno and sample order matches
meth_data <- meth_data[match(rownames(sample_data), rownames(meth_data)),]
cat("\nHave matched meth + sample rows.\n")

# Filter out if some GEO sets have already been used in training...
if (ext == "T") {
    if (loo == "F") { # Will deal with filtering for testing with LOO later
        sample_data <- sample_data[!((sample_data["GEO"] == "GSE40279") | (sample_data["GEO"] == "GSE42861")),]
        meth_data <- meth_data[!((sample_data["GEO"] == "GSE40279") | (sample_data["GEO"] == "GSE42861")),]
    }
}   

print(dim(sample_data))
# If cohort(s) has been selected for testing
if (test_cohort != "all") {
    sample_data <- sample_data[sample_data$GEO %in% test_cohort,]
    meth_data <- meth_data[sample_data$GEO %in% test_cohort,]
}
print(dim(sample_data))

# Scale per sample if needed
if (scale_sample == "T") {
    meth_data <- t(meth_data)
    meth_data <- scale(meth_data)
    meth_data <- t(meth_data)
    cat("\nHave scaled per sample as per user's choice.\n")
}

# Scale per CpG if needed
if (scale_cpg == "T") {
    for (set in unique(sample_data["GEO"])) { # Scale per GEO set
        if (sex_strat == "T") {
            samples_F <- rownames(sample_data[(sample_data["GEO"] == set) & (sample_data$Sex == "F"),])
            samples_M <- rownames(sample_data[(sample_data["GEO"] == set) & (sample_data$Sex == "M"),])
            sub_F <- scale(meth_data[samples_F,])
            sub_M <- scale(meth_data[samples_M,])
            rownames(sub_F) <- samples_F
            rownames(sub_M) <- samples_M
            meth_data[samples_F,] <- sub_F
            meth_data[samples_M,] <- sub_M
            rm(sub_F)
            rm(sub_M) 
        } else {
            samples <- rownames(sample_data[sample_data["GEO"] == set,])
            sub <- scale(meth_data[samples,])
            rownames(sub) <- samples
            meth_data[samples,] <- sub
            rm(sub)
        }
    }
    cat("\nHave scaled per CpG as per user's choice.\n")
}
    
# Impute NAs
meth_data <- na_mean(meth_data)

# Square methylation if needed
if ((cpg_s == "T") | (cpg_s_subset == "T")) {
    names <- c(colnames(meth_data),paste0(colnames(meth_data),'_2'))
    meth_data <- cbind(meth_data,meth_data^2)
    colnames(meth_data) <- names
    cat("\nSquared CpGs.\n")
}

# Import weights
if (loo == "T") {
    if (sex_strat == "F") {
        weights_list <- list()
        for (geo in geos) {
            w <- paste0(weights_dir, "elnet_coefficients_", name, "_loo_", geo, squared_f, squaredsubset_f, squared_cpg_f, random_f, lasso_f, ".tsv") 
            weights_list[[geo]] <- read.delim(w, row.names = 1)
        }
    } else {
        weights_list_F <- list()
        weights_list_M <- list()
        for (geo in geos) {
            w_F <- paste0(weights_dir, "elnet_coefficients_", name, "_loo_", geo, squared_f, squaredsubset_f, squared_cpg_f, random_f, lasso_f, "_F.tsv") 
            w_M <- paste0(weights_dir, "elnet_coefficients_", name, "_loo_", geo, squared_f, squaredsubset_f, squared_cpg_f, random_f, lasso_f, "_M.tsv") 
            weights_list_F[[geo]] <- read.delim(w_F, row.names = 1)
            weights_list_M[[geo]] <- read.delim(w_M, row.names = 1)
        }
    }
} else {
    if (sex_strat == "T") {
        w_F <- paste0(opt$weights_dir, "elnet_coefficients_", opt$name, squared_f, squaredsubset_f, squared_cpg_f, random_f, lasso_f, "_F.tsv") 
        w_M <- paste0(opt$weights_dir, "elnet_coefficients_", opt$name, squared_f, squaredsubset_f, squared_cpg_f, random_f, lasso_f, "_M.tsv") 
        weights_F <- read.delim(w_F, row.names = 1)
        weights_M <- read.delim(w_M, row.names = 1)
    } else {
        w <- paste0(opt$weights_dir, "elnet_coefficients_", opt$name, squared_f, squaredsubset_f, squared_cpg_f, random_f, lasso_f, ".tsv") 
        weights <- read.delim(w, row.names = 1)
    }
}
cat("\nHave imported weights for CpGs.\n")

# Only keep CpG data for selected features and prep x and y
if (loo == "T") {
    if (sex_strat == "T") {
        # Model parameters
        intercept_list_F <- list()
        intercept_list_M <- list()
        features_list_F <- list()
        features_list_M <- list()
        # Age (y) and meth (x)
        x_test_list_F <- list()
        x_test_list_M <- list()
        y_test_list_F <- list()
        y_test_list_M <- list()
        for (geo in geos) {
            # F
            intercept_list_F[[geo]] <- weights_list_F[[geo]]["Intercept", "Coefficient"] # Extract intercept
            features_list_F[[geo]] <- weights_list_F[[geo]][2:nrow(weights_list_F[[geo]]),"Coefficient",drop=FALSE] # Keep only CpGs
            features_list_F[[geo]] <- features_list_F[[geo]][rownames(features_list_F[[geo]]) %in% colnames(meth_data),,drop = FALSE] # Keep only features present in testing data
            
            # M
            intercept_list_M[[geo]] <- weights_list_M[[geo]]["Intercept", "Coefficient"] # Extract intercept
            features_list_M[[geo]] <- weights_list_M[[geo]][2:nrow(weights_list_M[[geo]]),"Coefficient",drop=FALSE] # Keep only CpGs
            features_list_M[[geo]] <- features_list_M[[geo]][rownames(features_list_M[[geo]]) %in% colnames(meth_data),,drop = FALSE] # Keep only features present in testing data
            
            # Filter methylation data just for females/males + GEO set
            geo_sample_data <- sample_data[(sample_data["GEO"] == geo),]
            geo_meth_data <- meth_data[(sample_data["GEO"] == geo),]
            fem <- rownames(geo_sample_data[geo_sample_data$Sex == "F",])
            mal <- rownames(geo_sample_data[geo_sample_data$Sex == "M",])
            meth_data_F <- geo_meth_data[rownames(geo_meth_data) %in% fem,]
            meth_data_M <- geo_meth_data[rownames(geo_meth_data) %in% mal,]
            
            # x
            x_test_list_F[[geo]] <- meth_data_F[,match(rownames(features_list_F[[geo]]), colnames(meth_data_F))]
            x_test_list_M[[geo]] <- meth_data_M[,match(rownames(features_list_M[[geo]]), colnames(meth_data_M))]
            x_test_list_F[[geo]] <- t(x_test_list_F[[geo]])
            x_test_list_M[[geo]] <- t(x_test_list_M[[geo]])

            # y
            y_test_list_F[[geo]] <- fem[,"Age",drop=FALSE]
            y_test_list_F[[geo]] <- y_test_list_F[[geo]][colnames(x_test_list_F[[geo]]),,drop=FALSE]
            y_test_list_M[[geo]] <- mal[,"Age",drop=FALSE]
            y_test_list_M[[geo]] <- y_test_list_M[[geo]][colnames(x_test_list_M[[geo]]),,drop=FALSE]
        }
    } else {
        # Model parameters
        intercept_list <- list()
        features_list <- list()
        # Age (y) and meth (x)
        x_test_list <- list()
        y_test_list <- list()
        for (geo in geos) {
            intercept_list[[geo]] <- weights_list[[geo]]["Intercept","Coefficient"] # Extract intercept
            features_list[[geo]] <- weights_list[[geo]][2:nrow(weights_list[[geo]]),"Coefficient",drop=FALSE] # Keep only CpGs
            features_list[[geo]] <- features_list[[geo]][rownames(features_list[[geo]]) %in% colnames(meth_data),,drop = FALSE] # Keep only features present in testing data
            # Meth data for geo
            geo_sample_data <- sample_data[(sample_data["GEO"] == geo),]
            geo_meth_data <- meth_data[(sample_data["GEO"] == geo),]
            # Keep just CpGs
            x_test_list[[geo]] <- geo_meth_data[,match(rownames(features_list[[geo]]), colnames(geo_meth_data))]
            # Match x and y
            x_test_list[[geo]] <- t(x_test_list[[geo]])
            y_test_list[[geo]] <- geo_sample_data[,"Age",drop=FALSE]
            y_test_list[[geo]] <- y_test_list[[geo]][colnames(x_test_list[[geo]]),,drop=FALSE]
        }
    }
} else {
    if (sex_strat == "T") {
        # F
        intercept_F <- weights_F["Intercept","Coefficient"] # Extract intercept
        features_F <- weights_F[2:nrow(weights_F),"Coefficient",drop=FALSE] # Keep only CpGs
        features_F <- features_F[rownames(features_F) %in% colnames(meth_data),,drop = FALSE] # Keep only features present in testing data
        
        # M
        intercept_M <- weights_M["Intercept","Coefficient"] # Extract intercept
        features_M <- weights_M[2:nrow(weights_M),"Coefficient",drop=FALSE] # Keep only CpGs
        features_M <- features_M[rownames(features_M) %in% colnames(meth_data),,drop = FALSE] # Keep only features present in testing data
        
        # Filter methylation data just for females/males
        fem <- rownames(sample_data[sample_data$Sex == "F",]) # 1544
        mal <- rownames(sample_data[sample_data$Sex == "M",]) # 1329
        meth_data_F <- meth_data[rownames(meth_data) %in% fem,] # 1445
        meth_data_M <- meth_data[rownames(meth_data) %in% mal,] # 1230
        
        # Keep just CpGs
        x_test_F <- meth_data_F[,match(rownames(features_F), colnames(meth_data_F))]
        x_test_M <- meth_data_M[,match(rownames(features_M), colnames(meth_data_M))]

        # Match x and y
        x_test_F <- t(x_test_F)
        y_test_F <- sample_data[,"Age",drop=FALSE]
        y_test_F <- y_test_F[rownames(y_test_F) %in% colnames(x_test_F),,drop=FALSE]
        x_test_M <- t(x_test_M)
        y_test_M <- sample_data[,"Age",drop=FALSE]
        y_test_M <- y_test_M[rownames(y_test_M) %in% colnames(x_test_M),,drop=FALSE]
    } else {
        intercept <- weights["Intercept","Coefficient"] # Extract intercept
        features <- weights[2:nrow(weights),"Coefficient",drop=FALSE] # Keep only CpGs
        features <- features[rownames(features) %in% colnames(meth_data),,drop = FALSE] # Keep only features present in testing data
        
        # Keep just CpGs
        x_test <- meth_data[,match(rownames(features), colnames(meth_data))]

        # Match x and y
        x_test <- t(x_test)
        y_test <- sample_data[,"Age",drop=FALSE]
        y_test <- y_test[colnames(x_test),,drop=FALSE]
    }    
}

cat("\nHave prepped + scaled per CpG as per user's choice.\n")


# Prediction
######################################################

if (loo == "T") {
    pred_df <- ""
    i <- 1
    if (sex_strat == "T") {
        for (geo in geos) {
            pred_geo_F <- x_test_list_F[[geo]] * features_list_F[[geo]][,"Coefficient"]
            pred_pp_geo_F <- colSums(pred_geo_F)
            pred_pp_geo_F <- pred_pp_geo_F + intercept_list_F[[geo]]
            pred_df_geo_F <- data.frame("Age_pred" = pred_pp_geo_F, "Age" = y_test_list_F[[geo]]["Age"], "GEO" = sample_data[rownames(y_test_list_F[[geo]]), "GEO"], "Sex" = sample_data[rownames(y_test_list_F[[geo]]), "Sex"])

            pred_geo_F <- x_test_list_F[geo] * features_list_F[[geo]][,"Coefficient"]
            pred_pp_geo_F <- colSums(pred_geo_F)
            pred_pp_geo_F <- pred_pp_geo_F + intercept_list_F[[geo]]
            pred_df_geo_F <- data.frame("Age_pred" = pred_pp_geo_F, "Age" = y_test_list_F[[geo]]["Age"], "GEO" = sample_data[rownames(y_test_list_F[[geo]]), "GEO"], "Sex" = sample_data[rownames(y_test_list_F[[geo]]), "Sex"])

            pred_df_geo <- rbind(pred_df_geo_F, pred_df_geo_M)
            if (i == 1) {
                pred_df <- pred_df_geo
            } else {
                pred_df <- rbind(pred_df, pred_df_geo)
            }
            i <- i + 1
        }
    } else {
        for (geo in geos) {
            pred_geo <- x_test_list[[geo]] * features_list[[geo]][,"Coefficient"]
            pred_pp_geo <- colSums(pred_geo)
            pred_pp_geo <- pred_pp_geo + intercept_list[[geo]]
            pred_df_geo <- data.frame("Age_pred" = pred_pp_geo, "Age" = y_test_list[[geo]]["Age"], "GEO" = sample_data[rownames(y_test_list[[geo]]), "GEO"], "Sex" = sample_data[rownames(y_test_list[[geo]]), "Sex"])
            
            if (i == 1) {
                pred_df <- pred_df_geo
            } else {
                pred_df <- rbind(pred_df, pred_df_geo)
            }
            i <- i + 1        
        }
    }
} else {
    if (sex_strat == "T") {
        pred_F <- x_test_F * features_F[,"Coefficient"]
        pred_pp_F <- colSums(pred_F)
        pred_pp_F <- pred_pp_F + intercept_F
        pred_df_F <- data.frame("Age_pred" = pred_pp_F, "Age" = y_test_F["Age"], "GEO" = sample_data[rownames(y_test_F), "GEO"], "Sex" = rep("F", length(pred_pp_F)))

        pred_M <- x_test_M * features_M[,"Coefficient"]
        pred_pp_M <- colSums(pred_M)
        pred_pp_M <- pred_pp_M + intercept_M
        pred_df_M <- data.frame("Age_pred" = pred_pp_M, "Age" = y_test_M["Age"], "GEO" = sample_data[rownames(y_test_M), "GEO"], "Sex" = rep("M", length(pred_pp_M)))

        pred_df <- rbind(pred_df_F, pred_df_M)
    } else {
        pred <- x_test * features[,"Coefficient"]
        pred_pp <- colSums(pred)
        pred_pp <- pred_pp + intercept
        pred_df <- data.frame("Age_pred" = pred_pp, "Age" = y_test["Age"], "GEO" = sample_data[rownames(y_test), "GEO"], "Sex" = sample_data[rownames(y_test), "Sex"])
    }    
}
cat("\nObtained predictions! Exporting...\n")

pred_df <- pred_df[!is.na(pred_df$Age),]
age <- 20
cat("\nAll people:\n")
cat(dim(pred_df))
cat("\n")
if (logage == "T") {
    if (logage20 == "F") {
        pred_df$Age_pred <- exp(pred_df$Age_pred)
    }
    if (logage20 == "Y") {
        old_file <- paste0(output_dir, "predictions", sex_f, squared_f, squaredsubset_f, squared_cpg_f, random_f, lasso_f, loo_f, lbc_f, test_cohort_f, "_logage20_O", ".tsv")
        old_pred <- read.delim(old_file, row.names = 1)
        pred_df$Age_pred <- exp(pred_df$Age_pred)
        pred_df <- pred_df[!(rownames(pred_df) %in% rownames(old_pred)),]
        cat("\nYoung people:\n")
        cat(dim(pred_df))
    }
    if (logage20 == "O") {
        pred_df <- pred_df[pred_df$Age_pred > age,]
        cat("\nOld people:\n")
        cat(dim(pred_df))
    }
}

# Output
######################################################

# Output table
write.table(data.frame("Sample" = rownames(pred_df), pred_df), file = paste0(output_dir, "predictions", sex_f, squared_f, squaredsubset_f, squared_cpg_f, random_f, lasso_f, loo_f, lbc_f, test_cohort_f, logage20_f, ".tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
