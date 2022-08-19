#!/usr/bin/env Rscript 
# -*- coding: utf-8 -*-
# By Elena

library("optparse")
library("glue")
library("imputeTS")
library("dplyr")
library("clusterSim")

# Get arguments
######################################################

option_list = list(
    make_option(c("-x", "--wave1"), type="character", default=NULL, 
              help="Methylation wave 1 file location", metavar="character"),
    make_option(c("-y", "--wave3"), type="character", default=NULL, 
              help="Methylation wave 3 file location", metavar="character"),
    make_option(c("-z", "--wave4"), type="character", default=NULL, 
              help="Methylation wave 4 file location", metavar="character"),
    make_option(c("-w", "--gs20k"), type="character", default=NULL, 
              help="Methylation GS20k file location", metavar="character"),
    make_option(c("-l", "--lbc"), type="character", default=NULL, 
              help="Methylation for LBC file location", metavar="character"),
    make_option(c("-g", "--geo"), type="character", default=NULL, 
              help="Methylation for GEO file location", metavar="character"),
    make_option(c("-c", "--cpgs"), type="character", default=NULL, 
              help="File with CpGs to limit to", metavar="character"),
    make_option(c("-p", "--pheno"), type="character", default=NULL, 
              help="Pheno file with fold info location", metavar="character"),
    make_option(c("-e", "--epic"), type="character", default=NULL, 
              help="EPIC file location", metavar="character"),
    make_option(c("-n", "--name"), type="character", default=NULL, 
              help="Run name", metavar="character"),
    make_option(c("-o", "--out"), type="character", default=NULL, 
              help="Output directory", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



# To test in terminal
######################################################

#w1 <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/norm_mvals_5087.rds")
#w3 <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/wave3_mvals.rds")
#w4 <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/wave4/w4-mvals.rds")
#gs20k <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/mvals.rds")
#lbc <- readRDS("/Cluster_Filespace/Marioni_Group/Elena/data/lbc_data/lbc_mvals_3489.rds")
#geo <- t(readRDS("/Cluster_Filespace/Marioni_Group/Elena/data/geo_data/test_methylation.RDS"))
#epic <- readRDS("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds")
#pheno <- read.delim("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/chronological_age/elasticnet_models/cv_folds/w1w3w4/random/gs_basic_folds_random.tsv", header = TRUE, row.names = 2)
#o_name_rds <- "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/chronological_age/data_prep/w1w3w4/methylation_training_noscale_gs20knorm.rds"



# Some functions
######################################################

m2beta <- function(m) { 
  beta <- 2^m/(2^m + 1)
  return(beta)
}

beta2m <- function(beta) { 
  m <- log2(beta/(1-beta))
  return(m)
}


# Read in the DNAm data 
######################################################

if (is.null(opt$gs20k)) {
    # Import data
    w1 <- readRDS(opt$wave1)
    w3 <- readRDS(opt$wave3)

    # Run with all waves?
    if (!is.null(opt$wave4)) {
    w4 <- readRDS(opt$wave4)
    }
} else {
    gs20k <- readRDS(opt$gs20k)
}

# External data in training?
if (!is.null(opt$lbc)) {
    lbc <- readRDS(opt$lbc)
    geo <- readRDS(opt$geo)
    geo <- t(geo)
    geo <- beta2m(geo)
}

# Fuse into list
if (!is.null(opt$wave4)) {
    if (!is.null(opt$lbc)) {
        train <- list("w1" = w1, "w3" = w3, "w4" = w4, "lbc" = lbc, "geo" = geo)
        rm(w1, w3, w4, lbc, geo) # Erase unnecessary
    } else {
        train <- list("w1" = w1, "w3" = w3, "w4" = w4)
        rm(w1, w3, w4) # Erase unnecessary
    }
} else if (!is.null(opt$gs20k)) {
    if (!is.null(opt$lbc)) {
        train <- list("gs20k" = gs20k, "lbc" = lbc, "geo" = geo)
        rm(gs20k, lbc, geo) # Erase unnecessary
    } else {
        train <- list("gs20k" = gs20k)
        rm(gs20k) # Erase unnecessary
    }    
} else {
    if (!is.null(opt$lbc)) {
        train <- list("w1" = w1, "w3" = w3, "lbc" = lbc, "geo" = geo)
        rm(w1, w3, lbc, geo) # Erase unnecessary
    } else {
        train <- list("w1" = w1, "w3" = w3)
        rm(w1, w3) # Erase unnecessary
    }
}

gc()
cat("\nImported methylation data.\n")



# Read in pheno + fold data
######################################################

if (!is.null(opt$lbc)){
    pheno <- read.delim(opt$pheno, header = TRUE, row.names = 1)
} else {
    pheno <- read.delim(opt$pheno, header = TRUE, row.names = 2)
}

# Remove NAs in pheno data for age
pheno <- pheno[!(is.na(pheno$age)),]
cat("\nImported age and fold data.\n")



# Clean up methylation data a lil'
######################################################

# Get common probes and keep common
if (is.null(opt$wave4)) {
    if (is.null(opt$gs20k)) {
        if (is.null(opt$lbc)) {
            common <- Reduce(intersect, list(rownames(train[["w1"]]), rownames(train[["w3"]])))
        } else {
            common <- Reduce(intersect, list(rownames(train[["w1"]]), rownames(train[["w3"]]), rownames(train[["lbc"]]), rownames(train[["geo"]])))
        }    
    } else {
        if (is.null(opt$lbc)) {
            common <- Reduce(intersect, list(rownames(train[["gs20k"]])))
        } else {
            common <- Reduce(intersect, list(rownames(train[["gs20k"]]), rownames(train[["lbc"]]), rownames(train[["geo"]])))
        }         
    }
} else {
    if (is.null(opt$lbc)) {
        common <- Reduce(intersect, list(rownames(train[["w1"]]), rownames(train[["w3"]]), rownames(train[["w4"]])))
    } else{
        common <- Reduce(intersect, list(rownames(train[["w1"]]), rownames(train[["w3"]]), rownames(train[["w4"]]), rownames(train[["lbc"]]), rownames(train[["geo"]])))
    }  
}

# Get EPIC array file and subset to probes common to 450k and EPIC array
epic <- readRDS(opt$epic)
common_array <- rownames(epic[which(epic$Methyl450_Loci == "TRUE"),])
rm(epic)

# Filter (also making sure sample in pheno file)
for (wave in names(train)) {
    train[[wave]] <- train[[wave]][which(rownames(train[[wave]]) %in% common),]
    train[[wave]] <- train[[wave]][which(rownames(train[[wave]]) %in% common_array),]
    train[[wave]] <- train[[wave]][,which(colnames(train[[wave]]) %in% rownames(pheno))]
}
cat("\nKept only samples with pheno (age) data, as well as probes in EPIC array and those common to all waves.\n")

# Print dimensions to make sure everything is alright
cat("\n#### Dimensions:\n")
for (wave in names(train)) {
    cat(paste(wave, "\n", sep = ""))
    cat(paste("\t", dim(train[[wave]])))
    cat("\n")
}
# w1:   386296  5087
# w3:   386296  4450
# lbc:  386296  3489
# geo:  386296  1345

cat("\nRAM clean up...\n\n")
gc()


# Obtain t() of each
for (wave in names(train)) {
    train[[wave]] <- t(train[[wave]])
    cat("\nRAM clean up...\n\n")
    gc()
}

cols <- colnames(train[[names(train)[1]]])

# Convert to beta values, impute NAs, match order
######################################################

for (wave in names(train)) {
    train[[wave]] <- m2beta(train[[wave]])
    train[[wave]] <- na_mean(train[[wave]])
    train[[wave]] <- train[[wave]][,cols]
    cat("\nRAM clean up...\n\n")
    gc()
}


# Fuse
######################################################

# Fuse all samples + methylation data into a single table
if (is.null(opt$wave4)) {
    if (is.null(opt$gs20k)) {
        if (is.null(opt$lbc)) {
            train_df <- rbind(train[["w1"]], train[["w3"]])
        } else {
            train_df <- rbind(train[["w1"]], train[["w3"]], train[["lbc"]], train[["geo"]])
        }
    } else {
        if (is.null(opt$lbc)) {
            train_df <- train[["gs20k"]]
        } else {
            train_df <- rbind(train[["gs20k"]], train[["lbc"]], train[["geo"]])
        }
    }
} else {
    if (is.null(opt$lbc)) {
        train_df <- rbind(train[["w1"]], train[["w3"]], train[["w4"]])
    } else{
        train_df <- rbind(train[["w1"]], train[["w3"]], train[["w4"]], train[["lbc"]], train[["geo"]])
    }  
}

cat("\n#### Dimensions after fusing waves:\n")
cat(paste("\t", dim(train_df), "\n"))
# 383320 
# 24940 

names(train_df) <- sub("^X", "", names(train_df))
rm(train)
cat("\nRAM clean up...\n\n")
gc()

train_df <- train_df[match(rownames(pheno), rownames(train_df)),] # Match order of methylation table (x) and phenotype table (y) 
cat("\nRAM clean up...\n\n")
gc()


# QC CpGs
######################################################

cpgs <- read.delim(opt$cpgs, header = FALSE)$V1
train_df <- train_df[,which(colnames(train_df) %in% cpgs)]
print(dim(train_df))


# Export data
######################################################

cat("\nExporting prepped data for models...\n")
o_name_rds <- paste0(opt$out, "methylation_training_", opt$name, ".rds")
saveRDS(train_df, o_name_rds, compress = FALSE)


# If external, export just external data without GS
######################################################

if (!is.null(opt$lbc)) {
    pheno_noGS <- pheno[!(pheno$cohort %in% c("W1", "W3", "W4")),]
    train_df <- train_df[rownames(train_df) %in% rownames(pheno_noGS),]
    
    cat("\nExporting prepped data for models (just external)...\n")
    o_name_rds <- paste0(opt$out, "methylation_training_", opt$name, "_noGS.rds")
    saveRDS(train_df, o_name_rds, compress = FALSE)
}
