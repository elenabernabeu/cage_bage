#!/usr/bin/env Rscript 
# -*- coding: utf-8 -*-
# By Elena

library("optparse")
library("data.table")


option_list = list(
    make_option(c("-p", "--pcaprep"), type="character", default=NULL, 
              help="PCA option", metavar="character"),
    make_option(c("-r", "--cpgprep"), type="character", default=NULL, 
              help="CpG trajectory across cohorts option", metavar="character"),
    make_option(c("-m", "--meth"), type="character", default=NULL, 
              help="Methylation file", metavar="character"),
    make_option(c("-t", "--target"), type="character", default=NULL, 
              help="Target file", metavar="character"),
    make_option(c("-c", "--cpg"), type="character", default=NULL, # ELOVL2
              help="CpG for which to obtain trajectory", metavar="character"),
    make_option(c("-n", "--name"), type="character", default="cg16867657",
              help="Name for output file", metavar="character"),
    make_option(c("-o", "--outputdir"), type="character", default=NULL,
              help="Path to output directory", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# For testing
######################################################

#cpg <- "cg16867657"
#name <- "noscale_gs20knorm_external"
#meth_file <- "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/chronological_age/data_prep/w1w3w4/methylation_training_noscale_gs20knorm_external.rds"
#target_file <- "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/chronological_age/elasticnet_models/cv_folds/w1w3w4/random/gs_basic_folds_random_withexternal.tsv"
#output_dir <- "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/chronological_age/data_prep/w1w3w4/"
#cpg_prep(cpg = cpg, meth_file = meth_file, target_file = target_file, name = "noscale_gs20knorm_external", output_dir = output_dir)


# Prepare ELOVL2 data across GS, LBC, and GEO cohorts
# cg16867657
######################################################


cpg_prep = function(cpg, meth_file, target_file, name, output_dir, scale = FALSE) {
    # Methylation
    meth <- readRDS(meth_file)
    meth <- meth[,cpg,drop=FALSE]
    
    # Target
    if (grepl("external", name, fixed = TRUE)) {
        target <- read.delim(target_file, header = TRUE, row.names = 1)
    } else {
        target <- read.delim(target_file, header = TRUE, row.names = 2)
    }
    target <- target[match(rownames(meth), rownames(target)),]
    
    # Fuse
    if (grepl("external", name, fixed = TRUE)) {
        df <- data.frame("ID" = rownames(meth), "CpG" = meth[,cpg], "Age" = target$age, "Cohort" = target$cohort, "Sex" = target$sex)
    } else {
        df <- data.frame("ID" = rownames(meth), "CpG" = meth[,cpg], "Age" = target$age, "Cohort" = target$Set, "Sex" = target$sex)
    }
    
    # Export
    write.table(df, paste0(output_dir, cpg, "_trajectory_", name, ".tsv"), row.names = FALSE, quote = FALSE, sep = "\t")
}

pca_prep = function(meth_file, name, output_dir) {
    # Methylation
    meth <- readRDS(meth_file)
    # Export to txt file so can be read by python
    fwrite(data.frame("ID" = rownames(meth), meth), paste0(output_dir, "methbetavals_", name, ".csv"), quote = FALSE, sep = ",", row.names = FALSE)
}

if (opt$pcaprep == TRUE) {
    pca_prep(meth_file = opt$meth, name = opt$name, output_dir = opt$outputdir)
}

if (opt$cpgprep == TRUE) {
    cpg_prep(cpg = opt$cpg, meth_file = opt$meth, target_file = opt$target, name = opt$name, output_dir = opt$outputdir)
}