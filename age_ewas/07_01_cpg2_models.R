#!/usr/bin/env Rscript 
# -*- coding: utf-8 -*-
# By Elena

library("optparse")

option_list = list(
    make_option(c("-c", "--chrom"), type="character", default=NULL, 
              help="Chromosome to iterate through", metavar="character"),
    make_option(c("-l", "--linear"), type="character", default="F", 
              help="Include only linear form of CpG in model", metavar="character")     
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
chrom <- opt$chrom
linear <- opt$linear

# M-vals (per chromosome)
df <- readRDS(paste0("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/Chromosomes/GS20k_chr", chrom, "_mvals.rds"))
df <- t(df)
gc()

# CpGs and samples to include in model
#cpgs <- read.delim("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_cpgs_common.txt", header = FALSE)$V1
cpgs <- read.delim("/Cluster_Filespace/Marioni_Group/Elena/gs_osca/data/cpgs_all_831733.txt", header = FALSE)$V1
cpgs <- cpgs[cpgs %in% colnames(df)]
target <- read.delim("/Cluster_Filespace/Marioni_Group/Elena/data/gs_data/gs20ktargets_fixedage.tsv", header = TRUE, row.names = 1)

# Filter to just CpGs and samples of interest
df <- df[,cpgs]
df <- df[rownames(target),]

# Scale CpG and CpG^2
df_scaled <- scale(df)
df_squared_scaled <- scale(df^2) 
print(df_scaled[1:5, 1:5])
print(df_squared_scaled[1:5, 1:5])
rm(df)
gc()

# Import covariates 
covars <- read.delim("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_allcovars.tsv", header = TRUE, row.names = 1)
covars <- covars[rownames(target),]

if (linear == "T") {
    # Make dataframe for results
    results_chrom <- data.frame("cpg_beta" = rep(NA, length(colnames(df_scaled))), "cpg_se" = rep(NA, length(colnames(df_scaled))), "cpg_p" = rep(NA, length(colnames(df_scaled))))
    rownames(results_chrom) <- colnames(df_scaled)

    # Iterate over all CpGs and run models
    i <- 1
    for (cpg in colnames(df_scaled)) {
        cat(paste0("CpG: ", cpg, "(", i, " of ", length(colnames(df_scaled)), ")...\n"))
        df_cpg <- data.frame("cpg" = df_scaled[,cpg], covars)
        cpg_model <- lm(age ~ cpg + as.factor(sex) + as.factor(batch) + as.factor(smoking_status) + pack_years +
                        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + 
                        CD8T + CD4T + NK + Bcell + Gran, data = df_cpg)
        results_chrom[cpg, "cpg_beta"] = summary(cpg_model)$coefficients["cpg", 1]
        results_chrom[cpg, "cpg_se"] = summary(cpg_model)$coefficients["cpg", 2]
        results_chrom[cpg, "cpg_p"] = summary(cpg_model)$coefficients["cpg", 4]
        i <- i + 1
    }
    write.table(data.frame("cpg" = rownames(results_chrom), results_chrom), paste0("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/cpg2_lm/w1w3w4/cpg_lm_chrom", chrom, "_scaled.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

} else {
    # Make dataframe for results
    results_chrom <- data.frame("cpg_beta" = rep(NA, length(colnames(df_scaled))), "cpg2_beta" = rep(NA, length(colnames(df_scaled))), "cpg_se" = rep(NA, length(colnames(df_scaled))), "cpg2_se" = rep(NA, length(colnames(df_scaled))), "cpg_p" = rep(NA, length(colnames(df_scaled))), "cpg2_p" = rep(NA, length(colnames(df_scaled))))
    rownames(results_chrom) <- colnames(df_scaled)

    # Iterate over all CpGs and run models
    i <- 1
    for (cpg in colnames(df_scaled)) {
        cat(paste0("CpG: ", cpg, "(", i, " of ", length(colnames(df_scaled)), ")...\n"))
        df_cpg <- data.frame("cpg" = df_scaled[,cpg], "cpg2" = (df_squared_scaled[,cpg]), covars)
        cpg_model <- lm(age ~ cpg + cpg2 + as.factor(sex) + as.factor(batch) + as.factor(smoking_status) + pack_years +
                        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + 
                        CD8T + CD4T + NK + Bcell + Gran, data = df_cpg)
        results_chrom[cpg, "cpg_beta"] = summary(cpg_model)$coefficients["cpg", 1]
        results_chrom[cpg, "cpg2_beta"] = summary(cpg_model)$coefficients["cpg2", 1]
        results_chrom[cpg, "cpg_se"] = summary(cpg_model)$coefficients["cpg", 2]
        results_chrom[cpg, "cpg2_se"] = summary(cpg_model)$coefficients["cpg2", 2]
        results_chrom[cpg, "cpg_p"] = summary(cpg_model)$coefficients["cpg", 4]
        results_chrom[cpg, "cpg2_p"] = summary(cpg_model)$coefficients["cpg2", 4]
        i <- i + 1
    }
    write.table(data.frame("cpg" = rownames(results_chrom), results_chrom), paste0("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/cpg2_lm/w1w3w4/cpg2_lm_chrom", chrom, "_scaled.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
