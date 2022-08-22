#!/usr/bin/env Rscript 
# -*- coding: utf-8 -*-
# By Elena

library(data.table)

# Import
df <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/mvals.rds")
target <- read.delim("/Cluster_Filespace/Marioni_Group/Elena/data/gs_data/gs20ktargets_fixedage.tsv", header = TRUE, row.names = 1)
cpgs <- read.delim("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/chronological_age/data_prep/w1w3w4/methylation_training_noscale_commonsnps.txt", header = FALSE)$V1
cpgs <- cpgs[which(cpgs %in% rownames(df))]

# Filter
df <- df[,rownames(target)] # Keep only 18413 individuals
df <- df[cpgs,]
df <- t(df)

# Export
fwrite(data.frame("FID" = rownames(df), "IID" = rownames(df), df), "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_mvals_fixedage_commonsnps.txt", quote = FALSE, row.names = FALSE, sep = " ")
write.csv(rownames(df), "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_samples.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

# Squared CpGs
df <- df^2
fwrite(data.frame("FID" = rownames(df), "IID" = rownames(df), df), "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_mvals_squared_fixedage_commonsnps.txt", quote = FALSE, row.names = FALSE, sep = " ")

# Then in shell
# nohup osca_Linux  --efile /Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_mvals_fixedage_commonsnps.txt --methylation-m --make-bod --out /Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_mvals_fixedage_commonsnps &
# nohup osca_Linux  --efile /Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_mvals_squared_fixedage_commonsnps.txt --methylation-m --make-bod --out /Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_mvals_squared_fixedage_commonsnps &