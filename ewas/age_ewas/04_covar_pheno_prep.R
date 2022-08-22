#!/usr/bin/env Rscript 
# -*- coding: utf-8 -*-
# By Elena

# Target file
######################################################

target <- read.delim("/Cluster_Filespace/Marioni_Group/Elena/data/gs_data/gs20ktargets_fixedage.tsv", header = TRUE, row.names = 1)

# PCs
######################################################

#pcs <- read.delim("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_PCA.tsv", header = TRUE, row.names = 1)
pcs <- read.delim("/Cluster_Filespace/Marioni_Group/Elena/gs_osca/data/norm20k_18413_752722_PCA.tsv", header = TRUE, row.names = 1)
pcs <- pcs[rownames(target),]
pcs <- pcs[,1:20]

# Immune cells
######################################################

cells  <- c("CD8T", "CD4T", "NK", "Bcell", "Gran")
immune_w1 <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/stradl-samples-5087.rds")
immune_w1["Sample_Sentrix_ID"] <- paste(immune_w1$Sentrix_ID, immune_w1$Sentrix_Position, sep = "_")
rownames(immune_w1) <- immune_w1$Sample_Sentrix_ID
immune_w1 <- immune_w1[, cells]

immune_w3 <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/wave3-final/samplesheet.final.csv")
rownames(immune_w3) <- immune_w3$Sample_Sentrix_ID
immune_w3 <- immune_w3[, cells]

immune_w4 <- read.table("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/wave4_cellcomp.tsv", row.names = 1, header = TRUE)
immune_w4 <- immune_w4[, cells]

cell_comp <- rbind(immune_w1, immune_w3, immune_w4)
cell_comp <- cell_comp[rownames(target),]

# Smoking status + years
######################################################

smoking <- read.csv("/Cluster_Filespace/Marioni_Group/Elena/data/gs_data/smoking.csv", header = TRUE, row.names = 1)
smoking <- smoking[target$Sample_Name,]
smoking[which(is.na(smoking$ever_smoke)),"ever_smoke"] <- 5
for(i in c(1,2,3)) {
    smoking[which(is.na(smoking$pack_years) & smoking$ever_smoke == i), "pack_years"] <- median(smoking[which(!is.na(smoking$pack_years) &  smoking$ever_smoke == i), "pack_years"])
}
smoking[which(smoking$ever_smoke == 5), "pack_years"] <- 0        

# Prep .opi files for OSCA with old info
######################################################

anno <- readRDS("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds")
osca_dat <- read.delim("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_mvals_squared_fixedage_commonsnps.opi", header = FALSE)
rownames(osca_dat) <- osca_dat$V2
opi <- anno[rownames(osca_dat),c("chr","Name", "pos","UCSC_RefGene_Name", "strand")]
opi$chr <- gsub("chr", "", opi$chr)
opi$chr <- gsub("X", "23", opi$chr)
opi$chr <- gsub("Y", "24", opi$chr)
opi$chr <- as.numeric(opi$chr)
opi$UCSC_RefGene_Name  <- as.factor(opi$UCSC_RefGene_Name )
opi$strand <- as.factor(opi$strand)
opi[which(opi$UCSC_RefGene_Name ==""), "UCSC_RefGene_Name"] <- NA
write.table(opi, file="/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_mvals_squared_fixedage_commonsnps.opi", col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')

# Prep qualitative covar files for OSCA
######################################################

covar <- data.frame("FID" = rownames(target), "IID" = rownames(target), "sex" = target$sex, "batch" = target$Batch, "wave" = target$Set, "smoking_status" = smoking$ever_smoke)
write.table(covar, "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_covar.txt", sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)

# With only batch
covar_onlybatch <- covar[c("FID", "IID", "sex", "batch", "smoking_status")]
write.table(covar_onlybatch, "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_covar_nowave.txt", col.names = FALSE, sep = " ", row.names = FALSE, quote = FALSE)

# Prepare quantitative covar files for OSCA
######################################################

target$agesc <- scale(target$age, scale = FALSE)
target$agesc2 <- (target$agesc)^2

qcovar <- data.frame("FID" = rownames(target), "IID" = rownames(target), "pack_years" = smoking$pack_years, pcs, cell_comp)
qcovar_age <- data.frame("FID" = rownames(target), "IID" = rownames(target), "age" = target$agesc, "pack_years" = smoking$pack_years, pcs, cell_comp)
write.table(qcovar, "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_qcovar.txt", sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(qcovar_age, "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_qcovar_withage.txt", sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Prep pheno files for OSCA
######################################################

age <- data.frame("FID" = rownames(target), "IID" = rownames(target), "age" = target$agesc)
age2 <- data.frame("FID" = rownames(target), "IID" = rownames(target), "age2" = target$agesc2)
write.table(age, "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_age.txt", sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(age2, "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_age2.txt", sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Prep all covars and age in single file to have
######################################################

all_df <- data.frame("ID" = rownames(target), "age" = target$agesc, "age2" = target$agesc2, "sex" = target$sex, "batch" = target$Batch, "wave" = target$Set, "smoking_status" = smoking$ever_smoke, "pack_years" = smoking$pack_years, pcs, cell_comp)
write.table(all_df, "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_allcovars.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
