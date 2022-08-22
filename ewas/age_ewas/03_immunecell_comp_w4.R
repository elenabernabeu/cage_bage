#!/usr/bin/env Rscript 
# -*- coding: utf-8 -*-
# By Elena

library("minfi")
library("stringr")
library("FlowSorted.Blood.EPIC")
library("FlowSorted.Blood.450k")

# Target info
target <- read.delim("/Cluster_Filespace/Marioni_Group/Elena/data/gs_data/gs20ktargets_fixedage.tsv", header = TRUE, row.names = 1)
w4_target <- target[target["Set"] == "W4",]
w4_target[,c("Slide", "Array")] <- str_split_fixed(rownames(w4_target), "_", 2)

# Idat location
idat <- "/Cluster_Filespace/Marioni_Group/GS/GS_methylation/wave4/idats/wave4-idats/"

# Targets dataframe for minfi
targets <- data.frame("Sample_Name" = w4_target$Sample_Name, "Slide" = w4_target$Slide, "Array" = w4_target$Array, "age" = w4_target$age, "sex" = w4_target$sex, "Basename" = paste0(idat, rownames(w4_target)))

# Separate into chunks
chunks <- 10
por <- ceiling(length(rownames(targets))/chunks)
chunk_list <- list()

for (i in 1:chunks) {
    if (i == 1) {
        t_chunk <- targets[1:(por-1),]
    } else if (i == chunks) {
        t_chunk <- targets[(por*(i-1)):length(rownames(targets)),]
    } else {
        t_chunk <- targets[(por*(i-1)):((por*i)-1),]
    }
    chunk_list <- append(chunk_list, list(t_chunk))
}

# Obtain cell compositions for those
j <- 1
count_list <- list()
for (chunk in chunk_list) {
    cat(paste0("\nWorking on chunk ", j, "...\n"))
    j <- j + 1

    # Create RGset
    RGset <- read.metharray.exp(targets = chunk, force = TRUE)
    RGset@annotation=c(array='IlluminaHumanMethylationEPIC', annotation='ilm10b2.hg19')
    cat("- Imported RGsets.\n")
    gc()

    # Get cell counts
    cell_counts <- estimateCellCounts(RGset, compositeCellType = "Blood",
                    processMethod = "auto", probeSelect = "auto",
                    cellTypes = c("CD8T","CD4T","NK","Bcell","Mono","Gran"),
                    returnAll = FALSE, meanPlot = FALSE, verbose = TRUE)
    cat("- Calculated cell counts.\n")
    gc()

    # Append
    count_list <- append(count_list, list(cell_counts))
}

# Export
composition <- do.call("rbind", count_list)
write.table(data.frame("ID" = rownames(composition), composition), "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/wave4_cellcomp.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Compare with DNAm predictions
cells <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
minfi <- read.table("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/wave4_cellcomp.tsv", header = TRUE, row.names = 1)
dnam <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS_w4_DNAmAge.csv", row.names = 1)
dnam <- dnam[, cells]
names(dnam) <- paste0("DNAm_", names(dnam))
comp <- merge(minfi, dnam, by="row.names")
for (cell in cells) {
    r <- cor(comp[cell], comp[paste0("DNAm_", cell)])
    cat(paste0("Correlation for cell ", cell, ": ", r, "\n"))
}