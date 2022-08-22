#!/usr/bin/env Rscript 
# -*- coding: utf-8 -*-
# By Elena

library("optparse")
library("survival")
library("coxme")

option_list = list(
    make_option(c("--chrom"), type="character", default=NULL, 
              help="Chromosome to iterate through", metavar="character"),    
    make_option(c("--meth"), type="character", default=NULL, 
              help="Location of methylation file", metavar="character"),
    make_option(c("--cpgs"), type="character", default=NULL, 
              help="File with cpgs to filter to", metavar="character"),  
    make_option(c("--covars"), type="character", default=NULL, 
              help="File with other covar info", metavar="character"),
    make_option(c("--death"), type="character", default=NULL, 
              help="File with survival info", metavar="character"),
    make_option(c("--kinship"), type="character", default=NULL, 
              help="File with kinship matrix", metavar="character"),
    make_option(c("--output_dir"), type="character", default=NULL, 
              help="Directory to output Martingale residuals", metavar="character")  
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

extract_coxme_table <- function (mod){
    beta <- mod$coefficients
    nvar <- length(beta)
    nfrail <- nrow(mod$var) - nvar
    se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
    z <- round(beta/se, 2)
    p <- signif(1 - pchisq((beta/se)^2, 1), 2)
    confint_95_low <- exp(beta) - 1.96*se
    confint_95_high <- exp(beta) + 1.96*se
    table=data.frame(cbind(beta,exp(beta),se,z,p,confint_95_low,confint_95_high))
    return(table)
}


chrom <- opt$chrom

# M-vals (per chromosome)
# df <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/Chromosomes/GS20k_chr22_mvals.rds") # 18216 20210
df <- readRDS(opt$meth)
df <- t(df)
gc()

# Import covariates 
# covars <- read.delim("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_allcovars.tsv", header = TRUE, row.names = 1)
covars <- read.delim(opt$covars, header = TRUE, row.names = 1)

# Death
# death <- read.delim("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/biological_age/data_prep/gs_w1w3w4_bage_variables_scaledgrim.tsv", sep = "\t", row.names = 1)
death <- read.csv(opt$death, sep = "\t", row.names = 1)
death <- death[rownames(covars),]

# CpGs
# cpgs <- read.delim("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_cpgs_common.txt", header = FALSE)$V1
cpgs <- read.delim(opt$cpgs, header = FALSE)$V1
cpgs <- cpgs[cpgs %in% colnames(df)]

# Kinship matrix
# kinship_matrix <- readRDS("/Cluster_Filespace/Marioni_Group/Elena/data/gs_data/gs_kinship_20k.RDS")
kinship_matrix <- readRDS(opt$kinship)

# Filter to just CpGs and samples of interest
df <- df[,cpgs]
df <- df[rownames(covars),]

# Make dataframe for results
results_chrom <- data.frame(matrix(ncol=7,nrow=ncol(df), dimnames=list(colnames(df), c("logHR", "HR", "SE", "Z", "p", "HR_CI95_Low", "HR_CI95_High"))))

# Iterate over all CpGs and run models
i <- 1
for (cpg in colnames(df)) {
    cat(paste0("CpG: ", cpg, "(", i, " of ", length(colnames(df)), ")...\n"))
    df_cpg <- data.frame("cpg" = df[,cpg], covars, death)

    # Cox with Kinship matrix
    cpg_model <- coxme(Surv(tte, dead) ~ as.factor(sex) + as.factor(batch) + as.factor(smoking_status) + age + cpg + pack_years +
                    PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + 
                    CD8T + CD4T + NK + Bcell + Gran + (1|Sample_Name), data = df_cpg, varlist = kinship_matrix*2)
    
    results <- extract_coxme_table(cpg_model)
    results_chrom[cpg,] <- results["cpg",]
    i <- i + 1
}

write.table(data.frame("cpg" = rownames(results_chrom), results_chrom), paste0(opt$output_dir, "coxme_survival_ewas_", chrom, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
