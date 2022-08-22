#!/usr/bin/env Rscript 
# -*- coding: utf-8 -*-
# By Elena

library("optparse")
library("survival")

option_list = list(
    make_option(c("--death"), type="character", default=NULL, 
              help="File with death info", metavar="character"),    
    make_option(c("--covars"), type="character", default=NULL, 
              help="File with other covar info", metavar="character"),
    make_option(c("--output_dir"), type="character", default=NULL, 
              help="Directory to output Martingale residuals", metavar="character")        
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# For terminal testing
# death <- read.csv("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/biological_age/data_prep/gs_w1w3w4_bage_variables_scaledgrim.tsv", sep = "\t", row.names = 1)
# covars <- read.csv("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_allcovars.tsv", sep = "\t", row.names = 1)
# output_dir <- "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/survival_ewas/coxph_osca/"

# Import stuff
death <- read.csv(opt$death, sep = "\t", row.names = 1) # Note that 48 have missing values for TTE, hence smaller subset in osca_df
covars <- read.csv(opt$covars, sep = "\t", row.names = 1)
covars <- covars[rownames(death),] # Match order
df <- cbind(death[,c("dead", "tte")], covars[,c("sex", "age")])

# Cox model
cox <- coxph(Surv(tte, dead) ~ age + factor(sex), data=df)
martin <- residuals(cox, type = "martingale")

# OSCA pheno export
osca_df <- data.frame("FID" = names(martin), "IID" = names(martin), "Martingale" = martin)
write.table(osca_df, paste0(output_dir, "coxph_survival_martingale.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")

# Since not all samples are accounted for, make file with samples with no NAs
write.table(osca_df[,c(1,2)], paste0(output_dir, "coxph_survival_samples.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")
