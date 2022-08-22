#!/usr/bin/env python
# -*- coding: utf-8 -*-
# By Elena

################################################################
# Script that coordinates steps to perform survival EWAS
################################################################

# Run on cluster
# For VSCode + cluster: export PATH=/home/ebernab3/anaconda3/bin:$PATH

import os
import numpy as np
import pandas as pd
import random
import subprocess
import time
from datetime import datetime
import argparse
from survival_ewas.plot_funcs import *

# USER INPUT
######################################################
######################################################

parser = argparse.ArgumentParser(description='Survival EWAS pipeline.')
parser.add_argument('--coxph_osca', help='Option to run EWAS first fitting CoxPH with age and sex, then OSCA on Martingale residuals.', action="store_true", default=False)
parser.add_argument('--coxph', help='Option to run EWAS one CpG at a time using CoxPH.', action="store_true", default=False)
parser.add_argument('--coxme', help='Option to run EWAS one CpG at a time using CoxME (accounting for relatedness).', action="store_true", default=False)
parser.add_argument('--cpgs', type=str, help='File location of CpGs to subset analysis to.', default=None)
parser.add_argument('--chrom', type=int, help='Chromosome to consider in analysis (for coxph and coxme options). If none chosen, all are run one after the other.', default = 0, choices=range(1,23))
args = parser.parse_args()

####### Parameters
# tf = {True: "T", False: "F"}
if args.chrom == 0:
    chrom = range(1,23)
else:
    chrom = [args.chrom]

# BASE INFO
######################################################
######################################################

mvals_base = "/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/Chromosomes/GS20k_chr%s_mvals.rds"
#cpgs = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_cpgs_common.txt"
#cpgs = "/Cluster_Filespace/Marioni_Group/Elena/gs_osca/data/cpgs_tokeep.txt"
cpgs = "/Cluster_Filespace/Marioni_Group/Elena/gs_osca/data/cpgs_all_831733.txt"
if args.cpgs != None:
    cpgs = args.cpgs # CoxPH Cpgs located here: /Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/survival_ewas/coxph/coxph_survival_cpgs_ews.txt
covars = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_allcovars.tsv"
death_20k = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/biological_age/data_prep/gs_w1w3w4_bage_variables_scaledgrim.tsv"
gs_20k_kinship = "/Cluster_Filespace/Marioni_Group/Elena/data/gs_data/gs_kinship_20k.RDS"
output_dir = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/survival_ewas/%s/"
cpgs_tokeep = "/Cluster_Filespace/Marioni_Group/Elena/gs_osca/data/cpgs_tokeep.txt"

# OSCA files
#osca_base = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_mvals_fixedage_commonsnps"
osca_base = "/Cluster_Filespace/Marioni_Group/Elena/gs_osca/data/mvals-norm20k-18413-831733"
osca_meth = osca_base + ".bod"
osca_samp = osca_base + ".oii"
osca_cpgs = osca_base + ".opi"
osca_qcovs = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_qcovar.txt" # Quantitative covars
osca_ccovs = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_covar_nowave_nosex.txt" # Qualitative covars


# RUN EWAS
######################################################
######################################################

now = datetime.now()
start = time.time()

if args.coxph_osca:
    output_dir = output_dir % "coxph_osca"

    # First run Cox-PH with age and sex as covars & extract Martingale residuals for OSCA
    ##################################################################################################################
    subprocess.call("survival_ewas/coxph_osca.R --death %s --covars %s --output_dir %s" % (death_20k, covars, output_dir), shell  = True)

    # Prepare OSCA files - Run OSCA on Martingale residuals, using remaining covariates + CpG as independent variables
    ##################################################################################################################
    subprocess.call("osca_Linux --befile %s --keep %s --pheno %s --qcovar %s --covar %s --linear --fast-linear --out %s" % (osca_base, output_dir + "coxph_survival_samples.txt", output_dir + "coxph_survival_martingale.txt", osca_qcovs, osca_ccovs, output_dir + "osca_survival_martingale"), shell  = True)

    # Explore results
    ##################################################################################################################
    
    # Extract info of interest
    survival_df = pd.read_csv(output_dir + "osca_survival_martingale.linear", delim_whitespace = True, index_col = 1)
    survival_df["-log10p"] = -np.log10(survival_df["p"])
    survival_df = survival_df.sort_values(by = ["-log10p", "b"], ascending = [False, False])
    survival_df.to_csv(output_dir + "osca_survival_ewas.tsv", sep = "\t", na_rep = "NA")
    ews = 3.6 * 10**(-8)
    survival_sig = survival_df[survival_df["p"] < ews] # 1737 signigicant CpGs
    survival_sig.to_csv(output_dir + "osca_survival_ewas_ews.tsv", sep = "\t", na_rep = "NA")
    
    # Export as list for coxme
    sig_cpgs = open(output_dir + "osca_survival_cpgs_ews.txt", "w")
    sig_cpgs.write("\n".join(survival_sig.index.values.tolist()))
    sig_cpgs.close()

    # Manhattan
    fig, axes = plt.subplots(1, 1, figsize = (10, 3))
    manhattan(values = survival_df["-log10p"], meta = survival_df, ax = axes, col_chr = "Chr", col_bp = "bp", ylabel = "-log10(P-Value)", colors = "custom2")
    axes.set_ylim([0, 60])
    axes.set_title("Survival")
    plt.tight_layout()
    fig.savefig(output_dir + "osca_survival_manhattan.pdf", dpi=300)
    plt.close(fig)


elif args.coxph:
    output_dir = output_dir % "coxph"

    # Run Cox-PH with all covariates per CpG
    ##################################################################################################################
    """
    for c in chrom:
        subprocess.call("survival_ewas/coxph.R --chrom %s --cpgs %s --meth %s --covars %s --death %s --output_dir %s" % (c, cpgs, mvals_base % c, covars, death_20k, output_dir), shell  = True)

    # Explore results (once all chroms have been run, if not comment out)
    ##################################################################################################################
    """
    file_base = "coxph_survival_ewas_%s.tsv"
    survival_df = pd.DataFrame()
    for chrom in range(1,23):
        df_chrom = pd.read_table(output_dir + file_base % chrom, index_col = 0)
        if chrom == 1:
            survival_df = df_chrom
        else:
            survival_df = pd.concat([survival_df, df_chrom])
    print(survival_df)

    # Filter to CpGs of interest
    cpgs_tokeep_df = [i.rstrip() for i in open(cpgs_tokeep, "r").readlines()]
    survival_df = survival_df.loc[survival_df.index.isin(cpgs_tokeep_df)]
    print(survival_df.shape)

    # Add CpG info
    cpgs = pd.read_table(osca_cpgs, delim_whitespace = True, index_col = 1, header = None, names = ["Chrom", "CpG", "Position", "Gene", "Strand"])
    cpgs = cpgs.loc[cpgs.index.isin(cpgs_tokeep_df),]
    survival_df = pd.concat([survival_df, cpgs], axis = 1)
    print(survival_df.shape)

    # Significant hits
    survival_df["-log10p"] = -np.log10(survival_df["p"])
    survival_df = survival_df.sort_values(by = ["-log10p", "HR"], ascending = [False, False])
    survival_df.to_csv(output_dir + "coxph_survival_ewas.tsv", sep = "\t", na_rep = "NA")
    ews = 3.6 * 10**(-8)
    survival_sig = survival_df[survival_df["p"] < ews] # 1182 signigicant CpGs
    print(len(survival_sig.index))
    survival_sig.to_csv(output_dir + "coxph_survival_ewas_ews.tsv", sep = "\t", na_rep = "NA")
    
    # Export as list for coxme
    sig_cpgs = open(output_dir + "coxph_survival_cpgs_ews.txt", "w")
    sig_cpgs.write("\n".join(survival_sig.index.values.tolist()))
    sig_cpgs.close()

    # Manhattan
    fig, axes = plt.subplots(1, 1, figsize = (10, 3))
    manhattan(values = survival_df["-log10p"], meta = survival_df, ax = axes, col_chr = "Chrom", col_bp = "Position", ylabel = "-log10(P-Value)", colors = "custom2")
    axes.set_title("Survival")
    axes.set_ylim([0, 60])
    axes.axhline(y=-np.log10(3.6*(10**-8)), color='dimgray', linestyle='--')
    #axes.set_ylim([0, 340])
    plt.tight_layout()
    fig.savefig(output_dir + "coxph_survival_manhattan.pdf", dpi=300)
    plt.close(fig)
    
elif args.coxme:
    output_dir = output_dir % "coxme"

    # Run Cox_ME with all covariates per CpG
    ##################################################################################################################

    for c in chrom:
        subprocess.call("survival_ewas/coxme.R --chrom %s --cpgs %s --meth %s --covars %s --death %s --kinship %s --output_dir %s" % (c, cpgs, mvals_base % c, covars, death_20k, gs_20k_kinship, output_dir), shell  = True)
    """
    # Explore results (once all chroms have been run, if not comment out)
    ##################################################################################################################
    file_base = "coxme_survival_ewas_%s.tsv"
    survival_df = pd.DataFrame()
    for chrom in range(1,23):
        df_chrom = pd.read_table(output_dir + file_base % chrom, index_col = 0)
        if chrom == 1:
            survival_df = df_chrom
        else:
            survival_df = pd.concat([survival_df, df_chrom])
    
    # Add CpG info
    cpgs = pd.read_table(osca_cpgs, delim_whitespace = True, index_col = 1, header = None, names = ["Chrom", "CpG", "Position", "Gene", "Strand"])
    cpgs = cpgs[cpgs.index.isin(survival_df.index)]
    survival_df = pd.concat([cpgs, survival_df], axis = 1)

    # Significant hits
    survival_df["-log10p"] = -np.log10(survival_df["p"])
    survival_df = survival_df.sort_values(by = ["-log10p", "HR"], ascending = [False, False])
    survival_df.to_csv(output_dir + "coxme_survival_ewas.tsv", sep = "\t", na_rep = "NA")
    ews = 3.6 * 10**(-8)
    survival_sig = survival_df[survival_df["p"] < ews] # 796 signigicant CpGs
    survival_sig.to_csv(output_dir + "coxme_survival_ewas_ews.tsv", sep = "\t", na_rep = "NA")
    """

print("\nDone!")
total_time = time.time() - start
print("\nTotal time: %s minutes (%s hours)." % (round(total_time/60, 3), round((total_time/60)/60, 3)))