#!/usr/bin/env python
# -*- coding: utf-8 -*-
# By Elena

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as mpatches
import matplotlib.colors as mc
import seaborn as sns
import colorsys
import datatable as dt
from plot_funcs import *
import math

# For VSCode + cluster: export PATH=/home/ebernab3/anaconda3/bin:$PATH

def m2beta(m):
    beta = (2**m)/(2**m+1) 
    return beta

# Import lm results
##########################################

age2_df = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/age2_lm/w1w3w4/age2_lm_allchrom.tsv", index_col = 0)
age_df = age2_df.sort_values(by = ["age_log10p", "age_beta"], ascending = [False, False])
age2_df.loc[age2_df["age2_log10p"] == float(inf), "age2_log10p"] = 320
age_df.loc[age_df["age_log10p"] == float(inf), "age_log10p"] = 320

# Just sig
##########################################

ews = 3.6 * 10**(-8) # Epigenome-wide significance
age2_sig = age2_df[age2_df["age2_p"] < ews] # 30114 sig CpGs
age_sig = age_df[age_df["age_p"] < ews] # 106442 sig CpGs 
age2_sig.to_csv("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/ewas_results/age2_lm_epigenomewidesig.tsv", sep = "\t", na_rep = "NA")
age_sig.to_csv("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/ewas_results/age_lm_epigenomewidesig.tsv", sep = "\t", na_rep = "NA")

# Sig in both
siginboth = list(set(age2_sig.index).intersection(set(age_sig.index))) # 20692

# Output lists of cpgs
cpgs_age2 = open("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/ewas_results/cpgs_age2_epigenomewidesig_lm.txt", "w")
cpgs_age2.write("\n".join(age2_sig.index))
cpgs_age2.close()
cpgs_age = open("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/ewas_results/cpgs_age_epigenomewidesig_lm.txt", "w")
cpgs_age.write("\n".join(age_sig.index))
cpgs_age.close()
#cpgs_both_age = open("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/ewas_results/cpgs_both_cpg2_age_epigenomewidesig_lm.txt", "w")
#cpgs_both_age.write("\n".join(siginboth_age))
#cpgs_both_age.close() 
#cpgs_both_age2 = open("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/ewas_results/cpgs_both_cpg2_age2_epigenomewidesig_lm.txt", "w")
#cpgs_both_age2.write("\n".join(siginboth_age2))
#cpgs_both_age2.close() 

# Subsets
##########################################

cpgs_for_subsets = [i.rstrip() for i in open("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_cpgs_common.txt", "r").readlines()]
output_dir = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/subsets/"
age2_sub = [100,200,300,400,500,600,700,800,900,1000,1250,1500,2000,3000,5000,7500,10000,12500,15000,20000]
age_sub = [1,2,3,4,5,10,15,20,35,50,75,100,200,300]

for i in age2_sub:
    sub = age2_df.loc[age2_df.index.isin(cpgs_for_subsets),]
    sub = sub.loc[sub.index[0:int(i)],]
    print(sub.shape)
    sub.to_csv(output_dir + "age2_lm_%s_gs20k.tsv" % i, sep = "\t", na_rep = "NA")

for i in age_sub:
    sub = age_df.loc[age_df.index.isin(cpgs_for_subsets),]
    sub = sub.loc[sub.index[0:int(i*1000)],]
    print(sub.shape)
    sub.to_csv(output_dir + "age_lm_%sK_gs20k.tsv" % i, sep = "\t", na_rep = "NA")

# Manhattan plot
##########################################

fig, axes = plt.subplots(1, 1, figsize = (10, 3))
manhattan(values = age2_df["age2_log10p"], meta = age2_df, ax = axes, col_chr = "chr", col_bp = "pos", ylabel = "-log10(P-Value)", colors = "custom2")
axes.axhline(y=-np.log10(3.6*(10**-8)), color='dimgray', linestyle='--')
axes.set_ylim([0, 340])
axes.set_title("CpG ~ Age^2")
plt.tight_layout()
fig.savefig("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/manhattan/manhattan_age2_gs20k_lm.pdf", dpi=300)
plt.close(fig)

fig, axes = plt.subplots(1, 1, figsize = (10, 3))
manhattan(values = age_df["age_log10p"], meta = age_df, ax = axes, col_chr = "chr", col_bp = "pos", ylabel = "-log10(P-Value)", colors = "custom2")
axes.axhline(y=-np.log10(3.6*(10**-8)), color='dimgray', linestyle='--')
axes.set_ylim([0, 340])
axes.set_title("CpG ~ Age")
plt.tight_layout()
fig.savefig("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/manhattan/manhattan_age_gs20k_lm.pdf", dpi=300)
plt.close(fig)


# CpG trajectories
##########################################

cpgs_age2_lm = ["cg13108341", "cg00329615", "cg17621438", "cg21184711"] # DNAH9, IGSF11, RNF180, CADPS2
cpgs_age_lm = ["cg01873886", "cg04875128", "cg07955995", "cg05404236", "cg16867657"] # BMP4, OTUD7A, KLF14, IRS2, ELOVL2
cpgs_lm = list(set(cpgs_age2_lm + cpgs_age_lm))
#df = dt.fread("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/cpg_trajectories/meth_sig.tsv", sep = "\t")
df = dt.fread("/Cluster_Filespace/Marioni_Group/Elena/gs_osca/data/mvals-norm20k-18413-831733.txt", sep = " ")
df = df[:,["IID"] + cpgs_lm]
df = df.to_pandas()
df = df.set_index("IID")
df_beta = m2beta(df)

# Target file for ages
target = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/data/gs_data/gs20ktargets_fixedage.tsv", index_col = 0)
output_dir = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/cpg_trajectories/"

for cpg in cpgs_age2_lm:
    print("##### Plots for CpGs associated to age quadratically")
    print("Making trajectory plot for %s..." % cpg)
    df_s = pd.concat([target["age"], df[cpg]], axis = 1)
    cpg_trajectory_nodens(df = df_s, cpg = cpg, gene_name = age2_sig.at[cpg, "UCSC_RefGene_Name"].split(";")[0], output_dir = output_dir, cpg_col = cpg, age_col = "age", color = "#07B1D8", wave = "gs20k_age2_lm")
    df_s_b = pd.concat([target["age"], df_beta[cpg]], axis = 1)
    cpg_trajectory_nodens(df = df_s_b, cpg = cpg, gene_name = age2_sig.at[cpg, "UCSC_RefGene_Name"].split(";")[0], output_dir = output_dir, cpg_col = cpg, age_col = "age", color = "#07B1D8", wave = "gs20k_age2_beta_lm")

for cpg in cpgs_age_lm:
    print("##### Plots for CpGs associated to age quadratically")
    print("Making trajectory plot for %s..." % cpg)
    df_s = pd.concat([target["age"], df[cpg]], axis = 1)
    cpg_trajectory_nodens(df = df_s, cpg = cpg, gene_name = age_sig.at[cpg, "UCSC_RefGene_Name"].split(";")[0], output_dir = output_dir, cpg_col = cpg, age_col = "age", color = "#07B1D8", wave = "gs20k_age_lm")
    df_s_b = pd.concat([target["age"], df_beta[cpg]], axis = 1)
    cpg_trajectory_nodens(df = df_s_b, cpg = cpg, gene_name = age_sig.at[cpg, "UCSC_RefGene_Name"].split(";")[0], output_dir = output_dir, cpg_col = cpg, age_col = "age", color = "#07B1D8", wave = "gs20k_age_beta_lm")
