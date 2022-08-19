#!/usr/bin/env python
# -*- coding: utf-8 -*-
# By Elena

# For VSCode + cluster: export PATH=/home/ebernab3/anaconda3/bin:$PATH

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

def m2beta(m):
    beta = (2**m)/(2**m+1) 
    return beta

# Import lm results
##########################################

cpg2_df_nonscaled = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/cpg2_lm/w1w3w4/cpg2_lm_allchrom.tsv", index_col = 0)
cpg2_df = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/cpg2_lm/w1w3w4/cpg2_lm_allchrom_scaled.tsv", index_col = 0)
cpg2_df["cpg_beta_abs"] = abs(cpg2_df["cpg_beta"])
cpg2_df["cpg2_beta_abs"] = abs(cpg2_df["cpg2_beta"])
cpg2_df_nonscaled["cpg_beta_abs"] = abs(cpg2_df_nonscaled["cpg_beta"])
cpg2_df_nonscaled["cpg2_beta_abs"] = abs(cpg2_df_nonscaled["cpg2_beta"])
cpg2_df = cpg2_df.sort_values(by = ["cpg2_log10p", "cpg2_beta_abs"], ascending = [False, False])
cpg2_df_nonscaled = cpg2_df_nonscaled.sort_values(by = ["cpg2_log10p", "cpg2_beta_abs"], ascending = [False, False])
cpg_df = cpg2_df.sort_values(by = ["cpg_log10p", "cpg_beta_abs"], ascending = [False, False])
cpg_df_nonscaled = cpg2_df_nonscaled.sort_values(by = ["cpg_log10p", "cpg_beta_abs"], ascending = [False, False])
cpg2_df.loc[cpg2_df["cpg2_log10p"] == float(inf), "cpg2_log10p"] = 320
cpg_df.loc[cpg_df["cpg_log10p"] == float(inf), "cpg_log10p"] = 320

# Just sig
##########################################

ews = 3.6 * 10**(-8) # Epigenome-wide significance
cpg2_sig = cpg2_df[cpg2_df["cpg2_p"] < ews] # 137915 sig CpGs
cpg_sig = cpg_df[cpg_df["cpg_p"] < ews] # 99832 sig CpGs
cpg2_sig.to_csv("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/ewas_results/cpg2_lm_scaled_epigenomewidesig.tsv", sep = "\t", na_rep = "NA")
cpg_sig.to_csv("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/ewas_results/cpg_lm_scaled_epigenomewidesig.tsv", sep = "\t", na_rep = "NA")

# Sort cpg2 hits by beta instead of pval and export
cpg2_df_beta = cpg2_df.sort_values(by = ["cpg2_beta_abs"], ascending = [False])

# Compare to age and age^2 EWAS OSCA
age_sig = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/ewas_results/age_ewas_epigenomewidesig.tsv", index_col = 0)
#age2_sig = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/ewas_results/age2_ewas_epigenomewidesig.tsv", index_col = 0)
#cpg2_osca_sig = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/ewas_results/cpg2_ewas_epigenomewidesig.tsv", index_col = 0)
siginboth_age = list(set(age_sig.index).intersection(set(cpg_sig.index))) # 48312
#siginboth_age2 = list(set(age2_sig.index).intersection(set(cpg2_sig.index))) # 16167
#siginboth_cpg2 = list(set(cpg2_osca_sig.index).intersection(set(cpg2_sig.index)))

# Sig in both
siginboth = list(set(cpg2_sig.index).intersection(set(cpg_sig.index))) # 48312

# Output lists of cpgs
cpgs_cpg2 = open("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/ewas_results/cpgs_cpg2_epigenomewidesig_lm.txt", "w")
cpgs_cpg2.write("\n".join(cpg2_sig.index))
cpgs_cpg2.close()
cpgs_cpg = open("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/ewas_results/cpgs_cpg_epigenomewidesig_lm.txt", "w")
cpgs_cpg.write("\n".join(cpg_sig.index))
cpgs_cpg.close()
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
cpg2_sub = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 2, 3, 4, 5, 10, 15, 20, 35, 50, 75, 100, 200, 300]
cpg_sub = [1,2,3,4,5,10,15,20,35,50,75,100,200,300]

for i in cpg2_sub:
    sub = cpg2_df.loc[cpg2_df.index.isin(cpgs_for_subsets),]
    print(sub.shape)
    sub = sub.loc[sub.index[0:int(i*1000)],]
    print(sub.shape)
    sub.to_csv(output_dir + "cpg2_lm_%sK_gs20k_scaled.tsv" % i, sep = "\t", na_rep = "NA")

for i in cpg2_sub:
    sub = cpg2_df_beta.loc[cpg2_df_beta.index.isin(cpgs_for_subsets),]
    print(sub.shape)
    sub = sub.loc[sub.index[0:int(i*1000)],]
    print(sub.shape)
    sub.to_csv(output_dir + "cpg2_lm_%sK_gs20k_scaled_betasort.tsv" % i, sep = "\t", na_rep = "NA")

for i in cpg_sub:
    sub = cpg_df.loc[cpg_df.index.isin(cpgs_for_subsets),]
    print(sub)
    sub = sub.loc[sub.index[0:int(i*1000)],]
    print(sub.shape)
    sub.to_csv(output_dir + "cpg_lm_%sK_scaled_gs20k.tsv" % i, sep = "\t", na_rep = "NA")

# Manhattan plot
##########################################

fig, axes = plt.subplots(1, 1, figsize = (10, 3))
manhattan(values = cpg2_df["cpg2_log10p"], meta = cpg2_df, ax = axes, col_chr = "chr", col_bp = "pos", ylabel = "-log10(P-Value)", colors = "custom2")
axes.axhline(y=-np.log10(3.6*(10**-8)), color='dimgray', linestyle='--')
axes.set_ylim([0, 340])
axes.set_title("Age ~ CpG^2")
plt.tight_layout()
fig.savefig("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/manhattan/manhattan_cpg2_gs20k_lm.pdf", dpi=300)
plt.close(fig)

fig, axes = plt.subplots(1, 1, figsize = (10, 3))
manhattan(values = cpg_df["cpg_log10p"], meta = cpg_df, ax = axes, col_chr = "chr", col_bp = "pos", ylabel = "-log10(P-Value)", colors = "custom2")
axes.axhline(y=-np.log10(3.6*(10**-8)), color='dimgray', linestyle='--')
axes.set_ylim([0, 340])
axes.set_title("Age ~ CpG")
plt.tight_layout()
fig.savefig("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/manhattan/manhattan_cpg_gs20k_lm.pdf", dpi=300)
plt.close(fig)

# CpG trajectories
##########################################

cpgs_cpg2_lm = ["cg11084334", "cg15996534", "cg23527621", "cg08935541"] # LHFPL4, LOC134466, ECE2, LINC00899
cpgs_cpg2_lm_betasort = ["cg14858786", "cg01660407", "cg12589308", "cg15600935"] # PLEKHG6, CHST10, ATP8B3, CUX1
# cpgs_cpg_lm = ["cg16867657", "cg22454769", "cg09018739", "cg01673856"] # ELOVL2, FHL2, CPNE2, TCF7L2
cpgs_cpg_lm = ["cg16867657", "cg08097417", "cg24724428", "cg12841266"] # ELOVL2, KLF14, ELOVL2, LHFPL4
cpgs_lm = list(set(cpgs_cpg2_lm + cpgs_cpg_lm + cpgs_cpg2_lm_betasort))
#df = dt.fread("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/cpg_trajectories/meth_sig.tsv", sep = "\t")
df = dt.fread("/Cluster_Filespace/Marioni_Group/Elena/gs_osca/data/mvals-norm20k-18413-831733.txt", sep = " ")
df = df[:,["IID"] + cpgs_lm]
df = df.to_pandas()
df = df.set_index("IID")
df_beta = m2beta(df)

# Target file for ages
target = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/data/gs_data/gs20ktargets_fixedage.tsv", index_col = 0)
output_dir = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/cpg_trajectories/"

print("##### Plots for CpGs associated to age quadratically")
for cpg in cpgs_cpg2_lm:
    print("Making trajectory plot for %s..." % cpg)
    df_s = pd.concat([target["age"], df[cpg]], axis = 1)
    cpg_trajectory_nodens(df = df_s, cpg = cpg, gene_name = cpg2_sig.at[cpg, "UCSC_RefGene_Name"].split(";")[0], output_dir = output_dir, cpg_col = cpg, age_col = "age", color = "#07B1D8", wave = "gs20k_cpg2_lm")
    df_s_b = pd.concat([target["age"], df_beta[cpg]], axis = 1)
    cpg_trajectory_nodens(df = df_s_b, cpg = cpg, gene_name = cpg2_sig.at[cpg, "UCSC_RefGene_Name"].split(";")[0], output_dir = output_dir, cpg_col = cpg, age_col = "age", color = "#07B1D8", wave = "gs20k_cpg2_beta_lm")

print("##### Plots for CpGs associated to age quadratically (beta sorted)")
for cpg in cpgs_cpg2_lm_betasort:
    print("Making trajectory plot for %s..." % cpg)
    df_s = pd.concat([target["age"], df[cpg]], axis = 1)
    cpg_trajectory_nodens(df = df_s, cpg = cpg, gene_name = cpg2_sig.at[cpg, "UCSC_RefGene_Name"].split(";")[0], output_dir = output_dir, cpg_col = cpg, age_col = "age", color = "#07B1D8", wave = "gs20k_cpg2_lm")
    df_s_b = pd.concat([target["age"], df_beta[cpg]], axis = 1)
    cpg_trajectory_nodens(df = df_s_b, cpg = cpg, gene_name = cpg2_sig.at[cpg, "UCSC_RefGene_Name"].split(";")[0], output_dir = output_dir, cpg_col = cpg, age_col = "age", color = "#07B1D8", wave = "gs20k_cpg2_beta_lm")

print("##### Plots for CpGs associated to age linearly")
for cpg in cpgs_cpg_lm:
    print("Making trajectory plot for %s..." % cpg)
    df_s = pd.concat([target["age"], df[cpg]], axis = 1)
    cpg_trajectory_nodens(df = df_s, cpg = cpg, gene_name = cpg_sig.at[cpg, "UCSC_RefGene_Name"].split(";")[0], output_dir = output_dir, cpg_col = cpg, age_col = "age", color = "#07B1D8", wave = "gs20k_cpg_lm")
    df_s_b = pd.concat([target["age"], df_beta[cpg]], axis = 1)
    cpg_trajectory_nodens(df = df_s_b, cpg = cpg, gene_name = cpg_sig.at[cpg, "UCSC_RefGene_Name"].split(";")[0], output_dir = output_dir, cpg_col = cpg, age_col = "age", color = "#07B1D8", wave = "gs20k_cpg_beta_lm")

