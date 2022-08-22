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


# OSCA results
##########################################

age_df = pd.read_csv("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/osca/w1w3w4/gs20k_age.linear", delim_whitespace = True, index_col = 1)
age2_df = pd.read_csv("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/osca/w1w3w4/gs20k_age2.linear", delim_whitespace = True, index_col = 1)
#cpg2_df = pd.read_csv("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/osca/w1w3w4/gs20k_cpg2.linear", delim_whitespace = True, index_col = 1)


# Filter to QC'd CpGs
##########################################

cpgs = [i.rstrip() for i in open("/Cluster_Filespace/Marioni_Group/Elena/gs_osca/data/cpgs_tokeep.txt", "r").readlines()]
age_df = age_df.loc[age_df.index.isin(cpgs),]
age2_df = age2_df.loc[age2_df.index.isin(cpgs),]


# Sort by p-val and export 
##########################################

age_df["-log10p"] = -np.log10(age_df["p"])
age2_df["-log10p"] = -np.log10(age2_df["p"])
#cpg2_df["-log10p"] = -np.log10(cpg2_df["p"])
age_df = age_df.sort_values(by = ["-log10p", "b"], ascending = [False, False])
age2_df = age2_df.sort_values(by = ["-log10p", "b"], ascending = [False, False])
#cpg2_df = cpg2_df.sort_values(by = ["-log10p", "b"], ascending = [False, False])
age_df.to_csv("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/ewas_results/age_ewas.tsv", sep = "\t", na_rep = "NA")
age2_df.to_csv("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/ewas_results/age2_ewas.tsv", sep = "\t", na_rep = "NA")
#cpg2_df.to_csv("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/ewas_results/cpg2_ewas.tsv", sep = "\t", na_rep = "NA")


# Just sig
##########################################

ews = 3.6 * 10**(-8) # Epigenome-wide significance

age_sig = age_df[age_df["p"] < ews] # 99832 sig CpGs
age2_sig = age2_df[age2_df["p"] < ews] # 30114 sig CpGs
#cpg2_sig = cpg2_df[cpg2_df["p"] < ews]
siginboth_age2 = list(set(age_sig.index).intersection(set(age2_sig.index))) # 18826 sig CpGs for both age and age^2
#siginboth_cpg2 = list(set(age_sig.index).intersection(set(cpg2_sig.index)))

# Output dfs
age_sig.to_csv("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/ewas_results/age_ewas_epigenomewidesig.tsv", sep = "\t", na_rep = "NA")
age2_sig.to_csv("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/ewas_results/age2_ewas_epigenomewidesig.tsv", sep = "\t", na_rep = "NA")
#cpg2_sig.to_csv("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/ewas_results/cpg2_ewas_epigenomewidesig.tsv", sep = "\t", na_rep = "NA")

# Output lists of CpGs
cpgs_age = open("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/ewas_results/cpgs_age_epigenomewidesig.txt", "w")
cpgs_age.write("\n".join(age_sig.index))
cpgs_age.close() 
cpgs_age2 = open("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/ewas_results/cpgs_age2_epigenomewidesig.txt", "w")
cpgs_age2.write("\n".join(age2_sig.index))
cpgs_age2.close()
#cpgs_cpg2 = open("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/ewas_results/cpgs_cpg2_epigenomewidesig.txt", "w")
#cpgs_cpg2.write("\n".join(cpg2_sig.index))
#cpgs_cpg2.close()
cpgs_both_age2 = open("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/ewas_results/cpgs_both_epigenomewidesig.txt", "w")
cpgs_both_age2.write("\n".join(siginboth_age2))
cpgs_both_age2.close() 
#cpgs_both_cpg2 = open("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/ewas_results/cpgs_both_epigenomewidesig_cpg2.txt", "w")
#cpgs_both_cpg2.write("\n".join(siginboth_cpg2))
#cpgs_both_cpg2.close() 


# Subsets
##########################################

output_dir = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/subsets/"
cpgs_for_subsets = [i.rstrip() for i in open("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_cpgs_common.txt", "r").readlines()]
age_sub = [1, 2, 3, 4, 5, 10, 15, 20, 35, 50, 75, 100, 200, 300]
age2_sub = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1250, 1500, 2000, 3000, 5000, 7500, 10000, 12500, 15000, 20000]
#cpg2_sub = [1, 2, 3, 4, 5, 10, 15, 20, 35, 50, 75, 100, 200, 300]

for limit in age_sub:
    sub = age_df.loc[age_df.index.isin(cpgs_for_subsets),]
    sub = age_df.loc[age_df.index[0:(limit*1000)],]
    sub.to_csv(output_dir + "age_ewas_%sK_gs20k.tsv" % limit, sep = "\t", na_rep = "NA")

for limit in age2_sub:
    sub = age2_df.loc[age2_df.index.isin(cpgs_for_subsets),]
    sub = age2_df.loc[age2_df.index[0:limit],]
    sub.to_csv(output_dir + "age2_ewas_%s_gs20k.tsv" % limit, sep = "\t", na_rep = "NA")
"""
for limit in cpg2_sub:
    sub = cpg2_df.loc[cpg2_df.index[0:limit],]
    sub.to_csv(output_dir + "cpg2_ewas_%s_gs20k.tsv" % limit, sep = "\t", na_rep = "NA")
"""

# Manhattan plot
##########################################

# Remove inf and set to 310 (approx -log10 of smallest number representable)
age_df.loc[age_df["-log10p"] == float(inf), "-log10p"] = 310
age2_df.loc[age2_df["-log10p"] == float(inf), "-log10p"] = 310
#cpg2_df.loc[cpg2_df["-log10p"] == float(inf), "-log10p"] = 310

fig, axes = plt.subplots(1, 1, figsize = (10, 3))
manhattan(values = age_df["-log10p"], meta = age_df, ax = axes, col_chr = "Chr", col_bp = "bp", ylabel = "-log10(P-Value)", colors = "custom2")
axes.set_title("Age")
axes.axhline(y=-np.log10(3.6*(10**-8)), color='dimgray', linestyle='--')
plt.tight_layout()
fig.savefig("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/manhattan/manhattan_age_gs20k.pdf", dpi=300)
plt.close(fig)

fig, axes = plt.subplots(1, 1, figsize = (10, 3))
manhattan(values = age2_df["-log10p"], meta = age2_df, ax = axes, col_chr = "Chr", col_bp = "bp", ylabel = "-log10(P-Value)", colors = "custom2")
axes.set_title("Age^2")
axes.axhline(y=-np.log10(3.6*(10**-8)), color='dimgray', linestyle='--')
plt.tight_layout()
fig.savefig("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/manhattan/manhattan_age2_gs20k.pdf", dpi=300)
plt.close(fig)
"""
fig, axes = plt.subplots(1, 1, figsize = (10, 3))
manhattan(values = cpg2_df["-log10p"], meta = cpg2_df, ax = axes, col_chr = "Chr", col_bp = "bp", ylabel = "-log10(P-Value)", colors = "custom2")
axes.set_title("Age")
plt.tight_layout()
fig.savefig("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/manhattan/manhattan_cpg2_gs20k.pdf", dpi=300)
plt.close(fig)
"""
"""
age_df_no23 = age_df[age_df["Chr"] != 23]
age2_df_no23 = age2_df[age2_df["Chr"] != 23]
age_df_no23.loc[age_df_no23["-log10p"] == float(inf), "-log10p"] = 310
age2_df_no23.loc[age2_df_no23["-log10p"] == float(inf), "-log10p"] = 310

fig, axes = plt.subplots(1, 1, figsize = (10, 3))
manhattan(values = age_df_no23["-log10p"], meta = age_df_no23, ax = axes, col_chr = "Chr", col_bp = "bp", ylabel = "-log10(P-Value)", colors = "custom2")
axes.set_title("Age")
plt.tight_layout()
fig.savefig("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/manhattan/manhattan_age_gs20k_no23.pdf", dpi=300)
plt.close(fig)

fig, axes = plt.subplots(1, 1, figsize = (10, 3))
manhattan(values = age2_df_no23["-log10p"], meta = age2_df_no23, ax = axes, col_chr = "Chr", col_bp = "bp", ylabel = "-log10(P-Value)", colors = "custom2")
axes.set_title("Age^2")
plt.tight_layout()
fig.savefig("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/manhattan/manhattan_age2_gs20k_no23.pdf", dpi=300)
plt.close(fig)

"""
# CpG trajectories
##########################################

# Some CpGs to investigate across ages
cpgs_age = ["cg21572722", "cg16717122", "cg06493994", "cg11071401", "cg07547549"] # Genes: ELOVL2, SCG3, SCGN, CACNA1G, SLC12A5
cpgs_age2 = ["cg07850154", "cg00329615", "cg21184711", "cg05412028", "cg21878650"] # Genes: RNF180, IGSF11, GRM2, CADPS2, ABCC4
#cpgs_cpg2 = ["cg16867657", "cg07850154", "cg07955995", "cg22353329"]

# First take OG m-val matrix and filter to sig CpGs for age and age^2
#cpgs = list(set(age_sig.index.values.tolist() + age2_sig.index.values.tolist()))
#df = dt.fread("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_mvals.txt", sep = " ")
#df = df.to_pandas()
#df = df.set_index("IID")
#df = df.drop("FID", axis = 1)
#df = df[cpgs]
#df.to_csv("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/cpg_trajectories/meth_sig.tsv", sep = "\t")
df = dt.fread("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/cpg_trajectories/meth_sig.tsv", sep = "\t")
df = df.to_pandas()
df = df.set_index("IID")
df_beta = m2beta(df)

# Target file for ages
target = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/data/gs_data/gs20ktargets_fixedage.tsv", index_col = 0)

output_dir = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/cpg_trajectories/"
for cpg in cpgs_age:
    print("##### Plots for CpGs associated to age linearly")
    print("Making trajectory plot for %s..." % cpg)
    df_s = pd.concat([target["age"], df[cpg]], axis = 1)
    cpg_trajectory(df = df_s, cpg = cpg, gene_name = age_sig.at[cpg, "Gene"].split(";")[0], output_dir = output_dir, cpg_col = cpg, age_col = "age", color = "#07B1D8", wave = "gs20k_age")
    df_s_b = pd.concat([target["age"], df_beta[cpg]], axis = 1)
    cpg_trajectory(df = df_s_b, cpg = cpg, gene_name = age_sig.at[cpg, "Gene"].split(";")[0], output_dir = output_dir, cpg_col = cpg, age_col = "age", color = "#07B1D8", wave = "gs20k_age_beta")


for cpg in cpgs_age2:
    print("##### Plots for CpGs associated to age quadratically")
    print("Making trajectory plot for %s..." % cpg)
    df_s = pd.concat([target["age"], df[cpg]], axis = 1)
    cpg_trajectory(df = df_s, cpg = cpg, gene_name = age2_sig.at[cpg, "Gene"].split(";")[0], output_dir = output_dir, cpg_col = cpg, age_col = "age", color = "#07B1D8", wave = "gs20k_age2")
    df_s_b = pd.concat([target["age"], df_beta[cpg]], axis = 1)
    cpg_trajectory(df = df_s_b, cpg = cpg, gene_name = age2_sig.at[cpg, "Gene"].split(";")[0], output_dir = output_dir, cpg_col = cpg, age_col = "age", color = "#07B1D8", wave = "gs20k_age2_beta")
"""
for cpg in cpgs_cpg2:
    print("##### Plots for CpGs associated to age quadratically")
    print("Making trajectory plot for %s..." % cpg)
    df_s = pd.concat([target["age"], df[cpg]], axis = 1)
    cpg_trajectory(df = df_s, cpg = cpg, gene_name = age2_sig.at[cpg, "Gene"].split(";")[0], output_dir = output_dir, cpg_col = cpg, age_col = "age", color = "#07B1D8", wave = "gs20k_cpg2")
    df_s_b = pd.concat([target["age"], df_beta[cpg]], axis = 1)
    cpg_trajectory(df = df_s_b, cpg = cpg, gene_name = age2_sig.at[cpg, "Gene"].split(";")[0], output_dir = output_dir, cpg_col = cpg, age_col = "age", color = "#07B1D8", wave = "gs20k_cpg2_beta")
"""
