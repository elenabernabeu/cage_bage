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
import scipy.stats

# Font config
font = "Galvji"
path = '/Cluster_Filespace/Marioni_Group/Elena/other/fonts/Galvji/Galvji-01.ttf'
fe = matplotlib.font_manager.FontEntry(
    fname= path,
    name=font)
matplotlib.font_manager.fontManager.ttflist.insert(0, fe)
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.family'] = font


##### Comparison of mortality EWAS vs Colicino's 2020 Aging 9 CpGs
##################################################################

ews = 3.6 * 10**(-8)
colicino_cpgs = ["cg07677157", "cg09615688", "cg18424841", "cg17086398", "cg14866069", "cg23666362", "cg12619262", "cg20045320", "cg07839457"]
coxph_ewas = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/survival_ewas/coxph/coxph_survival_ewas.tsv", index_col = 0)
coxph_ewas_sig = coxph_ewas[coxph_ewas["p"] < ews]

##### Overlap

def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

inter = intersection(colicino_cpgs, coxph_ewas_sig.index.values.tolist()) # 1

##### Comparison of mortality EWAS vs Colcino's 2020 Aging 257 CpGs (from non-complete model)
#############################################################################################

colicino_cpgs_basic = pd.read_csv("/Cluster_Filespace/Marioni_Group/Elena/data/random/colicino_aging2020_fdr_cpgs_basicmodel.csv", sep = ";", index_col = 0)
inter_basic = intersection(colicino_cpgs_basic.index, coxph_ewas_sig.index.values.tolist()) # 25
data4both = [i for i in coxph_ewas.index if i in colicino_cpgs_basic.index] # Only 200 of the 257 have data for them in our EWAS
colicino_cpgs_basic["logHR"] = np.log(colicino_cpgs_basic["HRb"])
colicino_cpgs_basic["SE"] = ([np.log(float(i.split(";")[1][:-1])) for i in colicino_cpgs_basic["95% CI"]] - colicino_cpgs_basic["logHR"])/1.96
colicino_cpgs_basic["Z"] = colicino_cpgs_basic["logHR"]/colicino_cpgs_basic["SE"]

##### Prep data for comparison
##############################

df = pd.concat([coxph_ewas[["Z", "p"]], colicino_cpgs_basic[["Z", "p"]]], join = "inner", axis = 1)
df.columns = ["Z", "p", "colicino_Z", "colicino_p"]
df["p"] = -np.log10(df["p"])
df["colicino_p"] = -np.log10(df["colicino_p"])

##### Correlations: overall and just for + and - Z vals
#######################################################

neg = df[df["colicino_Z"] < 0]
pos = df[df["colicino_Z"] > 0]
r_neg = scipy.stats.pearsonr(neg["Z"], neg["colicino_Z"])
r_pos = scipy.stats.pearsonr(pos["Z"], pos["colicino_Z"])

##### Plot comparison
#####################

r_p = scipy.stats.pearsonr(df["p"], df["colicino_p"])
r_z = scipy.stats.pearsonr(df["Z"], df["colicino_Z"])
output_dir = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/survival_ewas/coxph/"
fig, axes = plt.subplots(1, 2, figsize = (8, 4))
sns.scatterplot(ax = axes[0], data = df, x = "Z", y = "colicino_Z", s=15, color = "#07B1D8", label = "r = %s" % (round(r_z[0], 3)), linewidth=0)
sns.kdeplot(data = df, x = "Z", y = "colicino_Z", ax = axes[0], color = "black", levels = 4, thresh=.1, linewidth=0.5, alpha = 0.5)
axes[0].axhline(y=0, color='r', linestyle='-', linewidth = 1)
axes[0].axvline(x=0, color='r', linestyle='-', linewidth = 1)
axes[0].plot(df["Z"], df["Z"], linewidth=1.5, color = "grey")
axes[0].set_xlabel("Z")
axes[0].set_ylabel("Colicino Z")
axes[0].set_title("Z")
leg = axes[0].legend(frameon = False, handletextpad=0, handlelength=0)
for item in leg.legendHandles:
    item.set_visible(False)

sns.scatterplot(ax = axes[1], data = df, x = "p", y = "colicino_p", s=15, color = "#07B1D8", label = "r = %s" % (round(r_p[0], 3)), linewidth=0)
sns.kdeplot(data = df, x = "p", y = "colicino_p", ax = axes[1], color = "black", levels = 4, thresh=.1, linewidth=0.5, alpha = 0.5)
axes[1].plot(df["p"], df["p"], linewidth=1.5, color = "grey")
axes[1].set_ylabel("Colicino -log10(p)")
axes[1].set_xlabel("-log10(p)")
axes[1].set_title("-log10(p)")
leg = axes[1].legend(frameon = False, handletextpad=0, handlelength=0)
for item in leg.legendHandles:
    item.set_visible(False)

plt.suptitle("Colicino Mortality EWAS Comparison (200/257 CpGs)")
plt.tight_layout()
plt.savefig(output_dir + "colicino_mortalityewas_basicmodel_comp.pdf", dpi = 300)
plt.close(fig)

####### Just Z

r_z = scipy.stats.pearsonr(df["Z"], df["colicino_Z"])
output_dir = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/survival_ewas/coxph/"
fig, axes = plt.subplots(1, 1, figsize = (4, 4))
sns.scatterplot(ax = axes, data = df, x = "Z", y = "colicino_Z", s=15, color = "#07B1D8", label = "r = %s" % (round(r_z[0], 3)), linewidth=0)
sns.kdeplot(data = df, x = "Z", y = "colicino_Z", ax = axes, color = "black", levels = 4, thresh=.1, linewidth=0.5, alpha = 0.5)
axes.axhline(y=0, color='r', linestyle='-', linewidth = 1)
axes.axvline(x=0, color='r', linestyle='-', linewidth = 1)
axes.plot(df["Z"], df["Z"], linewidth=1.5, color = "grey")
axes.set_xlabel("Z")
axes.set_ylabel("Colicino Z")
axes.set_title("Colicino et al 2020 Comparison")
leg = axes.legend(frameon = False, handletextpad=0, handlelength=0)
for item in leg.legendHandles:
    item.set_visible(False)

plt.tight_layout()
plt.savefig(output_dir + "colicino_mortalityewas_basicmodel_comp_justZ.pdf", dpi = 300)
plt.close(fig)