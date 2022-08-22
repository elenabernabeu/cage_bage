#!/usr/bin/env python
# -*- coding: utf-8 -*-
# By Elena

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
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

# EWAS results
osca_ewas = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/survival_ewas/coxph_osca/osca_survival_ewas.tsv", index_col = 0)
coxph_ewas = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/survival_ewas/coxph/coxph_survival_ewas.tsv", index_col = 0)
coxme_ewas = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/survival_ewas/coxme/coxme_survival_ewas.tsv", index_col = 0)

# Output
output_dir = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/survival_ewas/model_comp/"

# Scatterplot, all CpGs
df = pd.concat([osca_ewas["-log10p"], coxph_ewas["-log10p"]], axis = 1)
df.columns = ["osca_-log10p", "coxph_-log10p"]
r = scipy.stats.pearsonr(df["osca_-log10p"], df["coxph_-log10p"])[0]
fig, axes = plt.subplots(1, 1, figsize = (5, 5))
sns.scatterplot(ax = axes, data = df, x = "osca_-log10p", y = "coxph_-log10p", color = "#07B1D8", rasterized = True, linewidth=0)
axes.plot(range(0, 60), range(0, 60), linewidth=1.5, color = "grey")
axes.set_xlabel("CoxPH + OSCA -log10(p)")
axes.set_ylabel("CoxPH -log10(p)")
axes.set_title("OSCA vs CoxPH Mortality EWAS\nr = %s" % round(r,3))
plt.tight_layout()
plt.savefig(output_dir + "osca_coxph_mortality_ewas_comparison.pdf", dpi = 300)
plt.close(fig)

# Scatterplot, CoxPH vs CoxME
ews = 3.6 * 10**(-8)
coxph_ewas_sig = coxph_ewas.loc[coxph_ewas["p"] < ews,]
df = pd.concat([coxph_ewas_sig["-log10p"], coxme_ewas["-log10p"]], axis = 1)
df.columns = ["coxph_-log10p", "coxme_-log10p"]
df.loc[df["coxme_-log10p"] == np.inf, "coxme_-log10p"] = np.nan
df = df[~df["coxme_-log10p"].isna()]
r = scipy.stats.pearsonr(df["coxph_-log10p"], df["coxme_-log10p"])[0]
fig, axes = plt.subplots(1, 1, figsize = (5, 5))
sns.scatterplot(ax = axes, data = df, x = "coxph_-log10p", y = "coxme_-log10p", color = "#07B1D8", rasterized = True, linewidth=0)
axes.plot(range(0, 25), range(0, 25), linewidth=1.5, color = "grey")
axes.set_xlabel("CoxPH -log10(p)")
axes.set_ylabel("CoxME -log10(p)")
axes.set_title("CoxPH vs CoxME Mortality EWAS\nr = %s" % round(r,3))
plt.tight_layout()
plt.savefig(output_dir + "coxph_coxme_mortality_ewas_comparison.pdf", dpi = 300)
plt.close(fig)