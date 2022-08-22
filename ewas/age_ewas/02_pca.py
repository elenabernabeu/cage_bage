#!/usr/bin/env python
# -*- coding: utf-8 -*-
# By Elena

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import datatable as dt
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mc
import seaborn as sns
import statsmodels.api as sm
import colorsys
import math

# Aesthetic stuff
################################################################

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

def lighten_color(color, amount=1.2):  
    # --------------------- SOURCE: @IanHincks ---------------------
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


def set_colors(N, name = "default"):
    colors = []
    if name == "default":
        colors = sns.color_palette("Set2", 6)
    elif name == "custom1":
        colors = ["#DD1155", "#7D82B8", "#FCBA04", "#63C7B2", "#45364B"] #Hotpink, lilacky-blue, teal, yellow, dark purple
    elif name == "custom2":
        colors = ["#47C3F4", "#FF9429", "#63C7B2", "#FFCE00", "#EF553C", "#2C819B", "#07B1D8", "#C5E3EA", "#326B64", "#FFD15F", "#EA4C15"] #Lightblue, orange, teal, yellow, red, dark blue, midtone blue, very light blue, 
    colors = (colors*int(math.ceil(N/len(colors))))[:N]
    return colors

# Import data
################################################################

df = dt.fread("/Cluster_Filespace/Marioni_Group/Elena/gs_osca/data/gs20k_mvals.txt", sep = " ")
df = df.to_pandas()
df = df.set_index("IID")
df = df.drop("FID", axis = 1)
names = df.index.values.tolist()

# Standardize
################################################################

df = StandardScaler().fit_transform(df)

# PCA
################################################################

pca = PCA(n_components=100)
pcs = pca.fit_transform(df)
pca_df = pd.DataFrame(data = pcs, columns = ["PC%s" % i for i in range(1,101)], index = names)
expvar_df = pd.DataFrame(data = pca.explained_variance_ratio_, index = ["PC%s" % i for i in range(1,101)], columns = ["ExpVar"])

# Export results
################################################################

pca_df.to_csv("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_PCA.tsv", sep = "\t", index_label = "ID")
expvar_df.to_csv("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_PCA_explainedvariance.tsv", sep = "\t", index_label = "PC")

# Plot coloring by wave
################################################################

pca_df = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/gs_osca/data/norm20k_18413_752722_PCA.tsv", index_col = 0)
expvar_df = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/gs_osca/data/norm20k_18413_752722_PCA_explainedvariance.tsv", index_col = 0)

# Add target info
target = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/data/gs_data/gs20ktargets.tsv", index_col = 0)
pca_df = pd.concat([pca_df, target], axis = 1)

# Plot
col = set_colors(3, "custom2")
sns.set_palette(sns.color_palette(col))
fig, axes = plt.subplots(2, 3, figsize = (12, 8))
# PC1 vs PC2
sns.scatterplot(ax = axes[0,0], data = pca_df, x = "PC1", y = "PC2", hue = "Set", rasterized = True)
sns.kdeplot(ax = axes[1,0], data = pca_df, x = "PC1", y = "PC2", hue = "Set", rasterized = True)
axes[0,0].set_ylabel("PC2 (%s%%)" % round(expvar_df.at["PC2", "ExpVar"]*100, 2))
axes[0,0].set_xlabel("PC1 (%s%%)" % round(expvar_df.at["PC1", "ExpVar"]*100, 2))
axes[0,0].set_title("PC1 vs PC2")
axes[1,0].set_ylabel("PC2 (%s%%)" % round(expvar_df.at["PC2", "ExpVar"]*100, 2))
axes[1,0].set_xlabel("PC1 (%s%%)" % round(expvar_df.at["PC1", "ExpVar"]*100, 2))
# PC1 vs PC3
sns.scatterplot(ax = axes[0,1], data = pca_df, x = "PC1", y = "PC3", hue = "Set", rasterized = True)
sns.kdeplot(ax = axes[1,1], data = pca_df, x = "PC1", y = "PC3", hue = "Set", rasterized = True)
axes[0,1].set_ylabel("PC3 (%s%%)" % round(expvar_df.at["PC3", "ExpVar"]*100, 2))
axes[0,1].set_xlabel("PC1 (%s%%)" % round(expvar_df.at["PC1", "ExpVar"]*100, 2))
axes[0,1].set_title("PC1 vs PC3")
axes[1,1].set_ylabel("PC3 (%s%%)" % round(expvar_df.at["PC3", "ExpVar"]*100, 2))
axes[1,1].set_xlabel("PC1 (%s%%)" % round(expvar_df.at["PC1", "ExpVar"]*100, 2))
# PC2 vs PC3
sns.scatterplot(ax = axes[0,2], data = pca_df, x = "PC2", y = "PC3", hue = "Set", rasterized = True)
sns.kdeplot(ax = axes[1,2], data = pca_df, x = "PC2", y = "PC3", hue = "Set", rasterized = True)
axes[0,2].set_ylabel("PC3 (%s%%)" % round(expvar_df.at["PC3", "ExpVar"]*100, 2))
axes[0,2].set_xlabel("PC2 (%s%%)" % round(expvar_df.at["PC2", "ExpVar"]*100, 2))
axes[0,2].set_title("PC2 vs PC3")
axes[1,2].set_ylabel("PC3 (%s%%)" % round(expvar_df.at["PC3", "ExpVar"]*100, 2))
axes[1,2].set_xlabel("PC2 (%s%%)" % round(expvar_df.at["PC2", "ExpVar"]*100, 2))
# Export
plt.tight_layout()
plt.savefig("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/pca/gsmeth20k_pca.pdf", dpi = 300)
plt.close(fig)









