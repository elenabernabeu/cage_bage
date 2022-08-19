#!/usr/bin/env python
# -*- coding: utf-8 -*-
# By Elena

import os
import numpy as np
import scipy.stats
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as mpatches
import matplotlib.colors as mc
import seaborn as sns
import statsmodels.api as sm
import colorsys
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import datatable as dt
import math

################################################################
# Functions to explore data prior to fitting elnet models
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

# Aesthetic functions
################################################################

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


# Age distribution
################################################################

def age_dist(target, output_dir, cohort_col, age_col, external, name):
    if external == True:
        df = pd.read_table(target, index_col = 0)
        df["New_Cohort"] = df[cohort_col]
    else:
        df = pd.read_table(target, index_col = 1)
    df = df[~df[age_col].isnull()].copy()

    if external == True:
        # Fuse GS observations
        df.loc[df[cohort_col].isin(["W1", "W3", "W4"]), "New_Cohort"] = "GS"

        # Color assignment
        colors = {"GS" : "#47C3F4"}
        remaining_colors = [ "#FF9429", "#63C7B2", "#FFCE00", "#EF553C", "#2C819B", "#07B1D8", "#C5E3EA", "#326B64", "#FFD15F", "#EA4C15"]
        geo_lbc = sorted([i for i in set(df[cohort_col]) if not i.startswith("W")], reverse = True)
        i = 0
        for gset in geo_lbc:
            colors[gset] = remaining_colors[i]
            i += 1

        # Density plot
        fig, axes = plt.subplots(1, 2, figsize = (10, 3), gridspec_kw=dict(width_ratios=[1.6,1]))
        for cohort in ["GS"] + geo_lbc:
            sns.kdeplot(ax = axes[0], data = df[df["New_Cohort"] == cohort], x = age_col, shade = True, label = cohort, color = colors[cohort], linewidth=2)
        axes[0].set_xlabel("Age")
        axes[0].set_ylabel("Density")

        # Box plot
        sns.boxplot(data = df, y = "New_Cohort", x = age_col, ax = axes[1], palette = colors, linewidth = 0.8, saturation = 1)
        #sns.stripplot(data = df, y = "New_Cohort", x = age_col, ax = axes[1], palette = colors, alpha = 0.3)
        axes[1].yaxis.tick_right()
        axes[1].set_xlabel("Age")
        axes[1].set(ylabel=None)
        
        # Final bits
        plt.suptitle("Age Distribution")
        #elements = [mpatches.Patch(edgecolor = colors[i], facecolor = colors[i], label = i) for i in (["GS"] + geo_lbc)]
        #fig.legend(handles=elements, loc = 3, ncol = 4, fancybox = False, frameon=False, prop={"size":9})
        plt.tight_layout()
        #fig.subplots_adjust(bottom=0.35) 
        plt.savefig(output_dir + "agedistribution_%s.pdf" % (name))
        plt.close(fig)
    
    else: 
        # Color assignment
        colors = {}
        remaining_colors = ["#47C3F4", "#FF9429", "#63C7B2"]
        new_names = {"W1" : "Set 1", "W3" : "Set 2", "W4" : "Set 3"}
        df = df.sort_values(by = "Set")
        for wave in set(df["Set"]):
            df.loc[df["Set"] == wave, "Set"] = new_names[wave]
        waves = sorted([i for i in set(df["Set"])])
        i = 0
        for wave in waves:
            colors[wave] = remaining_colors[i]
            i += 1

        # Density plot
        fig, axes = plt.subplots(1, 2, figsize = (10, 3), gridspec_kw=dict(width_ratios=[1.6,1]))
        for cohort in waves:
            sns.kdeplot(ax = axes[0], data = df[df["Set"] == cohort], x = age_col, shade = True, label = cohort, color = colors[cohort], linewidth=2)
        axes[0].set_xlabel("Age")
        axes[0].set_ylabel("Density")

        # Box plot
        sns.boxplot(data = df, y = "Set", x = age_col, ax = axes[1], palette = colors, linewidth = 1, saturation = 1)
        #sns.stripplot(data = df, y = "Set", x = age_col, ax = axes[1], palette = colors, alpha = 0.3)
        axes[1].yaxis.tick_right()
        axes[1].set_xlabel("Age")
        axes[1].set(ylabel=None)
        
        # Final bits
        plt.suptitle("Age Distribution")
        #elements = [mpatches.Patch(edgecolor = colors[i], facecolor = colors[i], label = i) for i in waves]
        #fig.legend(handles=elements, loc = 3, ncol = 4, fancybox = False, frameon=False, prop={"size":9})
        plt.tight_layout()
        #fig.subplots_adjust(bottom=0.3) 
        plt.savefig(output_dir + "agedistribution_%s.pdf" % (name))
        plt.close(fig)

# Sex distribution
################################################################

# target = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/chronological_age/elasticnet_models/cv_folds/w1w3w4/random/gs_basic_folds_random_withexternal.tsv"
# output_dir = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/chronological_age/data_explore/w1w3w4/"
# cohort_col = "cohort"
# age_col = "age"
# sex_col = "sex"
# external = True
# name = "noscale_external"

def sex_dist(target, output_dir, cohort_col, sex_col, external, name):
    palette = set_colors(2, "custom2")
    sns.set_palette(sns.color_palette(palette))
    
    if external == True:
        df = pd.read_table(target, index_col = 0)
        df.loc[df[cohort_col].isin(["W1", "W3", "W4"]), cohort_col] = "GS"
    else:
        df = pd.read_table(target, index_col = 1)
        cohort_col = "Set"
        new_names = {"W1" : "Set 1", "W3" : "Set 2", "W4" : "Set 3"}
        df = df.sort_values(by = cohort_col)
        for wave in set(df[cohort_col]):
            df.loc[df["Set"] == wave, cohort_col] = new_names[wave]
    df = df[~df[sex_col].isnull()].copy()
    df = df.rename(columns = {"sex":"Sex"})
    sex_col = "Sex"

    count = df.groupby(cohort_col)[sex_col].value_counts()
    count = count.rename('Count').reset_index()
    count['Percentage'] = count['Count'].div(count.groupby(cohort_col)['Count'].transform(lambda x: x.sum()))
    
    if external == True:
        fig, axes = plt.subplots(2, 1, figsize = (6, 4.5), sharex=True)
        sns.barplot(ax = axes[0], x = cohort_col, y = "Count", data = count, hue = sex_col, saturation = 1, linewidth = 0)
        sns.barplot(ax = axes[1], x = cohort_col, y = "Percentage", data = count, hue = sex_col, saturation = 1, linewidth = 0)
        axes[0].set_xlabel("")
        axes[1].set_xlabel("")
        axes[1].tick_params(axis='x', labelrotation=45)
        axes[1].get_legend().remove()
        plt.suptitle("Sex Distribution")
        plt.tight_layout()
        plt.savefig(output_dir + "sexdistribution_%s.pdf" % (name))
        plt.close(fig)
    else:
        fig, axes = plt.subplots(1, 2, figsize = (5.5, 2.5))
        sns.barplot(ax = axes[0], x = cohort_col, y = "Count", data = count, hue = sex_col, saturation = 1, linewidth = 0)
        sns.barplot(ax = axes[1], x = cohort_col, y = "Percentage", data = count, hue = sex_col, saturation = 1, linewidth = 0)
        axes[0].set_xlabel("")
        axes[1].set_xlabel("")
        #axes[1].tick_params(axis='x', labelrotation=45)
        axes[1].get_legend().remove()
        axes[0].legend(loc = 2)
        plt.suptitle("Sex Distribution")
        plt.tight_layout()
        plt.savefig(output_dir + "sexdistribution_%s.pdf" % (name))
        plt.close(fig)


# CpG across ages (ELOVL2)
################################################################

def cpg_trajectory_external(prepped_file, cpg, output_dir, name, cpg_col = "CpG", age_col = "Age", sex_col = "Sex", cohort_col = "Cohort"):
    # Import data
    df = pd.read_table(prepped_file, index_col = 0)
    
    # Fuse GS and LBC observations
    df["New_Cohort"] = df[cohort_col]
    df.loc[df[cohort_col].isin(["W1", "W3", "W4"]), "New_Cohort"] = "GS"
    df.loc[df[cohort_col].isin(["LBC21", "LBC36"]), "New_Cohort"] = "LBC"
    
    # Color assignment
    colors = {"GS" : "#47C3F4", "LBC": "#FF9429"}
    remaining_colors = ["#63C7B2", "#FFCE00", "#EF553C", "#2C819B", "#07B1D8", "#C5E3EA", "#326B64", "#FFD15F", "#EA4C15"]
    geo = [i for i in set(df[cohort_col]) if i.startswith("GSE")]
    i = 0
    for gset in geo:
        colors[gset] = remaining_colors[i]
        i += 1
    
    # Plot - all together
    fig, axes = plt.subplots(1, 2, figsize = (8, 3.5), sharey=True, gridspec_kw=dict(width_ratios=[4,1]))
    for cohort in (["GS", "LBC"] + geo):
        sns.regplot(ax = axes[0], data = df[df["New_Cohort"] == cohort], x = age_col, y = cpg_col, label = cohort, color = colors[cohort], marker = '.', scatter_kws={'alpha':1, 'rasterized':True}, line_kws={'color':lighten_color(colors[cohort])}, lowess = True)
        sns.kdeplot(ax = axes[1], data = df[df["New_Cohort"] == cohort], y = cpg_col, shade = True, label = cohort, color = colors[cohort], linewidth = 2)
    axes[0].set_ylabel("CpG Beta Value")
    plt.suptitle(cpg)
    elements = [mpatches.Patch(edgecolor = colors[i], facecolor = colors[i], label = i) for i in (["GS", "LBC"] + geo)]
    fig.legend(handles=elements, loc = 3, ncol = 5, fancybox = False, frameon=False)
    plt.tight_layout()
    fig.subplots_adjust(bottom=0.3)  
    plt.savefig(output_dir + "%s_trajectory%s.pdf" % (cpg, name), dpi = 300)
    plt.close(fig)

    # Plot - separately
    rc = [(i, j) for i in range(0,2) for j in range(0,5)]
    fig, axes = plt.subplots(2, 5, figsize = (15, 6))
    i = 0
    for cohort in (["GS", "LBC"] + geo):
        subset = df[df["New_Cohort"] == cohort]
        r = scipy.stats.pearsonr(subset[cpg_col], subset[age_col])
        sns.regplot(ax = axes[rc[i][0], rc[i][1]], data = subset, x = age_col, y = cpg_col, label = "r = %s (pv = %s)" % (round(r[0], 3), round(r[1], 3)), color = colors[cohort], marker = '.', scatter_kws={'alpha':1, 'rasterized':True}, line_kws={'color':lighten_color(colors[cohort])}, lowess = True)
        axes[rc[i][0], rc[i][1]].legend()
        axes[rc[i][0], rc[i][1]].set_xlim([0, 100])
        axes[rc[i][0], rc[i][1]].set_ylabel("CpG Beta Value")
        axes[rc[i][0], rc[i][1]].set_title(cohort)
        i += 1
    plt.suptitle(cpg)
    #elements = [mpatches.Patch(edgecolor = colors[i], facecolor = colors[i], label = i) for i in (["GS", "LBC"] + geo)]
    #fig.legend(handles=elements, loc = 3, ncol = 5, fancybox = False)
    plt.tight_layout()
    #fig.subplots_adjust(bottom=0.2)  
    plt.savefig(output_dir + "%s_trajectory%s_perdataset.pdf" % (cpg, name), dpi = 300)
    plt.close(fig)

    # Plot - separately (same axes)
    rc = [(i, j) for i in range(0,2) for j in range(0,5)]
    fig, axes = plt.subplots(2, 5, figsize = (15, 6))
    i = 0
    for cohort in (["GS", "LBC"] + geo):
        subset = df[df["New_Cohort"] == cohort]
        r = scipy.stats.pearsonr(subset[cpg_col], subset[age_col])
        sns.regplot(ax = axes[rc[i][0], rc[i][1]], data = subset, x = age_col, y = cpg_col, label = "r = %s (pv = %s)" % (round(r[0], 3), round(r[1], 3)), color = colors[cohort], marker = '.', scatter_kws={'alpha':1, 'rasterized':True}, line_kws={'color':lighten_color(colors[cohort])}, lowess = True)
        axes[rc[i][0], rc[i][1]].legend()
        axes[rc[i][0], rc[i][1]].set_xlim([0, 100])
        axes[rc[i][0], rc[i][1]].set_ylim([0, 1])
        axes[rc[i][0], rc[i][1]].set_ylabel("CpG Beta Value")
        axes[rc[i][0], rc[i][1]].set_title(cohort)
        i += 1
    plt.suptitle(cpg)
    #elements = [mpatches.Patch(edgecolor = colors[i], facecolor = colors[i], label = i) for i in (["GS", "LBC"] + geo)]
    #fig.legend(handles=elements, loc = 3, ncol = 5, fancybox = False)
    plt.tight_layout()
    #fig.subplots_adjust(bottom=0.2)  
    plt.savefig(output_dir + "%s_trajectory%s_perdataset_sameaxes.pdf" % (cpg, name), dpi = 300)
    plt.close(fig)

def cpg_trajectory(prepped_file, cpg, output_dir, name, cpg_col = "CpG", age_col = "Age", sex_col = "Sex", cohort_col = "Cohort"):
    # Import data
    df = pd.read_table(prepped_file, index_col = 0)
    
    # Color assignment
    new_names = {"W1" : "Set 1", "W3" : "Set 2", "W4" : "Set 3"}
    df = df.sort_values(by = "Cohort")
    for wave in set(df["Cohort"]):
        df.loc[df["Cohort"] == wave, "Cohort"] = new_names[wave]
    colors = {"Set 1" : "#47C3F4", "Set 2": "#FF9429", "Set 3": "#63C7B2"}

    # Plot - all together
    fig, axes = plt.subplots(1, 2, figsize = (8, 3.5), sharey=True, gridspec_kw=dict(width_ratios=[4,1]))
    for cohort in sorted(list(set(df["Cohort"]))):
        sns.regplot(ax = axes[0], data = df[df[cohort_col] == cohort], x = age_col, y = cpg_col, label = cohort, color = colors[cohort], marker = '.', scatter_kws={'alpha':1, 'rasterized':True}, line_kws={'color':lighten_color(colors[cohort])}, lowess = True)
        sns.kdeplot(ax = axes[1], data = df[df[cohort_col] == cohort], y = cpg_col, shade = True, label = cohort, color = colors[cohort], linewidth = 2)
    axes[0].set_ylabel("CpG Beta Value")
    plt.suptitle(cpg)
    elements = [mpatches.Patch(edgecolor = colors[i], facecolor = colors[i], label = i) for i in sorted(list(set(df["Cohort"])))]
    fig.legend(handles=elements, loc = 3, ncol = 5, fancybox = False, frameon=False)
    plt.tight_layout()
    fig.subplots_adjust(bottom=0.2)  
    plt.savefig(output_dir + "%s_trajectory%s.pdf" % (cpg, name), dpi = 300)
    plt.close(fig)

    # Plot - separately
    fig, axes = plt.subplots(1, 3, figsize = (9, 3.25))
    i = 0
    for cohort in sorted(list(set(df["Cohort"]))):
        subset = df[df[cohort_col] == cohort]
        r = scipy.stats.pearsonr(subset[cpg_col], subset[age_col])
        sns.regplot(ax = axes[i], data = subset, x = age_col, y = cpg_col, label = "r = %s (pv = %s)" % (round(r[0], 3), round(r[1], 3)), color = colors[cohort], marker = '.', scatter_kws={'alpha':1, 'rasterized':True}, line_kws={'color':lighten_color(colors[cohort])}, lowess = True)
        axes[i].legend()
        axes[i].set_xlim([0, 100])
        axes[i].set_ylabel("CpG Beta Value")
        axes[i].set_title(cohort)
        i += 1
    plt.suptitle(cpg)
    #elements = [mpatches.Patch(edgecolor = colors[i], facecolor = colors[i], label = i) for i in (["GS", "LBC"] + geo)]
    #fig.legend(handles=elements, loc = 3, ncol = 5, fancybox = False)
    plt.tight_layout()
    #fig.subplots_adjust(bottom=0.2)  
    plt.savefig(output_dir + "%s_trajectory%s_perdataset.pdf" % (cpg, name), dpi = 300)
    plt.close(fig)

    # Plot - separately (same axes)
    fig, axes = plt.subplots(1, 3, figsize = (9, 3.25))
    i = 0
    for cohort in sorted(list(set(df["Cohort"]))):
        subset = df[df[cohort_col] == cohort]
        r = scipy.stats.pearsonr(subset[cpg_col], subset[age_col])
        sns.regplot(ax = axes[i], data = subset, x = age_col, y = cpg_col, label = "r = %s (pv = %s)" % (round(r[0], 3), round(r[1], 3)), color = colors[cohort], marker = '.', scatter_kws={'alpha':1, 'rasterized':True}, line_kws={'color':lighten_color(colors[cohort])}, lowess = True)
        axes[i].legend()
        axes[i].set_xlim([0, 100])
        axes[i].set_ylim([0, 1])
        axes[i].set_ylabel("CpG Beta Value")
        axes[i].set_title(cohort)
        i += 1
    plt.suptitle(cpg)
    #elements = [mpatches.Patch(edgecolor = colors[i], facecolor = colors[i], label = i) for i in (["GS", "LBC"] + geo)]
    #fig.legend(handles=elements, loc = 3, ncol = 5, fancybox = False)
    plt.tight_layout()
    #fig.subplots_adjust(bottom=0.2)  
    plt.savefig(output_dir + "%s_trajectory%s_perdataset_sameaxes.pdf" % (cpg, name), dpi = 300)
    plt.close(fig)

# PCA
################################################################

def pca_plot(pca_df_f, expvar_df_f, target_f, output_dir, name, external):
    # Import
    pca_df = pd.read_table(pca_df_f, index_col = 0)
    expvar_df = pd.read_table(expvar_df_f, index_col = 0)
    if external == True:
        target = pd.read_table(target_f, index_col = 0)
        group_col = "cohort"
        target.loc[target[group_col].isin(["W1", "W3", "W4"]), group_col] = "GS"
    else:
        target = pd.read_table(target_f, index_col = 1)
        group_col = "Set"
        new_names = {"W1" : "Set 1", "W3" : "Set 2", "W4" : "Set 3"}
        for wave in set(target[group_col]):
            target.loc[target[group_col] == wave, group_col] = new_names[wave]
    pca_df = pd.concat([pca_df, target], axis = 1)
    # Colors
    if external == True:
        colors = {"GS" : "#63C7B2"}
        remaining_colors = ["#326B64", "#FF9429", "#FFCE00", "#EF553C", "#47C3F4", "#07B1D8", "#C5E3EA", "#2C819B", "#FFD15F", "#EA4C15"]
        geo_lbc = [i for i in set(pca_df[group_col]) if i != "GS"]
        i = 0
        for gset in geo_lbc:
            colors[gset] = remaining_colors[i]
            i += 1
    else:
        colors = {"Set 1" : "#47C3F4", "Set 2": "#FF9429", "Set 3": "#63C7B2"}
    # Plot
    fig, axes = plt.subplots(2, 3, figsize = (10, 6.5))
    # PC1 vs PC2
    sns.scatterplot(ax = axes[0,0], data = pca_df, x = "PC1", y = "PC2", hue = group_col, palette = colors, rasterized = True, legend = False)
    sns.kdeplot(ax = axes[1,0], data = pca_df, x = "PC1", y = "PC2", hue = group_col, palette = colors, legend = False)
    axes[0,0].set_ylabel("PC2 (%s%%)" % round(expvar_df.at["PC2", "ExpVar"]*100, 2))
    axes[0,0].set_xlabel("PC1 (%s%%)" % round(expvar_df.at["PC1", "ExpVar"]*100, 2))
    axes[0,0].set_title("PC1 vs PC2")
    axes[1,0].set_ylabel("PC2 (%s%%)" % round(expvar_df.at["PC2", "ExpVar"]*100, 2))
    axes[1,0].set_xlabel("PC1 (%s%%)" % round(expvar_df.at["PC1", "ExpVar"]*100, 2))
    # PC1 vs PC3
    sns.scatterplot(ax = axes[0,1], data = pca_df, x = "PC1", y = "PC3", hue = group_col, palette = colors, rasterized = True, legend = False)
    sns.kdeplot(ax = axes[1,1], data = pca_df, x = "PC1", y = "PC3", hue = group_col, palette = colors, legend = False)
    axes[0,1].set_ylabel("PC3 (%s%%)" % round(expvar_df.at["PC3", "ExpVar"]*100, 2))
    axes[0,1].set_xlabel("PC1 (%s%%)" % round(expvar_df.at["PC1", "ExpVar"]*100, 2))
    axes[0,1].set_title("PC1 vs PC3")
    axes[1,1].set_ylabel("PC3 (%s%%)" % round(expvar_df.at["PC3", "ExpVar"]*100, 2))
    axes[1,1].set_xlabel("PC1 (%s%%)" % round(expvar_df.at["PC1", "ExpVar"]*100, 2))
    # PC2 vs PC3
    sns.scatterplot(ax = axes[0,2], data = pca_df, x = "PC2", y = "PC3", hue = group_col, palette = colors, rasterized = True, legend = False)
    sns.kdeplot(ax = axes[1,2], data = pca_df, x = "PC2", y = "PC3", hue = group_col, palette = colors, legend = False)
    axes[0,2].set_ylabel("PC3 (%s%%)" % round(expvar_df.at["PC3", "ExpVar"]*100, 2))
    axes[0,2].set_xlabel("PC2 (%s%%)" % round(expvar_df.at["PC2", "ExpVar"]*100, 2))
    axes[0,2].set_title("PC2 vs PC3")
    axes[1,2].set_ylabel("PC3 (%s%%)" % round(expvar_df.at["PC3", "ExpVar"]*100, 2))
    axes[1,2].set_xlabel("PC2 (%s%%)" % round(expvar_df.at["PC2", "ExpVar"]*100, 2))
    # Add a global legend
    if external == True:
        elements = [mpatches.Patch(edgecolor = colors[i], facecolor = colors[i], label = i) for i in (["GS"] + sorted(geo_lbc))]
    else: 
        elements = [mpatches.Patch(edgecolor = colors[i], facecolor = colors[i], label = i) for i in ["Set 1", "Set 2", "Set 3"]]
    fig.legend(handles=elements, loc = 3, ncol = 6, fancybox = False, frameon=False)
    plt.tight_layout()
    fig.subplots_adjust(bottom=0.2)  
    # Export
    plt.savefig(output_dir + name + "_pca.pdf", dpi = 300)
    plt.close(fig)

# prepped_file = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/chronological_age/data_explore/w1w3w4/methbetavals_noscale_withexternal.csv"
# output_dir = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/chronological_age/data_explore/w1w3w4/"
# target = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/chronological_age/elasticnet_models/cv_folds/w1w3w4/random/gs_basic_folds_random_withexternal.tsv"
# name = "random_noscale_external"
# external = True
# plot = True

def pca(prepped_file, output_dir, target, name, external):
    # Import
    if external == True:
        target = pd.read_table(target, index_col = 0)
    else:
        target = pd.read_table(target, index_col = 1)
    df = dt.fread(prepped_file, sep = ",")
    df = df.to_pandas()
    df = df.set_index("ID")
    df = df[df.index.isin(target.index)]
    names = df.index.values.tolist()
    
    # Scale
    df = StandardScaler().fit_transform(df)
    
    # PCA
    pca = PCA(n_components=20)
    pcs = pca.fit_transform(df)
    
    # Result
    pca_df = pd.DataFrame(data = pcs, columns = ["PC%s" % i for i in range(1,21)], index = names)
    expvar_df = pd.DataFrame(data = pca.explained_variance_ratio_, index = ["PC%s" % i for i in range(1,21)], columns = ["ExpVar"])

    # Export
    pca_df.to_csv(output_dir + name + "_PCA.tsv", sep = "\t", index_label = "ID")
    expvar_df.to_csv(output_dir + name + "_PCA_explainedvariance.tsv", sep = "\t", index_label = "PC")