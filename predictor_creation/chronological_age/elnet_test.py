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
from scipy.stats import gaussian_kde
from chronological_age.elnet_folds import *

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


def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

def mae_(predictions, targets):
    return(np.median(abs(predictions - targets)))

def explore_test(prediction_file, pred_col, true_col, dataset_col, sex_col, sex_stratified, squared, squared_subset, squared_subset_n, squared_cpg2, squared_subset_n_cpg2, age2cpg2, external, random, lasso, loo, lbc, logage20, output_dir, palette):
    # Naming stuff
    squared_dic = {"T" : "_squared", "F" : ""}
    totallyrandom_dic = {"T" : "_randomizedfolds", "F" : ""}
    lasso_dic = {"T" : "_lasso", "F" : ""}
    loo_dic = {True : "_loo", False : ""}
    squared_subset_dic = {"T" : "_squaredsubset%s" % squared_subset_n, "F" : ""}
    if age2cpg2 == "T":
        squared_subset_dic = {"T" : "_squaredsubset%s_squaredsubsetcpg2lm%sK" % (squared_subset_n, squared_subset_n_cpg2), "F" : ""}
    lbc_dic = {"T" : "_lbc", "F" : ""}
    cpg2_dic = {"T" : "K_cpg2lm", "F" : ""}
    if age2cpg2 == "T":
        cpg2_dic = {"T" : "", "F" : ""}
    logage20_dic = {"T" : "_logage20", "F" : ""}
    
    # Import
    df = pd.read_table(prediction_file, index_col = 0)
    df = df.dropna(subset = ["Age"])
    df = df.rename(columns={dataset_col : "Set"})
    dataset_col = "Set"
    df = df.sort_values(by = dataset_col)

    # Correlation total
    r = scipy.stats.pearsonr(df[pred_col], df[true_col])
    # Root mean square error total
    res_mse = rmse(df[pred_col], df[true_col])
    # Median absolute error
    mae = mae_(df[pred_col], df[true_col])

    # Plot
    ##############################################

    col = set_colors(len(set(df[dataset_col])), palette)
    sns.set_palette(sns.color_palette(col))
    x = range(-10, 110)
    y = range(-10, 110)

    fig, axes = plt.subplots(1, 1, figsize = (5.5, 5.5))
    sns.scatterplot(data = df, x = pred_col, y = true_col, hue = dataset_col)
    plt.plot(x, y, linewidth=1.5, color = "grey")
    axes.set_ylabel("True Age")
    axes.set_xlabel("Predicted Age")
    axes.set_xlim([0, 100])
    axes.set_ylim([0, 100])
    axes.set_title("r = %s (pv = %s)\nRMSE = %s | MAE = %s" % (round(r[0], 3), round(r[1], 3), round(res_mse, 3), round(mae, 3)))
    if sex_stratified == "T":
        plt.savefig(output_dir + "prediction_sexstratified%s%s%s%s%s%s%s%s.pdf" % (lasso_dic[lasso], squared_dic[squared], squared_subset_dic[squared_subset], cpg2_dic[squared_cpg2], totallyrandom_dic[random], loo_dic[loo], lbc_dic[lbc], logage20_dic[logage20]))
    else:
        plt.savefig(output_dir + "prediction%s%s%s%s%s%s%s%s.pdf" % (lasso_dic[lasso], squared_dic[squared], squared_subset_dic[squared_subset], cpg2_dic[squared_cpg2], totallyrandom_dic[random], loo_dic[loo], lbc_dic[lbc], logage20_dic[logage20]))
    plt.close(fig)

    # Plot, separate sexes
    ##############################################
    
    # Correlation total
    r_F = scipy.stats.pearsonr(df.loc[df[sex_col] == "F", pred_col], df.loc[df[sex_col] == "F", true_col])
    r_M = scipy.stats.pearsonr(df.loc[df[sex_col] == "M", pred_col], df.loc[df[sex_col] == "M", true_col])
    # Root mean square error total
    res_mse_F = rmse(df.loc[df[sex_col] == "F", pred_col], df.loc[df[sex_col] == "F", true_col])
    res_mse_M = rmse(df.loc[df[sex_col] == "M", pred_col], df.loc[df[sex_col] == "M", true_col])
    # Median absolute error
    mae_F = mae_(df.loc[df[sex_col] == "F", pred_col], df.loc[df[sex_col] == "F", true_col])
    mae_M = mae_(df.loc[df[sex_col] == "M", pred_col], df.loc[df[sex_col] == "M", true_col])

    fig, axes = plt.subplots(1, 2, figsize = (11, 5.5))
    # Females
    sns.scatterplot(ax = axes[0], data = df[df["Sex"] == "F"], x = pred_col, y = true_col, hue = dataset_col)
    axes[0].plot(x, y, linewidth=1.5, color = "grey")
    axes[0].set_ylabel("True Age")
    axes[0].set_xlabel("Predicted Age")
    axes[0].set_xlim([0, 100])
    axes[0].set_ylim([0, 100])
    axes[0].set_title("Females\nr = %s (pv = %s)\nRMSE = %s | MAE = %s" % (round(r_F[0], 3), round(r_F[1], 3), round(res_mse_F, 3), round(mae_F, 3)))
    # Males
    sns.scatterplot(ax = axes[1], data = df[df["Sex"] == "M"], x = pred_col, y = true_col, hue = dataset_col)
    axes[1].plot(x, y, linewidth=1.5, color = "grey")
    axes[1].set_ylabel("True Age")
    axes[1].set_xlabel("Predicted Age")
    axes[1].set_xlim([0, 100])
    axes[1].set_ylim([0, 100])
    axes[1].set_title("Males\nr = %s (pv = %s)\nRMSE = %s | MAE = %s" % (round(r_M[0], 3), round(r_M[1], 3), round(res_mse_M, 3), round(mae_M, 3)))
    if sex_stratified == "T":
        plt.savefig(output_dir + "prediction_persex_sexstratified%s%s%s%s%s%s%s%s.pdf" % (lasso_dic[lasso], squared_dic[squared], squared_subset_dic[squared_subset], cpg2_dic[squared_cpg2], totallyrandom_dic[random], loo_dic[loo], lbc_dic[lbc], logage20_dic[logage20]))
    else:
        plt.savefig(output_dir + "prediction_persex%s%s%s%s%s%s%s%s.pdf" % (lasso_dic[lasso], squared_dic[squared], squared_subset_dic[squared_subset], cpg2_dic[squared_cpg2], totallyrandom_dic[random], loo_dic[loo], lbc_dic[lbc], logage20_dic[logage20]))
    plt.close(fig)

    # Now make plot dividing each dataset
    ##############################################

    if external == "F":
        if lbc == "T":
            fig, axes = plt.subplots(2, 5, figsize = (16, 6.5))
            rc = [(i, j) for i in range(0,2) for j in range(0,5)]
        else:
            fig, axes = plt.subplots(2, 4, figsize = (15, 7))
            rc = [(i, j) for i in range(0,2) for j in range(0,4)]
    elif (external == "T") & (loo == False):
        fig, axes = plt.subplots(2, 3, figsize = (11.25, 7))
        rc = [(i, j) for i in range(0,2) for j in range(0,3)]
    elif (external == "T") & (loo == True):
        if lbc == "T":
            fig, axes = plt.subplots(2, 5, figsize = (16, 6.5))
            rc = [(i, j) for i in range(0,2) for j in range(0,5)]
        else:
            fig, axes = plt.subplots(2, 4, figsize = (15, 7))
            rc = [(i, j) for i in range(0,2) for j in range(0,4)]
    i = 0
    used = set()
    unique = [x for x in df[dataset_col] if x not in used and (used.add(x) or True)]
    for dataset in unique:
        subset = df[df[dataset_col] == dataset]
        
        # Correlation subset
        r = scipy.stats.pearsonr(subset[pred_col], subset[true_col])
        # Root mean square error subset
        res_mse = rmse(subset[pred_col], subset[true_col])
        # Median absolute error
        mae = mae_(subset[pred_col], subset[true_col])

        sns.scatterplot(ax = axes[rc[i][0], rc[i][1]], data = subset, x = pred_col, y = true_col, color = col[i], label = "r = %s (pv = %s)\nRMSE = %s\nMAE = %s" % (round(r[0], 3), round(r[1], 3), round(res_mse, 3), round(mae, 3)))
        axes[rc[i][0], rc[i][1]].plot(x, y, linewidth=1.5, color = "grey")
        axes[rc[i][0], rc[i][1]].set_ylabel("True Age")
        axes[rc[i][0], rc[i][1]].set_xlabel("Predicted Age")
        axes[rc[i][0], rc[i][1]].set_title(dataset)
        axes[rc[i][0], rc[i][1]].set_ylim([0, 100])
        axes[rc[i][0], rc[i][1]].set_xlim([0, 100])
        i += 1

    plt.tight_layout()
    if sex_stratified == "T": 
        plt.savefig(output_dir + "prediction_perdataset_sexstratified%s%s%s%s%s%s%s%s.pdf" % (lasso_dic[lasso], squared_dic[squared], squared_subset_dic[squared_subset], cpg2_dic[squared_cpg2], totallyrandom_dic[random], loo_dic[loo], lbc_dic[lbc], logage20_dic[logage20]))
    else:
        plt.savefig(output_dir + "prediction_perdataset%s%s%s%s%s%s%s%s.pdf" % (lasso_dic[lasso], squared_dic[squared], squared_subset_dic[squared_subset], cpg2_dic[squared_cpg2], totallyrandom_dic[random], loo_dic[loo], lbc_dic[lbc], logage20_dic[logage20]))
    plt.close(fig)


    # Now make plot dividing each dataset, density plot, and other changes
    ######################################################################

    if external == "F":
        if lbc == "T":
            fig, axes = plt.subplots(2, 5, figsize = (12, 5))
            rc = [(i, j) for i in range(0,2) for j in range(0,5)]
        else:
            fig, axes = plt.subplots(2, 4, figsize = (11, 5))
            rc = [(i, j) for i in range(0,2) for j in range(0,4)]
    elif (external == "T") & (loo == False):
        fig, axes = plt.subplots(2, 3, figsize = (10, 5))
        rc = [(i, j) for i in range(0,2) for j in range(0,3)]
    elif (external == "T") & (loo == True):
        if lbc == "T":
            fig, axes = plt.subplots(2, 5, figsize = (12, 5))
            rc = [(i, j) for i in range(0,2) for j in range(0,5)]
        else:
            fig, axes = plt.subplots(2, 4, figsize = (12, 5))
            rc = [(i, j) for i in range(0,2) for j in range(0,4)]
    i = 0
    used = set()
    unique = [x for x in df[dataset_col] if x not in used and (used.add(x) or True)]
    for dataset in unique:
        subset = df[df[dataset_col] == dataset]
        
        cmap = mc.LinearSegmentedColormap.from_list("", ["#0084AD", "#47C3F4", "#FF9429", "#FFE3CC"])
        col = set_colors(len(set(df[dataset_col])), palette)
        sns.set_palette(sns.color_palette(col))
        x = range(-10, 110)
        y = range(-10, 110)

        # Calculate the point density
        xy = np.vstack([subset[true_col], subset[pred_col]])
        z = gaussian_kde(xy)(xy)

        # Correlation subset
        r = scipy.stats.pearsonr(subset[pred_col], subset[true_col])
        # Root mean square error subset
        res_mse = rmse(subset[pred_col], subset[true_col])
        # Median absolute error
        mae = mae_(subset[pred_col], subset[true_col])

        axes[rc[i][0], rc[i][1]].scatter(y = subset[pred_col], x = subset[true_col], c = z, s = 20, cmap = cmap, label = "r = %s\nRMSE = %s\nMAE = %s" % (round(r[0], 3), round(res_mse, 3), round(mae, 3)))
        axes[rc[i][0], rc[i][1]].plot(x, y, linewidth=1, color = "black")
        if rc[i][0] == 1:
            axes[rc[i][0], rc[i][1]].set_xlabel("Chronological Age")
        else:
            axes[rc[i][0], rc[i][1]].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        if rc[i][1] == 0:
            axes[rc[i][0], rc[i][1]].set_ylabel("Epigenetic Age")
        else:
            axes[rc[i][0], rc[i][1]].tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
        axes[rc[i][0], rc[i][1]].set_title(dataset)
        axes[rc[i][0], rc[i][1]].set_ylim([0, 100])
        axes[rc[i][0], rc[i][1]].set_xlim([0, 100])
        leg = axes[rc[i][0], rc[i][1]].legend(frameon = False, handletextpad=0, handlelength=0)
        for item in leg.legendHandles:
            item.set_visible(False)
        i += 1

    plt.tight_layout()
    if sex_stratified == "T": 
        plt.savefig(output_dir + "prediction_perdataset_sexstratified%s%s%s%s%s%s%s%s_density.pdf" % (lasso_dic[lasso], squared_dic[squared], squared_subset_dic[squared_subset], cpg2_dic[squared_cpg2], totallyrandom_dic[random], loo_dic[loo], lbc_dic[lbc], logage20_dic[logage20]))
    else:
        plt.savefig(output_dir + "prediction_perdataset%s%s%s%s%s%s%s%s_density.pdf" % (lasso_dic[lasso], squared_dic[squared], squared_subset_dic[squared_subset], cpg2_dic[squared_cpg2], totallyrandom_dic[random], loo_dic[loo], lbc_dic[lbc], logage20_dic[logage20]))
    plt.close(fig)


    # Now make plot dividing each dataset, per sex
    ##############################################

    if external == "F":
        if lbc == "T":
            fig, axes = plt.subplots(2, 5, figsize = (16, 6.5))
            rc = [(i, j) for i in range(0,2) for j in range(0,5)]
        else:
            fig, axes = plt.subplots(2, 4, figsize = (15, 7))
            rc = [(i, j) for i in range(0,2) for j in range(0,4)]
    elif (external == "T") & (loo == False):
        fig, axes = plt.subplots(2, 3, figsize = (11.25, 7))
        rc = [(i, j) for i in range(0,2) for j in range(0,3)]
    elif (external == "T") & (loo == True):
        if lbc == "T":
            fig, axes = plt.subplots(2, 5, figsize = (16, 6.5))
            rc = [(i, j) for i in range(0,2) for j in range(0,5)]
        else:
            fig, axes = plt.subplots(2, 4, figsize = (15, 7))
            rc = [(i, j) for i in range(0,2) for j in range(0,4)]
    i = 0
    used = set()
    unique = [x for x in df[dataset_col] if x not in used and (used.add(x) or True)]
    for dataset in unique:
        subset = df[df[dataset_col] == dataset]
        
        # Correlation subset
        r_F = scipy.stats.pearsonr(subset.loc[subset[sex_col] == "F", pred_col], subset.loc[subset[sex_col] == "F", true_col])
        r_M = scipy.stats.pearsonr(subset.loc[subset[sex_col] == "M", pred_col], subset.loc[subset[sex_col] == "M", true_col])
        # Root mean square error subset
        res_mse_F = rmse(subset.loc[subset[sex_col] == "F", pred_col], subset.loc[subset[sex_col] == "F", true_col])
        res_mse_M = rmse(subset.loc[subset[sex_col] == "M", pred_col], subset.loc[subset[sex_col] == "M", true_col])
        # Median absolute error
        mae_F = mae_(subset.loc[subset[sex_col] == "F", pred_col], subset.loc[subset[sex_col] == "F", true_col])
        mae_M = mae_(subset.loc[subset[sex_col] == "M", pred_col], subset.loc[subset[sex_col] == "M", true_col])
        
        subset.loc[subset[sex_col] == "F", "Sex"] = "Females\nr = %s\nRMSE = %s\nMAE = %s" % (round(r_F[0], 3), round(res_mse_F, 3), round(mae_F, 3)) 
        subset.loc[subset[sex_col] == "M", "Sex"] = "Males\nr = %s\nRMSE = %s\nMAE = %s" % (round(r_M[0], 3), round(res_mse_M, 3), round(mae_M, 3)) 
        subset = subset.sort_values(by = "Sex")
        
        sns.scatterplot(ax = axes[rc[i][0], rc[i][1]], data = subset, x = pred_col, y = true_col, style = "Sex", color = col[i])
        axes[rc[i][0], rc[i][1]].plot(x, y, linewidth=1.5, color = "grey")
        axes[rc[i][0], rc[i][1]].set_ylabel("True Age")
        axes[rc[i][0], rc[i][1]].set_xlabel("Predicted Age")
        axes[rc[i][0], rc[i][1]].set_title(dataset)
        axes[rc[i][0], rc[i][1]].set_ylim([0, 100])
        axes[rc[i][0], rc[i][1]].set_xlim([0, 100])
        axes[rc[i][0], rc[i][1]].legend(fontsize = 8)
        i += 1

    plt.tight_layout()
    if sex_stratified == "T": 
        plt.savefig(output_dir + "prediction_perdataset_persex_sexstratified%s%s%s%s%s%s%s%s.pdf" % (lasso_dic[lasso], squared_dic[squared], squared_subset_dic[squared_subset], cpg2_dic[squared_cpg2], totallyrandom_dic[random], loo_dic[loo], lbc_dic[lbc], logage20_dic[logage20]))
    else:
        plt.savefig(output_dir + "prediction_perdataset_persex%s%s%s%s%s%s%s%s.pdf" % (lasso_dic[lasso], squared_dic[squared], squared_subset_dic[squared_subset], cpg2_dic[squared_cpg2], totallyrandom_dic[random], loo_dic[loo], lbc_dic[lbc], logage20_dic[logage20]))
    plt.close(fig)



def explore_resid_combo(prediction_file, pred_col, true_col, dataset_col, sex_col, sex_stratified, squared, squared_subset, squared_subset_n, squared_cpg2, squared_subset_n_cpg2, age2cpg2, external, random, lasso, loo, lbc, logage20, output_dir, palette):
    # Naming stuff
    squared_dic = {"T" : "_squared", "F" : ""}
    totallyrandom_dic = {"T" : "_randomizedfolds", "F" : ""}
    lasso_dic = {"T" : "_lasso", "F" : ""}
    loo_dic = {True : "_loo", False : ""}
    squared_subset_dic = {"T" : "_squaredsubset%s" % squared_subset_n, "F" : ""}
    lbc_dic = {"T" : "_lbc", "F" : ""}
    cpg2_dic = {"T" : "K_cpg2lm", "F" : ""}
    logage20_dic = {"T" : "_logage20", "F" : ""}
    if age2cpg2 == "T":
        squared_subset_dic = {"T" : "_squaredsubset%s_squaredsubsetcpg2lm%sK" % (squared_subset_n, squared_subset_n_cpg2), "F" : ""}
        cpg2_dic = {"T" : "", "F" : ""}
    # Import
    df = pd.read_table(prediction_file, index_col = 0)
    df = df.dropna(subset = ["Age"])
    df = df.rename(columns={dataset_col : "Set"})
    dataset_col = "Set"
    df = df.sort_values(by = dataset_col)

    # Correlation total
    r = scipy.stats.pearsonr(df[pred_col], df[true_col])
    # Root mean square error total
    res_mse = rmse(df[pred_col], df[true_col])
    # Median absolute error
    mae = mae_(df[pred_col], df[true_col])

    # Plot
    ##############################################

    col = set_colors(len(set(df[dataset_col])), palette)
    sns.set_palette(sns.color_palette(col))
    x = range(-10, 110)
    y = range(-10, 110)

    fig, axes = plt.subplots(2, 1, figsize = (5.5, 6), sharex=True, gridspec_kw=dict(height_ratios=[4,1]))
    sns.scatterplot(data = df, y = pred_col, x = true_col, hue = dataset_col, ax = axes[0])
    axes[0].plot(x, y, linewidth=1.5, color = "grey")
    #axes[0].set_xlabel("True Age")
    axes[0].set_ylabel("Predicted Age")
    axes[0].set_xlim([0, 100])
    axes[0].set_ylim([0, 100])
    axes[0].set_title("r = %s (pv = %s)\nRMSE = %s | MAE = %s" % (round(r[0], 3), round(r[1], 3), round(res_mse, 3), round(mae, 3)))
    # Add residuals on side
    df["Absolute Error"] = abs(df["Age_pred"] - df["Age"])
    color = "#07B1D8"
    sns.regplot(ax = axes[1], data = df, x = "Age", y = "Absolute Error", color = color, marker = '.', scatter_kws={'alpha':1, 'rasterized':True}, line_kws={'color':'black'}, lowess = True)
    #axes.set_title("Residual Analysis")
    axes[1].set_ylabel("Absolute Error")
    axes[1].set_xlabel("True Age")
    plt.tight_layout()

    if sex_stratified == "T":
        plt.savefig(output_dir + "prediction_sexstratified%s%s%s%s%s%s%s%s_withresiduals.pdf" % (lasso_dic[lasso], squared_dic[squared], squared_subset_dic[squared_subset], cpg2_dic[squared_cpg2], totallyrandom_dic[random], loo_dic[loo], lbc_dic[lbc], logage20_dic[logage20]), dpi = 300)
    else:
        plt.savefig(output_dir + "prediction%s%s%s%s%s%s%s%s_withresiduals.pdf" % (lasso_dic[lasso], squared_dic[squared], squared_subset_dic[squared_subset], cpg2_dic[squared_cpg2], totallyrandom_dic[random], loo_dic[loo], lbc_dic[lbc], logage20_dic[logage20]), dpi = 300)
    plt.close(fig)

    # Contour plot - all together
    ##############################################
    
    cmap = mc.LinearSegmentedColormap.from_list("", ["#0084AD", "#47C3F4", "#FF9429", "#FFE3CC"])
    col = set_colors(len(set(df[dataset_col])), palette)
    sns.set_palette(sns.color_palette(col))
    x = range(-10, 110)
    y = range(-10, 110)

    fig, axes = plt.subplots(2, 1, figsize = (5.5, 6), sharex=True, gridspec_kw=dict(height_ratios=[4,1]))
    sns.kdeplot(data = df, y = pred_col, x = true_col, ax = axes[0], cmap = cmap, fill = True)
    axes[0].plot(x, y, linewidth=1, color = "black")
    #axes[0].set_xlabel("True Age")
    axes[0].set_ylabel("Predicted Age")
    axes[0].set_xlim([0, 100])
    axes[0].set_ylim([0, 100])
    axes[0].set_title("r = %s (pv = %s)\nRMSE = %s | MAE = %s" % (round(r[0], 3), round(r[1], 3), round(res_mse, 3), round(mae, 3)))
    # Add residuals on side
    df["Absolute Error"] = abs(df["Age_pred"] - df["Age"])
    #color = "#07B1D8"
    color = "#0084AD"
    sns.regplot(ax = axes[1], data = df, x = "Age", y = "Absolute Error", color = color, marker = '.', scatter_kws={'alpha':1, 'rasterized':True}, line_kws={'color':'black'}, lowess = True)
    #axes.set_title("Residual Analysis")
    axes[1].set_ylabel("Absolute Error")
    axes[1].set_xlabel("True Age")
    plt.tight_layout()

    if sex_stratified == "T":
        plt.savefig(output_dir + "prediction_sexstratified%s%s%s%s%s%s%s%s_withresiduals_contour.pdf" % (lasso_dic[lasso], squared_dic[squared], squared_subset_dic[squared_subset], cpg2_dic[squared_cpg2], totallyrandom_dic[random], loo_dic[loo], lbc_dic[lbc], logage20_dic[logage20]), dpi = 300)
    else:
        plt.savefig(output_dir + "prediction%s%s%s%s%s%s%s%s_withresiduals_contour.pdf" % (lasso_dic[lasso], squared_dic[squared], squared_subset_dic[squared_subset], cpg2_dic[squared_cpg2], totallyrandom_dic[random], loo_dic[loo], lbc_dic[lbc], logage20_dic[logage20]), dpi = 300)
    plt.close(fig)


    # Density - all together
    ##############################################

    cmap = mc.LinearSegmentedColormap.from_list("", ["#0084AD", "#47C3F4", "#FF9429", "#FFE3CC"])
    col = set_colors(len(set(df[dataset_col])), palette)
    sns.set_palette(sns.color_palette(col))
    x = range(-10, 110)
    y = range(-10, 110)
   
    # Calculate the point density
    xy = np.vstack([df[true_col], df[pred_col]])
    z = gaussian_kde(xy)(xy)
    print(z)

    fig, axes = plt.subplots(2, 1, figsize = (5.5, 6), sharex=True, gridspec_kw=dict(height_ratios=[4,1]))
    axes[0].scatter(y = df[pred_col], x = df[true_col], c = z, s = 15, cmap = cmap)
    axes[0].plot(x, y, linewidth=1, color = "black")
    #axes[0].set_xlabel("True Age")
    axes[0].set_ylabel("Epigenetic Age")
    axes[0].set_xlabel("Chronological Age")
    axes[0].set_xlim([0, 100])
    axes[0].set_ylim([0, 100])
    axes[0].set_title("r = %s\nRMSE = %s | MAE = %s" % (round(r[0], 3), round(res_mse, 3), round(mae, 3)))
    # Add residuals on side
    df["Absolute Error"] = abs(df["Age_pred"] - df["Age"])
    #color = "#07B1D8"
    color = "#0084AD"
    sns.regplot(ax = axes[1], data = df, x = "Age", y = "Absolute Error", color = color, marker = '.', scatter_kws={'alpha':1, 'rasterized':True}, line_kws={'color':'black'}, lowess = True)
    #axes.set_title("Residual Analysis")
    axes[1].set_ylabel("Absolute Error")
    axes[1].set_xlabel("Chronological Age")
    plt.tight_layout()

    if sex_stratified == "T":
        plt.savefig(output_dir + "prediction_sexstratified%s%s%s%s%s%s%s%s_withresiduals_density.pdf" % (lasso_dic[lasso], squared_dic[squared], squared_subset_dic[squared_subset], cpg2_dic[squared_cpg2], totallyrandom_dic[random], loo_dic[loo], lbc_dic[lbc], logage20_dic[logage20]), dpi = 300)
    else:
        plt.savefig(output_dir + "prediction%s%s%s%s%s%s%s%s_withresiduals_density.pdf" % (lasso_dic[lasso], squared_dic[squared], squared_subset_dic[squared_subset], cpg2_dic[squared_cpg2], totallyrandom_dic[random], loo_dic[loo], lbc_dic[lbc], logage20_dic[logage20]), dpi = 300)
    plt.close(fig)

def residual_plot(pred_file, output_dir, name):
    pred = pd.read_table(pred_file, index_col = 0, na_values = "NA")
    pred = pred.dropna(subset = ["Age"])
    pred["Absolute Error"] = abs(pred["Age_pred"] - pred["Age"])
    color = "#07B1D8"
    fig, axes = plt.subplots(1, 1, figsize = (7, 3.5))
    sns.regplot(ax = axes, data = pred, x = "Age", y = "Absolute Error", color = color, marker = '.', scatter_kws={'alpha':1, 'rasterized':True}, line_kws={'color':'black'}, lowess = True)
    axes.set_title("Residual Analysis")
    axes.set_ylabel("Absolute Error")
    axes.set_xlabel("Age")
    plt.tight_layout()
    plt.savefig(output_dir + "%s.pdf" % name, dpi = 300)
    plt.close(fig)
    