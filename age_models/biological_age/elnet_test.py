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
from chronological_age.elnet_folds import *
from lifelines import KaplanMeierFitter


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

def explore_test(prediction_file, bage_pred_col, grim_pred_col, true_col, dataset_col, output_dir, palette, tte = False, death_col = "Dead", name = False):
	df = pd.read_table(prediction_file, index_col = 0)
	df = df[df["TTE"].notna()]

	if tte == True:
		df = df[df[death_col] == 1]
		df = df[df["TTE"].notna()]

	# Correlation total
	r_bage = scipy.stats.pearsonr(df[bage_pred_col], df[true_col])
	r_grim = scipy.stats.pearsonr(df[grim_pred_col], df[true_col])
	# Residual mean square error total
	res_mse_bage = rmse(df[bage_pred_col], df[true_col])
	res_mse_grim = rmse(df[grim_pred_col], df[true_col])

	# Plot
	col = set_colors(8, palette)
	sns.set_palette(sns.color_palette(col))
	#x = range(-10, 110)
	#y = range(-10, 110)

	fig, axes = plt.subplots(1, 2, figsize = (11, 5.5))
	sns.scatterplot(ax = axes[0], data = df, x = bage_pred_col, y = true_col, hue = dataset_col)
	axes[0].set_xlabel("Predicted Age")
	axes[0].set_title("bAge\nr = %s (pv = %s)\nRMSE = %s" % (round(r_bage[0], 3), round(r_bage[1], 3), round(res_mse_bage, 3)))
	sns.scatterplot(ax = axes[1], data = df, x = grim_pred_col, y = true_col, hue = dataset_col)
	axes[1].set_xlabel("Predicted Age")
	axes[1].set_title("GrimAge\nr = %s (pv = %s)\nRMSE = %s" % (round(r_grim[0], 3), round(r_grim[1], 3), round(res_mse_grim, 3)))
	#plt.plot(x, y, linewidth=1.5, color = "grey")
	if tte == True:
		axes[0].set_ylabel("Time to Event")
		axes[1].set_ylabel("Time to Event")
	else:
		axes[0].set_ylabel("True Age")
		axes[1].set_ylabel("True Age")
		
	if tte == True:
		plt.savefig(output_dir + "prediction_tte.pdf")
	else:
		plt.savefig(output_dir + "prediction.pdf")
	plt.close(fig)

	# Now make plot dividing each dataset
	fig, axes = plt.subplots(2, 2, figsize = (7, 7))
	i = 0
	used = set()
	unique = [x for x in df[dataset_col] if x not in used and (used.add(x) or True)]
	for dataset in unique:
		subset = df[df[dataset_col] == dataset]
		
		# Correlation subset
		r_bage = scipy.stats.pearsonr(subset[bage_pred_col], subset[true_col])
		r_grim = scipy.stats.pearsonr(subset[grim_pred_col], subset[true_col])
		# Residual mean square error subset
		res_mse_bage = rmse(subset[bage_pred_col], subset[true_col])
		res_mse_grim = rmse(subset[grim_pred_col], subset[true_col])

		sns.scatterplot(ax = axes[i,0], data = subset, x = bage_pred_col, y = true_col, color = col[i], label = "r = %s (pv = %s)\nRMSE = %s" % (round(r_bage[0], 3), round(r_bage[1], 3), round(res_mse_bage, 3)))
		axes[i,0].set_xlabel("Predicted Age")
		axes[i,0].set_title("bAge | %s" % dataset)
		sns.scatterplot(ax = axes[i,1], data = subset, x = grim_pred_col, y = true_col, color = col[i], label = "r = %s (pv = %s)\nRMSE = %s" % (round(r_grim[0], 3), round(r_grim[1], 3), round(res_mse_grim, 3)))
		axes[i,1].set_xlabel("Predicted Age")
		axes[i,1].set_title("GrimAge | %s" % dataset)
		#axes[i].plot(x, y, linewidth=1.5, color = "grey")
		if tte == True:
			axes[i,0].set_ylabel("Time to Event")
			axes[i,1].set_ylabel("Time to Event")
		else:
			axes[i,0].set_ylabel("True Age")
			axes[i,1].set_ylabel("True Age")
		#axes[i].set_ylim([0, 100])
		#axes[i].set_xlim([0, 100])
		i += 1

	plt.tight_layout()
	if tte == True:
		plt.savefig(output_dir + "prediction_perdataset_tte.pdf")
	else:   
		plt.savefig(output_dir + "prediction_perdataset.pdf")
	plt.close(fig)


def correlation_predictors(prediction_file, bage_pred_col, grim_pred_col, dataset_col, output_dir, palette, name = False, cpg_bage_comp = False, prediction_file_bagecpg = None):
	df = pd.read_table(prediction_file, index_col = 0)
	df = df[df["TTE"].notna()]

	# Correlation total
	if cpg_bage_comp == True:
		cpg_bage = pd.read_table(prediction_file_bagecpg, index_col = 0)
		cpg_bage = cpg_bage[cpg_bage["TTE"].notna()]
		df = df[df.index.isin(cpg_bage.index)] # 11 samples with no bAge_cpg --> No methylation data?
		bage_pred_col_cpg = bage_pred_col + "_cpg"
		df[bage_pred_col_cpg] = cpg_bage[bage_pred_col]
		print(df)
		print(df[df["bAge_cpg"].isna()])
		r = scipy.stats.pearsonr(df[bage_pred_col], df[bage_pred_col_cpg])
	else:
		r = scipy.stats.pearsonr(df[bage_pred_col], df[grim_pred_col])

	# Plot
	col = set_colors(8, palette)
	sns.set_palette(sns.color_palette(col))

	fig, axes = plt.subplots(1, 1, figsize = (5, 5))
	if cpg_bage_comp == True:
		sns.scatterplot(ax = axes, data = df, x = bage_pred_col, y = bage_pred_col_cpg, hue = dataset_col)
		axes.set_ylabel("EWAS Age")
		axes.set_xlabel("bAge")
		axes.set_title("bAge vs EWAS Age\nr = %s (pv = %s)" % (round(r[0], 3), round(r[1], 3)))
		plt.tight_layout()
		plt.savefig(output_dir + "bagevscpgbage_prediction.pdf")
	else:
		sns.scatterplot(ax = axes, data = df, x = bage_pred_col, y = grim_pred_col, hue = dataset_col)
		axes.set_ylabel("GrimAge")
		axes.set_xlabel("bAge")
		axes.set_title("bAge vs GrimAge\nr = %s (pv = %s)" % (round(r[0], 3), round(r[1], 3)))
		plt.tight_layout()
		plt.savefig(output_dir + "bagevsgrimage_prediction.pdf")
	plt.close(fig)

	# Now make plot dividing each dataset
	fig, axes = plt.subplots(1, 2, figsize = (8, 4))
	used = set()
	unique = [x for x in df[dataset_col] if x not in used and (used.add(x) or True)]
	i = 0
	for dataset in unique:
		subset = df[df[dataset_col] == dataset]
		# Correlation subset
		if cpg_bage_comp == True:
			r = scipy.stats.pearsonr(subset[bage_pred_col], subset[bage_pred_col_cpg])
			sns.scatterplot(ax = axes[i], data = subset, x = bage_pred_col, y = bage_pred_col_cpg, color = col[i], label = "r = %s (pv = %s)" % (round(r[0], 3), round(r[1], 3)))
			axes[i].set_xlabel("bAge")
			axes[i].set_ylabel("EWAS Age")
			axes[i].set_title("bAge vs EWAS Age | %s" % dataset)
			plt.tight_layout()
			plt.savefig(output_dir + "bagevscpgbage_prediction_perdataset.pdf")
		else:
			r = scipy.stats.pearsonr(subset[bage_pred_col], subset[grim_pred_col])
			sns.scatterplot(ax = axes[i], data = subset, x = bage_pred_col, y = grim_pred_col, color = col[i], label = "r = %s (pv = %s)" % (round(r[0], 3), round(r[1], 3)))
			axes[i].set_xlabel("bAge")
			axes[i].set_ylabel("GrimAge")
			axes[i].set_title("bAge vs GrimAge | %s" % dataset)
			plt.tight_layout()
			plt.savefig(output_dir + "bagevsgrimage_prediction_perdataset.pdf")
		i += 1 
	plt.close(fig)


def forest_plot(cox_file, predictions_file, output_dir, palette, hr_col = "HR", ci_low_col = "HR_CI95_Low", ci_high_col = "HR_CI95_High", grim_row = "GrimAge_GrimAgeAccel", bage_row = "bAge_bAgeAccel", cage_row = "bAge_Age", coh = False, cohorts = ["LBC21", "LBC36"], nocage = True): # Makes forest plot comparing bAge estimator to GrimAge
	# Color
	col = set_colors(7, palette)[-1]
	
	# Data
	data = pd.read_table(cox_file, index_col = 0)
	pred = pd.read_table(predictions_file, index_col = 0)
	pred = pred[pred["TTE"].notna()]
	deaths = len(pred[pred["Dead"] == 1].index)
	samples = len(pred.index)

	# Rows and row names
	if nocage == False:
		rows = [cage_row, grim_row, bage_row]
		pred_names = ["cAge", "GrimAgeAccel", "bAgeAccel"]
	else:
		rows = [grim_row, bage_row]
		pred_names = ["GrimAgeAccel", "bAgeAccel"]
	if coh == True:
		rows = ["%s_%s" % (i,j) for i in cohorts for j in rows]
		pred_names = pred_names * len(cohorts)
	data = data.loc[rows,]
	data["Predictor"] = pred_names
	
	# Errors for errorbars
	data["ci_low"] = data[hr_col] - data[ci_low_col]
	data["ci_high"] = data[ci_high_col] - data[hr_col]

	# Labels
	labels = []
	for row in rows:
		labels.append("%s [%s, %s]" % (round(data.loc[row, hr_col], 2), round(data.loc[row, ci_low_col], 2), round(data.loc[row, ci_high_col], 2)))
	data["Labels"] = labels

	# Plot   
	if coh == False:
		if nocage == True:
			fig, axes = plt.subplots(1, 1, figsize = (6.5, 1.75))
		else:
			fig, axes = plt.subplots(1, 1, figsize = (6.5, 2))
		axes.errorbar(x = data[hr_col], y = data["Predictor"], xerr = [data["ci_low"], data["ci_high"]], fmt='o', color = col)
		axes.set_title("Mortality Predictors (N = %s, D = %s)" % (samples, deaths))
		axes.set_xlabel("Hazard Ratio (CI 95%)")
		axes.margins(0.3, 0.3)
		axes2 = axes.twinx()
		axes2.errorbar(x = data[hr_col], y = data["Labels"], xerr = [data["ci_low"], data["ci_high"]], fmt='o', color = col, alpha = 0)

		# Line at HR = 1
		axes.axvline(1, color = "grey", linestyle = "--", linewidth = 1)
		
		# Aesthetics
		axes.set_xlim([0.7, 2])
		axes2.set_xlim([0.7, 2])
		axes.spines['right'].set_visible(False)
		axes.spines['left'].set_visible(False)
		axes2.spines['right'].set_visible(False)
		axes2.spines['left'].set_visible(False)
		axes.tick_params(axis=u'both', which=u'both',length=0)
		axes2.tick_params(axis=u'both', which=u'both',length=0)
		axes2.set_yticklabels(data["Labels"])
		axes.margins(0.3, 0.3)
		axes2.margins(0.3, 0.3)     
		
		# Export
		plt.tight_layout()
		plt.savefig(output_dir + "forest_plot.pdf")
		plt.close(fig)
   
	else:
		if nocage == True:
			fig, axes = plt.subplots(len(cohorts), 1, figsize = (6.5, 1.75*len(cohorts)))
		else:
			fig, axes = plt.subplots(len(cohorts), 1, figsize = (6.5, 2*len(cohorts)))
		for i in range(0, len(cohorts)):
			data_c = data.loc[[j for j in data.index if j.startswith(cohorts[i])],]
			pred_c = pred[pred["Cohort"] == cohorts[i]]
			deaths_c = len(pred_c[pred_c["Dead"] == 1].index)
			samples_c = len(pred_c.index)

			axes[i].errorbar(x = data_c[hr_col], y = data_c["Predictor"], xerr = [data_c["ci_low"], data_c["ci_high"]], fmt='o', color = col)
			axes[i].set_title("%s (N = %s, D = %s)" % (cohorts[i], samples_c, deaths_c))
			axes[i].set_xlabel("Hazard Ratio (CI 95%)")
			axes[i].margins(0.3, 0.3)
			axes2 = axes[i].twinx()
			axes2.errorbar(x = data_c[hr_col], y = data_c["Labels"], xerr = [data_c["ci_low"], data_c["ci_high"]], fmt='o', color = col, alpha = 0)

			# Line at HR = 1
			axes[i].axvline(1, color = "grey", linestyle = "--", linewidth = 1)
			
			# Aesthetics
			axes[i].set_xlim([0.7, 2])
			axes2.set_xlim([0.7, 2])
			axes[i].spines['right'].set_visible(False)
			axes[i].spines['left'].set_visible(False)
			axes2.spines['right'].set_visible(False)
			axes2.spines['left'].set_visible(False)
			axes[i].tick_params(axis=u'both', which=u'both',length=0)
			axes2.tick_params(axis=u'both', which=u'both',length=0)
			axes2.set_yticklabels(data_c["Labels"])
			axes[i].margins(0.3, 0.3)
			axes2.margins(0.3, 0.3)     
		
		# Export
		plt.suptitle("Mortality Predictors")
		plt.tight_layout()
		plt.subplots_adjust(top=0.85)
		plt.savefig(output_dir + "forest_plot_perdataset.pdf")
		plt.close(fig)       


def kaplanmeier_plot(prediction_file, grim_col, bage_col, dataset_col, output_dir, palette, tte_col = "TTE", death_col = "Dead", coh = False):
	# Import df
	df = pd.read_table(prediction_file, index_col = 0)
	df = df[df[tte_col].notna()]
	
	# Datasets
	ds = set(df[dataset_col])
	
	# Columns of interest
	ci_col = "%s_lohi"
	ci_col_ds = "%s_lohi_%s"

	# Colors
	col = set_colors(2, palette)
	col_low = col[0]
	col_high = col[1]

	if coh == False: # Not stratified by cohort
		########## Plot
		fig, axes = plt.subplots(1, 2, figsize = (7, 3.5))

		# First plot is bAge
		kmf1 = KaplanMeierFitter()
		groups = df[ci_col % bage_col]   
		i1 = (groups == 0) # Bottom quartile
		i2 = (groups == 1) # Top quartile
		print(sum(df.loc[i1, tte_col].isna()))
		kmf1.fit(df.loc[i1, tte_col], df.loc[i1, death_col], label='Bottom quartile')
		kmf1.plot(ax = axes[0], color = col_low)
		## fit the model for 2nd cohort
		kmf1.fit(df.loc[i2, tte_col], df.loc[i2, death_col], label='Top quartile')
		kmf1.plot(ax = axes[0], color = col_high)
		axes[0].set_title("bAge")
		axes[0].set_xlabel("Years - Follow Up")
		axes[0].set_ylabel("Survival Probability")

		# Second plot is GrimAge
		kmf2 = KaplanMeierFitter()
		groups = df[ci_col % grim_col]   
		i1 = (groups == 0) # Bottom quartile
		i2 = (groups == 1) # Top quartile
		kmf2.fit(df.loc[i1, tte_col], df.loc[i1, death_col], label='Bottom quartile')
		kmf2.plot(ax = axes[1], color = col_low)
		## fit the model for 2nd cohort
		kmf2.fit(df.loc[i2, tte_col], df.loc[i2, death_col], label='Top quartile')
		kmf2.plot(ax = axes[1], color = col_high)
		axes[1].set_title("GrimAge")
		axes[1].set_xlabel("Years - Follow Up")
		
		plt.tight_layout()
		plt.savefig(output_dir + "kaplanmeier_plot.pdf")
		plt.close(fig)    

	elif coh == True:
		########## Plot: Now make plot dividing each dataset
		fig, axes = plt.subplots(2, 2, figsize = (7, 7))
		i = 0
		for dataset in ds:
			subset = df[df[dataset_col] == dataset]
			# First plot is bAge
			kmf1 = KaplanMeierFitter()
			groups = df[ci_col % bage_col]   
			i1 = (groups == 0) # Bottom quartile
			i2 = (groups == 1) # Top quartile
			kmf1.fit(df.loc[i1, tte_col], df.loc[i1, death_col], label='Bottom quartile')
			kmf1.plot(ax = axes[i, 0], color = col_low)
			## fit the model for 2nd cohort
			kmf1.fit(df.loc[i2, tte_col], df.loc[i2, death_col], label='Top quartile')
			kmf1.plot(ax = axes[i, 0], color = col_high)
			axes[i, 0].set_title("bAge | %s" % dataset)
			axes[i, 0].set_xlabel("Years - Follow Up")
			axes[i, 0].set_ylabel("Survival Probability")

			# Second plot is bAge
			kmf2 = KaplanMeierFitter()
			groups = df[ci_col % grim_col]   
			i1 = (groups == 0) # Bottom quartile
			i2 = (groups == 1) # Top quartile
			kmf2.fit(df.loc[i1, tte_col], df.loc[i1, death_col], label='Bottom quartile')
			kmf2.plot(ax = axes[i, 1], color = col_low)
			## fit the model for 2nd cohort
			kmf2.fit(df.loc[i2, tte_col], df.loc[i2, death_col], label='Top quartile')
			kmf2.plot(ax = axes[i, 1], color = col_high)
			axes[i, 1].set_title("GrimAge | %s" % dataset)
			axes[i, 1].set_xlabel("Years - Follow Up")
			i += 1

		plt.tight_layout()
		plt.savefig(output_dir + "kaplanmeier_plot_perdataset.pdf")
		plt.close(fig)    