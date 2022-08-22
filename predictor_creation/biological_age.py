#!/usr/bin/env python
# -*- coding: utf-8 -*-
# By Elena

################################################################
# Script that coordinates all data prep pre-running elastic net
# + fitting of elastic net model for bAge
################################################################

# Run on cluster
# For VSCode + cluster: export PATH=/home/ebernab3/anaconda3/bin:$PATH

import os
import numpy as np
import pandas as pd
import random
import subprocess
from chronological_age.elnet_folds import *
from biological_age.elnet_test import *
import argparse
import time
from datetime import datetime


# USER INPUT
######################################################
######################################################

parser = argparse.ArgumentParser(description='bAge analysis pipeline.')

# Toggles to run bits of code - needs at least one of these flags or pipeline won't do anything
parser.add_argument('--folds', help='Option to assign batches to folds. WARNING: Needs to be run at least once prior to training. No need if cAge analysis has already been run since same file is used!', action="store_true", default=False)
parser.add_argument('--train', help='Option to train models.', action="store_true", default=False)
parser.add_argument('--test', help='Option to test models.', action="store_true", default=False)

# Parameters - if nothing is stated will run on default
#parser.add_argument('--grimscale', help='Option to run with all of GS (waves 1, 3, and 4) instead of just waves 1 and 3 as base of models.', action="store_true", default=False)
parser.add_argument('--wave4', help='Option to run with all of GS (waves 1, 3, and 4) instead of just waves 1 and 3 as base of models.', action="store_true", default=False)
#parser.add_argument('--justw4', help='Option to run with all of GS (waves 1, 3, and 4) instead of just waves 1 and 3 as base of models.', action="store_true", default=False)
parser.add_argument('--gsnorm', help='Option to run with all of GS, where all waves have been normalized together.', action="store_true", default=False)
parser.add_argument('--deathage', help='Option to run with only deaths over certain age in training.', action="store_true", default=False) # Not yet implemented
parser.add_argument('--age', help='Age to consider deaths over.', type = int, default=60) # Not yet implemented
parser.add_argument('--deathcause', help='Option to run with only certain causes of death in training.', action="store_true", default=False) # Not yet implemented
parser.add_argument('--cause', help='Cause of death to consider in training.', type = str, default=None) # Not yet implemented
parser.add_argument('--kinship', help='Option to run with kinship matrix, obtain Martingale residuals, and then obtain predictor running elnet on that.', action="store_true", default=False) # Not yet implemented
parser.add_argument('--sex', help='Option to sex-stratify models.', action="store_true", default=False)
parser.add_argument('--squared', help='Option to also include features^2 in models.', action="store_true", default=False)
parser.add_argument('--nonrandom', help='Option to assign batches to folds in ascending order instead of randomly.', action="store_true", default=False)
parser.add_argument('--cpg_bage', help='Option to obtain mortality predictor from methylation data instead of episcores.', action="store_true", default=False)
parser.add_argument('--cpg_subset', help='Option to obtain portality predictor from methylation data, subsetting to certain number of CpGs from mortality EWAS.', action="store_true", default=False)
parser.add_argument('--cpg_subset_n', type=int, help='Number (integer) of top CpGs to include in training (default is 900).', default = 900, choices=[100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1250, 1500, 1750, 2000, 3000, 4000, 5000, 7500, 10000])
parser.add_argument('--combo_bage', help='Option to include both episcores and CpGs in training. Options for subsets same as before.', action="store_true", default=False)

args = parser.parse_args()
tf = {True: "T", False: "F"}

####### Waves considered in model
if args.wave4:
    waves = [1,3,4]
else:
    waves = [1,3]

####### Parameters
sex = tf[args.sex] # Sex stratify? Will fit model for males and females separately
squared = tf[args.squared] # Squared CpGs? Will fit model considering CpG^2 values as well in training models
random_plates = tf[not args.nonrandom] # Plate assignment to folds
deathage = tf[args.deathage]
deathage_age = "F"
if args.deathage:
    deathage_age = args.age
deathcause = tf[args.deathcause]
deathcause_cause = "F"
if (args.deathcause) & (args.cause != None):
    deathcause_cause = args.cause
gsnorm = args.gsnorm
kinship = tf[args.kinship]
cpg_bage = tf[args.cpg_bage]
cpg_subset = tf[args.cpg_subset]
cpg_subset_n = ""
if cpg_subset == "T":
    cpg_subset_n = args.cpg_subset_n
combo_bage = tf[args.combo_bage]
if combo_bage == "T":
    cpg_bage = "T"

####### Stuff to run
gFolds = args.folds # True if we have yet to assign folds - if false skip this step 
runELNET_train = args.train # True if want to run training models
runELNET_test = args.test # True to test models
pTraining = False
pTesting = False

# Dictionaries for file naming
random_dic = {"T": "random", "F" : "nonrandom"}
gsnorm_dic = {True : "_gs20knorm", False: ""}
kinship_dic = {"T" : "_kinship", "F" : ""}
deathage_dic = {"T" : "_agedeath%s" % deathage_age, "F" : ""}
deathcause_dic = {"T" : "_causedeath_%s" % deathcause_cause, "F" : ""}
cpgpred_dic = {"T" : "_cpgpredictor", "F" : ""}
cpgpred_subset_dic = {"T" : "_cpgsubset%s" % cpg_subset_n, "F" : ""}
squared_dic = {"T" : "_squared", "F": ""}
combo_dic = {"T" : "_comboepiscorecpg", "F" : ""}

# Stuff yet to be implemented in code (obtain EpiScore projections for training/testing and fusing with GrimAge projections obtained online)
# pTraining = False # True if we want to obtain the episcore projections for the training data set and merge them with GrimAge obtained online
# pTesting = False # True if we want to obtain the episcore projections for the testing data set and merge them with GrimAge obtained online


# BASE INFO
######################################################
######################################################

# Output paths/where things will be stored
output_dir_folds = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/biological_age/elasticnet_models/cv_folds/"
output_dir_train = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/biological_age/elasticnet_models/elnet_train/"
output_dir_test = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/biological_age/elasticnet_models/elnet_test/"

# Where to find testing data
testing_dir = "/Cluster_Filespace/Marioni_Group/Elena/data/lbc_data/"

# Where to find target info
gs_df_10k = "/Cluster_Filespace/Marioni_Group/Elena/data/gs_data/gs10ktargets.tsv" # File requires ID of sample, age, sex, wave & batch columns && set and batch must follow W1_1 format
gs_df_20k = "/Cluster_Filespace/Marioni_Group/Elena/data/gs_data/gs20ktargets_fixedage.tsv"

# Kinship matrix
gs_20k_kinship = "/Cluster_Filespace/Marioni_Group/Elena/data/gs_data/gs_kinship_20k.RDS"

# Where to find prepped methylation data (betas)
gs_20k_methylation = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/chronological_age/data_prep/w1w3w4/methylation_training_noscale.rds"
lbc_methylation = "/Cluster_Filespace/Marioni_Group/Elena/data/lbc_data/lbc_bvals_cagetesting.rds"

# epic_rds = "/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds"

# Subset file from survival EWAS
mortality_ewas_subset = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/survival_ewas/coxph/subsets/coxph_survival_ewas_%s.tsv"

# File with GrimAge and Episcore projections, as well as sex and cAge (if they are available, if not these will not be used)
projections_training_10k = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/biological_age/data_prep/gs_w1w3_bage_variables_scaledgrim.tsv"
projections_training_20k = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/biological_age/data_prep/gs_w1w3w4_bage_variables_scaledgrim.tsv"
projections_training_gsnorm20k = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/biological_age/data_prep/gs_20k_bage_variables.tsv"
projections_testing = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/biological_age/data_prep/lbc_21w1_36w1_bage_variables.rds"

# File with ages and causes for death
death_20k = "/Cluster_Filespace/Marioni_Group/Elena/data/gs_data/gs20k_reasons_death_10032022.tsv"

# In case Episcores have yet to been calculated, calculate them and fuse them with the GrimAge projections that can obly be obtained online (not currently in use -> for future implementations)
# grim_training = ""
# grim_testing = ""

# General info for folds + models
nindivs_d = {1: 5087, 3: 4450, 4: 9241} # Individuals in each methylation wave
plates_d = {1: list(range(1,32)), 3: list(range(1,30)) + [31], 4: list(range(1,61))} # W1 has 31 batches and W3 has 30 batches, but note that batch 30 doesn't exist and batch 31 does, W4 has 60 batches
platesperfold_d = {1:[6,6,6,6,7], 3:[6,6,6,6,6], 4:[6,6,6,6,6,6,6,6,6,6]} # Number of batches per fold originating in each wave, each fold has around 1K indivds

# Color palette
palette = "custom2"


# ANALYSIS
######################################################
######################################################

now = datetime.now()
run_name = "%s%s%s%s%s%s%s%s" % (random_dic[random_plates], cpgpred_dic[cpg_bage], combo_dic[combo_bage], cpgpred_subset_dic[cpg_subset], gsnorm_dic[gsnorm], kinship_dic[kinship], deathage_dic[deathage], deathcause_dic[deathcause])

print(now.strftime("%d/%m/%Y %H:%M:%S"))
print("\n###########################################################")
print("##################################  Starting bAge analysis!")
print("###########################################################\n")

start = time.time()

print("Will be running: ")

print("Run name: %s" % run_name)
if gFolds:
    print("- Fold assignment")
if runELNET_train:
    print("- Training of predictive bAge models")
if runELNET_test:
    print("- Testing of predictive bAge models")

print("\n~ ~ ~ ~ ~ ~ ~ ~ ~ ~")

print("\nUser options:")
print("- Run name: %s" % run_name)
if cpg_bage == "T":
    print("- Making predictor with CpG methylation data")
if cpg_subset == "T":
    print("- Subsetting CpGs to a subset of %s from mortality EWAS" % cpg_subset_n)
    mortality_ewas_subset = mortality_ewas_subset % cpg_subset_n

# Modify some variables based on input
######################################################

if random_plates == False:
    ran = "nonrandom"
    plate_dict = plates_d
else:
    ran = "random"
    plate_dict = {1: random.sample(plates_d[1], len(plates_d[1])), 3: random.sample(plates_d[3], len(plates_d[3]))}

if waves == [1,3]:
    output_dir_folds = output_dir_folds + "w1w3/%s/" % ran
    output_dir_train = output_dir_train + "w1w3/%s/" % run_name
    output_dir_test = output_dir_test + "w1w3/%s/" % run_name
    output_dir_folds = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/chronological_age/elasticnet_models/cv_folds/w1w3/%s/" % ran
elif waves == [1,3,4]:
    output_dir_folds = output_dir_folds + "w1w3w4/%s/" % ran
    output_dir_train = output_dir_train + "w1w3w4/%s/" % run_name
    output_dir_test = output_dir_test + "w1w3w4/%s/" % run_name
    output_dir_folds = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/chronological_age/elasticnet_models/cv_folds/w1w3w4/%s/" % ran

if not os.path.exists(output_dir_folds):
    os.makedirs(output_dir_folds)
if not os.path.exists(output_dir_train):
    os.makedirs(output_dir_train)
if not os.path.exists(output_dir_test):
    os.makedirs(output_dir_test)


# Get folds and explore
######################################################

if gFolds == True:
    # Import table
    df = pd.read_table(gs_df, index_col = "Sample_Name")
    df = df.drop("Numb", axis = 1)
    # Calculate folds
    folds = get_folds(df = df, plate_col = "Batch", index_col = "Sample_Name", waves = waves, plate_dictionary = plate_dict, platesperfold_dictionary = platesperfold_d, output_dir = output_dir_folds, name = ran)
    # Add folds to table
    df["Fold"] = folds
    # Export table with folds
    df.to_csv(output_dir_folds + "gs_basic_folds_%s.tsv" % ran, sep = "\t")
    # Explore
    explore_folds(df, ran, output_dir_folds)


# Elastic net model data prep - projections, and creation of input file 
#######################################################################

# if pTraining == True:
#     print("TBD.")


# Elastic net model training
#######################################################################

if runELNET_train == True:
    if pTraining == False: # Meaning that projection file has already been provided
        if waves == [1,3]:
            args = [projections_training_10k, cpg_bage, gs_20k_methylation, cpg_subset, mortality_ewas_subset, combo_bage, death_20k, output_dir_folds + "gs_basic_folds_%s.tsv" % ran, sex, squared, run_name, kinship, gs_20k_kinship, deathage_age, deathcause_cause, output_dir_train] # Arguments for R script
        elif waves == [1,3,4]:
            if gsnorm == False:
                args = [projections_training_20k, cpg_bage, gs_20k_methylation, cpg_subset, mortality_ewas_subset, combo_bage, death_20k, output_dir_folds + "gs_basic_folds_%s.tsv" % ran, sex, squared, run_name, kinship, gs_20k_kinship, deathage_age, deathcause_cause, output_dir_train] # Arguments for R script
            else:
                args = [projections_training_gsnorm20k, cpg_bage, gs_20k_methylation, cpg_subset, mortality_ewas_subset, combo_bage, death_20k, output_dir_folds + "gs_basic_folds_%s.tsv" % ran, sex, squared, run_name, kinship, gs_20k_kinship, deathage_age, deathcause_cause, output_dir_train] # Arguments for R script
    else:
        print("TBD.")
    subprocess.call("biological_age/elnet_models.R --input %s --cpg_bage %s --cpg_file %s --cpg_subset %s --cpg_subset_file %s --combo_bage %s --agecause %s --folds %s --sex %s --squared %s --name %s --kinship %s --kinship_file %s --deathage %s --deathcause %s --out %s" % tuple(args), shell  = True)


# Elastic net model data prep - projections for testing data
#######################################################################

# if pTesting == True:
#     print("TBD.")


# Elastic net model testing
#######################################################################

if runELNET_test == True:
    args = [projections_testing, cpg_bage, combo_bage, lbc_methylation, cpg_subset, testing_dir + "lbc_testing_info.tsv", output_dir_folds + "gs_basic_folds_%s.tsv" % ran, output_dir_train + "elnet_coefficients_%s.tsv" % run_name, run_name, output_dir_test] # Arguments for R script
    subprocess.call("biological_age/elnet_test.R --input %s --cpg_bage %s --combo_bage %s --cpg_file %s --cpg_subset %s --pheno %s --training %s --weights %s --name %s --out %s" % tuple(args), shell  = True)
    
    # Scatter plots comparing age and tte to prediction
    explore_test(prediction_file = output_dir_test + "predictions.tsv", bage_pred_col = "bAge", grim_pred_col = "GrimAge", true_col = "Age", dataset_col = "Cohort", output_dir = output_dir_test, palette = palette)
    explore_test(prediction_file = output_dir_test + "predictions.tsv", bage_pred_col = "bAge", grim_pred_col = "GrimAge", true_col = "TTE", dataset_col = "Cohort", output_dir = output_dir_test, palette = palette, tte = True)

    # Scatter plots comparing bAge to GrimAge
    correlation_predictors(prediction_file = output_dir_test + "predictions.tsv", bage_pred_col = "bAge", grim_pred_col = "GrimAge", dataset_col = "Cohort", output_dir = output_dir_test, palette = palette)
    if cpg_bage == "T":
        # Scatterplots comparing bAge to bAge based on CpGs
        episcore_bage = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/biological_age/elasticnet_models/elnet_test/w1w3w4/random/predictions.tsv"
        correlation_predictors(prediction_file = episcore_bage, bage_pred_col = "bAge", grim_pred_col = "GrimAge", dataset_col = "Cohort", output_dir = output_dir_test, palette = palette, cpg_bage_comp = True, prediction_file_bagecpg = output_dir_test + "predictions.tsv")

    # Make forest plot comparing to GrimAge, for all datasets and per subset
    forest_plot(cox_file = output_dir_test + "coxph_testing.tsv", predictions_file = output_dir_test + "predictions.tsv", output_dir = output_dir_test, palette = palette)
    forest_plot(cox_file = output_dir_test + "coxph_testing.tsv", predictions_file = output_dir_test + "predictions.tsv", output_dir = output_dir_test, palette = palette, coh = True)

    # Make Kaplan-Meier plots, for all datasets and per subset
    kaplanmeier_plot(prediction_file = output_dir_test + "predictions.tsv", grim_col = "GrimAgeAccel", bage_col = "bAgeAccel", dataset_col = "Cohort", output_dir = output_dir_test, palette = palette, death_col = "Dead", coh = False)
    kaplanmeier_plot(prediction_file = output_dir_test + "predictions.tsv", grim_col = "GrimAgeAccel", bage_col = "bAgeAccel", dataset_col = "Cohort", output_dir = output_dir_test, palette = palette, death_col = "Dead", coh = True)


print("\nDone!")
total_time = time.time() - start
print("\nTotal time: %s minutes (%s hours)." % (round(total_time/60, 3), round((total_time/60)/60, 3)))