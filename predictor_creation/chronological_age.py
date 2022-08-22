#!/usr/bin/env python
# -*- coding: utf-8 -*-
# By Elena

################################################################
# Script that coordinates all data prep pre-running elastic net
# + fitting of elastic net model for cAge
################################################################

# Run on cluster
# For VSCode + cluster: export PATH=/home/ebernab3/anaconda3/bin:$PATH
# Recommended way of running this: nohup python -u chronological_age.py [options] > logs/output_file.out &

import os
from os.path import exists
import numpy as np
import pandas as pd
import random
import subprocess
import time
from datetime import datetime
import argparse
from chronological_age.data_diagnostics import *
from chronological_age.elnet_folds import *
from chronological_age.elnet_test import *


# USER INPUT
######################################################
######################################################

parser = argparse.ArgumentParser(description='cAge analysis pipeline.')

# Toggles to run bits of code - needs at least one of these flags or pipeline won't do anything
parser.add_argument('--folds', help='Option to assign batches to folds. WARNING: Needs to be run at least once prior to training.', action="store_true", default=False)
parser.add_argument('--trainprep', help='Option to prepare training data prior to fitting models. WARNING: Needs to be run at least once prior to training.', action="store_true", default=False)
parser.add_argument('--diagnose', help='Option to diagnose data sets used in analysis prior to fitting models.', action="store_true", default=False)
parser.add_argument('--train', help='Option to train models.', action="store_true", default=False)
parser.add_argument('--test', help='Option to test models.', action="store_true", default=False)

# Parameters that start kicking in in fold assignment
parser.add_argument('--wave4', help='Option to run with all of GS (waves 1, 3, and 4) instead of just waves 1 and 3 as base of models.', action="store_true", default=False)
parser.add_argument('--external', help='Option to use external data to train (LBC21/36 and GEO datasets).', action="store_true", default=False)
parser.add_argument('--nonrandom', help='Option to assign batches to folds in ascending order instead of randomly.', action="store_true", default=False)

# Parameters that start kicking in data prep
parser.add_argument('--gsnorm', help='Option to use the 20k normalized GS methylation file instead of independent wave files.', action="store_true", default=False)

# Parameters that start kicking in diagnostics
parser.add_argument('--cpg', help='CpG that we want to obtain the trajectory of in our diagnostics, default is ELOVL2.', type=str, default="cg16867657")

# Parameters that start kicking in in training
parser.add_argument('--name', type=str, help='Name to add to generic run name.', default="")
parser.add_argument('--cpg_train', help='Option to scale per CpG in training, per CV fold.', action="store_true", default=False)
parser.add_argument('--sample_train', help='Option to scale per sample in training.', action="store_true", default=False)
parser.add_argument('--subset', help='Option to limit training to most associated CpGs to age.', action="store_true", default=False)
parser.add_argument('--subset_n', type=int, help='Number (integer) of Ks of top CpGs to include in training (default is 20, so top 20K CpGs).', default = 20, choices=[1,2,3,4,5,10,15,20,35,50,75,100,200,300])
parser.add_argument('--loo', help='Leave-one-out option when training and testing with GEO datasets (one model will be fit excluding the set to then be tested on, for all sets considered except GS and LBC).', action="store_true", default=False)
parser.add_argument('--sex', help='Option to sex-stratify models.', action="store_true", default=False)
parser.add_argument('--squared', help='Option to also include CpG^2 in models. Will include CpG^2 for all CpGs already in model.', action="store_true", default=False)
parser.add_argument('--squared_subset', help='Option to also include CpG^2 in models, selecting from CpGs most associated to age^2. If selected, it overrides --squared.', action="store_true", default=False)
parser.add_argument('--squared_subset_n', type=int, help='Number (integer) of top CpGs^2 for age to include in training (default is 900).', default = 900, choices=[100,200,300,400,500,600,700,800,900,1000,1250,1500,2000,3000,5000,7500,10000,12500,15000,20000])
parser.add_argument('--cpg2', help='Option for squared CpGs to be obtained from age + CpG2 EWAS, not from age^2 EWAS.', action="store_true", default=False)
parser.add_argument('--cpg2_n', type=float, help='Number (integer) of top CpGs^2 for age to include in training (default is 5K).', default = 5., choices=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.25,1.5,2,3,4,5,10,15,20,35,50,75,100,200,300])
parser.add_argument('--age2cpg2', help='Option to include both CpG^2 derived from CpG^2-age EWAS and from age^2 EWAS. Subsets will be determined by already described options.', action="store_true", default=False)
parser.add_argument('--trandom', help='Option to ignore fold assignment and use totally random folds in training (totally random).', action="store_true", default=False)
parser.add_argument('--lasso', help='Option to use lasso instead of elastic net in training of models.', action="store_true", default=False)
parser.add_argument('--logage', help='Option to train on log(age) instead of age.', action="store_true", default=False)

# Parameters that start kicking in testing
parser.add_argument('--sample_test', help='Option to scale per sample in testing.', action="store_true", default=False)
parser.add_argument('--cpg_test', help='Option to scale CpG in testing, per cohort/study.', action="store_true", default=False)
parser.add_argument('--lbc_test', help='Include LBC in testing.', action="store_true", default=False)
parser.add_argument('--logage20', help='Option to test on log(age) for under 18s (non-log for older individuals).', action="store_true", default=False)
parser.add_argument('--cohort', type = str, help='Option to select cohort to test in', default="all")

args = parser.parse_args()
tf = {True: "T", False: "F"}

####### Waves considered in model
if args.wave4:
    waves = [1,3,4]
else:
    waves = [1,3]

####### Parameters
gsnorm = args.gsnorm # Option to use combined 20k normalized m-val file (only relevant if running with all waves)
if args.wave4 == False: # Option only active if wave 4 is selected
    gsnorm = False
cpg_traject = args.cpg # CpG of interest for CpG trajectory plots
cpg = tf[args.cpg_train] # True to scale per CpG in training data
sample = tf[args.sample_train] # True to scale per sample in training data
cpg_test = tf[args.cpg_test] # True to scale per CpG in test data
sample_test = tf[args.sample_test] # True to scale per sample in test data
subset = tf[args.subset] # Subset CpGs instead of full circa 400K? 
subset_n = args.subset_n # Number of Ks to subset CpGs to
external = tf[args.external] # External data in training?
loo = args.loo # Leave-one-out model option - will fit in total one model per GEO set (8 GEOs), where said set is left out
if loo:
    external = "T"
sex = tf[args.sex] # Sex stratify? Will fit model for males and females separately
squared = tf[args.squared] # Squared CpGs? Will fit model considering CpG^2 values as well in training models
squared_subset = tf[args.squared_subset] # Squared CpGs from selected bunch most associated to age^2? Will fit model considering these CpG^2 as well in training models
squared_subset_n = args.squared_subset_n # Number of top associated CpGs to age^2 to include as CpG^2 in training model (analogous to subset_n)
squared_cpg2 = tf[args.cpg2]
squared_subset_n_cpg2 = args.cpg2_n
if squared_subset_n_cpg2.is_integer():
    squared_subset_n_cpg2 = int(squared_subset_n_cpg2)
age2cpg2 = tf[args.age2cpg2]
if squared_subset == "T":
    squared = "F" # Override this option in case both called
if (squared_cpg2 == "T") & (age2cpg2 == "F"):
    squared_subset_n = squared_subset_n_cpg2
elif (squared_cpg2 == "T") & (age2cpg2 == "T"):
    squared_subset_n_cpg2 = squared_subset_n_cpg2
coh = args.cohort
coh_f = "F"
if coh != "all":
    coh_f = "T"

random_plates = tf[not args.nonrandom] # Plate assignment to folds
totally_random = tf[args.trandom] # Totally random option? Will ignore fold assignment and run totally random folds
lasso = tf[args.lasso] # Lasso instead of elnet?
lbc = tf[args.lbc_test] # LBC in testing?
if (external == "T") & (loo == False):
    lbc = "F" # If using external in training LBC can't be used in testing as it's part of the training, unless LOO is being run
if args.name != "":
    if (squared_subset == "T") or (subset == "T"):
        name = "_%s_%s" % ("gs20k_lmewas_scaledewas", args.name)
    else:
        name = "_%s" % args.name
else:
    if (squared_subset == "T") or (subset == "T"):
        name = "_gs20k_lmewas_scaledewas"
    else:
        name = ""
logage = tf[args.logage]
logage20 = tf[args.logage20]
if args.logage20:
    logage = "T" # For file naming purposes


####### Stuff to run
diagnose = args.diagnose # True if we want to explore data before starting analysis
gFolds = args.folds # True if we have yet to assign folds - if false skip this step
runELNET_train_prep = args.trainprep # True if want to prepare data for training of elastic net models
runELNET_train = args.train # True if want to run training models
runELNET_test = args.test # True to test models


# BASE INFO
######################################################
######################################################

# Output paths/where things will be stored
output_dir_diagnose = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/chronological_age/data_explore/"
output_dir_prep = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/chronological_age/data_prep/"
output_dir_folds = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/chronological_age/elasticnet_models/cv_folds/"
output_dir_train = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/chronological_age/elasticnet_models/elnet_train/"
output_dir_test = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/chronological_age/elasticnet_models/elnet_test/"

# Where to find testing data
testing_dir = "/Cluster_Filespace/Marioni_Group/Elena/data/geo_data/"
lbc_testing_samp = "/Cluster_Filespace/Marioni_Group/Elena/data/lbc_data/lbc_targets_cagetesting.tsv"
lbc_testing_meth = "/Cluster_Filespace/Marioni_Group/Elena/data/lbc_data/lbc_bvals_cagetesting.rds"

# CpGs to limit analysis to
cpgs = "/Cluster_Filespace/Marioni_Group/Elena/gs_osca/data/cpgs_tokeep.txt"

# Where to find target info
gs_df_10k = "/Cluster_Filespace/Marioni_Group/Elena/data/gs_data/gs10ktargets.tsv" # File requires ID of sample, age, sex, wave & batch columns && set and batch must follow W1_1 format
gs_df_20k = "/Cluster_Filespace/Marioni_Group/Elena/data/gs_data/gs20ktargets_fixedage.tsv"
lbc_df = "/Cluster_Filespace/Marioni_Group/Elena/data/lbc_data/lbc_targets_3489.tsv" # File with sample, age, sex, wave in LBC if using non-GS data to train as well
geo_df = "/Cluster_Filespace/Marioni_Group/Elena/data/geo_data/test_samples.tsv" # File with age and sex

# Where to find subset of CpGs of interest
#subset_f = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/subsets/age_ewas_%sK_gs20k.tsv" # Format should be .tsv with first column with CpG names to subset
#squared_subset_f = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/subsets/age2_ewas_%s_gs20k.tsv" # Format same as above
subset_f = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/subsets/cpg_lm_%sK_scaled_gs20k.tsv"
squared_subset_f = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/subsets/age2_lm_%s_gs20k.tsv"
squared_subset_f_cpg2 = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_explore/w1w3w4/subsets/cpg2_lm_%sK_gs20k_scaled.tsv" # Format same as above

# Where to find methylation data
w1_methylation_rds = "/Cluster_Filespace/Marioni_Group/GS/GS_methylation/norm_mvals_5087.rds"
w3_methylation_rds = "/Cluster_Filespace/Marioni_Group/GS/GS_methylation/wave3_mvals.rds"
w4_methylation_rds = "/Cluster_Filespace/Marioni_Group/GS/GS_methylation/wave4/w4-mvals.rds"
gs_methylation_rds = "/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/mvals.rds"
epic_rds = "/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds"
lbc_methylation_rds = "/Cluster_Filespace/Marioni_Group/Elena/data/lbc_data/lbc_mvals_3489.rds"
geo_methylation_rds = "/Cluster_Filespace/Marioni_Group/Elena/data/geo_data/test_methylation.RDS"

# General info for folds + models
nindivs_d = {1: 5087, 3: 4450, 4: 9241} # Individuals in each methylation wave
plates_d = {1: list(range(1,32)), 3: list(range(1,30)) + [31], 4: list(range(1,61))} # W1 has 31 batches and W3 has 30 batches, but note that batch 30 doesn't exist and batch 31 does, W4 has 60 batches
platesperfold_d = {1:[6,6,6,6,7], 3:[6,6,6,6,6], 4:[6,6,6,6,6,6,6,6,6,6]} # Number of batches per fold originating in each wave, each fold has around 1K indivds

# External cohorts considered
cohorts = ["LBC21", "LBC36", "GSE40279", "GSE42861", "GSE72775", "GSE78874", "GSE41169", "GSE53740", "GSE72773", "GSE72777"]
ext_select = "LBC21,LBC36,GSE40279,GSE42861" # Default when using external data

# Color palette
palette = "custom2"

# Dictionaries for file naming
gsnorm_dic = {True : "_gs20knorm", False: ""}
random_dic = {"T": "random", "F" : "nonrandom"}
scalesample_dic = {"T": "_scalesample", "F": "_noscalesample"}
scalecpg_dic = {"T": "_scalecpg", "F": "_noscalecpg"}
squared_dic = {"T" : "_squared", "F" : ""}
external_dic = {"T" : "_withexternal", "F" : ""}
external_dic_ii = {"T" : "_external", "F" : ""}
totallyrandom_dic = {"T": "_randomizedfolds", "F" : ""}
lasso_dic = {"T": "_lasso", "F" : ""}
subset_dic = {"T": "_subset%sK" % subset_n , "F" : ""}
squared_subset_dic = {"T": "_squaredsubset%s" % squared_subset_n , "F" : ""}
if (squared_cpg2 == "T") & (age2cpg2 == "F"):
    squared_subset_dic = {"T": "_squaredsubset%sK" % squared_subset_n, "F" : ""}
if (squared_cpg2 == "T") & (age2cpg2 == "T"):
    squared_subset_dic = {"T": "_squaredsubset%s_squaredsubsetcpg2lm%sK" % (squared_subset_n, squared_subset_n_cpg2), "F" : ""}
sex_dic = {"T" : "_sexstratified", "F" : ""}
loo_dic = {True : "_loo", False: ""}
lbc_dic = {"T" : "_lbc", "F": ""}
cpg2_dic = {"T" : "_cpg2lm", "F" : ""}
if age2cpg2 == "T":
    cpg2_dic = {"T" : "", "F" : ""}
logage_dic = {"T" : "_logage", "F" : ""}
logage20_dic = {"T" : "_logage20", "F" : ""}
test_cohort_dic = {"T" : "_test%s" % coh, "F" : ""}

# ANALYSIS
######################################################
######################################################

now = datetime.now()
#time4file = now.strftime("%d-%m-%Y_%H-%M-%S")
run_name = "%s%s%s%s%s%s%s%s" % (random_dic[random_plates], gsnorm_dic[gsnorm], scalesample_dic[sample], scalecpg_dic[cpg], subset_dic[subset], external_dic_ii[external], logage_dic[logage], name)
#log_filename = 'logs/cage_%s_%s.out' % (run_name, time4file)
#os.makedirs(os.path.dirname(log_filename), exist_ok=True)

print(now.strftime("%d/%m/%Y %H:%M:%S"))
print("\n###########################################################")
print("##################################  Starting cAge analysis!")
print("###########################################################\n")

start = time.time()


# Modify some variables based on input
######################################################

print("Will be running: ")

if gFolds:
    print("- Fold assignment")
if runELNET_train_prep:
    print("- Data prep")
if diagnose:
    print("- Data exploration/diagnosis")
if runELNET_train:
    print("- Training of predictive cAge models")
if runELNET_test:
    print("- Testing of predictive cAge models")

print("\n~ ~ ~ ~ ~ ~ ~ ~ ~ ~")

print("\nUser options:")
print("- Run name: %s" % run_name)

if logage == "T":
    if logage20 == "T":
        print("- Running on log(age) instead of age for under 20s")
    elif logage20 == "F":
        print("- Running on log(age) instead of age")
if cpg == "T":
    print("- Scaling per CpG (per fold) in training data")
if sample == "T":
    print("- Scaling per sample in training data")
if cpg_test == "T":
    print("- Scaling per CpG (per fold) in testing data")
if sample_test == "T":
    print("- Scaling per sample in testing data")
if subset == "T":
    subset_f = subset_f % subset_n 
    print("- Subsetting to %sK CpGs" % subset_n)

# Random plates for fold assignment
if random_plates == False:
    print("- Batch folds: Non-random fold allocation")
    ran = "nonrandom"
    plate_dict = plates_d
else:
    print("- Batch folds:  Random fold allocation")
    ran = "random"
    plate_dict = {1: random.sample(plates_d[1], len(plates_d[1])), 3: random.sample(plates_d[3], len(plates_d[3])), 4: random.sample(plates_d[4], len(plates_d[4]))}

# Directory names based on training data
if waves == [1,3]:
    print("- GS: Training with GS waves 1 and 3 as base")
    gs_df = gs_df_10k
    output_dir_folds = output_dir_folds + "w1w3/%s/" % ran
    output_dir_prep = output_dir_prep + "w1w3/"
    output_dir_diagnose = output_dir_diagnose + "w1w3/"
    output_dir_train = output_dir_train + "w1w3/%s/" % run_name
    output_dir_test = output_dir_test + "w1w3/%s/" % run_name
elif waves == [1,3,4]:
    print("- GS: Training with GS waves 1, 3 and 4 as base")
    if gsnorm == True:
        print("- GS: Using normalized 20k methylation data")
    gs_df = gs_df_20k
    output_dir_folds = output_dir_folds + "w1w3w4/%s/" % ran
    output_dir_prep = output_dir_prep + "w1w3w4/"
    output_dir_diagnose = output_dir_diagnose + "w1w3w4/"
    output_dir_train = output_dir_train + "w1w3w4/%s/" % run_name
    output_dir_test = output_dir_test + "w1w3w4/%s/" % run_name

if not os.path.exists(output_dir_diagnose):
    os.makedirs(output_dir_diagnose)
if not os.path.exists(output_dir_prep):
    os.makedirs(output_dir_prep)
if not os.path.exists(output_dir_folds):
    os.makedirs(output_dir_folds)
if not os.path.exists(output_dir_train):
    os.makedirs(output_dir_train)
if not os.path.exists(output_dir_test):
    os.makedirs(output_dir_test)

if external == "T":
    print("- External: Training will also include LBC + GEO data")
if squared == "T":
    print("- Squared: Training will also include CpG^2 features")
if squared_subset == "T":
    squared_subset_f = squared_subset_f % squared_subset_n
    if squared_cpg2 == "F":
        print("- Squared: Training will also include a subset of %s CpG^2 features" % squared_subset_n)
        squared_subset_f_cpg2 = "NA"
    elif (squared_cpg2 == "T") & (age2cpg2 == "F"):
        squared_subset_f = squared_subset_f_cpg2 % squared_subset_n
        squared_subset_f_cpg2 = "NA"
        print("- Squared: Training will also include a subset of %sK CpG^2 features" % squared_subset_n)
        print("- Squared: CpG^2 selected from linear models for age + CpG^2 (not age^2 EWAS)")
    elif (squared_cpg2 == "T") & (age2cpg2 == "T"):
        squared_subset_f_cpg2 = squared_subset_f_cpg2 % squared_subset_n_cpg2
        print("- Squared: Training will also include a subset of %sK CpG^2 from CpG^2 lm models and %s CpG^2 features from Age^2 EWAS models" % (squared_subset_n_cpg2, squared_subset_n))

if sex == "T":
    print("- Sex: Training will be done sex-stratified")
if totally_random == "T":
    print("- Random: Training will ignore fold allocation and totally randomize samples to folds")
if lasso == "T":
    print("- Lasso: Training will include fitting of lasso models instead of elastic net")
if loo:
    print("- Leave-one-out: models will be fit leaving one GEO/LBC set out for testing + a final model w/ all")
if lbc == "T":
    print("- LBC: Including LBC21 and 36 in testing")



# Get folds + explore
######################################################

if gFolds == True:

    print("\n################### Getting folds\n")

    # Import table
    df = pd.read_table(gs_df, index_col = "Sample_Name")
    if waves == [1,3]:
        df = df.drop("Numb", axis = 1)
        index_col = "Sample_Name"
    elif waves == [1,3,4]:
        index_col = "ID"
    # Calculate folds
    folds = get_folds(df = df, plate_col = "Batch", index_col = index_col, waves = waves, plate_dictionary = plate_dict, platesperfold_dictionary = platesperfold_d, output_dir = output_dir_folds, name = ran)
    # Add folds to table
    df["Fold"] = folds
    # Export table with folds
    df.to_csv(output_dir_folds + "gs_basic_folds_%s.tsv" % ran, sep = "\t")

    print("\n################### Exploring folds\n")

    explore_folds(df, ran, output_dir_folds)


# Add extra info if also training with non-GS data
######################################################

if external == "T":

    if gFolds == True:
        
        print("\n################### Adding external info to fold file\n")

        coi = ["age", "sex", "cohort"]
        if waves == [1,3]:
            df = pd.read_table(output_dir_folds + "gs_basic_folds_%s.tsv" % ran, index_col = "Sample_Sentrix_ID")
        else:
            df = pd.read_table(output_dir_folds + "gs_basic_folds_%s.tsv" % ran, index_col = "ID")
        df = df.rename({"Set" : "cohort"}, axis = 1)
        lbc = pd.read_table(lbc_df, index_col = "Basename")
        geo = pd.read_table(geo_df, index_col = 0, header = 0, names = coi)

        # Empty cols
        lbc["Fold"] = np.nan
        geo["Fold"] = np.nan

        # Add folds in LBC
        max_fold = max(df["Fold"]) # Max fold
        for c in set(lbc["cohort"]):
            for w in set(lbc.loc[lbc["cohort"] == c, "WAVE"]):
                max_fold += 1
                lbc.loc[(lbc["cohort"] == c) & (lbc["WAVE"] == w), "Fold"] = max_fold

        # Add folds in GEO
        for g in set(geo["cohort"]):
            max_fold += 1
            geo.loc[geo["cohort"] == g, "Fold"] = max_fold

        # Fuse external targets
        external_df = pd.concat([df[coi + ["Fold"]], lbc[coi + ["Fold"]], geo[coi + ["Fold"]]], axis = 0)

        # Remove NAs in age column
        external_df = external_df[external_df["age"].notna()]
        external_df.to_csv(output_dir_folds + "gs_basic_folds_%s_withexternal.tsv" % ran, sep = "\t", index_label = "ID")

        # Explore folds
        explore_folds(external_df, ran, output_dir_folds, external = "T")


# Elastic net model data prep
######################################################

if runELNET_train_prep == True:
    
    print("\n################### Methylation data prep for training\n")

    if waves == [1,3]:
        if external == "F":
            args = [w1_methylation_rds, w3_methylation_rds, output_dir_folds + "gs_basic_folds_%s.tsv" % ran, epic_rds, "noscale", output_dir_prep] #Arguments for R script
            subprocess.call("chronological_age/elnet_dataprep.R --wave1 %s --wave3 %s --pheno %s --epic %s --name %s --out %s" % tuple(args), shell  = True)
        elif external == "T":
            args = [w1_methylation_rds, w3_methylation_rds, lbc_methylation_rds, geo_methylation_rds, output_dir_folds + "gs_basic_folds_%s_withexternal.tsv" % ran, epic_rds, "noscale_external", output_dir_prep] #Arguments for R script
            subprocess.call("chronological_age/elnet_dataprep.R --wave1 %s --wave3 %s --lbc %s --geo %s --pheno %s --epic %s --name %s --out %s" % tuple(args), shell  = True)
           
    elif waves == [1,3,4]:
        if external == "F":
            if gsnorm == True:
                args = [gs_methylation_rds, cpgs, output_dir_folds + "gs_basic_folds_%s.tsv" % ran, epic_rds, "noscale_gs20knorm", output_dir_prep] #Arguments for R script
                subprocess.call("chronological_age/elnet_dataprep.R --gs20k %s --cpgs %s --pheno %s --epic %s --name %s --out %s" % tuple(args), shell  = True)
            else:
                args = [w1_methylation_rds, w3_methylation_rds, w4_methylation_rds, cpgs, output_dir_folds + "gs_basic_folds_%s.tsv" % ran, epic_rds, "noscale", output_dir_prep] #Arguments for R script
                subprocess.call("chronological_age/elnet_dataprep.R --wave1 %s --wave3 %s --wave4 %s --cpgs %s --pheno %s --epic %s --name %s --out %s" % tuple(args), shell  = True)
        elif external == "T":
            if gsnorm == True:
                args = [gs_methylation_rds, lbc_methylation_rds, geo_methylation_rds, cpgs, output_dir_folds + "gs_basic_folds_%s_withexternal.tsv" % ran, epic_rds, "noscale_gs20knorm_external", output_dir_prep] #Arguments for R script
                subprocess.call("chronological_age/elnet_dataprep.R --gs20k %s --lbc %s --geo %s --cpgs %s --pheno %s --epic %s --name %s --out %s" % tuple(args), shell  = True)           
            else:
                args = [w1_methylation_rds, w3_methylation_rds, w4_methylation_rds, lbc_methylation_rds, geo_methylation_rds, cpgs, output_dir_folds + "gs_basic_folds_%s_withexternal.tsv" % ran, epic_rds, "noscale_external", output_dir_prep] #Arguments for R script
                subprocess.call("chronological_age/elnet_dataprep.R --wave1 %s --wave3 %s --wave4 %s --lbc %s --geo %s --cpgs %s --pheno %s --epic %s --name %s --out %s" % tuple(args), shell  = True)
           

# Data diagnostics
######################################################

if diagnose == True:

    print("\n################### Diagnostics of data\n")

    # CpG trajectory
    #subprocess.call("chronological_age/data_diagnostics.R --cpgprep TRUE --pcaprep FALSE --cpg %s --meth %s --target %s --name %s --outputdir %s" % (cpg_traject, output_dir_prep + "methylation_training_noscale%s%s.rds" % (gsnorm_dic[gsnorm], external_dic_ii[external]), output_dir_folds + "gs_basic_folds_%s%s.tsv" % (ran, external_dic[external]), "noscale%s%s" % (gsnorm_dic[gsnorm],external_dic_ii[external]), output_dir_prep), shell  = True)
    if external == "T":
        cpg_trajectory_external(prepped_file = output_dir_prep + "%s_trajectory_noscale%s%s.tsv" % (cpg_traject, gsnorm_dic[gsnorm], external_dic_ii[external]), cpg = cpg_traject, output_dir = output_dir_diagnose, name = "_noscale%s%s" % (gsnorm_dic[gsnorm], external_dic_ii[external]))
    else:
        cpg_trajectory(prepped_file = output_dir_prep + "%s_trajectory_noscale%s%s.tsv" % (cpg_traject, gsnorm_dic[gsnorm], external_dic_ii[external]), cpg = cpg_traject, output_dir = output_dir_diagnose, name = "_noscale%s%s" % (gsnorm_dic[gsnorm], external_dic_ii[external]))
    
    # Age distribution
    #age_dist(target = output_dir_folds + "gs_basic_folds_%s%s.tsv" % (ran, external_dic[external]), output_dir = output_dir_diagnose, cohort_col = "cohort", age_col = "age", external = args.external, name = "noscale%s%s" % (gsnorm_dic[gsnorm],external_dic_ii[external]))

    # Sex distribution
    #sex_dist(target = output_dir_folds + "gs_basic_folds_%s%s.tsv" % (ran, external_dic[external]), output_dir = output_dir_diagnose, cohort_col = "cohort", sex_col = "sex", external = args.external, name = "noscale%s%s" % (gsnorm_dic[gsnorm],external_dic_ii[external]))
    
    # PCA
    #subprocess.call("chronological_age/data_diagnostics.R --cpgprep FALSE --pcaprep TRUE --meth %s --name %s --outputdir %s" % (output_dir_prep + "methylation_training_noscale%s%s.rds" % (gsnorm_dic[gsnorm], external_dic_ii[external]), "noscale%s%s" % (gsnorm_dic[gsnorm],external_dic[external]), output_dir_diagnose), shell  = True)
    #pca(external = args.external, prepped_file = output_dir_diagnose + "methbetavals_noscale%s_withexternal.csv" % (gsnorm_dic[gsnorm]), target = output_dir_folds + "gs_basic_folds_%s%s.tsv" % (ran, external_dic[external]), name = "noscale%s%s" % (gsnorm_dic[gsnorm], external_dic_ii[external]), output_dir = output_dir_diagnose)
    #pca_plot(pca_df_f = output_dir_diagnose + "noscale%s%s_PCA.tsv" % (gsnorm_dic[gsnorm],external_dic_ii[external]), expvar_df_f = output_dir_diagnose + "noscale%s%s_PCA_explainedvariance.tsv" % (gsnorm_dic[gsnorm],external_dic_ii[external]), target_f = output_dir_folds + "gs_basic_folds_%s%s.tsv" % (ran, external_dic[external]), output_dir = output_dir_diagnose, name = "noscale%s%s" % (gsnorm_dic[gsnorm],external_dic_ii[external]), external = args.external)

# Elastic net model training
######################################################

if runELNET_train == True:

    print("\n################### Training of models\n")
    if loo == True:
        print("**Leave one out framework!**\n")
        i = 0
        for cohort in cohorts:
            i += 1
            print("Training models excluding %s (%s/%s)...\n" % (cohort, i, len(cohorts)+1))
            ext_select_loo = ",".join([i for i in cohorts if i != cohort])
            print("GEO/LBC sets in model:\n%s\n" % ext_select_loo)
            args = ["methylation_training_noscale%s%s.rds" % (gsnorm_dic[gsnorm], external_dic_ii[external]), output_dir_prep, output_dir_folds + "gs_basic_folds_%s%s.tsv" % (ran, external_dic[external]), run_name + "_loo_%s" % cohort, sex, squared, external, ext_select_loo, totally_random, lasso, cpg, sample, subset, subset_f, squared_subset, squared_subset_f, squared_subset_f_cpg2, squared_cpg2, age2cpg2, logage, output_dir_train]
            subprocess.call("chronological_age/elnet_models.R --meth %s --methdir %s --pheno %s --name %s --sex %s --squared %s --external %s --ext_select %s --random %s --lasso %s --cpg %s --sample %s --subset %s --subfil %s --squaredsubset %s --squaredsubsetfil %s --squaredsubsetfil_cpg2 %s --squaredcpg2 %s --age2cpg2 %s --logage %s --out %s" % tuple(args), shell  = True)
        # Final model with all cohorts
        print("Training final model with all cohorts (%s/%s)...\n" % (len(cohorts)+1, len(cohorts)+1))
        ext_select_loo = ",".join(cohorts)
        args = ["methylation_training_noscale%s%s.rds" % (gsnorm_dic[gsnorm], external_dic_ii[external]), output_dir_prep, output_dir_folds + "gs_basic_folds_%s%s.tsv" % (ran, external_dic[external]), run_name + "_loo_all", sex, squared, external, ext_select_loo, totally_random, lasso, cpg, sample, subset, subset_f, squared_subset, squared_subset_f, squared_subset_f_cpg2, squared_cpg2, age2cpg2, logage, output_dir_train]
        subprocess.call("chronological_age/elnet_models.R --meth %s --methdir %s --pheno %s --name %s --sex %s --squared %s --external %s --ext_select %s --random %s --lasso %s --cpg %s --sample %s --subset %s --subfil %s --squaredsubset %s --squaredsubsetfil %s --squaredsubsetfil_cpg2 %s --squaredcpg2 %s --age2cpg2 %s --logage %s --out %s" % tuple(args), shell  = True)
    else:
        args = ["methylation_training_noscale%s%s.rds" % (gsnorm_dic[gsnorm], external_dic_ii[external]), output_dir_prep, output_dir_folds + "gs_basic_folds_%s%s.tsv" % (ran, external_dic[external]), run_name, sex, squared, external, ext_select, totally_random, lasso, cpg, sample, subset, subset_f, squared_subset, squared_subset_f, squared_subset_f_cpg2, squared_cpg2, age2cpg2, logage, output_dir_train]
        subprocess.call("chronological_age/elnet_models.R --meth %s --methdir %s --pheno %s --name %s --sex %s --squared %s --external %s --ext_select %s --random %s --lasso %s --cpg %s --sample %s --subset %s --subfil %s --squaredsubset %s --squaredsubsetfil %s --squaredsubsetfil_cpg2 %s --squaredcpg2 %s --age2cpg2 %s --logage %s --out %s" % tuple(args), shell  = True)


# Elastic net model testing
######################################################

if runELNET_test == True:
    
    print("\n################### Model testing\n")
    """
    if logage20 == "T":
        logage20_old = "O"
        logage20_young = "Y"
        
        # Prediction for PREDICTED > 20 
        args = [testing_dir + "test_methylation.RDS", testing_dir + "test_samples.RDS", output_dir_train.replace(logage_dic[logage], ""), run_name.replace(logage_dic[logage], ""), output_dir_test, cpg_test, sample_test, coh, sex, squared, squared_subset, squared_subset_n, squared_cpg2, squared_subset_n_cpg2, age2cpg2, external, totally_random, lasso, logage, logage20_old, tf[loo], lbc, lbc_testing_samp, lbc_testing_meth] # Arguments for R script
        subprocess.call("chronological_age/elnet_test.R --meth %s --pheno %s --weights_dir %s --name %s --out %s --cpg %s --sample %s --test_cohort %s --sex %s --squared %s --squaredsubset %s --squaredsubset_n %s --squaredcpg2 %s --squaredsubset_n_cpg2 %s --age2cpg2 %s --external %s --random %s --lasso %s --logage %s --logage20 %s --loo %s --lbc %s --lbc_df %s --lbc_meth %s" % tuple(args), shell  = True)
        
        # Prediction for PREDICTED < 20
        args = [testing_dir + "test_methylation.RDS", testing_dir + "test_samples.RDS", output_dir_train, run_name, output_dir_test, cpg_test, sample_test, coh, sex, squared, squared_subset, squared_subset_n, squared_cpg2, squared_subset_n_cpg2, age2cpg2, external, totally_random, lasso, logage, logage20_young, tf[loo], lbc, lbc_testing_samp, lbc_testing_meth] # Arguments for R script
        subprocess.call("chronological_age/elnet_test.R --meth %s --pheno %s --weights_dir %s --name %s --out %s --cpg %s --sample %s --test_cohort %s --sex %s --squared %s --squaredsubset %s --squaredsubset_n %s --squaredcpg2 %s --squaredsubset_n_cpg2 %s --age2cpg2 %s --external %s --random %s --lasso %s --logage %s --logage20 %s --loo %s --lbc %s --lbc_df %s --lbc_meth %s" % tuple(args), shell  = True)

        # Fuse together
        pred_old = pd.read_table(output_dir_test + "predictions%s%s%s%s%s%s%s%s%s_logage20_O.tsv" % (sex_dic[sex], squared_dic[squared], squared_subset_dic[squared_subset], cpg2_dic[squared_cpg2], totallyrandom_dic[totally_random], lasso_dic[lasso], loo_dic[loo], lbc_dic[lbc], test_cohort_dic[coh_f]), index_col = 0)
        if exists(output_dir_test + "predictions%s%s%s%s%s%s%s%s%s_logage20_Y.tsv" % (sex_dic[sex], squared_dic[squared], squared_subset_dic[squared_subset], cpg2_dic[squared_cpg2], totallyrandom_dic[totally_random], lasso_dic[lasso], loo_dic[loo], lbc_dic[lbc], test_cohort_dic[coh_f])):
            pred_young = pd.read_table(output_dir_test + "predictions%s%s%s%s%s%s%s%s%s_logage20_Y.tsv" % (sex_dic[sex], squared_dic[squared], squared_subset_dic[squared_subset], cpg2_dic[squared_cpg2], totallyrandom_dic[totally_random], lasso_dic[lasso], loo_dic[loo], lbc_dic[lbc], test_cohort_dic[coh_f]), index_col = 0)
            pred_df = pd.concat([pred_young, pred_old], axis = 0)
        else:
            pred_df = pred_old
        pred_df.to_csv(output_dir_test + "predictions%s%s%s%s%s%s%s%s%s_logage20.tsv" % (sex_dic[sex], squared_dic[squared], squared_subset_dic[squared_subset], cpg2_dic[squared_cpg2], totallyrandom_dic[totally_random], lasso_dic[lasso], loo_dic[loo], lbc_dic[lbc], test_cohort_dic[coh_f]), sep = "\t")
    else:
        args = [testing_dir + "test_methylation.RDS", testing_dir + "test_samples.RDS", output_dir_train, run_name, output_dir_test, cpg_test, sample_test, coh, sex, squared, squared_subset, squared_subset_n, squared_cpg2, squared_subset_n_cpg2, age2cpg2, external, totally_random, lasso, logage, logage20, tf[loo], lbc, lbc_testing_samp, lbc_testing_meth] # Arguments for R script
        subprocess.call("chronological_age/elnet_test.R --meth %s --pheno %s --weights_dir %s --name %s --out %s --cpg %s --sample %s --test_cohort %s --sex %s --squared %s --squaredsubset %s --squaredsubset_n %s --squaredcpg2 %s --squaredsubset_n_cpg2 %s --age2cpg2 %s --external %s --random %s --lasso %s --logage %s --logage20 %s --loo %s --lbc %s --lbc_df %s --lbc_meth %s" % tuple(args), shell  = True)
    """
    # Now explore a bit
    explore_test(prediction_file = output_dir_test + "predictions%s%s%s%s%s%s%s%s%s%s.tsv" % (sex_dic[sex], squared_dic[squared], squared_subset_dic[squared_subset], cpg2_dic[squared_cpg2], totallyrandom_dic[totally_random], lasso_dic[lasso], loo_dic[loo], lbc_dic[lbc], test_cohort_dic[coh_f], logage20_dic[logage20]), pred_col = "Age_pred", true_col = "Age", dataset_col = "GEO", sex_col = "Sex", sex_stratified = sex, squared = squared, squared_subset = squared_subset, squared_subset_n = squared_subset_n, squared_cpg2 = squared_cpg2, squared_subset_n_cpg2 = squared_subset_n_cpg2, age2cpg2 = age2cpg2, external = external, random = totally_random, lasso = lasso, loo = loo, lbc = lbc, logage20 = logage20, output_dir = output_dir_test, palette = palette)


# Diagnostics
######################################################

if runELNET_test == True:

    print("\n################### Model diagnostics\n")
    explore_resid_combo(prediction_file = output_dir_test + "predictions%s%s%s%s%s%s%s%s%s%s.tsv" % (sex_dic[sex], squared_dic[squared], squared_subset_dic[squared_subset], cpg2_dic[squared_cpg2], totallyrandom_dic[totally_random], lasso_dic[lasso], loo_dic[loo], lbc_dic[lbc], test_cohort_dic[coh_f], logage20_dic[logage20]), pred_col = "Age_pred", true_col = "Age", dataset_col = "GEO", sex_col = "Sex", sex_stratified = sex, squared = squared, squared_subset = squared_subset, squared_subset_n = squared_subset_n, squared_cpg2 = squared_cpg2, squared_subset_n_cpg2 = squared_subset_n_cpg2, age2cpg2 = age2cpg2, external = external, random = totally_random, lasso = lasso, loo = loo, lbc = lbc, logage20 = logage20, output_dir = output_dir_test, palette = palette)
    residual_plot(pred_file = output_dir_test + "predictions%s%s%s%s%s%s%s%s%s%s.tsv" % (sex_dic[sex], squared_dic[squared], squared_subset_dic[squared_subset], cpg2_dic[squared_cpg2], totallyrandom_dic[totally_random], lasso_dic[lasso], loo_dic[loo], lbc_dic[lbc], test_cohort_dic[coh_f], logage20_dic[logage20]), output_dir = output_dir_test, name = "residualanalysis%s%s%s%s%s%s%s" % (sex_dic[sex], squared_dic[squared], squared_subset_dic[squared_subset], totallyrandom_dic[totally_random], lasso_dic[lasso], loo_dic[loo], lbc_dic[lbc]))


print("\nDone!")
total_time = time.time() - start
print("\nTotal time: %s minutes (%s hours)." % (round(total_time/60, 3), round((total_time/60)/60, 3)))