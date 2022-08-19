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

# Script to make subsets from mortality EWAS results
coxph = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/survival_ewas/coxph/coxph_survival_ewas.tsv", index_col = 0)

# Filter to those in all datasets considered in study
cpgs = [i.rstrip() for i in open("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_cpgs_common.txt", "r").readlines()]
coxph = coxph.loc[coxph.index.isin(cpgs),]

# Subset
subsets = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1250, 1500, 1750, 2000, 3000, 4000, 5000, 7500, 10000]

for subset in subsets:
    df = coxph.loc[coxph.index[0:subset],]
    df.to_csv("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/survival_ewas/coxph/subsets/coxph_survival_ewas_%s.tsv" % subset, index_label = "cpg", sep = "\t", na_rep = "NA")