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

coxph_sig = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/survival_ewas/coxph/coxph_survival_ewas_ews.tsv", index_col = 0)

# Get unique genes for FUMA
unique_genes = [i for i in set(coxph_sig["Gene"]) if not pd.isna(i)]
unique_genes = [gene.split(";")[0] for gene in unique_genes] # Keep only first gene
unique_genes = list(set(unique_genes)) # And filter out after that

new_file = open("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/survival_ewas/coxph/unique_genes.txt", "w")
new_file.write("\n".join(unique_genes))
new_file.close()

