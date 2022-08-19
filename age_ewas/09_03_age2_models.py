#!/usr/bin/env python
# -*- coding: utf-8 -*-
# By Elena

import pandas as pd
import pyreadr
import numpy as np

output_dir = "/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/age2_lm/w1w3w4/"
file_base = "age2_lm_chrom%s.tsv"

df = pd.DataFrame()
for chrom in range(1,23):
    df_chrom = pd.read_table(output_dir + file_base % chrom, index_col = 0)
    if chrom == 1:
        df = df_chrom
    else:
        df = pd.concat([df, df_chrom])

# Filter to QC´d CpGs
# cpgs = [line.strip() for line in open("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/data_prep/w1w3w4/gs20k_cpgs_common.txt", 'r')]
cpgs = [i.rstrip() for i in open("/Cluster_Filespace/Marioni_Group/Elena/gs_osca/data/cpgs_tokeep.txt", "r").readlines()]
df = df[df.index.isin(cpgs)]

# Add in location and chromosome info from EPIC
anno = pyreadr.read_r("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds")
anno = anno[None]
anno = anno.loc[anno.index.isin(df.index),["chr","pos","UCSC_RefGene_Name"]]
anno["chr"] = [i[3:] for i in anno["chr"].values.tolist()]
df = pd.concat([anno, df], axis = 1)

# Sort by p-val
df["cpg_log10p"] = -np.log10(df["cpg_p"])
df["cpg2_log10p"] = -np.log10(df["cpg2_p"])
df = df.sort_values(by = ["cpg2_log10p", "cpg2_beta"], ascending = [False, False])
df = df.fillna("NA")

# Export
df.to_csv(output_dir + "age2_lm_allchrom.tsv", sep = "\t", index_label = "cpg")

