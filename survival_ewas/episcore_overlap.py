#!/usr/bin/env python
# -*- coding: utf-8 -*-
# By Elena

import os
import numpy as np
import pandas as pd
import random
import subprocess
import argparse
import time
from datetime import datetime

# For VSCode + cluster: export PATH=/home/ebernab3/anaconda3/bin:$PATH

episcore_df = pd.read_csv("/Cluster_Filespace/Marioni_Group/Elena/data/episcores/protein_episcore_cpgs.csv", sep = ";")
mortality_ewas = pd.read_table("/Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/survival_ewas/coxph/coxph_survival_ewas_ews.tsv", index_col = 0) # 796

# CpGs in episcores
cpgs = set(episcore_df["CpG Site"]) # 9101
genes = set(episcore_df["Identifier"]) # 109

# Sig CpGs from mortality ewas in episcore CpGs
overlap = [cpg for cpg in mortality_ewas.index.values.tolist() if cpg in cpgs] # 148