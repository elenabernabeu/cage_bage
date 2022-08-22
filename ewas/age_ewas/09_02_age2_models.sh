#!/usr/bin/env bash
# -*- coding: utf-8 -*-
# By Elena

export PATH=/home/ebernab3/anaconda3/bin:$PATH
cd /Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/scripts/age_ewas/w1w3w4/

#N=10
#(
for CHROM in {1..22}
do
#    ((i=i%N)); ((i++==0)) && wait
    nohup Rscript 09_01_age2_models.R --chrom $CHROM > /Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/age2_lm/w1w3w4/logs/age2_chrom$CHROM.out & 
done
#)