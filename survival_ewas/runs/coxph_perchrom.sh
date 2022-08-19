#!/usr/bin/env bash
# -*- coding: utf-8 -*-
# By Elena

export PATH=/home/ebernab3/anaconda3/bin:$PATH
cd /Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/scripts/survival_ewas/

#N=10
#(
for CHROM in {1..22}
do
#    ((i=i%N)); ((i++==0)) && wait
    nohup python survival_ewas.py --coxph --chrom $CHROM > /Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/scripts/survival_ewas/logs/06062022_coxphsurvival_$CHROM.out & 
done
#)