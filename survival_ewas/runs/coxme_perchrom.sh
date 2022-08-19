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
    nohup python survival_ewas.py --coxme --chrom $CHROM --cpgs /Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/survival_ewas/coxph/coxph_survival_cpgs_ews.txt > /Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/scripts/survival_ewas/logs/10062022_coxmesurvival_$CHROM.out & 
done
#)