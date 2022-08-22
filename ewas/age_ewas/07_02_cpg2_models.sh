#!/usr/bin/env bash
# -*- coding: utf-8 -*-
# By Elena

export PATH=/home/ebernab3/anaconda3/bin:$PATH
cd /Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/scripts/age_ewas/w1w3w4/

#N=15
#(
for CHROM in {1,2,3,5,6,7,8}
do
#    ((i=i%N)); ((i++==0)) && wait
    nohup Rscript 07_01_cpg2_models.R --chrom $CHROM --linear T > /Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/cpg2_lm/w1w3w4/logs/cpg_chrom_scaled_$CHROM.out & 
done
#)

#N=15
#(
#for CHROM in {1..22}
#do
#    ((i=i%N)); ((i++==0)) && wait
#    nohup Rscript 07_01_cpg2_models.R --chrom $CHROM --linear F > /Cluster_Filespace/Marioni_Group/Elena/epigenetic_clocks/age_ewas/cpg2_lm/w1w3w4/logs/cpg2_chrom_scaled_$CHROM.out & 
#done
#)