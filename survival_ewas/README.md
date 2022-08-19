# Survival EWAS

Scripts detailing survival EWAS.

`survival_ewas.py` takes multiple arguments, which indicate whether to run EWAS using coxph, coxme, or using coxph ahead of OSCA. Scripts for the main analyses contained in `survival_ewas/` directory. 

~~~
usage: survival_ewas.py [-h] [--coxph_osca] [--coxph] [--coxme] [--cpgs CPGS] [--chrom {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}]

Survival EWAS pipeline.

optional arguments:
  -h, --help            show this help message and exit
  --coxph_osca          Option to run EWAS first fitting CoxPH with age and sex, then OSCA on Martingale residuals.
  --coxph               Option to run EWAS one CpG at a time using CoxPH.
  --coxme               Option to run EWAS one CpG at a time using CoxME (accounting for relatedness).
  --cpgs CPGS           File location of CpGs to subset analysis to.
  --chrom {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}
                        Chromosome to consider in analysis (for coxph and coxme options). If none chosen, all are run one after the other.
~~~


