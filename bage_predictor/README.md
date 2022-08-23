# bAge predictor

## Requirements
To test bAge you need:
1. Methylation data (in the form of M-values or beta-values), preferably rows being CpGs and columns being samples.
2. Phenotype data: chronological age, sex (coded as 1 = females, 0 = males), time-to-event (death), and death status (1 = dead, 0 = alive).
3. GrimAge results: these must be obtained from http://dnamage.genetics.ucla.edu/. Guidance can be found here: http://dnamage.genetics.ucla.edu/important-hints. Some further notes:
    - You need to make an account
    - You will need to upload your DNAm file as well as a phenotype file with basic info like sample, age, sex, etc.
    - It has a file size limit so if data is large you may need to submit in batches
    - Results will be e-mailed to you

## Instructions

Download the scripts contained in this subdirectory and place them in your working directory. 

Once files are ready, you need to edit the top of the `bage_external_testing.R` script and assign the paths to your methylation, phenotype, and GrimAge files as well as their formats (rds, tsv, or csv). You also need to specify the column names used for age, sex, TTE, and death status. 

Once that's done you can run the script and that will output predictions as well as evaluate bAge's predictive ability and compare it to GrimAge.

Files that will be outputted by script:
- `projections_episcores.tsv`: Episcore projections for samples included in phenotype file/for which methylation data was provided. 
- `predictions.tsv`: bAge predictions.
- `coxph_testing.tsv`: Cox-PH model results assessing predictive ability of bAge. 



