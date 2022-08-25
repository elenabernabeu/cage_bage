# cAge predictor

## Requirements

To obtain a cAge prediction you need methylation data (in the form of M-values or beta-values), preferably rows being CpGs and columns being samples (but not required).

## Instructions

Download the scripts contained in this subdirectory and place them in your working directory. 

Once files are ready, you need to edit the top of the `cage_predictor.R` script and assign the paths to your methylation file, as well as its format (rds, tsv, or csv). 

Once that's done you can run the script and that will output predictions.

Files that will be outputted by script:
- `cage_predictions.tsv`: bAge predictions.
