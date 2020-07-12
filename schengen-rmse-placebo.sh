#!/bin/bash
#SBATCH --job-name="schengen-placebo"
#SBATCH -e schengen-placebo.err
#SBATCH -p common
#SBATCH -c 6
#SBATCH --mem=30G

module load R/3.6.0 

# Run 
R --no-save < schengen-rmse-placebo.R