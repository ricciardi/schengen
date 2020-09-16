#!/bin/bash
#SBATCH --job-name="schengen-rmse-placebo"
#SBATCH -e schengen-rmse-placebo.err
#SBATCH -p common
#SBATCH -c 6
#SBATCH --mem=30G

module load R/3.6.3

# Run 
R --no-save < schengen-rmse-placebo.R