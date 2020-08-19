#!/bin/bash
#SBATCH --job-name="schengen-estimates-covars"
#SBATCH -e schengen-estimates-covars.err
#SBATCH -p common
#SBATCH -c 12
#SBATCH --mem=30G

module load R/3.6.0

# Run 
R --no-save < Schengen_MCM_covars.R