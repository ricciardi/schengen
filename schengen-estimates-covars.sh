#!/bin/bash
#SBATCH --job-name="schengen-estimates-covars"
#SBATCH -e schengen-estimates-covars.err
#SBATCH -p scavenger
#SBATCH -c 24
#SBATCH --mem=50G

module load R/3.6.3

# Run 
R --no-save < Schengen_MCM_covars.R