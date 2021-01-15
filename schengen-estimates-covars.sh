#!/bin/bash
#SBATCH --job-name="schengen-estimates-covars"
#SBATCH -e schengen-estimates-covars.err
#SBATCH -p scavenger-gpu --gres=gpu:1
#SBATCH -c 6
#SBATCH --mem=30G

module load R/3.6.3

# Run 
R --no-save < Schengen_MCM_covars.R