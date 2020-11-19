#!/bin/bash
#SBATCH --job-name="schengen-placebo"
#SBATCH -e schengen-placebo.err
#SBATCH -p scavenger
#SBATCH -c 6
#SBATCH --mem=30G

module load R/3.6.3

# Run 
R --no-save < mc-placebo.R