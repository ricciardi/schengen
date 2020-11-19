#!/bin/bash
#SBATCH --job-name="schengen-estimates"
#SBATCH -e schengen-estimates.err
#SBATCH -p scavenger
#SBATCH -c 6
#SBATCH --mem=30G

module load R/3.6.3

# Run 
R --no-save < Schengen_MCM.R