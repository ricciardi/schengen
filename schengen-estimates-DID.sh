#!/bin/bash
#SBATCH --job-name="schengen-estimates-DID"
#SBATCH -e schengen-estimates-DID.err
#SBATCH -p scavenger
#SBATCH -c 6
#SBATCH --mem=30G

module load R/3.6.3

# Run 
R --no-save < Schengen_DID.R