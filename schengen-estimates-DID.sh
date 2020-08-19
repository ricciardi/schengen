#!/bin/bash
#SBATCH --job-name="schengen-estimates-DID"
#SBATCH -e schengen-estimates-DID.err
#SBATCH -p common
#SBATCH -c 12
#SBATCH --mem=30G

module load R/3.6.0

# Run 
R --no-save < Schengen_DID.R