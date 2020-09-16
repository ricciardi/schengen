#!/bin/bash
#SBATCH --job-name="schengen-data"
#SBATCH -e schengen-data.err
#SBATCH -p common
#SBATCH -c 6
#SBATCH --mem=25G

module load R/3.6.3

# Run 
R --no-save < schengen_MCM_data.R