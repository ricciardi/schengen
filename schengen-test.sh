#!/bin/bash
#SBATCH --job-name="schengen-test"
#SBATCH -e schengen-test.err
#SBATCH -p common
#SBATCH -c 12
#SBATCH --mem=30G

module load R/3.6.0

# Run 
R --no-save < test_MCPanel.R