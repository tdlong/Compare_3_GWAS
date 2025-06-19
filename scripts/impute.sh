#!/bin/bash
#SBATCH --job-name=impute
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=8

module load R/3.6.2
Rscript scripts/impute.R $1 $2 $3 $4

