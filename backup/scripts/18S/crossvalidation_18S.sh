#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 64G
#SBATCH -c 4
#SBATCH -t 7-00:00:00

Rscript scripts/18S/crossvalidation_18S.r
