#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 4G
#SBATCH -c 16
#SBATCH -t 7-00:00:00

Rscript scripts/18S/crossvalidation_18S.r
