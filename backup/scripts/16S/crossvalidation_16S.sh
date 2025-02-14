#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 64G
#SBATCH -c 16
#SBATCH -t 7-00:00:00

Rscript scripts/16S/crossvalidation_16S.r
