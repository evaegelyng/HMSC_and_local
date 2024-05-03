#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 12G
#SBATCH -c 4
#SBATCH -t 4-00:00:00

Rscript scripts/VP/crossvalidation_18S.r
