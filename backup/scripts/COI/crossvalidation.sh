#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 128G
#SBATCH -c 4
#SBATCH -t 7-00:00:00

Rscript scripts/VP/crossvalidation.r
