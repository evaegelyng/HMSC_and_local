#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 256G
#SBATCH -c 1
#SBATCH -t 7-00:00:00

Rscript scripts/16S/hmsc_16S_rich_water.r