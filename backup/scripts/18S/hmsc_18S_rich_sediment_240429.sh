#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 128G
#SBATCH -c 1
#SBATCH -t 7-00:00:00

Rscript scripts/18S/hmsc_18S_rich_sediment_240429.r