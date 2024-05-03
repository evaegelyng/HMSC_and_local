#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 128G
#SBATCH -c 1
#SBATCH -t 7-00:00:00

Rscript scripts/COI/hmsc_COI_rich_water_240429.r
