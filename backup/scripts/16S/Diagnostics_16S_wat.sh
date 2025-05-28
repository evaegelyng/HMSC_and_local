#!/bin/bash
#SBATCH --account edna
#SBATCH --partition normal
#SBATCH --mem-per-cpu 4G
#SBATCH -c 8
#SBATCH -t 7-00:00:00

Rscript ../scripts/16S/Diagnostics_wat.r