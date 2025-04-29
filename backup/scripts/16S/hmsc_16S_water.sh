#!/bin/bash
#SBATCH --account edna
#SBATCH -c 2
#SBATCH --mem 128g
#SBATCH --partition gpu
#SBATCH --gres=gpu:1
#SBATCH --time 7-00:00:00
#SBATCH --array=0-3

SAM=250
THIN=1000
input_path="results/Water/16S/init_file.rds"
output_path=$(printf "results/Water/16S/post_file%.2d.rds" $SLURM_ARRAY_TASK_ID)

srun python3 -m hmsc.run_gibbs_sampler --input $input_path --output $output_path --samples $SAM --transient $(($SAM*$THIN)) --thin $THIN --verbose 100 --chains $SLURM_ARRAY_TASK_ID