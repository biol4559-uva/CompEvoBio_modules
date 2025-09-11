#!/usr/bin/env bash

#SBATCH -J runFASTQC # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0:00:15 ### 15 seconds
#SBATCH --mem 10G
#SBATCH -o /scratch/COMPUTE_ID/logs/demo_1.%A_%a.out # Standard output
#SBATCH -e /scratch/COMPUTE_ID/logs/demo_1.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account biol4559-aob2x

### run as: sbatch --array=1-5 PATH_TO_THIS_FILE
### sacct -j XXXXXXXXX
### cat /scratch/COMPUTE_ID/logs/demo_1.*.err

echo "this is job: "
echo ${SLURM_ARRAY_TASK_ID}

bioproj="SRP002024"
# SLURM_ARRAY_TASK_ID=1
file=$( ls -d /standard/BerglandTeach/data/fastq/${bioproj}/* | tr '\t' '\n' | sed "${SLURM_ARRAY_TASK_ID}q;d" )
echo $file
