#!/usr/bin/env bash

#SBATCH -J runFASTQC # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 2:00:00 ### 6 hours
#SBATCH --mem 10G
#SBATCH -o /scratch/COMPUTE_ID/logs/runFASTQC.%A_%a.out # Standard output
#SBATCH -e /scratch/COMPUTE_ID/logs/runFASTQC.%A_%a.err # Standard error
#SBATCH -p instructional
#SBATCH --account biol4559-aob2x

### run as: sbatch --array=1-NUMBER_OF_FILES PATH_TO_THIS_FILE
### sacct -j XXXXXXXXX
### cat /scratch/COMPUTE_ID/logs/runFASTQC.*.err

module load fastqc

proj=SRP002024 ### YOU'LL NEED TO REPLACE THIS with your file
file=$( ls -d /project/biol4559-aob2x/data/fastq/${bioproj}/* | tr '\t' '\n' | sed "${SLURM_ARRAY_TASK_ID}q;d" )

if [ ! -d /scratch/COMPUTEID/fastq_QC_Out/ ]; then
  mkdir   /scratch/COMPUTEID/fastq_QC_Out/
fi

echo ${file}
fastqc -o /scratch/COMPUTEID/fastq_QC_Out/ --noextract ${file}
