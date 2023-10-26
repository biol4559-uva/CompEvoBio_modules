#!/usr/bin/env bash
#
#SBATCH -J fst # A single job name for the array
#SBATCH -c 10 ### 10 cores
#SBATCH -N 1 # on one node
#SBATCH -t 0:15:00
#SBATCH --mem 5G
#SBATCH -o /scratch/aob2x/logs/fst.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/logs/fst.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account biol4559-aob2x

### run as: sbatch --array=1-NUMBER_OF_FILES PATH_TO_THIS_FILE
### sacct -j XXXXXXXXX
### cat /scratch/COMPUTE_ID/logs/fst.*.err


### modules
  module load intel/18.0 intelmpi/18.0 R/4.0.3

### run window

  Rscript --vanilla FULL_PATH_TO-outFLANK_Fst.R ${SLURM_ARRAY_TASK_ID}
