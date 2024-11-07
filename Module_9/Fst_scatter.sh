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

### run as: sbatch --array=1-NUMBER_OF_WINDOWS PATH_TO_THIS_FILE
### sacct -j XXXXXXXXX
### cat /scratch/COMPUTE_ID/logs/fst.*.err


### modules
  module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1; echo "R_LIBS_USER=~/R/goolf/4.3" > ~/.Renviron

### run window

  Rscript --vanilla FULL_PATH_TO-Fst_template.R ${SLURM_ARRAY_TASK_ID}
