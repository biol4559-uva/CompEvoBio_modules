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
#SBATCH --account biol4020-aob2x

### run as: sbatch --array=1-1085 ~/CompEvoBio_modules/Module_9/Fst_scatter.sh
### sacct -j 22497591
### cat /scratch/COMPUTE_ID/logs/fst.*.err


### modules
  module load gcc/14.2.0 openmpi/5.0.7 R/4.6.0

### run window

  Rscript --vanilla ~/CompEvoBio_modules/Module_9/2025_fstwindow.R ${SLURM_ARRAY_TASK_ID}
