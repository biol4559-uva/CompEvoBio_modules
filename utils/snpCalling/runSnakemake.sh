#!/bin/bash
#
#SBATCH -J runSnakemake # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 80:00:00 ### 1 hours
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/DEST/logs/runSnakemake.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/DEST/logs/runSnakemake.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

#module load gcc/9.2.0 openmpi/3.1.6 python/3.7.7 snakemake/6.0.5
module load gcc/11.4.0  openmpi/4.1.4 python/3.11.4 snakemake/9.8.1

cd ~/CompEvoBio_modules/utils/snpCalling
snakemake --profile ~/CompEvoBio_modules/utils/snpCallingslurm
