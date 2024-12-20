#!/bin/bash
#
#SBATCH -J runSnakemake # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 80:00:00 ### 1 hours
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/compBio_SNP_28Sept2024/logs/runSnakemake.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/compBio_SNP_28Sept2024/logs/runSnakemake.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account biol4020-aob2x

### cat /scratch/aob2x/DESTv2_output_26April2023/logs/runSnakemake.53518248*.err



#module load gcc/9.2.0 openmpi/3.1.6 python/3.7.7 snakemake/6.0.5
module load gcc/11.4.0  openmpi/4.1.4 python/3.11.4 snakemake/7.24.2

cd ~/CompEvoBio_modules/utils/snpCalling/
snakemake --profile ~/CompEvoBio_modules/utils/snpCalling/slurm --ri
