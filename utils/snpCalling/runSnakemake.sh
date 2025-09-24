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

### sbatch ~/CompEvoBio_modules/utils/snpCalling/runSnakemake.sh
### sacct -j 3969679
# cluster: "sbatch --account berglandlab -p standard -o {resources.log_dir}/{rule}.%A.out -e {resources.log_dir}/{rule}.%A.err -J {rule} --ntasks-per-node={resources.ntasks_per_node} -N 1 -t {resources.time_limit} --mem {resources.memory_limit}G --parsable"
### cat /scratch/aob2x/DEST/logs/runSnakemake/*3959458*.err

#module load gcc/9.2.0 openmpi/3.1.6 python/3.7.7 snakemake/6.0.5
#module load gcc/11.4.0  openmpi/4.1.4 python/3.11.4 snakemake/9.8.1

module load miniforge/24.11.3-py3.12
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake=7.24.2
conda activate snakemake


cd ~/CompEvoBio_modules/utils/snpCalling
snakemake --profile ~/CompEvoBio_modules/utils/snpCalling/slurm
