#!/usr/bin/env bash
#
#SBATCH -J index # A single job name for the array
#SBATCH -c 10 ### 10 cores
#SBATCH -N 1 # on one node
#SBATCH -t 00:10:00
#SBATCH --mem 10G
#SBATCH -o /scratch/aob2x/logs/samtoolsindex.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/logs/samtoolsindex.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account biol4020-aob2x

### sbatch --array=1-$( cat /scratch/aob2x/small_file.csv | wc -l ) /scratch/aob2x/indexBam.sh

### modules
module load samtools

###################################
# Part  1. Get Sample information #
###################################
#SLURM_ARRAY_TASK_ID=1

sampleId=$( cat /scratch/aob2x/small_file.csv | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f3 -d',' )
srr=$( cat /scratch/aob2x/small_file.csv | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f1 -d',' )
numFlies=$( cat /scratch/aob2x/small_file.csv | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f2 -d',' )


echo ${sampleId}
echo ${srr}
echo ${numFlies}

samtools index -@ 10 /standard/BerglandTeach/mapping_output/${sampleId}/${sampleId}.original.bam
