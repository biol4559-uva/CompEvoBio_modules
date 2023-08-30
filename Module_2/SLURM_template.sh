#!/usr/bin/env bash

#SBATCH -J runFASTQC # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 2:00:00 ### 6 hours
#SBATCH --mem 10G
#SBATCH -o /scratch/aob2x/logs/runFASTQC.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/logs/runFASTQC.%A_%a.err # Standard error
#SBATCH -p instructional
#SBATCH --account biol4559-aob2x

### run as: sbatch /project/biol4559-aob2x/repos/CompEvoBio_modules/utils/makeSingularityDocker.sh
### sacct -j XXXXXXXXX
### cat /scratch/aob2x/logs/runFASTQC.*.err

module load fastqc
module load multiqc

proj=SRP002024
srr=SRR036932 #srr=SRR036934


if [ ! -d ~/fastq_QC_Out/ ]; then
  mkdir ~/fastq_QC_Out/
fi


### single call
fastqc -o ~/fastq_QC_Out/ --noextract /project/biol4559-aob2x/data/fastq/${proj}/${srr}.fastq.gz

### loop
for file in /project/biol4559-aob2x/data/fastq/${proj}/*.fastq.gz; do
  echo ${file}
  fastqc -o ~/fastq_QC_Out/ --noextract ${file}
done


if [ ! -d "/path/to/dir" ]; then
  mkdir ~/multiqc_Out/
fi

multiqc -d ~/fastq_QC_Out/ -o ~/multiqc_Out/
