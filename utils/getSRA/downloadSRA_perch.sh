#!/usr/bin/env bash
#
#SBATCH -J download_SRA # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6:00:00 ### 6 hours
#SBATCH --mem 10G
#SBATCH -o /scratch/aob2x/compBio/logs/prefetch.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/compBio/logs/prefetch.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

wd=/scratch/aob2x/euroPerch

### run as: sbatch --array=1-$( wc -l < ~/misc/1000G/__todo/EuroPerch/sras.csv )%20 ~/misc/1000G/__todo/EuroPerch/downloadSRA.sh
### sacct -j 9735165
### cat /scratch/aob2x/compBio/logs/prefetch.52222298_*.out | grep -B1 "do not"
### cat /scratch/aob2x/compBio/logs/prefetch.9735165_3.out

module load gcc/11.4.0 sratoolkit/3.1.1 aspera-connect/4.2.8

#SLURM_ARRAY_TASK_ID=2
# cat /home/aob2x/CompEvoBio_modules/data/runs.csv | nl | grep "SRR1988514"
# SLURM_ARRAY_TASK_ID=1

sranum=$( sed "${SLURM_ARRAY_TASK_ID}q;d" /home/aob2x/misc/1000G/__todo/EuroPerch/sras.csv | cut -f2 -d',' )
sampName=$( sed "${SLURM_ARRAY_TASK_ID}q;d" /home/aob2x/misc/1000G/__todo/EuroPerch/sras.csv | cut -f3 -d',' )

echo $sampName " / " $sranum

### sranum=SRR1184609; proj=PRJNA194129


if ls ${wd}/fastq/${sranum}*fastq.gz 1> /dev/null 2>&1; then
    echo "files do exist"
else
  echo "files do not exist"

  echo "force re-download"
  prefetch \
  -o ${wd}/sra/${sranum}.sra \
  -p \
  ${sranum}

  fasterq-dump \
  --split-files \
  --split-3 \
  --outfile ${wd}/fastq/${sranum} \
  -e 10 \
  -p \
  --temp /scratch/aob2x/tmp \
  ${wd}/sra/${sranum}.sra

  ls -lh ${wd}/fastq/${sranum}*

fi

if [ -f "/scratch/aob2x/compBio/fastq/${proj}/${sranum}_1.fastq" ]; then
  gzip ${wd}/fastq//${sranum}_1.fastq
  gzip ${wd}/fastq//${sranum}_2.fastq
fi

if [ -f "/scratch/aob2x/compBio/fastq/${proj}/${sranum}" ]; then
  gzip -c ${wd}/fastq/${sranum} > ${wd}/fastq/${sranum}.fastq.gz
  rm /scratch/aob2x/compBio/fastq/${proj}/${sranum}
fi

#rm /scratch/aob2x/fastq/${sranum}.sra
