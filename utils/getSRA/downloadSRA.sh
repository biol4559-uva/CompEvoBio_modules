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

wd=/scratch/aob2x/compBio
### run as: sbatch --array=1-$( wc -l < ~/CompEvoBio_modules/data/runs.csv )%10 ~/CompEvoBio_modules/utils/getSRA/downloadSRA.sh
### sacct -j 52141033
### cat /scratch/aob2x/compBio/logs/prefetch.52096862_1.out

module load sratoolkit/2.10.5

#SLURM_ARRAY_TASK_ID=1

sranum=$( sed "${SLURM_ARRAY_TASK_ID}q;d" /home/aob2x/CompEvoBio_modules/data/runs.csv | cut -f2 -d' ' )
sampName=$( sed "${SLURM_ARRAY_TASK_ID}q;d" /home/aob2x/CompEvoBio_modules/data/runs.csv | cut -f2 -d' ' )
proj=$( sed "${SLURM_ARRAY_TASK_ID}q;d" /home/aob2x/CompEvoBio_modules/data/runs.csv | cut -f1 -d' ' )

echo $sampName " / " $sranum

prefetch \
-o /scratch/aob2x/compBio/sra/${sranum}.sra \
-p \
${sranum}

### sranum=SRR1184609; proj=PRJNA194129

if [ ! -d "/scratch/aob2x/compBio/fastq/${proj}" ]; then
  mkdir /scratch/aob2x/compBio/fastq/${proj}
fi


if ls /scratch/aob2x/compBio/fastq/${proj}/${sranum}*fastq.gz 1> /dev/null 2>&1; then
    echo "files do exist"
else
  fasterq-dump \
  --split-files \
  --split-3 \
  --outfile /scratch/aob2x/compBio/fastq/${proj}/${sranum} \
  -e 10 \
  -p \
  /scratch/aob2x/compBio/sra/${sranum}.sra
fi

if [ -f "/scratch/aob2x/compBio/fastq/${proj}/${sranum}_1.fastq" ]; then
  gzip /scratch/aob2x/compBio/fastq/${proj}/${sranum}_1.fastq
  gzip /scratch/aob2x/compBio/fastq/${proj}/${sranum}_2.fastq
fi

if [ -f "/scratch/aob2x/compBio/fastq/${proj}/${sranum}" ]; then
  gzip -c /scratch/aob2x/compBio/fastq/${proj}/${sranum} > /scratch/aob2x/compBio/fastq/${proj}/${sranum}.fastq.gz
  rm /scratch/aob2x/compBio/fastq/${proj}/${sranum}
fi

#rm /scratch/aob2x/fastq/${sranum}.sra
