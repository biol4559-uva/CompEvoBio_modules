#!/usr/bin/env bash
#
#SBATCH -J dockerMap # A single job name for the array
#SBATCH -c 11 ### 11 cores
#SBATCH -N 1 # on one node
#SBATCH -t 72:00:00
#SBATCH --mem 90G
#SBATCH -o /scratch/aob2x/logs/RunDest.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/logs/RunDest.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account biol4020-aob2x

### modules
  module load gcc/11.4.0 sratoolkit/3.1.1 aspera-connect/4.2.8
  module load apptainer

###################################
# Part  1. Get Sample information #
###################################
  #SLURM_ARRAY_TASK_ID=1

  sampleId=$( cat $4 | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f3 -d',' )
  srr=$( cat $4 | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f1 -d',' )
  numFlies=$( cat $4 | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f2 -d',' )
  proj=$( cat $4 | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f3 -d',' | cut -f2 -d'_' )

  echo ${sampleId}
  echo ${srr}
  echo ${numFlies}
  echo ${proj}


#####
# download missing data

  if [ ! -d "/scratch/aob2x/compBio/fastq/${proj}" ]; then
    mkdir ${2}/${proj}
  fi


  if ls /scratch/aob2x/compBio/fastq/${proj}/${sranum}*fastq.gz 1> /dev/null 2>&1; then
      echo "files do exist"
  else
    echo "files do not exist"

    echo "force re-download"
    prefetch \
    -o /scratch/aob2x/compBio/sra/${sranum}.sra \
    -p \
    ${sranum}

    fasterq-dump \
    --split-files \
    --split-3 \
    --outfile ${2}/${proj}/${sranum} \
    -e 10 \
    -p \
    --temp /scratch/aob2x/tmp \
    /scratch/aob2x/compBio/sra/${sranum}.sra

    ls -lh /scratch/aob2x/compBio/fastq/${proj}/${sranum}*

  fi

  if [ -f "${2}/${proj}/${sranum}_1.fastq" ]; then
    gzip ${2}/${proj}/${sranum}_1.fastq
    gzip ${2}/${proj}/${sranum}_2.fastq
  fi

  if [ -f "/scratch/aob2x/compBio/fastq/${proj}/${sranum}" ]; then
    gzip -c ${2}/${proj}/${sranum} > ${2}/${proj}/${sranum}.fastq.gz
    rm ${2}/${proj}/${sranum}
  fi

  #rm /scratch/aob2x/fastq/${sranum}.sra
  #cat /home/aob2x/CompEvoBio_modules/data/runs.csv | nl | grep "SRR12463313"



###################################
# Part  2. Run the Container      #
###################################
### If your reads are Paired End use this version (and delete the other)

  singularity run \
  $1 \
  $2/${proj}/${srr}_1.fastq.gz \
  $2/${proj}/${srr}_2.fastq.gz \
  ${sampleId} \
  $3 \
  --cores $SLURM_CPUS_PER_TASK \
  --max-cov 0.95 \
  --min-cov 4 \
  --base-quality-threshold 25 \
  --num-flies ${numFlies} \
  --do_poolsnp
