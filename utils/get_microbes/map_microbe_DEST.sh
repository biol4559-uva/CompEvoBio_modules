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

# sbatch --array=1-12 ~/CompEvoBio_modules/utils/get_microbes/map_microbe_DEST.sh
# sacct -j 22499092

module load apptainer/1.4.5
#apptainer build -F /scratch/aob2x/dest_v3.sif docker://alanbergland/dest_v3:latest


###################################
# Part  1. Get Sample information #
###################################
  #SLURM_ARRAY_TASK_ID=1


    #tr -cd '\11\12\40-\176' < ~/CompEvoBio_modules/utils/get_microbes/small_file.a_pomoroum.csv | sed 's/ _/_/g' | sed -e '$a\' >  /scratch/aob2x/small_file.a_pomoroum.clean.csv

  sampleName=$( cat /scratch/aob2x/small_file.a_pomoroum.clean.csv | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f1 -d',' )
  srr=$(        cat /scratch/aob2x/small_file.a_pomoroum.clean.csv | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f3 -d',' )
  numFlies=$(   cat /scratch/aob2x/small_file.a_pomoroum.clean.csv | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f2 -d',' )


  echo ${sampleName}
  echo ${srr}
  echo ${numFlies}


  # grep ">" /project/berglandlab/Dmel_genomic_resources/References/DESTv3_dmelholo/holo.dmel_6.54.dsim_3.1.dest3.fa | grep "pomorum" | sed 's/>//g' | sed 's/^/A_pomorum,/g' | sed 's/$/,A/g'  > /scratch/aob2x/apomoroum.focalFile

    singularity run \
    /project/berglandlab/multispecies_endemism/dest_v3.sif  \
    --fastq_1 /standard/BerglandTeach/misc/data/fastq/new/PRJEB5713/${srr}_1.fastq.gz \
    --fastq_2 /standard/BerglandTeach/misc/data/fastq/new/PRJEB5713/${srr}_2.fastq.gz \
    --sample_name ${sampleName} \
    --output_directory /scratch/aob2x/martins_microbes/ \
    --threads 20 \
    --max-cov 0.95 \
    --min-cov 4 \
    --base-quality-threshold 25 \
    --num-flies ${numFlies} \
    --reference_genome /project/berglandlab/Dmel_genomic_resources/References/DESTv3_dmelholo/holo.dmel_6.54.dsim_3.1.dest3.fa \
    --focal_file /scratch/aob2x/apomoroum.focalFile \
    --prep_reference 0 \
    --do_map 1 \
    --do_pileup 1 \
    --do_poolsnp 1 \
    --do_snape 0 \
    --do_cleanup 0


