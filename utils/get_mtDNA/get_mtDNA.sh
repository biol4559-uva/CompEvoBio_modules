#!/usr/bin/env bash
#
#SBATCH -J getMtDNA # A single job name for the array
#SBATCH -c 20 ### 20 cores
#SBATCH -N 1 # on one node
#SBATCH -t 1:00:00
#SBATCH --mem 10G
#SBATCH -o /scratch/aob2x/logs/RunDest.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/logs/RunDest.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### ijob -A berglandlab -c20 -p standard --mem=40G
### run as: sbatch --array=1-816 /home/aob2x/CompEvoBio_modules/utils/get_mtDNA/get_mtDNA.sh
### sacct -j 4918194
### cat /scratch/aob2x/logs/RunDest.4841903*.err

### modules
  module load samtools
  module load bcftools
  module load gcc/11.4.0 openmpi/4.1.4 python/3.11.4

### get bam file list
  # SLURM_ARRAY_TASK_ID=2
  bam=$( ls -d /project/berglandlab/DEST/dest_mapped/*/*/*.mel.bam | tr ' ' '\n' | sed "${SLURM_ARRAY_TASK_ID}q;d" )
  echo $bam

  sample=$( echo ${bam} | rev | cut -d'/' -f1 | rev | sed 's/bam//g' )

### get reads mapping to mitochondria

  #inputFile=/project/berglandlab/DEST/dest_mapped/Cville/US_Vir_Cha_1_2018-11-29/US_Vir_Cha_1_2018-11-29.original.bam
  #samtools idxstats ${inputFile} | grep -E "mitochondrion_genome" | cut -f1,2 | awk '{print $1"\t"1"\t"$2}' > /standard/BerglandTeach/mtDNA/mtDNA.bed
  #sed -i '$d' /scratch/aob2x/DESTv2_unmapped_reads/mitoGenome.bed

  samtools view -@ 20 -L /standard/BerglandTeach/mtDNA/mtDNA.bed ${bam} -b > \
  /standard/BerglandTeach/mtDNA/bams/${sample}mt.bam

  samtools mpileup \
  --reference /scratch/aob2x/tmpRef/holo_dmel_6.12.fa \
  /standard/BerglandTeach/mtDNA/bams/${sample}mt.bam > \
  /standard/BerglandTeach/mtDNA/pileup/${sample}mt.pileup

  python ~/DESTv3/mappingPipeline/scripts/Mpileup2Sync.py \
  --mpileup /standard/BerglandTeach/mtDNA/pileup/${sample}mt.pileup \
  --base-quality-threshold 25 \
  --minIndel 5 \
  --coding 1.8 \
  --ref /scratch/aob2x/tmpRef/holo_dmel_6.12.mel.fa.ref \
  --output /standard/BerglandTeach/mtDNA/sync/${sample}mt.sync
