#!/bin/bash
#
#SBATCH -J manual_annotate # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 2:00:00 ### 1 hours
#SBATCH --mem 20G
#SBATCH -o /scratch/aob2x/compBio_SNP_25Sept2023/logs/manual_annotate.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/compBio_SNP_25Sept2023/logs/manual_annotate.%A_%a.err # Standard error
#SBATCH -p instructional
#SBATCH --account biol4559-aob2x

### cat /scratch/aob2x/DESTv2_output_SNAPE/logs/runSnakemake.49369837*.err

### sbatch /scratch/aob2x/CompEvoBio_modules/utils/snpCalling/mergeVCF.sh
### sacct -j 54019500
### cat /scratch/aob2x/DESTv2_output_26April2023/logs/manual_annotate.49572492*.out

module purge
module load  htslib/1.10.2 bcftools/1.9 intel/18.0 intelmpi/18.0 parallel/20200322 R/3.6.3 samtools vcftools

module load gcc/9.2.0 bedtools/2.29.2

bcftools merge -0 --threads 1 \
-o /scratch/aob2x/compBio_SNP_25Sept2023/dest.expevo.PoolSNP.001.50.11Oct2023.norep.ann.vcf.gz \
/scratch/aob2x/compBio_SNP_25Sept2023/dest.PoolSeq.PoolSNP.001.50.25Sept2023.norep.vcf.gz \
/scratch/aob2x/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.vcf.gz
