#!/usr/bin/env bash
#
#SBATCH -J vcf2gds # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0:04:00  ### 48 hours
#SBATCH --mem 24G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/vcf2gds.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/vcf2gds.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

### sbatch /scratch/aob2x/DESTv2/snpCalling/scatter_gather_annotate/vcf2gds.sh
### sacct -j 50241880
### cat /scratch/aob2x/dest/slurmOutput/vcf2gds.22867938

module load htslib/1.17 bcftools/1.17 parallel/20200322 gcc/11.4.0 openmpi/4.1.4 R/4.3.1 samtools vcftools

Rscript --vanilla ~/CompEvoBio_modules/utils/snpCalling/scatter_gather_annotate/vcf2gds.R \
/scratch/aob2x/compBio_SNP_29Sept2025/dest.PoolSeq.PoolSNP.001.50.29Sept2025_ExpEvo.norep.ann.vcf.gz
