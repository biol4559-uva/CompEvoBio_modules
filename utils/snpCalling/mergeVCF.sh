#!/bin/bash
#
#SBATCH -J manual_annotate # A single job name for the array
#SBATCH --ntasks-per-node=20 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 20:00:00 ### 1 hours
#SBATCH --mem 80G
#SBATCH -o /scratch/aob2x/compBio_SNP_28Sept2024/manual_annotate.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/compBio_SNP_28Sept2024/manual_annotate.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account biol4020-aob2x

### cat /scratch/aob2x/DESTv2_output_SNAPE/logs/runSnakemake.49369837*.err

### sbatch ~/CompEvoBio_modules/utils/snpCalling/mergeVCF.sh
### sacct -j 64525372
### cat /scratch/aob2x/compBio_SNP_28Sept2024/manual_annotate.64525372

module load gcc/11.4.0  openmpi/4.1.4 python/3.11.4

module load htslib/1.17  bcftools/1.17 parallel/20200322 gcc/11.4.0 openmpi/4.1.4 python/3.11.4 vcftools/0.1.16 R/4.3.1
module load bedtools/2.30.0


### annotate
popSet=PoolSeq
method=PoolSNP
maf=001
mac=50
version=28Sept2024_ExpEvo
wd=/scratch/aob2x/compBio_SNP_28Sept2024

# bcftools view \
# --threads 20 \
# ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz | \
# java -jar ~/snpEff/snpEff.jar \
# eff \
# BDGP6.86 - > \
# ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf
#
#
# echo "bgzip & tabix"
#   bgzip -@20 -c ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf > ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf.gz
#   tabix -f -p vcf ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf.gz
#

bcftools merge -0 --threads 20 \
-o /scratch/aob2x/compBio_SNP_25Sept2023/dest.expevo.PoolSeq.PoolSNP.001.50.28Sept2024_ExpEvo.norep.vcf.gz
${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf.gz \
/project/berglandlab/DEST/dest.all.PoolSNP.001.50.3May2024.ann.vcf.gz



Rscript --vanilla ~/CompEvoBio_modules/utils/snpCalling/scatter_gather_annotate/vcf2gds.R \
/scratch/aob2x/compBio_SNP_25Sept2023/dest.expevo.PoolSeq.PoolSNP.001.50.28Sept2024_ExpEvo.norep.vcf.gz
