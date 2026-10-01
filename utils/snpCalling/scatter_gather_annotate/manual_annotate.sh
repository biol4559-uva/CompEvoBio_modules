#!/bin/bash
#
#SBATCH -J manual_annotate # A single job name for the array
#SBATCH --ntasks-per-node=48 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 12:00:00 ### 1 hours
#SBATCH --mem 40G
#SBATCH -o /scratch/aob2x/29Sept2025_ExpEvo/logs/manual_annotate.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/29Sept2025_ExpEvo/logs/manual_annotate.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### cat /scratch/aob2x/DESTv2_output_SNAPE/logs/runSnakemake.49369837*.err

### sbatch ~/CompEvoBio_modules/utils/snpCalling/scatter_gather_annotate/manual_annotate.sh
### sacct -j 4425471
### cat /scratch/aob2x/29Sept2025_ExpEvo/logs/manual_annotate.4425471*.err
# # ijob -A biol4559-aob2x -c10 -p largemem --mem=40G

module purge

module load htslib/1.17 bcftools/1.17 parallel/20200322 gcc/11.4.0 openmpi/4.1.4 R/4.3.1 samtools vcftools bedtools/2.30.0


popSet=all
method=PoolSNP
maf=001
mac=50
version=22Sept2026_ExpEvo
wd=/scratch/aob2x/compBio_SNP_22Sept2026
script_dir=~/CompEvoBio_modules/utils/snpCalling/
pipeline_output=/project/berglandlab/DEST/dest_mapped/

snpEffPath=~/snpEff

cd ${wd}

 echo "concat"
   ls -d ${wd}/sub_bcf/dest.*.${popSet}.${method}.${maf}.${mac}.${version}.norep.eff.vcf.gz  > \
   ${wd}/sub_bcf/vcf_order.genome

   bcftools concat \
   -f ${wd}/sub_bcf/vcf_order.genome \
   -O z \
   --threads 10 \
   -o ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.eff.vcf.gz

   tabix -p vcf ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz

Rscript --vanilla ~/DESTv3/snpCalling_dev/scatter_gather_annotate/vcf2gds.R \
${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz \
10
