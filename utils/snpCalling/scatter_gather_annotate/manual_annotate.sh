#!/bin/bash
#
#SBATCH -J manual_annotate # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6:00:00 ### 1 hours
#SBATCH --mem 40G
#SBATCH -o /scratch/aob2x/compBio_SNP_25Sept2023/logs/manual_annotate.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/compBio_SNP_25Sept2023/logs/manual_annotate.%A_%a.err # Standard error
#SBATCH -p largemem
#SBATCH --account biol4559-aob2x

### cat /scratch/aob2x/DESTv2_output_SNAPE/logs/runSnakemake.49369837*.err

### sbatch /scratch/aob2x/CompEvoBio_modules/utils/snpCalling/scatter_gather_annotate/manual_annotate.sh
### sacct -j 53998649
### cat /scratch/aob2x/DESTv2_output_26April2023/logs/manual_annotate.49572492*.out
# # ijob -A biol4559-aob2x -c10 -p largemem --mem=40G

module purge

module load  htslib/1.10.2 bcftools/1.9 intel/18.0 intelmpi/18.0 parallel/20200322 R/3.6.3 samtools vcftools


popSet=PoolSeq
method=PoolSNP
maf=001
mac=50
version=11Oct2023
wd=/scratch/aob2x/compBio_SNP_25Sept2023

snpEffPath=~/snpEff

cd ${wd}


# echo "index"
#   bcftools index -f ${wd}/sub_bcf/dest.2L.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz
#   bcftools index -f ${wd}/sub_bcf/dest.2R.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz
#   bcftools index -f ${wd}/sub_bcf/dest.3L.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz
#   bcftools index -f ${wd}/sub_bcf/dest.3R.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz
#   bcftools index -f  ${wd}/sub_bcf/dest.X.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz
#   bcftools index -f  ${wd}/sub_bcf/dest.Y.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz
#   bcftools index -f  ${wd}/sub_bcf/dest.4.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz
#   bcftools index -f ${wd}/sub_bcf/dest.mitochondrion_genome.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz
#
#
 echo "concat"
   # bcftools concat \
   # -O z \
   # ${wd}/sub_bcf/dest.*.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz \
   # -o ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.new.vcf.gz

   ls -d ${wd}/sub_bcf/dest.*.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz | grep -E "2L|2R|3L|3R|X" > \
   ${wd}/sub_bcf/vcf_order.genome

   bcftools concat \
   -f ${wd}/sub_bcf/vcf_order.genome \
   -O z \
   -o ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz

   # vcf-concat \
   # -f ${wd}/sub_bcf/vcf_order.genome \
   # -s \
   # |  \
   # bgzip -c > ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz

   tabix -p vcf ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz

   # ijob -A berglandlab -c20 -p standard --mem=60G

 echo "convert to vcf & annotate"
   bcftools view \
   --threads 20 \
   ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz | \
   java -jar ~/snpEff/snpEff.jar \
   eff \
   BDGP6.86 - > \
   ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf


# echo "fix header" #this is now fixed in PoolSNP.py
# #  sed -i '0,/CHROM/{s/AF,Number=1/AF,Number=A/}' ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf
# #  sed -i '0,/CHROM/{s/AC,Number=1/AC,Number=A/}' ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf
# #  sed -i '0,/CHROM/{s/AD,Number=1/AD,Number=A/}' ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf
# #  sed -i '0,/CHROM/{s/FREQ,Number=1/FREQ,Number=A/}' ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf
#
# #  bcftools view -h ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf > ${wd}/tmp.header
# #
# #  bcftools reheader --threads 10 -h ${wd}/tmp.header -o ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.header.bcf ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.bcf
#
# echo "make GDS"
   Rscript --vanilla /scratch/aob2x/DESTv2/snpCalling/scatter_gather_annotate/vcf2gds.R ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf

echo "bgzip & tabix"
  bgzip -@20 -c ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf > ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf.gz
  tabix -p vcf ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf.gz

#  bgzip ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf.gz

  #bcftools sort \
  #-o /scratch/aob2x/DESTv2_output_26April2023/dest.all.PoolSNP.001.50.26April2023.norep.ann.sort.vcf.gz \
  #-O z \
  #/scratch/aob2x/DESTv2_output_26April2023/dest.all.PoolSNP.001.50.26April2023.norep.ann.vcf.gz

  #bgzip -d -c -@20 /scratch/aob2x/DESTv2_output_26April2023/dest.all.PoolSNP.001.50.26April2023.norep.ann.vcf.gz > /scratch/aob2x/DESTv2_output_26April2023/dest.all.PoolSNP.001.50.26April2023.norep.ann.vcf

  #cat /scratch/aob2x/DESTv2_output_26April2023/dest.all.PoolSNP.001.50.26April2023.norep.ann.vcf | vcf-sort > /scratch/aob2x/DESTv2_output_26April2023/dest.all.PoolSNP.001.50.26April2023.norep.ann.sort.vcf

  #bgzip -c -@20 /scratch/aob2x/DESTv2_output_26April2023/dest.all.PoolSNP.001.50.26April2023.norep.ann.sort.vcf > /scratch/aob2x/DESTv2_output_26April2023/dest.all.PoolSNP.001.50.26April2023.norep.ann.sort.vcf.gz

  # bcftools reheader \
  # -o /scratch/aob2x/DESTv2_output_26April2023/dest.all.PoolSNP.001.50.26April2023.norep.ann.reheader.vcf.gz \
  # --threads 20 \
  # --fai /scratch/aob2x/DESTv2/snpCalling/scatter_gather_annotate/holo_dmel_6.12.fa.fai \
  # /scratch/aob2x/DESTv2_output_26April2023/dest.all.PoolSNP.001.50.26April2023.norep.ann.vcf.gz
  #

#  Rscript --vanilla /scratch/aob2x/DESTv2/snpCalling/scatter_gather_annotate/gds2vcf.R ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.gds
#  tabix -p vcf ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.new.vcf.gz


# rm ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf
# rm ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf
