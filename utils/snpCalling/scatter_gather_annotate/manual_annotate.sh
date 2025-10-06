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
### sacct -j 4308447
### cat /scratch/aob2x/29Sept2025_ExpEvo/logs/manual_annotate.4308447*.err
# # ijob -A biol4559-aob2x -c10 -p largemem --mem=40G

module purge

module load htslib/1.17 bcftools/1.17 parallel/20200322 gcc/11.4.0 openmpi/4.1.4 R/4.3.1 samtools vcftools bedtools/2.30.0



popSet=all
method=PoolSNP
maf=001
mac=50
version=29Sept2025_ExpEvo
wd=/scratch/aob2x/compBio_SNP_29Sept2025
script_dir=~/CompEvoBio_modules/utils/snpCalling/
pipeline_output=/project/berglandlab/DEST/dest_mapped/

snpEffPath=~/snpEff

cd ${wd}

echo "no rep & index"

  noRepIndex () {

    popSet=all
    method=PoolSNP
    maf=001
    mac=50
    version=29Sept2025_ExpEvo
    wd=/scratch/aob2x/compBio_SNP_29Sept2025
    script_dir=~/CompEvoBio_modules/utils/snpCalling/
    pipeline_output=/project/berglandlab/DEST/dest_mapped/
    chr=${1}


    bedtools intersect -sorted -v -header \
    -b ${script_dir}/scatter_gather_annotate/repeat_bed/repeats.sort.bed.gz \
    -a $bcf_outdir/dest.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz |
    bgzip -c > \
    $bcf_outdir/dest.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz

    bcftools index -f $bcf_outdir/dest.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz
  }
  export -f noRepIndex

  parallel -j5 noRepIndex ::: 2L 2R 3L 3R X

 echo "concat"
   ls -d ${wd}/sub_bcf/dest.*.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz | grep -E "2L|2R|3L|3R|X" > \
   ${wd}/sub_bcf/vcf_order.genome

   bcftools concat \
   -f ${wd}/sub_bcf/vcf_order.genome \
   -O z \
   --threads 10 \
   -o ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz

   tabix -p vcf ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz


 echo "convert to vcf & annotate"
   bcftools view \
   --threads 48 \
   ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz | \
   java -jar ~/snpEff/snpEff.jar \
   eff \
   BDGP6.86 - > \
   ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf

echo "make GDS"
   Rscript --vanilla ~/CompEvoBio_modules/utils/snpCalling/scatter_gather_annotate/gds2vcf.R ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf

echo "bgzip & tabix"
  bgzip -@10 -c ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf > ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf.gz
  tabix -p vcf ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf.gz
