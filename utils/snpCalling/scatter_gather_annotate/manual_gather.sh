#!/bin/bash
#
#SBATCH -J manual_gather # A single job name for the array
#SBATCH --ntasks-per-node=48 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 14:00:00 ### 1 hours
#SBATCH --mem 20G
#SBATCH -o /scratch/aob2x/29Sept2025_ExpEvo/manual_gather.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/29Sept2025_ExpEvo/logs/manual_gather.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch ~/CompEvoBio_modules/utils/snpCalling/scatter_gather_annotate/manual_gather.sh
### sacct -j 4393980
### cat /scratch/aob2x/29Sept2025_ExpEvo/logs/manual_gather.4393980_*.err
### cat /scratch/aob2x/compBio_SNP_25Sept2023/logs/manual_gather
### cd /scratch/aob2x/compBio_SNP_25Sept2023

module load htslib/1.17  bcftools/1.17 parallel/20200322 gcc/11.4.0 openmpi/4.1.4 python/3.11.4 perl/5.40.2 vcftools/0.1.16

concatVCF() {


  popSet=all
  method=PoolSNP
  maf=001
  mac=50
  version=29Sept2025_ExpEvo
  wd=/scratch/aob2x/compBio_SNP_29Sept2025
  script_dir=~/CompEvoBio_modules/utils/snpCalling/
  pipeline_output=/project/berglandlab/DEST/dest_mapped/


  # chr=2L

  chr=${1}

  echo "Chromosome: $chr"

  bcf_outdir="${wd}/sub_bcf"
  if [ ! -d $bcf_outdir ]; then
      mkdir $bcf_outdir
  fi

  outdir=$wd/sub_vcfs
  cd ${wd}

  echo "generate list"
  #ls -d *.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz | grep '^${chr}_' | sort -t"_" -k2n,2 -k4g,4 \
  #> $outdir/vcfs_order.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.sort


  ls -d ${outdir}/*.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz | \
  rev | cut -f1 -d '/' |rev | grep -E "^${chr}_" | sort -t"_" -k2n,2 -k4g,4 | \
  sed "s|^|$outdir/|g" > $outdir/vcfs_order.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.sort

  # less -S $outdir/vcfs_order.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.sort


  echo "Concatenating"

  bcftools concat \
  -f $outdir/vcfs_order.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.sort \
  -O z \
  -n \
  --threads 48 \
  -o $bcf_outdir/dest.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz

  # vcf-concat \
  # -f $outdir/vcfs_order.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.sort \
  # -s | \
  # bgzip -c > $bcf_outdir/dest.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz

  tabix -p vcf $bcf_outdir/dest.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz

  bedtools intersect -sorted -v -header \
  -b ${script_dir}/scatter_gather_annotate/repeat_bed/repeats.sort.bed.gz \
  -a $bcf_outdir/dest.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz |
  bgzip -c > \
  $bcf_outdir/dest.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz

}
export -f concatVCF

parallel -j1 concatVCF ::: 2L 2R 3L 3R 4 mitochondrion X Y
#parallel -j8 concatVCF ::: 3L
