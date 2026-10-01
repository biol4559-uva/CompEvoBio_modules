#!/bin/bash
#
#SBATCH -J manual_gather # A single job name for the array
#SBATCH --ntasks-per-node=48 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 1:00:00 ### 1 hours
#SBATCH --mem 20G
#SBATCH -o /scratch/aob2x/compBio_SNP_22Sept2026/logs/manual_gather.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/compBio_SNP_22Sept2026/logs/manual_gather.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch ~/CompEvoBio_modules/utils/snpCalling/scatter_gather_annotate/manual_gather.sh
### sacct -j 20672612
### cat /scratch/aob2x/29Sept2025_ExpEvo/logs/manual_gather.4393980_*.err
### cat /scratch/aob2x/compBio_SNP_25Sept2023/logs/manual_gather
### cd /scratch/aob2x/compBio_SNP_25Sept2023

# # ijob -A biol4020-aob2x -c10 -p standard --mem=40G


module purge
module load gcc/14.2.0 xz/5.8.3 htslib/1.23 bcftools/1.23 parallel/20250722
concatVCF() {


  popSet=all
  method=PoolSNP
  maf=001
  mac=50
  version=22Sept2026_ExpEvo
  wd=/scratch/aob2x/compBio_SNP_22Sept2026
  script_dir=~/CompEvoBio_modules/utils/snpCalling/
  pipeline_output=/project/berglandlab/DEST/dest_mapped/



  # chr=3L

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


  ls -d ${outdir}/*.${popSet}.${method}.${maf}.${mac}.${version}.norep.eff.vcf.gz | \
  rev | cut -f1 -d '/' |rev | grep -E "^${chr}_" | sort -t"_" -k2n,2 -k4g,4 | \
  sed "s|^|$outdir/|g" > $outdir/vcfs_order.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.sort

  # less -S $outdir/vcfs_order.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.sort
  #sed -i '$d' $outdir/vcfs_order.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.sort | tail

  echo "Concatenating"

  bcftools concat \
  -f $outdir/vcfs_order.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.sort \
  -O z \
  --threads 10 \
  -o $bcf_outdir/dest.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.norep.eff.vcf.gz

 
  tabix -p vcf $bcf_outdir/dest.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.norep.eff.vcf.gz

}
export -f concatVCF

parallel -j1 concatVCF ::: 2L 2R 3L 3R 4 mitochondrion_genome X Y
#parallel -j8 concatVCF ::: 3L
