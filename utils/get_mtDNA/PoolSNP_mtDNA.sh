#!/usr/bin/env bash

#SBATCH -J manual_scatter # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 12:00:00 ### 1 hours
#SBATCH --mem 64G
#SBATCH -o /scratch/aob2x/29Sept2025_ExpEvo/logs/manual_gather.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/29Sept2025_ExpEvo/logs/manual_gather.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard


# ijob -A berglandlab -c10 -p standard --mem=50G
# sbatch ~/CompEvoBio_modules/utils/get_mtDNA/PoolSNP_mtDNA.sh
# sacct -j 5024727 | grep -v "COMPLE"
# cat /scratch/aob2x/29Sept2025_ExpEvo/logs/manual_gather.5024727_*.err


  module purge
  trap 'rm -rf ${tmpdir}' EXIT

  #module load htslib bcftools parallel intel/18.0 intelmpi/18.0 mvapich2/2.3.1 R/3.6.3 python/3.6.6 vcftools/0.1.16
  #module load htslib/1.10.2 bcftools/1.9 parallel/20200322 intel/18.0 intelmpi/18.0 R/3.6.3 python/3.6.6 vcftools/0.1.16
  module load htslib/1.17  bcftools/1.17 parallel/20200322 gcc/11.4.0 openmpi/4.1.4 python/3.11.4 vcftools/0.1.16 R/4.3.1
  module load bedtools/2.30.0

  cat /scratch/aob2x/allpops.mtDNA.sites | python ~/CompEvoBio_modules/utils/snpCalling/PoolSNP/PoolSnp.py \
  --sync - \
  --min-cov 4 \
  --max-cov 0.999 \
  --min-count 4 \
  --min-freq 0.001 \
  --miss-frac 0.5 \
  --names $( cat /standard/BerglandTeach/mtDNA/sync/allpops.mtDNA.names |  tr '\n' ',' | sed 's/,$//g' ) > /standard/BerglandTeach/mtDNA/mtDNA.vcf

  cat /standard/BerglandTeach/mtDNA/mtDNA.vcf | sed 's/mitochondrion_genome/dmel_mitochondrion_genome/g' | \
  java -jar ~/snpEff/snpEff.jar \
  eff \
  BDGP6.86 - > \
  /standard/BerglandTeach/mtDNA/mtDNA.ann.vcf

  echo "make GDS"
  Rscript --vanilla ~/CompEvoBio_modules/utils/snpCalling/scatter_gather_annotate/gds2vcf.R /standard/BerglandTeach/mtDNA/mtDNA.ann.vcf
