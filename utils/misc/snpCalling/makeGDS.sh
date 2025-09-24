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

### sbatch /scratch/aob2x/CompEvoBio_modules/utils/snpCalling/makeGDS.sh
### sacct -j 54019500
### cat /scratch/aob2x/DESTv2_output_26April2023/logs/manual_annotate.49572492*.out



Rscript --vanilla /scratch/aob2x/CompEvoBio_modules/utils/snpCalling/scatter_gather_annotate/vcf2gds.R \
/scratch/aob2x/compBio_SNP_25Sept2023/dest.expevo.PoolSNP.001.50.11Oct2023.norep.ann.vcf.gz
