#!/bin/bash
#
#SBATCH -J manual_annotate # A single job name for the array
#SBATCH --ntasks-per-node=5 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 12:00:00 ### 1 hours
#SBATCH --mem 40G
#SBATCH -o /scratch/aob2x/29Sept2025_ExpEvo/logs/manual_annotate.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/29Sept2025_ExpEvo/logs/manual_annotate.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### cat /scratch/aob2x/DESTv2_output_SNAPE/logs/runSnakemake.49369837*.err

### sbatch ~/CompEvoBio_modules/utils/snpCalling/scatter_gather_annotate/manual_annotate_slices.sh
### sacct -j 4425471
### cat /scratch/aob2x/29Sept2025_ExpEvo/logs/manual_annotate.4425471*.err
# # ijob -A biol4559-aob2x -c10 -p largemem --mem=40G

module purge

module load gcc/14.2.0 xz/5.8.3 htslib/1.23 bcftools/1.23

popSet=all
method=PoolSNP
maf=001
mac=50
version=22Sept2026_ExpEvo
wd=/scratch/aob2x/compBio_SNP_22Sept2026
script_dir=~/CompEvoBio_modules/utils/snpCalling/
pipeline_output=/project/berglandlab/DEST/dest_mapped/
job=${SLURM_ARRAY_TASK_ID}    # job=1

snpEffPath=~/snpEff

cd ${wd}/sub_vcfs

input_file=$( ls -d *norep.vcf.gz | sed -n "${job}p" )
output_file=$( echo ${input_file} | sed 's/norep/norep.eff/g' | sed 's/\.gz//g' )
echo ${input_file}
echo ${output_file}

echo "convert to vcf & annotate"
   bcftools view \
   --threads 5 \
   ${input_file} | \
   java -jar ~/snpEff/snpEff.jar \
   eff \
   BDGP6.86 - | head -n 100 > \
   ${output_file}

echo "bgzip & tabix"
  bgzip -@5 -c ${output_file} > ${output_file}.gz
  tabix -p vcf ${output_file}.gz
