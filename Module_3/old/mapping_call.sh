sbatch --array=1-$( cat /scratch/aob2x/small_file.csv  | wc -l  ) \
~/CompEvoBio_modules/Module_3/old/mapping_script.aob2x.sh \
/scratch/aob2x/dest_v2.6.1_latest.sif \
/scratch/aob2x/compBio/fastq/SRP002024 \
/standard/BerglandTeach/mapping_output \
/scratch/aob2x/small_file.csv
