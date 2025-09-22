sbatch --array=1-$( cat /standard/BerglandTeach/misc/aob_small_file.csv | wc -l  ) \
~/CompEvoBio_modules/utils/mapping_help/mapping_script.aob2x.sh \
/standard/BerglandTeach/dest_freeze2.6.1_latest.sif \
/scratch/aob2x/compBio/fastq/ \
/standard/BerglandTeach/mapping_output \
/standard/BerglandTeach/misc/aob_small_file.csv
