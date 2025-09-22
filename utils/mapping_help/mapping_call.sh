sbatch --array=1-$( cat /standard/BerglandTeach/misc/aob_small_file.csv | wc -l  ) \
~/CompEvoBio_modules/utils/mapping_help/mapping_script.aob2x.sh \
/standard/BerglandTeach/dest_freeze2.6.1_latest.sif \
/scratch/aob2x/compBio/fastq \
/standard/BerglandTeach/mapping_output \
/standard/BerglandTeach/misc/aob_small_file.csv

sacct -j 3920688

cat /scratch/aob2x/logs/RunDest.3920688_1.out

### try redownloading hte data
sbatch --array=1-31,42-48 \
~/CompEvoBio_modules/utils/mapping_help/mapping_script.aob2x.sh \
/standard/BerglandTeach/dest_freeze2.6.1_latest.sif \
/scratch/aob2x/compBio/fastq \
/standard/BerglandTeach/mapping_output \
/standard/BerglandTeach/misc/aob_small_file.csv
