
### try redownloading hte data
# sbatch --array=1-31,42-48 \
# ~/CompEvoBio_modules/utils/mapping_help/mapping_script.aob2x.sh \
# /standard/BerglandTeach/dest_freeze2.6.1_latest.sif \
# /scratch/aob2x/compBio/fastq \
# /standard/BerglandTeach/mapping_output \
# /standard/BerglandTeach/misc/aob_small_file.csv
#
# sacct -j 3921031
# cat /scratch/aob2x/logs/RunDest.3921031_17.err


### PRJNA315172
#   sbatch --array=1-6 \
#   ~/CompEvoBio_modules/utils/mapping_help/mapping_script.aob2x.sh \
#   /standard/BerglandTeach/dest_freeze2.6.1_latest.sif \
#   /scratch/aob2x/compBio/fastq \
#   /standard/BerglandTeach/mapping_output \
#   /standard/BerglandTeach/misc/aob_small_file.csv
#
#   sacct -j 3947077
#   cat /scratch/aob2x/logs/RunDest.3947077_2.err

## PRJNA170455

## PRJNA657615
sbatch --array=18-31 \
~/CompEvoBio_modules/utils/mapping_help/mapping_script.aob2x.sh \
/standard/BerglandTeach/dest_freeze2.6.1_latest.sif \
/scratch/aob2x/compBio/fastq \
/standard/BerglandTeach/mapping_output \
/standard/BerglandTeach/misc/aob_small_file.csv

sacct -j 3954329
cat /scratch/aob2x/logs/RunDest.3947144_19.err

cd ExpEvo_PRJNA657615_H_pop3_gen12_1999-10-01


PRJNA194129
PRJNA315172
