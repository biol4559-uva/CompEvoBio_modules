sbatch --array=1-$( cat /project/biol4559-aob2x/mapping_scripts/aob2x/data/sras.txt | wc -l  ) \
/project/biol4559-aob2x/mapping_scripts/aob2x/scripts/mapping_script.sh \
/project/biol4559-aob2x/singularity \
/project/biol4559-aob2x/data/fastq/PROJECTID \
/project/biol4559-aob2x/mapping_output \
PATH_to_your_small_file
