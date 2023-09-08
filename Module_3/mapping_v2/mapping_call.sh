sbatch --array=1-$( cat /project/biol4559-aob2x/mapping_scripts/COMPUTEID/data/sras.txt | wc -l  ) \
/project/biol4559-aob2x/mapping_scripts/COMPUTEID/scripts/mapping_script.sh \
/project/biol4559-aob2x/singularity \
/project/biol4559-aob2x/data/fastq/PROJECTID \
/project/biol4559-aob2x/mapping_output \
/project/biol4559-aob2x/mapping_scripts/COMPUTEID/data/sras.txt
