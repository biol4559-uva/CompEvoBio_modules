# ijob -A berglandlab -c1 -p standard --mem=6G


module load apptainer


apptainer -v pull /standard/BerglandTeach/multiqc.sif docker://multiqc/multiqc:latest
apptainer exec /standard/BerglandTeach/multiqc.sif multiqc .