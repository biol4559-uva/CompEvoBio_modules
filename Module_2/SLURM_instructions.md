# **SLURM**

## Overview
Up to this point, we have been using Rivanna in an interactive way. That is great when we are analyzing data and learning to code. But, sometimes, we need to let Rivanna work for hours and we want it to run in the background. To submit jobs to Rivanna, we use the job handler SLURM.

Goals:
1. learn about the structure of a SLURM file
2. learn about modules on Rivanna

Objective:
1. Run a quality control tool on FASTQ files for your paper using a slurm script

## Instructions
1. Open a new Rstudio job. For this exercise, we are only going to RStudio to work as a script editor and to use GitHub. Once you have Rstudio open, load your GitHub Repo project, and create a new script. Copy the contents of the contents of `SLURM_template.sh` into the script.

2. Lines 1-11 are mandatory for a SLURM script. What do each of the lines mean?

3. Lines 13-15 are how we will ultimately run the script

4. Lines 17-30 load modules, and run a small example. Run these lines to run fastqc. Download the output of that program to your computer using the Files menu on OpenOnDemand. What is the basic pattern of change in quality across the sequence read?

<p align="center">
  <img src="/Module_2/images/download_fastqc.jpeg" width="1000"/>
</p>

5. Once your paper is assigned, you will be running FASTQC on all of the samples associated with your paper. Download the multiqc report and describe one of the figures. To modify the slurm example script to run well, you need to do the following: <br>
• In terminal run this line `mkdir /scratch/aob2x/logs/` <br>
• Replace "COMPUTE_ID" with your actual computeID<br>
• Get the correct PATH to file on line 13<br>
• Comment out line 30<br>
• Using the terminal, copy everything on line 13 starting at `sbatch`<br>
• If your job is submitted, you will be returned a job ID. Copy that and paste it in the XXXX'd out section on line 14. Run that line (without the # symbol). Is it running successfully? If not, look at your log files. What do they say?<br><br>

6. To complete this assignment take screen shots of the General Stats table, and the Sequence Quality Histograms figures and upload them to Canvas. Write 3-4 sentences to describe these figures.
