# **SLURM**

## Overview
Up to this point, we have been using Rivanna in an interactive way. That is great when we are analyzing data and learning to code. But, sometimes, we need to let Rivanna work for hours and we want it to run in the background. To submit jobs to Rivanna, we use the job handler SLURM.

Goals:
1. learn about the structure of a SLURM file
2. learn about modules on Rivanna
3. learn about batch processing

Objective:
1. Run a quality control tool on FASTQ files for your paper using a slurm script.

## Instructions -
### Basics of a SLURM script
1. Open a new Rstudio job. For this exercise, we are only going to RStudio to work as a script editor and to use GitHub. Once you have Rstudio open, load your GitHub Repo project, and create a new script. Copy the contents of the contents of `SLURM_demo_1.sh` into the script. Save the file in your GitHub folder in a new folder called "SLURM_assignment". Make sure to save the file name as "SLURM_demo_1.sh"

2. Open a terminal window and make a new "log" directory in your scratch directory. To do this, nagivate to your scratch directory (`/scratch/COMPUTE_ID`). Then, make a new directory (`mkdir logs`)

3. In RStudio, edit the script text to point to your scratch directory (lines 8 & 9).

4. Modify line 13 give the full and path name to the SLURM_demo_1.sh script. Save everything.

5. Copy the text on line 13 after the "run as: " bit. It should be something like `sbatch ~/biol4559/aob2x/SLURM_assignment/SLURM_demo_1.sh`. Go back to your terminal, and paste that in the shell.

6. Copy the job job number and paste it back on line 14, where the XXXXXXX is. Now, copy the `sacct -j XXXXXXXXX` line and paste it back in the terminal. What is the status of your job?

7. Look at the contents of the log files. What do they say?


### Basics of batch processing
1. We want to be able to run fastqc on all of the samples associated with our paper. To do this, we will utilize batch processing which is a bit like a foreach loop in R. We have a list of samples we want to process, and each sample is processed independently. For SLURM, this is called an array.

2. To run a slurm script as an array, we will type:

```
sbatch --array=1-5 SLURM_demo_2.sh
```

In this way, our slurm script is a bit like a function that only takes one parameter, a number between X and Y in this case, between 1-5. Within the slurm script, that parameter gets caputred by the `${SLURM_ARRAY_TASK_ID}` variable. To run this script, you will need to copy it into R, and modify the paths, and save on Rivanna. When you run it, what files do you look at to ensure it works?

3. Let's say that we have a list of sample names, and we want to run `fastqc` on each of them. In order to use the `${SLURM_ARRAY_TASK_ID}` variable, we need to link it to a premade table of possible values. Then, we can look up the ith line and use that information as our parameter values that we actually need. Copy `SLURM_demo_3.sh` and modify as appropriate.

### Run FASTQC on all of the samples from your paper
1. Copy and modify `SLURM_demo_4.sh` so that it fits with your path names, as well as your paper `proj` variable. Run it, and check that it is running.

### Run MULTIQC.
1. Combine your FASTQC output files into a single comparative file using MULTIQC in an interactive terminal:
```
module load multiqc
multiqc -o ~/multiQC_output/ /scratch/COMPUTEID/fastq_QC_Out/
```

2. Download the two contents of that foler and open the html file.

### To complete this assignment
6. To complete this assignment take screen shots of the General Stats table, and the Sequence Quality Histograms figures and upload them to Canvas. Write 3-4 sentences to describe these figures. Watch the [instructional video](https://www.youtube.com/watch?v=fcSIM1hkIBo). Are there any samples that you are concerned about? Why?
