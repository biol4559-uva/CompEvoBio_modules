# **Short Read Mapping & Pipelines**

## Overview
Well, we are going to have to try this again. In this second version of the mapping script, we are going to be working out of the `/project` directory. That means that all of your scripts, input files, and error logs will be put on `/project` in  folders that I have made for you. The benefit of this approach is that Connor and I can see all of your files easily in one place and can help debug things if they arise. I'm fairly certain that most of the issues that people are having are with path names to files plus some oddball formatting issues.

Please follow the instructions as carefully as you can! If you run into issues, drop me a line and I will try to debug.

## Instructions
1. Start a new RStudio job. Close any open scripts. We are going to change Rstudio's working directory to the new directory that I made for you. To do this, click on the button on the right hand side with the 3 periods. The arrow in the immage below shows you where that button is. Type in the following path into the box that opens. <b>CHANGE `aob2x` to your usename</b>

<p align="center">
  <img src="/Module_3/mapping_v2/images/step1.jpg" width="1000"/>
</p>

2. Create a new text file by going to the File menu, New File, Text File. The arrow points to it below.
<p align="center">
  <img src="/Module_3/mapping_v2/images/step4.jpg" width="1000"/>
</p>

3. Copy the contents of the ![NEW MAPPING SCRIPT](/Module_3/mapping_v2/mapping_instructions_v2.md) and paste it in the text file. Change your COMPUTEID to your computing ID
<p align="center">
  <img src="/Module_3/mapping_v2/images/step6.jpg" width="1000"/>
</p>




Using the [csv table](/Module_3/sras.txt)that you and your group constructed, divide up the samples roughly evenly. Copy your lines to a new file and save in your assignment folder.

2. Copy the `mapping_script.sh` to a new file using Rstudio. Modify that file to have the proper paths and decide if you need PE or SE. Copy the lines below to a dummy text file, and modify `PATH_to_your_small_file` to point it to the small csv file you made in the previous step. Also modify the `PROJECTID` to include the proper SRA project ID. Record the job number.

```
sbatch --array=1-$( PATH_to_your_small_file  | wc -l  ) \
PATH_to_your_version_of_mapping_script.sh \ # The script
/project/biol4559-aob2x/singularity \ # Argument 1: Where the SIF file is located
/project/biol4559-aob2x/data/fastq/PROJECTID \ # Argument 2: Where the reads are located
/project/biol4559-aob2x/mapping_output \ # Argument 3: Output folder
PATH_to_your_small_file
```

3. Confirm your job is successfully running by running the command `sacct -j JOBID`

4. To complete this assignment, upload your JOBID to Cavnas.
