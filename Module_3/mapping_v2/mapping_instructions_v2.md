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

3. Copy the contents of the ![NEW MAPPING SCRIPT](/Module_3/mapping_v2/mapping_script.sh) and paste it in the text file. Change your COMPUTEID to your computing ID. Determine if your samples are Single Ended or Paired End reads. Delete the lines of code (found in lines 43-70) that you DO NOT NEED.
<p align="center">
  <img src="/Module_3/mapping_v2/images/step6.jpg" width="1000"/>
</p>

4. Save the script in your scripts folder. In the image below, you can see it pointing to my `aob2x` folder. Make sure that you are saving to YOUR folder.
<p align="center">
  <img src="/Module_3/mapping_v2/images/step7.jpg" width="1000"/>
</p>

5. Save the file as `mapping_script.sh`
<p align="center">
  <img src="/Module_3/mapping_v2/images/step8.jpg" width="1000"/>
</p>

6. Create a new text file
<p align="center">
  <img src="/Module_3/mapping_v2/images/step4.jpg" width="1000"/>
</p>

7. Copy the contents of ["mapping_call.sh"](/Module_3/mapping_v2/mapping_call.sh) to the new file. Note, that this is a different version compared compared ot the first time. Replace the `COMPUTEID` with your computing ID. Replace the `PROJECTID` with the BioProject ID. If you have forgotten what the BioProject Id for your paper is, you can look it up in this [Excel file](/data/SRA_accessions_v2.xlsx):

<p align="center">
  <img src="/Module_3/mapping_v2/images/step9.jpg" width="1000"/>
</p>

8. Save that file as `mapping_call.sh` in the scripts folder.

9. Create a new text file
<p align="center">
  <img src="/Module_3/mapping_v2/images/step4.jpg" width="1000"/>
</p>

10. Find your "small csv file" wherever you have it. Copy the contents of that file to the new text file.

11. Save the file in the `data` folder that is inside your new directory. You'll see below that I am saving my file in my data folder. Notice the path and that it points to `aob2x`. Yours should point to your directory.
<p align="center">
  <img src="/Module_3/mapping_v2/images/step10.jpg" width="1000"/>
</p>

12. Once you have created and saved all three files (mapping_script.sh, sras.txt, mapping_call.sh), return to the mapping_call.sh script. Select all of the lines and copy them.
<p align="center">
  <img src="/Module_3/mapping_v2/images/step11.jpg" width="1000"/>
</p>

13. Open a new terminal window, and paste the lines in terminal. Record the SLURM job ID.
<p align="center">
  <img src="/Module_3/mapping_v2/images/step12.jpg" width="1000"/>
</p>

14. It might take a few moments (or maybe an hour) for the job to start. In the meantime, submit your job ID to Canvas.

15. have a good weekend!
