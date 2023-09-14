# **Mapping pipeline output**

## Overview
The mapping pipeline that you ran (or are working to get running) spits out a lot of files. We will examine some of the important ones and discuss how our FASTQ data turned into this output.

Goals:
1. Learn about BAM files
2. Learn how to use samtools programs (`view`, `tview`)

Objectives:
1. Assemble the estimates of the PCR duplicate rates for your study & estimate the frequency of the Adh F/S polymorphism in your data

## Instructions
1. Open a new Rstudio session. We will just be using it as a text editor.

2. Open a new terminal screen and navigate to one of your mapping output folders. If you don't have one, you can get oriented by using any folder in the directory. That directory should look something like this:
<p align="center">
  <img src="/Module_4/images/mapping_directory.jpeg" width="1000"/>
</p>

3. First, we are going to look at the `*bam` files. These files are where the raw information about where the FASTQ reads map to. These files are compressed using a special algorithm, and we need to use a special tool to unpack them. To do this, we need to load the `samtools` module using the command:
```
module load samtools
```

4. To unpack the `bam` files, we will use the `samtools view` command. Be careful! If you don't include the `| less -S` your terminal screen will become flooded! You can use the arrow keys to pan left or right, up or down.

```
samtools view ExpEvo_SRP002024_CO_4_1975-MM-DD.original.bam | less -S
```

It should look something like this.
<p align="center">
  <img src="/Module_4/images/bam.jpeg" width="1000"/>
</p>

Some of these field should be familiar. What is the first column? The columns with "2L" and numbers represent where the read mapping starts. The "54M" (or whatever you have) represents how the read maps to the genome. You can learn to decode the info in your file using information [here](https://jef.works/blog/2017/03/28/CIGAR-strings-for-dummies/). "54M" means that the read Matches at all 54 nucleotides in it.

5. To make better sense of this file, we can use a more visual approach. We will use the program `tview` in samtools to accomplish this. Modify the command below for your file:
```
samtools tview ExpEvo_SRP002024_CO_4_1975-MM-DD.original.bam /project/biol4559-aob2x/data/refGenome/holo_dmel_6.12.fa
```

If your sample has low coverage, like mine, the screen will look quite blank. This is because we are looking at the very tip of the chromosome and mapping is difficult there. If your sample has a ton of reads, you might see something.  We are going to navigate to our famous Adh F/S polymorphism that we discussed on the first day of class. To do this, type `g` in the `tview` window and then type: `2L:14617051`. Your cursor shoudl be pointed right over the site. Is there a polymorphism at your site? What is its frequency? How many reads cover the site (to exit tview press `q`)
<p align="center">
  <img src="/Module_4/images/adh.jpeg" width="500"/>
</p>

6. One important aspect of your data is the read depth. This is especially true for Pool-Seq data (as we will discuss in Class 5). One factor that can influence read depth is PCR duplication. During the library prep stage, PCR is performed. PCR is an important step for preparring the actual DNA fragment for sequencing by attaching adapter sequence to the 5' ends. It also copies fragments of DNA, and amplifies the library to a sufficient concentration. If you overdo PCR, then some specific fragments can amplify many many times and can represent a large fraction of your sequencing library. During the mapping pipeline, we identified these duplicated fragments and flagged them in the BAM file. Our downstream analysis discards these duplicated reads. This [post on StackExchange](https://bioinformatics.stackexchange.com/questions/2866/how-do-pcr-duplicates-arise-and-why-is-it-important-to-remove-them-for-ngs-analy) describes the issue well:
<p align="center">
  <img src="/Module_4/images/seqAnswer_pcrDup.jpeg" width="1000"/>
</p>

7. Fortunately, all we have to do to assess PCR duplicate rate is look at the report. Create an Excel (or Google Sheet) file with the PCR duplicate rates for all of the samples in your study. You can do this manually by navigating to each file, and typing this command:
```
tail -n4 ExpEvo_SRP002024_CO_4_1975-MM-DD.mark_duplicates_report.txt | head -n2
```

You shoudl get something like this; the 9th column should be the duplication rate (3.6% in my case):
<p align="center">
  <img src="/Module_4/images/dupRate.jpeg" width="1000"/>
</p>

8. To complete this assignment, upload your estimates to Canvas and determine if the authors of your study reported performing PCR de-duplication.
