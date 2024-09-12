# **FASTQ files**

## Overview
In this exercise we will learn about the basic nature of short-read data. Short-read data is stored in `fastq` files. These files are simple text files and we can look directly at them. A single read is represented by four lines. Lines 1 and 3 are headers that give information about the sample. In this example, we learn that the sample name is `@SRR036932.24`. The second column of data reflects the machine that produced the data, as well as the physical location on the slide where the read was read. The last column tells you that this read has 54 bases.

<p align="center">
  <img src="/Module_2/images/fastq.jpeg" width="1000"/>
</p>

The Wikipedia page for FASTQ is very informative:
https://en.wikipedia.org/wiki/FASTQ_format

Goals:
1. Practice using bash commands to look at FASTQ files
2. Learn about quality encodings
3. Learn about single-end vs paired-end reads

Objective:
1. Determine the Illumina version that produced two sample data files

## Instructions
0. Open a terminal window using OpenOnDemand

1. Navigate to `/standard/BerglandTeach/data/qual_encoding`.

2. Look at the contents of the two sample files. Wow! That is a lot of lines. What are strategies for not overloading your terminal window with TMI? Here are some commands to try:

```
cat
less -S
head
```

3. What version of the Illumina platform generated these samples? Here is a table from the wikipedia page that you might find helpful.

<p align="center">
  <img src="/Module_2/images/fastq_encoding.jpeg" width="1000"/>
</p>

3. Most modern Illumina data comes as "paired-end" data as described in the video that we watched. Paried end reads represent DNA sequence at either end of the DNA insert. Paired end data is represented by two files, often with file names that include `_1` or `_2` to designate the first and second of the paired read. When the reads are mapped back to the reference genome, both pairs of a read are treated as a single unit and the mapping programs need to be able to know which reads belong together. How do we tell the program that the reads are linked?

4. One of your friends generated some short-read sequence data and wants your help analyzing it. But, they messed up by naming their files poorly. They need your help. Can you figure out which samples are the proper `Read1` and `Read2` pairs? What clues did you use? You can find your friend's data here: `/project/biol4559-aob2x/data/mixed_up_pairs/`

4. To complete this assignment submit your answers on Canvas.
