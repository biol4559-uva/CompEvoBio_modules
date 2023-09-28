# **Data shuffling**

## Overview
Have you been getting the emails about "project" storage? I have. The situation will get resolved and no data seems to be lost. But, for the next few weeks (months?) `/project` storage will be a bit slow to access. They are upgrading the hardware that this storage option uses, and I think in the process things are slowing down. So, the folks at Rivanna suggest that we work off of `/scratch`.

The reason why this is important moving forward is that we have a single ~11Gb compressed file that you will be using for most of rest of the class. This file contains the allele frequency estimates + a bunch of other info for your samples, plus about 750 Pool-Seq samples from across the world.

This big file is something that we can easily access using R. You might be interested to learn that when you typically use R, you are loading an entire file into RAM. This is a limitation of R, especially when you are using your laptop. The 11Gb file uncompressed is something like 50Gb. So, we do not want to load the whole file in. Also, loading the whole file in would take hours.

The nice thing about GDS files is that you can access specific SNPs and samples and other bits and boodles of information quickly from the harddrive, and then load them into R. That means that you need to have very fast read-speeds from the drive. Clearly this will be an issue if we are trying to work off of `/project`.

So, we are going to transfer the files you need, including a new folder of packges to your scratch. We are going to use the Globus interface. This interface allows for very fast transfer between drives on Rivanna.

Rivanna has more detailed information here

Goals:
1. Learn about Globus

Objective:
1. Transfer data to `/scratch/COMPUTEID`

## Instructions
1. Log into (Globus)[https://globus.org].
<p align="center">
  <img src="/Module_6/images/globus1.jpeg" width="1000"/>
</p>

2. Click in the text field next to "Collection" and type UVA Standard
<p align="center">
  <img src="/Module_6/images/globus2.jpeg" width="1000"/>
</p>

3. You'll see three or four listings. I'm a little curious, actually, what you can see. We'll either be pulling data from `/standard` or `/project` to our `/scratch`. Click on the diagonal arrow button.
<p align="center">
 <img src="/Module_6/images/globus3.jpeg" width="1000"/>
</p>

4. We are going to navigate our way to the folder that contains `module_6` and copy that to our `/scratch`. Click the check box beside module_6, and click start.
<p align="center">
 <img src="/Module_6/images/globus4.jpeg" width="1000"/>
</p>
