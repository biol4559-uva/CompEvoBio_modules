# **Navigating around Rivanna exercise instructions**

## Overview
In this exercise you will explore Rivanna's file structure and gain basic skills in how to use `bash`. `Bash` is a a program that lets you interact and script with unix computers, and using `bash` you invoke other programs, or execute scripts wiht a string of commands, or submit jobs to Rivanna's job handler. Mastering `bash` can take a while and in this course we will be covering the basics.

Goals:
1. Practice using bash.

Objective:
1. Explore the file structure on Rivanna using `bash`
2. Learn about manual pages, flags.

To learn more click here:
https://ucsbcarpentry.github.io/2021-01-21-SWC-Bash-online/02-filedir/

## Instructions
1. Open up the shell terminal (the `bash` terminal):
<p align="center">
<img src="/Module_1/images/create_new_repo8.jpeg" width="500"/>
</p>

2. One of the first steps in using `bash` is to learn how to navigate the file structure. Just like your computer, unix has a nested file structure and directory structures. One basic command is `pwd`. This gives you the **Present Working Directory**; this is the directory that you are currently in. You can see the files and folders of the directory that you are in by typing `ls`. This **lists* the files and folders. You can change the directory that you are in by typing `cd ___` where the blank is filled into the directory that you want. For instance, you can move to your scratch directory by typing `cd /scratch/aob2x`

3. Programs and commands like `ls` take arguments. For instance `ls -l` reports the output as a list. `ls -lh` reports the file sizes in human readable numbers. `ls -lha` reports all hidden files and folders too. Try `ls -lha ~/`. The files and folders starting wiht a `.` are hidden but you can still access them. They usually contain the settings for various programs.

4. You can read about the arguments for any program by typing `man ____`, e.g. `man lh`.

5. Here are some other useful `bash` programs that let you see the contents of files:<br>
• `cat`: prints out the contents of a file <br>
• `less -S`:  lets you scroll through a file. The -S flag prevents line-wrapping<br>
• `wc -l`: How many lines are in an text input stream<br>



6. You have three main places to go on Rivanna:
• `~/`: this is your home directory. How much space do you have on it?
• `/scratch/USERNAME`: this is your scratch directory. How much space do you have on it?
• `/project/biol4559-aob2x`: this is our shared project directory. Are there hidden files there? If so, what does the hidden file say?
 
