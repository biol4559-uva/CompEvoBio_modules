# **GitHub assignment instructions**

## Overview
Once you have created a GitHub account and shared your username with me, you will be added to a GitHub Organization called [biol4559-uva](https://github.com/biol4559-uva). As you have already seen, this organization contains repositories such as the (CompEvoBio_modules)[https://github.com/biol4559-uva/CompEvoBio_modules] which contains scripts that we will use scripts that we will use throughout the course.

As a member of this GitHub Organization, you can create your own repositories. In today's class, you will create your own repository and use it to store scripts that you write for assignments. This will be a private repository that you own and that you share with me and the TA. Down the road, you might share your repository with your team during the analysis of the data-set that you are working on. Using GitHub will allow us to more easily share code and also to get feedback and help when things go awry.

Objective:
1. Learn how to use Github for cloning, pushing and pulling

Goals:
1. Creat an `ssh-key`
2. **To complete the GitHub assignment today, you will be upload your figure and script for your HWE assignment to a repository that you have made***

## Steps
1. Create a new repository on the [biol4559-uva](https://github.com/biol4559-uva) organization. To create a repository, navigate to the [biol4559-uva](https://github.com/biol4559-uva) and click the "New" button.

2. Use your UVA computing ID as the repository-name and put your full name in the Description. Select "Private", and leave the rest as is. Click "Create Repo"
<p align="center">
<img src="/Module_1/images/create_new_repo1.jpeg" width="1000"/>
</p>


3. Once your repository is created, click on it. Then click on "Manage Access". You might be asked for your GitHub password again. Click "Add People" and type in Alan's id ("alanbergland") and Robert's id (________). Give us "Write" privileges.

4. Next you will need to set up an `ssh-key`. This is a way that GitHub can verify your identity when you are pushing and pulling scripts from your computer or Rivanna. You will be setting up your `ssh-key` on Rivanna, although you might also want to set it up on your own computer. Follow the steps listed [here](https://jdblischak.github.io/2014-09-18-chicago/novice/git/05-sshkeys.html). You will be using OpenOnDemand's shell access to complete these steps:
<p align="center">
<img src="/Module_1/images/create_new_repo8.jpeg" width="500"/>
</p>


5. Once you have generated your `ssh-key` clone your repository onto Rivanna. There are different ways to clone your repository. One way is using the command line:
```git clone git@github.com:biol4559-uva/aob-test.git``` The other is using a program like RStudio or GitHub desktop. Today, we are going to try using RStudio. Using OpenOnDemand, launch an RStudio job. Go to the File menu and click "New Project". Click on Version Control:
<p align="center">
<img src="/Module_1/images/create_new_repo2.jpeg" width="500"/>
</p>
<br>
And then on Git:
<p align="center">
<img src="/Module_1/images/create_new_repo3.jpeg" width="500"/>
</p>
<br>
And then enter the name of your repository and the directory where it belongs on Rivanna:
<p align="center">
<img src="/Module_1/images/create_new_repo4.jpeg" width="500"/>
</p>
<br>
You'll need to enter the password that you set when you generated an `ssh-key`.
<p align="center">
<img src="/Module_1/images/create_new_repo5.jpeg" width="500"/>
</p>
<br>

6. From Rstudio, you will save your HWE script in your repository. Please keep things tidy. So, for example, create a directory with the name `Module_1_HWE`, and save your HWE script and PDF figure there.
<p align="center">
<img src="/Module_1/images/create_new_repo6.jpeg" width="500"/>
</p>
<br>

7. Next, you will push your script and figure to GitHub. Click on on the `Git` tab, and then click the `Staged` boxes for every file. Then click commit and follo the instrucitons.
<p align="center">
<img src="/Module_1/images/create_new_repo7.jpeg" width="500"/>
</p>
<br>

8. There is a lot more to Git than that, and some of those other aspects might come up later on in the course. If you are curious, check out this lesson on Git from [Software Carpentry](https://swcarpentry.github.io/git-novice/)
