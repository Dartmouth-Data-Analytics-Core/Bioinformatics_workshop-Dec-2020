# Welcome to the DAC Fundamentals of Bioinformatics workshop #

Before you attend the workshop there are a couple of things we would like you to do to get setup for using the tools that will be required during the course of the workshop.  

For those of you that indicated that you did not have an account on *discovery* you should have received an email from me explaing how to set that up, please make sure this is done and you are able to log into your account **BEFORE** the workshop begins. 

---

## Downloading the data ##

The commands that you will be following can be found in markdown `(.md)` files where there is a brief description of the command and how it is applied to the data and what it does followed by an example command that you can copy and paste into the terminal window. The majority of day 1 and 2 will be using the terminal window on your local machine, with an open `ssh` connection to discovery7, as we will be running `bash` code. For some of day 2 and most of day 3 you will be using RStudio on your local machine to run the commands in the  markdown files (`.md`) located in this GitHub repo. 

To start make sure that you are able to use a terminal emulator, select one of the following based on your operating system, download it and open it up. 

Operating system| Terminal emulators
---|---
Mac| Terminal (comes pre-installed)
Windows| [MobaXterm](https://mobaxterm.mobatek.net/download.html) <br> [PuTTY](https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html)
Linux| Konsole, Terminal, etc. (should be pre-installed but depends on the desktop environment you are running)


In your terminal window navigate to where you want to download the files needed for this workshop onto your local machine. Then execute the following command:

```bash
git clone https://github.com/Dartmouth-Data-Analytics-Core/Bioinformatics_workshop_setup/
```

**On Monday** before you log onto the first zoom session we will make the workshop materials public and you should download those to your local machine (prefereably in the same location as you downloaded the setup materials) with the folloing command: 

```bash
git clone https://github.com/Dartmouth-Data-Analytics-Core/Bioinformatics_workshop/
```

---

## Install the Integrative Genomics Viewer (IGV)

We will be using the [Integrative Genomics Viewer (IGV)](http://software.broadinstitute.org/software/igv/), a genome browser produced by researchers at the Broad Institute, to explore and visualize genomics data. 

<img src="figures/igv.png" height="100" width="100"/>

You will need to download and install the IGV desktop application for your operating system before the workshop begins. The latest versions of IGV can be found at their [downloads page](http://software.broadinstitute.org/software/igv/download). After installing IGV, try opening the application on your computer to confirm the installation was successful. 

---

## Setting up a Conda Environment ## 

Conda is a package management system that helps you find, install, and organize groups of packages needed for a particular task. Conda environments are really useful when working on high performance computing (HPC) environments like Dartmouth's Discovery system because you can install packages locally without needing administrator permission. Conda environments are also useful for project continuity, the versions of the packages that you install in your environment and all of their dependencies will remain the same (unless you update them). We will be using a conda environment to make sure we all have the same version of many different bioinformatics software programs available to us. 

Before you begin using conda environments on discovery you will need to ensure that you have run the source code to enable the conda commands. Log into discovery and run the following command:

```bash
source /optnfs/common/miniconda3/etc/profile.d/conda.sh
```

We recommend that you add the above line of code to your `.bashrc` file in your home directory, otherwise you will need to run this command each time you start a new session on discovery. To do this use the nano text editor on discovery to copy the line above into your `.bashrc` file.

```bash
# open the file with the text editor
nano .bashrc

# copy this line to the file : source /optnfs/common/miniconda3/etc/profile.d/conda.sh

# use ctrl+x to exit the editor and save the changes you made
```


Next you will have to run the following command to create a .conda/ directory in your home drive to store all of your personal conda environments. You only have to run this command once to make this directory, so it does not need to be added to your .bashrc file.

```bash
cd ~
mkdir -p .conda/pkgs/cache .conda/envs
```

Now you will need to create the conda environment that we will be using for the course into your personal directory of accessible conda environments. This takes about 15 minutes to execute and you will see all of the packages that are loaded into this environment. The number of packages should indicate why conda environments are so useful, imagine having to load all of these packages individually it is much easier to load them with a single command in a conda environment.

```bash
conda env create -f /scratch/fund_of_bioinfo/environment.yml
```

When you are ready activate the conda environment, which you will need for the work we are doing for days 1 and 2 of the workshop you can use the following command. 

```bash
conda activate bioinfo/
```

You will see that the activate command has worked when it reads (fund_of_bioinfo) rather than (base) to the left of the prompt. When you are finished using a conda environment it is good practice to deactivate your session with the following command.

```bash
conda deactivate
```

---
## Installing an SFTP client ##

This is optional but for those of you that are new to the command line this might be an easier way to move files between the HPC environment and your local machine. An SFTP client stands for secure file transfer protocol and will enable you to drag and drop files as you might in a finder window between your local machine and a remote location. I use FileZilla, which I believe works on Mac, Windows, and linux operating systems. You can download [FileZilla](https://filezilla-project.org/download.php?show_all=1) by following the link and selecting the version that is correct for your OS, then open the program to ensure that you have downloaded it successfully. 

---

## Setting up an R project ##

We will be using R-Studio to explore and analyze genomics data on day 2 and 3, therefore we ask that you have R and R-Studio installed prior to attending the workshop. You will need to be running at least version 3.6 to ensure all of the packages needed will run smoothly. The latest versions for R and R-Studio can be found [here](https://cran.r-project.org) and [here](https://rstudio.com/products/rstudio/download/).

Next you will need to set up a new project in R-Studio for this workshop. Projects in R are like containers for various jobs that you will perform, the history of the project will be loaded when you open a new project. By using a project you can install all of the tools ahead of time and they will be there for you when you need them during the workshop. In R-Studio under File select New directory and then select New Project and name the project something you will remember (bioinfo_workshop).

Now that you have your project loaded, run the following code to install all of the packages we will be using during the workshop: 
```r
if (!any(rownames(installed.packages()) == "biomaRt")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("biomaRt")
}
library(biomaRt)

if (!any(rownames(installed.packages()) == "IRanges")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("IRanges")
}
library(IRanges)

if (!any(rownames(installed.packages()) == "GenomicRanges")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("GenomicRanges")
}
library(biomaRt)

if (!any(rownames(installed.packages()) == "GViz")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("GViz")
}
library(GViz)

if (!any(rownames(installed.packages()) == "org.Hs.eg.db")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("org.Hs.eg.db")
}
library(org.Hs.eg.db)

if (!any(rownames(installed.packages()) == "EnsDb.Hsapiens.v86")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("EnsDb.Hsapiens.v86")
}
library(EnsDb.Hsapiens.v86)

if (!any(rownames(installed.packages()) == "GenomicFeatures")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("GenomicFeatures")
}
library(GenomicFeatures)

if (!any(rownames(installed.packages()) == "VariantAnnotation")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("VariantAnnotation")
}
library(VariantAnnotation)

if (!any(rownames(installed.packages()) == "TxDb.Hsapiens.UCSC.hg38.knownGene")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
}
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

if (!any(rownames(installed.packages()) == "TxDb.Mmusculus.UCSC.mm10.knownGene")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
}
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

if (!any(rownames(installed.packages()) == "BSgenome")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("BSgenome")
}
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

if (!any(rownames(installed.packages()) == "ChIPseeker")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("ChIPseeker")
}
library(ChIPseeker)

sessionInfo()
```

If you have issues with any part of the installation and setup, please reach out to us directly (contact details are in the workshop ReadMe page) or come to bioinformatics help hours **December 11, 2020 1-2PM** the link is in your welcome email. 

