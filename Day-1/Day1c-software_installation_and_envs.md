# Installing and managing bioinformatics software  

## Introduction

Bioinformatics software can be installed and managed in a number of ways. 

Important to be able to keep track of software verisons that you are using, and which versions where used for a particular analysis (for example, when you go to publish, you should include version numbers of critical software programs). 

In this lesson, we will introduce the major ways bioinformatics packages can be installed and managed, and provide some recommendations for how to use them in your own work. It is useful to note that other ways of managing bioinformatics software that are not discussed here do exist, however we think these are the most important ones to know when starting out.  

---






## Installing Bioinformatics Software Using a CLI

Bioinformatics software can be obtained in several ways.  Depending on the software to be installed, it may be available in one of the following formats:
 - Pre-installed on your system (eg. modules on Discovery)
 - Language-specific package managers (eg. R/Bioconductor, Python/pip)
 - Full package and environment management tools (eg. Conda)
 - Pre-compiled binary executable
 - Source code to be compiled
 - Virtual machine images (eg. Docker, Singularity)





## Software pre-installed on the system
Linux systems will have many core utilities for navigating the file system, creating, editing and removing files, downloading and uploading files, compiling code, submitting jobs to a cluster, and many more.  These utilities are commonly found in `/usr/bin`.  On Discovery, software is also made available to users via system-wide modules.  See modules section.

## Language-specific package managers
Package managers for specific programming languages aim to make the installation of packages or libraries more simple, and from a central location. This allows software to be installed using a single command, rather than having to search the internet for each piece of software and download/install it separately. For R, packages are available from two major sources:  
- [*CRAN*](https://cran.r-project.org/web/packages/available_packages_by_name.html) - A large diverse collection of R packages currently approaching 17,000 in total
- [*Bioconductor*](https://www.bioconductor.org/) - a specific collection of packages specifically geared toward facilitating bioinformatic data analysis in R

To install R packages from CRAN (within R): 
```R
# Install ggplot2 from CRAN
install.packages('ggplot2')
```

To install R packages from Bioconductor (within R): 
```
# Install DESeq2 from Bioconductor
install.packages("BiocManager") # Get Bioconductor, if not installed already
BiocManager::install("DESeq2")
```

In Python, packages are available in PyPI. To install Python packages from PyPI (within bash shell):
```shell
# Install matplotlib from PyPI
pip install matplotlib
```

## Full package and environment management tools (eg. Conda)



```bash
conda create -n python=3.7.1
```



```bash
conda activate x
```

![](../figures/conda.png)

```bash
conda deactivate
```

# channels 


![](../figures/conda_screenshot.png)





## Pre-compiled binary executable
Some developers will pre-compile releases of their software for several operating systems and make them available for download.  If a pre-compiled executable is available for the Linux system we are using (for Discovery, this is CentOS 7), this can be a painless way to install software. It only requires downloading the executable to a directory and running it.

Programs written in Java are frequently distributed as JAR files, which are similar to pre-compiled binaries, in that only a single file is required to download and install the software.  The JAR file is then run using the `java -jar` command.
```shell
java -jar some_software.jar
```

## Source code to be compiled 
If software is not available via a package manager, or via a pre-compiled executable, it must be compiled from source code.  For Bioinformatics software, this will usually be C or C++ code, and will be distributed with a "makefile", which can be compiled with the following commands.  The `--prefix="/path/to/install"` defines the directory where the software will be installed.
```shell
./configure --prefix="/path/to/install"
make
make install
```







## Environment Modules


ised on HPC systems to make software versions available 
usually installed by the system admin, so you cannot add or modify them yourself 




Each person who uses the HPCs at Dartmouth has a different set of tasks and data that they need to work on, and as such we do not need all the same software loaded to complete the tasks that you are interested. Instead discovery/polaris/andes have modules that contain pre-loaded software that you can load into your current environment so that they are available for you to use. In order to see the modules you currently have loaded in your environment use the command `modules list`. To see the breadth of software available for you to load use the command `module avail`.

<p align="center">
  <img src="../figures/modules_avail.png" title="xxxx" alt="context"
	width="100%" height="100%" />
 </p>
 </p>

You can see that not only is there a lot of software available for you, there are multiple versions of the same software (java 1.6, 1.7, & 1.8) in case you need a newer or older version of a specific program. For exmaple, let's suppose we need to use an older version of R than the version loaded by default on discovery (3.6.0). 

First, lets start running R interactively on discovery, and confirm the version currently loaded. 
```bash 
# start running R interactively
R

# (from within R) check the version
R.version

# quite running R interactively to return to the bash terminal 
q()
```

Now go ahead and load the module for the older version of R that you need (3.3.1). 
```bash
# Load a module to the current environment
module load R/3.3.1

# List your currently loaded modules to check that the module was loaded
module list
```

There may also be circumstances where you want to remove a loaded module as the software is interfering with another process that you would like to run or you want to load a different version of the same software. Let's remove the module we just loaded.
```bash

# remove the module
module rm R/3.3.1

# confirm that the module was removed
module list
```

There are some pieces of software that you will want to make sure are loaded each time that you log onto the HPC. It would be annoying if you had to use the command `module load` with each piece of software you always want loaded every time you log in. Instead you can edit a hidden file that is located in each of your home directories called `.bash_profile`. This file has preferences for setting up your environment and is executed each time that you log into the HPC. To edit a file from your command line interface we will use the tool `nano`. First let's look at the current contents of your `.bash_profile`

```bash

# Look at the contents of your .bash_profile
cat .bash_profile

```

You can see there is a location that says `# put your own module loads here` under this line is where we will add the commands for the modules that we would like to load. Use the `nano` command to add the latest R version to your environment each time you log in by adding the line `module load R/3.3.1` under the line `# put your own module loads here`. Then use **ctrl+X** to exit the `nano` program. You will be asked if you would like to save the changes you made, type `Y`, and if the changes should be saved to the file name .bash_profile hit return. Let's check that the changes we made have been saved.

```bash
# Look at the contents of your modified .bash_profile
cat .bash_profile
``` 



remove all of your modules using purge command 
```bash
module purge
```






## Virtual machine images (eg. Docker, Singularity)
Virtual machine images allow software to be distributed along with an entire linux environment.  This ensures that anyone running the software will be able to, regardless of software installed or environment variables.  However, Docker is not currently available on Discovery, as it requires administrator permissions to run.  Some images may be run with Singularity.

![](../figures/containers.png)
