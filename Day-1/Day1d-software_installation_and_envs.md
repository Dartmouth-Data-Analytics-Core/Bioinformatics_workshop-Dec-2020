# Customizing your working environment on discovery/polaris/andes 

## Modules

Each person who uses the HPCs at Dartmouth has a different set of tasks and data that they need to work on, and as such we do not need all the same software loaded to complete the tasks that you are interested. Instead discovery/polaris/andes have modules that contain pre-loaded software that you can load into your current environment so that they are available for you to use. In order to see the modules you currently have loaded in your environment use the command `modules list`. To see the breadth of software available for you to load use the command `module avail`.

<p align="center">
  <img src="../figures/modules_avail.png" title="xxxx" alt="context"
	width="100%" height="100%" />
 </p>
 </p>

You can see that not only is there a lot of software available for you, there are multiple versions of the same software (java 1.6, 1.7, & 1.8) in case a circumstance calls for a specific version of that software. Let's add the latest version of R (R/3.3.1) to your current environment.

```bash

# Load a module to the current environment
module load R/3.3.1

# List your module to check that the module was loaded
module list
```
There may also be circumstances where you want to remove a loaded module as the software is interfering with another process that you would like to run or you want to load a different version of the same software. Let's remove the module we just loaded.

```bash

# Remove the module
module rm R/3.3.1

# Check that the module was removed
module list

```

There are some pieces of software that you will want to make sure are loaded each time that you log onto the HPC. It would be annoying if you had to use the command `module load` with each piece of software you always want loaded every time you log in. Instead you can edit a hidden file that is located in each of your home directories called `.bash_profile`. This file has preferences for setting up your environment and is executed each time that you log into the HPC. To edit a file from your command line interface we will use the tool `nano`. First let's look at the current contents of your `.bash_profile`

```bash

# Look at the contents of your .bash_profile
cat .bash_profile

```

You can see there is a location that says `# put your own module loads here` under this line is where we will add the commands for the modules that we would like to load. Use the `nano` command to add the latest R version to your environment each time you log in by adding the line `module load R/3.3.1` under the line `# put your own module loads here`. Then use **ctrl+X** to exit the `nano` program. You will be asked if you would like to save the changes you made, type `Y`, and if the changes should be saved to the file name .bash_profile hit return. Let's check that the changes we made have been saved.

```bash

#Look at the contents of your modified .bash_profile
cat .bash_profile

``` 
## Installing Bioinformatics Software Using a CLI

Bioinformatics software can be obtained in several ways.  Depending on the software to be installed, it may be available in one of the following formats:
 - Pre-installed on your system (eg. modules on Discovery)
 - Language-specific package managers (eg. R/Bioconductor, Python/pip)
 - Full package and environment management tools (eg. Conda)
 - Pre-compiled binary executable
 - Source code to be compiled
 - Virtual machine images (eg. Docker, Singularity)

## Environments
A command line environment is the collection of all programs and environment variables available to a user on a specific system.  Because all systems are different, it may be necessary to edit the environment before using certain software, or even create separate environments for using specific software. The `env` command will show all environment variables available in the current linux shell.  Some important variables are `$PATH`, `$PWD`, `$HOME`, and `$TMPDIR`.

### The $PATH environment variable
The `$PATH` environment variable stores a list of directories where programs are available.  The list is stored as strings separated by colons, so that many directories can be defined.  If the list contains multiple directories, it will be searched from left to right, until the command is found.  
```shell
echo $PATH
/usr/local/bin:/usr/lib64/qt-3.3/bin:/opt/java/latest/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/opt/bin:/opt/slurm/bin:/optnfs/moab/active/bin:/dartfs-hpc/admin/etc/cluster_files/scripts:/optnfs/mkerberos/bin:/dartfs-hpc/admin/local/bin:/dartfs-hpc/rc/home/b/f00527b/bin

#Using 'tr' to swap the colons for newlines, so a human may read the list
echo $PATH| tr ":" "\n"
/usr/local/bin
/usr/lib64/qt-3.3/bin
/opt/java/latest/bin
/usr/local/bin
/usr/bin
/usr/local/sbin
/usr/sbin
/opt/bin
/opt/slurm/bin
/optnfs/moab/active/bin
/dartfs-hpc/admin/etc/cluster_files/scripts
/optnfs/mkerberos/bin
/dartfs-hpc/admin/local/bin
/dartfs-hpc/rc/home/b/f00527b/bin
```

A command for finding where a program lives in the $PATH is the `which` command.  This can be useful for debugging environment issues, if one loses track of where a program is installed:
```shell
which echo
/usr/bin/echo
```


## Software pre-installed on the system
Linux systems will have many core utilities for navigating the file system, creating, editing and removing files, downloading and uploading files, compiling code, submitting jobs to a cluster, and many more.  These utilities are commonly found in `/usr/bin`.  On Discovery, software is also made available to users via system-wide modules.  See modules section.

## Language-specific package managers
Package manages for specific programming languages aim to make the installation of packages or libraries more simple, and from a central location.  This allows software to be installed using a single command, rather than having to search the internet for each piece of software and download/install it separately.  For R, packages are available in CRAN or in Bioconductor.  For Python, packages are available in PyPI.

To install R packages from CRAN or Bioconductor (within R shell):
```R
#Install ggplot2 from CRAN
install.packages('ggplot2')
#Install DESeq2 from Bioconductor
install.packages("BiocManager") #Get Bioconductor, if not installed already
BiocManager::install("DESeq2")
```

To install Python packages from PyPI (within bash shell):
```shell
#Install matplotlib from PyPI
pip install matplotlib
```

## Full package and environment management tools (eg. Conda)


## Pre-compiled binary executable
Some developers will pre-compile releases of their software for several operating systems and make them available for download.  If a pre-compiled executable is available for the Linux system we are using (for Discovery, this is CentOS 7), this can be a painless way to install software.  It only requires downloading the executable to a directory and running it.

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

## Virtual machine images (eg. Docker, Singularity)
Virtual machine images allow software to be distributed along with an entire linux environment.  This ensures that anyone running the software will be able to, regardless of software installed or environment variables.  However, Docker is not currently available on Discovery, as it requires administrator permissions to run.  Some images may be run with Singularity.

