# High Performance Computing

High performance computing (HPC) is a group of interconnected computers, called **nodes**, that form what is commonly called a **cluster**. An HPC system is constructed from several different types of nodes, the **head node** is the node that you will be directed to upon logging onto the cluster. The head node is the middle man between the **compute nodes** that do most of the data processing and the outside network. Data processing should not be performed on the head node as it is not designed for this purpose.

Similar to your personal computer each node is made up of **cores**, **memory**, and **disk space**. Again just as in your personal computed the **core** processes the commands that are handed to it, the **memory** is for storing information needed by the core to execute the process, and the **disk space** is for longer term storage. Most personal computers are made up of multiple cores (quadcore, 6-core) and so data processing is often **parallelized**, that is parts of the task that do not depend on each other are preformed simultaneously by different cores **threading - explain what threading is and how you can tell if a task can be parallelized** . HPCs generally have hundreds or thousands of cores and thus are optimized for parallel computing of tasks that require more memory or disk space than is generally present on a personal computer. Even if a task could be performed locally, on your personal computer, it is generally performed faster and more efficiently on an HPC system. 

<p align="center">
  <img src="../figures/HPC_PC.png" title="xxxx" alt="context"
	width="70%" height="70%" />
 </p>
 </p>
    
When you interact with an HPC system you do so using a terminal program (Mac/Linux - terminal, PC - MobaXterm, Cmder,etc.) rather than a graphics user interface (GUI). Learning to use the terminal to interact with files can be a steep learning curve, but once you get the hang of it you may find yourself interacting with your local files using the terminal! 

<p align="center">
  <img src="../figures/terminal.png" title="xxxx" alt="context"
	width="100%" height="100%" />
 </p>
 </p>

Dartmouth's HPC is maintained by [Research Computing](https://rc.dartmouth.edu)(RC) and is made up of 3000+ compute cores, 120,000 graphics processing cores, 12+ TB of memory, and 1PB of storage space. You can request an HPC account that comes with 50 GB of network file storage in [DartFS](https://rc.dartmouth.edu/index.php/dartfs/) for your personal account and 1TB of storage for a lab account. RC takes snapshots of your directories daily (kept for a week), weekly (kept for a month), and monthly(kept for a year), so that if you accidentaly delete or modify a file you have a repository of previous versions of that file. Furthmore they offer an option for web sharing, so that links to data can be shared with collaborators on or off campus (public_html).

## Discovery

`ssh netID@discovery7.dartmouth.edu`

Discovery is a linux cluster with 2016 cores, 14.5TB of memory, and 1.8PB of disk space. The compute nodes are managed by a [job scheduler](https://rc.dartmouth.edu/index.php/using-discovery/scheduling-jobs/scheduler-policy/) when you log onto discovery you are automatically directed to the head node where data processing should not be executed. Interactive processing on discovery should only be performed on the test node x01. Once you are satisfied that the commands work as you expect they can be submitted to the scheduler using a PBS script.   

To use the test node simply log onto the test node interactively `mksub -I -l node=x01 -l walltime=3:00:00` and execute the commands that you want to test. 

To submit a job to the scheduler it is best practice to denote the submission options along with the data processing commands from within a PBS script and then submit that script to the scheduler `mksub my_script.pbs`. 

An example PBS script:
```bash
#!/bin/bash -l  
#declare a name for this job to be my_serial_job  
#it is recommended that this name be kept to 16 characters or less  
#PBS -N my_job  

#request the queue (enter the possible names, if omitted, default is the default)  
#this job is going to use the default  
#PBS -q default  

#request 1 node and 16  processors on that node  
#PBS -l nodes=1:ppn=16  

#request 0 hours and 15 minutes of wall time  
#(Default is 1 hour without this directive)  
#PBS -l walltime=0:15:00  

#mail is sent to you when the job starts and when it terminates (e) or aborts (a) - can also have an email sent when the job begins (b)   
#PBS -m ea  

#specify your email address   
#PBS -M my_email@dartmouth.edu  

#By default, PBS scripts execute in your home directory, not the  
#directory from which they were submitted. The following line  
#places the job in the directory from which the job was submitted.   
cd $PBS_O_WORKDIR  

some command that I would like to submit to the scheduler  
```

Though it is possible to submit a command to the scheduler without a pbs script `mksub -N my_job -l nodes=1;ppn:16 -l walltime=0:15:00 -M my_email@dartmouth.edu -m ea my command goes here` saving the command and the submission settings the form of a PBS script enables you to easily resubmit the job if you run into an error, or edit and resubmit th job using different parameters. 

## Andes/Polaris

`ssh netID@andes.dartmouth.edu`  

`ssh netID@polaris.dartmouth.edu`

Andes and polaris are shared memory computers which run jobs that require a lot of memory or scratch space (temp files that are created during processing but discarded later). Andes has 60 cores, 512 GB of memory, and 5TB of scratch space. Polaris has 40 cores, 1TB of memory, and 5TB of scratch space. You will notice that there are far fewer cores on andes/polaris, these clusters do not use a job scheduler and jobs are executed interactively on these HPCs. If you feel more comfortable executing jobs interactively this is where you should work. 

You will also notice there is a lot more memory on polaris than andes, jobs that require a lot of memory and scratch space should be executed interactively on polaris as discovery may not have the scratch space or memory available to execute these types of jobs.

## Logging onto an hpc
Logging onto an hpc -ssh  `ssh netID@discovery7.dartmouth.edu`

creating an alias for logging on  `alias discovery="ssh netID@discovery.dartmouth.edu"`

Interacting with files on your hpc from your local computer - SSH clients (filezilla, cyberduck)

## customizing your account on discovery/polaris/andes 
modules   
List the modules you currently have loaded `module list`. 
List the modules you could potentially load `module avail`. 
Load a module from the list of available modules for the current login session you are running `module load`  
Remove a module from the list that are loaded `module rm`  
conda environments -bioconda vs. yml file   
.bash_profile - what is it and how do you use it  
