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

Discovery is a linux cluster with 2016 cores, 14.5TB of memory, and 1.8PB of disk space. The compute nodes are managed by a [job scheduler](https://rc.dartmouth.edu/index.php/using-discovery/scheduling-jobs/scheduler-policy/) when you log onto discovery you are automatically directed to the head node where data processing should **not be executed**. Let's log onto the discovery cluster now. We will use a secure shell command `ssh` to log onto the discovery cluster. 

```bash

# Establish the secure shell connection
ssh netID@discovery7.dartmouth.edu

# Enter your password at the prompt (when you type no characters will show up to preserve privacy)
netID@discovery7.dartmouth.edu's password:

# You're in! 
(base) [netID@discovery7 ~]$

```
It is always useful to orient yourself when you're working on an HPC so that you know where the output of all of the commands you run will end up. Lets run our first command to get your location. 

```bash

# Check your location on the cluster
pwd

```

You should see something like `/dartfs-hpc/rc/home/h/netID` displayed in response to your command. Initially when you log on you will always be directed to your home directory (the address or path listed above). It is a good idea when working on projects on an HPC to stay organized so lets start by making a folder, or directory to store all of the work you do today we will call it `fundamentals_of_bioinformatics`. You will notice that I chose a title that has no spaces in it, this is because the space is a special character, special characters need to be *escaped* with the `\` and so `funadmentals_of_bioinformatics` would look like `fundamentals\ of\ bioinformatics` with the escape characters. You can see that file names with spaces become unweildy to type out so most programmers will replace spaces with `_`, `.`, or `-` in their filenames to keep everything neat.

```bash

# Create a directory 
mkdir fundamentals_of_bioinformatics

# enter the directory
cd fundamentals_of_bioinformatics

# check your location on the cluster
pwd 

# list the contents of your directory
ls

```
You can see that you have now created a directroy and entered it such that your path has an extended to include your new directory as part of your location on the cluster. You can also see that the new directory that you created is empty there are no files. Lets quickly navigate back to our home directory and check the contents of the home directory.

```bash

# navigate to the home directory
cd ~

# Check the contents of the home directory
ls

```

The `~` stands for the path you saw when we typed `pwd` in your home directory and you can navigate to the home directory quickly from anywhere else on the cluster using that command. Another way we could have navigated to the home directory would be to go up one level with the command `cd ../` since the home directory was only one level up in our path from the `fundamentals_of_bioinformatics` directory where we were. 

Especially if you type quickly typos and long filenames can trip you up, a good tip to remember is that you can use the `tab` key to autofinish a filename as long as it is unique. If we had two directories one calles `fundamentals_of_bioinformatics` and one called `fundamentals_of_chess` after you typed `fund` and used the tab key you would see `fundamentals_of_` with the prompt asking you to type more so that the correct directory can be entered. If you then typed `fundamentals_of_b` and hit the tab the rest of the name would finish for you. This seems like a little thing but it can be a real time saver if you have many hundreds of similarly named files with long complicated names. 

## Working interactively

Interactive processing on discovery should only be performed on the *test node x01*. Once you are satisfied that the commands you're testing work as you expect they can be submitted to the scheduler using a PBS script. To use the test node simply log onto the test node interactively and execute the commands that you want to test. 

```bash

# Log onto to the test node x01 interactively (the -I flag) for up to 3 hours, the -l flags indicate the details of the submission that you are asking for 
mksub -I -l nodes=x01 -l walltime=3:00:00 

```
#### Andes/Polaris

Alternatively you can work interactively on Andes/Polaris clusters. Lets take a minute to log onto Polaris. 

```bash

# Exit discovery
exit

# Ssh onto Polaris, notice the command is the same except that the location that you are logging into has changed from discovery7 to polaris
ssh netID@polaris.dartmouth.edu

# Enter your password at the prompt (the password is the same as the one you used for discovery)
netID@polaris.dartmouth.edu's password:

# You're in! 
(base) [netID@polaris ~]$

#check your location on polaris
pwd

```
You will notice the the inital location you are logged into on polaris is the same home directory you found yourself in after your initial login to discovery. You can even see that all of the same files are there with the `ls` command. The content of your directories does not change whether you are logged onto discovery or polaris, what changes are the capabilities of the computing nodes that you are hosted on. 

Andes and polaris are shared memory computers which run jobs that require a lot of memory or scratch space (temp files that are created during processing but discarded later). Andes has 60 cores, 512 GB of memory, and 5TB of scratch space. Polaris has 40 cores, 1TB of memory, and 5TB of scratch space. You will notice that there are far fewer cores on andes/polaris, these clusters do not use a job scheduler and jobs are executed interactively on these HPCs. If you feel more comfortable executing jobs interactively this is where you should work. 

You will also notice there is a lot more memory on polaris than andes, jobs that require a lot of memory and scratch space should be executed interactively on polaris as discovery may not have the scratch space or memory available to execute these types of jobs.

## SSH clients (filezilla, cyberduck)

Sometimes you will want to move files from your account on the cluster to your local computer, an SSH client is useful for this. 

## Submitting a job to the cluster

For many tasks that you will want to run the task will require lots of memory and perhaps even threading to spread parts of the job onto multiple cores to get the job completed more efficiently (think of a car being built in a factory it makes sense to build the headlights and doors in separate locations at the same time). We do this by submitting a job to the scheduler, in your submission you will tell the scheduler the name you want to give your job, the number and types of nodes you need, the time you think the job should take, where to send notification when the job begins/ends or runs into an error, and what directory the output files from the job should be placed in. It is possible to do this all in one command but often the best practice is to use a script with all of this information endcoded and submit that script to the scheduler. In this way you have a written record of the parameters you used to submit the job and if anything goes wrong you can easily modify the script and resubmit it rather than typing the whole command out again. 

Let's take a look at an example of a job submission script below: 

```bash
#!/bin/bash -l  
#declare a name for this job to be my_serial_job  
#it is recommended that this name be kept to 16 characters or less  
#PBS -N my_job  

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

The lines that start with `#PBS` are the lines that are denoting the settings we would like for the scheduler to use for our command.

- `#PBS -N` declares the name that you would like to use for your job
- `#PBS -l nodes=1:ppn=16` tells the scheduler you need 1 node with at least 16 cores to run this job
- `#PBS -l walltime=0:15:00` tells the scheduler you need a maximum of 15 minutes to complete this job, the job will be killed after 15 minutes so it is always best to err on the side of more time 
- `#PBS -m bea` tells the scheduler that you would like to be notified when the job starts, ends, and if any errors are encountered while running the job
- `#PBS -M` tells the scheduler how to contact you with the information you requested in the `-m` argument
- `cd $PBS_O_WORKDIR` tells the scheduler to put outfiles, log files and error files in the directory that you submitted your code from. Leaving this out will cause these files to be placed in your home directory. 





## customizing your account on discovery/polaris/andes 
modules   
List the modules you currently have loaded `module list`. 
List the modules you could potentially load `module avail`. 
Load a module from the list of available modules for the current login session you are running `module load`  
Remove a module from the list that are loaded `module rm`  
conda environments -bioconda vs. yml file   
.bash_profile - what is it and how do you use it  
