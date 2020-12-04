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

## Moving files between computers 

Sometimes you will want to move files from your account on the cluster to your local computer, or visa versa. You can do this is with the commands `scp` and `rsync`. The `scp` and `rsync` command have a similar syntax, for simplicity let's focus on the scp command. The `scp` command takes two arguments each will be a path to copy from/to and the order does not matter (though in `rsync` the source comes first and the destination comes second - so lest stick with that syntax here). 

`scp source_path destination_path`

Let's use these commands to move the file `all_counts.txt` which should be in the Day1 directory you copied from the github repository, into your `fundamentals_of_bioinformatics` directory on discovery. To start open a new terminal window, this window should default to your local home directory (you can check this with `pwd`), navigate to the location you copied the github repo to so that you can see the `all_counts.txt` file when you use the `ls` command. From this new terminal window enter the following command:

```bash

# Move the file from the cluster to your local machine
scp all_conuts.txt netID@discovery7.dartmouth.edu:/dartfs-hpc/rc/home/h/netID/fundamentals_of_bioinformatics/ 

```

You will be prompted for your password to ensure you have permissions to access the file on discovery, and you should see the `all_counts.txt` file in your fundamentals_of_bioinformatics directory (you can check with the `ls` command).

Now lets copy the `fundamentals_of_bioinformatics` directory to your local directory, using the `scp` command.

```bash

# Move the file from the cluster to your local machine
scp -r netID@discovery7.dartmouth.edu:/dartfs-hpc/rc/home/h/netID/fundamentals_of_bioinformatics/ ./

```

You will notice that we used the `-r` option with this command this stands for recursive, when copying a directory this will enable you to copy that directory and all of its contents organized exactly as they are in the location you are copying from. You will also notice I used the `./` to indicate the destination, his is shorthand to mean the directory that I am currently in.

#### SFTP clients (FileZilla, CyberDuck, WinSCP)

Another way to do this that can be a little less daunting when you are new to the command line is to use an SFTP client. I use [FileZilla](https://filezilla-project.org) for this but there are many other programs available to visualize the files that you have locally and on a remote site of your choosing. 

You will notice on the left is a list of the files in my local directory and on the right is a list of files in my remote directory. Once we have a connection established double-clicking a filename will transfer the file from one location to the other. If the filename already exists in the receiving directory a warning will pop up asking if you would really like to replace the current version of the file. In order to make sure the directories in both locations are showing the most current version of their contents it is often useful to use the **Refresh** button on the top menu. 

<p align="center">
  <img src="../figures/Intro_SFTP_client.png" title="xxxx" alt="context"
	width="100%" height="100%" />
 </p>
 </p>

To establish a connection to a new site you would click on the **Site Manager** image on the top left of the screen. From there you would select the **SFTP** option under the protocol drop down menu, then enter the address of the remote host you would like to connect to, your user name (netID in this case), and password. Click connect and you should see a list displayed on the right pane of the window showing the contents of your home directory. 

<p align="center">
  <img src="../figures/logging_in2.png" title="xxxx" alt="context"
	width="100%" height="100%" />
 </p>
 </p>



## TO DO: 
Storage - 50GB for personal accounts and 1TB for lab accounts
Scratch drive



