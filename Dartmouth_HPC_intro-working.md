# High Performance Computing
what is an HPC, why is it better than working on your local computer?
What are the limitations of an hpc?

## Discovery
What are the details of this cluster?
What types of jobs are best for running on this cluster?
What is the test node- how do you use it?

## Polaris/Andes
What are the details of this cluster?
What types of jobs are best for running on this cluster?

## Logging onto an hpc
Logging onto an hpc -ssh
creating an alias for logging on
Interacting with files on your hpc from your local computer - SSH clients vs mount points

## customizing your account on discovery/polaris/andes 
modules 
conda environments -bioconda vs. yml file 
.bash_profile

## Interacting with files from your hpc
Editing files on the hpc - nano
navigating the prompt - ctrl+a/ctrl+e
figuring out where you are - pwd
creating directories and organizing files - cd, mv, cp
navigating to directories - space limitations and the scratch directory
copying files from one location to another - recursive copying, rsync, cp, scp
Downloading data from external sources (NCBI, ensembl, SRA, etc.) - rsync, curl, wget
Looping commands

## Large jobs on the hpc
building a command - reading the manual!!! 
writing a bash script by linking commands
writing a pbs script to submit a job
checking the traffic on discovery- pbsmons
submitting a job - mksub
checking the status of a submitted job- qstat, myjobs 
nohup

## Error mitigation
stack exchange - what is it and how do yuo use it?
Error messages, where can I find them and how do I know what they mean?
Problem set
