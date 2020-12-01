
## R on the hpc 
submitting an R script
loading R libraries 
running R interactively 

## Large jobs on the hpc
building a command - reading the manual!!!   
checking the traffic on discovery- pbsmons  
submitting a job - mksub  
checking the status of a submitted job- qstat, myjobs  

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
 


## Error mitigation  
errors with job submission - exit status what does this tell you about what happened, look in your error/output file - does it look like you would expect?
stack exchange - what is it and how do yuo use it?  
Error messages, where can I find them and how do I know what they mean?  
Problem set  
