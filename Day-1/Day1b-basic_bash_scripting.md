## Basic bash coding

You have already learned some useful commands in the last section (`pwd`, `cd`, `ls`, `scp`, `cat`, & tab to autocomplete) there are a couple more commands that will come in handy as you begin to work with more files using the command line interface (CLI). 

#### Viewing the contents of files

We had previously used the `cat` command to look at the contents of our .bash_profile. The `cat` command will list the entire contents of a file to your screen, this is not a big deal when the file is small like the .bash_profile we were looking at, however this can be a little difficult to use for very large files or if you only wanted to see the first line in a file. There are a couple of other commands that enable you to view the contents of a file with a little more control. The `more` command shows you as much of the file as can be shown in the size of the terminal screen you have open, you can continue to "scroll" through the rest of the file by using the space bar which will proceed through the file page by page (a page being the size of the terminal window you have open), or the return key which will "scroll" through the file line by line. The `head` command will show you the first ten lines of a file, conversely the `tail` command shows you the last ten lines of a file. Either of these commands can be modified with the flag `-n` to show a different number of lines. 

```bash

# Show the first 20 lines of the all_counts.txt file
head -n 20 all_counts.txt

# Show the last 50 lines of the all_counts.txt file
tail -n 50 all_counts.txt

```

#### Copying and renaming files 

Sometimes you will want to copy a file into a new directory, this command is very similar to the `scp` command that we used in the previous lesson. Let's copy the all_counts.txt file from the fundamentals_of_bioinformatics directory to your home directory.

```bash

# Copy the all_counts.txt file to your home directory
cp all_counts.txt ~/all_counts.txt

```
Another useful tool for your toolbelt is renaming a file using the `mv` command. Let's rename the copy of the all_counts.txt file that we just created.

```bash

# Rename the copied all_counts.txt file
mv ~/all_counts.txt ~/all_counts.copy.txt

```
You can also use the `mv` command to move a file rather than copy a file

creating directories and organizing files -  mv, cp  
manipulating file contents - cut, paste  
piping commands to link them together

grep - looking for patterns in your files  
sed - editing batches of files  
redirecting output - >   
regex special characters - \* ,\t, \n, ., ~, ^, $  
navigating the prompt - ctrl+a/ctrl+e, ctrl+c   


**NEEDS TO BE IN DAY1A** navigating to directories - space limitations and the scratch directory  



Downloading data from external sources (NCBI, ensembl, SRA, etc.) - rsync, curl, wget  


## For & while loops 

Loops allow us repeat operations over a defined variable or set of files. Essentially, you need to tell Bash what you want to loop over, and what operation you want it to do to each item. 

Notice that the variable ***i*** set in the conditions for our loop is used to reference all the elements to be looped over in the operation using the term ***$i*** in this **for*** loop example: 
```bash 
# loop over numbers 1:10, printing them as we go
for i in {1..10}; do \
   echo "$i"; \
done
``` 


Alternatively, if you do not know how many times you might need to run a loop, using a ***while*** loop may be useful, as it will continue the loop until the boolean (logical) specified in the first line evaluates to `false`. An example would be looping over all of the files in your directory to perform a specific task. e.g. 
```bash
ls *.fastq.gz | while read x; do \
   # tell me what the shell is doing 
   echo $x is being processed...; 
   # provide an empty line for ease of viewing 
   yes '' | sed 1q;  \
   # unzip w/ zcat and print head of file
   zcat $x | head -n 4;  \
   # print 3 lines to for ease of viewing 
   yes '' | sed 3q ;
done
```

Perhaps we wanted to check how many reads contain the start codon `ATG`. We can do this by searching for matches and counting how many times it was found, and repeating this process for each sample using a for loop. 
```bash
ls *.fastq.gz | while read x; do \
   echo $x
   zcat $x | sed -n '2~4p' | head -4 | grep -o "ATG" | wc -l
done
```

We could use one of these loops to perform the nucleotide counting task that we performed on a single sample above. 
```bash
ls *.fastq.gz | while read x; do \
   yes '' | sed 1q 
   echo processing sample $x 
   zcat $x | sed -n '2~4p' | sed -n '1,10000p' | grep -o . | sort | grep 'C\|G' | uniq -c ;
done
```
## Scripting in bash 

So loops are pretty useful, but what if we wanted to make it even simpler to run. Maybe we even want to share the program we just wrote with other lab members so that they can execute it on their own FASTQ files. One way to do this would be to write this series of commands into a Bash script, that can be executed at the command line, passing the files you would like to be operated on to the script. 

To generate the script (suffix `.sh`) we could use the `nano` editor: 
```bash 
nano count_GC_content.sh
```

Add our program to the script, using a shebang `#!/bin/bash` at the top of our script to let the shell know this is a bash script. As in the loops we use the `$` to specify the input variable to the script. `$1` represents the variable that we want to be used in the first argument of the script. Here, we only need to provide the file name, so we only have 1 `$`, but if we wanted to create more variables to expand the functionality of our script, we would do this using `$2`, `$3`, etc. 
```bash 
#!/bin/bash
echo processing sample "$1"; zcat $1 | sed -n '2~4p' | sed -n '1,10000p' | grep -o . | sort | grep 'C\|G' | uniq -c
```

Now run the script, specifying the a FASTQ file as variable 1 (`$1`)
```bash
# have a quick look at our script 
cat count_GC_content.sh

# now run it with bash 
bash count_GC_content.sh SRR1039508_1.chr20.fastq.gz
```

Now we can use our while loop again to do this for all the FASTQs in our directory 
```bash
ls *.fastq.gz | while read x; do \
   bash count_GC_content.sh $x
done
```

What if we wanted to write the output into a file instead of printing to the screen? We could save the output to a *Standard output* (stout) file that we can look at, save to review later, and document our findings. The `1>>` redirects the output that would print to the screen to a file.
```bash
# create the text file you want to write to
touch stout.txt

# run the loop 
ls *.fastq.gz | while read x; do \
   bash count_GC_content.sh $x 1>> stout.txt 
done

# view the file 
cat stout.txt
```

These example programs run fairly quickly, but stringing together mutiple commands in a bash script is common and these programs take much longer to run. In these cases we might want to close our computer and go and do some other stuff while our program is running. We can do this using `nohup` which allows us to run a series of commands in the background, but disconnects the process from the shell you initally submit it through, so you are free to close this shell and the process will continue to run until completion. e.g. 
```bash
nohup bash count_GC_content.sh SRR1039508_1.chr20.fastq.gz &

# show the result 
cat nohup.out 
```
