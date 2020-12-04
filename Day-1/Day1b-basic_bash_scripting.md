# Basic bash coding

You have already learned some useful commands in the last section (`pwd`, `cd`, `ls`, `scp`, `cat`, & tab to autocomplete) there are a couple more commands that will come in handy as you begin to work with more files using the command line interface (CLI). 

## Viewing the contents of files

We had previously used the `cat` command to look at the contents of our .bash_profile. The `cat` command will list the entire contents of a file to your screen, this is not a big deal when the file is small like the .bash_profile we were looking at, however this can be a little difficult to use for very large files or if you only wanted to see the first line in a file. There are a couple of other commands that enable you to view the contents of a file with a little more control. The `more` command shows you as much of the file as can be shown in the size of the terminal screen you have open, you can continue to "scroll" through the rest of the file by using the space bar which will proceed through the file page by page (a page being the size of the terminal window you have open), or the return key which will "scroll" through the file line by line. The `head` command will show you the first ten lines of a file, conversely the `tail` command shows you the last ten lines of a file. Either of these commands can be modified with the flag `-n` to show a different number of lines. 

```bash

# Show the first 20 lines of the all_counts.txt file
head -n 20 all_counts.txt

# Show the last 50 lines of the all_counts.txt file
tail -n 50 all_counts.txt

```

## Copying, renaming, and removing files 

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
You can also use the `mv` command to move a file to a new location. Let's move the alll_counts.copy.txt from your home directory into your fundamentals_of_bioinformatics directory.

```bash

# Move the all_counts.copy.txt into your fundamentals_of_bioinformatics directory
mv ~/all_counts.copy.txt fundamentals_of_bioinformatics/

#check the contents of your fundamentals_of_bioinformatics directory
ls

```

Copying the all_counts.copy.txt file was just an exercise to show you how the tools works, in practice you will want to keep your directories as neat as possible as you accumulate a lot of files in a short time when processing next-gen sequencing files. Let's remove the all_counts.copy.txt file with the `rm` command.

```bash

# Remove the all_counts.copy.txt file
rm all_counts.copy.txt

```
You will notice that before the file was deleted you were asked if you were sure you wanted this file deleted, you want to be careful not to remove files that you did not create. If you do accidentally remove a file that you realize you need later discovery uses snapshots to take yearly, monthly, and weekly records of the files in your directory and you can go into the `.snapshot/` directory to recover the file. If you only recently signed up for a discovery account you will not have much in your `.snapshot/` directory, but it's something that might come in useful to you someday.

## Manipulating file contents

There are times that you will only be interested in a subset of your data, for example, the all_counts.txt file is a counts matrix with the first column as the gene names and the next several columns list the number of reads mapping to each gene for each sample, with the sample name at the top of the column. We might be interested in pulling out the counts for only a certain subset of our samples. Let's first look at the list of samples (the first line in the file) 

```bash

# List the column names in all_counts.txt
head -n 1 all_counts.txt

```

Let's say that we would like to subset this count matrix so that we are only looking at the counts for Samples SRR1039508, SRR1039509, SRR1039510, and SRR1039511 (the first four samples). To do this we can use the `cut` command to subset our data, besides including the counts for these samples we will want to see them displayed next to the gene names (first column), so ultimately we would like to view only columns 1, 2, 3, 4, and 5 from our all_counts.txt file. 

```bash

# Look at only the counts from the first four samples
cut -f 1,2,3,4,5 all_counts.txt

```

Teh `cut` command automatically cuts on the tab character and our columns in this file happen to be tab delimited so we only need to use the command `cut` and the argument for the fields we are interestsed in keeping. Now lets say that we are interested in looking at samples SRR1039508 and SRR1039523 (the first sample and the sixteenth sample). 

```bash

# Look at only the counts from SRR1039508 and SRR1039523
cut -f 1,2,17 all_counts.txt

```

This file is pretty long, let's see how many lines this file has with the word count `wc` command using the lines `-l` argument.

```bash

# Count the lines in the all_counts.txt file
wc -l all_counts.txt

```

60,623 lines is unweild to get a feel for how counts are different between two samples that we are interested in, but we can combine our `cut` command with the `head` command that we previously learned to look at only the first 100 lines of the file and get a feel for how different our samples are using the `|` to join the commands. The `|` uses the output from the first command as the input for the second command.

```bash

# List only the first 100 lines of only samples SRR1039508 and SRR1039523
cut -f 1,2,17 all_counts.txt|head -n 100

# List only the first 100 lines of only samples SRR1039508 and SRR1039523
head -n 100 all_counts.txt| cut -f 1,2,17

```

You can see that changing the order of the commands doesn't affect the output, the output from each command is identical. Now that we have had a chance to really look at the data by capturing the first hundred lines we might decide that we need to generate a new file containing only samples SRR1039508 and SRR1039523. To do this we will use the redirect command `>` to force the output of the command into a new file rather than printing it out to the screen. 

```bash

# Print the counts from SRR1039508 and SRR1039523 to a new file 
cut -f 1,2,17 all_counts.txt > 8_23_counts.txt

# list the contents of your directory
ls

```

The redirect command comes in two flavors, the `>` version if the filename you indicate does not already exist will create a new file and add the output from your command OR if the filename does exist it will write over the existing file with the output of your command. The other version is often called the append `>>` again if the filename you indicate does not exist it will create a new file and write the output from your command, however if the filename already exists it will append the output of your command to the bottom of the existing file. We will see an example of this version used later on.  

Another way that we can append content to files is to use the command `paste`, this command will add to an existing file, or combine two files but rather than adding to the bottom of the file like the append command `>>` does `paste` uses a tab delimiter to separate the contents of two files as multiple columns. One way we might want to use this is if we wanted to add a column to the 8_23_counts.txt file we created with the information from sample SRR1039509. 

```bash

## Add a column to the 8_23_counts.txt file with the info from sample SRR1039509

#create a new file with the counts for sample SRR1039509
cut -f 3 all_counts.txt > 9_counts.txt

# Combine the contents of 8_23_counts.txt with 9_counts.txt and save to a new file 
paste 8_23_counts.txt 9_counts.txt > 8_23_9_counts.txt

# Check that the command behaved as expected
head 8_23_9_counts.txt

```


## Looking for patterns in your files with Grep

Often we will want to pull a specific piece of information from a large file, lets say that we were interested in the read counts for a specific gene, ALDH3B1, the aldehyde dehydrogenase 3 gene family B1 which plays a major role in the detoxification of aldehydes generated by alcohol metabolism and lipid peroxidation. First we know that in our all_counts.txt file the gene IDs are listed using their ensembl identifiers so we need to look for the ensembl identifier that corresponds to the ALDH3B1 gene, you can find that [here](https://uswest.ensembl.org/index.html) by selecting human in the drop down menu and typing the gene name in the field below. Once we have the ensemble ID we can use the tool `grep` to find that in our data. 

```bash

# Get the count data for ENSG00000006534 (ALDH3B1) from all_counts.txt
grep "ENSG00000006534" all_counts.txt

```

What was returned are the counts for this gene across all of your samples. `grep` is a pattern recognition tool and while it is useful to pull out pieces of information that you know it can be even more powerful to pull out pieces of information that all fit into a certain category. This is the real power of the `grep` tool, but in order to use the tool to its full potential we first need to discuss patterns, or regular expressions (regex). Regular expressions are special characters that stand for useful patterns, see the table below. 

Regex| Pattern
---|---
* | wildcard stands for any number of anything
^ | start of the line
$ | end of the line
[0-9] or \d| any number (0123456789)
[a-z]| any lowercase letter
[A-Z]| any uppercase letter
\t | a tab
\n | a newline
\s | any white space (tab, newline, space)
\S | non-white space (the opposite of \s)

These regular expressions can be used with any of the tools that you have learned thus far, so if we wanted to list all of the files in our directory that end in .txt we could use the following command.

```bash

# List all files that end in .txt
ls *.txt

```

We can even enhance the power of these regular expressions by specifying how many times we expect to see the regular expression with quantifiers.

Quantifier| Operation
---|---
X* | 0 or more repetitions of X
X+ | 1 or more repetitions of X
X? | 0 or 1 instances of X
X{*m*} | exactly *m* instances of X
X{*m*,} | at least *m* instances of X
X{*m*,*n*} | between *m* and *n* instances of X

Now lets use some of these regular expressions in a `grep` command  to see their utility. Let's use regular expressions to see how many genes have no reads expressed for the first four samples.  

```bash
# Count the number of genes with no reads in the first four samples
grep "^ENSG[0-9]*\t0\t0\t0\t0\t" all_counts.txt| wc -l

# Count the number of genes with no reads expressed in any of the samples
grep "^ENSG[0-9]*\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0$" all_counts.txt| wc -l

```


## For & while loops 

Here we are going to stretch your muscles a little we are going to switch from our counts file to a fastq file. Fastq files are the basic file format that most next-gen data are stored in their raw format. Fastq files store 
##### Left off here 
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
