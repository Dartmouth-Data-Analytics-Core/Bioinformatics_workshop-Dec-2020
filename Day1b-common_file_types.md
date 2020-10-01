## FASTQ file format

FASTQ files are arguably the workhorse format of bioinformatics. FASTQs are used to store sequence reads generated in next-generatoon sequencing (NGS) experiments. Similarly to FASTA files, FASTQ files contain a herder line, followed by the sequence read, however individual quality of base calls from the sequencer are included for each record in a FASTQ file. 

Here is what a the first record of an example FASTQ file looks like
```
@SRR1039508.1 HWI-ST177:290:C0TECACXX:1:1101:1225:2130 length=63
CATTGCTGATACCAANNNNNNNNGCATTCCTCAAGGTCTTCCTCCTTCCCTTACGGAATTACA
+
HJJJJJJJJJJJJJJ########00?GHIJJJJJJJIJJJJJJJJJJJJJJJJJHHHFFFFFD
```

**Four rows exist for each record in a FASTQ file:**
- **Row 1:** Header line that stores information about the read (always starts with an `@`), such as the *instrument ID*, *flowcell ID*, *lane on flowcell*, *file number*, *cluster coordinates*, *sample barcode*, etc.
- **Row 2:** The sequence of bases called
- **Row 3:** Usually just a `+` and sometimes followed by the read information in line 1
- **Row 4:** Individual base qualities (must be same length as line 2)

Quality scores, also known as **Phred scores**, in row 4 represent the probability that the associated base call is incorrect, which are defined by the below formula for current Illumina machines:
```
Q = -10 x log10(P), where Q = base quality, P = probability of incorrect base call
```
or 
```
P = 10^-Q/10
```

Intuitively, this means that a base with a Phred score of `10` has a `1 in 10` chance of being an incorrectly called base, or *90%*. Likewise, a score of `20` has a `1 in 100` chance (99% accuracy), `30` a `1 in 1000` chance (99.9%) and `40` a `1 in 10,000` chance (99.99%). 

However, we can clearly see that these are not probabilities. Instead, quality scores are encoded by a character that is associated with an *ASCII (American Standard Code for Information Interchange)* characters. *ASCII* codes provide a convenient way of representing a number with a character. In FASTQ files, Q-score is linked to a specific ASCII character by **adding 33 to the Phred-score**, and matching the resulting number with its *ASCII* character according to the standard code. You can see the full table used for ASCII character to Phred-score conversion [here](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm). The reason for doing it this way is so that quality scores only take up **1 byte per value** in the FASTQ file.

For example, the first base call in our sequence example above, the `C` has a quality score encoded by an `H`, which corresponds to a Q-score of 39, meaning this is a good quality base call. 

Generally, you can see this would be a good quality read if not for the strech of `#`s indicating a Q-score of 2. Looking at the FASTQ record, you can see these correspond to a string of `N` calls, which are bases that the sequencer was not able to make a base call for. Streches of Ns' are generally not useful for your analysis. 

**Paired-end reads:**  

If you sequenced paired-end reads, you will have two FASTQ files:  
**..._R1.fastq** - contains the forward reads  
**..._R2.fastq**- contains the reverse reads  

Most downstream analysis tools will recognize that such files are paired-end, and the reads in the forward file correspond to the reads in the reverse file, although you often have to specify the names of both files to these tools. 

It is critical that the R1 and R2 files have **the same number of records in both files**. If one has more records than the other, which can sometimes happen if there was an issue in the demultiplexing process, you will experience problems using these files as paired-end reads in downstream analyses. 

## Working with FASTQ files 

### Basic operations 

While you don't normally need to go looking within an individual FASTQ file, it is very important to be able to manipulate FASTQ files if you are going to be doing more involved bioinformatics. There are a lot of operations we can do with a FASTQ file to gain more information about our experiment, and being able to interact with FASTQ files can be useful for troubleshooting problems that might come up in your analyses. 

Due to their large size, we often perform gzip copmpression of FASTQ files so that they take up less space, however this means we have to unzip them if we want to look inside them and perform operations on them. We can do this with the `zcat` command and a pipe (|). The pipe command is a way of linking commands, the pipe sends the output from the first command to the second command. `zcat` lists the contents of a zipped file to your screen, and head limits the output to the first ten lines. 

Lets use `zcat` and `head` to have a look at the first few records in our FASTQ file. 
```bash
# unzip and view first few lines of FASTQ file 
zcat SRR1039508_1.chr20.fastq.gz | head
zcat SRR1039508_2.chr20.fastq.gz | head
```

How many records do we have in total? (don't forget to divide by 4..) 
```bash
zcat SRR1039508_1.chr20.fastq.gz | wc -l
zcat SRR1039508_2.chr20.fastq.gz | wc -l
```
Paired-end reads should have the same number of records! 

What if we want to count how many unique barcodes exist in the FASTQ file. To do this, we would need to print all the sequence lines of each FASTQ entry, then search those for the barcode by specifying a regular expression. To print all the sequence lines (2nd line) of each FASTQ entry, we can use a command called ***sed***, short for ***stream editor*** which allows you to streamline edits to text that are redirected to the command. You can find a tutorial on using **sed** [here](https://www.digitalocean.com/community/tutorials/the-basics-of-using-the-sed-stream-editor-to-manipulate-text-in-linux). 

First we can use sed with with the `'p'` argument to tell it that we want the output to be printed, and the `-n` option to tell sed we want to suppress automatic printing (so we don't get the results printed 2x). Piping this to `head` we can get the first line of the first 10 options in the FASTQ file (the header line). We specify `'1-4p'` as we want sed tp *print 1 line, then skip forward 4*. 
```bash
zcat SRR1039508_1.chr20.fastq.gz | sed -n '1~4p' | head -10
```

Using this same approach, we can print the second line for the first 10000 entires of the FASTQ file, and use the ***grep*** command to search for regular expressions in the output. Using the `-o` option for grep, we tell the command that we want it to print lines that match the character string. 
```bash
# print the first 10 lines to confirm we are getting bthe sequence lines 
zcat SRR1039508_1.chr20.fastq.gz | sed -n '2~4p' | head -10

# pipe the sequence line from the first 10000 FASTQ records to grep to search for our (pretend) adapter sequence
zcat SRR1039508_1.chr20.fastq.gz | sed -n '2~4p' | head -10000 | grep -o "ATGGGA"
```

This is a bit much to count by each, so lets count the how many lines were printed by grep using the ***wc*** (word count) command with the `-l` option specified for lines.
```bash
zcat SRR1039508_1.chr20.fastq.gz | sed -n '2~4p' | head -10000 | grep -o "ATGGGA" | wc -l
```

Using a similar approach, we could count up all of the instances of individual DNA bases (C,G) called by the sequencer in this sample. Here we use the ***sort*** command to sort the bases printed by grep, grep again to just get the bases we are interested in, then using the ***uniq*** command with the `-c` option to count up the unique elements. 
```bash
zcat SRR1039508_1.chr20.fastq.gz | sed -n '2~4p' | head -10000 | grep -o . | sort | grep 'C\|G' | uniq -c 
```

Now we have the number of each nuleotide across the reads from the first 10000 records. A quick and easy program to get GC content. GC content is used in basic quality control of sequence from FASTQs to check for potential contamination of the sequencing library. We just used this code to check 1 sample, but what if we want to know for our 4 samples?

## SAM/BAM/CRAM file format

SAM format is a common way of representing sequenced reads, especially after reads have been aligned or mapped to a reference genome. BAM is a binary (non-human readable) format of SAM that takes up less space. CRAM files are compressed versions of BAM files - these take up the least space and are recommended for longer term storage of alignment files. SAM/BAM/CRAM files can be converted back and forth with the tool **samtools view**.

```bash

#convert sam to bam
samtools view -b sample.sam > sample.bam

#convert bam to cram
samtools view -C sample.bam > sample.cram
```

A SAM file is made up of two basic parts, the header and the alignment. All header lines will start with the **@** symbol. The mandatory flag **@HD** will come first in the file and should only occur once, this flag has the meta-data that pertains to the SAM file and will either have a **GO** field indicating that reads are grouped but not sorted or a **SO** field indicating that reads are sorted. If the reads have been mapped there will be a series of **@SQ** flags. Additional optional flags are **@RG** which denotes the read groups, **@PG** which denotes the programs used, and **@CO** which is used for additional comments.

![example SAM header](/figures/sam_header_example.png)

To view a SAM/BAM/CRAM file you can use the **samtools view** tool with the **-H** flag:

```bash

samtool view -H sample.bam
```

The alignment field has eleven mandatory fields for each read. 

![table of alignment fields](/figures/sam_alignment_fields.png)

- **QNAME** denotes the query name, if there are muleiple alignment lines in this flag it indicates multimapping or chimeric reads

- **FLAG** a combination of bitwise flags that describe the alignment properties of each segment of the sequence

![table of bitwise FLAGS](/figures/sam_flag-bit-decoder.png)

- **RNAME** the name of the reference sequence aligned to the read in this field

- **POS** the left most position of the first CIGAR operation that "consumes" a reference base

- **MAPQ** mapping quality, 255 means no mapping quality is available

- **CIGAR** represents the type of match between the query and reference

![cigar table](/figures/sam_cigar_key.png)

- **QUAL** Phred scaled base error probability

You can learn more about the SAM file format [here](https://samtools.github.io/hts-specs/SAMv1.pdf).

## FASTA file format
