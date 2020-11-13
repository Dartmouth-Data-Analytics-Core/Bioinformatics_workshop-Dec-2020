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


## Basic sequence quailty control check 

- FastQC   
- MultiQC  

### Discuss basics of adapter trimming 
### 



----


## Principles of read alignment for RNA-seq
The goal of aliginging reads to a reference genome is to find the ***most likely location in that reference genome where the read originated from***. In the context of RNA-seq, this means we try to find the most likely gene in the reference genome that transcribed the RNA fragment which ultimately ended up in our cDNA library. Doing this for millions of reads allows us to quantify how many RNA fragments were transcribed from a given gene, so we are using the read alignment to measure gene expression. 

Although we won't go into the theory here, alignming reads to reference genomes involves ***mapping*** to identify the most likely position of the read in the reference genome, followed by the ***alignment***, which describes the base-by-base relationship between the read and the reference. Alignments are often imperfect, and are associated with quality scores (***MAPQ scores***) that describe the quality of the alignment. 

**Challenges of aligining millions of short reads to a refence genome involve:**
- Mismatches introduced by **genetic variation** and **sequencing errors**
- **Repeitive sequences** in genomes (e.g. start and end of chromosomes)
- Presence of **intron** in reference genomes, meaning aligners must be able to consider **splice-junctions**

It is important when selecting an aligner to use for your dataset that it is appropriate for your experiment, as numerous aligners exist and make different assumptions and have different strengths/weaknesses. Importantly, some aligners are ***splice-aware*** while others are not. ***Splice-aware*** aligners can generate alignments to a reference genome that span the intronic regions and therefore account for splicing, e.g. `STAR` and `HISAT2`. If your dataset is prokaryotic (non-splicosomal) you would **not** want to use a splice-aware aligner, and instead using an aligner that is not designed to map across intronic regions such as `bwa-mem` or `bowtie2`. 

It is also worth noting here that generating alignments is the most time consuming step of the analytical pipeline for most NGS analyses, and for RNA-seq, the field is moving toward quantification of gene expression using novel *'lightweight'* alignment tools, that are **extremely fast**, e.g. [Kallisto](https://pachterlab.github.io/kallisto/about), [Salfish](https://www.nature.com/articles/nbt.2862), and [Salmon](https://combine-lab.github.io/salmon/). These algorithms avoid generation of typical base-by-base alignments, and instead generate ***psuedo-alignments***, which have been shown to produce very accurate estimates of transcript abundances in RNA-seq data. 

![Read alignment](../figures/read_alignment.png)
  

### Concepts for read alignment 

**Read clipping**  
Aligners are capable of 'clipping' reads from sequence ends if they do not improve the quality of an alignment that exists for the rest of the sequence.  

There are two type of clipping:  
- *Soft-clipping*: bases at 5' and 3' ends of the read will be kept in the read sequence in the BAM file, but are NOT part of the alignment
- *Hard-clipping*: bases at 5' and 3' ends of the read will be removed from the BAM file altogether and are NOT part of the alignment 

Such clipping is commonly used by aligners to get rid of sequence contamination, e.g. ***adapter sequences*** or ***polyA tails*** from mRNAs, so that it does not affect the alignment. At least for RNA-seq, this is why you do not necessairily need to be very aggressive in read trimming and pre-processing steps. 

Clipping can be very advantageous, but also can potentially cause some issues, read more [here](https://sequencing.qcfail.com/articles/soft-clipping-of-reads-may-add-potentially-unwanted-alignments-to-repetitive-regions/). 

**Splicing**  
As discussed above, numerous aligners exist, consisting of both ***splie-aware*** and ***splice-unaware*** aligners. Splice-aware aligners, such as `STAR` and `HISAT2` will produce alignments spanning splice junctions, which is obviously an important characteristic of RNA-seq data that the aligner needs to be able to account for. Furthermore, if you provide coordinates of splice-junctions to aligners like `STAR`, it can improve the mapping over spliced regions and improve detection of novel splice-functions. 

**Genome vs transcriptome mapping?**  
While there are times when one may want to map to a transcriptome, there are issues with this approach.  
- If your annotated transcriptome is not complete, you may fail to map some reads simply because the sequences aren't in your reference, which would not be an issue if you mapped to the genome. 
- With multiple splice isoforms it is difficult to disambiguate which splice isoform a read should be aligned to in transcriptome mapping. 
- You cannot identify novel transcripts this way.

**What input do I need for an alignment?**  
At miniumum:  
- `FASTQ` file(s)
- A reference genome (`.fasta`)

Optional:   
- `.gtf` file for the reference genome that species the genomic feature annotation. As mentioned above, if you know where the splice-junctions in your genome are, you can give this to aligners such as STAR and they will use this information to improve the quality of mapping in these regions. 

![](../figures/gtf.png)


## SAM/BAM/CRAM file format

**Alignment file formats**  

Read alignments are stored in the ***SAM (.sam)*** and )***BAM (.bam))*** file format. )***SAM)*** stands for )***Sequence Alignment/Map)*** format and is in tab-delimited text format, making it a human readable file (should you dare to look inside, these file are huge). )***BAM)*** files are the **compressed, indexed, binary version** of SAM files and are **NOT** human readable, but are much faster to parse and do complex downstream operations on. You can read all about the SAM/BAM file format specification in the documentation [here](https://samtools.github.io/hts-specs/SAMv1.pdf). While you may never need to actually look inside of a SAM/BAM file, it is important to have an understanding of what information is stored in one. 

Both formats contain a number of slots for each read alignment that describe key information about the alignment. 11 slots are mandatory, while others are optional and depend on the aligner used, and the settings used in that alignment.

![SAM file](../figures/sam-file.png)
The image for the example BAM file is take from the [SAM/BAM file format documentation](https://samtools.github.io/hts-specs/SAMv1.pdf)

#### Notes on select fields:

**FLAG**:  
Encodes important information about the read, for example, is it a ***primary***, ***secondary***, or ***supplementary*** alignment. Since a single read will likely have a number of properties that we want to ***'flag'***, SAM files use a special way of encoding the FLAG field to pack as much information as possible into a single number. While we won't go into detail on this here, SAM/BAM file use a *bit-wise* system to combine information across flags into a **single integer**. I encourage you to go read more about FLAGs and how they are specified in the SAM/BAM documentation. 

The Broad institute provides an [excellent tool](https://broadinstitute.github.io/picard/explain-flags.html) for decomposing SAM flags into the proprties of the read that make up a specific `FLAG` value. 

This command will provide basic information on FLAGs from samtools.
```bash 
samtools flags
```
The values shown here relate the the [hexadecimal system](https://www.electronics-tutorials.ws/binary/bin_3.html)

**MAPQ**:   
Corresponds to the quality of the mapping. These are calculated in the same way as the Phred scores `Q = -10 x log10(P)`, although are generally considered to be a best guess form the aligner. A MAPQ of 255 is used where mapping quality is not available. Some aligners also use specific values to represent certain types of alignments, which may affect use of downstream tools, so it is worth understanding those that are specific to your aligner. 

**CIGAR**  
An alphanumerical string that tells you information about the alignment. For relatively short reads, these are nice, but for long reads, they are a headache. Numbers correspond to number of bases, and letters correspond to features of those bases.  

Letter key for CIGAR strings: 
M = match or mismatch  
S = soft clip  
H = hard clip  
I = insertion  
D = deletion  
N = skipping  

So for example, alignment in row 3 of our SAM file example above (`5S6M`) would describe an alignment where 5 bases are soft-clipped, followed by 6 matching bases. 


SAM format is a common way of representing sequenced reads, especially after reads have been aligned or mapped to a reference genome. BAM is a binary (non-human readable) format of SAM that takes up less space. CRAM files are compressed versions of BAM files - these take up the least space and are recommended for longer term storage of alignment files. SAM/BAM/CRAM files can be converted back and forth with the tool **samtools view**.

```bash

#convert sam to bam
samtools view -b sample.sam > sample.bam

#convert bam to cram
samtools view -C sample.bam > sample.cram
```

A SAM file is made up of two basic parts, the header and the alignment. 

### SAM file header field

All header lines will start with the **@** symbol. The mandatory flag **@HD** will come first in the file and should only occur once, this flag has the meta-data that pertains to the SAM file and will either have a **GO** field indicating that reads are grouped but not sorted or a **SO** field indicating that reads are sorted. If the reads have been mapped there will be a series of **@SQ** flags. Additional optional flags are **@RG** which denotes the read groups, **@PG** which denotes the programs used, and **@CO** which is used for additional comments.

To view a SAM/BAM/CRAM file you can use the **samtools view** tool with the **-H** flag:

```bash

samtool view -H sample.bam
```

![example SAM header](/figures/sam_header_example.png)

### SAM file alignment field

<p align="center">
  <img src="/figures/sam_alignment_fields.png" title="xxxx" alt="context"
	width="80%" height="80%" />
 </p>
 </p>
 
 The alignment field has eleven mandatory fields for each read, outlined in the table above. Some of the major ones are detailed below. 
 
- **QNAME** denotes the query name, if there are muleiple alignment lines in this flag it indicates multimapping or chimeric reads

- **FLAG** a combination of bitwise flags that describe the alignment properties of each segment of the sequence

<p align="center">
  <img src="/figures/sam_flag-bit-decoder.png" title="xxxx" alt="context"
	width="80%" height="80%" />
 </p>
 </p>

- **RNAME** the name of the reference sequence aligned to the read in this field

- **POS** the left most position of the first CIGAR operation that "consumes" a reference base

- **MAPQ** mapping quality, 255 means no mapping quality is available

- **CIGAR** represents the type of match between the query and reference

<p align="center">
  <img src="/figures/sam_cigar_key.png" title="xxxx" alt="context"
	width="80%" height="80%" />
 </p>
 </p>
 
- **QUAL** Phred scaled base error probability

You can learn more about the SAM file format [here](https://samtools.github.io/hts-specs/SAMv1.pdf).


## QUantification

<p align="center">
<img src="../figures/genoic-content.png" title="xxxx" alt="context"
	width="100%" height="100%" />
</p>
</p>

***Fig. 1.** Genomic context of aligned reads can fall into several categories dictated by the annotation used. *



---- 


cover zwero and 1 -based coord systems 



## TO DO: 
- Add basic concepts of alignment, assembly, quantification, peak calling as file types are introduced  



