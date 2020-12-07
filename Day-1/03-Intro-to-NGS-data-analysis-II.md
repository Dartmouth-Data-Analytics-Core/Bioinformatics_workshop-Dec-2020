# Working with NGS data Part II

## Aligned read files, BAM/SAM/CRAM formats

### Principles of read alignment for RNA-seq

The goal of aligning reads to a reference genome is to find the ***most likely location in that reference genome where the read originated from***.

Although we won't go into the theory here, aligning reads to reference genomes involves ***mapping*** to identify the most likely position of the read in the reference genome, followed by the ***alignment***, which describes the base-by-base relationship between the read and the reference. Alignments are often imperfect, and are associated with quality scores (***MAPQ scores***) that describe the quality of the alignment.

**Challenges of aligining millions of short reads to a refence genome involve:**
- Mismatches introduced by **genetic variation** and **sequencing errors**
- **Repeitive sequences** in genomes (e.g. start and end of chromosomes)
- For Eukaryotic genomes the presence of **introns** in reference genomes, meaning aligners must be able to consider **splice-junctions**
- For Prokaryotic genomes the presence of **mobile genetic elements** or **recombination hotspots** in reference genomes

It is important when selecting an aligner to use for your dataset that it is appropriate for your experiment, as numerous aligners exist and make different assumptions and have different strengths/weaknesses. Importantly, some aligners are ***splice-aware*** while others are not. ***Splice-aware*** aligners can generate alignments to a reference genome that span the intronic regions and therefore account for splicing, e.g. `STAR` and `HISAT2`. If your dataset is prokaryotic (non-splicosomal) you would **not** want to use a splice-aware aligner, and instead using an aligner that is not designed to map across intronic regions such as `bwa-mem` or `bowtie2`.

### Concepts for read alignment

**Read clipping**  
Aligners are capable of 'clipping' reads from sequence ends if they do not improve the quality of an alignment that exists for the rest of the sequence.  

There are two type of clipping:  
- *Soft-clipping*: bases at 5' and 3' ends of the read will be kept in the read sequence in the BAM file, but are NOT part of the alignment
- *Hard-clipping*: bases at 5' and 3' ends of the read will be removed from the BAM file altogether and are NOT part of the alignment

Such clipping is commonly used by aligners to get rid of sequence contamination, e.g. ***adapter sequences*** or ***polyA tails*** from mRNAs, so that it does not affect the alignment. This is why you do not necessairily need to be very aggressive in read trimming and pre-processing steps.

Clipping can be very advantageous, but also can potentially cause some issues, read more [here](https://sequencing.qcfail.com/articles/soft-clipping-of-reads-may-add-potentially-unwanted-alignments-to-repetitive-regions/).

**Splicing**  
As discussed above, numerous aligners exist, consisting of both ***splie-aware*** and ***splice-unaware*** aligners. Splice-aware aligners, such as `STAR` and `HISAT2` will produce alignments spanning splice junctions, which is obviously an important characteristic of RNA-seq data that the aligner needs to be able to account for. Furthermore, if you provide coordinates of splice-junctions to aligners like `STAR`, it can improve the mapping over spliced regions and improve detection of novel splice-functions.

**What input do I need for an alignment?**  
At miniumum:  
- `FASTQ` file(s)
- A reference genome (`.fasta`)

Optional:   
- `.gtf` file for the reference genome that species the genomic feature annotation. As mentioned above, if you know where the splice-junctions in your genome are, you can give this to aligners such as STAR and they will use this information to improve the quality of mapping in these regions.

![](../figures/gtf.png)








### Basic concepts of a reference genome

Reference genomes are a concept, not a reality. Reference genomes are a idealized representation of a standard genome from a particular organism.
Reference genomes are created through organizing sequence reads through a process referred to as [*genome assembly*](https://link.springer.com/referenceworkentry/10.1007%2F978-0-387-09766-4_402). 

### Sources of reference genomes

Most reference genomes are hosted on ftp sites where files can copied to a location of your choosing with the `rsync` command. 

```bash

# example rsync command with the UCSC ftp site
# the -a option denotes archive mode, the -P option indicates you want to show the progress as files are downloaded 
rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeDir/ ./

```
 Here are some examples of commonly used sites where reference genomes are hosted.
 
[*RefSeq*](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/)
[*JGI*](https://genome.jgi.doe.gov/portal/help/download.jsf)
[*UCSC*](https://hgdownload.soe.ucsc.edu/downloads.html)
[*Ensembl*]()
[*GENCODE*]()

#### Reference genome annotation

**RefSeq:**
NX_XXXX represents manually curated transcripts.
XM_XXXX represents automatically* annotated transcripts

**Ensembl/GENCODE:**
ENST_XXXX represents manually curated transcripts.



#### Limitations of reference genomes
genetic variation etc.
including [suggestions that we fundamentally change how we use reference genomes](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1774-4).












## SAM/BAM/CRAM file format

**Alignment file formats**  

Read alignments are stored in the ***SAM (.sam)***, ***BAM (.bam))***, and ***CRAM (.cram)*** file formats. ***SAM*** stands for ***Sequence Alignment/Map*** format and is in tab-delimited text format, making it a human readable file (should you dare to look inside, these file are huge). ***BAM*** files are the **compressed, indexed, binary version** of SAM files and are **NOT** human readable, but are much faster to parse and do complex downstream operations on. ***CRAM*** files are compressed versions of the ***BAM*** format, and are not human readable, they are generally only used for storage purposes. You can read all about the SAM/BAM/CRAM file format specification in the documentation [here](https://samtools.github.io/hts-specs/SAMv1.pdf). While you may never need to actually look inside of a SAM/BAM file, it is important to have an understanding of what information is stored in one.

Both formats contain a number of slots for each read alignment that describe key information about the alignment. 11 slots are mandatory, while others are optional and depend on the aligner used, and the settings used in that alignment.

![SAM file](../figures/sam_file.png)
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

![example SAM header](../figures/sam_header_example.png)

### SAM file alignment field

<p align="center">
  <img src="../figures/sam_alignment_fields.png" title="xxxx" alt="context"
	width="80%" height="80%" />
 </p>
 </p>

 The alignment field has eleven mandatory fields for each read, outlined in the table above. Some of the major ones are detailed below.

- **QNAME** denotes the query name, if there are muleiple alignment lines in this flag it indicates multimapping or chimeric reads

- **FLAG** a combination of bitwise flags that describe the alignment properties of each segment of the sequence

<p align="center">
  <img src="../figures/sam_flag-bit-decoder.png" title="xxxx" alt="context"
	width="80%" height="80%" />
 </p>
 </p>

- **RNAME** the name of the reference sequence aligned to the read in this field

- **POS** the left most position of the first CIGAR operation that "consumes" a reference base

- **MAPQ** mapping quality, 255 means no mapping quality is available

- **CIGAR** represents the type of match between the query and reference

<p align="center">
  <img src="../figures/sam_cigar_key.png" title="xxxx" alt="context"
	width="80%" height="80%" />
 </p>
 </p>

- **QUAL** Phred scaled base error probability

You can learn more about the SAM file format [here](https://samtools.github.io/hts-specs/SAMv1.pdf).

# Add coding exercises for building an alignment with a couple of different programs



