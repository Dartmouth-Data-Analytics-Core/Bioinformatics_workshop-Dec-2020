# Working with NGS data Part II

## Alignment files (BAM/SAM/CRAM formats)

### Principles of read mapping for RNA-seq

For NGS experiments, the goal of read mapping is to find an alignment that describes the **most likely location in the reference genome where that read originated from**. This is generally determined by reducing the information in the reference and query (read to be mapped) into smaller strings and looking for the position in the reference with the highest number of matching smaller strings. This process is also used in *de novo* genome assembly, local alignments (BLAST or BLAT), and global alignments.

Although we won't go into the theory here, aligning reads to reference genomes involves **mapping** to identify the most likely position of the read in the reference genome, followed by the **alignment**, which describes the base-by-base relationship between the read and the reference.

<p align="center">
<img src="../figures/read_alignment.png" title="xxxx" alt="context"
	width="100%" height="100%" />
</p>

Challenges of aligning millions of short reads to a reference genome involve:
- Mismatches introduced by genetic variation and sequencing errors
- Repetitive sequences in genomes (e.g. start and end of chromosomes)
- For Eukaryotic genomes the presence of introns in reference genomes, meaning aligners must be able to consider splice-junctions
- For Prokaryotic genomes the presence of mobile genetic elements or recombination hotspots in reference genomes

It is important when selecting an aligner to select one appropriate for your experiment. Different aligners exist and generally have different properties and applications. For example, some aligners are **splice-aware** while others are not. Splice-aware aligners can generate alignments that span intronic regions and therefore account for splicing, e.g. `STAR` and `HISAT2`. If your dataset is prokaryotic (non-splicosomal) you would not want to use a splice-aware aligner, and instead using an aligner that is not designed to map across intronic regions such as `bwa-mem` or `bowtie2`.













### What is a reference genome?

Reference genomes are a concept, not a reality. Reference genomes are a idealized representation of a standard genome from a particular organism.
Reference genomes are created through organizing sequence reads through a process referred to as [*genome assembly*](https://link.springer.com/referenceworkentry/10.1007%2F978-0-387-09766-4_402).

### Sources of reference genomes

Most reference genomes are hosted on ftp sites where files can copied to a location of your choosing with the `rsync` command.

*Do not run this command - this is only an example*
```bash

# example rsync command with the UCSC ftp site
# the -a option denotes archive mode, the -P option indicates you want to show the progress as files are downloaded
rsync -a -P rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeDir/ ./

```
 Here are some examples of commonly used sites where reference genomes are hosted.

- [*RefSeq*](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/)
- [*JGI*](https://genome.jgi.doe.gov/portal/help/download.jsf)
- [*UCSC*](https://hgdownload.soe.ucsc.edu/downloads.html)
- [*Ensembl*](http://ftp.ensembl.org/pub/)
- [*GENCODE*](http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/)

# Reference genome annotation

**RefSeq:**
NX_XXXX represents manually curated transcripts.
XM_XXXX represents automatically* annotated transcripts

**Ensembl/GENCODE:**
ENST_XXXX represents manually curated transcripts.

#### Limitations of reference genomes
By definition a constant property of genomes is genetic variation between individuals, this feature holds true for both single celled organisms where variation is introduced through mutation and recombination and multicellular sexual organisms where mutation and sexual recombination constantly stir the gene pool to introduce variation. It is therefore impossible to define a reference genome for a given species that encompasses all of the variation within that species.

However, Current reference genomes work as the foundation for all genomic data and databases. They provide a scaffold for genome assembly, variant calling, RNA or other sequencing read alignment, gene annotation, and functional analysis. Genes are referred to by their loci, with their base positions defined by reference genome coordinates. Variants and alleles are labeled as such when compared to the reference (i.e., reference (REF) versus alternative (ALT)). Related (strains, diploid, or personal) genomes are assembled using the reference as a scaffold, and RNA-seq reads are typically mapped to the reference genome.

At this point we need a baseline for all of the genomics tools that we use, however many scientists are considering [solutions that could fundamentally change how we think of and build reference genomes](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1774-4).















### General concepts for read alignment

**Read clipping:**  
Aligners are capable of 'clipping' reads from sequence ends if they do not improve the quality of an alignment that exists for the rest of the sequence.  

There are two type of clipping:  
- *Soft-clipping*: bases at 5' and 3' ends of the read will be kept in the read sequence in the BAM file, but are NOT part of the alignment
- *Hard-clipping*: bases at 5' and 3' ends of the read will be removed from the BAM file altogether and are NOT part of the alignment

Such clipping is commonly used by aligners to get rid of sequence contamination, e.g. adapter sequences or polyA tails from mRNAs, so that it does not affect the alignment. This is why you do not necessarily need to be very aggressive in read trimming and pre-processing steps. Clipping can be very advantageous, but also can potentially cause some issues, read more [here](https://sequencing.qcfail.com/articles/soft-clipping-of-reads-may-add-potentially-unwanted-alignments-to-repetitive-regions/).

**Splicing:**  
As discussed above, numerous aligners exist, consisting of both ***splice-aware*** and ***splice-unaware*** aligners. Splice-aware aligners, such as `STAR` and `HISAT2` will produce alignments spanning splice junctions, which is obviously an important characteristic of RNA-seq data that the aligner needs to be able to account for. Furthermore, if you provide coordinates of splice-junctions to aligners like `STAR`, it can improve the mapping over spliced regions and improve detection of novel splice-functions.


**What input do I need for an alignment?**  
At miniumum:  
- `FASTQ` file(s)
- A reference genome (`.fasta`)

Optional:   
- `.gtf` file for the reference genome that species the genomic feature annotation. As mentioned above, if you know where the splice-junctions in your genome are, you can give this to aligners such as STAR and they will use this information to improve the quality of mapping in these regions.

![](../figures/gtf.png)

**Alignment file formats**  

Read alignments for NGS data are stored in three major file formats: *SAM (.sam)*, *BAM (.bam)*, and *CRAM (.cram)*.

- **SAM (Sequence Alignment/Map)** format - tab-delimited text format, so is human readable file (should you dare to look inside)
- **BAM** files are compressed, binary version of SAM files and are NOT human readable, but are much faster to parse and do complex operations with.
- **CRAM** files are compressed versions of the BAM format, and are not human readable, they are generally only used for storage purposes.

You can read all about the SAM/BAM/CRAM file format specification in the documentation [here](https://samtools.github.io/hts-specs/SAMv1.pdf). While you may never need to actually look inside of a SAM/BAM file, it is important to have an understanding of what information they contain.

Alignment files are composed of two basic sections:
- the header
- the alignments

All header lines start with the `@` symbol. The mandatory flag `@HD` will come first in the header and can be followed by a number of additional flags that represent features of the alignments in the file (e.g. `SO`, indicating reads are sorted by coordinate).

The alignment section contains a number of 'slots' for each read alignment that describe key information about the alignment. 11 slots are mandatory, while others are optional and depend on the aligner used, and the settings used for mapping.

<p align="center">
<img src="../figures/sam-file.png" title="xxxx" alt="context"
	width="95%" height="95%" />
</p>

SAM/BAM/CRAM files can be viewed, queried, and maniuplated using the [Samtools software suite](http://www.htslib.org/), which we will explore the usage of in more detail later in this lesson.


#### Notes on select SAM fields:

**FLAG:**  
Encodes important information about the read, for example, is it a *primary*, *secondary*, or *supplementary* alignment. Since a single read will likely have a number of properties that we want to *'flag'*, SAM files use a special way of encoding the FLAG field to pack as much information as possible into a single number. While we won't go into detail on this here, SAM/BAM file use a bit-wise system to combine information across flags into a single integer.

I encourage you to go read more about FLAGs and how they are specified in the SAM/BAM documentation.The Broad institute also provides an [excellent tool](https://broadinstitute.github.io/picard/explain-flags.html) for decomposing SAM flags into the properties of the read that make up a specific `FLAG` value.

This command will provide basic information on FLAGs from samtools.
```bash
samtools flags
```
The values shown here relate the the [hexadecimal system](https://www.electronics-tutorials.ws/binary/bin_3.html)

**MAPQ:**  
Corresponds to the quality of the mapping. These are calculated in the same way as the Phred scores `Q = -10 x log10(P)`, although are generally considered to be a best guess form the aligner. A MAPQ of 255 is used where mapping quality is not available. Some aligners also use specific values to represent certain types of alignments, which may affect use of downstream tools, so it is worth understanding those that are specific to your aligner.

**CIGAR:**  
An alphanumerical string that tells you information about the alignment. For relatively short reads, these are nice, but for long reads, they are a headache. Numbers correspond to number of bases, and letters correspond to features of those bases.  

Letter key for CIGAR strings:
`M` = match or mismatch  
`S` = soft clip  
`H` = hard clip  
`I` = insertion  
`D` = deletion  
`N` = skipping  

So for example, alignment in row 3 of our SAM file example above (`5S6M`) would describe an alignment where 5 bases are soft-clipped, followed by 6 matching bases.














### Perform an alignment

For this dataset we are using STAR because the reads were generated from a eukaryotic RNA-seq experiment which require a splice aware aligner to generate appropriate alignments.




index




```bash
# create a directory to store aligned files
mkdir ../aligned/

# enter that directory
cd ../aligned

#run splice aware alignment
STAR --genomeDir /scratch/fund_of_bioinfo/ref/hg38_chr20_index \
  --readFilesIn ../trim/SRR1039508_1.trim.chr20.fastq.gz ../trim/SRR1039508_2.trim.chr20.fastq.gz \
  --readFilesCommand zcat \
  --sjdbGTFfile /scratch/fund_of_bioinfo/ref/Homo_sapiens.GRCh38.97.chr20.gtf \
  --runThreadN 1 \
  --outSAMtype BAM SortedByCoordinate \
  --outFilterType BySJout \
  --outFileNamePrefix SRR1039508.
```

Option details:
- `--genomeDir`: the path to the directory with genome indices
- `--readFilesIn`: read files to map to reference alignment
- `--readFilesCommand`: uncompression command to apply to read files
- `--sjdbGTFfile`: the path to the annotation file that includes cooordinates of splice-junctions
- `--runThreadN`: number of threads to use in the run
- `--outSAMtype`: (SAM/BAM unsorted/ BAM SortedByCoordinate)
- `--outFilterType`: how mapped reads will be filtered (normal/BySJout)
- `--outFileNamePrefix`: prefix for outfiles generated in the run






The predominant metric of interest to detrmine if our aignment was successful is the % of uniquely mappng reads, what value means the experiment was successful depends on the data type




alignments can be viewed in igv, we will do this on day-2.


<p align="center">
<img src="../figures/igv-03.png" title="xxxx" alt="context"
	width="95%" height="95%" />
</p>












```bash
#convert bam to sam
samtools view -H sample.bam > sample.sam

#convert sam to bam
samtools view -b sample.sam > sample.bam

#convert bam to cram
samtools view -C sample.bam > sample.cram
```

header

```bash
samtools view -H sample.bam
```






### Run alignments for multiple samples
```bash
ls ../trim/*_1.trim.chr20.fastq.gz | while read x; do

  # save the file name
  sample=`echo "$x"`
  # get everything in file name before "/" (to remove '../trim/')
  sample=`echo "$sample" | cut -d"/" -f3`
  # get everything in file name before "_" e.g. "SRR1039508"
  sample=`echo "$sample" | cut -d"_" -f1`
  echo processing "$sample"

  # run STAR for each sample
  STAR --genomeDir /dartfs-hpc/scratch/rnaseq1/refs/hg38_chr20_index \
    --readFilesIn ../trim/${sample}_1.trim.chr20.fastq.gz ../trim/${sample}_2.trim.chr20.fastq.gz \
    --readFilesCommand zcat \
    --sjdbGTFfile /dartfs-hpc/scratch/rnaseq1/refs/Homo_sapiens.GRCh38.97.chr20.gtf \
    --runThreadN 4 \
    --outSAMtype BAM SortedByCoordinate \
    --outFilterType BySJout \
    --outFileNamePrefix ${sample}.
done

# index the BAMs
samtools index *sortedByCoord.out.bam
```
Note that I change `--outSAMtype` to `BAM sortedByCoord` so that we dont have to convert SAM to BAM and run `sort`.

View the reports quickly:
```bash
ls *Log.final.out | while read x; do
   yes '' | sed 4q
   echo Printing $x
   yes '' | sed 1q
   cat $x
done
```






> Important Note: After generating alignments, it is critical to inspect their quality.





## Break out room exercises

- Look at the SAM file that you created

- Convert the SAM file to a BAM file

- Convert the SAM file to a CRAM file

- How much space does each file take up?

- What is the best way to store the aligned file to minimize the space constraints?
