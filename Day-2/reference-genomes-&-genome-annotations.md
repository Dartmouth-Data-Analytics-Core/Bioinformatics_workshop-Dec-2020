# Reference Genomes and Genome annotation



### Basic concepts of a reference genome 

Reference genomes are a concept, not a reality. Reference genomes are a idealized representation of a standard genome from a particular organism. 

Reference genomes are created through organizing sequence reads through a process referred to as *genome assembly*.  


N50
Contig number
Gaps
Patches 
Alternate-loci/-contigs

Coordinates are the same within *patches*, but are completely different between *builds*. e.g. the nucleotide at position 12,332,681 on chr12 of *hg19* may be in a very different location of chromosome 12 in *hg38*. If you ever need to convert coordiates between genome builds, there are various *'liftOver'* tools that exist that are capable of this.  


### Sources of reference genomes

[*RefSeq*](https://www.ncbi.nlm.nih.gov/refseq/)
*UCSC* 
*Ensembl* 
*GENCODE* 



[Genome Reference Consortium (GRC)](https://www.ncbi.nlm.nih.gov/grc). 

Lets go and explore a refence genome from the **GRC** on their [website](https://www.ncbi.nlm.nih.gov/grc). 

[Genome in a Bottle Consortium (GIAB)](https://www.nist.gov/programs-projects/genome-bottle). 



#### Which reference genome should I use? 



#### Limitations of reference genomes

including [suggestions that we fundamentally change how we use reference genomes](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1774-4).


### Genome annotation 


discuss annotation pipelines


#### Transcript annotation 

**RefSeq:** 
NX_XXXX represents manually curated transcripts. 
XM_XXXX represents automattically* annotated transcripts

Used a lot in clinical annptation projects (e.g. 

**Ensembl/GENCODE:** 
ENST_XXXX represents manually curated transcripts. 

Used by several major projects such as GTEx, ICGC, 1000Genomes 

