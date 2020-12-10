
### Genome reference sequences in Bioconductor

Beyond providing access to extensive annotation data in R, Bioconductor also provides functionality to obtain and maniuplate the complete reference sequences for commonly used genomes. Specifically, the [BSgenome](https://bioconductor.org/packages/release/bioc/html/BSgenome.html) family of Bioconductor packages provides an efficient way to obtain, query, and maniuplate genomic sequence data from reference genomes. You can return a vector of the currently available genomes to your console by printing `available.genomes()` after loading the `BSgenome` package.

```r
library(BSgenome)

available.genomes()
```

The available genomes are preominantly based around NCBI and UCSC genomes, however functionality does exist for [forging a *BSGenome*](https://bioconductor.org/packages/release/bioc/html/BSgenome.html) package, allowing you to leverage the *BioStrings* framework for genomes that are not part of the available set from `BSgenome`.

If your genome is not included in the available genomes but you would still like to leverage the `BioStrings` and `BSGenome` framework, you can [forge a BSGenome package](https://bioconductor.org/packages/release/bioc/html/BSgenome.html) following instructions available at the BioConductor website.

BSgenome packages are all heavily dependent on another Bioconductor package, [BioStrings](http://bioconductor.org/packages/release/bioc/html/Biostrings.html) , which defines a set of methods and object classes for storing and accessing the sequence data, and is loaded automatically when you loaded `BSgenome`. We will introduce the basic object classes and methods introduced by `BioStrings` to demonstrate how they form the basis for `BSgenome` packages, and what genomic sequence-based operations they allow you to perform on reference genome data (or any sequence you define).

Analyzing genomic sequence data directly can be used for a number of common research tasks, all possible in the *BSGenome/BioStrings* framework, for example:  
* Extracting DNA/RNA/protein sequences of specific genes or gene regions of interest (e.g. DNA sequence flanking a ChIP-seq peak)
* Calculating nucleotide frequencies among specific genomic features
* Searching for matching sequences of interest (e.g. barcode matching)

#### Basic object-classes and methods in the *BioStrings* package

Before working with a complete reference genome sequence from *BSGenome*, lets discuss the basic object-classes implemented in the *BioStrings* package to store sequence data, as well as the methods use to parse these sequences.

The most basic object class in *BioStrings* is the *XString* class, which is technically a '*virtual class*' (meaning it cannot actually store objects itself, but can be used to set rules for a a group of classes) encompassing object classes for DNA, RNA, and protein sequences: `DNAString`, `RNAString`, and `AAString`. Today we will focus on DNA sequences using *DNAString* class objects.

Lets start by creating a **really** simple *DNAString* object and looking at some basic features of it:
```r
# use the DNAString contructor function to create a 10 letter DNA sequence
seq <- DNAString(x="AGCT", start=1, nchar=NA)
seq

# how long is it
length(seq)

# show the structure of the DNAString object
str(dna.st)
```

This sequence is a bit short and unrealistic. Lets make a longer sequence using the pre-stored *BioStrings* object `DNA_ALPHABET` and some functions from base R.
```r
# print DNA alphabet to see what it returns
DNA_ALPHABET
```

You can see the 4 standard DNA bases are returned as the first 4 elements of this character string. The remaining elements represent ambiguous bases or specific combinations/relationships using something called the *extended The International Union of Pure and Applied Chemistry (IUPAC) genetic alphabet.* **BioStrings** and its object classes use the extended IUPAC genetic alphabet to describe nucelic acid sequences, therefore we will nbriefly cover the basics of the extended IUPAC alphabet now.

The extended IUPAC code is a specific nomenclature designed to describe incompletely specified nucleic acids. The standard IUPAC code uses 16 characters to specify single bases (A, G, C, T, U) in nucleic acid sequences, or various possible states that a specific nucleic acid may exist as in a sequence.

The standard IUPAC code is used by numerous bioinformatics tools and softwares in order to represent complex sequences of nucleic acids for which we may not be confident in some individual base identities (e.g. complex genomic regions that are challenging to sequence using short read approaches).

**Table 1. Standard IUPAC genetic alphabet.**

|Symbol |	Mnemonic| Translation  |
|---|---|---|
| A	|	A | (adenine) |  
| C	|	C | (cytosine)  |
| G	|	G	| (guanine)  |
| T	|		T	| (thymine)  |
| U	|		U	| (uracil)  |
| R	|	pu**R**ine		| A or G  |
| Y		| p**Y**rimidine		| C or T/U  |
| S		| **S**trong interaction	|	C or G  |
| W		| **W**eak interaction		| A or T/U  |
| M		| a**M**ino group		| A or C  |
| K		| **K**eto group		| G or T/U |   
| H		| not G		| A, C or T/U |   
| B		| not A		| C, G or T/U |  
| V		| not T/U		| A, C or G |  
| D		| not C		| A, G or T/U |  
| N		| a**N**y		| A, C, G or T/U |  
| - | none | Gap  

This table was adapted from [Johnson, 2010, *Bioinformatics*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2865858/#B1).

An extended IUPAC genetic alphabet was also described in 2010 [(Johnson, 2010, *Bioinformatics*)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2865858/#B1). The extended code uses additional characters, underlining, and bolding as well as the original 16 character code (all meanings maintained) to denote all possible combinations or relationships between bases. Among other uses, this has been valuable for representing genetic varaition in DNA sequences. You can explore the details on the extended code in **Tables 2 & 3** of [(Johnson, 2010)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2865858/#B1).

 If you're working with *BioStrings* objects, and need a reminder of the basic characters of the extended code, you can just type `IUPAC_CODE_MAP` when you have the *BioStrings* package loaded into your R session.

For example, entering the command below into your command-line..
```
IUPAC_CODE_MAP
```
...will return the following to your console:
```r
     A      C      G      T      M      R      W      S      Y      K      V      H      D      B  N
   "A"    "C"    "G"    "T"   "AC"   "AG"   "AT"   "CG"   "CT"   "GT"  "ACG"  "ACT"  "AGT"  "CGT" "ACGT"
```

For now, lets use the standard, unambiguous DNA bases (A, G, C, T).
```{r}
# randomly sample from specific characters in DNA_ALPHABET to create a longer sequence
seq = sample(DNA_ALPHABET[c(1:4, 16)], size=100, replace=TRUE)
seq

# use the paste command to collapse the characters into a single string
seq = paste(seq, collapse="")
seq

# use the DNAString constructor function to turn this sequence into a DNAString class object
seq.dnastring <- DNAString(seq)
seq.dnastring
```

Now collect some basic information on your sequence.
```r
# confirm how long it is
length(seq.dnastring)

# what is the frequency of each base in your sequence
alphabetFrequency(seq.dnastring, baseOnly=TRUE, as.prob=TRUE)

# what is the frequency of your favourite base (which is obviously Adenine)
letterFrequency(seq.dnastring, "A", as.prob=TRUE)

# return the frequency of dinucleotide pairs
dinucleotideFrequency(seq.dnastring, as.prob=TRUE)

# or even trinucleotides
trinucleotideFrequency(seq.dnastring, as.prob=TRUE)
```

We can also perform some basic manipulations of our sequence using *BioStrings* functions.
```r
# subset the sequence for 10 specific bases of interest  
seq.dnastring[10:19]

# this can also be done with the `subseq()` function  
subseq(seq.dnastring, 10, 19)

# get the reverse of our sequence
reverse(seq.dnastring)

# get the reverse COMPLEMENT sequence
reverseComplement(seq.dnastring)

# get the reverse complement of the first 10 bases in your sequence
reverseComplement(subseq(seq.dnastring, 1, 10))

# translate our DNA sequence
translate(seq.dnastring)
```

Again, our example is a little impractical since we are usually working with a set of sequences, for example the chromosomes in a reference genome. This is where the `DNAStringSet` object class becomes useful. `DNAStringSet` allows you to store, name, and manipulate multiple sequences in one *BioStrings* object.
```r
# remove the old single sequence from our global R environment
rm(seq)

# create a new variable, and fill it with individual sequences created as we did above
seq.dnass <- NULL
seq.dnass[1] = paste(sample(DNA_ALPHABET[c(1:4)], size=50, replace=TRUE), collapse="")
seq.dnass[2] = paste(sample(DNA_ALPHABET[c(1:4)], size=50, replace=TRUE), collapse="")
seq.dnass[3] = paste(sample(DNA_ALPHABET[c(1:4)], size=50, replace=TRUE), collapse="")
seq.dnass[4] = paste(sample(DNA_ALPHABET[c(1:4)], size=50, replace=TRUE), collapse="")
seq.dnass[5] = paste(sample(DNA_ALPHABET[c(1:4)], size=50, replace=TRUE), collapse="")
seq.dnass

# how long is this object
length(seq.dnass)

# use the constructor function DNAStringSet to make the DNAStringSet object with your sequences
dna.st.set = DNAStringSet(seq.dnass)

# how long is the DNAStringSet object
length(seq.dnass)

# name all your sequences
names(seq.dnass) = paste("barcode-", 1:5, sep="")
seq.dnass
```

Like *XString*, a virtual class exists for `DNAStringSet` class objects called *XStringSet*, which also contains object classes for storing RNA and AA sequences (`RNAStringSet` and `AAStringSet`).

#### Working with *BSGenome* reference genomes

Now that we understand the major classes implemented by *BioStrings*, lets load a complete reference genome and start exploring it.
```r
# assign the genome to a variable using getBSgenome() (you need to have the package for the BSgenome you are trying to load already installed)
genome <- getBSgenome("BSgenome.Mmusculus.UCSC.mm10")
genome

# check the structure
str(genome)

# print the metadata for the genome
metadata(genome)
```

By default, the *BSGenomes* come with no sequence masking of, for example known repetitive regions that you may wish to ignore in your analyses. To obtain a masked genome, you should set `masked=TRUE` in the `getBSgenome()` function. This will load a genome in which specific sequences have been masked in a hierachical fashion using the following criteria:  
1. Gaps in the genome assembly
2. Sequences with intra-contig ambiguities
3. regions flagged by [*RepeatMasker*](http://www.repeatmasker.org/)
4. regions flagged by [*Tandem Repeat Finder*](https://tandem.bu.edu/trf/trf.html)

Lets load in the masked reference and compare to the unmasked version.
```r
genome.m <- getBSgenome("BSgenome.Mmusculus.UCSC.mm10.masked")
class(genome.m)
genome.m

# unmaksed genome
class(genome)
genome
```

The masked genomes utilized the `MaskedXString` class implemented by *BioStrings* to denote the masked sequences. If you print a specific sequence to the console from the masked genome, you will also get a high-level summary of masking for that sequence.
```r
# return basic sequence information summary
seqinfo(genome.m)

# print chromosome 1
genome.m$chr1

# unmasked genome
seqinfo(genome)
genome$chr1
```

Lets move forward with the masked genome for today. Remove the `genome` variable from your working environment and replace it with the masked genome for convienience.
```r
rm(genome)

genome <- getBSgenome("BSgenome.Mmusculus.UCSC.mm10.masked")
```

We can perform all the same *BioStrings* based-methods on the sequences stored in our *BSGenome* object. For example:
```r
# assign chr 1
chr1 <- genome$chr1

# confirm how long it is
length(chr1)

# subset it
chr1[1:100]
chr1[100498:100598]

# what is the frequency of each base in your sequence
alphabetFrequency(chr1, baseOnly=TRUE, as.prob=TRUE)

# what is the frequency of your favourite base
letterFrequency(chr1, "A", as.prob=TRUE)
```

Beyond the basic *BioStrings* based methods, there is one very important method implemented by the *BSGenome* using the `getSeq()` function, that allows you to extract specific sequences at request from a *BSgenome* or *XStringSet* class object. We will use `getSeq()` functionality in our next example to illustrate how you might use BSGenome packages in a standard NGS analysis.

#### Example usage of a BSGenome package: Extracting peak flanking sequences from ChIP-seq data

Once peak regions have been identified to describe the potential binding locations of a particular transcription factor (TF), a common task in the analysis of ChIP-seq data is to scan the sequences immediately surrounding these peaks in order to identify sequences enriched over these peak regions that may represent the binding motif for that TF. In order to achieve this, we need to obtain the sequences for these peaks from the reference genome that the samples were aligned to (mm10). The cartoon below depicts this overall workflow.

<img src="../figures/motif-example.png" height="550" width="900"/>

As an example, we will continue our mouse forebrain theme, using ChIP-seq data from the developing mouse forebrain that was performed using an antibody specific for the CTCF transcription factor (TF), a critical TF for diverse cellular processes that performs a plethora of transcriptional activation/repression functions at a genome-wide level. Called CTCF peaks for this experiment were downloaded from the ENCODE website [here](https://www.encodeproject.org/experiments/ENCSR677HXC/).

Lets read in the BED file as a *GRanges* object using *rtracklayer* function `import()` as we have done previously. We can then use the `getSeq()` function in Bioconductor to return the sequences from our previously assigned *BSGenome* (UCSC - mm10, assigned to `genome`) that cover the regions specified in the *GRanges* object.
```r
# read in peaks
bed <- import("data/CTCF-mouse-forebrain-mm10.bed", format="BED")
bed

# extract sequences for peak regions and print to console
ctcf_seqs <- getSeq(genome, bed)
ctcf_seqs
```

Since the object returned by `getSeq()` is a DNAStringSet class object, we can use *BioStrings* based methods to perform operations on the sequences directly. For example, we might be interested in checking the nucleotide frequencies across all peaks.
```r
# calculate nucleotide freqs.
nt_freqs <- alphabetFrequency(ctcf_seqs, baseOnly=TRUE, as.prob=TRUE)

# calculate mean value for nucleotide freqs across all peaks
round(apply(nt_freqs, 2, mean), digits=2)
```

We might also be interested in visualizaing the distribution of the peak width, to get an idea of how much they vary. We can use the `width` accessor function to extract the width of each peak, and base R functions for plotting.
```r
hist(width(ctcf_seqs),
     col = "darkgray",
     xlab = "Peak width (bp)",
     main = "CTCF peak width distribution")
```

We could now export these sequences to a FASTA file (using `writeXStringSet()`) however several motif discovery softwares require that peaks be of the same size (width). To do this in a meaningful way for our ChIP-seq data, we will need to find the center of each peak, and then restrict to the a certain number of bases flanking either side of the center position. We will need to go back to the ranges from our original BED file to resize the peaks to the desired width around the center, then re-extract the sequences fopr those regions.
```r
# resize the regions from the BED file
bed_centered <- resize(bed, width = 400, fix = "center")
bed_centered

# check their with
width(bed_centered)

# extract sequences again
ctcf_seqs_cent <- getSeq(genome, bed_centered)
ctcf_seqs_cent
```

Now we are ready to export these sequences in FASTA file format, which is used as the default format as input to many motif discovery algorithms. As mentioned above, we can do this for DNAStringSet objects with the function `writeXStringSet()`.
```r
# export peaks to FASTA file
writeXStringSet(ctcf_seqs, file=paste0(peak_dir, "CTCF-peaks-resized.fa"))
```

After you write the file, go to your the bash command line and have a look at your FASTA file to confirm it looks correct.

**Note:** As we have discussed before, there are several other ways you could have performed this task outside of the functionality implemented by *BioStrings* and *BSGenome*. The major advantage of performing this analysis in R/Bioconductor is that you do not need to host any large reference genome files locally except for installing the *BSGenome* packages, functionality from other Bioconductor packages can be utilized (such as how we used the *GRanges* `resize()` function above), and you can easily leverage the other functionality for sequence operations available in *BioStrings* (of which there are many).

If you did have direct access to the reference genome locally and other functionality in Bioconductor wasn't a priority for you, you could perform this analysis at the command line in bash with [*bedtools*](https://bedtools.readthedocs.io/en/latest/) and its `getfasta` tool, which allows you to extract sequences from a BED/GTF/VCF file and export them to a FASTA file.

---

#### Using *BioStrings* without *BSGenome*

It is also worth noting that *BioStrings* can be used independently from *BSGenome* with any set of sequences you are able to define in your R environment as an *XString* or *XStringSet* class object. For example, perhaps you are studying the *Amphimedon queenslandica*, a marine sponge organism native to the Great Barrier Reef, and want to explore some basic features of its coding sequences.

<img src="../figures/Amphimedon_queenslandica_adult.png" height="450" width="650"/>

Image source: [Wikipedia](https://en.wikipedia.org/wiki/Amphimedon_queenslandica)

We can retrieve a FASTA file for the coding sequences  (13 overall) from NCBI [(RefSeq ID: NC_008944.1)](https://www.ncbi.nlm.nih.gov/genome/2698) and read the FASTA file into R as a DNAStringSet object using the `readDNAStringSet()` function.
```{r}
fasta.file <- "../data/a.queenslandica.fasta"
a.queen <- readDNAStringSet(fasta.file, "fasta")
a.queen
```

Just as we have done earlier in this lesson, we can again use the *BioStrings* functions to perform basic operations on these sequences. For example:
```r
# confirm how long it is
length(a.queen)

# what is the frequency of each base in your sequence
base.freqs <- alphabetFrequency(a.queen, baseOnly=TRUE, as.prob=TRUE)
base.freqs

# what is the frequency of your favourite base
a.freqs <- letterFrequency(a.queen, "A", as.prob=TRUE)
a.freqs
```

*BioStrings* also implements extensive functionality for **pattern matching**, allowing you to search sequences for specific patterns of interest. For example, we may want to confirm that each of our coding sequences begins with an `ATG` start codon. We can do this using the *BioStrings* functions `matchPattern()` and `countPattern()`.
```r
# return all matches in a DNAString subject sequence to a query pattern
matchPattern("ATG", a.queen[[1]])

# only return how many counts were found
countPattern("ATG", a.queen[[1]])

# what happens if we remove the indexing of the DNAStringSet object 'a.queen'? Why?
matchPattern("ATG", a.queen)

# match a query sequence against multiple sequences that you want to search in  
vmatch <- vmatchPattern("ATG", a.queen)
vmatch

# look at the structure of this object
str(vmatch)

# extract the IRanges for matches in the first subject sequence
vmatch[[1]]
```

We may also have several patterns that we want to search for in each our coding sequences. For example, perhaps we want to search for standard stop codons (`TAG`, `TAA`, `TGA`) in the *A.queenslandica* coding sequences. *BioStrings* functions `matchPDict()` and `vmatchPDict()` provide functionality for such tasks. e.g.
```r
# create a DNAStringSet of the stop codons
stop.codons <- DNAStringSet(c("TAG", "TAA", "TGA"))

# create a dictionary of patterns (PDict class object) that you want to search your sequences for
stop.codons.dict <- PDict(stop.codons)

# search in the first coding sequence
match1 <- matchPDict(stop.codons.dict, a.queen[[1]])
match1
match1[[3]]

# use a loop to search for stop codons in all our coding sequences
matches <- list()
for(i in 1:length(gc.content)){
  matches[[i]] <- matchPDict(stop.codons.dict, a.queen[[i]])
}
length(matches)
str(matches)
matches[[4]]
matches[[4]][[3]]
```

#### Other functionality in *BioStrings*

*BioStrings* also provides functionality for a number of other analytical tasks that you may want to perform on a set of sequences stored using the *XString* and *XStringSet* method, for example:  
* trimming sequence ends based on pattern matching using `trimLRPatterns()`
* local and global alignment problems using `pairwiseAlignment()`
* read in multiple sequence alignments using `readDNAMultipleAlignment()`
* motif searches with a Position Weight Matrix (PWM) using `matchPWM()` (commonly done in ChIP-seq & ATAC-seq)
* palindrome searching using findPalindromes `findPalindromes()`
* computing edit distances between sets of sequences using `stringDist()`  

An excellent *BioStrings* tutorial is available [here](https://bioconductor.org/help/course-materials/2011/BioC2011/LabStuff/BiostringsBSgenomeOverview.pdf) from one of the *BioStrings* creators, that covers much of the same material as we have above, but in more detail and with more complex examples.

> Similar tasks (and many other things) can be performed in python using tools like [*biopython*](https://biopython.org/). The major advantage of the *BioStrings* package is its interoperability with other R/BioConductor packages, such as the *BSGenome* package.
