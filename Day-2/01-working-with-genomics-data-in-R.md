# Part 1 - Wokring with Genomics Data in R

### Learning objectives: 
- Understand why we might want to work with genomics data in R, and what packages are available to do this
- Learn the basics behind the GenomicRanges format and how it can be used to store information on genomic regions
- Learn how to work with reference genomes in R from BioConductor packages 
- Learn how to read in and manipulate a refernece genome in R 


### The BioStrings R-package

[BioStrings](http://bioconductor.org/packages/release/bioc/html/Biostrings.html) is a powerful R-package that lies at the core of several other BioConductor packages, and facilitates the storing sequence data in R so that it can be efficienty and quickly viewed, searched, and maniuplated. BioStrings achieves this through implementing several unique object classes for storing sequences, and asets of sequences, in the R-environment. 



### The GenomicRanges R-package

[GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) is an extremely useful package from the BioConductor project for working with genomic regions and coordinates in R, and lies at the core of numerous other BioConductor packages such as [BSgenome](https://bioconductor.org/packages/release/bioc/html/BSgenome.html), [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html) and [VariantAnnotation](https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html). 

The goal of this lesson is to introduce you to the major object classes used by *GenomicRanges*, how they work, and together explore some examples of how you could use the package in your own work. This is not mean't to be a comprehensive introduction to the complete functionality of *GenomicRanges*, or replace the excellent tutorials or Vigenttes available on the [BioConductor website](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)

Lets load the GenomicRanges package: 
```{r}
library(GenomicRanges)
```



### The BSGenome R-package

The [BSgenome](https://bioconductor.org/packages/release/bioc/html/BSgenome.html) package provides BioConductor-based access to reference genome sequences that utilize the object classes of the *BioStrings* package discussed above (e.g. *DNAString* and *DNAStringSet*). 

Check which genomes are currently available from BSGenome. 
```{r}
available.genomes()
```

Lets use the `getBSgenome()` function to retrieve the human genome (build *hg38* from UCSC). 
```{r}
genome <- getBSgenome("BSgenome.Hsapiens.UCSC.hg38")

# check the structure 
str(genome) 
```

**Usage note:** In the `getBSgenome()` function, the `masked` argument is set to `FALSE` by default. This means that repetitive or otherwise masked sequences that can be problematic in a number of applications (e.g. peak calling) will not be masked (represented as `N` nucleotides). If you would like these sequences to be masked, you should use `masked=TRUE`. 




If your genome is not included in the available genomes but you would still like to leverage the BioStrings and BSGenome framework, you can [forge a BSGenome package](https://bioconductor.org/packages/release/bioc/html/BSgenome.html) following instructions available at the BioConductor website. 




