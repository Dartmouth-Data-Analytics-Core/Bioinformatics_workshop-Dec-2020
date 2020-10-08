# Bioinformatics_workshop

## Day1

- Introduction 
  - Intro to workshop, DAC, and plug the DAC service 
  - How to do this workshop (e.g. managing windows etc. can be reused from that iontro lecture I gave for the RNA-seq workshop) 
   - Personal introductions 
- Introduction to bioinformatics (what is it and how is it done) + intro. to HPC & the Discovery cluster (including having them log in and navigate to files we will be working with) 
- Fred sequencing lecture (45 min with time for discussion?)
- Common file formats for bioinformatics
  - FASTA and FASTQ (exercises for manipulating these files)
  - Principles of alignments and BAM/SAM/CRAM 
  - Genomic region files (GTF/GFF, BED, BED+ examples (e.g. narrowPeak format for ChIP-seq))
- Basic scripting (basic bash scripting, submitting jobs to HPC)
- Installing bioinformatics software using a CLI 
  - Conda and environments 

## Day2

- Reference genomes & genome annotation
  - what are they and where do you find them
  - include exercises of downloading one and then looking at the sequences headers, extracting sequences for an individual chromosome/contig
  - Could mention principles of things like RepeatMasking and n50 in context of microbial genomes?
- IGV tutorial
  - Looking at alignment files
  - Loading in region files (BED files or similar)
  - Do some sort of independent exercise for 20-ish-min (perhaps review evidence for evidence for peak calls from ChIP-seq opr ATAC-seq data)
- Introduction to R (I think we should probably plan to acnowledge this will be a mega short cross course in R, and encourage people as much as possible to do this as much as possible as a prereq and provide many links to more comprhensive R tutorials) 
  - Basic object classes
  - Basic data structures (vector, matrix, list, dataframe)
  - Reading data into R
  - Simple data vis.
  - Installing Packages (CRAN & Bioconductor)
  - Writing data files and figures to files
  - More complicated data struictures anbd specific classes (e.g. S4 objects, summarizedExperiment class as example?)


## Day 3
- Doing bioinformatics with R 
  - GenomicRanges and operating on genomic intervals in R 
  - Importing genomic file types into R (example with a BED file into a GenomicRanges object)
  - Extracting sequences from a reference genome using GenomicRanges objects  
- Statistics and bioinformatics
  - Not sure how to approach this in the context of this workshop, but perhaps this can be a short section introducing basic concepts of how to approach basic statiatsics in R, 
- Using R on an HPC 
  - Running jobs interactively while using a text editor 
  - Submitting R jobs on Discovery 
- Wrap -up/closing remarks
  - Plug coming RNA-seq workshops 


## Particpant number?
I was thinkning we could perhaps go as high as 35 with this one, given that we will have Tim and can recuit more TAs this time around. This also ensures we generqate more revnue from the NOv workshop


## Propsed dates:

**Week of Dec 14th** 

## examples of similar courses (placeholder section)
- [UC Davis bioinfo prereqs workshop](https://ucdavis-bioinformatics-training.github.io/2020-Bioinformatics_Prerequisites_Workshop/). They do 5 full day of 8 hours a day qnd cover a lot of material. 
- 







## Fundamentals of Bioinformatics, December 2020

This workshop will be delivered on xxxx 2020 by the Data Analytics Core (DAC) of the [Center for Quantitative Biology at Dartmouth](https://sites.dartmouth.edu/cqb/). 

The DAC aims to facilitate advanced bioinformatic, computational, and statistical analysis of complex genomics data for the Dartmouth research community. 

If you have questions about this workshop, or would like to discuss data analysis services available from the Data Analytics Core, please visit out [website](https://sites.dartmouth.edu/cqb/projects-and-cores/data-analytics-core/), or email us at: DataAnalyticsCore@groups.dartmouth.edu

<img src="figures/logo (1).jpg" width="250" height="140" >

### Workshop goals: 
- Develop a working understanding what Bioinformatic data analysis involves, how it is done, and what skills it requires
- Gain an appreciation for how next-generation sequencing data is generated (NGS) and how the information generated is stored
- Learn the major file-types used in bioinformatic data analysis and how to manipulate them
- Learn how to install standard bioinformatic software using *Conda*
- Understand the concepts of reference genomes and genome annotations and where to find them 
- Learn how to leverage the *Integrative Genomics Viewer (IGV)* for exploring genomics data 
- Gain a working knowledge of basic programming in R and how it can be used for Bioinformatics 
- Learn how to leverage high performance computing systems (HPCs) to perform Bioinformatic data-analysis 

### Workshop Contacts: 
- Shannon Soucy (Shannon.Margaret.Soucy@Dartmouth.edu)
- Owen Wilkins (omw@Dartmouth.edu)



