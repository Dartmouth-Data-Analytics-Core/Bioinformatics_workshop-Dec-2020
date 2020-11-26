
# Introduction to the Integrative Genomics Viewer (IGV) 

### What is IGV?

The **Integrative Genomics Viewer (IGV)** is a very powerful piece of genomics software produced and distributed by researchers from the Broad Institute at MIT. IGV provides a fast and easy way to explore and visualize genomics data stored in various formats. 

<img src="../figures/igv.png" height="50" width="100"/>

<img src="../figures/igv_example.png" height="50" width="100"/>

File types suppoprted by IGV include:  
* .BAM - alignments  
* .GTF/GFF - genomic features  
* .VCF - variant call format  
* .BED - genmoic regions   
* .BIGWIG - signal tracks  

The range of file formats supported by IGV means it is able to facilitate exploration and visualization of virtually all types of genomics data generated from diverse experimental procedures, for example: 
* Gene expression data (RNA-seq)  
* SNP, mutation and copy number data (WES, WGS)  
* Transcription factor binding sites (ChIP-seq)  
* Hostone modifications (ChIP-seq)  
* Chromatin accessibility (ATAC-seq)  
* Methylation (Bisulfite sequencing)  

The IGV server also hosts a number of reference genomes and annotations, meaning you do not need to load your own genome from a file for many model organisms. You can view the list of hosted genomes on their website [here](http://software.broadinstitute.org/software/igv/Genomes). IGV also provide access to data from large consortia-scale projects such as [*ENCODE*](https://www.encodeproject.org/), [*1000 Genomes*](https://www.internationalgenome.org/home), and [*The Cancer Genome Atlas (TCGA)*](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga). 

If you use IGV in your publications, you should at cite at least the original publication [(found here)](https://www.nature.com/articles/nbt.1754). 

### How do we use IGV?

IGV can be installed and run locally on MacOS, Linux and Windows as a Java desktop application (which is how we will use it today). You should have all downloaded and installed the Desktop version of IGV for the operating system you are working on. 

There is now also an [IGV web-app](https://igv.org/app/) that does not use Java and only needs an internet browser. 

### Learning objectives: 

In todays lesson, we will:  
* Overview of the IGV user interface and basic navigation
* Loading genomes and data into an IGV session
* Discuss and explore how IGV represents data stored across distinct file types 
* Reading in a custom genome
* Saving and restoring an IGV session 

**Note:** This is by no means a comprehensive guide to IGV. Far more functionality exists than we have discussed here, which can be explored in more detail on their website and using the [IGV User Guide](https://software.broadinstitute.org/software/igv/UserGuide). 

### The basic IGV interface and basic navigation

exploring the gene track (+/- strand)


### Loading genomes and data into an IGV session


![](../figures/igv-ex.png)



#### Working BAM files (alignments) in IGV


viewing alignments in igv (.bai indexes)

Example data from 1000G 
Hovering over reads 
gray colors in reads vs colored bases for mapping 
empty reads for bad mapping quality 
paired end reads vs single end reads 
colors for bases different to reference 
insertions and deletions 


#### Working with genomic coordinate data (.BED) and signal tracks (.BIGWIG)




OWEN - NEXT
- pull in anything relevant from rnaseq workshop 
- make all screenshots of things your want to write up and cover 
- fill in text sections describing how that is done 



<img src="../figures/igv_example.png" height="50" width="100"/>






<img src="../figures/igv_example.png" height="50" width="100"/>
<img src="../figures/igv_example.png" height="50" width="100"/>
<img src="../figures/igv_example.png" height="50" width="100"/>
<img src="../figures/igv_example.png" height="50" width="100"/>
<img src="../figures/igv_example.png" height="50" width="100"/>
<img src="../figures/igv_example.png" height="50" width="100"/>







owen -  pull previous igv stuff from rnaseq workshop









comparing tracks in igv, need for normalization 
maybe also subset rnaseq chr 20 bam 


loading in signal track data from chip seq (bigwig)
load in bed files 
defining regions of interest 


loading reference genomes
reading in custom reference genomes 

saving and restoring a session 




what are we not covering - advanced variant evaluation and review in igv 




# create subset vcf file 
bgzip -c ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf > chr20.vcf.gz
tabix -p vcf chr20.vcf.gz
tabix -h chr20.vcf.gz 20:29,875,852-30,548,378 > chr20.sub.vcf

bgzip -c chr20.sub.vcf > chr20.sub.vcf.gz
tabix -p vcf chr20.sub.vcf.gz
tabix chr20.sub.vcf.gz 20:30000000-30001000 | cut -f1-10

# create subset bam file 
samtools view -b HG00099.chrom20.ILLUMINA.bwa.GBR.low_coverage.20130415.bam 20:30096440-30555551 > HG00099.chrom20-sub.low_coverage.bam
samtools index HG00099.chrom20-sub.low_coverage.bam



