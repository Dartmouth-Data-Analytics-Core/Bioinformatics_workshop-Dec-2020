# Working with NGS data Part III

After aligning reads to a reference genome, there are a couple of different ways we can continue with the analysis, which are dependent on the data type.

**Read quantification:**
- When working with RNA-seq data for example, we are often interested in counting how many reads overlap each gene. This process, often referred to as *quantification*, allows us to infer the expression levels of individual genes.

**Peak calling:**
- Alternatively, in experiments where we have performed some sort of enrichment for genomic regions of interest (e.g. DNAase-seq, ChIP-seq, ATAC-seq), we are usually interested in identifying regions containing signal above some background level. This is referred to as *peak calling* and allows us to confidently identify the regions containing enriched signals, e.g. transcription-factor binding sites in a ChIP-seq experiment.

In this lesson, we will briefly explore the fundamental concepts of *read quantification* and *peak calling*, while introducing useful software and relevant file formats for each.

---

## Read count quantification

For most downstream analyses in RNA-seq, especially differential expression, we care about how many reads aligned to a specific gene, as this tells us about the genes expression level, which we can then compare to other samples. Inherently, this means that we want to make these data count-based, so that we can use statistical models to compare these counts between experimental conditions of interest.

<p align="center">
<img src="../figures/genoic-content.png" title="xxxx" alt="context"
	width="100%" height="100%" />
</p>

> Although RNA-seq is the most common scenario where read counting is performed, read counting is relevant in the analysis of other genomic data types. For example, in ChIP-seq we often want to perform a differential binding analysis between two or more conditions, which requires us to count how many reads overlap each of our called peaks.

Read quantification methods generally require two inputs:  
- an alignment file (.bam)
- a set of features over which to count (e.g. GTF/GFF).

Recall that a GTF/GFF file is used to store genome annotation data, therefore contains the coordinates over all of the exons that we want to count reads.

![](../figures/gtf.png)

The most simplistic methods (e.g. [htseq-count](https://htseq.readthedocs.io/en/release_0.11.1/count.html), [featureCounts](http://subread.sourceforge.net/)) use a specific set of rules to count the number of reads overlapping specific features. These are a good choice if your data is less complex, e.g. 3'-end data. More complex methods such as [RSEM](https://deweylab.github.io/RSEM/)), determines the probability that a read should be counted for a particular feature.

As an example, lets use [htseq-count](https://htseq.readthedocs.io/en/release_0.11.1/count.html) to quantify reads for an alignment we created in the previous lesson. Some important options in *htseq-count* include:

**Feature type (`-t`):**  
Specifies the feature in your GTF file you want to count over (r3rd column). The default is **exon**. However, this can be changed to any feature in your GTF file, so theoretically can be used to count any feature you have annotated.

**Strandedness (`-s`):**  
Specifies if reads in your experiment come from a stranded (`yes`) or unstranded (`no`) library type. It is critical to set this correctly, as incorrect selection will result in needlessesly throwing away 50% of your reads.  

```r
htseq-count \
	-f bam \
	-s no \
	-r pos \
	../alignment/SRR1039508.Aligned.sortedByCoord.out.bam \
	../Homo_sapiens.GRCh38.97.chr20.gtf > SRR1039508.htseq-counts
```

Have a look at the resulting file.
```bash
# how many lines
wc -l SRR1039508.htseq-counts

# first few rows
head SRR1039508.htseq-counts

# importantly, lets check the last few rows as these contain some important info
tail -n 12 SRR1039508.htseq-counts
```

This process can be repeated for each sample in your dataset, and the resulting files compiled to generate a matrix of raw read counts that serve as input to downstream analysis (e.g. differential expression or binding analysis).


## Peak calling

In contrast to RNA-seq data, other types of genomics data exist for which we wish to first identify regions of the reference genome containing increased signal above background. For example, enriched regions in a ChIP-seq experiment are indicative of

While peak calling is typically associated with ChIP-seq, it is utilized in a growing number of genomic analysis workflows, especially in more recent years as the number of technologies being designed to profile various genomic features grows rapidly.

In order to accurately identify peaks in signal intensity that truly represent the signal of interest, we need some sort of statistical model that can compare the distribution of reads between samples and determine where signal truly exceeds the background in a robust and quantitative way.



<p align="center">
<img src="../figures/peak-calling.png" title="xxxx" alt="context"
	width="100%" height="100%" />
</p>








<p align="center">
<img src="../figures/bed.png" title="xxxx" alt="context"
	width="100%" height="100%" />
</p>




bed extension formats (narrowpeak, broadpeak)



Read quantification can then be performed on these regions and be used as input to a differential binding analysis (ChIP-seq) or a differential accessibility analysis (ATAC-seq).




In analyses where you have identified a set of called peaks with signal significantly above the background level, it can be useful to visualize extent of the signal in the peak regions, in order to gain an idea of how enriched above background the signal in those regions was.

The best way to achieve this is to convert our alignment files into a **bigWig* file. bigWig files are an indexed binary file format used to store continuous data (i.e. signal) over a set of genomic positions (e.g chromosomes).


<p align="center">
<img src="../figures/bigwig.png" title="xxxx" alt="context"
	width="100%" height="100%" />
</p>

As bigwig files contain signal enrichment data across all positions in our reference genome, we can use them to evaluate signal across genomic loci simultaneously. One common example is to evaluate signal in the regions directly upstream and downstream of transcriptional start sites (TSSs).

<p align="center">
<img src="../figures/bigwig.png" title="xxxx" alt="context"
	width="100%" height="100%" />
</p>

Again, the signal is not limited to ChIP-seq and TFBS, and could instead represent Tn5 insertions in an ATAC-seq experiment for example.

> For the purposes of this workshop, we only need understand the idea behind of bigwig files and what they are used for. We hope to address their generation and use in more detail in future workshops.
