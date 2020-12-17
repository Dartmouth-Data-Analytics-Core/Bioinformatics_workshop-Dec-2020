# Working with genomics data in R/Bioconductor - Part II

## Genome annotation

Commonly in genomic data analysis, we want to provide additional context to our data that helps us addrees our scientific hypothesis. We often achieve this through integrating results from an analysis with publicly available annotation data. Bioconductor provides functionality to directly interface with many popular annotation databases (NCBI, Ensembl, GenBank, UniProt) and import these data into your R environment as Bioconductor type objects (e.g. GRanges).

Examples of common annotation tasks include:  
* Mapping unique gene identifiers (e.g. ENSEMBL or NCBI IDs) to gene symbols in an RNA-seq experiment
* Identifying coding variants from a WGS/WES dataset based on transcriptional context (e.g. coding variants)
* Annotation of genomic context for peak set from ChIP- or ATAC-seq experiment (e.g. promoter, intron, exon, intergenic)

In this lesson, we will introduce the major Bioconductor packages for genome annotation and how you might use them to achieve common tasks in NGS data analysis.

#### Key annotation packages in Bioconductor:

| **Package/package-family** | **Contents & uses**                                                |
|------------------------|-------------------------------------------------------------|
| *AnnotationDbi*      | Methods for accessing data from SQLite-based annotation packages |
| *GenomicFeatures*   | Methods for storing and manipulating transcriptome annotations using the *TxDb* object class  |
| *Org.X.db*               | Gene-based annotation for current genome versions, useful for mapping IDs, symbols and identifiers |
| *EnsDb.X.vX*   | Access to the most up-to-date transcriptome annotations directrly from Ensembl       |
| *biomaRT*   | Release-specific transcriptome annotations from Ensembl for any version or organism     |
| *BS.genome*              | Full sequences for common reference genomes                 |
| *genomation*             | annotation of genomic context and basic data visualization  |

These packages can be broadly categorized into **annotation-centric** packages and **method-centric** packages.

**Annotation-centric** packages such as *Org.X.Db*, *EnsDb.X.vX*, *biomaRT*, and *BS.genome* are designed to provide access to specific annotations, e.g. Ensembl annotations for all organisms from a specific release, or access to legacy Ensembl genome annotations.

**Method-centric** packages such as *AnnotationDbi* and *GenomicFeatures* provide functionality for convienient and efficient access to multiple databases, and do not focus on providing access to any one annotation resource alone. For example, *Org.X.DB*, *EnsDb.X.vX*, and *biomaRT* objects all inherit methods from *AnnotationDbi*, meaning that we can use these common methods to access data from different annotation packages (as we will see in this lesson).

**Note:** Another method-centric package that we won't discuss here is [*AnnotationHub*](https://www.bioconductor.org/packages/release/bioc/html/AnnotationHub.html), which provides methods to query annotation data from a very large range of databases.

### Mapping gene identifiers with *Org.X.DB*

OrgDb represents a family of Bioconductor packages that store gene identifiers from a numerous annotation databases (e.g. gene ontologies) for a wide range or organisms. For example, `org.Hs.eg.db` loads the current annotations for human genes, `org.Mm.eg.db` for mouse, `org.Sc.sgd.db` for yeast, and so on.

OrgDb packages also provide access to some basic additional annotation data such as membership in gene ontology groups. Just as we did in the last lesson, you can navigate through Bioconductor's package index to view all of the existing Org.X.db packages [here](https://www.bioconductor.org/packages/release/BiocViews.html#___OrgDb).

**Note:**OrgDb packages contain annotation data for the most recent genome annotation for that organism, and are therefore not build specific. If you need to annotation for a specific genome build (e.g. h19 vs hg38) you may want to use a different annotation package.


Load the org.db package for the human genome:
```r
library(org.Hs.eg.db)
```

Once loaded, OrgDB allows access to the annotation data for the package you loaded through hidden objects that can be accessed using a set of common Org.Db functions. For example, we can access data from the `org.Hs.eg.db` package after loading it by using the hidden object *org.Hs.eg.db*.
```r
org.Hs.eg.db

# what class is it
class(org.Hs.eg.db)

# what types of data can be extracted?
keytypes(org.Hs.eg.db)
```

OrgDb objects use the `keytypes()` function to access specific types of data from the annotation source. We can ask OrgDb to return a specific keytype that we are interested in.
```r
# obtain all ENSEMBL IDs
entrez.ids <- keys(org.Hs.eg.db, keytype="ENSEMBL")
head(entrez.ids)

# how many are there
length(entrez.ids)
```

The situation we usually end up in however, is when we want to return a set of keytypes given one set of keytypes, for example, returning the gene symbol, entrez ID, and RefSeq ID from a list of Ensembl IDs. For this we can use the `select()` method or `mapIds()` method.
```r
# use ensembl id of first 6 in entrez.ids to get desired keytypes
select(org.Hs.eg.db, keys = head(entrez.ids), columns = c("SYMBOL","ENTREZID", "REFSEQ"), keytype="ENSEMBL")

# using mapIds but only to get gene symbol
mapIds(org.Hs.eg.db, keys = head(entrez.ids), column = c("SYMBOL"), keytype="ENSEMBL")
```

### RNA-seq results annotation using OrgDb  

Lets apply this approach to annotate some real RNA-seq differential expression results. Start by reading in the data, currently stored in `.csv` format.
```r
# read in data
results <- read.csv("diff-exp-results.csv", stringsAsFactors = F, row.names = "ensembl")

# check the first few lines  
head(results)
```

Now use mapIDs to obtain 1:1 mappings between the Ensembl ID and the gene symbol.
```r
# using mapIds but only to get gene symbol
gene.symbols <- mapIds(org.Hs.eg.db, keys = rownames(results), column = c("SYMBOL"), keytype="ENSEMBL")

# add a symbols column to results with the gene symbols we just retreived
results$symbol <- gene.symbols

# check the new results table
head(results)

# make sure there are no NAs in the symbols column we just created
table(is.na(results$symbol))
```

Uh Oh! There are lots of NAs, meaning many genes didn't have a symbol mapped to them... Turns out Org.Db is most built around **entrez IDs** from NCBI and does not contain the annotations for many Ensembl genes, which includes a lot of non-coding RNAs like lincRNAs. Instead, we can use an Ensembl package, `EnsDb.Hsapiens.v86` to pull data directly from Ensembl.

```r
library(EnsDb.Hsapiens.v86)
```

Now lets use `EnsDb.Hsapiens.v86` to retrieve annotation data for our genes and see how many missing genes occur.

```r
# using mapIds but only to get gene symbol
gene.symbols.2 <- mapIds(EnsDb.Hsapiens.v86, keys = head(entrez.ids), column = c("SYMBOL"), keytype="ENSEMBL")

# how long is it
length(gene.symbols.2)

# how many NAs
table(is.na(gene.symbols.2))
```

Many fewer NAs are identified, meaning we were able to annotate more of the genes in our dataset with gene symbols. There are still a few missing though, why might this be?

To ensure we annotate all possible genes, we need to make sure we are using annotation data from the genome annotation used to in the read count quantification process for these data (think back to the `GTF` file we used during alignment and quantification).

These data were annotated using Ensembl verion **97** (which explains why the R-package based off of Ensembl v86 was not able to find matching symbols for all our Ensembl IDs) therefore we could read the GTF file directly into R and manually link ensembl IDs to gene symbols. However, the GTF file is very large and may exceed available memory on our local machines if we load it into R.

Alternatively, we can download annotation data for the human Ensembl annotation releases **97**, and all other Ensembl genomes, using the [BioMart annotation database](https://www.ensembl.org/biomart/martview/6f49181a012124dc5569c6e9d5f71720). BioMart is an easy to use web-based tool that allows you to efficiently obtain annotation data for Ensembl genomes. **Use BioMart to obtain annotation data for the human genome version hg38 from annotation release v97.**

Now read this file into R:
```r
anno <- read.table("GRCh38.p12_ensembl-97.txt", sep="\t", header=TRUE, stringsAsFactors = F)

# check the first few rows and dimensions
head(anno)
dim(anno)

# check how many Ensembl IDs overlap with our results
table(anno$Gene.stable.ID %in% rownames(results))
table(rownames(results) %in% anno$Gene.stable.ID)

# lets rename the ensembl ID column in both datasets so that we can merge them together based on those IDs
colnames(anno)[1] <- "ensembl"
results$ensembl <- rownames(results)

results_merge <- merge(results, anno, by="ensembl")
head(results_merge)
table(is.na(results_merge$Gene.name))
```

Great! We now have gene symbols for all the genes in our dataset, and some additional annotation data integrated directly with our results. Save these data and send it to your PI!
```r
write.csv(results_merge, file = "diff-exp-results-annotated.csv")
```

As we have seen, while the R-packages discussed above can present powerful and quick ways to access lots of annotation data (e.g. gene ontology etc.), there are some obvious limitations which are important to understand when you are annotating your own datasets.

Using BioMart is also valuable if you need annotation data for a model organism that doesn't have an EnsDb or OrgDb R-package availble for it.

### OPTIONAL SECTION: Access Biomart data from within R

If you want to access BioMart from within R, you can use the `BioMart` package to directly interface with the database. This can be a little slower than the approach described above, but can allow more flexibility depending on what you need annotation data for.
```r
library(biomaRt)

# check available martslistMarts(), 10)
head(listMarts())
```

You can see the same marts are listed as were available from the BiomaRt website. We need to choose one of these marts in order to see the available annotation datasets that can be accessed from that mart.
```r
# use 'ensembl to select the 'ENSEMBL_MART_ENSEMBL' mart
ensembl <- useMart("ensembl")

# show available datasets to pull annotation data from
head(listDatasets(ensembl), 10)
tail(listDatasets(ensembl), 10)
nrow(listDatasets(ensembl))

# check for human dataset
table(listDatasets(ensembl)$dataset %in% "hsapiens_gene_ensembl")
```

Now select the `hsapiens_gene_ensembl` dataset from the mart and view the available *attributes* (values that can be returned, e.g. gene symbol) for that dataset.
```r

# pick the ensembl mart for humans  
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# list the attributes for this dataset
head(listAttributes(ensembl), 10)
tail(listAttributes(ensembl), 10)
nrow(listAttributes(ensembl))
```

The flagship function of the BiomaRt package is `getBM()` (for get BiomaRt presumably) which allows us to obtain specific data (attributes) given a set of values that we provide (e.g. Ensembl IDs). The process is very similar to how we used the `select()` and `mapIDs()` functions from OrgDb. Lets use `getBM()` to return annotation data from BiomaRt for our RNA-seq data.
```r
# submit query for 1st 2000 genes (to save time in class) in our RNAseq results
anno_bm <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand"),
                 filters = "ensembl_gene_id",
                 values = head(results$ensembl, 2000),
                 mart = ensembl,
		 useCache = FALSE)
head(anno_bm, 10)
```

You can see we now have a `data.frame` stored in anno_bm in our environment that looks similar to the text file that we downloaded directly from biomaRt and read into R. You could follow a similar protocol to that which we performed to merge the BiomaRt data downloaded in the text file with our RNA-seq results.

---

### Transcript-specific annotation data with Bioconductor

Another core Bioconductor package is the **GenomicFeatures** package, which implements the *TxDb* object class, and provides a convenient and efficient way to store and access transcript specific data from a genome annotation. TxDb objects store a wide-range of transcript-specific information including coordinates and sequences for promoters, exons, introns, and untranslated regions (UTRs).

TxDb objects for common genome annotations can be loaded directly by calling the corresponding annotation package. We can view the available TxDb packages by going to the Bioconductor website and using the search tool. Lets start by loading the TxDb package for the human genome and having a look at the contents.
```r
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# assign to txdb variable
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
txdb
```

You can see this available TxDb object is for gene annotations generated under the UCSC annotation pipeline. What if your genome is not included in the available TxDb packages, for example, if we wanted to continue using an Ensembl annotation? Or perhaps there is no pre-constructed TxDb object available for your organism. GenomicFeatures provides a number of functions specifically to address these issues, for example:  
* `makeTxDbFromEnsembl()` - Construct a TxDb object from Ensembl
* `makeTxDbPackageFromBiomaRt()` - Construct a TxDb object from Biomart
* `makeTxDbPackageFromGFF()` - Construct a TxDb object from a GFF/GTF file

Lets construct a TxDb object from the latest release for human genes from Ensembl. We won't actually build it from scratch right now as it takes a bit of time, but we have a pre-made TxDb ready for you to read into your environment and work with.
```r
#### DO NOT RUN ####
txdb <- makeTxDbFromEnsembl("homo_sapiens", release = 101)

Fetch transcripts and genes from Ensembl ... OK
Fetch exons and CDS from Ensembl ... OK
Fetch chromosome names and lengths from Ensembl ...OK
Gather the metadata ... OK
Make the TxDb object ... OK

#### DO RUN ####
txdb <- loadDb("data/TxDb.Hsapiens.Ensembl.101.db")
txdb
```

Printing the object to the console tells us some basic information about the annotation. For example, you can see the data include hundreds of thousands of rows for unique transcripts, exons, and coding sequences. We can access this information with some basic accessor functions provided by the GenomicFeatures package.
```r
# retrieve all transcript level info
txs <- transcripts(txdb)
txs

# what class is it?
class(txs)

# how long is it
length(txs)

# what is the distribution of transcripts across strands
table(strand(txs))
```

The `transcripts()` function conveniently returns a GRanges class object. This means we can apply all the same methods and accessor functions we used in the previous lesson to the transcript data here (e.g. `seqnames()`, `strand()`, `findOverlaps()`, etc.). There are also several other useful accessor functions that we can use to return specific subsets of the data in our TxDb object.
```r
# retireve all the exon ranges
ex <- exons(txdb)

# retireve all the gene ranges
ge <- genes(txdb)

# promoter ranges for a specified width around TSS
prom <- promoters(txdb, upstream=1000, downstream=0)

# non-overlapping introns or exons
exonicParts(txdb)
intronicParts(txdb)
```

Some of these ranges might be more useful if they were organized by their relation to a specific transcript or gene. There are several accessor functions that provide functionality to achieve this, and return a GRangesList class object rather than ordinary Granges objects.
```r
# return all transcript ranges organized by gene
txs_by_gene <- transcriptsBy(txdb, by = "gene")
txs_by_gene

# index by gene.id of interest to get all transcripts annotated to that gene
txs_by_gene["ENSG00000000419"]

# index by exons by transcript (to identify unique exons)
ex_by_gene <- exonsBy(txdb, by = "tx")
ex_by_gene
```

Equivalent functions exist to return organized GRangesLists for specific features, including:  
* `exonsBy()` - exons by feature
* `intronsByTranscript()` - introns by transcript
* `exonsByTranscript()` - exons by transcript
* `threeUTRsByTranscript()` - 3'UTRs by transcript
* `fiveUTRsByTranscript()` - 5'-UTRs by transcript

Data can also be accessed from a TxDb object using the `select()` method with the `columns` and `keytypes` arguments just as we did for *OrgDBb* objects. This convenient approach is made possible by the fact that *TxDb* objects inheret from *AnnotationDbi* objects, just as *OrgDb* objects do. Using `select` in this way allows us to return data for a large list of features, or a specific subset that we request using the `keys` argument. For example, we might wish to return transcript to gene mapping for specific gene IDs, or we may want to obtain all the exon IDs and their genomic location info for a specific set of transcripts.
```r
# look at the columns avaialble to be returned in the Txdb
columns(txdb)

# return the transcripts annotated to a specific gene of interest
gene_to_tx <- select(txdb, keys = "ENSG00000273696", columns="TXNAME", keytype="GENEID")
gene_to_tx

# return tx to gene mapping for top 500 RNA-seq diff. exp. results
gene_to_tx <- select(txdb, keys = head(rownames(results), 500) ,
                     columns="TXNAME",
                     keytype="GENEID")
head(gene_to_tx)
dim(gene_to_tx)

# check for duplicate entries
table(duplicated(gene_to_tx$GENEID))
table(duplicated(gene_to_tx$TXNAME))

# return exons IDs, their coordinates, and strand for top 10 transcripts from RNA-seq results
tx_to_exon <- select(txdb, keys = head(gene_to_tx, 10)$TXNAME ,
                      columns=c("EXONCHROM", "EXONNAME", "EXONSTART",
                      "EXONEND", "EXONSTRAND", "GENEID"),
                      keytype="TXNAME")

# again, check for duplicate entries
table(duplicated(tx_to_exon$TXNAME))
```

### Example application 1: Variant annotation

Transcript annotation data can be used in many ways. One common usage example is in the annotation of variant calls, where we need to identify the transcriptional context of a variant set (e.g. promoter-associated, exon, intron, untranslated regions, etc.).

To demonstrate how we could achieve this in Bioconductor, we will also use the `VariantAnnotation` package that uses TxDb objects directly to annotate a set of variants. An example set of variant calls is provided, representing variants present over multiple cancer types, identified as part of The Cancer Genome Atlas (TCGA) [*Pan-Cancer Analysis of Whole Genomes (PCAWG) project*](https://www.nature.com/articles/s41586-020-1969-6). Genomic coordinates (hg38) for all identified variants present on chromosome 17 are included in the file `../data/TCGA.pcawg.chr17.bed`.

> Note that *'variant annotation'* also commonly includes annotating variants with their functional consequence

To demonstrate how we could use our TxDb object created above to annotate variants, we will leverage functionality from another BioConductor package, `VariantAnnotation` that uses TxDb objects directly to annotate a set of variants (that can be in GRanges format).

```r
library(VariantAnnotation)

# import the variant locations in bed file format
bed <- import("data/TCGA.pcawg.chr17.bed", format="BED")
bed

# annotate the variants based on our Ensembl Txdb
vars <- locateVariants(bed, txdb, AllVariants())
vars
```

As you can see by printing this object to the console, we now have variants annotated by their transcriptional context, as it relates to the human Ensembl annotation release 101. We can perform some simple operations on this object to explore it further and answer some basic questions, such as how many variants are annotated in each group variuant class.

```r
# sum up variants in each group
sum.tab <- table(vars$LOCATION)
sum.tab

# calculate a quick proprtions table
round(prop.table(sum.tab), digits = 2)

# quick visualization
barplot(round(prop.table(table(coding$LOCATION)), digits = 2))
```

It would also be nice to have the gene symbols included in the TxDb object. We can add gene symbols in the same way we did using above using annotation data downloaded from Biomart. For these data, we need annotation release 101, which has been provided for you in the `Day-3/data/` directory.
```r
#
anno <- read.table("data/GRCh38.p12_ensembl-101.txt", sep="\t", header=TRUE, stringsAsFactors = F)

# return indicies of ENSEMBL geneIDs from variants annotation in the Ensembl v101 annotation data
indicies_of_matches <- match(vars$GENEID, anno$Gene.stable.ID)

# add gene symbols to vars object
vars$GENE.SYMBOL <- anno$Gene.name[indicies_of_matches]
```

Adding gene symbols allows us to easily search for genes of interest, by their transcript ID, gene ID, or gene symbol. We demonstrate this below by restricting to variants identified in the *CD97B* gene.
```r
# exmaple gene of interest:
vars_cd79b <- vars[vars$GENE.SYMBOL %in% "CD97B",]
vars_cd79b

# check how many of each variant type
table(vars_cd79b$LOCATION)
```

We could also use the visualization approaches we learn't in the last lesson to plot the variants in this region using the `Gviz` package.
```r
# required to set expectation for format of chromosome names ('chr17' vs '17')
options(ucscChromosomeNames=FALSE)

# set gene region track from our txdb
txTr <- GeneRegionTrack(txdb,
                        chromosome = "17",
                        start = (min(start(vars_cd79b)) - 500),  
                        end =  (max(start(vars_cd79b) + 500)),
                        name = "Ensembl v101")

# create the annotation track for the variants of interest
track1 <- AnnotationTrack(granges(vars_cd79b), name = "TCGA variants",
                          col.line = "red", fill = "red")

# add the genome axis for scale
gtrack <- GenomeAxisTrack()

# generate the plot
plotTracks(list(gtrack, txTr, track1), main="CD97B variants")
```

<p align="center">
<img src="../figures/cd97b-variants.png" title="xxxx" alt="context"
	width="80%" height="80%" />
</p>

### Example application 2: Peak annotation

Another example usage of how you might use a TxDb object is in the annotation of peak regions from a ChIP-seq experiment. To demonstrate how we could approach this task, we will return to the ChIP-seq data from [Gorkin *et al*, *Nature*, 2020](https://www.nature.com/articles/s41586-020-2093-3) used in the previous lesson, which describes the dynamic chromatin landscape of the developing mouse.

Start by reading the peak regions back in from the narrowpeak files:
```r
# we will use a pre-loaded txdb for mm10 in this example
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

# set txdb to variable
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# set extracols for reading in narrowpeak data
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

# forbrain H3K27ac ChIP-seq peaks
fr_h3k27ac <- rtracklayer::import("forebrain_E15.5_H3K27ac.bed",
                                  format = "BED",
                                  extraCols = extraCols_narrowPeak,
                                  genome = "mm10")

# heart H3K27ac ChIP-seq peaks
ht_h3k27ac <- rtracklayer::import("heart_E15.5_H3K27ac.bed",
                                  format = "BED",
                                  extraCols = extraCols_narrowPeak,
                                  genome = "mm10")

# forbrain H3K9ac ChIP-seq peaks
fr_h3k9ac <- rtracklayer::import("forebrain_E15.5_H3K9ac.bed",
                                 format = "BED",
                                 extraCols = extraCols_narrowPeak,
                                 genome = "mm10")

# heart H3K9ac ChIP-seq peaks
ht_h3k9ac <- rtracklayer::import("heart_E15.5_H3K9ac.bed",
                                 format = "BED",
                                 extraCols = extraCols_narrowPeak,
                                 genome = "mm10")

# combine with H3K27ac peak sets to make GrangesList objects
fr <- GRangesList("h3K27ac" = fr_h3k27ac, "h3K9ac" = fr_h3k9ac)
ht <- GRangesList("h3K27ac" = ht_h3k27ac, "h3K9ac" = ht_h3k9ac)
```

To annotate the genomic context of the ChIP peaks, we will use functionality from the Bioconductor package [*ChIPseeker*](https://bioconductor.org/packages/devel/bioc/manuals/ChIPseeker/man/ChIPseeker.pdf) which provides object classes and methods for ChIP-seq peak annotation and visualization.

The specific function we will use to perform the annotation is the `annotatePeak` function, which accepts a *TxDb* class object directly to define the regions that the peaks should be annotated based on. Lets `annotatePeak` on the forebrain H3K27ac peak set.
```r
# load the chipseeker package
library('ChIPseeker')

# run annotatePeak
fr_h3K27ac_anno <- annotatePeak(fr$h3K27ac, tssRegion=c(-2000, 1000), TxDb = txdb)
fr_h3K27ac_anno

# extract and print the annotation data
fr_h3K27ac_anno <- fr_h3K27ac_anno@anno
fr_h3K27ac_anno

# what class is it
class(fr_h3K27ac_anno)
```

It would be useful if we could run `annotatePeak()` on all samples in one line. We can achieve this using `lapply()`:
```r
annolist <- lapply(list(fr$h3K27ac, ht$h3K27ac, fr$h3K9ac, ht$h3K9ac),
                   annotatePeak,
                   TxDb=txdb,
                   tssRegion=c(-2000, 1000), verbose=FALSE)

# set the names for each element of the list
names(annolist) <- c('Forebrain_H3K27ac', 'Heart_H3K27ac',
                     'Forebrain_H3K9ac', 'Heart_H3K9ac')

annolist
annolist[[1]]
annolist$Forebrain_H3K27ac
```

One way to explore the annotations and compare them across peak sets is to use the `plotAnnoBar()` function from *ChIPseeker*, which plots the proportion of peaks falling into each of the annotation categories.
```r
plotAnnoBar(annolist)
```

<p align="center">
<img src="../figures/chip-anno-example.png" title="xxxx" alt="context"
	width="80%" height="80%" />
</p>

While the proprotion of H3K27ac peaks distributed across the various annotation groups seem relatively stable between the forebrain and heart peak sets, there seems to be a substantially larger proportion of promoter-associated peaks in the H3K9ac peak set from heart tissue compared to that of the forebrain. Perhaps this suggests more transcriptional activity in the heart tissue.

If we were interested in specifically exploring the promoter-associated peaks further on their own, we could subset them.

```r
#extract annotation data for heart h3k9ac
ht_h3K9ac_anno <- annolist$Heart_H3K9ac@anno

# subset for promoter-associated peaks
ht_h3K9ac_anno_promoter <- ht_h3K9ac_anno[ht_h3K9ac_anno$annotation=="Promoter (<=1kb)" |
                                              ht_h3K9ac_anno$annotation=="Promoter (1-2kb)"]
ht_h3K9ac_anno_promoter
```

If we wanted to create a flat file storing the annotated peaks in a sharable file type, we could do this simply by converting the GRanges object to a data frame, and writing that dataframe to a `.csv` file.
```r
# convert GRanges to dataframe
df1 <- as.data.frame(fr_h3K27ac_anno)

# write to csv
write.csv(df1, file = "forebrain_h3K27ac_peaks_annotated_mm10.csv")
```

A far more comprehensive tutorial and description of the ChIPseeker package is available online at the [Bioconductor website](http://bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html).
