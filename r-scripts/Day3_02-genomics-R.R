# 02-genomics-R

################################################
library(org.Hs.eg.db)

################################################
org.Hs.eg.db

# what class is it
class(org.Hs.eg.db)

# what types of data can be extracted?
keytypes(org.Hs.eg.db)

################################################
# obtain all ENSEMBL IDs
entrez.ids <- keys(org.Hs.eg.db, keytype="ENSEMBL")
head(entrez.ids)

# how many are there
length(entrez.ids)

################################################
# use ensembl id of first 6 in entrez.ids to get desired keytypes
select(org.Hs.eg.db, keys = head(entrez.ids), columns = c("SYMBOL","ENTREZID", "REFSEQ"), keytype="ENSEMBL")

# using mapIds but only to get gene symbol
mapIds(org.Hs.eg.db, keys = head(entrez.ids), column = c("SYMBOL"), keytype="ENSEMBL")

################################################
# read in data
results <- read.csv("diff-exp-results.csv", stringsAsFactors = F, row.names = "ensembl")

# check the first few lines  
head(results)

################################################
# using mapIds but only to get gene symbol
gene.symbols <- mapIds(org.Hs.eg.db, keys = rownames(results), column = c("SYMBOL"), keytype="ENSEMBL")
#'select()' returned 1:many mapping between keys and columns

# add a symbols column to results with the gene symbols we just retreived
results$symbol <- gene.symbols

# check the new results table
head(results)

# make sure there are no NAs in the symbols column we just created
table(is.na(results$symbol))

################################################
library(EnsDb.Hsapiens.v86)

################################################
# using mapIds but only to get gene symbol
gene.symbols.2 <- mapIds(EnsDb.Hsapiens.v86, keys = head(entrez.ids), column = c("SYMBOL"), keytype="GENEID")

# how long is it
length(gene.symbols.2)

# how many NAs
table(is.na(gene.symbols.2))

################################################
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

################################################
write.csv(results_merge, file = "diff-exp-results-annotated.csv")

################################################
library(biomaRt)

# check available martslistMarts(), 10)
head(listMarts())

################################################
# use 'ensembl to select the 'ENSEMBL_MART_ENSEMBL' mart
ensembl <- useMart("ensembl")

# show available datasets to pull annotation data from
head(listDatasets(ensembl), 10)
tail(listDatasets(ensembl), 10)
nrow(listDatasets(ensembl))

# check for human dataset
table(listDatasets(ensembl)$dataset %in% "hsapiens_gene_ensembl")

################################################
# pick the ensembl mart for humans  
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# list the attributes for this dataset
head(listAttributes(ensembl), 10)
tail(listAttributes(ensembl), 10)
nrow(listAttributes(ensembl))

################################################ Left off here - an error with the filter in the getBM command
# submit query for 1st 2000 genes (to save time in class) in our RNAseq results
anno_bm <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand"),
                 filters = "ensembl_gene_id",
                 values = head(results$ensembl, 2000),
                 mart = ensembl,
                 useCache = FALSE)
head(anno_bm, 10)

################################################
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# assign to txdb variable
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
txdb

################################################
# retrieve all transcript level info
txs <- transcripts(txdb)
txs

# what class is it?
class(txs)

# how long is it
length(txs)

# what is the distribution of transcripts across strands
table(strand(txs))

################################################
# retireve all the exon ranges
ex <- exons(txdb)

# retireve all the gene ranges
ge <- genes(txdb)

# promoter ranges for a specified width around TSS
prom <- promoters(txdb, upstream=1000, downstream=0)

# non-overlapping introns or exons
exonicParts(txdb)
intronicParts(txdb)

################################################
# return all transcript ranges organized by gene
txs_by_gene <- transcriptsBy(txdb, by = "gene")
txs_by_gene

# index by gene.id of interest to get all transcripts annotated to that gene - Error: subscript contains invalid names
txs_by_gene["ENSG00000000419"]

# index by exons by transcript (to identify unique exons)
ex_by_gene <- exonsBy(txdb, by = "tx")
ex_by_gene

################################################ - Left off here 
# look at the columns avaialble to be returned in the Txdb
columns(txdb)

# return the transcripts annotated to a specific gene of interest
gene_to_tx <- select(txdb, keys = "ENSG00000273696", columns="TXNAME", keytype="GENEID")
gene_to_tx

# return tx to gene mapping for top 500 RNA-seq diff. exp. results
gene_to_tx <- select(txdb, keys = head(rownames(results), 500) , columns="TXNAME", keytype="GENEID")
head(gene_to_tx)
dim(gene_to_tx)

# check for duplicate entries
table(duplicated(gene_to_tx$GENEID))
table(duplicated(gene_to_tx$TXNAME))

# return exons IDs, their coordinates, and strand for top 10 transcripts from RNA-seq results
tx_to_exon <- select(txdb, keys = head(gene_to_tx, 10)$TXNAME ,
                     columns=c("EXONCHROM", "EXONNAME", "EXONSTART", "EXONEND", "EXONSTRAND", "GENEID"), keytype="TXNAME")
tx_to_exon

# again, check for duplicate entries
table(duplicated(tx_to_exon$TXNAME))

################################################
library(VariantAnnotation)

# import the variant locations in bed file format
bed <- import("data/TCGA.pcawg.chr17.bed", format="BED")
bed

# annotate the variants based on our Ensembl Txdb
vars <- locateVariants(bed, txdb, AllVariants())
vars

################################################
# sum up variants in each group
sum.tab <- table(vars$LOCATION)
sum.tab

# calculate a quick proprtions table
round(prop.table(sum.tab), digits = 2)

# quick visualization
barplot(round(prop.table(table(coding$LOCATION)), digits = 2))

################################################
#
anno <- read.table("data/GRCh38.p12_ensembl-101.txt", sep="\t", header=TRUE, stringsAsFactors = F)

# return indicies of ENSEMBL geneIDs from variants annotation in the Ensembl v101 annotation data
indicies_of_matches <- match(vars$GENEID, anno$Gene.stable.ID)

# add gene symbols to vars object
vars$GENE.SYMBOL <- anno$Gene.name[indicies_of_matches]

# exmaple gene of interest:
vars_erbb2 <- vars[vars$GENE.SYMBOL %in% "ERBB2",]
vars_erbb2

# check how many of each variant type
table(vars_erbb2$LOCATION)

################################################
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

################################################
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

################################################
# run annotatePeak
fr_h3K27ac_anno <- annotatePeak(fr$h3K27ac, tssRegion=c(-2000, 1000), TxDb = txdb)
fr_h3K27ac_anno

# extract and print the annotation data
fr_h3K27ac_anno <- fr_h3K27ac_anno@anno
fr_h3K27ac_anno

# what class is it
class(fr_h3K27ac_anno)

################################################
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

################################################
plotAnnoBar(annolist)

################################################
#extract annotation data for heart h3k9ac
ht_h3K9ac_anno <- annolist$Heart_H3K9ac@anno

# subset for promoter-associated peaks
ht_h3K9ac_anno_promoter <- ht_h3K9ac_anno[ht_h3K9ac_anno$annotation=="Promoter (<=1kb)" |
                                            ht_h3K9ac_anno$annotation=="Promoter (1-2kb)"]
ht_h3K9ac_anno_promoter

################################################
# convert GRanges to dataframe
df1 <- as.data.frame(fr_h3K27ac_anno)

# write to csv
write.csv(df1, file = "forebrain_h3K27ac_peaks_annotated_mm10.csv")
