library(IRanges); library(GenomicRanges)
############################################################################################
# 1st region
IRanges(start = c(1), width = 4)

# 2nd region 
IRanges(start = c(11), width = 3)
############################################################################################
ir <- IRanges(start = c(1,11), width = c(4, 3))
ir
############################################################################################
# shift all of the regions by a specified offset 
shift(ir, 2)

# resize all regions to only the integer at the center of each region 
resize(ir, fix="center", width=1)

############################################################################################
ir <- IRanges(start = c(1,2,3,3,5,6,7,7,8,11), 
              width = c(4,4,4,4,3,3,3,3,3,3))
ir

############################################################################################
gr <- GRanges(
  seqnames = rep("chr1", 10),
  ranges = IRanges(start = c(1,2,3,3,5,6,7,7,8,11), width = c(4,4,4,4,3,3,3,3,3,3)),
  names = paste0("r", "-", seq(1,10)),
  strand = c(rep("+", 2), rep("-", 3), rep("*", 3), rep("+", 2)),
  score = rnorm(10,5,2))
gr

# return the region ranges only 
granges(gr)

# return the strand info for all regions 
strand(gr)

# return the region names  
names(gr)

# extract the metadata columns 
mcols(gr)

############################################################################################
# calculate coverage of each base over this genomic region 
coverage(gr)

# sum across all regions to get the total coverage for this genomic region 
sum(coverage(gr))

# perhaps we are only interested in the regions on the + strand 
sum(coverage(gr[strand(gr)=="+"]))

############################################################################################
# index GRange object for specific elements 
gr[1]
gr[[1]]
###### Error in getListElement(x, i, ...) : GRanges objects don't support [[, as.list(), lapply(), or unlist() at the moment
gr[1][["r-1"]]
##### Error in getListElement(x, i, ...) : GRanges objects don't support [[, as.list(), lapply(), or unlist() at the moment

# view the top X regions of interest 
head(gr, n=5)

# view the top X regions with scores greater than a value of interest 
head(gr[score(gr)>4], n=5)

############################################################################################
# shift all regions 5bp 
shift(gr, 5)
##### Error in shift(gr, 5) : type 'S4' passed to shift(). Must be a vector, list, data.frame or data.table

# resize all regions by requiring them to be 5bp wide 
resize(gr, 5)

# reduce the regions to one simplified set of non-overlapping regions 
reduce(gr)

############################################################################################
# we need to establish a vector describing what the extra extended BED columns are
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")
# Note: if we had a reglar BED file (no extended fields) we could ignore the extraCols argument

# use the import() function to read in the peaks called in the forebrain H3K27ac ChIP-seq
fr_h3k27ac <- rtracklayer::import("forebrain_E15.5_H3K27ac.bed", 
                                  format = "BED", 
                                  extraCols = extraCols_narrowPeak,
                                  genome = "mm10")


# do the same for the heart H3K27ac ChIP-seq peaks 
ht_h3k27ac <- rtracklayer::import("heart_E15.5_H3K27ac.bed", 
                                  format = "BED", 
                                  extraCols = extraCols_narrowPeak, 
                                  genome = "mm10")

# print both GRanges objects to get an idea for their contents 
fr_h3k27ac
ht_h3k27ac

# check their lengths 
length(fr_h3k27ac)
length(ht_h3k27ac)

############################################################################################
# use findOverlaps() to return matches of genomic ranges between a 'query' and a 'subject'
overlaps <- findOverlaps(query = fr_h3k27ac, subject = ht_h3k27ac)
overlaps

# subset the forebrain GRanges object for H3K27ac peaks that overlap with peaks in heart
fr_h3k27ac_ov1 <- fr_h3k27ac[queryHits(overlaps)]
fr_h3k27ac_ov1

# and vice versa for heart with forebrain
ht_h3k27ac_ov1 <- ht_h3k27ac[subjectHits(overlaps)]
ht_h3k27ac_ov1

# use these objects to calculate the % of overlapping peaks between
length(fr_h3k27ac_ov1)/length(fr_h3k27ac)*100
length(ht_h3k27ac_ov1)/length(ht_h3k27ac)*100

# we could directly subset for the overlapping peaks using subsetByOverlaps()
subsetByOverlaps(fr_h3k27ac, ht_h3k27ac)

# alternatively, we could get the H3K27ac peaks that are unique to each tissue  
#### forebrain
fr_h3k27ac_uniq1 <- fr_h3k27ac[-queryHits(overlaps)]
fr_h3k27ac_uniq1

#### heart
fr_h3k27ac_uniq1 <- fr_h3k27ac[-queryHits(overlaps)]
fr_h3k27ac_uniq1

############################################################################################

# forebrain H3K9ac ChIP-seq peaks
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

# have a look at them
fr
ht

# check their length
length(fr)
length(ht)

# explore individual elements of the list
fr[[1]]
fr[[2]]
length(fr[[1]])
length(fr[[2]])

############################################################################################
# subset for overlapping regions within the forebrain peaks, across both histone marks
fr_overlaps <- findOverlaps(query = fr$h3K27ac, subject = fr$h3K9ac)
fr_overlaps

# subset the forebrain H3K27ac GRanges for peaks overlapping with firebrain H3K9ac peaks
fr_h3k27ac_ov_h3K9ac <- fr$h3K27ac[queryHits(fr_overlaps)]

# calculate % overlapping peaks based on all forebrain H3K27ac peaks
length(fr_h3k27ac_ov_h3K9ac)/length(fr$h3K27ac)*100

# do the same for heart
ht_overlaps <- findOverlaps(query = ht$h3K27ac, subject = ht$h3K9ac)
ht_h3k27ac_ov_h3K9ac <- ht$h3K27ac[queryHits(ht_overlaps)]
length(ht_h3k27ac_ov_h3K9ac)/length(ht$h3K27ac)*100

############################################################################################
fr_ov2 <- subsetByOverlaps(fr$h3K27ac, fr$h3K9ac)
fr_ov2

ht_ov2 <- subsetByOverlaps(ht$h3K27ac, ht$h3K9ac)
ht_ov2

############################################################################################
# create an annotation track from the Granges object for H3K27ac - this takes ~1 minute
fr_h3k27ac_track <- AnnotationTrack(fr$h3K27ac, chromosome = "chr17", start = 9e6, end = 10e6,
                                    name = "Forebrain - H3K27ac", stacking = "dense", col = "indianred")

# do the same for heart H3K27ac - this takes ~ 1 minute
hr_h3k27ac_track <- AnnotationTrack(ht$h3K27ac, chromosome = "chr17", start = 9e6, end = 10e6,
                                    name = "Heart - H3K27ac", stacking = "dense", col = "cornflowblue")

# create a genomic axis object to add to plot
gtrack <- GenomeAxisTrack()

# plot the tracks for this region
plotTracks(list(gtrack, fr_h3k27ac_track, hr_h3k27ac_track), from = 9e6, to = 10e6)

