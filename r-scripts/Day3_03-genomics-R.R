## 03-genomics-w-R

########################################
library(BSgenome)

available.genomes()

########################################
# use the DNAString contructor function to create a 10 letter DNA sequence
seq <- DNAString(x="AGCT", start=1, nchar=NA)
seq

# how long is it
length(seq)

# show the structure of the DNAString object
str(seq)

########################################
# print DNA alphabet to see what it returns
DNA_ALPHABET

########################################
IUPAC_CODE_MAP

########################################
# randomly sample from specific characters in DNA_ALPHABET to create a longer sequence
seq = sample(DNA_ALPHABET[c(1:4, 16)], size=100, replace=TRUE)
seq

# use the paste command to collapse the characters into a single string
seq = paste(seq, collapse="")
seq

# use the DNAString constructor function to turn this sequence into a DNAString class object
seq.dnastring <- DNAString(seq)
seq.dnastring

########################################
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

########################################
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
# Error in .Call2("DNAStringSet_translate", x, skip_code, dna_codes[codon_alphabet],  : 
# not a base at pos 3

########################################
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

########################################
# assign the genome to a variable using getBSgenome() 
install("BSgenome.Mmusculus.UCSC.mm10")
genome <- getBSgenome("BSgenome.Mmusculus.UCSC.mm10")
genome

# check the structure
str(genome)

# print the metadata for the genome
metadata(genome)
#Error in (function (classes, fdef, mtable)  : 
#unable to find an inherited method for function ‘metadata’ for signature ‘"BSgenome"’

########################################
install("BSgenome.Mmusculus.UCSC.mm10.masked")
genome.m <- getBSgenome("BSgenome.Mmusculus.UCSC.mm10.masked")
class(genome.m)
genome.m

# unmasked genome
class(genome)
genome

########################################
# return basic sequence information summary
seqinfo(genome.m)

# print chromosome 1
genome.m$chr1

# unmasked genome
seqinfo(genome)
genome$chr1

########################################
rm(genome)

genome <- getBSgenome("BSgenome.Mmusculus.UCSC.mm10.masked")

########################################
# assign chr 1
chr1 <- genome$chr1

# confirm how long it is
length(chr1)

# subset it
chr1[1:100]
#Error in chr1[1:100] : object of type 'S4' is not subsettable
chr1[100498:100598]

# what is the frequency of each base in your sequence
alphabetFrequency(chr1, baseOnly=TRUE, as.prob=TRUE)

# what is the frequency of your favourite base
letterFrequency(chr1, "A", as.prob=TRUE)

######################################## 
# we need to establish a vector describing what the extra extended BED columns are
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

# read in peaks
bed <- import("CTCF-forebrain-mm10.bed", extraCols = extraCols_narrowPeak, format="BED")
bed

# extract sequences for peak regions and print to console -first command take a minute or two
ctcf_seqs <- getSeq(genome, bed)
ctcf_seqs

########################################
# calculate nucleotide freqs.
nt_freqs <- alphabetFrequency(ctcf_seqs, baseOnly=TRUE, as.prob=TRUE)

# calculate mean value for nucleotide freqs across all peaks
round(apply(nt_freqs, 2, mean), digits=2)

########################################
hist(width(ctcf_seqs),
     col = "darkgray",
     xlab = "Peak width (bp)",
     main = "CTCF peak width distribution")

########################################
# resize the regions from the BED file
bed_centered <- resize(bed, width = 400, fix = "center")
bed_centered

# check their with
width(bed_centered)

# extract sequences again
ctcf_seqs_cent <- getSeq(genome, bed_centered)
ctcf_seqs_cent

########################################
# export peaks to FASTA file
writeXStringSet(ctcf_seqs, file="CTCF-peaks-resized.fa")

######################################## Started here again 
fasta.file <- "a.queenslandica.fasta"
a.queen <- readDNAStringSet(fasta.file, "fasta")
a.queen

########################################
# confirm how long it is
length(a.queen)

# what is the frequency of each base in your sequence
base.freqs <- alphabetFrequency(a.queen, baseOnly=TRUE, as.prob=TRUE)
base.freqs

# what is the frequency of your favourite base
a.freqs <- letterFrequency(a.queen, "A", as.prob=TRUE)
a.freqs

########################################
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

########################################
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
for(i in 1:length(a.queen)){
  matches[[i]] <- matchPDict(stop.codons.dict, a.queen[[i]])
}
length(matches)
str(matches)
matches[[4]]
matches[[4]][[3]]
