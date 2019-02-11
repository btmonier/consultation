# ------------------------------------------------------
# Merritt Burch
#
# mbb262@cornell.edu
#
# 2018-11-07
#
# Script to parse three genome BLAST results
# for summary statstics,
# This version has updated gene IDs that include
# the start and stop positions
#
# Input: blast tabular output from alignment of
#        B73 CDSs with Brachy, sorghum, and rice CDSs
#
# Output: .txt file containing most conserved CDSs
#         .bed file containig testable ref-ranges for PHG
#
# Flavor: Positions of ATG-start codon and terminator with
#         merged overlaps
# --------------------------------------------------------

# Set working directory
setwd("~/Box Sync/Cornell_PhD/PHG/blast_b73cds_refranges/data/blast4Genomes")

# Windows
setwd("C:/Users/merri/Box Sync/Cornell_PhD/PHG/blast4Genomes")

# Load file containg positions of start codon - stop codon
blast <- read.table("matches_tabular_option6_2019_02_01.txt", header=FALSE)


# ------------------
#  Parse out file
# ------------------

#Separate first column into general gene, chromosome, start,stop
# Not pretty but this will work
blast$gene = blast$V1
blast$chr = blast$V1
blast$start = blast$V1
blast$stop = blast$V1

# Format
blast$gene = gsub("-[0-9]:.*", "", blast$gene)
blast$gene = gsub("-[0-9][0-9]:.*", "", blast$gene)
#blast$gene = gsub("-B73V4_.*", "", blast$gene)

# Remove rows of mitochondrial and chloroplast genes
blast <- blast[!grepl("-Mt:", blast$gene),]
blast <- blast[!grepl("-Pt:", blast$gene),]
blast <- blast[!grepl("B73V4_ctg", blast$gene),]

# More general formatting
blast$transcript = blast$gene
blast$gene = gsub("_T[0-9][0-9][0-9]", "", blast$gene)
blast$gene = gsub("_T[0-9][0-9]", "", blast$gene)
blast$chr = gsub(".*[0-9]-", "", blast$chr)
blast$chr = gsub(":.*", "", blast$chr)

# Format start position
#blast$start = gsub(".*_ctg[0-9][0-9][0-9]:", "", blast$start)
#blast$start = gsub(".*_ctg[0-9][0-9]:", "", blast$start)
#blast$start = gsub(".*_ctg[0-9]:", "", blast$start)
blast$start = gsub(".*-[0-9]:", "", blast$start)
blast$start = gsub(".*-10:", "", blast$start)
blast$start = gsub(":.*", "", blast$start)

# Format stop position
blast$stop = gsub(":-1$", "", blast$stop)
blast$stop = gsub(":1$", "", blast$stop)
blast$stop = gsub(".*:", "", blast$stop)

# Rearrange
blast <- blast[,c(11,15,12,13,14,6,7,8,9,10,4)]

# Rename columns
colnames(blast) <- c("gene", "trans", "chr", "five_prime_UTR", "three_prime_UTR",
                     "evalue", "bitScore", "rawValue",
                     "percentMatch", "other", "proteinMatch")

# ---------------------------------------------
# Find what species it mapped to by protein ID
# ---------------------------------------------
# KQK matches both
# Make new column from protein matches to simplify into species
blast$speciesMatch <- blast$proteinMatch

# Rearrange again
blast <- blast[,c(1:2, 12,3:11)]

# Gsub out protein ID matches with species
# Rice only
blast$speciesMatch = gsub("Os.*", "rice", blast$speciesMatch) #rice
blast$speciesMatch = gsub("BAC.*", "rice", blast$speciesMatch) #rice
blast$speciesMatch = gsub("CAA.*", "rice", blast$speciesMatch) #rice

# Sorghum only
blast$speciesMatch = gsub("EER.*", "sorghum", blast$speciesMatch) #only sorghum
blast$speciesMatch = gsub("EES.*", "sorghum", blast$speciesMatch) # only sorghum
blast$speciesMatch = gsub("OQU.*", "sorghum", blast$speciesMatch)
blast$speciesMatch = gsub("KXG.*", "sorghum", blast$speciesMatch)

# bracy only
blast$speciesMatch = gsub("KQJ.*", "brachypodium", blast$speciesMatch)
blast$speciesMatch = gsub("PNT.*", "brachypodium", blast$speciesMatch)
blast$speciesMatch = gsub("KQK10.*", "brachypodium", blast$speciesMatch)
blast$speciesMatch = gsub("KQK0.*", "brachypodium", blast$speciesMatch)
blast$speciesMatch = gsub("KQK2.*", "brachypodium", blast$speciesMatch)
blast$speciesMatch = gsub("KQK1.*", "brachypodium", blast$speciesMatch)
blast$speciesMatch = gsub("PNS.*", "brachypodium", blast$speciesMatch)

# setaria only
blast$speciesMatch = gsub("KQL.*", "setaria", blast$speciesMatch)
blast$speciesMatch = gsub("KQK9.*", "setaria", blast$speciesMatch)
blast$speciesMatch = gsub("KQK8.*", "setaria", blast$speciesMatch)

# Save this file and write to csv
write.csv(blast, file="blast4Genomes_parsed.csv", row.names = FALSE)


# ------------------------------
# Read in parsed blast file
# ------------------------------

# Load in file
blastpar <- read.csv("blast4Genomes_parsed.csv")

# Load in required packages
library(dplyr)
library(data.table)

# Count number of species matches for each transcript
tgene <- as.data.table(blastpar)[, count := length(unique(speciesMatch)), by = trans][]

# Add count column to blastpar
blastpar$noSpeciesMatchesPerTranscript = tgene$count



t <- test %>%
    group_by(trans, speciesMatch) %>%
    filter(percentMatch==max(percentMatch) &
               bitScore==max(bitScore) &
               rawValue==max(rawValue))

# ----
# Try to subset things the M&B way
# ----
tgene <- as.data.table(blast)[, count := length(unique(speciesMatch)), by = trans][]
t <- test %>%
    group_by(trans, speciesMatch) %>%
    filter(percentMatch==max(percentMatch) &
               bitScore==max(bitScore) &
               rawValue==max(rawValue))
t <- t %>%
    group_by(gene) %>%
    filter(percentMatch==max(percentMatch) &
               bitScore==max(bitScore) &
               rawValue==max(rawValue))


# -------------------------------
# Brandon look at this snippet
# -------------------------------

# Trying on 1k genes
onek <- tibble::as_tibble(blastpar[1:1000,])
onek$count <- onek[, count := length(unique(speciesMatch)), by = onek$trans]


onek <- as.data.table(onek)[, count := length(unique(speciesMatch)), by = trans][]

# Specs for later comparisons
onek = onek[which(onek$count==2),]

# Make a play dataset to test subsetting functions
test <- rbind((subset(blastpar, gene== "Zm00001d027450")),(subset(blastpar, gene== "Zm00001d027450")))

# Takes all top transcripts
onek2 <- test %>%
    group_by(trans) %>%
    filter(percentMatch==max(percentMatch) &
               bitScore==max(bitScore) &
               noSpeciesMatchesPerTranscript==max(noSpeciesMatchesPerTranscript))

# Doesn't work, takes only 1 genes with all transcripts,
#`    not top transcript for each gene
onek2 <- onek2 %>%
    dplyr::group_by(gene) %>%
    dplyr::filter(percentMatch==max(percentMatch) &&
                      bitScore==max(bitScore))


# === Brandon - iterate =============================================

# Load in required packages
library(dplyr)
library(data.table)
library(tibble)

# Count number of species matches for each transcript
blastpar <- read.csv("blast4Genomes_parsed.csv", stringsAsFactors = FALSE)
tgene <- as.data.table(blastpar)[, count := length(unique(speciesMatch)), by = trans][]

# Add count column to blastpar
blastpar$noSpeciesMatchesPerTranscript = tgene$count
onek <- tibble::as_tibble(blastpar[1:1000, ])

## Iterate through genes (unique) - using `onek` object
blast_ls <- list()
for (i in seq_along(unique(onek$gene))) {
    tmp <-  onek[which(onek$gene == unique(onek$gene)[i]), ]
    blast_ls[[i]] <- tmp[which(
        max(tmp$noSpeciesMatchesPerTranscript) &&
            max(tmp$percentMatch) &&
            max(tmp$bitScore)
    ), ]
}
blast_df <- data.frame(
    matrix(
        unlist(blast_ls),
        nrow = length(unique(onek$gene)),
        byrow = T
    )
)
colnames(blast_df) <- colnames(onek)








# --------------------------------------------------------
m = onek2
m$col4 = names(m)[max.col(m,ties.method="first")]

#The columns we want to take the max of
cols <- c("percentMatch","bitScore", "rawValue", "gene")
m$col4 <- names(m)[names(m) %in% cols][max.col(m[,cols],ties.method="first")]


onek3 <- onek2 %>%
    group_by(gene) %>%
    filter(percentMatch==max(percentMatch) &
               bitScore==max(bitScore) &
               rawValue==max(rawValue))
refRangegenes2 <- onek2[!duplicated(onek2$gene),]

length(unique(test$gene))
length(unique(onek2$gene))


# For things with only one speciesMatch - FINISHED
tgene <- as.data.table(blastpar)[, count := length(unique(speciesMatch)), by = gene][]
topblast = tgene[which(tgene$count==1),]

topblast <- topblast %>%
    group_by(gene) %>%
    filter(percentMatch==max(percentMatch) &
               bitScore==max(bitScore) &
               rawValue==max(rawValue))
library(plyr)
gdfg = ddply(topblast,.(gene),nrow)


# -------------------------------
# Count number of genes matched
# -------------------------------

# Counts number of maize genes that successfuly mapped
# to 1-4 species

library(data.table)
table2 <- as.data.table(blast)[, count := length(unique(speciesMatch)), by = gene][]
tef = table2[,c(1,13)]
tef = subset(tef,!duplicated(tef$gene))
hist(tef$count)

sum(tef$count==4)
sum(tef$count==3)
sum(tef$count==2)
sum(tef$count==1)

library(ggplot2)
count_matches = c(sum(tef$count==4),sum(tef$count==3),
                  sum(tef$count==2),sum(tef$count==1))
combined_match_values = cbind(4:1, count_matches)
colnames(combined_match_values) = c("speciesMatched", "numberMatched")
combined_match_values=as.data.frame(combined_match_values)
combined_match_values$speciesMatched = as.factor(combined_match_values$speciesMatched)

ggplot(combined_match_values, aes(x = speciesMatched, y = numberMatched)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label=numberMatched), vjust=-0.3, size=4) +
    scale_x_discrete(name = "Number of BLAST hits to other queried species") +
    scale_y_discrete(name = "Number of maize CDSs matched to queried species") +
    theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))

# -------------------
# Filter out results
# -------------------

# Take top BLAST hit with highest bit score & percent match, lowest evalue
# Subsets by general gene
library(dplyr)
refRangegenes <- blast %>%
    group_by(gene) %>%
    filter(evalue==min(evalue) & bitScore==max(bitScore))

# Sort data by gene ID and transcript ID
refRangegenes <- refRangegenes[order(refRangegenes$gene, refRangegenes$trans),]

# Subset out only unique genes (because many transcripts look
#   to be duplicated)
refRangegenes2 <- refRangegenes[!duplicated(refRangegenes$gene),]

# Export only transcript IDs
transcIds <- as.data.frame(refRangegenes2$trans)
write.table(transcIds, file = "B73_BLAST_ConservedTranscriptIDs.txt",
            row.names = FALSE, col.names = FALSE)


# --------------------------------
# Look for overlapping ref-ranges
# --------------------------------

# Create testing data.frame
ranges = refRangegenes2

# Use genomic ranges to merge overlapping ranges
library(GenomicRanges)

# Make genomic range format from data.frame
refRangegenes2 = makeGRangesFromDataFrame(refRangegenes2,
                                          start.field = "five_prime_UTR_start",
                                          end.field = "three_prime_UTR_stop",
                                          seqnames.field = "chr",
                                          keep.extra.columns = TRUE)

# Reduce ref-ranges (merge overlapping ref-ranges)
refRangegenes2= reduce(refRangegenes2)

# Turn this into a data frame and extract relavent columns
refRangegenes2 = as.data.frame(refRangegenes2)


# ------------------------------------------------------
# Get CDS start/stop positions for most conserved genes
# ------------------------------------------------------

# Format to Lynn's .bed format
bed <- refRangegenes2[,1:3]
bed$start <- as.numeric(bed$start)
bed$start <- bed$start - 1

# Sort data again
bed <- bed[order(bed$chr),]

# Write tab delim file
write.table(bed, file = "mbb262_UTRsIncluded_refRanges_mergedOverlapping.bed",
            row.names = FALSE, col.names = FALSE, quote = FALSE,
            sep = "\t")


