#!/usr/bin/env Rscript

#--------------------------------------------------------------------
# Script Name:   subsetting_blast.R
# Description:   Subset BLAST results based on criteria
# Author:        Brandon Monier
# Created:       2019-02-06 10:05:18
# Last Modified:
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to take the transcript
#    (trans) of each gene that has the most `speciesMatches`, and of
#    each of those `speciesMatches` has the:
#
#       1. Smallest `evalues`
#       2. Largest bit score
#       3. rawValue
#       4. percentMatch
#
#    In the end, results file should be very much reduced, with only
#    the best transcript for each gene and its associated metadata.
#--------------------------------------------------------------------

# === Preamble ======================================================

## Load packages
library(tidyverse)

## Import data
blast <- readr::read_csv("subsetting_blast.csv")

## Inspect data
blast



# === Subsetting ====================================================

## Get unique genes
genes <- unique(blast$gene)

## Find all transcripts for each gene
a1 <- blast[grep(pattern = genes[1], x = blast$trans), ]
a2 <- table(a1$trans)
a3 <- names(a2)[which.max(a2)]
a4 <- a1[grep(pattern = a3, x = blast$trans), ]



# === Merritt workflow ==============================================

# ------------------------------
# Read in parsed blast file
# ------------------------------

# Make a play dataset to test subsetting functions
test <- subset(blast, gene=="Zm00001d048284")

library(dplyr)

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



# Trying on 1k genes
onek =blastpar[1:1000,]
onek <- as.data.table(onek)[, count := length(unique(speciesMatch)), by = trans][]

# Specs for later comparisons
#onek = onek[which(onek$count==2),]

onek2 <- onek %>% 
    group_by(gene, trans) %>% 
    filter(percentMatch==max(percentMatch) &
               count==max(count) &
               bitScore==max(bitScore) & 
               rawValue==max(rawValue))

onek2 <- onek %>% 
    group_by(gene) %>% 
    filter(percentMatch==max(percentMatch) &
               count==max(count) &
               bitScore==max(bitScore) & 
               rawValue==max(rawValue))


# For things with only one speciesMatch - FINISHED
tgene <- as.data.table(blast)[, count := length(unique(speciesMatch)), by = gene][]
topblast = tgene[which(tgene$count==1),]
topblast <- topblast %>% 
    group_by(gene) %>% 
    filter(percentMatch==max(percentMatch) & 
               bitScore==max(bitScore) & 
               rawValue==max(rawValue))

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













