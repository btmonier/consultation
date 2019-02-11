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
library(dplyr)
library(data.table)
library(ggplot2)
library(microbenchmark)
library(tibble)

## Load data
onek_blast <- read.csv(
    file = "data/blast4Genomes_parsed.csv",
    stringsAsFactors = FALSE
)



# ==== Data munging =================================================

## Get species count (why use `data.table`?)
tgene <- as.data.table(onek_blast)[, count := length(unique(speciesMatch)), by = trans][]

## Add count column to onek_blast
onek_blast$noSpeciesMatchesPerTranscript = tgene$count



# === Deployment (all genes) ========================================

## Iterate through genes (unique) - using `onek` object
blast_ls <- list()
for (i in seq_along(unique(onek_blast$gene))) {
    tmp <- onek_blast[which(onek_blast$gene == unique(onek_blast$gene)[i]), ]
    blast_ls[[i]] <- tmp[which(
        max(tmp$noSpeciesMatchesPerTranscript) &&
            max(tmp$percentMatch) &&
            max(tmp$bitScore)
    ), ]
}
blast_df <- data.frame(
    matrix(
        unlist(blast_ls),
        nrow = length(unique(onek_blast$gene)),
        byrow = T
    )
)
colnames(blast_df) <- colnames(onek_blast)



# ==== Make pipeline more efficient =================================

## Use `lapply()` instead of a for loop
blast_ls <- lapply(X = seq_along(unique(onek_blast$gene)), FUN = function(i) {
    tmp <- onek_blast[which(onek_blast$gene == unique(onek_blast$gene)[i]), ]
    blast_ls <- tmp[which(
        max(tmp$noSpeciesMatchesPerTranscript) &&
            max(tmp$percentMatch) &&
            max(tmp$bitScore)
    ), ]
    return(blast_ls)
})
blast_df <- data.frame(
    matrix(
        unlist(blast_ls),
        nrow = length(unique(onek_blast$gene)),
        byrow = T
    )
)
colnames(blast_df) <- colnames(onek_blast)



# === Time it =======================================================

## Microbenchmark
times <- microbenchmark::microbenchmark(
    lapply = function() {
        blast_ls <- lapply(X = seq_along(unique(onek_blast$gene)), FUN = function(i) {
        tmp <- onek_blast[which(onek_blast$gene == unique(onek_blast$gene)[i]), ]
        blast_ls <- tmp[which(
            max(tmp$noSpeciesMatchesPerTranscript) &&
                max(tmp$percentMatch) &&
                max(tmp$bitScore)
        ), ]
        return(blast_ls)
        })
    },
    for_loop = function() {
        blast_ls <- list()
        for (i in seq_along(unique(onek_blast$gene))) {
            tmp <- onek_blast[which(onek_blast$gene == unique(onek_blast$gene)[i]), ]
            blast_ls[[i]] <- tmp[which(
                max(tmp$noSpeciesMatchesPerTranscript) &&
                    max(tmp$percentMatch) &&
                    max(tmp$bitScore)
            ), ]
        }
    },
    times = 1000L
)

## Plot the class with ggplot2
ggplot2::autoplot(times)
