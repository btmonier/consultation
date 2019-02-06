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


