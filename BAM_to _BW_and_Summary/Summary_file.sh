#!/bin/bash

# Join all the summaries of each chromosome.

# Output file name
outputFile="summary.tsv"

# Output file header
head -n 1 Coverage-summary_chr1.tsv > "$outputFile"

# Concatenate file contents
tail -n +2 -q Coverage-summary_chr*.tsv >> "$outputFile"
