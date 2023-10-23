#!/bin/bash

# Restrict the analysis to the region where >= 90% of samples have coverage of depth >=10, for cases and controls.

samples_coverage="samples.coverage10x.bed.gz"
controls_coverage="gnomAD.coverage10x_chr.bed.gz"
outputFile="intersect.coverage10x.bed.gz"

# bedtools intersect and bedtools merge
bedtools intersect -sorted -a "$samples_coverage" -b "$controls_coverage" | bedtools merge -i stdin | gzip > "$outputFile"
