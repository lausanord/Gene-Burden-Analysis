#!/bin/bash

# Use bwToCoverageSummary.py script to calculate coverage summaries from bw files.

set -eu -o pipefail

main() {

    bedFile="sample-input/P02_capture_targets_liftOverTohg19.bed"
    samtools="utils/samtools"
    referenceBuild="GRCh37"
 
    bigwigFileList="sample-outputssh/bwFiles.txt"
    outputPrefix="sample-outputssh/Coverage-summary"

    python3 bwToCoverageSummary.py -bw ${bigwigFileList} -bed ${bedFile} -out ${outputPrefix}


}

main "$@"
