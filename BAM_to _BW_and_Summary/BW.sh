#!/bin/bash

# Use the bamtoBW.sh script to convert BAM to BW with the same prefix.

main() {

    bamFileDir="/mnt/beegfs/home/lausanord/cocorv/bams/BAMS"
    outputDir="/mnt/beegfs/home/lausanord/cocorv/sample-output"
    bedFile="/mnt/beegfs/home/lausanord/cocorv/bams/BAMS/P02_capture_targets_liftOverTohg19.bed"
    samtools="/mnt/beegfs/home/lausanord/cocorv/utils/samtools"
    referenceBuild="GRCh37"

    for bamFile in ${bamFileDir}/*.bam; do
        baseFilename=$(basename "${bamFile}")
        baseFilenameWithoutExt="${baseFilename%.bam}"
        outputPrefix="${outputDir}/${baseFilenameWithoutExt}"

        sbatch bamToBW.sh "${bamFile}" "${outputPrefix}" "${bedFile}" "${samtools}" "${referenceBuild}"
    done
}

main "$@"
