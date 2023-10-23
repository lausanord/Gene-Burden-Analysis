#!/bin/bash

set -u -o pipefail

main() {
  # converts samtools depth output to bedGragh and then bedGRaph to bigWig 

  bamFile=$1          # input bam file
  outputPrefix=$2     # the output prefix, bigwig output will be 
                      # ${outputPrefix}.bw
  bedFile=$3          # the bed file that specifies the region to calculate the 
                      # covearge
  SAMTOOLS=$4         # the path of samtools executable
  referenceBuild=$5   # reference build, could be GRCh37 or GRCh38
  
  filename=$(basename -- "$bamFile")
  extension="${filename##*.}"
  
  cramOption=""
  if [[ ${extension} == "cram" ]]; then
    cramOption="--reference ${referenceBuild}"
  fi

  minBaseQuality=10
  minMappingQuality=20

  outputBedGraph=${outputPrefix}.bedgraph
  outputBigWig=${outputPrefix}.bw

  SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
  depthToBedGraphTool="${SCRIPT_DIR}/utils/depth2bedgraph.awk"
  bedGraphToBigWigTool="${SCRIPT_DIR}/utils//bedGraphToBigWig"
  chromSizesGRCh38="${SCRIPT_DIR}/utils/hg38.chrom.sizes"
  chromSizesGRCh37="${SCRIPT_DIR}/utils/hg19.chrom.sizes"

  ${SAMTOOLS} depth ${cramOption} \
    -q ${minBaseQuality} -Q ${minMappingQuality} -b ${bedFile} ${bamFile} \
    | awk -f ${depthToBedGraphTool} > ${outputBedGraph}

  if grep -qwr '^chr1' ${outputBedGraph}; then
    grep -wr '^chr1' ${outputBedGraph} > ${outputBedGraph}.chr1
  else 
    grep -wr '^1' ${outputBedGraph} > ${outputBedGraph}.chr1
  fi

  if grep -qwr '^chr2' ${outputBedGraph}; then
    grep -wr '^chr2' ${outputBedGraph} > ${outputBedGraph}.chr2
  else
    grep -wr '^2' ${outputBedGraph} > ${outputBedGraph}.chr2
  fi

  if grep -qwr '^chr3' ${outputBedGraph}; then
    grep -wr '^chr3' ${outputBedGraph} > ${outputBedGraph}.chr3
  else
    grep -wr '^3' ${outputBedGraph} > ${outputBedGraph}.chr3
  fi

  if grep -qwr '^chr4' ${outputBedGraph}; then
    grep -wr '^chr4' ${outputBedGraph} > ${outputBedGraph}.chr4
  else
    grep -wr '^4' ${outputBedGraph} > ${outputBedGraph}.chr4
  fi

  if grep -qwr '^chr5' ${outputBedGraph}; then
    grep -wr '^chr5' ${outputBedGraph} > ${outputBedGraph}.chr5
  else
    grep -wr '^5' ${outputBedGraph} > ${outputBedGraph}.chr5
  fi

  if grep -qwr '^chr6' ${outputBedGraph}; then
    grep -wr '^chr6' ${outputBedGraph} > ${outputBedGraph}.chr6
  else
    grep -wr '^6' ${outputBedGraph} > ${outputBedGraph}.chr6
  fi

  if grep -qwr '^chr7' ${outputBedGraph}; then
    grep -wr '^chr7' ${outputBedGraph} > ${outputBedGraph}.chr7
  else
    grep -wr '^7' ${outputBedGraph} > ${outputBedGraph}.chr7
  fi

  if grep -qwr '^chr8' ${outputBedGraph}; then
    grep -wr '^chr8' ${outputBedGraph} > ${outputBedGraph}.chr8
  else
    grep -wr '^8' ${outputBedGraph} > ${outputBedGraph}.chr8
  fi

  if grep -qwr '^chr9' ${outputBedGraph}; then
    grep -wr '^chr9' ${outputBedGraph} > ${outputBedGraph}.chr9
  else
    grep -wr '^9' ${outputBedGraph} > ${outputBedGraph}.chr9
  fi

  if grep -qwr '^chr10' ${outputBedGraph}; then
    grep -wr '^chr10' ${outputBedGraph} > ${outputBedGraph}.chr10
  else
    grep -wr '^10' ${outputBedGraph} > ${outputBedGraph}.chr10
  fi

  if grep -qwr '^chr11' ${outputBedGraph}; then
    grep -wr '^chr11' ${outputBedGraph} > ${outputBedGraph}.chr11
  else
    grep -wr '^11' ${outputBedGraph} > ${outputBedGraph}.chr11
  fi

  if grep -qwr '^chr12' ${outputBedGraph}; then
    grep -wr '^chr12' ${outputBedGraph} > ${outputBedGraph}.chr12
  else
    grep -wr '^12' ${outputBedGraph} > ${outputBedGraph}.chr12
  fi

  if grep -qwr '^chr13' ${outputBedGraph}; then
    grep -wr '^chr13' ${outputBedGraph} > ${outputBedGraph}.chr13
  else
    grep -wr '^13' ${outputBedGraph} > ${outputBedGraph}.chr13
  fi

  if grep -qwr '^chr14' ${outputBedGraph}; then
    grep -wr '^chr14' ${outputBedGraph} > ${outputBedGraph}.chr14
  else
    grep -wr '^14' ${outputBedGraph} > ${outputBedGraph}.chr14
  fi

  if grep -qwr '^chr15' ${outputBedGraph}; then
    grep -wr '^chr15' ${outputBedGraph} > ${outputBedGraph}.chr15
  else
    grep -wr '^15' ${outputBedGraph} > ${outputBedGraph}.chr15
  fi

  if grep -qwr '^chr16' ${outputBedGraph}; then
    grep -wr '^chr16' ${outputBedGraph} > ${outputBedGraph}.chr16
  else
    grep -wr '^16' ${outputBedGraph} > ${outputBedGraph}.chr16
  fi

  if grep -qwr '^chr17' ${outputBedGraph}; then
    grep -wr '^chr17' ${outputBedGraph} > ${outputBedGraph}.chr17
  else
    grep -wr '^17' ${outputBedGraph} > ${outputBedGraph}.chr17
  fi

  if grep -qwr '^chr18' ${outputBedGraph}; then
    grep -wr '^chr18' ${outputBedGraph} > ${outputBedGraph}.chr18
  else
    grep -wr '^18' ${outputBedGraph} > ${outputBedGraph}.chr18
  fi

  if grep -qwr '^chr19' ${outputBedGraph}; then
    grep -wr '^chr19' ${outputBedGraph} > ${outputBedGraph}.chr19
  else
    grep -wr '^19' ${outputBedGraph} > ${outputBedGraph}.chr19
  fi

  if grep -qwr '^chr20' ${outputBedGraph}; then
    grep -wr '^chr20' ${outputBedGraph} > ${outputBedGraph}.chr20
  else
    grep -wr '^20' ${outputBedGraph} > ${outputBedGraph}.chr20
  fi

  if grep -qwr '^chr21' ${outputBedGraph}; then
    grep -wr '^chr21' ${outputBedGraph} > ${outputBedGraph}.chr21
  else
    grep -wr '^21' ${outputBedGraph} > ${outputBedGraph}.chr21
  fi

  if grep -qwr '^chr22' ${outputBedGraph}; then
    grep -wr '^chr22' ${outputBedGraph} > ${outputBedGraph}.chr22
  else
    grep -wr '^22' ${outputBedGraph} > ${outputBedGraph}.chr22
  fi

  if grep -qwr '^chrX' ${outputBedGraph}; then
    grep -wr '^chrX' ${outputBedGraph} > ${outputBedGraph}.chrX
  else
    grep -wr '^X' ${outputBedGraph} > ${outputBedGraph}.chrX
  fi

  if grep -qwr '^chrY' ${outputBedGraph}; then
    grep -wr '^chrY' ${outputBedGraph} > ${outputBedGraph}.chrY
  else
    grep -wr '^Y' ${outputBedGraph} > ${outputBedGraph}.chrY
  fi

  for i in 1 {10..19} 2 {20..22} {3..9} X Y;
    do cat ${outputBedGraph}.chr${i} >> ${outputBedGraph}.sorted
  done

  if [ ${referenceBuild} == "GRCh38" ]; then
    ${bedGraphToBigWigTool} ${outputBedGraph}.sorted ${chromSizesGRCh38} ${outputBigWig}
  elif [ ${referenceBuild} == "GRCh37" ]; then
    ${bedGraphToBigWigTool} ${outputBedGraph}.sorted ${chromSizesGRCh37} ${outputBigWig}
  fi

  gzip ${outputBedGraph}.sorted
  rm ${outputBedGraph}.chr*
  rm ${outputBedGraph}
}

main "$@"

