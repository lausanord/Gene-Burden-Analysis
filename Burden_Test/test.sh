
set -eu
test() {
  controlCount=./example/1KG/concatenated.gnomad.VEP.fix2.vcf.gz.gds
  caseCount=./example/1KG/samples.annotated.biallelic.leftnorm.ABCheck.noCHR.VEP.fix2.vcf.gz.gds
  sampleList=./example/1KG/samples.txt
  intersectBed=./example/1KG/intersect.coverage10x.bed.gz
  variantExclude=./example/1KG/gnomAD.exclude.allow.segdup.lcr.v3.txt.gz
  AFMax=1e-4
  maxAFPopmax=1
  variantMissing=1
  ACANConfig=./example/1KG/stratified_config_gnomad.txt
  ancestryFile=./example/1KG/ethnicity.txt
  variantGroup="LOF_Alphamissense"
  pLDControl=0.05
  REVELThreshold=0.65
  highLDVariantFile=./example/1KG/full_vs_gnomAD.p0.05.OR1.ignoreEthnicityInLD.rds
  outputPrefix=./example/1KG/result/AFMax${AFMax}.${variantGroup}.${REVELThreshold}.excludeV3_Control.NFE.ARREGLADO.LDv2
  annotationList=./example/1KG/annotationList.txt
  variantGroupCustom=./example/1KG/LOF_missenseCADD.R
  mkdir -p ./example/1KG/result

  Rscript utilities/CoCoRV_wrapper.R \
    --sampleList ${sampleList} \
    --outputPrefix ${outputPrefix} \
    --variantGroupCustom ${variantGroupCustom} \
    --AFMax ${AFMax} \
    --annotationList ${annotationList} \
    --maxAFPopmax ${maxAFPopmax} \
    --AFUse "control" \
    --minREVEL ${REVELThreshold} \
    --bed ${intersectBed} \
    --variantMissing ${variantMissing} \
    --variantGroup ${variantGroup} \
    --ACANConfig ${ACANConfig} \
    --variantExcludeFile ${variantExclude}  \
    --checkHighLDInControl \
    --pLDControl ${pLDControl} \
    --highLDVariantFile ${highLDVariantFile} \
    ${controlCount} \
    ${caseCount}

  outputFile=${outputPrefix}.association.tsv

  # plot the dominant model with lambda estimation
  Rscript utilities/QQPlotAndFDR.R ${outputFile} \
             ${outputFile}.dominant.nRep1000 --setID gene \
             --outColumns gene --n 1000 \
        --pattern "case.*Mutation.*_DOM$|control.*Mutation.*_DOM$" --FDR

}

main() {
  test
}

main "$@"
