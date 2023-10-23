set -eu
test() {
  controlCount=./example/1KG/concatenated.gnomad.VEP.fix2.vcf.gz.gds
  caseCount=./example/1KG/samples.annotated.biallelic.leftnorm.ABCheck.noCHR.VEP.fix2.vcf.gz.gds
  sampleList=./example/1KG/samples.txt
  intersectBed=./example/1KG/intersect.coverage10x.bed.gz
  variantExclude=./example/1KG/gnomAD.exclude.allow.segdup.lcr.v3.txt.gz
  AFMax=5e-4
  maxAFPopmax=1
  variantMissing=1
  ACANConfig=./example/1KG/stratified_config_gnomad.txt
  ancestryFile=./example/1KG/ethnicity.txt
  variantGroup="annovar_pathogenic" 
  REVELThreshold=0.65
  pLDControl=0.05
  highLDVariantFile=./example/1KG/full_vs_gnomAD.p0.05.OR1.ignoreEthnicityInLD.rds
  outputPrefix=./example/1KG/result/AFMax${AFMax}.${variantGroup}.${REVELThreshold}.excludeV3.LDv2


  mkdir -p ./example/1KG/result/

  Rscript utilities/CoCoRV_wrapper.R \
    --sampleList ${sampleList} \
    --outputPrefix ${outputPrefix} \
    --AFMax ${AFMax} \
    --maxAFPopmax ${maxAFPopmax} \
    --bed ${intersectBed} \
    --variantMissing ${variantMissing} \
    --variantGroup ${variantGroup} \
    --removeStar \
    --ACANConfig ${ACANConfig} \
    --minREVEL ${REVELThreshold} \
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
