#!/bin/bash

# Annotate with annovar the VCF of controls

annovarFolder="annovar"
refbuild="hg19"
gnomADVCF="gnomad.exomes.r2.1.1.sites.vcf.gz"
outputPrefix="gnomAD.annotated"

# Ejecutar ANNOVAR
bash annotate.sh ${gnomADVCF} ${annovarFolder} ${refbuild} ${outputPrefix}
