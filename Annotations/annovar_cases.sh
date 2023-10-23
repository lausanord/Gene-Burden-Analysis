#!/bin/bash

# Annotate with annovar the VCF of samples

annovarFolder="annovar"
refbuild="hg38"
vcfFile="samples.vcf.gz"
outputPrefix="samples.annotated"

# Ejecutar ANNOVAR
bash annotate.sh ${vcfFile} ${annovarFolder} ${refbuild} ${outputPrefix}

# Combinar genotipos con anotaciones
vcfAnnotated="samples.annotated.GT.vcf.gz"
bcftools annotate -a ${outputPrefix}.vcf.gz -c INFO -Oz -o ${vcfAnnotated} ${vcfFile}
