#!/bin/bash -ue
mkdir -p AY.99.2

cd AY.99.2

# Run snpEff to annotate the input VCF file using the provided snpeff_db
snpEff ann -v -c snpEff.config -csvStats annotated_vcf_ann.csv -s annotated_vcf_ann.html NC_045512.2 ../AY.99.2.vcf > AY.99.2.snpEff.vcf
