#!/bin/bash -ue
mkdir -p AY.46.3

cd AY.46.3

# Run snpEff to annotate the input VCF file using the provided snpeff_db
snpEff ann -v -c snpEff.config -csvStats annotated_vcf_ann.csv -s annotated_vcf_ann.html NC_045512.2 ../AY.46.3.vcf > AY.46.3.snpEff.vcf
