#!/bin/bash -ue
mkdir -p BA.2

cd BA.2

# Run snpEff to annotate the input VCF file using the provided snpeff_db
snpEff ann -v -c snpEff.config -csvStats annotated_vcf_ann.csv -s annotated_vcf_ann.html NC_045512.2 ../BA.2.vcf > BA.2.snpEff.vcf
