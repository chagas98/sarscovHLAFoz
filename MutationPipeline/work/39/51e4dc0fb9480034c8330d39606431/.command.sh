#!/bin/bash -ue
mkdir -p BA.5.3

cd BA.5.3

# Run snpEff to annotate the input VCF file using the provided snpeff_db
snpEff ann -v -c snpEff.config -csvStats annotated_vcf_ann.csv -s annotated_vcf_ann.html NC_045512.2 ../BA.5.3.vcf > BA.5.3.snpEff.vcf
