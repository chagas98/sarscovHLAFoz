#!/bin/bash -ue
mkdir -p BE.10

cd BE.10

# Run snpEff to annotate the input VCF file using the provided snpeff_db
snpEff ann -v -c snpEff.config -csvStats annotated_vcf_ann.csv -s annotated_vcf_ann.html NC_045512.2 ../BE.10.vcf > BE.10.snpEff.vcf
