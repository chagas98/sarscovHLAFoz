#!/bin/bash -ue
mkdir -p BE.1

cd BE.1

# Run snpEff to annotate the input VCF file using the provided snpeff_db
snpEff ann -v -c snpEff.config -csvStats annotated_vcf_ann.csv -s annotated_vcf_ann.html NC_045512.2 ../BE.1.vcf > BE.1.snpEff.vcf
