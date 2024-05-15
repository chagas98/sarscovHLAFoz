#!/bin/bash -ue
mkdir -p DL.1

cd DL.1

# Run snpEff to annotate the input VCF file using the provided snpeff_db
snpEff ann -v -c snpEff.config -csvStats annotated_vcf_ann.csv -s annotated_vcf_ann.html NC_045512.2 ../DL.1.vcf > DL.1.snpEff.vcf
