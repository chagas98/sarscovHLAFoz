#!/bin/bash -ue
mkdir -p P.7

cd P.7

# Run snpEff to annotate the input VCF file using the provided snpeff_db
snpEff ann -v -c snpEff.config -csvStats annotated_vcf_ann.csv -s annotated_vcf_ann.html NC_045512.2 ../P.7.vcf > P.7.snpEff.vcf
