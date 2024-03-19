#!/bin/bash -ue
mkdir -p P.2

cd P.2

# Run snpEff to annotate the input VCF file using the provided snpeff_db
snpEff ann -v -c snpEff.config -csvStats annotated_vcf_ann.csv -s annotated_vcf_ann.html NC_045512.2 ../P.2.vcf > P.2.snpEff.vcf
