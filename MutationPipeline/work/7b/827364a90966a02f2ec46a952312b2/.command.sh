#!/bin/bash -ue
mkdir -p BQ.1.1

cd BQ.1.1

# Run snpEff to annotate the input VCF file using the provided snpeff_db
snpEff ann -v -c snpEff.config -csvStats annotated_vcf_ann.csv -s annotated_vcf_ann.html NC_045512.2 ../BQ.1.1.vcf > BQ.1.1.snpEff.vcf
