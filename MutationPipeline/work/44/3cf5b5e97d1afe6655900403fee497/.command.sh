#!/bin/bash -ue
mkdir -p B.1.1

cd B.1.1

# Run snpEff to annotate the input VCF file using the provided snpeff_db
snpEff ann -v -c snpEff.config -csvStats annotated_vcf_ann.csv -s annotated_vcf_ann.html NC_045512.2 ../B.1.1.vcf > B.1.1.snpEff.vcf
