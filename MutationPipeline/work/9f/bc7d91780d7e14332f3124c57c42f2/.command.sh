#!/bin/bash -ue
mkdir -p BA.5.2.1

cd BA.5.2.1

# Run snpEff to annotate the input VCF file using the provided snpeff_db
snpEff ann -v -c snpEff.config -csvStats annotated_vcf_ann.csv -s annotated_vcf_ann.html NC_045512.2 ../BA.5.2.1.vcf > BA.5.2.1.snpEff.vcf
