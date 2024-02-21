#!/bin/bash -ue
mkdir -p AY.101

cd AY.101

# Run snpEff to annotate the input VCF file using the provided snpeff_db
snpEff ann -v -c snpEff.config -csvStats annotated_vcf_ann.csv -s annotated_vcf_ann.html NC_045512.2 ../AY.101.vcf > AY.101.snpEff.vcf
