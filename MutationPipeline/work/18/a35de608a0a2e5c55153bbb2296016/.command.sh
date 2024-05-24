#!/bin/bash -ue
mkdir -p AY.122

cd AY.122

# Run snpEff to annotate the input VCF file using the provided snpeff_db
snpEff ann -v -c snpEff.config -csvStats annotated_vcf_ann.csv -s annotated_vcf_ann.html NC_045512.2 ../AY.122.vcf > AY.122.snpEff.vcf
