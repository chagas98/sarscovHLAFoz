#!/bin/bash -ue
mkdir -p AY.122

cd AY.122

# Run bcftools to perform variant calling on the input file using the reference file

bcftools mpileup --max-depth 3000 -f ../NC_045512.2.fasta ../AY.122.bam | bcftools call -c --ploidy 1 -Ov -o AY.122.vcf
