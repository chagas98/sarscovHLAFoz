#!/bin/bash -ue
variant=AY.101

mkdir -p ${variant}

cd ${variant}

#mkdir -p prediction

#cd prediction

export PYTHONWARNINGS="ignore::RuntimeWarning"
export PYTHONWARNINGS="ignore::FutureWarning"

# Run epaa_mod.py to predict epitope for the input VCF file    
epaa_mod.py --identifier ${variant} --alleles "B*07:02;B*44:03" --tools "netmhcpan-4.1" --max_length 12 --min_length 10 --versions ../versions.csv --variant_lineage ${variant} --somatic_mutation ../AY.101.snpEff.vcf
