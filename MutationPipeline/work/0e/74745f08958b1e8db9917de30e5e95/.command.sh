#!/bin/bash -ue
# Filter the input fasta file based on quality and number of sequences

filtering.py --input_fasta GISAID_sequences.fasta --perc_n 5 --minimum 1
