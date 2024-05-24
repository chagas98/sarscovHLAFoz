#!/bin/bash -ue
# Download the fasta file from NCBI using the provided sequence ID

wget -O NC_045512.2.fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?db=nuccore&report=fasta&id=NC_045512.2"
