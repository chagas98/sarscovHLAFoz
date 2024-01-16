#!/usr/bin/env python3

from Bio import SeqIO
from collections import defaultdict
import os
import argparse

def calculate_n_percentage(sequence):
    
    n_count = sequence.count('N')
    total_bases = len(sequence)
    return (n_count / total_bases) * 100 if total_bases > 0 else 0


def mkdir_lineage(lineage):

    path = os.path.join(lineage)
    
    # Check if the subdirectory already exists
    if not os.path.exists(path):
        os.mkdir(path)
        print(f"Subdirectory '{lineage}' created.")

    else:
        print(f"Subdirectory '{lineage}' already exists.")
    
    return path

def find_representatives(input_fasta, percentage):
    lineages = defaultdict(list)
    
    try:
        # Read the input FASTA file and organize sequences by lineage
        for record in SeqIO.parse(input_fasta, "fasta"):
            lineage = record.description.split('|')[1]
            lineages[lineage].append(record)
        
        if not lineages:
            raise ValueError("No sequences found in the input FASTA file.")

        # Write representatives to the output FASTA file
        
        list_of_output=[]
        for lineage, sequences in lineages.items():
            # Find the sequence with the smallest number of 'N' nucleotides
            # and lower than percentage% 'N' nucleotides
            representatives=[]

            for seq in sequences:
                n_percentage = calculate_n_percentage(seq.seq)
                
                if n_percentage < percentage:
                    representatives.append(seq)

            # Check if representatives (more than 2) was found and print its name
            if not representatives or len(representatives) < 2:
                print("No representatives found:", lineage)
            
            else: 
                print("Representatives found:", lineage)
                dir_lineage = mkdir_lineage(lineage)
                fasta_path= os.path.join(dir_lineage, f"{lineage}.fasta")

                # Write the representative to the output file
                SeqIO.write(representatives, fasta_path, "fasta")

            list_of_output.append(f"{lineage}.fasta")


    except FileNotFoundError:
        print(f"Error: The input file '{input_fasta}' was not found.")
    
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


if __name__ == "__main__":
    try:

        parser = argparse.ArgumentParser(description="Separa as sequências por linhagens e exclui sequências com %<perc_n .")
        parser.add_argument("--input_fasta", default="input.fasta", help="Input FASTA file")
        parser.add_argument("--perc_n",  type=float, default=5, help="N percentage")
        args = parser.parse_args()

        input_fasta = args.input_fasta
        N_perc_threshold = args.perc_n

        result = find_representatives(input_fasta, 
                                      N_perc_threshold)

    except KeyboardInterrupt:
        print("\nProcess interrupted by the user.")
