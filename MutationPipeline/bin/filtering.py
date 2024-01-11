#!/usr/bin/env python3

from Bio import SeqIO
from collections import defaultdict
import os
import argparse

def calculate_n_percentage(sequence):
    n_count = sequence.count('N')
    total_bases = len(sequence)
    return (n_count / total_bases) * 100 if total_bases > 0 else 0


def mkdir_lineage(lineage, main_directory):

    path = os.path.join(main_directory, lineage)
    
    # Check if the subdirectory already exists
    if not os.path.exists(path):
        os.mkdir(path)
        print(f"Subdirectory '{lineage}' created.")

    else:
        print(f"Subdirectory '{lineage}' already exists.")
    
    return path

def find_representatives(input_fasta, percentage, variants_dir):
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
            if not representatives or len(representatives) < 1:
                print("No representatives found:", lineage)
            
            else: 
                print("Representatives found:", lineage)
                dir_lineage = mkdir_lineage(lineage, variants_dir)
                fasta_path= os.path.join(dir_lineage, lineage)

                # Write the representative to the output file
                SeqIO.write(representatives, fasta_path, "fasta")

            list_of_output.append(fasta_path)

        return list_of_output 

    except FileNotFoundError:
        print(f"Error: The input file '{input_fasta}' was not found.")
    
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


if __name__ == "__main__":
    try:

        parser = argparse.ArgumentParser(description="Description of your script.")
        parser.add_argument("--input_fasta", default="test.fasta", help="Input FASTA file")
        args = parser.parse_args()

        input_fasta = args.input_fasta
        output_fasta = "representatives.fasta"
        variants_dir = "variants"
        N_perc_threshold = 5

        if not os.path.exists(variants_dir):
            os.mkdir(variants_dir)
            print(f"Directory '{variants_dir}' created.")

        else:
            print(f"Directory '{variants_dir}' already exists.")
    

        result = find_representatives(input_fasta, 
                                      N_perc_threshold, 
                                      variants_dir)

        print(result)

    except KeyboardInterrupt:
        print("\nProcess interrupted by the user.")
