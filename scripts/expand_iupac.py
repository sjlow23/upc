#!/usr/bin/env python 

import sys
from itertools import product
from Bio.Seq import Seq
from Bio.SeqUtils import IUPACData

# Function to expand an IUPAC code to its nucleotides
def expand_iupac(code):
	return IUPACData.ambiguous_dna_values[code]

# Function to generate all possible combinations of expanded sequences (primers+probes)
def generate_combinations_all(forward, reverse, probe):
	expanded_forward = [''.join(seq) for seq in product(*map(expand_iupac, forward))]
	expanded_reverse = [''.join(seq) for seq in product(*map(expand_iupac, reverse))]
	expanded_probe = [''.join(seq) for seq in product(*map(expand_iupac, probe))]
	return product(expanded_forward, expanded_reverse, expanded_probe)

# Function to generate all possible combinations of expanded primers
def generate_combinations_primers(forward, reverse):
	expanded_forward = [''.join(seq) for seq in product(*map(expand_iupac, forward))]
	expanded_reverse = [''.join(seq) for seq in product(*map(expand_iupac, reverse))]
	return product(expanded_forward, expanded_reverse)

# Function to generate all possible combinations of expanded probes
def generate_combinations_probes(probe):
	expanded_probe = [''.join(seq) for seq in product(*map(expand_iupac, probe))]
	return product(expanded_probe)

# Check if the correct number of arguments is provided
if len(sys.argv) != 4:
	print("Usage: python {sys.argv[0]} input_file output_file")
	sys.exit(1)


input_file = sys.argv[1]
output_file = sys.argv[2]
mode = sys.argv[3]

# Read input file and process each line
with open(input_file, 'r') as file:
	lines = file.readlines()

# Write output to a new file
if mode == "both":
	with open(output_file, 'w') as file:
		for line in lines:
			name, forward, reverse, probe = line.strip().split('\t')
			combinations = generate_combinations_all(forward, reverse, probe)
			name_counter = 1
			for comb in combinations:
				expanded_name = f"{name}_{name_counter}"
				file.write(f"{name}\t{expanded_name}\t{comb[0]}\t{comb[1]}\t{comb[2]}\n")
				name_counter += 1

# Write output to a new file- primers only
if mode == "primers":
	with open(output_file, 'w') as file:
		for line in lines:
			name, forward, reverse = line.strip().split('\t')
			combinations = generate_combinations_primers(forward, reverse)
			name_counter = 1
			for comb in combinations:
				expanded_name = f"{name}_{name_counter}"
				file.write(f"{name}\t{expanded_name}\t{comb[0]}\t{comb[1]}\n")
				name_counter += 1


# Write output to a new file- probes only
if mode == "probes":
	with open(output_file, 'w') as file:
		for line in lines:
			name, probe = line.strip().split('\t')
			combinations = generate_combinations_probes(probe)
			name_counter = 1
			for comb in combinations:
				expanded_name = f"{name}_probe{name_counter}"
				file.write(f"{name}\t{expanded_name}\t{comb[0]}\n")
				name_counter += 1