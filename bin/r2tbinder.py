#!/usr/bin/env python

import glob
from collections import defaultdict

def read_fasta(file_path):
    """Reads a FASTA file and returns a dictionary of header and sequence pairs."""
    with open(file_path, 'r') as file:
        sequence_data = defaultdict(list)  # Use a list to handle multiple sequences with the same identifier
        current_header = ""
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                current_header = line  # Store the full header
                sequence_data[current_header] = []  # Start a new sequence list for this header
            else:
                sequence_data[current_header].append(line)  # Append sequence lines to the current header
        return sequence_data

def bin_fasta_files(fasta_files_pattern):
    """Bins sequences from multiple FASTA files based on their common identifier in the headers."""
    all_sequences = defaultdict(list)
    
    # Read all fasta files and accumulate sequences by header's common identifier
    for fasta_file in glob.glob(fasta_files_pattern):
        sequences = read_fasta(fasta_file)
        for header, sequence_lines in sequences.items():
            identifier = header.split('|')[0].strip()  # Extract the common identifier
            all_sequences[identifier].append((header, sequence_lines))
    
    # Write binned sequences to separate files, each sequence retains its original header
    for identifier, headers_and_sequences in all_sequences.items():
        safe_identifier = identifier.replace(">", "").replace("|", "_").replace(" ", "_")[:20]  # Create a file-safe identifier
        file_name = f"{safe_identifier}.fasta"
        with open(file_name, 'w') as file:
            for header, seq_lines in headers_and_sequences:
                file.write(f"{header}\n")
                file.write('\n'.join(seq_lines) + '\n')
# Example usage
bin_fasta_files("r2t_markers/*.fa")
