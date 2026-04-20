#!/usr/bin/env python3

import sys

def clean_fasta(input_file, output_file):
    """
    Cleans a non-standard FASTA file:
    - Keeps only first token of header (sequence ID)
    - Joins multi-line sequences
    - Removes empty lines and whitespace
    """

    with open(input_file, "r") as fin, open(output_file, "w") as fout:
        seq_id = None
        seq_parts = []

        for line in fin:
            line = line.strip()

            if not line:
                continue  # skip empty lines

            if line.startswith(">"):
                # write previous record if exists
                if seq_id is not None:
                    fout.write(f">{seq_id}\n")
                    fout.write("".join(seq_parts) + "\n")

                # parse header: take only first token after ">"
                header_fields = line[1:].split()
                seq_id = header_fields[0]
                seq_parts = []

            else:
                # sequence line
                seq_parts.append(line)

        # write last record
        if seq_id is not None:
            fout.write(f">{seq_id}\n")
            fout.write("".join(seq_parts) + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: clean_fasta.py <input.fasta> <output.fasta>")
        sys.exit(1)

    clean_fasta(sys.argv[1], sys.argv[2])
