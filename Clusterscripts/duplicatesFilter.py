import sys
import pandas

def df_to_fasta(df):
    fasta_lines = []
    for _, row in df.iterrows():
        uniprot_id = str(row['uniprot_id']).strip()
        sequence = str(row['sequence']).replace(",", "").replace('"', "").strip()
        #p_interface = str(row['p_interface']).replace(",", "").replace('"', "").strip()
        fasta_lines.append(">{}".format(uniprot_id))
        fasta_lines.append(sequence)
        #fasta_lines.append(p_interface)
    return "\n".join(fasta_lines)

def fasta_to_df(file_path):
    ids = []
    sequences = []

    with open(file_path, 'r') as f:
        current_id = None
        current_seq = []

        for line in f:
            line = line.strip()

            if line.startswith(">"):
                if current_id is not None:
                    ids.append(current_id)
                    sequences.append("".join(current_seq))

                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)

        # last entry
        if current_id is not None:
            ids.append(current_id)
            sequences.append("".join(current_seq))

    return pandas.DataFrame({
        "uniprot_id": ids,
        "sequence": sequences
    })

if __name__ == "__main__":
    # Ensure the script is called with the right number of arguments
    if len(sys.argv) != 3:
        print("Usage: python duplicatesFilter.py <parsed_input.csv> <output.fasta>")
        sys.exit(1)

    # Parse the arguments
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Call the filter function
    df = fasta_to_df(input_file)
    df_unique = df.drop_duplicates(subset='sequence')
    fasta_lines = df_to_fasta(df_unique)
    with open(output_file, 'w') as f:
        f.write(fasta_lines)


