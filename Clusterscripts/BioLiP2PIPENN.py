# CLustered BioLiP.*.txt file convert to PIPENN format
# PDB ID-chain, sequence, Rlength, normalized_length, n_interface
# Additional: Uniprot_id, catalytic site residues, Go-terms

# Credits to Bart Vonk for inspo


# argv = argument vector
# A list of strings that contains everything you typed after python script.py in the terminal
# What to run in windows cmd:
# python script.py input.fasta output.csv
from sys import argv
import os 

def read_input(ClusteredFile):

    # Store unique sequences
    sequence_data = {}

    minRlength = 100
    maxRlength = 0

    with open(ClusteredFile, 'r') as file:

        for line in file:

            parts = line.strip().split("\t")

            if len(parts) < 21:
                continue

            header = parts[0] + '-' + parts[1]

            raw_sequence = parts[20]
            sequence = '"' + ','.join(raw_sequence) + '"'

            length = len(raw_sequence)

            # Track min/max length
            if length < minRlength:
                minRlength = length

            if length > maxRlength:
                maxRlength = length

            # Binding sites
            binding = parts[8].split()

            # Sequence already exists
            if raw_sequence in sequence_data:

                # Existing binding sites
                existing_binding = sequence_data[raw_sequence]["binding"]

                # Append only new residues; prevent duplicate bindingsites in list
                for site in binding:
                    if site not in existing_binding:
                        existing_binding.append(site)

            else:
                # New sequence
                sequence_data[raw_sequence] = {
                    "header": header,
                    "sequence": sequence,
                    "length": length,
                    "binding": binding
                }

    # Convert dictionary back to lists
    headers = ["uniprot_id"]
    sequences = ["sequence"]
    lengths = ["Rlength"]
    bindingsites = ["Bindingsites"]

    for data in sequence_data.values():

        headers.append(data["header"])
        sequences.append(data["sequence"])
        lengths.append(data["length"])
        bindingsites.append(data["binding"])

    return headers, sequences, lengths, bindingsites, minRlength, maxRlength

# Function calculate normalized length
def calc_norm(minRlength, maxRlength, lengths):
    # this function calculates the normalized length for every sequence length in the list 'lengths'
    normalized = ["normalized_length"]
    for seqLength in lengths[1:]:
        # list starts one index later due to the fact the first index is the column header
        aaNormList = []
        for i in range(seqLength):
            #aaNormList.append((seqLength - minRlength) / float(maxRlength - minRlength))
            aaNormList.append((seqLength - minRlength) / float(2050 - minRlength))
        thing = ','.join([str(x) for x in aaNormList])
        use = '"' + thing + '"'
        normalized.append(use)
    return normalized



# Prediction filled in with zeros; length equal to sequence length
def interface(lengths, bindingsites):
    interface_0 = ["n_interface"]
    
    for seq_idx, seqLength in enumerate(lengths):
        try:
            seqLength = int(seqLength) # Turn the value into a number, fails for Rlength, so excludes header
        except ValueError:
            continue
            
        prediction = [0] * seqLength # Initialize list of zeros equal to sequence length; output without bindingsites

        for site in bindingsites[seq_idx]:
            try:
                pos = int(site[1:])-1 # Only position without aminoacid

                if 0<= pos <seqLength: # Makes sure the position fits inside the sequence
                    prediction[pos] = 1 # replace 0 with 1
                else:
                    print(f"Warning: position {pos+1} out of range for sequence length {seqLength}")

            except ValueError:
                print(f"Invalid binding site format: {site}")

        # Convert list to string
        thing = ','.join([str(x) for x in prediction])
        use = f'"{thing}"'
        
        interface_0.append(use)
    return interface_0

    
def write_csv(OutData, csvFile):
    # function for rewriting FASTA data into csv format
    # all parsed information from the input file is found in 'fastaData', while 'csvFile' contains the name of the new csv file
    with open(csvFile, 'w') as file:
        for i in range(len(OutData[0])):
            row = [str(OutData[j][i]) for j in range(len(OutData))]
            file.write(','.join(row) + '\n')

def main():
    # If no input and/or output files are mentioned in command line
    if len(argv) < 2:
        print("Usage: python script.py ClusBio*.txt")
        exit()    
    # Sets variables from command line, input and output files
    ClusteredFile = argv[1]
    csvFile = os.path.splitext(ClusteredFile)[0] + "2.csv"
    # collects all data from the input file
    print("collecting data from input file")

    headers, sequences, lengths, bindingsites, minRlength, maxRlength = read_input(ClusteredFile)
    # Only keep list columns for write_csv
    OutData = [headers, sequences, lengths]

    # adds column of normalized lengths
    OutData.append(calc_norm(minRlength, maxRlength, OutData[2]))

    # adds column for now only zeros equal to the sequence length
    OutData.append(interface(OutData[2], bindingsites))
    
    # writes the collected and formatted data to csv
    print("writing output data in PIPENN format")
    write_csv(OutData, csvFile)


# Run code
main()


