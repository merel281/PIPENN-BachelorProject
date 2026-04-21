# Match clustered fasta file back to BioLiP.*.txt
# Output clustered BioLip.*.txt file
from Bio import SeqIO

# Identify PDB IDs and chain from clustered fasta file 
# Store them in a set
with open("clusteredNucleotides.fasta", 'r') as file: # 'r' open the file in read mode  
    peptide_set = set()

    for rec in SeqIO.parse("clusteredNucleotides.fasta", "fasta"):
        pdb_id = rec.id[:4].lower()
        chain = rec.id[4:] if len(rec.id) > 4 else ""

        peptide_set.add((pdb_id, chain)) # e.g., peptide_set = {('1a1a', 'B'), ('2abc', 'A')}

merged_data = {}

# Add this line for nucleotides, because combination of DNA and RNA
for nucleotides in [ "BioLiP_DNA.txt", "BioLiP_RNA.txt"]:
    
    with open(nucleotides, 'r') as file_in:
            for line in file_in:
                parts = line.strip().split("\t")
                if len(parts) < 2:
                    continue
    
                pdb_id2 = parts[0].lower() # Column 1: PDB ID
                chain2 = parts[1] # Column 2: chain
                key = (pdb_id2, chain2)
    
                # Only consider entries that are in your peptide_set
                if key not in peptide_set:
                    continue
                
                if key not in merged_data:
                    merged_data[key] = parts
                else:
                    # Append parts[8] to existing entry (avoid duplicates)
                    existing = merged_data[key][8]
                    new_value = parts[8]
    
                    if new_value not in existing.split(","):
                        merged_data[key][8] += " " + new_value
with open("ClusBioNucleotides.txt","w") as test_out:
    for parts in merged_data.values():
        test_out.write("\t".join(parts) + "\n") # write full line to new file




            