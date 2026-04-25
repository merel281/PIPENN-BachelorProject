# Match clustered fasta file back to BioLiP.*.txt
# Output clustered BioLip.*.txt file
from Bio import SeqIO
from sys import argv 

# Identify PDB IDs and chain from clustered fasta file 
# Store them in a set
with open(argv[1], 'r') as file: # 'r' open the file in read mode  
    peptide_set = set()

    for rec in SeqIO.parse(argv[1], "fasta"):
        pdb_id = rec.id[:4].lower()
        chain = rec.id[4:] if len(rec.id) > 4 else ""

        peptide_set.add((pdb_id, chain)) # e.g., peptide_set = {('1a1a', 'B'), ('2abc', 'A')}

merged_data = {}

# BioLiP file to match back to 
#biolip = argv[2] 
# For nucleotides use:
biolip_files = [argv[2], "../input/BioLiP_RNA.txt"]
for biolip in biolip_files: # Whatch out for correct indent
    
  with open(biolip, 'r') as file_in:
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
                      
out = argv[1].rsplit('.',1)[0] + ".csv"

with open(out,"w") as outbio:
    for parts in merged_data.values():
        outbio.write("\t".join(parts) + "\n") # write full line to new file




            