# Defide the protein fasta database into three subsets
# Peptide
# Nucleotides (DNA + RNA)
# Small molecules (Metal)

#Import BioPython
from Bio import SeqIO



#Identify PDB_IDs from BioLiP file (first column) and chain (second column) 
# and stor them in a set (less data storage than a list)
# PDB ID + chain peptide
with open("BioLiP_peptide.txt") as f1:
    peptide_set = set(
        (parts[0].lower(), parts[1])
        for line in f1 if (parts := line.split()) and len(parts) >= 2
    )

# PDB ID + chain metal
with open("BioLiP_metal.txt") as f2:
    metal_set = set(
        (parts[0].lower(), parts[1])
        for line in f2 if (parts := line.split()) and len(parts) >= 2
    )
# PDB ID + chain nucleotiden (DNA + RNA)
with open("BioLiP_DNA.txt") as f3, open("BioLiP_RNA.txt") as f4:
         nucleotide_set = set(
             (parts[0].lower(), parts[1])
             for f in (f3, f4)
             for line in f
             if (parts := line.split()) and len(parts) >= 2
         )
# This way does not include duplicates from BioLiP.*.txt files
# Less entries in *_set than in BioLiP.*.txt files 



# Counters
peptide_count = 0
metal_count = 0
nucleotide_count = 0

peptide_not_found = set(peptide_set) # Start assuming all are missing
metal_not_found = set(metal_set)
nucleotide_not_found = set(nucleotide_set)



# Open output FASTA files
# "w" removes all content from the fill, start with an empty file
with open("peptide.fasta", "w") as peptide_out, \
     open("metal.fasta", "w") as metal_out, \
     open("nucleotide.fasta", "w") as nucleotide_out:

    # Stream through protein.fasta one sequence at a time
    for rec in SeqIO.parse("protein.fasta", "fasta"):
        # Extract PDB ID from FASTA header
        #rec. id = the FASTA header up to first space
        rec_pdb_id = rec.id[:4].lower() # PDB IDs are always 4 characters and in lowercase in this database; avoids false positives from substring matching
        rec_chain = rec.id[4:] if len(rec.id) > 4 else "" # Chain of protein works for A, B, C but also A0, AA, BL, B10 etc

        # PDB ID + chain 
        key = (rec_pdb_id, rec_chain)

        # Check peptide
        if key in peptide_set:
            SeqIO.write(rec, peptide_out, "fasta") # Add fasta of found PDB IDs to output
            peptide_count += 1
            peptide_not_found.discard(key) # remove found IDs from missing set

        # Check metal
        if key in metal_set:
            SeqIO.write(rec, metal_out, "fasta")
            metal_count += 1
            metal_not_found.discard(key)

        # Check nucleotides
        if key in nucleotide_set:
            SeqIO.write(rec, nucleotide_out, "fasta")
            nucleotide_count += 1
            nucleotide_not_found.discard(key)



# Print stats
print("=== PEPTIDE ===")
print(f"Total PDB IDs in input: {len(peptide_set)}")
print(f"Total entries saved: {peptide_count}")
print(f"Total entries not found: {len(peptide_not_found)}")

print("\n=== METAL ===")
print(f"Total PDB IDs in input: {len(metal_set)}")
print(f"Total entries saved: {metal_count}")
print(f"Total entries not found: {len(metal_not_found)}")

print("\n=== NUCLEOTIDE ===")
print(f"Total PDB IDs in input: {len(nucleotide_set)}")
print(f"Total entries saved: {nucleotide_count}")
print(f"Total entries not found: {len(nucleotide_not_found)}")

# Print missing entries
print("\nMissing peptide entries:")
for pdb, chain in sorted(peptide_not_found):
    print(pdb, chain)

print("\nMissing metal entries:")
for pdb, chain in sorted(metal_not_found):
    print(pdb, chain)

print("\nMissing nucleotide entries:")
for pdb, chain in sorted(nucleotide_not_found):
    print(pdb, chain)