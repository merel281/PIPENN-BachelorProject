from Bio import SeqIO

# Read BioLiP peptide PDB-chain pairs
with open("BioLiP_peptide.txt") as f1:
    peptide_set = set(
        (parts[0].lower(), parts[1])
        for line in f1 if (parts := line.split()) and len(parts) >= 2
    )

# Read BioLiP metal PDB-chain pairs
with open("BioLiP_metal.txt") as f2:
    metal_set = set(
        (parts[0].lower(), parts[1])
        for line in f2 if (parts := line.split()) and len(parts) >= 2
    )
# Read BioLiP DNA and RNA PDB-chain pairs
with open("BioLiP_DNA.txt") as f3, open("BioLiP_RNA.txt") as f4:
         nucleotide_set = set(
             (parts[0].lower(), parts[1])
             for f in (f3, f4)
             for line in f
             if (parts := line.split()) and len(parts) >= 2
         )


# Initialize counters
peptide_count = 0
metal_count = 0
nucleotide_count = 0

peptide_not_found = set(peptide_set)
metal_not_found = set(metal_set)
nucleotide_not_found = set(nucleotide_set)



# Open output FASTA files
with open("peptide.fasta", "w") as peptide_out, \
     open("metal.fasta", "w") as metal_out, \
     open("nucleotide.fasta", "w") as nucleotide_out:

    # Stream through protein.fasta
    for rec in SeqIO.parse("protein.fasta", "fasta"):
        # EXACT jouw manier van parsen behouden
        rec_pdb_id = rec.id[:4].lower()
        rec_chain = rec.id[4:] if len(rec.id) > 4 else ""

        key = (rec_pdb_id, rec_chain)

        # Check peptide
        if key in peptide_set:
            SeqIO.write(rec, peptide_out, "fasta")
            peptide_count += 1
            peptide_not_found.discard(key)

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