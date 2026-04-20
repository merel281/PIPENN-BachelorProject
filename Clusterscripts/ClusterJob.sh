#!/bin/bash
# CLusters fasta file and removes ZK-488 proteins from clustered database
# input file is in fasta format
# line 81-102 only if ZK-488 proteins are removed from database
# line 84 make sure filename corresponds with input file

# e -> exit immediately if any command fails
# u -> treat unset variables as errors
# 0 pipefail -> if any command in a pipeline fails, the whole pipeline fails
set -euo pipefail

# Argument check; did the user provide exactly 1 argument; The input file consisting a fasta file that needs to be clustered.
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

# absolute path of the script's folder
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Where the user executed the script from
CALL_DIR="$PWD"

# Creates data/ folder, -p prevents errors if folder already exists
# Files from in between steps will be stored here + output file
mkdir -p data
# Create output file
mkdir -p output
# Creates temporary working directory inside call location
mkdir -p "$CALL_DIR/data/tmp"


echo "Cleaning input format"
# runs python script clean_fasta.py; $1 = first argument (input fasta file); output goes to data/RAW.fasta
# clean_fasta.py cleans a non standard fasta file (removes extra info from headline)
python "$SCRIPT_DIR/clean_fasta.py" "$1" "$CALL_DIR/data/RAW.fasta"

echo "Filtering sequences"
# Input: cleaned fasta; output: filtered fasta; sequences between 26 and 2000 amino acids
python "$SCRIPT_DIR/generalFilter.py" \
    "$CALL_DIR/data/RAW.fasta" \
    "$CALL_DIR/data/filtered.fasta" \
    26 2000

echo "Removing duplicates"
# remove identical sequences; input: filtered fasta; output: deduplicated.fasta
python "$SCRIPT_DIR/duplicatesFilter.py" \
    "$CALL_DIR/data/filtered.fasta" \
    "$CALL_DIR/data/deduplicated.fasta"

echo 'Clustering data'
# filtered fasta to clustered fasta
# Copy file; perserve the original file and create a working file for the next step
cp "$CALL_DIR/data/deduplicated.fasta" "$CALL_DIR/data/precluster.fasta"

# Create clean working directory for MMSeqs
# Deletes the entire tmp directory and everything inside it
rm -rf "$CALL_DIR/data/tmp" # rm = remove; -r=  recursive (delete folder + all contents); -f = force (no prompts, ignore missing files) 
# Create a fresh, empty tmp directory
mkdir -p "$CALL_DIR/data/tmp" # mkdir = make directory; -p= create parent dirs if needed, don't error if it exists

# Cluster data with the use of MMSeqs
# Create a subshell with (...) -> module load stay local
(
    # load software modules
    module load 2021
    module load netsurfp
    
    export MMSEQS_CALL_DEPTH=0

    # CLustering of protein sequences; sequences must be >25% identical to cluster; alignment must cover query sequence; 90% coverage required
    # Sequences grouped into clusters; similar proteins grouped together
    mmseqs easy-linclust \
        "$CALL_DIR/data/precluster.fasta" \
        "$CALL_DIR/data/clusterRes" \
        "$CALL_DIR/data/tmp" \
        --min-seq-id 0.25 --cov-mode 1 -c 0.9
    
    echo "MMseqs output files 1.0"
    ls -lh "$CALL_DIR/data" | grep clusterRes
        
    #!This part only if ZK-448 proteins are removed (nucleotides and metals)
    echo "removing mapped clusters"
    # extract 2nd column of ZK-488_*.csv and skip the header
    tail -n +2 ZK-488_metals.csv | cut -d',' -f2 > pdb_ZK488.txt
    python "$SCRIPT_DIR/eliminateClusters.py" \
     "$CALL_DIR/data/clusterRes_all_seqs.fasta" \
     pdb_ZK488.txt\
     "$CALL_DIR/data/filtered_mapped.fasta" 

    echo 'Fetching representatives'
    # Re-cluster the cleaned dataset; Produce representative clusters quickly without making an unnecessary big dataset; 
    # Output clusterRes2_* = final clustered dataset
    mmseqs easy-linclust \
        "$CALL_DIR/data/filtered_mapped.fasta" \
        "$CALL_DIR/data/clusterRes2" \
        "$CALL_DIR/data/tmp" \
        --min-seq-id 0.25 --cov-mode 1 -c 0.9
        
    echo "MMseqs output files 2.0"
    ls -lh "$CALL_DIR/data" | grep clusterRes
        
    #! End of ZK-488 protein removal

    # Avoids conflicts later
    module unload netsurfp
)

echo "Finalizing dataset"

ls -lh "$CALL_DIR/data"| grep clusterRes2

# All sequences grouped into final clusters
cp "$CALL_DIR/data/clusterRes2_all_seqs.fasta" "$CALL_DIR/output/finalclus3.0.fasta"
# Contains one representative per cluster -> non-redundant dataset
cp "$CALL_DIR/data/clusterRes2_rep_seq.fasta" "$CALL_DIR/output/final3.0.fasta"
