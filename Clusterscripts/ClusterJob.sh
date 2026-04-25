#!/bin/bash
# CLusters fasta file and removes ZK-488 proteins from clustered database
# input file is in fasta format
# line 90-111 only if ZK-488 proteins are removed from database
# line 93 make sure filename corresponds with input file


# e -> exit immediately if any command fails
# u -> treat unset variables as errors
# 0 pipefail -> if any command in a pipeline fails, the whole pipeline fails
set -euo pipefail

# Instal bio python
which python
python -V
python -c "import Bio; print(Bio.__file__)"

# Argument check; did the user provide exactly 1 argument; The input file consisting a fasta file that needs to be clustered.
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_file> <BioLiP_*.txt>"
    exit 1
fi

# absolute path of the script's folder
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"

# Create files and directories
# -p prevents errors if folder already exists
# Files from in between steps will be stored here + output file
DATA_DIR="$BASE_DIR/data"
IN_DIR="$BASE_DIR/input"
OUT_DIR="$BASE_DIR/output"
TMP_DIR="$DATA_DIR/tmp"

mkdir -p "$DATA_DIR"
mkdir -p "$OUT_DIR"
mkdir -p "$TMP_DIR"


echo "Cleaning input format"
# runs python script clean_fasta.py; $1 = first argument (input fasta file); output goes to data/RAW.fasta
# clean_fasta.py cleans a non standard fasta file (removes extra info from headline)
python "$SCRIPT_DIR/clean_fasta.py" "$IN_DIR/$1" "$DATA_DIR/RAW.fasta"

echo "Filtering sequences"
# Input: cleaned fasta; output: filtered fasta; sequences between 26 and 2000 amino acids
python "$SCRIPT_DIR/generalFilter.py" \
    "$DATA_DIR/RAW.fasta" \
    "$DATA_DIR/filtered.fasta" \
    26 2000

echo "Removing duplicates"
# remove identical sequences; input: filtered fasta; output: deduplicated.fasta
python "$SCRIPT_DIR/duplicatesFilter.py" \
    "$DATA_DIR/filtered.fasta" \
    "$DATA_DIR/deduplicated.fasta"

echo 'Clustering data'
# filtered fasta to clustered fasta
# Copy file; perserve the original file and create a working file for the next step
cp "$DATA_DIR/deduplicated.fasta" "$DATA_DIR/precluster.fasta"

# Create clean working directory for MMSeqs
# Deletes the entire tmp directory and everything inside it
rm -rf "$TMP_DIR" # rm = remove; -r=  recursive (delete folder + all contents); -f = force (no prompts, ignore missing files) 
# Create a fresh, empty tmp directory
mkdir -p "$TMPDIR" # mkdir = make directory; -p= create parent dirs if needed, don't error if it exists

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
        "$DATA_DIR/precluster.fasta" \
        "$DATA_DIR/clusterRes" \
        "$TMP_DIR" \
        --min-seq-id 0.25 --cov-mode 1 -c 0.9
    
    echo "MMseqs output files 1.0"
    ls -lh "$DATA_DIR" | grep clusterRes
        
    #!This part only if ZK-448 proteins are removed (nucleotides and metals)
    echo "removing mapped clusters"
    # extract 2nd column of ZK-488_*.csv and skip the header
    tail -n +2 "$IN_DIR/ZK-448/ZK-448_nucleotides.csv" | cut -d',' -f2 > pdb_ZK448.txt
    python "$SCRIPT_DIR/eliminateClusters.py" \
     "$DATA_DIR/clusterRes_all_seqs.fasta" \
     pdb_ZK488.txt\
     "$DATA_DIR/filtered_mapped.fasta" 

    echo 'Fetching representatives'
    # Re-cluster the cleaned dataset; Produce representative clusters quickly without making an unnecessary big dataset; 
    # Output clusterRes2_* = final clustered dataset
    mmseqs easy-linclust \
        "$DATA_DIR/filtered_mapped.fasta" \
       "$DATA_DIR/clusterRes2" \
        "$DATA_DIR/tmp" \
        --min-seq-id 0.25 --cov-mode 1 -c 0.9
        
    echo "MMseqs output files 2.0"
    ls -lh "$DATA_DIR" | grep clusterRes2
        
    #! End of ZK-488 protein removal

    # Avoids conflicts later
    module unload netsurfp
)

echo "Finalizing dataset"
# If clusterRes2 doesn't exist, fall back to clusterRes
if [ -f "$DATA_DIR/clusterRes2_all_seqs.fasta" ]; then
    FINAL_ALL="$DATA_DIR/clusterRes2_all_seqs.fasta"
    FINAL_REP="$DATA_DIR/clusterRes2_rep_seq.fasta"
else
    FINAL_ALL="$DATA_DIR/clusterRes_all_seqs.fasta"
    FINAL_REP="$DATA_DIR/clusterRes_rep_seq.fasta"
fi

# All sequences grouped into final clusters
cp "$FINAL_ALL" "$OUT_DIR/finalclus.fasta"
# Contains one representative per cluster -> non-redundant dataset
cp "$FINAL_REP" "$OUT_DIR/final.fasta"

echo "Set fasta to BioLip"
python "$SCRIPT_DIR/fasta2BioLiP.py" \
  "$OUT_DIR/final.fasta" \
  "$IN_DIR/$2"

# Output is final2.0.csv
echo "Set dataset to PIPENN format"
python "$SCRIPT_DIR/BioLiP2PIPENN.py" \
  "$OUT_DIR/final.csv"

echo "Splits in training and test set"
python "$SCRIPT_DIR/split7030.py" "$OUT_DIR/final2.csv"

echo "Copy output to Reza file"
mkdir -p "$DATA_DIR/Reza"

cp "$OUT_DIR/final2_training.csv" "$OUT_DIR/Reza/prepared_biolip_win_n_training.csv"
cp "$OUT_DIR/final2_testing.csv" "$OUT_DIR/Reza/prepared_biolip_win_n_testing.csv"

