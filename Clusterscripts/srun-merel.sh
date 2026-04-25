# Runs CLusterJob.sh <Input.fasta> <BioLiP_*.txt>
# Runs in ADA defq

srun -N 1 -p defq ./ClusterJob.sh peptide.fasta BioLiP_peptide.txt
#srun -N 1 -p binf  ./merel.sh
