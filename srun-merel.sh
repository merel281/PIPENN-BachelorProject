# Runs CLusterJob.sh <Input.fasta> 
# Runs in ADA defq

srun -N 1 -p defq ./ClusterJob.sh metal.fasta 
#srun -N 1 -p binf  ./merel.sh