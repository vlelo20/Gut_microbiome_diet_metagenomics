#!/bin/bash
#SBATCH --job-name=bracken
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --account=def-cottenie
#SBATCH --output=%x_%j_output.txt
#SBATCH --error=%x_%j_error.txt

module load bracken/3.0

cd /scratch/$USER/assign3

mkdir -p bracken_output

DB=/scratch/$USER/assign3/kraken_db

SAMPLES=(SRR8146951 SRR8146952 SRR8146944 SRR8146936 SRR8146935 SRR8146938)

for SRR in "${SAMPLES[@]}"; do
    echo "Running Bracken on $SRR..."
    bracken \
        -d $DB \
        -i kraken_output/${SRR}.report \
        -o bracken_output/${SRR}.bracken \
        -r 150 \
        -l S \
        -t 10
    echo "Done: $SRR"
done

echo "All Bracken complete!"
