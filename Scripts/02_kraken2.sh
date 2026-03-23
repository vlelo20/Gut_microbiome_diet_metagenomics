#!/bin/bash
#SBATCH --job-name=kraken2_classify
#SBATCH --time=24:00:00
#SBATCH --mem=200G
#SBATCH --cpus-per-task=32
#SBATCH --account=def-cottenie
#SBATCH --output=%x_%j_output.txt
#SBATCH --error=%x_%j_error.txt

module load kraken2/2.1.6

cd /scratch/$USER/assign3
mkdir -p kraken_output

DB=/scratch/$USER/assign3/kraken_db

SAMPLES=(SRR8146951 SRR8146952 SRR8146944 SRR8146936 SRR8146935 SRR8146938)

for SRR in "${SAMPLES[@]}"; do
    echo "Classifying $SRR..."
    kraken2 \
        --db $DB \
        --confidence 0.15 \
        --threads 32 \
        --output kraken_output/${SRR}.kraken \
        --report kraken_output/${SRR}.report \
        --paired \
        qc_output/${SRR}_1_clean.fastq.gz \
        qc_output/${SRR}_2_clean.fastq.gz
    echo "Done: $SRR"
    rm kraken_output/${SRR}.kraken
done

echo "All classification complete!"
