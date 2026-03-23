#!/bin/bash
#SBATCH --job-name=qc_fastp
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --account=def-cottenie
#SBATCH --output=%x_%j_output.txt
#SBATCH --error=%x_%j_error.txt

module load fastp/0.24.0

cd /scratch/$USER/assign3

mkdir -p qc_output

# Sample metadata
declare -A SAMPLES
SAMPLES=(
    [SRR8146951]="vegan"
    [SRR8146952]="vegan"
    [SRR8146944]="vegan"
    [SRR8146936]="omnivore"
    [SRR8146935]="omnivore"
    [SRR8146938]="omnivore"
)

for SRR in "${!SAMPLES[@]}"; do
    echo "Running QC on $SRR (${SAMPLES[$SRR]})..."
    fastp \
        --in1 data/${SRR}_1.fastq.gz \
        --in2 data/${SRR}_2.fastq.gz \
        --out1 qc_output/${SRR}_1_clean.fastq.gz \
        --out2 qc_output/${SRR}_2_clean.fastq.gz \
        --html qc_output/${SRR}_report.html \
        --json qc_output/${SRR}_report.json \
        --thread 8 \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --length_required 50
    echo "Done QC: $SRR"
done

echo "All QC complete!"
