#!/bin/bash
#SBATCH --job-name=krona
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --account=def-cottenie
#SBATCH --output=%x_%j_output.txt
#SBATCH --error=%x_%j_error.txt

module load kraken2/2.1.6
module load kronatools/2.8.1

cd /scratch/$USER/assign3
mkdir -p krona_output

SAMPLES=(SRR8146951 SRR8146952 SRR8146944 SRR8146936 SRR8146935 SRR8146938)

for SRR in "${SAMPLES[@]}"; do
    echo "Creating Krona plot for $SRR..."
    python kreport2krona.py \
        -r kraken_output/${SRR}.report \
        -o krona_output/${SRR}.krona
    echo "Done: $SRR"
done

# Combine all samples into one interactive plot
ktImportText \
    krona_output/*.krona \
    -o krona_output/all_samples_krona.html

echo "Krona plots complete!"
