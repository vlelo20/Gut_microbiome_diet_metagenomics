#!/bin/bash
#SBATCH --job-name=sra_download
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --account=def-cottenie

module load sra-toolkit/3.0.9

cd /scratch/$USER/assign3

# Configure sra-toolkit to not use external services
vdb-config --set /repository/user/main/public/root=/scratch/$USER/sra_cache

SAMPLES=(SRR8146951 SRR8146952 SRR8146944 SRR8146936 SRR8146935 SRR8146938)

for SRR in "${SAMPLES[@]}"; do
    echo "Downloading $SRR..."
    fasterq-dump $SRR --split-files --threads 4 --skip-technical
    gzip ${SRR}_*.fastq
    echo "Done with $SRR"
done

echo "All downloads complete!"
