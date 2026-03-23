#!/bin/bash
#SBATCH --job-name=kraken_db
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --account=def-cottenie

mkdir -p /scratch/$USER/assign3/kraken_db
cd /scratch/$USER/assign3/kraken_db

echo "Downloading Kraken2 Standard-8 database..."
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08_GB_20260226.tar.gz

echo "Extracting..."
tar -xzf k2_standard_08_GB_20260226.tar.gz

echo "Done!"
