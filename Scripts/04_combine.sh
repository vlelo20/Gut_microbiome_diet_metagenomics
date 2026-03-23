#!/bin/bash
#SBATCH --job-name=combine
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --account=def-cottenie
#SBATCH --output=%x_%j_output.txt
#SBATCH --error=%x_%j_error.txt

module load bracken/3.0
cd /scratch/$USER/assign3

# Combine all bracken reports into one BIOM-style table
python3 combine_bracken_outputs.py \
    --files bracken_output/*.bracken \
    --output combined_bracken.txt

echo "Combined table created!"
