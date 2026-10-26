#!/bin/bash
#SBATCH --job-name=scaffolds
#SBATCH --output=logs/scaffolds_%j.out
#SBATCH --error=logs/scaffolds_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=48:00:00
# --- Cluster specific: uncomment and fill in if your cluster requires them ---
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>

# Number of unique Bemis-Murcko scaffolds of a Parquet library:
#   mkdir -p logs
#   sbatch run_calculate_unique_scaffolds_database.sh /path/to/Enamine_REAL_65B_parquet results/scaffolds/enamine_65B
#   sbatch run_calculate_unique_scaffolds_database.sh /path/to/Enamine_REAL_5B_parquet  results/scaffolds/enamine_5B
#

INPUT_DIR=$1
OUTPUT_DIR=$2

# --- Cluster specific: load the Python environment (see environment.yml) ---
module purge
module load anaconda
conda activate alscreen-analysis

python -u calculate_unique_scaffolds_database.py \
    --input_dir "$INPUT_DIR" \
    --output_dir "$OUTPUT_DIR" \
    --smiles_column SMILES \
    --workers "$SLURM_CPUS_PER_TASK"
