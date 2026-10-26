#!/bin/bash
#SBATCH --job-name=mw_db
#SBATCH --output=logs/mw_%A_%a.out
#SBATCH --error=logs/mw_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=48:00:00
# --- Cluster specific: uncomment and fill in if your cluster requires them ---
#SBATCH --account=<your_account>
#SBATCH --qos=<your_qos>

# Molecular weights of a Parquet library. Each array task is one node that
# processes every N-th file, e.g. for Enamine REAL 65B split over 48 nodes:
#   mkdir -p logs
#   sbatch --array=0-47 run_calculate_mw_database.sh /path/to/Enamine_REAL_65B_parquet results/mw/enamine_65B
#   sbatch --array=0-6  run_calculate_mw_database.sh /path/to/Enamine_REAL_5B_parquet  results/mw/enamine_5B
#

INPUT_DIR=$1
OUTPUT_DIR=$2

# --- Cluster specific: load the Python environment (see environment.yml) ---
module purge
module load anaconda
conda activate alscreen-analysis

python -u calculate_mw_database.py \
    --input_dir "$INPUT_DIR" \
    --output_dir "$OUTPUT_DIR" \
    --smiles_column SMILES \
    --workers "$SLURM_CPUS_PER_TASK" \
    --task_id "$SLURM_ARRAY_TASK_ID" \
    --num_tasks "$SLURM_ARRAY_TASK_COUNT"
