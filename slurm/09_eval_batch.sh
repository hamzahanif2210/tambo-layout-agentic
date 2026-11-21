#!/bin/bash
#SBATCH -J evalbatch
#SBATCH -p serial_requeue
#SBATCH -c 1
#SBATCH --mem=20g
#SBATCH -t 1:00:00
#SBATCH -o /n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout/logs/eval/%A_%a.out
#SBATCH --open-mode=append
#SBATCH --requeue
set -euo pipefail
W=/n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout
read -r LIBSET PART NPARTS <<< "$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASKS")"
julia --project=/n/home04/hhanif/tambo-org/LlhamaFit.jl/examples $W/scripts/09_eval_batch.jl "$W" "$BATCH" "$LIBSET" "${SPLIT:-train}" "${PART:-1}" "${NPARTS:-1}"
