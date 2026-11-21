#!/bin/bash
#SBATCH -J weights
#SBATCH -p serial_requeue
#SBATCH -c 1
#SBATCH --mem=16g
#SBATCH -t 3:00:00
#SBATCH -o /n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout/logs/weights/%A_%a.out
#SBATCH --open-mode=append
#SBATCH --requeue
set -euo pipefail
W=/n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout
LIBSET=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASKS")
julia --project=/n/home04/hhanif/tambo-org/LlhamaFit.jl/examples $W/scripts/05_event_weights.jl "$LIBSET" "$W"
