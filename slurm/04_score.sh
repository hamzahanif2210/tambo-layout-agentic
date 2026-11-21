#!/bin/bash
#SBATCH -J score
#SBATCH -p serial_requeue
#SBATCH -c 1
#SBATCH --mem=16g
#SBATCH -t 6:00:00
#SBATCH -o /n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout/logs/score/%A_%a.out
#SBATCH --open-mode=append
#SBATCH --requeue
# Reads the netscratch library (+ geometry from official_simulations, read-only); writes results/ only.
set -euo pipefail
W=/n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout
read LAYOUT LIBSET MAPPER < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASKS")
echo "task=$SLURM_ARRAY_TASK_ID $LAYOUT $LIBSET $MAPPER host=$(hostname)"
julia --project=/n/home04/hhanif/tambo-org/LlhamaFit.jl/examples $W/scripts/04_score_layout.jl "$LAYOUT" "$LIBSET" "$W" "$MAPPER"
