#!/bin/bash
#SBATCH -J lib_full
#SBATCH -p serial_requeue
#SBATCH -c 1
#SBATCH --mem=16g
#SBATCH -t 2:00:00
#SBATCH -o /n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout/logs/lib_full/%A_%a.out
#SBATCH --open-mode=append
#SBATCH --requeue
# Read-only on official_simulations; writes only under $W/library/full.
set -euo pipefail
W=/n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout
read SIMSET RUN START SIZE < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASKS")
echo "task=$SLURM_ARRAY_TASK_ID $SIMSET $RUN $START $SIZE host=$(hostname)"
julia --project=/n/home04/hhanif/tambo-org/LlhamaFit.jl/examples $W/scripts/03_build_chunk.jl "$SIMSET" "$RUN" "$START" "$SIZE" "$W"
