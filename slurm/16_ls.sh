#!/bin/bash
#SBATCH -J ls_driver
#SBATCH -p serial_requeue
#SBATCH --requeue
#SBATCH -c 1
#SBATCH --mem=4g
#SBATCH -t 2-00:00:00
#SBATCH -o /n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout/logs/%x-%j.out
#SBATCH --open-mode=append
# Local-search orchestrator (resumable). Env: NAME, INIT, LAMBDA, GENS, SEED, TASKFILE, NPLANPART, TAMBO_MAX_SLOPE
set -euo pipefail
W=/n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout
julia --project=/n/home04/hhanif/tambo-org/LlhamaFit.jl/examples $W/scripts/16_local_search.jl "$W" "$NAME" "$INIT" "${LAMBDA:-20}" "${GENS:-30}" "${SEED:-1}"
