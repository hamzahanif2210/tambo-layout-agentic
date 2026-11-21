#!/bin/bash
#SBATCH -J de_driver
#SBATCH -p serial_requeue
#SBATCH --requeue
#SBATCH -c 1
#SBATCH --mem=4g
#SBATCH -t 2-00:00:00
#SBATCH -o /n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout/logs/%x-%j.out
#SBATCH --open-mode=append
# Long-lived orchestrator: writes plans, submits eval arrays + reward jobs, waits.
# Resumable from search/de_b<budget>_s<seed>/state.jls. Not on serial_requeue (driver must not be preempted mid-wait... it would resume anyway).
set -euo pipefail
W=/n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout
julia --project=/n/home04/hhanif/tambo-org/LlhamaFit.jl/examples $W/scripts/12_de_driver.jl "$W" "$BUDGET" "${NP:-20}" "${GENS:-25}" "${SEED:-1}"
