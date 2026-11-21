#!/bin/bash
#SBATCH -J tambo_feas
#SBATCH -p serial_requeue
#SBATCH -c 1
#SBATCH --mem=32g
#SBATCH -t 4:00:00
#SBATCH -o /n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout/logs/%x-%j.out
#SBATCH --open-mode=append
#SBATCH --requeue
# Read-only on official_simulations; writes only this log.
set -euo pipefail
W=/n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout
SIMSET=${SIMSET:-/n/holylfs05/LABS/arguelles_delgado_lab/Lab/TAMBO/official_simulations/datasets/v1/volcano_valley/nu_tau/00019}
RUN=${RUN:-00000}
MAXEV=${MAXEV:-1000}
time julia --project=/n/home04/hhanif/tambo-org/LlhamaFit.jl/examples $W/scripts/01_feasibility_timing.jl "$SIMSET" "$RUN" "$MAXEV"
