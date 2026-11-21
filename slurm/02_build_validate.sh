#!/bin/bash
#SBATCH -p serial_requeue
#SBATCH -c 1
#SBATCH --mem=32g
#SBATCH -t 6:00:00
#SBATCH -o /n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout/logs/%x-%j.out
#SBATCH --open-mode=append
#SBATCH --requeue
# Read-only on official_simulations; writes only under $W (library + log).
set -euo pipefail
W=/n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout
time julia --project=/n/home04/hhanif/tambo-org/LlhamaFit.jl/examples \
  $W/scripts/02_build_and_validate.jl "$SIMSET" "$RUN" "$MAXEV" "$W"
