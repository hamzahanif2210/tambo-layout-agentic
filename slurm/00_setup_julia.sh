#!/bin/bash
#SBATCH -J tambo_jl_setup
#SBATCH -p serial_requeue
#SBATCH -c 1
#SBATCH --mem=16g
#SBATCH -t 3:00:00
#SBATCH -o /n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout/logs/%x-%j.out
#SBATCH --open-mode=append
#SBATCH --requeue
# Instantiate + precompile LlhamaFit's examples env (pulls TamboSim v1.0.1 etc. into ~/.julia).
# Touches only ~/.julia and the LlhamaFit clone's Manifest; never official_simulations.
set -euo pipefail
export JULIA_PKG_USE_CLI_GIT=true   # TamboSim source is git@github.com (SSH)
cd /n/home04/hhanif/tambo-org/LlhamaFit.jl
julia --version
time julia --project=examples -e 'using Pkg; Pkg.instantiate(); Pkg.precompile(); using TamboSim, LlhamaFit; println("OK: TamboSim + LlhamaFit load")'
