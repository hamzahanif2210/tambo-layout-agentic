#!/bin/bash
#SBATCH -J analyze
#SBATCH -p serial_requeue
#SBATCH -c 1
#SBATCH --mem=8g
#SBATCH -t 1:00:00
#SBATCH -o /n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout/logs/%x-%j.out
set -euo pipefail
W=/n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout
for m in exact fast; do
  [ -d $W/results/$m ] && julia --project=/n/home04/hhanif/tambo-org/LlhamaFit.jl/examples $W/scripts/06_analyze.jl "$W" "$m"
done
