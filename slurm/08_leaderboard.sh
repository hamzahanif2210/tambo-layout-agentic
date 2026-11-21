#!/bin/bash
#SBATCH -J leaderboard
#SBATCH -p serial_requeue
#SBATCH -c 1
#SBATCH --mem=16g
#SBATCH -t 1:00:00
#SBATCH -o /n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout/logs/%x-%j.out
set -euo pipefail
W=/n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout
julia --project=/n/home04/hhanif/tambo-org/LlhamaFit.jl/examples $W/scripts/08_leaderboard.jl "$W" exact official125
