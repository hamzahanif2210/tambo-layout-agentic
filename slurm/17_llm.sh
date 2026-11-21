#!/bin/bash
#SBATCH -J llm
#SBATCH -p serial_requeue
#SBATCH -c 1
#SBATCH --mem=2g
#SBATCH -t 1-00:00:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH -o /n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout/logs/%x-%j.out
# LLM plan-writing loop driver (calls NRP endpoint, submits eval/reward jobs). Resumable.
set -euo pipefail
W=/n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout
python3 -u $W/scripts/17_llm_driver.py "$W" "$BUDGET" "${ROUNDS:-20}" "${K:-40}" "${MODEL:-glm-5}"
