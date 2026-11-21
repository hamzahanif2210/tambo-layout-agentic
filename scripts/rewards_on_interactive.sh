#!/bin/bash
# Run pending reward/llm_reward jobs on an interactive allocation ($1), then cancel the queued copies so drivers proceed.
IJ=$1; W=/n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout
for j in $(squeue -u $USER -h -t PD -n reward,llm_reward -o %i); do
  cmd=$(scontrol show job $j -o | grep -oP 'SubmitLine=.*' | grep -oP '\-\-wrap \K.*' | sed 's/^"//; s/"$//')
  echo "$(date +%H:%M) job $j: $cmd" >> $W/logs/rewards_interactive.log
  srun --jobid=$IJ --overlap -n1 bash -c "$cmd" >> $W/logs/rewards_interactive.log 2>&1 && scancel $j
done
