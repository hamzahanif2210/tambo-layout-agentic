#!/bin/bash
# restart_after_gen.sh <driver_jobid> <run_dir_name> <job_name> <export_list>
# Waits until the driver logs one more completed generation, then cancels it and any eval/reward
# jobs it already submitted for the next generation, and resubmits it (resumes from state.jls).
set -u
W=/n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout; cd $W
jid=$1; run=$2; name=$3; exp=$4; log=search/$run/log.csv
last(){ [ -f $log ] && tail -1 $log | cut -d, -f1 || echo -1; }
g0=$(last)
while [ "$(last)" = "$g0" ]; do squeue -j $jid -h >/dev/null 2>&1 && [ -n "$(squeue -j $jid -h)" ] || { echo "$name finished/gone"; exit 0; }; sleep 20; done
scancel $jid; sleep 3
for j in $(squeue -u $USER -h -o "%A" | sort -u); do sacct -j $j -X -n -o SubmitLine%400 | head -1 | grep -q "search/$run/" && scancel $j && echo "cancelled $j"; done
sbatch --parsable -J $name --export=ALL,$exp slurm/12_de.sh
