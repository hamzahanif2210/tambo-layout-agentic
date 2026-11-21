#!/bin/bash
# Preliminary exact scoring + select bootstrap of the current best LLM plans (2 per budget).
set -eu
W=/n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout; cd $W
J=/n/home04/hhanif/tambo-org/LlhamaFit.jl/examples
P=plans/llm_prelim; mkdir -p $P
for b in 917 300; do
  d=search/llm_b${b}_glm-5
  tail -n +2 $d/log.csv | awk -F, '$5==0 && $8=="llm"' | sort -t, -k3 -gr | awk "NR<=2" | while IFS=, read r i rest; do
    rd=$(printf "round_%03d" $r); cp $d/$rd/$i.toml $P/llm${b}_r$(printf "%02d" $r)_$i.toml
  done
done
ls $P
lid=$(sbatch --parsable -J mklayouts6 -p serial_requeue --requeue -c 1 --mem=16g -t 1:00:00 -o $W/logs/%x-%j.out \
  --wrap "TAMBO_MAX_SLOPE=50 SITE=site_s50 julia --project=$J $W/scripts/07_site_and_placer.jl $W $(ls $W/$P/*.toml | tr '\n' ' ')")
T=slurm/tasks_score_v3_llm.txt; : > $T
for n in $(ls $P/*.toml | xargs -n1 basename | sed 's/.toml$//'); do
  for s in $(awk '{print $2}' slurm/tasks_score_v3_official.txt | sort -u); do echo "file:$W/layouts/$n.jls $s exact" >> $T; done
done
sid=$(sbatch --parsable --dependency=afterok:$lid --array=1-$(wc -l < $T) --mem=20g --export=ALL,TASKS=$W/$T slurm/04_score.sh)
bid=$(sbatch --parsable -J boot6 --dependency=afterok:$sid -p serial_requeue -c 1 --mem=16g -t 3:00:00 --requeue -o $W/logs/%x-%j.out \
  --wrap "CHANNEL=forced EST=raw julia --project=$J $W/scripts/14_bootstrap.jl $W 90 0.1 1e9 1000 train,select")
echo "layouts=$lid score=$sid boot=$bid tasks=$(wc -l < $T)"
