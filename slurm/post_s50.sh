#!/bin/bash
# After the steep-aware (≤50°) runs: top-3 leak-free DE plans per budget + local-search bests
# -> layouts -> exact scoring -> forced/raw paired bootstrap (train, select).
set -eu
W=/n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout; cd $W
J=/n/home04/hhanif/tambo-org/LlhamaFit.jl/examples
P=plans/de5_s50_winners; mkdir -p $P
for b in 300 917; do
  d=search/de_b${b}_s5_forced_raw_s50
  tail -n +2 $d/log.csv | sort -t, -k3 -gr | awk -F, '$5==0' | awk "NR<=3" | while IFS=, read g i rest; do
    gd=$(printf "gen_%03d" $g); [ -d $d/${gd}_r ] && gd=${gd}_r
    cp $d/$gd/$(printf "ind_%02d.toml" $i) $P/de5s50_${b}_$(printf "g%03d_i%02d" $g $i).toml
  done
done
cp search/ls_off125_s50/best_units.jls layouts/ls_off125_s50.jls
cp search/ls_sub300_s50/best_units.jls layouts/ls_sub300_s50.jls
ls $P
lid=$(sbatch --parsable -J mklayouts5 -p serial_requeue --requeue -c 1 --mem=16g -t 1:00:00 -o $W/logs/%x-%j.out \
  --wrap "TAMBO_MAX_SLOPE=50 SITE=site_s50 julia --project=$J $W/scripts/07_site_and_placer.jl $W $(ls $W/$P/*.toml | tr '\n' ' ')")
T=slurm/tasks_score_v3_s50.txt; : > $T
for n in $(ls $P/*.toml | xargs -n1 basename | sed 's/.toml$//') ls_off125_s50 ls_sub300_s50; do
  for s in $(awk '{print $2}' slurm/tasks_score_v3_official.txt | sort -u); do echo "file:$W/layouts/$n.jls $s exact" >> $T; done
done
sid=$(sbatch --parsable --dependency=afterok:$lid --array=1-$(wc -l < $T) --export=ALL,TASKS=$W/$T slurm/04_score.sh)
bid=$(sbatch --parsable -J boot5 --dependency=afterok:$sid -p serial_requeue -c 1 --mem=16g -t 3:00:00 --requeue -o $W/logs/%x-%j.out \
  --wrap "CHANNEL=forced EST=raw julia --project=$J $W/scripts/14_bootstrap.jl $W 90 0.1 1e9 1000 train,select")
echo "layouts=$lid score=$sid boot=$bid tasks=$(wc -l < $T)"
