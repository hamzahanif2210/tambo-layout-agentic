#!/bin/bash
# After both v3 DE drivers finish: top-4 trial plans per budget -> layouts -> exact scoring -> forced/raw bootstrap.
set -eu
W=/n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout; cd $W
J=/n/home04/hhanif/tambo-org/LlhamaFit.jl/examples
P=plans/de3_winners; mkdir -p $P
for b in 300 917; do
  d=search/de_b${b}_s1_forced_raw
  tail -n +2 $d/log.csv | sort -t, -k3 -gr | awk -F, '$5==0' | awk "NR<=4" | while IFS=, read g i rest; do
    gd=$(printf "gen_%03d" $g); [ -d $d/${gd}_r ] && gd=${gd}_r
    cp $d/$gd/$(printf "ind_%02d.toml" $i) $P/de3_${b}_$(printf "g%03d_i%02d" $g $i).toml
  done
done
ls $P
lid=$(sbatch --parsable -J mklayouts3 -p serial_requeue -c 1 --mem=16g -t 1:00:00 -o $W/logs/%x-%j.out --wrap "julia --project=$J $W/scripts/07_site_and_placer.jl $W $(ls $W/$P/*.toml | tr '\n' ' ')")
T=slurm/tasks_score_v3_de3.txt; : > $T
for p in $P/*.toml; do n=$(basename $p .toml); for s in $(awk '{print $2}' slurm/tasks_score_v3_official.txt | sort -u); do echo "file:$W/layouts/$n.jls $s exact" >> $T; done; done
sid=$(sbatch --parsable --dependency=afterok:$lid --array=1-$(wc -l < $T) --export=ALL,TASKS=$W/$T slurm/04_score.sh)
for est in raw; do
  sbatch --parsable -J boot3 --dependency=afterok:$sid -p serial_requeue -c 1 --mem=16g -t 3:00:00 --requeue -o $W/logs/%x-%j.out \
    --wrap "CHANNEL=forced EST=$est julia --project=$J $W/scripts/14_bootstrap.jl $W 90 0.1 1e9 1000 train,select"
done
echo "layouts=$lid score=$sid"
