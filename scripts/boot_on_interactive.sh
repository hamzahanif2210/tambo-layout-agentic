#!/bin/bash
W=/n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout
J=/n/home04/hhanif/tambo-org/LlhamaFit.jl/examples
while squeue -j 62583565 -h 2>/dev/null | grep -q .; do sleep 60; done
srun --jobid=62583857 --overlap -n1 bash -c "cd $W; CHANNEL=forced EST=raw julia --project=$J scripts/14_bootstrap.jl $W 90 0.1 1e9 1000 train,select" > logs/boot6_interactive2.log 2>&1
echo DONE >> logs/boot6_interactive2.log
