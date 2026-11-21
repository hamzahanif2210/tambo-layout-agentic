# Build one library chunk: events [start, start+size-1] of one official run.
#
#   julia --project=<env> 03_build_chunk.jl <simset_dir> <run> <start> <size> <outdir>
#
# Output (netscratch only): <outdir>/library/full/<particle>/<set>/<run>_<start>-<stop>.jls
# plus a one-line .ok summary. Idempotent: an existing chunk is skipped, so a
# requeued (serial_requeue) task is safe. A chunk past the run's last event
# writes <...>.empty and exits.
# If the run has official stored hits, the chunk validates itself: the exact
# mapper on the official 125 m layout must reproduce them (mismatch count logged).

include(joinpath(@__DIR__, "LayoutScore.jl"))
using .LayoutScore, TamboSim, LlhamaFit, Printf

simset, run = ARGS[1], ARGS[2]
start, size = parse(Int, ARGS[3]), parse(Int, ARGS[4])
outdir = ARGS[5]
startswith(abspath(outdir), LayoutScore.OFFICIAL) && error("refusing to write inside official_simulations")

stop = start + size - 1
tag = joinpath(basename(dirname(simset)), basename(simset))      # nu_tau/00011
base = joinpath(outdir, "library", "full", tag, @sprintf("%s_%05d-%05d", run, start, stop))
if isfile(base * ".jls") || isfile(base * ".empty")
    println("exists, skipping: ", base); exit(0)
end

t = @elapsed lib = LayoutScore.build_library(simset, run; evrange = start:stop)
if lib.nev <= 0
    mkpath(dirname(base)); touch(base * ".empty")
    println("past end of run (n_events_in_run=$(lib.n_events_in_run)); wrote .empty"); exit(0)
end
LayoutScore.save_library(base * ".jls", lib)
@printf("built %s events=%d (run has %d) raw=%d selected=%d build_s=%.0f MB=%.0f\n",
        base, lib.nev, lib.n_events_in_run, lib.n_raw, length(lib.x), t, filesize(base * ".jls") / 1e6)

mism = -1; nhit = -1
if lib.has_stored
    geo = LayoutScore.load_geometry(simset)
    d = geo.d_frame; units = d["detector_units"]; bvh = d["detector_unit_bvh"]
    hx = LayoutScore.map_hits_exact(lib, units, bvh, geo.cs, geo.topo_bvh)
    mx = LayoutScore.event_modules(lib, hx)
    mism = count(i -> mx[i] != lib.stored[i], eachindex(mx)); nhit = length(hx)
    @printf("validation: exact_vs_stored mismatched_events=%d/%d hits=%d\n", mism, lib.nev, nhit)
end
open(base * ".ok", "w") do io
    println(io, "set=$tag run=$run start=$start stop=$stop nev=$(lib.nev) n_events_in_run=$(lib.n_events_in_run) ",
                "selected=$(length(lib.x)) has_stored=$(lib.has_stored) mismatch=$mism hits=$nhit")
end
println("DONE")
