# Per-event OneWeights for every event of one library set (layout-independent).
#
#   julia --project=<env> 05_event_weights.jl <library_set_dir> <outdir>
#
# Writes <outdir>/weights/<set>.csv:
#   run,event_index,energy_GeV,theta,phi,oneweight_GeV_m2_sr,n_runs_in_set
# oneweight is per single run; divide by n_runs_in_set when pooling the set.

include(joinpath(@__DIR__, "LayoutScore.jl"))
using .LayoutScore, Printf

libset, outdir = ARGS[1], ARGS[2]
startswith(abspath(outdir), LayoutScore.OFFICIAL) && error("refusing to write inside official_simulations")
setname = replace(relpath(libset, joinpath(outdir, "library", "full")), "/" => "_")
chunks = sort(filter(f -> endswith(f, ".jls"), readdir(libset; join = true)))
runs = sort(unique([first(split(basename(c), "_")) for c in chunks]))
mkpath(joinpath(outdir, "weights"))
out = joinpath(outdir, "weights", "$setname.csv")

function write_weights(chunks, runs, path)
    n = 0; nzero = 0; s = 0.0
    open(path * ".tmp", "w") do io
        println(io, "run,event_index,energy_GeV,theta,phi,oneweight_GeV_m2_sr,n_runs_in_set")
        for c in chunks
            lib = LayoutScore.load_library(c)
            ow = LayoutScore.oneweights(lib)
            for i in 1:lib.nev
                println(io, join((lib.run, first(lib.evrange) + i - 1, lib.truthE[i], lib.truthθ[i],
                                  lib.truthϕ[i], ow[i], length(runs)), ","))
            end
            n += lib.nev; nzero += count(iszero, ow); s += sum(ow)
        end
    end
    mv(path * ".tmp", path; force = true)
    return n, nzero, s
end
n, nzero, s = write_weights(chunks, runs, out)
@printf("set=%s runs=%d events=%d zero_weight=%d sum_OW_per_run=%.4g GeV m^2 sr -> %s\n",
        setname, length(runs), n, nzero, s / length(runs), out)
