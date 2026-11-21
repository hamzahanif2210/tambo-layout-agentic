# Tag every event of one library set with its injection channel.
#
#   julia --project=<env> 15_point_types.jl <library_set_dir> <outdir>
#
# Writes <outdir>/weights/<set>.ptype.csv: run,event_index,ptype
#   forced   = ForcedNeutrinoInteractionPoint (ν reached the region, interaction forced; OW carries P_int)
#   upstream = UpstreamNeutrinoInteractionPoint (ν→τ converted naturally in Earth transit; OW = 1/surface pdf)
#   surface  = SurfaceCRPoint (cosmic rays)          none = no phase-space point
include(joinpath(@__DIR__, "LayoutScore.jl"))
using .LayoutScore, Printf
libset, outdir = ARGS[1], ARGS[2]
startswith(abspath(outdir), LayoutScore.OFFICIAL) && error("refusing to write inside official_simulations")
setname = replace(relpath(libset, joinpath(outdir, "library", "full")), "/" => "_")
out = joinpath(outdir, "weights", "$setname.ptype.csv")
tag(p) = isnothing(p) ? "none" : p isa LayoutScore.TamboSim.ForcedNeutrinoInteractionPoint ? "forced" :
         p isa LayoutScore.TamboSim.UpstreamNeutrinoInteractionPoint ? "upstream" : "surface"
cnt = Dict{String,Int}()
open(out * ".tmp", "w") do io
    println(io, "run,event_index,ptype")
    for c in sort(filter(f -> endswith(f, ".jls"), readdir(libset; join = true)))
        lib = LayoutScore.load_library(c)
        for i in 1:lib.nev
            t = tag(lib.psp[i]); cnt[t] = get(cnt, t, 0) + 1
            println(io, lib.run, ",", first(lib.evrange) + i - 1, ",", t)
        end
    end
end
mv(out * ".tmp", out; force = true)
println(setname, " ", cnt)
