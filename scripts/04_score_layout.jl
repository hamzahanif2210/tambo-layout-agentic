# Exact (final-tier) scoring of one layout over every library chunk of one set.
#
#   julia --project=<env> 04_score_layout.jl <layout_spec> <library_set_dir> <outdir> [mapper]
#
# layout_spec:
#   official125           official volcano_valley_d_125m_large array
#   official150           official volcano_valley_d_150m_large array
#   hex:<spacing_m>       TamboSim.place_detector_units at that spacing (35 deg slope cut)
#   random:<frac>:<seed>  random subset of official125
# mapper: exact (default, final tier) | fast (search tier, R = 5 m)
#
# Writes (netscratch only): <outdir>/results_v3/<mapper>/<layout>/<set>.csv  (one row per event)
#                           <outdir>/results_v3/<mapper>/<layout>/<set>.summary
# Idempotent: skips if the summary exists.

include(joinpath(@__DIR__, "LayoutScore.jl"))
using .LayoutScore, TamboSim, LlhamaFit, AstroParticleUnits, TOML, Random, Statistics, Printf, Serialization

spec, libset, outdir = ARGS[1], ARGS[2], ARGS[3]
mapper = length(ARGS) >= 4 ? ARGS[4] : "exact"
startswith(abspath(outdir), LayoutScore.OFFICIAL) && error("refusing to write inside official_simulations")

setname = replace(relpath(libset, joinpath(outdir, "library", "full")), "/" => "_")   # nu_tau_00011
lname = startswith(spec, "file:") ? "plan_" * splitext(basename(spec[6:end]))[1] : replace(spec, ":" => "_")
res = joinpath(outdir, "results_v3", mapper, lname); mkpath(res)
csvp, sump = joinpath(res, "$setname.csv"), joinpath(res, "$setname.summary")
isfile(sump) && (println("exists: ", sump); exit(0))

chunks = sort(filter(f -> endswith(f, ".jls"), readdir(libset; join = true)))
isempty(chunks) && error("no library chunks in $libset")
simset = LayoutScore.load_library(chunks[1]).simset
geo = LayoutScore.load_geometry(simset)
d, cs, topo = geo.d_frame, geo.cs, geo.topo_bvh
units125 = copy(d["detector_units"])

const GCD = joinpath(LayoutScore.OFFICIAL, "gcd", "volcano_valley")
units = if spec == "official125"
    units125
elseif spec == "official150"
    tf150 = load_frames([geo.gc_path, joinpath(GCD, "volcano_valley_d_150m_large.jld2")])
    fr = collect(tf150)
    copy(fr[findlast(f -> f.stream == 'D' && haskey(f.data, "detector_units"), fr)]["detector_units"])
elseif startswith(spec, "hex:")
    sp = parse(Float64, split(spec, ":")[2])
    g = only(filter(f -> f.stream == 'G', collect(geo.tf)))
    obbs, _ = TamboSim.place_detector_units(g, d; spacing = sp * u"m")
    obbs
elseif startswith(spec, "file:")                  # Vector{OBB} written by 07_site_and_placer.jl
    open(deserialize, spec[6:end])
elseif startswith(spec, "random:")
    _, fr_, sd = split(spec, ":")
    rng = MersenneTwister(parse(Int, sd))
    n = round(Int, parse(Float64, fr_) * length(units125))
    units125[sort(randperm(rng, length(units125))[1:n])]
else
    error("unknown layout spec $spec")
end
bvh = TamboSim.BVHTree(units)
d["detector_units"] = units; d["detector_unit_bvh"] = bvh
geom = LlhamaFit.DetectorGeometry(d)
println("layout=$spec modules=$(length(units)) set=$setname chunks=$(length(chunks)) mapper=$mapper")

# In a function: top-level loop counters would hit Julia's soft-scope trap.
function score_all(chunks, units, bvh, cs, topo, geom, mapper, path)
    ntot = ntrig = nfit = 0; dψall = Float64[]; tmap = tfit = 0.0
    open(path, "w") do io
        println(io, "run,chunk_start,event_index,event_id,energy_GeV,theta,phi,n_hit_modules,triggered,fitted,dpsi_deg,reco_theta,reco_phi,sigma_psi_deg,u3")
        for cp in chunks
            lib = LayoutScore.load_library(cp)
            tmap += @elapsed hits = mapper == "exact" ?
                LayoutScore.map_hits_exact(lib, units, bvh, cs, topo) :
                first(LayoutScore.map_hits_fast(lib, LayoutScore.GridIndex(lib), units, bvh, cs, topo; R = 5.0))
            tabs = LayoutScore.event_tables(lib, hits, units)
            tfit += @elapsed s = LayoutScore.score(lib, tabs, geom)
            dψ = Dict(zip(s.ids, s.dψ)); rec = Dict(zip(s.ids, zip(s.rθ, s.rϕ, s.sψ, s.u3)))
            for i in 1:lib.nev
                println(io, join((lib.run, first(lib.evrange), first(lib.evrange) + i - 1, lib.event_id[i],
                                  lib.truthE[i], lib.truthθ[i], lib.truthϕ[i], length(tabs[i]),
                                  Int(s.trig[i]), Int(haskey(dψ, i)), get(dψ, i, NaN), get(rec, i, (NaN, NaN, NaN, NaN))...), ","))
            end
            ntot += lib.nev; ntrig += s.ntrig; nfit += s.nfit; append!(dψall, s.dψ)
        end
    end
    return (; ntot, ntrig, nfit, dψall, tmap, tfit)
end
(; ntot, ntrig, nfit, dψall, tmap, tfit) = score_all(chunks, units, bvh, cs, topo, geom, mapper, csvp * ".tmp")
# Guard against silent partial writes (e.g. a full quota left NUL-padded, truncated
# CSVs while the summary still got written): verify before publishing.
let nlines = countlines(csvp * ".tmp"), nnul = count(==(0x00), read(csvp * ".tmp"))
    (nlines == ntot + 1 && nnul == 0) ||
        error("CSV write check failed: $(nlines) lines for $(ntot) events, $(nnul) NUL bytes — not publishing")
end
mv(csvp * ".tmp", csvp; force = true)
q(p) = isempty(dψall) ? NaN : quantile(dψall, p)
summary = @sprintf("layout=%s modules=%d set=%s mapper=%s events=%d triggered=%d trig_frac=%.5f fitted=%d dpsi_median_deg=%.3f dpsi_68_deg=%.3f map_s=%.0f fit_s=%.0f",
                   spec, length(units), setname, mapper, ntot, ntrig, ntrig / max(ntot, 1), nfit, q(0.5), q(0.68), tmap, tfit)
println(summary)
open(io -> println(io, summary), sump, "w")
