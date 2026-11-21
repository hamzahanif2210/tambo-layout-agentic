# Search-tier evaluation of a BATCH of plans on ONE library set (fast mapper,
# one split only). Chunk-outer / plan-inner: each library chunk is loaded once
# and every plan in the batch is scored on it.
#
#   julia --project=<env> 09_eval_batch.jl <outdir> <batch_dir> <library_set_dir> [split=train] [part=1] [nparts=1]
#
# With nparts > 1 only every nparts-th chunk (starting at part) of the split is
# scored and the output is <set>.p<part>of<nparts>.csv; 10_reward.jl concatenates
# the parts (rates do not depend on the partition).
# <batch_dir>/*.toml are plans. Writes <batch_dir>/eval/<plan>/<set>.csv with one
# row per event of the split:
#   energy_GeV,rate_per_yr,triggered,fitted,reco_theta_deg,dpsi_deg,truth_theta_deg,u3,sigma_psi_deg,forced
# forced = 0 for UpstreamNeutrinoInteractionPoint events (ν→τ converted in Earth), else 1.
# rate_per_yr is already rescaled to the full set (1 / split fraction) and uses
# the reference flux for the set kind (flux_nu for nu_tau, flux_cr otherwise).
# Idempotent per (plan, set): skips plans whose CSV exists.
include(joinpath(@__DIR__, "LayoutScore.jl"))
include(joinpath(@__DIR__, "Placer.jl"))
using .LayoutScore, .Placer, TamboSim, LlhamaFit, AstroParticleUnits, TOML, Printf, Serialization

outdir, batch, libset = ARGS[1], ARGS[2], ARGS[3]
spl = length(ARGS) >= 4 ? Symbol(ARGS[4]) : :train
part = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 1
nparts = length(ARGS) >= 6 ? parse(Int, ARGS[6]) : 1
1 <= part <= nparts || error("bad part $part of $nparts")
startswith(abspath(outdir), LayoutScore.OFFICIAL) && error("refusing to write inside official_simulations")
const TAGS = Dict(:train => (0, 1, 2), :select => (3,), :test => (4,))
tags = TAGS[spl]; frac = length(tags) / 5
setname = replace(relpath(libset, joinpath(outdir, "library", "full")), "/" => "_")
signal = startswith(setname, "nu_tau")
flux = signal ? LayoutScore.flux_nu : LayoutScore.flux_cr

chunks = sort(filter(f -> endswith(f, ".jls"), readdir(libset; join = true)))
runs = unique([first(split(basename(c), "_")) for c in chunks]); by_run = length(runs) > 1
nruns = length(runs)
# run-split: drop whole chunks of non-split runs up front (no load needed)
by_run && (chunks = filter(c -> mod(parse(Int, first(split(basename(c), "_"))), 5) in tags, chunks))
chunks = chunks[part:nparts:end]                       # strided partition (balances cost)
outname = nparts == 1 ? "$setname.csv" : "$setname.p$(part)of$(nparts).csv"

# plans: *.toml (placed by Placer) or *.jls (a ready Vector{OBB}, e.g. from the local search)
plans = sort(filter(f -> endswith(f, ".toml") || endswith(f, ".jls"), readdir(batch; join = true)))
# PLANPART/NPLANPART env: score only every NPLANPART-th plan (parallelise over plans too)
let pp = parse(Int, get(ENV, "PLANPART", "1")), npp = parse(Int, get(ENV, "NPLANPART", "1"))
    global plans = plans[pp:npp:end]
end
geo = LayoutScore.load_geometry(LayoutScore.load_library(chunks[1]).simset)
d, cs, topo = geo.d_frame, geo.cs, geo.topo_bvh
region = Placer.Region(d, cs)

todo = []
for p in plans
    name = splitext(basename(p))[1]
    out = joinpath(batch, "eval", name, outname)
    isfile(out) && continue
    units, rep = endswith(p, ".jls") ? (open(deserialize, p), nothing) : Placer.place(region, p)
    isempty(units) && (mkpath(dirname(out)); write(out, "energy_GeV,rate_per_yr,triggered,fitted,reco_theta_deg,dpsi_deg,truth_theta_deg,u3,sigma_psi_deg,forced\n"); continue)
    bvh = TamboSim.BVHTree(units)
    d["detector_units"] = units; d["detector_unit_bvh"] = bvh
    push!(todo, (; name, out, units, bvh, geom = LlhamaFit.DetectorGeometry(d), io = IOBuffer()))
end
println("set=$setname part=$part/$nparts split=$spl chunks=$(length(chunks)) plans_to_eval=$(length(todo))")

function run_batch(chunks, todo)
    t0 = time()
    for (ci, c) in enumerate(chunks)
        lib = LayoutScore.subset_events(LayoutScore.load_library(c), tags; by_run)
        lib.nev == 0 && continue
        ow = LayoutScore.oneweights(lib)
        rate = [flux(lib.truthE[i]) * ow[i] / nruns / frac * LayoutScore.SECONDS_PER_YEAR for i in 1:lib.nev]
        grid = LayoutScore.GridIndex(lib)
        forced = [lib.psp[i] isa TamboSim.UpstreamNeutrinoInteractionPoint ? 0 : 1 for i in 1:lib.nev]
        for p in todo
            hits, _ = LayoutScore.map_hits_fast(lib, grid, p.units, p.bvh, cs, topo; R = 5.0)
            tabs = LayoutScore.event_tables(lib, hits, p.units)
            s = LayoutScore.score(lib, tabs, p.geom)
            rec = Dict(zip(s.ids, zip(s.rθ, s.dψ, s.u3, s.sψ)))
            for i in 1:lib.nev
                θ, dψ, u3, sψ = get(rec, i, (NaN, NaN, NaN, NaN))
                println(p.io, lib.truthE[i], ",", rate[i], ",", Int(s.trig[i]), ",", Int(haskey(rec, i)), ",",
                        rad2deg(θ), ",", dψ, ",", rad2deg(lib.truthθ[i]), ",", u3, ",", sψ, ",", forced[i])
            end
        end
        @printf("chunk %d/%d done, %.0f s elapsed\n", ci, length(chunks), time() - t0)
    end
end
run_batch(chunks, todo)
for p in todo
    mkpath(dirname(p.out))
    open(p.out * ".tmp", "w") do f
        println(f, "energy_GeV,rate_per_yr,triggered,fitted,reco_theta_deg,dpsi_deg,truth_theta_deg,u3,sigma_psi_deg,forced"); write(f, take!(p.io))
    end
    mv(p.out * ".tmp", p.out; force = true)
end
println("DONE")
