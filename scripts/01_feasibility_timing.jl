# Milestone 1: can we re-score arbitrary layouts from existing CORSIKA surface
# particles, and how fast?
#
# READ-ONLY w.r.t. official_simulations: only load_frames / read_corsika_mesh
# are called on those paths. Nothing is saved there (no save_frames at all).
# Output: one TOML-ish summary printed to stdout (captured by the SLURM log).
#
#   julia --project=/n/home04/hhanif/tambo-org/LlhamaFit.jl/examples \
#         01_feasibility_timing.jl <simset_dir> <run_number> [max_events]

using TamboSim, TOML, Random, Statistics
using LlhamaFit
using AstroParticleUnits

const OFFICIAL = "/n/holylfs05/LABS/arguelles_delgado_lab/Lab/TAMBO/official_simulations"

simset = length(ARGS) >= 1 ? ARGS[1] : joinpath(OFFICIAL, "datasets/v1/volcano_valley/nu_tau/00019")
run    = length(ARGS) >= 2 ? ARGS[2] : "00000"
maxev  = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : typemax(Int)

gc_path  = TOML.parsefile(joinpath(simset, "injection.toml"))["geometry"]["gc_path"]
d_path   = TOML.parsefile(joinpath(simset, "hits_125m_large.toml"))["detector"]["d_path"]
hit_file = joinpath(simset, "hits_125m_large", "$run.jld2")
println("simset   = ", simset, "\nrun      = ", run)

t_load = @elapsed tf = load_frames([gc_path, d_path, hit_file])
qs = tf.q_frames[1:min(end, maxev)]
d_frame = qs[1].d_frame
g_frame = qs[1].m_frame.g_frame
cs, topo_bvh = g_frame["cs"], g_frame["bvh"]
units0 = copy(d_frame["detector_units"])
println("load_s   = ", round(t_load; digits = 1), "   events = ", length(qs),
        "   modules(125m_large) = ", length(units0))

# ---- 1. cache surface particles once (the expensive, layout-independent part)
stored_hits = [haskey(q.data, "corsika_hits") ? q["corsika_hits"] : nothing for q in qs]
cache = Vector{Vector{Any}}(undef, length(qs))
t_cache = @elapsed for (i, q) in enumerate(qs)
    evs = Any[]
    if haskey(q, "corsika_directories")
        for d in unique(dirname.(q["corsika_directories"]))
            isdir(d) || continue
            try
                append!(evs, TamboSim.read_corsika_mesh(d, cs))
            catch err
                @warn "skip" dir = d exception = err
            end
        end
    end
    cache[i] = evs
end
n_particles = sum(length, cache)
println("cache_s  = ", round(t_cache; digits = 1), "   surface_particles = ", n_particles,
        "   events_with_particles = ", count(!isempty, cache))

# ---- 2. per-layout scoring
function score_layout!(units, label)
    bvh = TamboSim.BVHTree(units)
    d_frame["detector_units"] = units
    d_frame["detector_unit_bvh"] = bvh
    geom = LlhamaFit.DetectorGeometry(d_frame)      # fresh, bypasses the per-frame cache

    t_map = @elapsed for (i, q) in enumerate(qs)
        hs = NamedTuple{(:particle, :module_index, :weight),
                        Tuple{TamboSim.Particle{Float64}, Int, Float64}}[]
        for ev in cache[i]
            h = TamboSim._corsika_module_hit(ev, bvh, topo_bvh)
            isnothing(h) || push!(hs, h)
        end
        q["corsika_hits"] = hs                       # in-memory only, never saved
    end

    local tabs
    t_tab = @elapsed tabs = [LlhamaFit.hits_table(q, geom.units) for q in qs]
    trig = [LlhamaFit.TRIGGER_5x1(h) for h in tabs]
    ntrig = count(trig)

    dψ = Float64[]; nvalid = 0
    t_fit = @elapsed for (q, h, ok) in zip(qs, tabs, trig)
        ok || continue
        r = try reconstruct(h, geom, PlaneTiming()) catch; nothing end
        isnothing(r) && continue
        r.valid && (nvalid += 1)
        inj = q["injection_initial_state"]
        push!(dψ, Float64(great_circle_distance(inj.direction.θ, inj.direction.ϕ, r.θ, r.ϕ) / °))
    end
    med = isempty(dψ) ? NaN : median(dψ)
    println("[$label] modules=$(length(units)) map_s=$(round(t_map; digits=2)) ",
            "table_s=$(round(t_tab; digits=2)) triggered=$ntrig/$(length(qs)) ",
            "fit_s=$(round(t_fit; digits=2)) fit_s_per_event=$(round(t_fit/max(ntrig,1); digits=4)) ",
            "valid=$nvalid median_dpsi_deg=$(round(med; digits=2))")
    return tabs
end

# full official layout: must reproduce the stored hits exactly (sanity check)
score_layout!(units0, "full_125m")
let ncmp = 0, mism = 0, nh_new = 0, nh_old = 0
    for (q, sh) in zip(qs, stored_hits)
        isnothing(sh) && continue
        ncmp += 1
        a = sort([h.module_index for h in q["corsika_hits"]])
        b = sort([h.module_index for h in sh])
        nh_new += length(a); nh_old += length(b)
        a == b || (mism += 1)
    end
    println("sanity: compared_events=$ncmp mismatched_events=$mism ",
            "hits_recomputed=$nh_new hits_stored=$nh_old")
end

rng = MersenneTwister(1)
for frac in (0.5, 0.25)
    sub = units0[sort(randperm(rng, length(units0))[1:round(Int, frac * length(units0))])]
    score_layout!(sub, "random_$(Int(100frac))pct")
end
println("DONE")
