# LayoutScore: re-score arbitrary TAMBO detector layouts from existing CORSIKA
# surface particles (TamboSim ObservationMesh output).
#
# READ-ONLY w.r.t. official_simulations. Libraries are written only to the
# netscratch work dir passed in by the caller.
#
# Pieces:
#   build_library(simset, run; maxev)  -> NamedTuple of plain column arrays
#                                         (selected particles only + event truth
#                                          + stored official selected hits)
#   save_library / load_library        -> Serialization (internal cache format)
#   GridIndex(lib; cell)               -> horizontal (cs x,y) spatial hash
#   map_hits_fast(lib, grid, units, bvh, topo_bvh; R) -> per-particle hits
#   map_hits_exact(lib, units, bvh, topo_bvh)          -> brute force, reference
#   event_tables(lib, hits, units)     -> LlhamaFit-compatible hits table per event
#   score(...)                         -> trigger + fit summary
module LayoutScore

using TamboSim, LlhamaFit, AstroParticleUnits, TOML, Serialization, Statistics

const OFFICIAL = "/n/holylfs05/LABS/arguelles_delgado_lab/Lab/TAMBO/official_simulations"
const SELECT = LlhamaFit.DEFAULT_SELECTION          # e±/μ± with E >= 3 MeV

# ---------------------------------------------------------------- geometry ---
"""Load the shared G/C + official D frames (read-only). Returns (tf, d_frame, cs, topo_bvh)."""
function load_geometry(simset)
    gc_path = TOML.parsefile(joinpath(simset, "injection.toml"))["geometry"]["gc_path"]
    d_path  = TOML.parsefile(joinpath(simset, "hits_125m_large.toml"))["detector"]["d_path"]
    tf = load_frames([gc_path, d_path])
    # The GC bundle carries its own (unit-less) D frame; the array file adds the
    # D frame with detector_units. Take the last D that has units — the same
    # one Q frames resolve to via q.d_frame. `cs` / `bvh` resolve up its parents.
    fr = collect(tf)
    d = fr[findlast(f -> f.stream == 'D' && haskey(f.data, "detector_units"), fr)]
    return (; tf, d_frame = d, cs = d["cs"], topo_bvh = d["bvh"], gc_path, d_path)
end

# ----------------------------------------------------------------- library ---
"""
Read one official run once and keep only what scoring needs: selected particles
(position, direction, time, speed, energy, pdg, weight, event index) as Float64
columns, event truth, and the officially stored selected hits per event (for
validation). Streams event by event so peak memory is one event's particles.
"""
function build_library(simset, run; maxev = typemax(Int), evrange = nothing)
    # Event source: the hits file when it exists (carries official hits for
    # validation), else the planned-showers file (same Q frames + corsika dirs).
    hit_file = joinpath(simset, "hits_125m_large", "$run.jld2")
    src = isfile(hit_file) ? hit_file : joinpath(simset, "showers", "planned_showers_$run.jld2")
    geo = load_geometry(simset)
    tf = load_frames([geo.gc_path, geo.d_path, src])
    nall = length(tf.q_frames)
    rng = isnothing(evrange) ? (1:min(nall, maxev)) : (first(evrange):min(last(evrange), nall))
    qs = tf.q_frames[rng]
    cs = geo.cs
    injection_cfg = deepcopy(tf.m_frames[1]["injection"])    # emin/emax/gamma/nevent/angles/xs
    psp = Vector{Any}(undef, length(qs))                     # per-event phase_space_point

    ev = Int32[]; pdg = Int32[]; E = Float64[]; t = Float64[]; v = Float64[]; w = Float64[]
    x = Float64[]; y = Float64[]; z = Float64[]; dx = Float64[]; dy = Float64[]; dz = Float64[]
    nev = length(qs)
    truthE = fill(NaN, nev); truthθ = fill(NaN, nev); truthϕ = fill(NaN, nev)
    event_id = Vector{Any}(undef, nev)
    stored = Vector{Vector{Int}}(undef, nev)      # official selected-hit module indices
    n_raw = 0

    for (i, q) in enumerate(qs)
        inj = q["injection_initial_state"]
        truthE[i] = ustrip(u"GeV", inj.energy)
        truthθ[i] = Float64(inj.direction.θ); truthϕ[i] = Float64(inj.direction.ϕ)
        event_id[i] = get(q.data, "event_id", i)
        psp[i] = get(q.data, "phase_space_point", nothing)
        stored[i] = haskey(q.data, "corsika_hits") ?
            sort!([h.module_index for h in q["corsika_hits"] if SELECT(h)]) : Int[]
        haskey(q, "corsika_directories") || continue
        for d in unique(dirname.(q["corsika_directories"]))
            isdir(d) || continue
            evs = try TamboSim.read_corsika_mesh(d, cs) catch err
                @warn "skip" dir = d exception = err; continue end
            for e in evs
                n_raw += 1
                p = e.particle
                SELECT((; particle = p)) || continue
                push!(ev, i); push!(pdg, Int32(Int(p.pdg)))
                push!(E, ustrip(u"GeV", p.energy)); push!(t, ustrip(u"s", p.time))
                push!(v, ustrip(u"m/s", p.speed)); push!(w, e.weight)
                P = p.position.point; D = p.direction.point
                push!(x, ustrip(u"m", P[1])); push!(y, ustrip(u"m", P[2])); push!(z, ustrip(u"m", P[3]))
                push!(dx, D[1]); push!(dy, D[2]); push!(dz, D[3])
            end
        end
    end
    return (; simset, run, nev, n_raw, ev, pdg, E, t, v, w, x, y, z, dx, dy, dz,
            truthE, truthθ, truthϕ, event_id, stored,
            source = src, n_events_in_run = nall, evrange = rng,
            has_stored = isfile(hit_file), injection_cfg, psp)
end

# Atomic: a job preempted mid-write (serial_requeue) leaves only a .tmp, never a
# truncated library that a requeued job would then try to load.
function save_library(path, lib)
    mkpath(dirname(path))
    tmp = path * ".tmp.$(getpid())"
    open(io -> serialize(io, lib), tmp, "w")
    mv(tmp, path; force = true)
end
load_library(path) = open(deserialize, path)

"""Rebuild the exact TamboSim Particle for library row k (all Float64, same cs)."""
@inline function particle(lib, k, cs)
    pos = TamboSim.Coordinate([lib.x[k] * u"m", lib.y[k] * u"m", lib.z[k] * u"m"], cs)
    dir = TamboSim.Direction([lib.dx[k], lib.dy[k], lib.dz[k]], cs)
    return TamboSim.Particle(TamboSim.ParticleType(lib.pdg[k]), lib.E[k] * u"GeV", pos, dir,
                             lib.t[k] * u"s", lib.v[k] * u"m/s")
end

@inline function hit_one(lib, k, cs, bvh, topo_bvh)
    ce = TamboSim.CorsikaEvent(particle(lib, k, cs), lib.w[k])
    return TamboSim._corsika_module_hit(ce, bvh, topo_bvh)
end

# -------------------------------------------------------------------- grid ---
struct GridIndex
    cell::Float64
    keys::Dict{Tuple{Int,Int},UnitRange{Int}}   # cell -> range into `order`
    order::Vector{Int32}                         # particle rows sorted by cell
end

function GridIndex(lib; cell = 4.0)
    n = length(lib.x)
    ck = Vector{Tuple{Int,Int}}(undef, n)
    @inbounds for k in 1:n
        ck[k] = (floor(Int, lib.x[k] / cell), floor(Int, lib.y[k] / cell))
    end
    order = Int32.(sortperm(ck))
    keys = Dict{Tuple{Int,Int},UnitRange{Int}}()
    s = 1
    @inbounds for j in 2:n+1
        if j > n || ck[order[j]] != ck[order[s]]
            keys[ck[order[s]]] = s:j-1; s = j
        end
    end
    return GridIndex(cell, keys, order)
end

const HitRec = NamedTuple{(:row, :module_index, :time, :weight),Tuple{Int32,Int,Float64,Float64}}

"""Candidates = particles whose recorded surface crossing lies within R (horizontal) of a module centre; then the exact TamboSim hit test."""
function map_hits_fast(lib, grid::GridIndex, units, bvh, cs, topo_bvh; R = 10.0)
    cand = Set{Int32}()
    r = ceil(Int, R / grid.cell)
    for u in units
        c = u.center.point
        cx = floor(Int, ustrip(u"m", c[1]) / grid.cell); cy = floor(Int, ustrip(u"m", c[2]) / grid.cell)
        for ix in cx-r:cx+r, iy in cy-r:cy+r
            rg = get(grid.keys, (ix, iy), nothing)
            isnothing(rg) && continue
            for j in rg
                push!(cand, grid.order[j])
            end
        end
    end
    hits = HitRec[]
    for k in cand
        h = hit_one(lib, k, cs, bvh, topo_bvh)
        isnothing(h) || push!(hits, (row = k, module_index = h.module_index,
                                     time = ustrip(u"s", h.particle.time), weight = h.weight))
    end
    return hits, length(cand)
end

"""Brute force over every library particle (reference)."""
function map_hits_exact(lib, units, bvh, cs, topo_bvh)
    hits = HitRec[]
    for k in eachindex(lib.x)
        h = hit_one(lib, Int32(k), cs, bvh, topo_bvh)
        isnothing(h) || push!(hits, (row = Int32(k), module_index = h.module_index,
                                     time = ustrip(u"s", h.particle.time), weight = h.weight))
    end
    return hits
end

"""Per event: sorted unique hit module indices (for validation against `lib.stored`)."""
function event_modules(lib, hits)
    out = [Int[] for _ in 1:lib.nev]
    for h in hits
        push!(out[lib.ev[h.row]], h.module_index)
    end
    return [sort!(m) for m in out]
end

"""LlhamaFit-compatible hits tables (same columns/semantics as LlhamaFit.hits_table)."""
function event_tables(lib, hits, units)
    acc = [Dict{Int,Tuple{Float64,Float64}}() for _ in 1:lib.nev]
    for h in hits
        a = acc[lib.ev[h.row]]
        prev = get(a, h.module_index, nothing)
        a[h.module_index] = isnothing(prev) ? (h.weight, h.time) :
                            (prev[1] + h.weight, min(prev[2], h.time))
    end
    return map(acc) do a
        idx = sort!(collect(keys(a)))
        LlhamaFit.StructArray(
            module_index = idx,
            module_pos   = [units[i].center.point for i in idx],
            hit_time     = [uconvert(μs, a[i][2] * u"s") for i in idx],
            nhits        = [a[i][1] for i in idx])
    end
end

"""Trigger + PlaneTiming fit summary for one layout."""
function score(lib, tabs, geom; trigger = LlhamaFit.TRIGGER_5x1, method = LlhamaFit.PlaneTiming())
    trig = [trigger(h) for h in tabs]
    dψ = Float64[]; ids = Int[]; rθ = Float64[]; rϕ = Float64[]; sψ = Float64[]; u3 = Float64[]; nmod = Int[]
    for i in findall(trig)
        r = try LlhamaFit.reconstruct(tabs[i], geom, method) catch; nothing end
        (isnothing(r) || !r.valid) && continue
        push!(dψ, Float64(LlhamaFit.great_circle_distance(lib.truthθ[i], lib.truthϕ[i], r.θ, r.ϕ) / °))
        push!(ids, i); push!(rθ, Float64(r.θ)); push!(rϕ, Float64(r.ϕ))
        push!(sψ, try Float64(r.σψ / °) catch; NaN end)
        push!(u3, try Float64(r.u3) catch; NaN end); push!(nmod, length(tabs[i]))
    end
    return (; ntrig = count(trig), nfit = length(dψ),
            median_dpsi = isempty(dψ) ? NaN : median(dψ), dψ, ids, trig, rθ, rϕ, sψ, u3, nmod)
end

# ----------------------------------------------------------------- weights ---
# Mirrors TamboSim.build_phase_space(m) (weighting/phase_space.jl) but from the
# injection config dict stored in the library, so no M Frame is needed.
function phase_space(cfg)
    args = (UInt(cfg["geometry_hash"]), Int(cfg["pdg"]),
            Float64(cfg["emin"]) * u"GeV", Float64(cfg["emax"]) * u"GeV", Float64(cfg["gamma"]),
            deg2rad(Float64(cfg["thetamin"])), deg2rad(Float64(cfg["thetamax"])),
            deg2rad(Float64(cfg["phimin"])), deg2rad(Float64(cfg["phimax"])), Int(cfg["nevent"]))
    s = cfg["strategy"]
    s in ("NeutrinoInjection", "MuonNeutrinoInjection") && return TamboSim.NeutrinoInjectionPS(args...)
    if s == "CosmicRayInjection"
        # TamboSim v1.0.1 (pinned by LlhamaFit) has no `altitude` field; newer
        # versions add it (used only for campaign-compatibility, not the pdf).
        fieldcount(TamboSim.CosmicRayInjectionPS) == length(args) &&
            return TamboSim.CosmicRayInjectionPS(args...)
        return TamboSim.CosmicRayInjectionPS(args..., Float64(get(cfg, "altitude", 112.0)) * u"km")
    end
    error("unknown injection strategy $s")
end

"""
Per-event OneWeight in GeV·m²·sr for a SINGLE run (TamboSim convention:
1 / (pdf(pt) · nevent)). When pooling N runs of the same campaign, divide by N
(TamboSim's oneweight sums pdf·nevent over every compatible M frame).
Rate for a flux Φ [GeV⁻¹ m⁻² s⁻¹ sr⁻¹]: Σ_triggered Φ(E_i) · OW_i / N_runs.
"""
function oneweights(lib)
    ps = phase_space(lib.injection_cfg)
    return map(lib.psp) do pt
        isnothing(pt) && return 0.0
        dens = ps(pt).val                      # GeV⁻¹ m⁻² sr⁻¹ (functor uconverts)
        dens == 0 ? 0.0 : 1 / (dens * ps.nevent)
    end
end

# Collaboration reference fluxes (TamboMakie.jl src/references/flux_models.jl),
# both returned in GeV⁻¹ m⁻² s⁻¹ sr⁻¹ with E in GeV.
"Per-flavour astrophysical ν+ν̄: E^-2.52, 1.8e-18 GeV⁻¹cm⁻²s⁻¹sr⁻¹ at 100 TeV."
flux_nu(E) = 1.8e-18 * 1e4 * (E / 1e5)^(-2.52)
"All-particle cosmic rays: E^-2.7, 1.8e4 GeV⁻¹m⁻²s⁻¹sr⁻¹ at 1 GeV."
flux_cr(E) = 1.8e4 * E^(-2.7)
const SECONDS_PER_YEAR = 3.15576e7

"""
Keep only events whose split tag is in `tags` (tag = run mod 5 for multi-run
sets, event_index mod 5 for single-run sets — same rule as Objective.jl).
Returns a library with the same fields, events renumbered 1..k.
"""
function subset_events(lib, tags; by_run::Bool)
    evidx = [first(lib.evrange) + i - 1 for i in 1:lib.nev]
    keep = by_run ? (mod(parse(Int, lib.run), 5) in tags ? collect(1:lib.nev) : Int[]) :
                    [i for i in 1:lib.nev if mod(evidx[i], 5) in tags]
    newid = zeros(Int32, lib.nev); newid[keep] .= 1:length(keep)
    rows = findall(k -> newid[lib.ev[k]] > 0, eachindex(lib.ev))
    pcol(c) = c[rows]
    return merge(lib, (; nev = length(keep), ev = newid[lib.ev[rows]],
        pdg = pcol(lib.pdg), E = pcol(lib.E), t = pcol(lib.t), v = pcol(lib.v), w = pcol(lib.w),
        x = pcol(lib.x), y = pcol(lib.y), z = pcol(lib.z), dx = pcol(lib.dx), dy = pcol(lib.dy), dz = pcol(lib.dz),
        truthE = lib.truthE[keep], truthθ = lib.truthθ[keep], truthϕ = lib.truthϕ[keep],
        event_id = lib.event_id[keep], stored = lib.stored[keep], psp = lib.psp[keep],
        event_index = evidx[keep]))
end

end # module
