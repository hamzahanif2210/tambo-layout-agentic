# Placer: named mountain zones + deterministic plan -> detector units (OBBs).
#
# Region  = the D frame's detector_bvh triangles (2161 tris, 12.9 km²) — the same
#           surface CORSIKA records particles on; nothing outside it can be hit.
# Zones   = region triangles with slope <= max_slope, grouped by the direction
#           the slope faces (4 quadrants) x elevation tercile (low/mid/high).
# Plan    = TOML (see example_plan.toml):
#   budget = 300
#   [rules]  min_spacing_m = 60, max_slope_deg = 35, seed = 1
#   [[zones]] name = "east_facing_mid", n = 120, pattern = "tri_grid",
#             spacing_m = 110, orient_deg = 20
# Patterns: tri_grid | contour | cluster | random
#   tri_grid: triangular lattice (spacing, orientation), the n points nearest the zone centre
#   contour : rows along elevation contours every row_step_m (vertical, default 40), along-row step = spacing
#   cluster : Gaussian around zone centre (radius_m), thinned to min spacing
#   random  : uniform over the zone, thinned to min spacing
#   ring    : n_rings concentric horizontal circles (radius_m, ring_step_m) about the zone
#             centre or [center_x, center_y]; points every spacing_m along each circle
# Zone name "all" = every placeable triangle (lets one pattern span zones).
# Detector boxes are built exactly as TamboSim.place_detector_units does
# (vertical ray onto the region, OBB aligned to the triangle normal, 2x2x0.25 m).
module Placer

using TamboSim, AstroParticleUnits, TOML, Random, Statistics, LinearAlgebra, Printf

const HALF = [1.0u"m", 1.0u"m", 0.125u"m"]
const FACES = ("north_facing", "east_facing", "south_facing", "west_facing")
const BANDS = ("low", "mid", "high")

struct Region
    bvh::Any
    cs::Any
    cx::Vector{Float64}; cy::Vector{Float64}; cz::Vector{Float64}
    area::Vector{Float64}; slope::Vector{Float64}; aspect::Vector{Float64}   # deg; aspect = azimuth (from +y, clockwise) the slope faces
    zone::Vector{String}                                                      # "" if too steep
end

# max_slope default from env TAMBO_MAX_SLOPE (35°). Official arrays place units up to ~48°,
# so steep-aware runs use 50 (2026-11-19).
function Region(d_frame, cs; max_slope = parse(Float64, get(ENV, "TAMBO_MAX_SLOPE", "35")))
    bvh = d_frame["detector_bvh"]; tris = bvh.triangles
    c(t, i) = ustrip(u"m", (t.v1.point[i] + t.v2.point[i] + t.v3.point[i]) / 3)
    cx = [c(t, 1) for t in tris]; cy = [c(t, 2) for t in tris]; cz = [c(t, 3) for t in tris]
    area = [ustrip(u"m^2", TamboSim.area(t)) for t in tris]
    n = [TamboSim.normal(t).point for t in tris]
    n = [v[3] < 0 ? -v : v for v in n]                       # upward normal
    slope = [rad2deg(acos(clamp(v[3], -1, 1))) for v in n]
    aspect = [mod(rad2deg(atan(v[1], v[2])), 360) for v in n] # horizontal part of the normal points downhill
    ok = slope .<= max_slope
    zq = quantile(cz[ok], [1/3, 2/3])
    zone = map(eachindex(tris)) do i
        ok[i] || return ""
        f = FACES[mod(round(Int, aspect[i] / 90), 4) + 1]
        b = cz[i] < zq[1] ? "low" : cz[i] < zq[2] ? "mid" : "high"
        "$(f)_$(b)"
    end
    # The site is one dominant (east-facing) slope; split any zone > 1 km² at its
    # area-weighted median y into _north / _south halves for finer plan control.
    for z in unique(filter(!isempty, zone))
        idx = findall(==(z), zone)
        sum(area[idx]) > 1e6 || continue
        o = sortperm(cy[idx]); cum = cumsum(area[idx][o]); ymed = cy[idx][o][findfirst(>=(cum[end] / 2), cum)]
        for i in idx
            zone[i] = cy[i] < ymed ? "$(z)_south" : "$(z)_north"
        end
    end
    return Region(bvh, cs, cx, cy, cz, area, slope, aspect, zone)
end

zone_names(r::Region) = sort(unique(filter(!isempty, r.zone)))

"""Per-zone summary (the TAMBO analogue of the PCB project's board_context)."""
function zone_table(r::Region)
    map(zone_names(r)) do z
        idx = findall(==(z), r.zone); a = r.area[idx]
        (; name = z, area_km2 = sum(a) / 1e6, n_tri = length(idx),
           cx = sum(r.cx[idx] .* a) / sum(a), cy = sum(r.cy[idx] .* a) / sum(a),
           zmin = minimum(r.cz[idx]), zmax = maximum(r.cz[idx]),
           mean_slope = sum(r.slope[idx] .* a) / sum(a))
    end
end

function site_context(r::Region)
    io = IOBuffer()
    tot = sum(r.area) / 1e6; okA = sum(r.area[.!isempty.(r.zone)]) / 1e6
    println(io, "# TAMBO site context (volcano_valley)\n")
    @printf(io, "Detector region: %d mesh triangles, %.2f km^2; %.2f km^2 at slope <= 35 deg (placeable).\n",
            length(r.cx), tot, okA)
    @printf(io, "Coordinates: local x, y (m, horizontal), z = elevation %.0f-%.0f m.\n", extrema(r.cz)...)
    println(io, "Face names assume local +x = east, +y = north (unverified; face = direction the slope looks toward).")
    println(io, "Detector unit: 2 m x 2 m x 0.25 m slab on the surface. Shower particles are recorded on this surface only.\n")
    println(io, "| zone | area km^2 | centre x,y (m) | elevation (m) | mean slope (deg) |")
    println(io, "|---|---|---|---|---|")
    for z in zone_table(r)
        @printf(io, "| %s | %.2f | %.0f, %.0f | %.0f-%.0f | %.0f |\n", z.name, z.area_km2, z.cx, z.cy, z.zmin, z.zmax, z.mean_slope)
    end
    return String(take!(io))
end

# ----------------------------------------------------------- projection ---
const UP = Ref{Any}(nothing)
"""Vertical ray at (x, y) onto the region. Returns (tri_index, point, normal) or nothing."""
function project(r::Region, x, y)
    up = TamboSim.Direction([0.0, 0.0, 1.0], r.cs)
    ray = TamboSim.Ray(TamboSim.Coordinate(x * u"m", y * u"m", 0.0u"m", r.cs), up)
    ixs = TamboSim.intersect_all(r.bvh, ray)
    isempty(ixs) && return nothing
    ix = ixs[1]
    return (ix.index, ix.point, ix.normal)
end

"""OBB exactly as TamboSim.place_detector_units builds it."""
function make_obb(r::Region, point, nrm)
    up = TamboSim.Direction([0.0, 0.0, 1.0], r.cs)
    ψ = acos(clamp(dot(nrm, up), -1.0, 1.0))
    rot = ψ < eps(Float64) ? TamboSim.AngleAxis(0.0, 1.0, 0.0, 0.0) :
                             TamboSim.AngleAxis(ψ, cross(up, nrm)...)
    return TamboSim.OBB(point, rot, HALF)
end

# ------------------------------------------------------------- patterns ---
function lattice(cx, cy, spacing, orient_deg, R)
    θ = deg2rad(orient_deg); e1 = (cos(θ), sin(θ)); e2 = (cos(θ + π/3), sin(θ + π/3))
    k = ceil(Int, R / spacing) + 1
    pts = Tuple{Float64,Float64}[]
    for i in -k:k, j in -k:k
        push!(pts, (cx + spacing * (i * e1[1] + j * e2[1]), cy + spacing * (i * e1[2] + j * e2[2])))
    end
    return pts
end

zone_radius(r, idx, cx, cy) = maximum(hypot.(r.cx[idx] .- cx, r.cy[idx] .- cy)) + 200

"""Candidate (x,y) points for one zone entry, ordered by preference."""
function candidates(r::Region, z, entry, rng)
    idx = z.name == "all" ? findall(!isempty, r.zone) : findall(==(z.name), r.zone)
    pat = get(entry, "pattern", "tri_grid"); sp = Float64(get(entry, "spacing_m", 100.0))
    if pat == "tri_grid"
        pts = lattice(z.cx, z.cy, sp, Float64(get(entry, "orient_deg", 0.0)), zone_radius(r, idx, z.cx, z.cy))
        sort!(pts; by = p -> hypot(p[1] - z.cx, p[2] - z.cy))
    elseif pat == "contour"
        rs = Float64(get(entry, "row_step_m", 40.0))                                 # vertical step between rows
        pts = lattice(z.cx, z.cy, min(sp, rs) / 3, 0.0, zone_radius(r, idx, z.cx, z.cy))   # dense sampling, filtered below
        sort!(pts; by = p -> hypot(p[1] - z.cx, p[2] - z.cy))
    elseif pat == "cluster"
        rad = Float64(get(entry, "radius_m", 300.0))
        pts = [(z.cx + rad * randn(rng), z.cy + rad * randn(rng)) for _ in 1:20 * entry["n"]]
    elseif pat == "ring"
        c0x = Float64(get(entry, "center_x", z.cx)); c0y = Float64(get(entry, "center_y", z.cy))
        R0 = Float64(get(entry, "radius_m", 500.0)); dR = Float64(get(entry, "ring_step_m", 150.0))
        pts = Tuple{Float64,Float64}[]
        for k in 0:Int(get(entry, "n_rings", 1)) - 1
            R = R0 + k * dR; m = max(3, round(Int, 2π * R / sp))
            ϕ0 = deg2rad(Float64(get(entry, "orient_deg", 0.0)))
            append!(pts, [(c0x + R * cos(ϕ0 + 2π * j / m), c0y + R * sin(ϕ0 + 2π * j / m)) for j in 0:m-1])
        end
    elseif pat == "random"
        R = zone_radius(r, idx, z.cx, z.cy)
        pts = [(z.cx + R * (2rand(rng) - 1), z.cy + R * (2rand(rng) - 1)) for _ in 1:50 * entry["n"]]
    else
        error("unknown pattern $pat")
    end
    return pts, pat, sp
end

"""
    place(r, plan::Dict) -> (units::Vector{OBB}, report::NamedTuple)

Deterministic for a given plan (incl. rules.seed). Never fills a zone's shortfall
from elsewhere: shortfalls and budget mismatches are reported (the reward
penalises them) rather than hidden.
"""
function place(r::Region, plan::Dict)
    rules = get(plan, "rules", Dict{String,Any}())
    minsp = Float64(get(rules, "min_spacing_m", 50.0))
    rng = MersenneTwister(Int(get(rules, "seed", 1)))
    zt = Dict{String,Any}(z.name => z for z in zone_table(r))
    okA = r.area .* .!isempty.(r.zone)
    zt["all"] = (; name = "all", cx = sum(r.cx .* okA) / sum(okA), cy = sum(r.cy .* okA) / sum(okA))
    units = TamboSim.OBB{Float64}[]; placed_xyz = NTuple{3,Float64}[]
    zone_report = []
    for entry in get(plan, "zones", Any[])
        name = entry["name"]; n = Int(entry["n"])
        haskey(zt, name) || (push!(zone_report, (; name, requested = n, placed = 0, filled = 0, note = "unknown zone")); continue)
        pts, pat, sp = candidates(r, zt[name], entry, rng)
        got = 0; rows = Dict{Int,Vector{NTuple{3,Float64}}}()
        for (x, y) in pts
            got >= n && break
            hit = project(r, x, y); isnothing(hit) && continue
            ti, p, nrm = hit
            (name == "all" ? !isempty(r.zone[ti]) : r.zone[ti] == name) || continue   # stay inside the zone (also enforces slope)
            P = ustrip.(u"m", p.point); q = (P[1], P[2], P[3])
            if pat == "contour"                                  # keep points near elevation rows
                rs = Float64(get(entry, "row_step_m", 40.0))
                row = round(Int, q[3] / rs)
                abs(q[3] - row * rs) <= rs / 6 || continue
                any(o -> hypot(o[1] - q[1], o[2] - q[2]) < sp, get(rows, row, [])) && continue
                push!(get!(rows, row, NTuple{3,Float64}[]), q)
            end
            any(o -> norm((o[1] - q[1], o[2] - q[2], o[3] - q[3])) < minsp, placed_xyz) && continue
            push!(units, make_obb(r, p, nrm)); push!(placed_xyz, q); got += 1
        end
        # Fill step (rules.fill_shortfall, default true): if the pattern ran out of
        # positions, top the zone up with a dense lattice at ~min spacing, nearest
        # the zone centre first. Reported as `filled` so it is never hidden.
        filled = 0
        if got < n && get(rules, "fill_shortfall", true)
            z = zt[name]
            idxz = name == "all" ? findall(!isempty, r.zone) : findall(==(name), r.zone)
            fill_pts = lattice(z.cx, z.cy, 1.05minsp, 0.0, zone_radius(r, idxz, z.cx, z.cy))
            sort!(fill_pts; by = p -> hypot(p[1] - z.cx, p[2] - z.cy))
            for (x, y) in fill_pts
                got >= n && break
                hit = project(r, x, y); isnothing(hit) && continue
                ti, p, nrm = hit
                (name == "all" ? !isempty(r.zone[ti]) : r.zone[ti] == name) || continue
                P = ustrip.(u"m", p.point); q = (P[1], P[2], P[3])
                any(o -> norm((o[1] - q[1], o[2] - q[2], o[3] - q[3])) < minsp, placed_xyz) && continue
                push!(units, make_obb(r, p, nrm)); push!(placed_xyz, q); got += 1; filled += 1
            end
        end
        push!(zone_report, (; name, requested = n, placed = got, filled,
                              note = got < n ? "shortfall" : filled > 0 ? "filled $filled" : ""))
    end
    budget = Int(get(plan, "budget", length(units)))
    report = (; budget, placed = length(units), zones = zone_report,
              shortfall = sum(z -> z.requested - z.placed, zone_report; init = 0),
              over_budget = max(0, sum(z -> z.requested, zone_report; init = 0) - budget))
    return units, report
end

place(r::Region, path::AbstractString) = place(r, TOML.parsefile(path))

end # module
