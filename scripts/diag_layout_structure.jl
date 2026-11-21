# What distinguishes official125 from uniform grids and DE plans? Per-zone
# detector counts/densities, spacing statistics, elevation and slope profile.
#
#   julia --project=<env> diag_layout_structure.jl <outdir>
# Writes <outdir>/analysis/layout_structure.md
include(joinpath(@__DIR__, "LayoutScore.jl"))
include(joinpath(@__DIR__, "Placer.jl"))
using .LayoutScore, .Placer, TamboSim, AstroParticleUnits, Serialization, Statistics, Printf

W = ARGS[1]
geo = LayoutScore.load_geometry(joinpath(LayoutScore.OFFICIAL, "datasets/v1/volcano_valley/nu_tau/00010"))
d, cs = geo.d_frame, geo.cs
r = Placer.Region(d, cs)
GCD = joinpath(LayoutScore.OFFICIAL, "gcd", "volcano_valley")
fr150 = collect(load_frames([geo.gc_path, joinpath(GCD, "volcano_valley_d_150m_large.jld2")]))
g = only(filter(f -> f.stream == 'G', collect(geo.tf)))
L = Dict{String,Any}(
    "official125" => copy(d["detector_units"]),
    "official150" => copy(fr150[findlast(f -> f.stream == 'D' && haskey(f.data, "detector_units"), fr150)]["detector_units"]),
    "hex_100" => first(TamboSim.place_detector_units(g, d; spacing = 100u"m")),
    "hex_125" => first(TamboSim.place_detector_units(g, d; spacing = 125u"m")))
for f in ("de3_917_g020_i17", "de3_917_g007_i18", "de3_300_g014_i18")
    L["DE:" * f] = open(deserialize, joinpath(W, "layouts", "$f.jls"))
end
zt = Placer.zone_table(r); zarea = Dict(z.name => z.area_km2 for z in zt); znames = sort([z.name for z in zt])

function profile(units)
    xyz = [ustrip.(u"m", u.center.point) for u in units]
    zone = String[]; slope = Float64[]
    for c in xyz
        h = Placer.project(r, c[1], c[2])
        push!(zone, isnothing(h) ? "?" : (isempty(r.zone[h[1]]) ? "steep" : r.zone[h[1]])); push!(slope, isnothing(h) ? NaN : r.slope[h[1]])
    end
    P = reduce(hcat, xyz)
    nn = [minimum(sqrt(sum(abs2, P[1:2, i] .- P[1:2, j])) for j in axes(P, 2) if j != i) for i in axes(P, 2)]
    return (; n = length(units), zone, slope, z = P[3, :], nn)
end
prof = Dict(k => profile(v) for (k, v) in L)
order = ["official125", "official150", "hex_100", "hex_125", "DE:de3_917_g020_i17", "DE:de3_917_g007_i18", "DE:de3_300_g014_i18"]

open(joinpath(W, "analysis", "layout_structure.md"), "w") do io
    println(io, "# Layout structure\n")
    println(io, "| layout | units | NN dist median (m) | NN p10–p90 (m) | elevation median (m) | elev p10–p90 | slope median (deg) | off-zone (steep/none) |")
    println(io, "|---|---|---|---|---|---|---|---|")
    for k in order
        p = prof[k]
        @printf(io, "| %s | %d | %.0f | %.0f–%.0f | %.0f | %.0f–%.0f | %.1f | %d |\n", k, p.n, median(p.nn), quantile(p.nn, 0.1), quantile(p.nn, 0.9),
                median(p.z), quantile(p.z, 0.1), quantile(p.z, 0.9), median(filter(isfinite, p.slope)), count(z -> z in ("steep", "?"), p.zone))
    end
    println(io, "\n## Detectors per zone (density per km² in brackets)\n")
    println(io, "| zone | km² | " * join(order, " | ") * " |"); println(io, "|---|---|" * repeat("---|", length(order)))
    for zn in znames
        cells = [@sprintf("%d (%.0f)", count(==(zn), prof[k].zone), count(==(zn), prof[k].zone) / zarea[zn]) for k in order]
        @printf(io, "| %s | %.2f | %s |\n", zn, zarea[zn], join(cells, " | "))
    end
    # elevation profile: fraction of units per 250 m elevation band
    println(io, "\n## Fraction of detectors per elevation band\n")
    edges = 2400:250:3900
    println(io, "| band (m) | " * join(order, " | ") * " |"); println(io, "|---|" * repeat("---|", length(order)))
    for i in 1:length(edges)-1
        cells = [@sprintf("%.2f", count(z -> edges[i] <= z < edges[i+1], prof[k].z) / prof[k].n) for k in order]
        println(io, "| $(edges[i])–$(edges[i+1]) | " * join(cells, " | ") * " |")
    end
end
print(read(joinpath(W, "analysis", "layout_structure.md"), String))
