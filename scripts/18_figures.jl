# Figures for the layout-optimization talk (talk/figs/*.pdf).
#
#   TAMBO_MAX_SLOPE=50 julia --project=<env> 18_figures.jl <outdir>
#
# Reads only netscratch outputs (layouts/, search/, analysis/) and the read-only
# official geometry. Palette: reference categorical slots (blue = slope ≤35° runs,
# orange = ≤50° runs, aqua = LLM), gray for reference layouts.
include(joinpath(@__DIR__, "LayoutScore.jl"))
include(joinpath(@__DIR__, "Placer.jl"))
using .LayoutScore, .Placer, TamboSim, AstroParticleUnits, Serialization, Statistics, Printf, CairoMakie, TOML

const W = ARGS[1]
const FD = joinpath(W, "talk", "figs"); mkpath(FD)
const BLUE, ORANGE, AQUA = colorant"#2a78d6", colorant"#eb6834", colorant"#1baf7a"
const INK, INK2, GRID = colorant"#0b0b0b", colorant"#52514e", colorant"#e4e3df"
const REF = colorant"#7a7974"
set_theme!(Theme(fontsize = 15, textcolor = INK,
    Axis = (xgridcolor = GRID, ygridcolor = GRID, leftspinecolor = INK2, bottomspinecolor = INK2,
            topspinevisible = false, rightspinevisible = false, xtickcolor = INK2, ytickcolor = INK2,
            xlabelcolor = INK, ylabelcolor = INK, titlealign = :left, titlefont = :regular),
    Legend = (framevisible = false,)))

# ------------------------------------------------------------------ geometry / layouts
geo = LayoutScore.load_geometry(joinpath(LayoutScore.OFFICIAL, "datasets/v1/volcano_valley/nu_tau/00010"))
d, cs = geo.d_frame, geo.cs
R = Placer.Region(d, cs; max_slope = 50.0)
tris = R.bvh.triangles
tpoly(t) = Point2f[(ustrip(u"m", v.point[1]), ustrip(u"m", v.point[2])) for v in (t.v1, t.v2, t.v3)]
polys = [tpoly(t) for t in tris]
allowed = .!isempty.(R.zone)
xy(units) = [Point2f(ustrip.(u"m", u.center.point)[1:2]...) for u in units]
g = only(filter(f -> f.stream == 'G', collect(geo.tf)))
LAY = Dict(
    "official125" => copy(d["detector_units"]),
    "hex_125" => first(TamboSim.place_detector_units(g, d; spacing = 125u"m")),
    "de50_917" => open(deserialize, joinpath(W, "layouts", "de5s50_917_g015_i27.jls")),
    "de35_917" => open(deserialize, joinpath(W, "layouts", "de3_917_g020_i17.jls")),
    "ls_917" => open(deserialize, joinpath(W, "layouts", "ls_off125_s50.jls")),
    "de50_300" => open(deserialize, joinpath(W, "layouts", "de5s50_300_g017_i38.jls")))
unit_slope(u) = (c = ustrip.(u"m", u.center.point); h = Placer.project(R, c[1], c[2]); isnothing(h) ? NaN : R.slope[h[1]])
zones = TOML.parsefile(joinpath(W, "site_s50", "zones.toml"))["zones"]

# terrain background: slope class per triangle
tcol = [!allowed[i] ? colorant"#f1f0ec" : (R.slope[i] <= 35 ? colorant"#d9e6f5" : colorant"#fbd9c9") for i in eachindex(tris)]
tgray = [!allowed[i] ? colorant"#f4f3f0" : (R.slope[i] <= 35 ? colorant"#e3e2dd" : colorant"#c9c8c2") for i in eachindex(tris)]
function terrain!(ax; labels = false, gray = false)
    poly!(ax, polys; color = gray ? tgray : tcol, strokewidth = 0)
    if labels
        for z in zones
            short = replace(z["name"], "east_facing_" => "E ", "north_facing_" => "N ", "south_facing_" => "S ", "west_facing_" => "W ", "_" => " ")
            text!(ax, z["cx"], z["cy"]; text = short, fontsize = 10, color = INK2, align = (:center, :center))
        end
    end
    ax.aspect = DataAspect(); hidedecorations!(ax); hidespines!(ax)
end

# 1. site map with official125
let f = Figure(size = (700, 620))
    ax = Axis(f[1, 1])
    terrain!(ax; labels = true)
    u = LAY["official125"]; s = unit_slope.(u); P = xy(u)
    scatter!(ax, P[s .<= 35]; color = INK, markersize = 4, label = "official125 unit, slope ≤ 35°")
    scatter!(ax, P[s .> 35]; color = ORANGE, markersize = 5, label = "official125 unit, slope 35–48°")
    Legend(f[2, 1], [PolyElement(color = colorant"#d9e6f5"), PolyElement(color = colorant"#fbd9c9"),
                     MarkerElement(color = INK, marker = :circle, markersize = 8), MarkerElement(color = ORANGE, marker = :circle, markersize = 8)],
           ["terrain ≤ 35° (8.9 km²)", "terrain 35–50° (+4.0 km²)", "official125, ≤ 35° (659)", "official125, > 35° (258)"];
           orientation = :horizontal, nbanks = 2, labelsize = 13)
    save(joinpath(FD, "site_map.pdf"), f)
end

# 2. layout gallery
let f = Figure(size = (900, 640))
    for (k, (key, title, col)) in enumerate([("official125", "official125 (917)", INK), ("hex_125", "hex_125: TamboSim grid, ≤35° (659)", REF),
                                             ("de35_917", "DE, slope ≤ 35° (917)", BLUE), ("de50_917", "DE, slope ≤ 50° (917)", ORANGE),
                                             ("ls_917", "Local search from official125 (917)", AQUA), ("de50_300", "DE, slope ≤ 50° (300)", ORANGE)])
        ax = Axis(f[(k - 1) ÷ 3 + 1, (k - 1) % 3 + 1]; title = title, titlesize = 14)
        terrain!(ax; gray = true)
        scatter!(ax, xy(LAY[key]); color = col, markersize = 5)
    end
    Legend(f[3, 1:3], [PolyElement(color = colorant"#e3e2dd"), PolyElement(color = colorant"#c9c8c2")], ["terrain ≤ 35°", "terrain 35–50°"]; orientation = :horizontal, labelsize = 13)
    save(joinpath(FD, "layouts.pdf"), f)
end

# 3. slope distribution of official125 units
let f = Figure(size = (620, 360))
    s = filter(isfinite, unit_slope.(LAY["official125"]))
    ax = Axis(f[1, 1]; xlabel = "terrain slope at unit (deg)", ylabel = "units")
    hist!(ax, s[s .<= 35]; bins = 0:2.5:50, color = (BLUE, 0.85), strokewidth = 1, strokecolor = :white)
    hist!(ax, s[s .> 35]; bins = 0:2.5:50, color = (ORANGE, 0.9), strokewidth = 1, strokecolor = :white)
    vlines!(ax, [35]; color = INK2, linestyle = :dash, linewidth = 1.5)
    text!(ax, 34, 95; text = "TamboSim default\nmax slope 35°", align = (:right, :top), fontsize = 12, color = INK2)
    text!(ax, 44, 70; text = "258 units\n(28%)", align = (:center, :center), fontsize = 12, color = INK2)
    @printf("official125: %d of %d units above 35°\n", count(>(35), s), length(s))
    save(joinpath(FD, "slope_hist.pdf"), f)
end

# ------------------------------------------------------------------ search logs
readcsv(p) = (L = readlines(p); h = split(L[1], ","); [Dict(zip(h, split(l, ","))) for l in L[2:end]])
num(x) = (v = tryparse(Float64, x); isnothing(v) ? -Inf : v)
"best-so-far reward after each DE generation"
function de_curve(dir)
    p = joinpath(W, "search", dir, "log.csv"); isfile(p) || return nothing
    rows = readcsv(p); gens = sort(unique(parse.(Int, getindex.(rows, "gen"))))
    best = -Inf; out = Float64[]
    for gg in gens
        best = max(best, maximum(num(r["reward"]) for r in rows if parse(Int, r["gen"]) == gg)); push!(out, best)
    end
    return gens, out
end

# 4. DE convergence, ≤35° vs ≤50°
let f = Figure(size = (1000, 400))
    for (k, B) in enumerate((917, 300))
        ax = Axis(f[1, k]; title = "budget $B units", xlabel = "generation (40 plans each)", ylabel = k == 1 ? "best training reward" : "")
        for s in 1:4
            c = de_curve(s == 1 ? "de_b$(B)_s1_forced_raw" : "de_b$(B)_s$(s)_forced_raw_useed"); isnothing(c) && continue
            ok = isfinite.(c[2]); stairs!(ax, c[1][ok], c[2][ok]; step = :post, color = BLUE, linewidth = 2, label = s == 1 ? "DE, slope ≤ 35°" : nothing)
        end
        for s in 5:8
            c = de_curve("de_b$(B)_s$(s)_forced_raw_s50"); isnothing(c) && continue
            ok = isfinite.(c[2]); stairs!(ax, c[1][ok], c[2][ok]; step = :post, color = ORANGE, linewidth = 2, label = s == 5 ? "DE, slope ≤ 50°" : nothing)
        end
        B == 917 && hlines!(ax, [21.21]; color = INK, linestyle = :dash, linewidth = 1.5, label = "official125")
        axislegend(ax; position = :rb, labelsize = 13)
    end
    save(joinpath(FD, "de_convergence.pdf"), f)
end

# 5. local search trajectory
let f = Figure(size = (620, 380))
    ax = Axis(f[1, 1]; xlabel = "local-search generation (20 mutants each)", ylabel = "training reward")
    rows = readcsv(joinpath(W, "search", "ls_off125_s50", "log.csv"))
    gens = sort(unique(parse.(Int, getindex.(rows, "gen"))))
    best = [maximum(num(r["accepted_best"]) for r in rows if parse(Int, r["gen"]) == gg) for gg in gens]
    lines!(ax, gens, best; color = AQUA, linewidth = 2, label = "local search from official125")
    hlines!(ax, [best[1]]; color = INK, linestyle = :dash, linewidth = 1.5, label = "official125 (start)")
    axislegend(ax; position = :rb, labelsize = 13)
    save(joinpath(FD, "ls_trajectory.pdf"), f)
end

# ------------------------------------------------------------------ exact-tier bootstrap tables
function boot(sp)
    D = Dict{String,NTuple{3,Float64}}()
    for l in readlines(joinpath(W, "analysis", "bootstrap_$(sp)_forced_raw.md"))
        startswith(l, "| ") || continue
        c = strip.(split(l, "|")); length(c) < 8 && continue
        v = tryparse(Float64, c[6]); isnothing(v) && continue
        lo, hi = parse.(Float64, split(strip(c[7], ['[', ']']), ","))
        D[c[2]] = (v, lo, hi)
    end
    return D
end
TR, SE = boot("train"), boot("select")
NAMES = [("official125", "official125", INK), ("official150", "official150", REF), ("hex_100", "hex_100", REF), ("hex_125", "hex_125", REF),
         ("plan_de3_917_g020_i17", "DE ≤35° (917)", BLUE), ("plan_de5s50_917_g015_i27", "DE ≤50° (917)", ORANGE),
         ("plan_ls_off125_s50", "local search (917)", AQUA),
         ("plan_de3_300_g014_i18", "DE ≤35° (300)", BLUE), ("plan_de5s50_300_g017_i38", "DE ≤50° (300)", ORANGE),
         ("plan_ls_sub300_s50", "local search (300)", AQUA)]

# 6. forest plot: select split (filled) with train (open) markers
let f = Figure(size = (760, 470))
    sel = filter(n -> haskey(SE, n[1]) && haskey(TR, n[1]), NAMES); ys = length(sel):-1:1
    ax = Axis(f[1, 1]; xlabel = "ν_τ signal rate passing cuts (events / yr, forced channel)", yticks = (collect(ys), [n[2] for n in sel]))
    vlines!(ax, [SE["official125"][1]]; color = GRID, linewidth = 6)
    for (y, (k, _, col)) in zip(ys, sel)
        v, lo, hi = SE[k]; rangebars!(ax, [y], [lo], [hi]; direction = :x, color = col, linewidth = 2, whiskerwidth = 6)
        scatter!(ax, [v], [y]; color = col, markersize = 11, strokecolor = :white, strokewidth = 1.5)
        scatter!(ax, [TR[k][1]], [y + 0.25]; color = :white, strokecolor = col, strokewidth = 1.5, marker = :diamond, markersize = 10)
    end
    Legend(f[2, 1], [MarkerElement(color = INK2, marker = :circle, markersize = 11), MarkerElement(color = :white, strokecolor = INK2, strokewidth = 1.5, marker = :diamond, markersize = 10)],
           ["select split (held out), 95% bootstrap CI", "training split (searched on)"]; orientation = :horizontal, labelsize = 13)
    save(joinpath(FD, "forest.pdf"), f)
end

# 7. generalization: train vs select relative to official125
let f = Figure(size = (560, 500))
    ax = Axis(f[1, 1]; xlabel = "training split, rate / official125", ylabel = "select split, rate / official125")
    lines!(ax, [0.2, 1.4], [0.2, 1.4]; color = GRID, linewidth = 2)
    t0, s0 = TR["official125"][1], SE["official125"][1]
    for (k, lab, col) in NAMES
        (haskey(TR, k) && haskey(SE, k)) || continue
        x, y = TR[k][1] / t0, SE[k][1] / s0
        scatter!(ax, [x], [y]; color = col, markersize = 12, strokecolor = :white, strokewidth = 1.5)
        text!(ax, x, y; text = lab, fontsize = 11, color = INK2, offset = (7, -4), align = (:left, :top))
    end
    xlims!(ax, 0.3, 1.75); ylims!(ax, 0.25, 1.15)
    save(joinpath(FD, "train_vs_select.pdf"), f)
end

# 8. upstream channel: share of MC events vs share of predicted rate
let f = Figure(size = (640, 360))
    sets = ["00010", "00011", "00012", "00013", "00014", "00019"]
    nf = [32028, 31805, 31849, 3097, 3040, 160943]; nu = [2, 4, 186, 52, 130, 0]; rate = [0.90, 0.50, 0.58, 0.12, 0.12, 0.0]
    ax = Axis(f[1, 1]; xticks = (1:6, sets), ylabel = "upstream channel share (%)", xlabel = "ν_τ set (energy ↑ from 00010 to 00014)")
    barplot!(ax, (1:6) .- 0.18, 100 .* nu ./ (nu .+ nf); width = 0.34, color = REF, label = "share of MC events")
    barplot!(ax, (1:6) .+ 0.18, 100 .* rate; width = 0.34, color = ORANGE, label = "share of predicted rate")
    axislegend(ax; position = :rt, labelsize = 13)
    save(joinpath(FD, "upstream.pdf"), f)
end

# 9. LLM loop vs DE (best-so-far vs number of evaluated plans), if available
let f = Figure(size = (1000, 400)), any_llm = false
    for (k, B) in enumerate((917, 300))
        ax = Axis(f[1, k]; title = "budget $B units", xlabel = "plans evaluated", ylabel = k == 1 ? "best training reward" : "")
        for s in 5:8
            c = de_curve("de_b$(B)_s$(s)_forced_raw_s50"); isnothing(c) && continue
            ok = isfinite.(c[2]); stairs!(ax, 40 .* (c[1][ok] .+ 1), c[2][ok]; step = :post, color = ORANGE, linewidth = 2, label = s == 5 ? "DE, slope ≤ 50°" : nothing)
        end
        p = joinpath(W, "search", "llm_b$(B)_glm-5", "log.csv")
        if isfile(p)
            rows = readcsv(p); rs = sort(unique(parse.(Int, getindex.(rows, "round"))))
            best = -Inf; bb = Float64[]
            for r in rs; best = max(best, maximum(num(x["reward"]) for x in rows if parse(Int, x["round"]) == r)); push!(bb, best); end
            scatterlines!(ax, 40 .* (rs .+ 1), bb; color = AQUA, linewidth = 2.5, markersize = 10, label = "LLM (glm-5), best so far")
            scatter!(ax, 40 .* (parse.(Int, getindex.(rows, "round")) .+ 1) .+ 6 .* randn(length(rows)), num.(getindex.(rows, "reward"));
                     color = (AQUA, 0.35), markersize = 5)
            any_llm = true
        end
        B == 917 && hlines!(ax, [21.21]; color = INK, linestyle = :dash, linewidth = 1.5, label = "official125")
        ylims!(ax, B == 917 ? 10 : 4, nothing)
        axislegend(ax; position = :rb, labelsize = 13)
    end
    save(joinpath(FD, "llm_vs_de.pdf"), f)
    println("llm data present: ", any_llm)
end
println("figures -> $FD: ", join(readdir(FD), ", "))

# 10. terrain with official125 / best DE / best LLM, one row (917 units)
let f = Figure(size = (1000, 400))
    llm = open(deserialize, joinpath(W, "layouts", "llm917_r12_ind_05.jls"))
    for (k, (units, title)) in enumerate([(LAY["official125"], "official125"), (LAY["de50_917"], "best DE (≤ 50°)"), (llm, "best LLM")])
        ax = Axis(f[1, k]; title = title, titlesize = 15)
        terrain!(ax)
        s = unit_slope.(units); P = xy(units)
        scatter!(ax, P[s .<= 35]; color = INK, markersize = 4)
        scatter!(ax, P[s .> 35]; color = ORANGE, markersize = 4)
        ax.subtitle = @sprintf("%d units, %.0f%% on > 35°", length(units), 100 * count(>(35), filter(!isnan, s)) / length(units))
    end
    Legend(f[2, 1:3], [PolyElement(color = colorant"#d9e6f5"), PolyElement(color = colorant"#fbd9c9"),
                       MarkerElement(color = INK, marker = :circle, markersize = 8), MarkerElement(color = ORANGE, marker = :circle, markersize = 8)],
           ["terrain ≤ 35°", "terrain 35–50°", "unit on ≤ 35°", "unit on > 35°"]; orientation = :horizontal, labelsize = 13)
    save(joinpath(FD, "terrain_three.pdf"), f); save(joinpath(FD, "terrain_three.png"), f, px_per_unit = 2)
end
