# Milestone 2: build a compact particle library for one official run and prove
# the fast (grid-indexed) mapper is exact.
#
#   julia --project=<LlhamaFit examples env> 02_build_and_validate.jl <simset_dir> <run> <maxev> <outdir>
#
# Reads official_simulations (read-only). Writes only <outdir>/library/....

include(joinpath(@__DIR__, "LayoutScore.jl"))
using .LayoutScore, TamboSim, LlhamaFit, AstroParticleUnits, Random, Statistics, Printf

simset = ARGS[1]; run = ARGS[2]; maxev = parse(Int, ARGS[3]); outdir = ARGS[4]
startswith(abspath(outdir), LayoutScore.OFFICIAL) && error("refusing to write inside official_simulations")
tag = joinpath(basename(dirname(simset)), basename(simset))           # e.g. nu_tau/00010
libpath = joinpath(outdir, "library", tag, "$(run)_max$(maxev).jls")

println("simset=$simset run=$run maxev=$maxev\nlibrary=$libpath")

# ---- build (or reuse) the library
lib = if isfile(libpath)
    t = @elapsed l = LayoutScore.load_library(libpath); @printf("load_library_s=%.1f\n", t); l
else
    t = @elapsed l = LayoutScore.build_library(simset, run; maxev)
    @printf("build_library_s=%.1f\n", t)
    ts = @elapsed LayoutScore.save_library(libpath, l); @printf("save_library_s=%.1f\n", ts); l
end
@printf("events=%d raw_particles=%d selected_particles=%d (%.1f%%) file_MB=%.0f\n",
        lib.nev, lib.n_raw, length(lib.x), 100 * length(lib.x) / max(lib.n_raw, 1), filesize(libpath) / 1e6)

geo = LayoutScore.load_geometry(simset)
d, cs, topo = geo.d_frame, geo.cs, geo.topo_bvh
units0 = copy(d["detector_units"])

tg = @elapsed grid = LayoutScore.GridIndex(lib; cell = 4.0)
@printf("grid_build_s=%.1f cells=%d\n", tg, length(grid.keys))

function set_layout!(units)
    bvh = TamboSim.BVHTree(units)
    d["detector_units"] = units; d["detector_unit_bvh"] = bvh
    return bvh, LlhamaFit.DetectorGeometry(d)
end
nmismatch(a, b) = count(i -> a[i] != b[i], eachindex(a))

# ---- 1. official layout: exact vs stored, fast(R) vs stored
bvh, geom = set_layout!(units0)
te = @elapsed hx = LayoutScore.map_hits_exact(lib, units0, bvh, cs, topo)
mx = LayoutScore.event_modules(lib, hx)
@printf("[official exact] map_s=%.1f hits=%d mismatch_vs_stored=%d/%d\n",
        te, length(hx), nmismatch(mx, lib.stored), lib.nev)

# distance (horizontal) from recorded surface crossing to the struck module centre
dist = [hypot(lib.x[h.row] - ustrip(u"m", units0[h.module_index].center.point[1]),
              lib.y[h.row] - ustrip(u"m", units0[h.module_index].center.point[2])) for h in hx]
isempty(dist) || @printf("hit_distance_m: median=%.2f p99=%.2f p99.99=%.2f max=%.2f\n",
                         median(dist), quantile(dist, 0.99), quantile(dist, 0.9999), maximum(dist))

for R in (3.0, 5.0, 10.0)
    tf_ = @elapsed (hf, ncand) = LayoutScore.map_hits_fast(lib, grid, units0, bvh, cs, topo; R)
    mf = LayoutScore.event_modules(lib, hf)
    @printf("[official fast R=%4.1f] map_s=%.2f candidates=%d hits=%d mismatch_vs_stored=%d mismatch_vs_exact=%d speedup=%.0fx\n",
            R, tf_, ncand, length(hf), nmismatch(mf, lib.stored), nmismatch(mf, mx), te / tf_)
end

# ---- 2. other layouts: fast(R=10) vs exact
rng = MersenneTwister(1)
layouts = Pair{String,Vector}[]
for frac in (0.5, 0.25)
    push!(layouts, "random_$(Int(100frac))pct" =>
          units0[sort(randperm(rng, length(units0))[1:round(Int, frac * length(units0))])])
end
for (name, units) in layouts
    b, _ = set_layout!(units)
    te2 = @elapsed hx2 = LayoutScore.map_hits_exact(lib, units, b, cs, topo)
    tf2 = @elapsed (hf2, _) = LayoutScore.map_hits_fast(lib, grid, units, b, cs, topo; R = 10.0)
    @printf("[%s] modules=%d exact_s=%.1f fast_s=%.2f mismatch_fast_vs_exact=%d/%d\n", name, length(units),
            te2, tf2, nmismatch(LayoutScore.event_modules(lib, hf2), LayoutScore.event_modules(lib, hx2)), lib.nev)
end

# ---- 3. full score for the official layout from fast hits
bvh, geom = set_layout!(units0)
(hf, _) = LayoutScore.map_hits_fast(lib, grid, units0, bvh, cs, topo; R = 10.0)
tt = @elapsed tabs = LayoutScore.event_tables(lib, hf, units0)
ts = @elapsed s = LayoutScore.score(lib, tabs, geom)
@printf("[official score] tables_s=%.2f score_s=%.2f triggered=%d fitted=%d median_dpsi_deg=%.2f\n",
        tt, ts, s.ntrig, s.nfit, s.median_dpsi)
println("DONE")
