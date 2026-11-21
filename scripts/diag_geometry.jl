# What surfaces exist (topography, detector region), and where do the official arrays sit?
include(joinpath(@__DIR__, "LayoutScore.jl"))
using .LayoutScore, TamboSim, AstroParticleUnits, Statistics, Printf
const S = joinpath(LayoutScore.OFFICIAL, "datasets/v1/volcano_valley/nu_tau/00010")
geo = LayoutScore.load_geometry(S)
fr = collect(geo.tf)
for f in fr
    println("frame ", f.stream, " keys=", sort(collect(keys(f.data))))
end
d = geo.d_frame
for k in ("bvh", "detector_bvh")
    b = d[k]
    tris = b.triangles
    cz = [ustrip(u"m", (t.v1.point[3] + t.v2.point[3] + t.v3.point[3]) / 3) for t in tris]
    cx = [ustrip(u"m", (t.v1.point[1] + t.v2.point[1] + t.v3.point[1]) / 3) for t in tris]
    cy = [ustrip(u"m", (t.v1.point[2] + t.v2.point[2] + t.v3.point[2]) / 3) for t in tris]
    up = [0.0, 0.0, 1.0]
    sl = [rad2deg(acos(clamp(abs(TamboSim.normal(t).point[3]), 0, 1))) for t in tris]
    ar = [ustrip(u"m^2", TamboSim.area(t)) for t in tris]
    @printf("%-13s triangles=%d area=%.2f km2  x[%.0f,%.0f] y[%.0f,%.0f] z[%.0f,%.0f]  slope<=35deg area frac=%.2f\n",
            k, length(tris), sum(ar) / 1e6, extrema(cx)..., extrema(cy)..., extrema(cz)...,
            sum(ar[sl .<= 35]) / sum(ar))
end
u = d["detector_units"]
c = [ustrip.(u"m", x.center.point) for x in u]
@printf("official125 centres: n=%d x[%.0f,%.0f] y[%.0f,%.0f] z[%.0f,%.0f]\n", length(c),
        extrema(first.(c))..., extrema(getindex.(c, 2))..., extrema(last.(c))...)
println("OBB fields: ", fieldnames(typeof(u[1])), "  half_lengths=", u[1].half_extents)
println("Rotations binding in TamboSim: ", isdefined(TamboSim, :AngleAxis))
