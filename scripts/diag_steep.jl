include(joinpath(@__DIR__, "LayoutScore.jl")); include(joinpath(@__DIR__, "Placer.jl"))
using .LayoutScore, .Placer, TamboSim, AstroParticleUnits, Statistics, Printf
geo = LayoutScore.load_geometry(joinpath(LayoutScore.OFFICIAL, "datasets/v1/volcano_valley/nu_tau/00010"))
r = Placer.Region(geo.d_frame, geo.cs)
u = geo.d_frame["detector_units"]; s = Float64[]
for x in u; c = ustrip.(u"m", x.center.point); h = Placer.project(r, c[1], c[2]); push!(s, isnothing(h) ? NaN : r.slope[h[1]]); end
@printf("official125 slope: median %.1f  p90 %.1f  max %.1f ; >35: %d  >45: %d  >55: %d\n", median(s), quantile(s, .9), maximum(s), count(>(35), s), count(>(45), s), count(>(55), s))
A = r.area; @printf("region area by slope: <=35 %.2f km2, 35-45 %.2f, 45-55 %.2f, >55 %.2f\n", sum(A[r.slope .<= 35]) / 1e6, sum(A[35 .< r.slope .<= 45]) / 1e6, sum(A[45 .< r.slope .<= 55]) / 1e6, sum(A[r.slope .> 55]) / 1e6)
