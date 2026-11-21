# Build zones + site context, then expand a plan into detector units.
#
#   julia --project=<env> 07_site_and_placer.jl <outdir> [plan.toml ...]
#
# Writes <outdir>/site/site_context.md, <outdir>/site/zones.toml, and for each
# plan <outdir>/layouts/<plan>.jls (Vector{OBB}) + <plan>.report.toml.
include(joinpath(@__DIR__, "LayoutScore.jl"))
include(joinpath(@__DIR__, "Placer.jl"))
using .LayoutScore, .Placer, TamboSim, AstroParticleUnits, TOML, Serialization, Printf

outdir = ARGS[1]; plans = ARGS[2:end]
startswith(abspath(outdir), LayoutScore.OFFICIAL) && error("refusing to write inside official_simulations")
geo = LayoutScore.load_geometry(joinpath(LayoutScore.OFFICIAL, "datasets/v1/volcano_valley/nu_tau/00010"))
t = @elapsed r = Placer.Region(geo.d_frame, geo.cs)
@printf("region built in %.1f s: %d triangles, %d zones\n", t, length(r.cx), length(Placer.zone_names(r)))

sd = joinpath(outdir, get(ENV, "SITE", "site")); mkpath(sd)   # SITE env: e.g. site_s50 for TAMBO_MAX_SLOPE=50
ctx = Placer.site_context(r); write(joinpath(sd, "site_context.md"), ctx); print(ctx)
open(joinpath(sd, "zones.toml"), "w") do io
    TOML.print(io, Dict("zones" => [Dict(string(k) => v for (k, v) in pairs(z)) for z in Placer.zone_table(r)]))
end

# make_obb must reproduce TamboSim's own boxes: re-place official 125 m units at their (x, y)
let off = geo.d_frame["detector_units"], nbad = 0, nmiss = 0, nchk = 0
    for u in off
        c = ustrip.(u"m", u.center.point)
        hit = Placer.project(r, c[1], c[2])
        isnothing(hit) && (nmiss += 1; continue)
        o = Placer.make_obb(r, hit[2], hit[3]); nchk += 1
        dc = maximum(abs.(ustrip.(u"m", o.center.point) .- c)); da = maximum(abs.(o.axes .- u.axes))
        (dc > 1e-6 || da > 1e-9) && (nbad += 1)
    end
    @printf("make_obb check vs official125: checked=%d differ=%d not_on_region=%d\n", nchk, nbad, nmiss)
end

ld = joinpath(outdir, "layouts"); mkpath(ld)
for p in plans
    name = splitext(basename(p))[1]
    t = @elapsed (units, rep) = Placer.place(r, p)
    open(io -> serialize(io, units), joinpath(ld, "$name.jls"), "w")
    open(joinpath(ld, "$name.report.toml"), "w") do io
        TOML.print(io, Dict("budget" => rep.budget, "placed" => rep.placed, "shortfall" => rep.shortfall,
                            "over_budget" => rep.over_budget,
                            "zones" => [Dict(string(k) => v for (k, v) in pairs(z)) for z in rep.zones]))
    end
    @printf("plan %s: placed %d / budget %d (shortfall %d, over_budget %d) in %.1f s\n",
            name, rep.placed, rep.budget, rep.shortfall, rep.over_budget, t)
    for z in rep.zones
        @printf("   %-20s requested=%4d placed=%4d %s\n", z.name, z.requested, z.placed, z.note)
    end
end
