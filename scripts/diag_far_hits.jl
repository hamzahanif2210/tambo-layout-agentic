# Which particles hit a module far (> 3 m horizontally) from their recorded
# surface crossing? Characterise by direction vs. the struck module's normal.
include(joinpath(@__DIR__, "LayoutScore.jl"))
using .LayoutScore, TamboSim, LlhamaFit, AstroParticleUnits, Statistics, Printf

simset, libpath = ARGS[1], ARGS[2]
lib = LayoutScore.load_library(libpath)
geo = LayoutScore.load_geometry(simset)
d, cs, topo = geo.d_frame, geo.cs, geo.topo_bvh
units = d["detector_units"]; bvh = d["detector_unit_bvh"]
hx = LayoutScore.map_hits_exact(lib, units, bvh, cs, topo)

# module normal = OBB local z axis (as LlhamaFit's DetectorGeometry does)
nrm(u) = u.axes * [0.0, 0.0, 1.0]
rows = map(hx) do h
    k = h.row; u = units[h.module_index]; c = ustrip.(u"m", u.center.point)
    dist = hypot(lib.x[k] - c[1], lib.y[k] - c[2])
    dvec = [lib.dx[k], lib.dy[k], lib.dz[k]]
    (; dist, dz = lib.dz[k], dn = abs(sum(dvec .* nrm(u))), pdg = lib.pdg[k], E = lib.E[k])
end
far  = filter(r -> r.dist > 3, rows); near = filter(r -> r.dist <= 3, rows)
@printf("hits=%d near(<=3m)=%d far(>3m)=%d\n", length(rows), length(near), length(far))
for (name, s) in (("near", near), ("far", far))
    isempty(s) && continue
    @printf("%-4s |d.n| min=%.3f median=%.3f max=%.3f   |dz| min=%.3f median=%.3f max=%.3f\n", name,
            minimum(r -> r.dn, s), median([r.dn for r in s]), maximum(r -> r.dn, s),
            minimum(r -> abs(r.dz), s), median([abs(r.dz) for r in s]), maximum(r -> abs(r.dz), s))
end
for r in sort(far; by = r -> -r.dist)[1:min(end, 15)]
    @printf("  far: dist=%6.1f m  |d.n|=%.3f  dz=%+.3f  pdg=%d  E=%.3g GeV\n", r.dist, r.dn, r.dz, r.pdg, r.E)
end
# fraction of ALL library particles with |dz| below candidate thresholds
for thr in (0.1, 0.2, 0.3, 0.5)
    @printf("library fraction with |dz| < %.1f : %.3f\n", thr, mean(abs.(lib.dz) .< thr))
end
