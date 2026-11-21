# Why is Z ≈ 0? Per-set S/B contributions, weight concentration, and the
# background events that leak into the signal region (official125, exact).
include(joinpath(@__DIR__, "Objective.jl"))
using .Objective, DelimitedFiles, Statistics, Printf
W = ARGS[1]; lay = length(ARGS) >= 2 ? ARGS[2] : "official125"; θcut = 82.5
root = joinpath(W, "results_v2", "exact", lay)
for f in sort(readdir(root))
    endswith(f, ".csv") || continue
    set = f[1:end-4]; wp = joinpath(W, "weights", f); isfile(wp) || continue
    ev = Objective.events(joinpath(root, f), wp, set)
    r = Objective.readcsv(joinpath(root, f))
    nh = Int.(r["n_hit_modules"]); sψ = Float64.(r["sigma_psi_deg"]); tθ = rad2deg.(Float64.(r["theta"]))
    trig = ev.trig; fit = ev.fit; inreg = fit .& (ev.rθ .<= θcut)
    w = ev.rate                     # full-set rate per event [1/yr]
    tot = sum(w[trig]; init = 0.0)
    top = sort(w[trig]; rev = true)
    @printf("%-14s %-10s trig=%5d (%.3g/yr) top1=%4.0f%% top5=%4.0f%% | fit=%5d inreg=%4d (%.3g/yr) | trig truthθ<90: %d\n",
            set, string(Objective.kind(set)), count(trig), tot,
            isempty(top) ? 0 : 100top[1] / tot, isempty(top) ? 0 : 100sum(top[1:min(5, end)]) / tot,
            count(fit), count(inreg), sum(w[inreg]; init = 0.0), count(trig .& (tθ .< 90)))
    if Objective.kind(set) != :signal
        for i in findall(inreg)
            @printf("    leak: E=%.3g GeV truthθ=%.1f recoθ=%.1f dψ=%.1f σψ=%.2f nmod=%d rate=%.3g/yr\n",
                    ev.E[i], tθ[i], ev.rθ[i], ev.dpsi[i], sψ[i], nh[i], w[i])
        end
    else
        g = inreg
        any(g) && @printf("    signal in region: median dψ=%.1f, frac σψ<2°: %.2f, median nmod=%.0f, truthθ range %.0f–%.0f\n",
                          median(ev.dpsi[g]), mean(sψ[g] .< 2), median(nh[g]), extrema(tθ[g])...)
    end
end
