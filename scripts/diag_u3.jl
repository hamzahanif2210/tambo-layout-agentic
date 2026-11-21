# Does the "arrived from the air side" requirement (u3 > c) remove the
# mirror-solution background leakage while keeping signal? (results_v3)
include(joinpath(@__DIR__, "Objective.jl"))
using .Objective, Statistics, Printf
W = ARGS[1]; lay = length(ARGS) >= 2 ? ARGS[2] : "official125"
root = joinpath(W, "results_v3", "exact", lay)
sig = Dict{String,Any}(); bkg = Dict{String,Any}()
for f in sort(readdir(root))
    endswith(f, ".csv") || continue
    set = f[1:end-4]; wp = joinpath(W, "weights", f); isfile(wp) || continue
    ev = Objective.events(joinpath(root, f), wp, set)
    r = Objective.readcsv(joinpath(root, f))
    x = (; ev, u3 = Float64.(r["u3"]), sψ = Float64.(r["sigma_psi_deg"]), tθ = rad2deg.(Float64.(r["theta"])))
    Objective.kind(set) == :signal ? (sig[set] = x) : (bkg[set] = x)
end
println("u3 of LEAKING background (fit, θ̂<=90) vs non-leaking:")
for (s, x) in sort(collect(bkg))
    lk = x.ev.fit .& (x.ev.rθ .<= 90); ok = x.ev.fit .& .!lk
    @printf("  %-14s leak n=%3d u3 min/med/max = %s | non-leak n=%5d u3 min/med = %s\n", s, count(lk),
            count(lk) > 0 ? @sprintf("%.3f/%.3f/%.3f", minimum(x.u3[lk]), median(x.u3[lk]), maximum(x.u3[lk])) : "-",
            count(ok), count(ok) > 0 ? @sprintf("%.3f/%.3f", minimum(x.u3[ok]), median(x.u3[ok])) : "-")
end
println("\nscan: cut = fit & θ̂<=θcut & u3>c (all splits pooled, full-set rates)")
for θcut in (85.0, 90.0, 95.0), c in (-1.0, 0.0, 0.05, 0.1, 0.2, 0.3)
    S = 0.0; nS = 0; B = 0.0; nB = 0
    for (_, x) in sig
        m = x.ev.fit .& (x.ev.rθ .<= θcut) .& (x.u3 .> c); S += sum(x.ev.rate[m]; init = 0.0); nS += count(m)
    end
    for (s, x) in bkg
        Objective.kind(s) == :background || continue
        m = x.ev.fit .& (x.ev.rθ .<= θcut) .& (x.u3 .> c); B += sum(x.ev.rate[m]; init = 0.0); nB += count(m)
    end
    @printf("  θcut=%5.1f u3>%5.2f : S=%.4g/yr (n=%4d)  B=%.4g/yr (n=%d)\n", θcut, c, S, nS, B, nB)
end
