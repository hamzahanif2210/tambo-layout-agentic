# Validate the binned estimator against the unbiased raw weighted sum (exact tier).
# For each fully scored layout and split: S_raw, S_binned at several settings.
# Pass criterion: mean(binned/raw) over layouts×splits ≈ 1 (no bias) and the
# split-to-split spread of binned is smaller than that of raw (less noise).
include(joinpath(@__DIR__, "Objective.jl"))
using .Objective, Statistics, Printf
W = ARGS[1]; θcut = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 90.0
u3min = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : -1.0
q = Objective.QualityCut(θcut = θcut, u3min = u3min)
root = joinpath(W, "results_v3", "exact")
settings = [(0.5, 30), (0.5, 60), (1.0, 30), (0.25, 30)]
ratios = Dict(s => Float64[] for s in settings)
spread_raw = Float64[]; spread_bin = Dict(s => Float64[] for s in settings)
@printf("quality cut: θ̂<=%.1f u3>%.2f\n", θcut, u3min)
for lay in sort(readdir(root))
    tabs = Dict{String,Any}()
    for f in readdir(joinpath(root, lay))
        endswith(f, ".csv") || continue
        set = f[1:end-4]; Objective.kind(set) == :signal || continue
        wp = joinpath(W, "weights", f); isfile(wp) || continue
        tabs[set] = Objective.table_exact(joinpath(root, lay, f), wp, set)
    end
    length(tabs) < 6 && continue
    raws = Float64[]; bins = Dict(s => Float64[] for s in settings)
    line = @sprintf("%-28s", lay)
    for sp in (:train, :select, :test)
        sub = Dict(k => Objective.restrict(t, sp) for (k, t) in tabs)
        r = sum(sum(t.rate[Objective.passes(t, q)]; init = 0.0) for t in values(sub))
        push!(raws, r); line *= @sprintf(" | %s raw=%.4f", sp, r)
        for s in settings
            b = sum(Objective.binned_rate(t, Objective.passes(t, q); dlogE = s[1], nmin = s[2]) for t in values(sub))
            push!(bins[s], b); r > 0 && push!(ratios[s], b / r)
            s == settings[1] && (line *= @sprintf(" bin=%.4f", b))
        end
    end
    println(line)
    cv(x) = std(x) / mean(x)
    push!(spread_raw, cv(raws)); for s in settings; push!(spread_bin[s], cv(bins[s])); end
end
println("\nsplit-to-split coefficient of variation (median over layouts): raw = ", round(median(spread_raw); digits = 3))
for s in settings
    @printf("  dlogE=%.2f nmin=%3d : binned/raw mean=%.3f median=%.3f  |  CV(binned) median=%.3f\n",
            s[1], s[2], mean(ratios[s]), median(ratios[s]), median(spread_bin[s]))
end
