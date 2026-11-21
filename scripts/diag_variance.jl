# Which (set, zenith band, energy group) of the binned S_good estimator drives
# its variance? Groups are fixed from the full split (same rule as
# Objective.binned_rate); bootstrap replicates resample events per set and
# recompute each group's contribution R_g · n_pass/n.
#
#   julia --project=<env> diag_variance.jl <outdir> [B=300] layout ...
include(joinpath(@__DIR__, "Objective.jl"))
using .Objective, Random, Statistics, Printf

W = ARGS[1]; B = parse(Int, ARGS[2]); LAYS = ARGS[3:end]
q = Objective.QualityCut(θcut = 90.0, u3min = 0.1, sψmax = 1e9)
root = joinpath(W, "results_v3", "exact")
SIG = ["nu_tau_00010", "nu_tau_00011", "nu_tau_00012", "nu_tau_00013", "nu_tau_00014", "nu_tau_00019"]
EDGES = (0.0, 85.0, 95.0, 180.0)

"Group label per event (band, first energy bin of group), same merge rule as binned_rate."
function groups(t; dlogE = 0.5, nmin = 30)
    lab = fill((0, 0, 0), length(t.E))
    for band in 1:3
        idx = findall(i -> EDGES[band] <= t.tθ[i] < EDGES[band+1], eachindex(t.E))
        isempty(idx) && continue
        eb = [floor(Int, log10(t.E[i]) / dlogE) for i in idx]
        gs = Vector{Vector{Int}}(); cur = Int[]; lo = Int[]; curlo = typemin(Int)
        for b in sort(unique(eb))
            isempty(cur) && (curlo = b)
            append!(cur, idx[eb .== b])
            length(cur) >= nmin && (push!(gs, cur); push!(lo, curlo); cur = Int[])
        end
        !isempty(cur) && (isempty(gs) ? (push!(gs, cur); push!(lo, curlo)) : append!(gs[end], cur))
        for (k, g) in enumerate(gs), i in g
            lab[i] = (band, lo[k], k)
        end
    end
    return lab
end
contrib(rate, pass, m) = (n = count(m); n == 0 ? 0.0 : sum(rate[m]) * count(pass .& m) / n)

for sp in (:train, :select)
    println("\n#################### split = $sp")
    for lay in LAYS
        tabs = Dict(s => Objective.restrict(Objective.table_exact(joinpath(root, lay, "$s.csv"), joinpath(W, "weights", "$s.csv"), s), sp) for s in SIG)
        rows = []; rng = MersenneTwister(7)
        reps = Dict{Any,Vector{Float64}}(); tot = zeros(B)
        for s in SIG
            t = tabs[s]; p = Objective.passes(t, q); lab = groups(t); n = length(t.E)
            G = unique(lab)
            for gl in G
                m = lab .== Ref(gl)
                gl[1] == 0 && continue
                w = t.rate[m]; ws = sort(w; rev = true)
                push!(rows, (; s, band = gl[1], logE = gl[2] * 0.5, n = count(m), np = count(p .& m),
                            R = sum(w), top1 = ws[1] / sum(w), c = contrib(t.rate, p, m),
                            craw = sum(t.rate[p .& m]; init = 0.0), key = (s, gl)))
                reps[(s, gl)] = zeros(B)
            end
            for b in 1:B
                i = rand(rng, 1:n, n); r = t.rate[i]; pp = p[i]; lb = lab[i]
                for gl in G
                    gl[1] == 0 && continue
                    v = contrib(r, pp, lb .== Ref(gl)); reps[(s, gl)][b] = v; tot[b] += v
                end
            end
        end
        S = sum(r.c for r in rows); vt = var(tot)
        @printf("\n== %s  S_binned=%.4f  S_raw=%.4f  boot rel.SD=%.0f%%\n", lay, S, sum(r.craw for r in rows), 100sqrt(vt) / mean(tot))
        println("  set           band  logE≥  n_ev  n_pass   R_group    top1/R   contrib  raw_contrib  var share")
        for r in sort(rows; by = r -> -var(reps[r.key]))[1:8]
            @printf("  %-13s %4d %6.1f %5d %6d %10.4f %8.2f %9.4f %11.4f %8.0f%%\n", r.s, r.band, r.logE, r.n, r.np,
                    r.R, r.top1, r.c, r.craw, 100var(reps[r.key]) / vt)
        end
    end
end
