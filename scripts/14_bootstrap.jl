# Paired event bootstrap of S_good (objective v2) for every exactly-scored layout.
#
#   julia --project=<env> 14_bootstrap.jl <outdir> <θcut> <u3min> <sψmax> [B=1000] [splits=train,select]
#
# All layouts are scored on the same library events, so each replicate resamples
# event indices (with replacement, per signal set) once and applies them to every
# layout: differences between layouts are paired. Reports per layout the S_good
# (binned, as in the reward) and S_raw spread, and ΔS vs the reference layouts
# with 95% intervals and P(Δ > 0). Background leaks are not resampled (counts are
# reported in the leaderboard). The test split is never touched.
# Writes <outdir>/analysis/bootstrap_<split>.md
include(joinpath(@__DIR__, "Objective.jl"))
using .Objective, Random, Statistics, Printf

W = ARGS[1]
q = Objective.QualityCut(θcut = parse(Float64, ARGS[2]), u3min = parse(Float64, ARGS[3]),
                         sψmax = parse(Float64, ARGS[4]))
B = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 1000
SPL = Symbol.(split(length(ARGS) >= 6 ? ARGS[6] : "train,select", ","))
:test in SPL && error("test split stays blind")
root = joinpath(W, "results_v3", "exact")
const REFS = ("official125", "official150")
# CHANNEL env: all (default) | forced | upstream — restrict ν_τ events by injection
# channel (weights/<set>.ptype.csv from 15_point_types.jl)
const CHANNEL = get(ENV, "CHANNEL", "all")
ptypes(set) = (d = Objective.readcsv(joinpath(W, "weights", "$set.ptype.csv"));
               Dict((Int(d["run"][i]), Int(d["event_index"][i])) => String(d["ptype"][i]) for i in eachindex(d["run"])))
const PT = CHANNEL == "all" ? nothing : Dict{String,Any}()
function chmask(resultcsv, set)
    r = Objective.readcsv(resultcsv); pt = get!(() -> ptypes(set), PT, set)
    return [pt[(Int(r["run"][i]), Int(r["event_index"][i]))] == CHANNEL for i in eachindex(r["run"])]
end
filt(t, m) = merge(t, (; E = t.E[m], rate = t.rate[m], tag = t.tag[m], tθ = t.tθ[m], trig = t.trig[m],
                       fit = t.fit[m], rθ = t.rθ[m], u3 = t.u3[m], sψ = t.sψ[m], dψ = t.dψ[m]))

function load(sp)
    L = Dict{String,Dict{String,Any}}()
    for lay in sort(readdir(root))
        tabs = Dict{String,Any}()
        for f in readdir(joinpath(root, lay))
            endswith(f, ".csv") || continue
            set = f[1:end-4]; Objective.kind(set) == :signal || continue
            wp = joinpath(W, "weights", f); isfile(wp) || continue
            t = Objective.table_exact(joinpath(root, lay, f), wp, set)
            CHANNEL == "all" || (t = filt(t, chmask(joinpath(root, lay, f), set)))
            tabs[set] = Objective.restrict(t, sp)
        end
        length(tabs) == 6 && (L[lay] = tabs)
    end
    # pairing requires identical event lists across layouts
    ref = L[REFS[1]]
    for (lay, tabs) in L, (set, t) in tabs
        r = ref[set]
        (t.E == r.E && t.tθ == r.tθ) || error("$lay/$set events not aligned with $(REFS[1])")
    end
    return L
end

sub(t, i) = merge(t, (; E = t.E[i], rate = t.rate[i], tθ = t.tθ[i], fit = t.fit[i], rθ = t.rθ[i],
                       u3 = t.u3[i], sψ = t.sψ[i]))
function S_of(tabs, idx)
    sb = 0.0; sr = 0.0
    for (set, t) in tabs
        s = idx === nothing ? t : sub(t, idx[set]); p = Objective.passes(s, q)
        sb += Objective.binned_rate(s, p); sr += sum(s.rate[p]; init = 0.0)
    end
    return sb, sr
end

q95(x) = (quantile(x, 0.025), quantile(x, 0.975))

for sp in SPL
    L = load(sp); lays = sort(collect(keys(L)))
    sets = sort(collect(keys(L[REFS[1]]))); n = Dict(s => length(L[REFS[1]][s].E) for s in sets)
    rng = MersenneTwister(2026)
    Sb = Dict(l => zeros(B) for l in lays); Sr = Dict(l => zeros(B) for l in lays)
    for b in 1:B
        idx = Dict(s => (n[s] == 0 ? Int[] : rand(rng, 1:n[s], n[s])) for s in sets)
        for l in lays
            Sb[l][b], Sr[l][b] = S_of(L[l], idx)
        end
    end
    pt = Dict(l => S_of(L[l], nothing) for l in lays)
    # EST env: binned (default) | raw — which estimator the paired Δ columns and the ordering use
    est = get(ENV, "EST", "binned"); k = est == "raw" ? 2 : 1; SD = est == "raw" ? Sr : Sb
    out = joinpath(W, "analysis", (CHANNEL == "all" ? "bootstrap_$(sp)" : "bootstrap_$(sp)_$(CHANNEL)") * (est == "raw" ? "_raw.md" : ".md"))
    open(out, "w") do io
        @printf(io, "# Paired event bootstrap, %s split, ν_τ channel = %s, B = %d (exact tier)\n\n", sp, CHANNEL, B)
        println(io, "Paired Δ columns use the $est estimator.\n")
        println(io, "Signal events per set: ", join(["$s=$(n[s])" for s in sets], ", "), "\n")
        println(io, "| layout | S binned | 95% CI | rel. SD | S raw | 95% CI | ΔS vs official125 [95% CI] | P(Δ>0) | ΔS vs official150 [95% CI] | P(Δ>0) |")
        println(io, "|---|---|---|---|---|---|---|---|---|---|")
        for l in sort(lays; by = l -> -pt[l][k])
            lo, hi = q95(Sb[l]); rlo, rhi = q95(Sr[l])
            cells = String[]
            for r in REFS
                d = SD[l] .- SD[r]; dl, dh = q95(d)
                push!(cells, l == r ? "—" : @sprintf("%+.4f [%+.4f, %+.4f]", pt[l][k] - pt[r][k], dl, dh))
                push!(cells, l == r ? "—" : @sprintf("%.2f", mean(d .> 0)))
            end
            @printf(io, "| %s | %.4f | [%.4f, %.4f] | %.0f%% | %.4f | [%.4f, %.4f] | %s |\n", l, pt[l][1], lo, hi,
                    100std(Sb[l]) / mean(Sb[l]), pt[l][2], rlo, rhi, join(cells, " | "))
        end
    end
    println(read(out, String))
end
