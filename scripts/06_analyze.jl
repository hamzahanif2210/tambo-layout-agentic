# Join per-event layout results with per-event weights -> physical metrics.
#
#   julia --project=<env> 06_analyze.jl <outdir> [mapper=exact]
#
# For every results/<mapper>/<layout>/<set>.csv with a matching weights/<set>.csv:
#   signal (nu_tau_*): ν_τ triggered rate /yr under flux_nu (per-flavour ν+ν̄,
#                      applied to ν_τ-only MC as a proxy for ν+ν̄), fitted rate,
#                      aperture Σ OW / ΔE per decade bin [m² sr], median/68% Δψ
#   background (proton_*, Fe-56_*): triggered rate /yr under the all-particle
#                      flux_cr (each species set weighted as if it were the whole
#                      flux — an upper bound per set, NOT a composition model)
# Writes <outdir>/analysis/<mapper>_leaderboard.md and _bins.csv.

include(joinpath(@__DIR__, "LayoutScore.jl"))
using .LayoutScore, DelimitedFiles, Statistics, Printf

outdir = ARGS[1]; mapper = length(ARGS) >= 2 ? ARGS[2] : "exact"
startswith(abspath(outdir), LayoutScore.OFFICIAL) && error("refusing to write inside official_simulations")

function readcsv(p)
    m, h = readdlm(p, ','; header = true)
    return Dict(String(h[j]) => m[:, j] for j in axes(m, 2))
end

function analyze(outdir, mapper)
    resroot = joinpath(outdir, "results", mapper)
    rows = String[]; bins = String[]
    push!(bins, "layout,set,log10E_lo,log10E_hi,n_events,n_trig,aperture_m2sr")
    for layout in sort(readdir(resroot)), f in sort(readdir(joinpath(resroot, layout)))
        endswith(f, ".csv") || continue
        set = f[1:end-4]
        wp = joinpath(outdir, "weights", f); isfile(wp) || continue
        r = readcsv(joinpath(resroot, layout, f)); w = readcsv(wp)
        key(run, ev) = (string(run), Int(ev))
        ow = Dict(key(w["run"][i], w["event_index"][i]) => (Float64(w["oneweight_GeV_m2_sr"][i]),
                  Int(w["n_runs_in_set"][i])) for i in eachindex(w["run"]))
        isbg = !startswith(set, "nu_tau")
        flux = isbg ? LayoutScore.flux_cr : LayoutScore.flux_nu
        rate = 0.0; ratefit = 0.0; nmiss = 0
        E = Float64.(r["energy_GeV"]); trig = Int.(r["triggered"]) .== 1; fit = Int.(r["fitted"]) .== 1
        owv = zeros(length(E))
        for i in eachindex(E)
            k = key(r["run"][i], r["event_index"][i])
            haskey(ow, k) || (nmiss += 1; continue)
            o, nr = ow[k]; owv[i] = o / nr
            trig[i] && (rate += flux(E[i]) * owv[i])
            fit[i] && (ratefit += flux(E[i]) * owv[i])
        end
        rate *= LayoutScore.SECONDS_PER_YEAR; ratefit *= LayoutScore.SECONDS_PER_YEAR
        dψ = Float64.(r["dpsi_deg"])[fit]
        med = isempty(dψ) ? NaN : median(dψ); q68 = isempty(dψ) ? NaN : quantile(dψ, 0.68)
        push!(rows, @sprintf("| %s | %s | %d | %d | %.3g | %.3g | %.2f | %.2f | %d |", layout, set, length(E),
                             count(trig), rate, ratefit, med, q68, nmiss))
        # aperture per half-decade of true energy
        lE = log10.(E)
        for lo in floor(minimum(lE) * 2) / 2:0.5:maximum(lE)
            sel = (lE .>= lo) .& (lE .< lo + 0.5)
            any(sel) || continue
            ΔE = 10^(lo + 0.5) - 10^lo
            push!(bins, @sprintf("%s,%s,%.1f,%.1f,%d,%d,%.6g", layout, set, lo, lo + 0.5,
                                 count(sel), count(sel .& trig), sum(owv[sel .& trig]) / ΔE))
        end
    end
    adir = joinpath(outdir, "analysis"); mkpath(adir)
    open(joinpath(adir, "$(mapper)_leaderboard.md"), "w") do io
        println(io, "| layout | set | events | triggered | trig rate /yr | fitted rate /yr | median Δψ (°) | 68% Δψ (°) | unweighted |")
        println(io, "|---|---|---|---|---|---|---|---|---|")
        foreach(l -> println(io, l), rows)
        println(io, "\nSignal: flux_nu (per-flavour ν+ν̄, E^-2.52) applied to ν_τ MC. ",
                    "Background: all-particle flux_cr (E^-2.7) applied to each species set — per-set upper bound, not a composition model.")
    end
    open(io -> foreach(l -> println(io, l), bins), joinpath(adir, "$(mapper)_bins.csv"), "w")
    return rows
end
foreach(println, analyze(outdir, mapper))
