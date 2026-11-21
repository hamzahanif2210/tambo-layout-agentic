# Objective: per-event layout results + weights -> S, B, significance, resolution.
#
# Shared by the leaderboard (08_leaderboard.jl) and the search reward, so both
# use one definition.
#
# Conventions
#   θ (truth and reco) = zenith of the direction of TRAVEL (TamboSim): cosmic rays
#   travel downward (θ > 90°); Earth-skimming ν_τ products travel up / near-horizontal.
#   Signal region: fitted event with reco θ̂ <= θcut.
#
# Splits (fixed, deterministic)
#   multi-run sets (ν_τ): by run number   — mod 5 ∈ {0,1,2} train, 3 select, 4 test
#   single-run sets (cosmic rays): by event index — same mod-5 rule
#   Weights are rescaled by 1 / (fraction of the set's MC in the split), so each
#   split is an unbiased rate estimate.
#
# Rates use the collaboration reference fluxes (LayoutScore.flux_nu / flux_cr).
# Background = proton sets weighted with the ALL-PARTICLE flux (all CRs treated as
# protons); Fe-56 is reported separately as a cross-check, never added.
module Objective

using DelimitedFiles, Statistics

const FLUX_NU(E) = 1.8e-18 * 1e4 * (E / 1e5)^(-2.52)   # GeV⁻¹ m⁻² s⁻¹ sr⁻¹ (TamboMakie flux_nu)
const FLUX_CR(E) = 1.8e4 * E^(-2.7)                     # GeV⁻¹ m⁻² s⁻¹ sr⁻¹ (TamboMakie flux_CR)
const YEAR = 3.15576e7
const SPLITS = (train = (0, 1, 2), select = (3,), test = (4,))

function readcsv(p)
    m, h = readdlm(p, ','; header = true)
    return Dict(String(h[j]) => m[:, j] for j in axes(m, 2))
end

kind(set) = startswith(set, "nu_tau") ? :signal : startswith(set, "proton") ? :background :
            startswith(set, "Fe-56") ? :crosscheck : :other

"""
Per-event table for one (layout, set): joins results with weights, assigns
split, and returns weighted-rate weights w [1/yr per event] for every split.
"""
function events(resultcsv, weightcsv, set)
    r = readcsv(resultcsv); w = readcsv(weightcsv)
    key(run, ev) = (lpad(string(Int(run)), 5, '0'), Int(ev))
    wd = Dict(key(w["run"][i], w["event_index"][i]) => (Float64(w["oneweight_GeV_m2_sr"][i]),
              Int(w["n_runs_in_set"][i])) for i in eachindex(w["run"]))
    n = length(r["run"])
    runs = [lpad(string(Int(x)), 5, '0') for x in r["run"]]
    evi = Int.(r["event_index"])
    nruns = length(unique(runs))
    by_run = nruns > 1
    tag = [by_run ? mod(parse(Int, runs[i]), 5) : mod(evi[i], 5) for i in 1:n]
    E = Float64.(r["energy_GeV"])
    flux = kind(set) == :signal ? FLUX_NU : FLUX_CR
    ow = [get(wd, (runs[i], evi[i]), (NaN, 1))[1] for i in 1:n]
    nset = [get(wd, (runs[i], evi[i]), (NaN, 1))[2] for i in 1:n]
    rate = flux.(E) .* ow ./ nset .* YEAR                  # 1/yr if the whole set were used
    frac = Dict(s => (count(t -> t in m, tag) / n) for (s, m) in pairs(SPLITS))
    return (; set, E, tag, rate, frac,
            trig = Int.(r["triggered"]) .== 1, fit = Int.(r["fitted"]) .== 1,
            dpsi = Float64.(r["dpsi_deg"]), rθ = rad2deg.(Float64.(r["reco_theta"])),
            nmissing = count(isnan, ow))
end

"""Weighted rates [1/yr] for one split: triggered, fitted-in-region, plus Δψ of in-region signal."""
function split_rates(ev, split::Symbol, θcut)
    m = [t in SPLITS[split] for t in ev.tag]
    f = ev.frac[split]; f > 0 || return nothing
    inreg = ev.fit .& (ev.rθ .<= θcut)
    w = ev.rate ./ f
    sel = m .& inreg
    return (; trig = sum(w[m .& ev.trig]; init = 0.0), region = sum(w[sel]; init = 0.0),
            n_region = count(sel), w_region = w[sel], dpsi = ev.dpsi[sel], w_trig = w[m .& ev.trig])
end

"Asimov discovery significance for s signal over b background (counts)."
function zasimov(s, b)
    b <= 0 && return s > 0 ? Inf : 0.0
    r = s / b
    r < 1e-4 && return s / sqrt(b)                       # small-signal limit; avoids round-off
    return sqrt(max(0.0, 2 * ((s + b) * log1p(r) - s)))
end

"weighted median"
function wmedian(x, w)
    isempty(x) && return NaN
    o = sortperm(x); c = cumsum(w[o]); x[o][findfirst(>=(c[end] / 2), c)]
end

"""
Layout metrics on one split from a Dict(set => events(...)).
Background floor: if no proton event survives the cut, B = 2.3 × mean triggered
proton rate-weight (≈ 90% CL upper bound of one MC event's worth), flagged.
"""
function metrics(evs::Dict, split::Symbol, θcut; years = 10.0)
    S = 0.0; B = 0.0; Bfe = 0.0; dps = Float64[]; wps = Float64[]; nS = 0; nB = 0; wtrigB = Float64[]
    for (set, ev) in evs
        r = split_rates(ev, split, θcut); isnothing(r) && continue
        k = kind(set)
        if k == :signal
            S += r.region; nS += r.n_region; append!(dps, r.dpsi); append!(wps, r.w_region)
        elseif k == :background
            B += r.region; nB += r.n_region; append!(wtrigB, r.w_trig)
        elseif k == :crosscheck
            Bfe += r.region
        end
    end
    floor_used = nB == 0
    Beff = floor_used ? (isempty(wtrigB) ? 0.0 : 2.3 * mean(wtrigB)) : B
    return (; S, B, Beff, Bfe, nS, nB, floor_used,
            Z = zasimov(S * years, Beff * years), years,
            median_dpsi = wmedian(dps, wps))
end

# =========================================================== objective v2 ===
# Approved 2026-11-18 after the Z≈0 diagnosis (background set by a handful of
# mirror-solution misfits, MC far too small to measure 1e-11 rejection):
#   maximise  S_good = rate of well-reconstructed upgoing ν_τ
#   gate      zero simulated background (proton) events pass the same cuts
#   estimator binned efficiency × generated rate (low variance), raw sum reported too
#
# Quality cut (parameters; u3min / sψmax to be frozen from diag_u3):
#   fitted & θ̂ <= θcut & u3 > u3min (arrived from the air side) & σψ <= sψmax

Base.@kwdef struct QualityCut
    θcut::Float64 = 90.0
    u3min::Float64 = 0.0
    sψmax::Float64 = Inf
end

"""
Uniform per-set table. From an exact results_v3 CSV + weights CSV (`table_exact`)
or a search-batch CSV (`table_batch`). Fields (vectors, one per event):
E, rate (1/yr, full set), tag (split), tθ (truth θ, deg), trig, fit, rθ, u3, sψ, dψ;
plus frac (split => fraction) and set.
"""
function table_exact(resultcsv, weightcsv, set)
    ev = events(resultcsv, weightcsv, set); r = readcsv(resultcsv)
    return (; set, E = ev.E, rate = ev.rate, tag = ev.tag, frac = ev.frac,
            tθ = rad2deg.(Float64.(r["theta"])), trig = ev.trig, fit = ev.fit, rθ = ev.rθ,
            u3 = Float64.(r["u3"]), sψ = Float64.(r["sigma_psi_deg"]), dψ = ev.dpsi)
end

function table_batch(csv, set)
    r = readcsv(csv); n = length(r["energy_GeV"])
    return (; set, E = Float64.(r["energy_GeV"]), rate = Float64.(r["rate_per_yr"]),
            tag = fill(-1, n), frac = Dict(:all => 1.0),           # batch CSVs are already one split
            tθ = Float64.(r["truth_theta_deg"]), trig = Int.(r["triggered"]) .== 1,
            fit = Int.(r["fitted"]) .== 1, rθ = Float64.(r["reco_theta_deg"]),
            u3 = Float64.(r["u3"]), sψ = Float64.(r["sigma_psi_deg"]), dψ = Float64.(r["dpsi_deg"]),
            forced = haskey(r, "forced") ? Int.(r["forced"]) .== 1 : trues(n))
end

"Keep only events where mask m is true (all per-event fields)."
keep(t, m) = merge(t, (; (k => getfield(t, k)[m] for k in (:E, :rate, :tag, :tθ, :trig, :fit, :rθ, :u3, :sψ, :dψ, :forced) if hasproperty(t, k))...))

"Restrict a table to one split, rescaling rates by 1/fraction. :all = no restriction."
function restrict(t, split::Symbol)
    split == :all && return t
    m = [g in SPLITS[split] for g in t.tag]; f = t.frac[split]
    sub(v) = v[m]
    return merge(t, (; E = sub(t.E), rate = sub(t.rate) ./ f, tag = sub(t.tag), tθ = sub(t.tθ),
                     trig = sub(t.trig), fit = sub(t.fit), rθ = sub(t.rθ), u3 = sub(t.u3),
                     sψ = sub(t.sψ), dψ = sub(t.dψ)))
end

passes(t, q::QualityCut) = t.fit .& (t.rθ .<= q.θcut) .& (t.u3 .> q.u3min) .& (t.sψ .<= q.sψmax)

"""
Binned estimator: Σ_bins R_bin · (n_pass / n_events), bins = 0.25 dex in true E ×
truth-θ bands. R_bin = Σ rate over all library events in the bin. Removes the
single-heavy-event noise of the raw weighted sum; assumes detection is
uncorrelated with the remaining weight spread inside a bin (checked by
comparing against the raw sum on exact results).
"""
function binned_rate(t, pass; dlogE = 0.5, θedges = (0.0, 85.0, 95.0, 180.0), nmin = 30)
    # v2 (2026-11-18): 0.25 dex × 7 bands gave sparse bins where one passing event
    # was multiplied by the whole bin rate (binned up to 7× off the raw sum).
    # Now 0.5 dex × 3 zenith bands, and within each band adjacent energy bins are
    # merged (low → high) until a group holds >= nmin events; a short tail joins
    # the last group.
    total = 0.0
    edges = collect(θedges)
    for band in 1:length(edges)-1
        idx = findall(i -> edges[band] <= t.tθ[i] < edges[band+1], eachindex(t.E))
        isempty(idx) && continue
        eb = [floor(Int, log10(t.E[i]) / dlogE) for i in idx]
        groups = Vector{Vector{Int}}(); cur = Int[]
        for b in sort(unique(eb))
            append!(cur, idx[eb .== b])
            length(cur) >= nmin && (push!(groups, cur); cur = Int[])
        end
        if !isempty(cur)
            isempty(groups) ? push!(groups, cur) : append!(groups[end], cur)
        end
        for g in groups
            n = length(g); np = count(pass[g])
            np == 0 && continue
            total += sum(t.rate[g]) * np / n
        end
    end
    return total
end

"""
Objective v2 for one layout. `tables`: Dict(set => table), already restricted
to the wanted split. Returns S (binned) and S_raw [1/yr], signal MC count,
background leak count (proton = gate; Fe-56 reported only), weighted median Δψ.
"""
function metrics_v2(tables::Dict, q::QualityCut)
    S = 0.0; Sraw = 0.0; nS = 0; leak = 0; leakfe = 0; Bleak = 0.0
    dps = Float64[]; wps = Float64[]; per = []
    for (set, t) in tables
        p = passes(t, q); k = kind(set)
        if k == :signal
            sb = binned_rate(t, p); sr = sum(t.rate[p]; init = 0.0)
            S += sb; Sraw += sr; nS += count(p)
            append!(dps, t.dψ[p]); append!(wps, t.rate[p])
            push!(per, (; set, kind = k, S_binned = sb, S_raw = sr, n_pass = count(p), n_trig = count(t.trig)))
        elseif k == :background
            leak += count(p); Bleak += sum(t.rate[p]; init = 0.0)
            push!(per, (; set, kind = k, n_leak = count(p), n_trig = count(t.trig)))
        elseif k == :crosscheck
            leakfe += count(p)
            push!(per, (; set, kind = k, n_leak = count(p), n_trig = count(t.trig)))
        end
    end
    # gate uses the UNWEIGHTED median: the rate-weighted one is set by a few heavy
    # high-energy events (0.7–13° for the same layout across splits); kept for reporting.
    return (; S, Sraw, nS, leak, leakfe, Bleak, median_dpsi = isempty(dps) ? NaN : median(dps),
            median_dpsi_weighted = wmedian(dps, wps), per)
end

"""
Search reward from metrics_v2 (gated, tiered):
  R = 1000·S  − 10·leak  − 2·max(0, median Δψ − 5°)
S in 1/yr (≈0.01–0.3 for current layouts → R ≈ 10–300); each leaking background
MC event costs 10. Hard placement failures are handled by the caller.
"""
# estimator = :raw until the binned estimator is validated against the raw sum
# on exact results (diag_binned.jl); then switch to :binned.
const REWARD_ESTIMATOR = Ref(:binned)   # validated 2026-11-18: binned/raw mean 0.98, split CV 0.48 vs 0.73 (diag_binned)
reward_v2(m; dpsi_gate = 5.0, estimator = REWARD_ESTIMATOR[]) =
    1000 * (estimator == :binned ? m.S : m.Sraw) - 10 * m.leak -
    (isfinite(m.median_dpsi) ? 2 * max(0.0, m.median_dpsi - dpsi_gate) : 0.0)

end # module
