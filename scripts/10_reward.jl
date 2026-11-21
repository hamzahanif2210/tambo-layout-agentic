# Reward (objective v2) for every plan in a batch, after 09_eval_batch.jl ran on
# all signal + proton sets (Fe-56 optional, reported only).
#
#   julia --project=<env> 10_reward.jl <outdir> <batch_dir> <θcut> <u3min> <sψmax> [channel=all|forced] [estimator=binned|raw]
#
# Objective v3 (2026-11-19): channel=forced estimator=raw — ν_τ upstream-channel events
# (2–186 MC events/set, weights ~1e7× larger) are excluded from the search objective.
#
# Gated, tiered (mirrors the PCB project's compute_rl_reward):
#   invalid plan (over budget / unknown zone) -> -1000 - over_budget - 10·unknown
#   shortfall (zone could not fit its count)   -> scored as placed, −0.1 per missing detector
#   else R = 1000·S_good − 10·proton_leak − 2·max(0, median Δψ − 5°)        (Objective.reward_v2)
#   S_good = binned rate [1/yr] of ν_τ passing: fitted & θ̂<=θcut & u3>u3min & σψ<=sψmax
# Writes <batch_dir>/rewards.toml and <batch_dir>/feedback.md (LLM-readable).
include(joinpath(@__DIR__, "LayoutScore.jl"))
include(joinpath(@__DIR__, "Placer.jl"))
include(joinpath(@__DIR__, "Objective.jl"))
using .LayoutScore, .Placer, .Objective, TOML, Printf, Serialization

W, batch = ARGS[1], ARGS[2]
q = Objective.QualityCut(θcut = parse(Float64, ARGS[3]), u3min = parse(Float64, ARGS[4]),
                         sψmax = parse(Float64, ARGS[5]))
const CHANNEL = length(ARGS) >= 6 ? ARGS[6] : "all"
const ESTIMATOR = Symbol(length(ARGS) >= 7 ? ARGS[7] : "binned")
geo = LayoutScore.load_geometry(joinpath(LayoutScore.OFFICIAL, "datasets/v1/volcano_valley/nu_tau/00010"))
region = Placer.Region(geo.d_frame, geo.cs)

function plan_reward(name)
    jf = joinpath(batch, "$name.jls")
    units, rep = isfile(jf) ? (u = open(deserialize, jf); (u, (; budget = length(u), placed = length(u), over_budget = 0, shortfall = 0, zones = NamedTuple[]))) :
                              Placer.place(region, joinpath(batch, "$name.toml"))
    unknown = count(z -> z.note == "unknown zone", rep.zones)
    # Hard gate only for real plan errors (over budget / unknown zone). A shortfall
    # (zone could not fit its count) is scored on the layout actually placed, with a
    # soft −0.1 per missing detector — fewer detectors already lowers S.
    if rep.over_budget > 0 || unknown > 0
        return (; name, reward = -1000.0 - rep.over_budget - 10unknown, valid = false, rep, m = nothing)
    end
    ed = joinpath(batch, "eval", name)
    tabs = Dict{String,Any}(); parts = Dict{String,Vector{Any}}(); nexp = Dict{String,Int}()
    for f in readdir(ed)
        endswith(f, ".csv") || continue
        mt = match(r"^(.*)\.p(\d+)of(\d+)\.csv$", f)
        set = mt === nothing ? f[1:end-4] : mt[1]
        Objective.kind(set) == :other && continue
        haskey(tabs, set) && error("$ed: $set present both whole and in parts")
        if mt === nothing
            haskey(parts, set) && error("$ed: $set present both whole and in parts")
            tabs[set] = Objective.table_batch(joinpath(ed, f), set)
        else
            push!(get!(parts, set, Any[]), Objective.table_batch(joinpath(ed, f), set)); nexp[set] = parse(Int, mt[3])
        end
    end
    for (set, ps) in parts                                # concatenate event-chunk parts
        length(ps) == nexp[set] || error("$ed: $set has $(length(ps)) of $(nexp[set]) parts")
        cat(k) = reduce(vcat, [getfield(p, k) for p in ps])
        tabs[set] = merge(ps[1], (; E = cat(:E), rate = cat(:rate), tag = cat(:tag), tθ = cat(:tθ), trig = cat(:trig),
                                  fit = cat(:fit), rθ = cat(:rθ), u3 = cat(:u3), sψ = cat(:sψ), dψ = cat(:dψ),
                                  forced = cat(:forced)))
    end
    if CHANNEL == "forced"
        for (s, t) in tabs
            Objective.kind(s) == :signal && (tabs[s] = Objective.keep(t, t.forced))
        end
    end
    m = Objective.metrics_v2(tabs, q)
    return (; name, reward = Objective.reward_v2(m; estimator = ESTIMATOR) - 0.1 * rep.shortfall, valid = true, rep, m)
end

plans = sort(unique([splitext(f)[1] for f in readdir(batch) if (endswith(f, ".toml") && f != "rewards.toml") || endswith(f, ".jls")]))
res = sort([plan_reward(p) for p in plans]; by = r -> -r.reward)

open(joinpath(batch, "rewards.toml"), "w") do io
    TOML.print(io, Dict("quality_cut" => Dict("theta_cut_deg" => q.θcut, "u3_min" => q.u3min, "sigma_psi_max_deg" => q.sψmax),
        "plans" => [Dict("name" => r.name, "reward" => r.reward, "valid" => r.valid,
                         "S_good_per_yr" => r.valid ? r.m.S : NaN, "S_raw_per_yr" => r.valid ? r.m.Sraw : NaN,
                         "proton_leak" => r.valid ? r.m.leak : -1, "median_dpsi_deg" => r.valid ? r.m.median_dpsi : NaN,
                         "placed" => r.rep.placed) for r in res]))
end
open(joinpath(batch, "feedback.md"), "w") do io
    @printf(io, "# Batch feedback (train split, fast tier)\nQuality cut: fitted, θ̂ <= %.1f°, u3 > %.2f, σψ <= %.1f°. ν_τ channel: %s, estimator: %s.\n", q.θcut, q.u3min, q.sψmax, CHANNEL, ESTIMATOR)
    println(io, "Reward = 1000·S_good(raw) − 10·(background leaks) − 2·max(0, median Δψ − 5°). Higher is better.\n")
    for r in res
        @printf(io, "## %s — reward %.2f%s\n", r.name, r.reward, r.valid ? "" : " (INVALID PLAN)")
        if !r.valid
            for z in r.rep.zones
                @printf(io, "- zone %s: requested %d, placed %d %s\n", z.name, z.requested, z.placed, z.note)
            end
            @printf(io, "- budget %d, placed %d, over budget by %d\n\n", r.rep.budget, r.rep.placed, r.rep.over_budget)
            continue
        end
        m = r.m
        r.rep.shortfall > 0 && @printf(io, "- shortfall: %d of %d detectors could not be placed (−0.1 each)\n", r.rep.shortfall, r.rep.budget)
        @printf(io, "- S_good = %.4g /yr (raw weighted sum; binned estimate %.4g, not yet validated), %d signal MC events pass; background leaks: %d proton, %d Fe-56 (x-check); median Δψ = %.2f°\n",
                m.Sraw, m.S, m.nS, m.leak, m.leakfe, m.median_dpsi)
        for p in sort(m.per; by = p -> p.set)
            if p.kind == :signal
                @printf(io, "  - %s: %.4g /yr good (%d MC pass of %d triggered)\n", p.set, p.S_raw, p.n_pass, p.n_trig)
            else
                @printf(io, "  - %s: %d triggered, %d leak\n", p.set, p.n_trig, p.n_leak)
            end
        end
        println(io)
    end
end
for r in res
    r.valid ? @printf("%-28s reward=%8.2f S=%.4g Sraw=%.4g nS=%d leak=%d dpsi=%.2f\n", r.name, r.reward, r.m.S, r.m.Sraw, r.m.nS, r.m.leak, r.m.median_dpsi) :
              @printf("%-28s reward=%8.2f INVALID\n", r.name, r.reward)
end
