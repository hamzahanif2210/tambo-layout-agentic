# Leaderboard on objective v2 (exact tier, results_v3).
#
#   julia --project=<env> 11_leaderboard_v2.jl <outdir> <θcut> <u3min> <sψmax>
#
# Per layout × split: S_good (binned and raw, 1/yr), signal MC count, proton
# leak count (gate) and Fe-56 leak count (cross-check), weighted median Δψ,
# search reward. Layouts need all signal + proton sets to be listed.
# Writes <outdir>/analysis/v2obj_leaderboard.md
include(joinpath(@__DIR__, "Objective.jl"))
using .Objective, Printf

W = ARGS[1]
q = Objective.QualityCut(θcut = parse(Float64, ARGS[2]), u3min = parse(Float64, ARGS[3]),
                         sψmax = parse(Float64, ARGS[4]))
root = joinpath(W, "results_v3", "exact")
const NEED = 10                              # 6 signal + 4 proton sets
# SPLITS env (default "train,select"): keep the test split blind until the final evaluation
const SPLITS_RUN = Symbol.(split(get(ENV, "SPLITS", "train,select"), ","))

function main(W, root, q)
    rows = String[]
    for lay in sort(readdir(root))
        tabs = Dict{String,Any}()
        for f in readdir(joinpath(root, lay))
            endswith(f, ".csv") || continue
            set = f[1:end-4]; wp = joinpath(W, "weights", f)
            isfile(wp) && (tabs[set] = Objective.table_exact(joinpath(root, lay, f), wp, set))
        end
        n = count(s -> Objective.kind(s) in (:signal, :background), keys(tabs))
        n < NEED && (push!(rows, "| $lay | (incomplete: $n/$NEED sets) | | | | | | |"); continue)
        for sp in SPLITS_RUN
            m = Objective.metrics_v2(Dict(s => Objective.restrict(t, sp) for (s, t) in tabs), q)
            push!(rows, @sprintf("| %s | %s | %.4g | %.4g | %d | %d | %d | %.2f | %.1f |", lay, sp,
                  m.S, m.Sraw, m.nS, m.leak, m.leakfe, m.median_dpsi, Objective.reward_v2(m)))
        end
    end
    adir = joinpath(W, "analysis"); mkpath(adir)
    open(joinpath(adir, get(ENV, "LBNAME", "v2obj_leaderboard.md")), "w") do io
        @printf(io, "Quality cut: fitted & θ̂ <= %.1f° & u3 > %.2f & σψ <= %.1f°. Exact tier.\n", q.θcut, q.u3min, q.sψmax)
        println(io, "S binned = efficiency × generated rate (used in reward, validated: mean ratio to raw 0.98, lower split variance); S_raw = plain weighted sum (unbiased reference). Gate: proton leak = 0.\n")
        println(io, "| layout | split | S binned /yr | S_raw /yr | signal MC | proton leak | Fe-56 leak | median Δψ (°) | reward |")
        println(io, "|---|---|---|---|---|---|---|---|---|")
        foreach(l -> println(io, l), rows)
    end
    return rows
end
foreach(println, main(W, root, q))
