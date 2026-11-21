# Leaderboard on the significance objective.
#
#   julia --project=<env> 08_leaderboard.jl <outdir> [mapper=exact] [reference=official125]
#
# 1. θcut: scanned on the SELECTION split for the reference layout (maximises Z),
#    then frozen for every layout and split.
# 2. For each layout with results on >=1 signal and >=1 background set:
#    S, B (proton, all-particle flux), Fe-56 cross-check, Z (10 yr, Asimov),
#    weighted median Δψ of in-region signal — on train / select / test.
# Writes <outdir>/analysis/v2_<mapper>_leaderboard.md and v2_<mapper>_thetacut.txt.
include(joinpath(@__DIR__, "Objective.jl"))
using .Objective, Printf

outdir = ARGS[1]; mapper = length(ARGS) >= 2 ? ARGS[2] : "exact"
ref = length(ARGS) >= 3 ? ARGS[3] : "official125"
root = joinpath(outdir, "results_v2", mapper)

function load_layout(dir)
    evs = Dict{String,Any}()
    for f in readdir(dir)
        endswith(f, ".csv") || continue
        set = f[1:end-4]; wp = joinpath(outdir, "weights", f)
        isfile(wp) && (evs[set] = Objective.events(joinpath(dir, f), wp, set))
    end
    return evs
end

function main(root, ref)
    refev = load_layout(joinpath(root, ref))
    scan = [(θ, Objective.metrics(refev, :select, θ)) for θ in 60.0:2.5:120.0]
    θcut, best = scan[argmax([isfinite(m.Z) ? m.Z : -1.0 for (_, m) in scan])]
    io = IOBuffer()
    println(io, "θcut scan on SELECT split, reference layout = $ref (θ = reco zenith of travel direction)")
    for (θ, m) in scan
        @printf(io, "  θcut=%5.1f  S=%.4g/yr  B=%.4g/yr (nB=%d%s)  Z10=%.3f\n", θ, m.S, m.Beff, m.nB,
                m.floor_used ? ", floor" : "", m.Z)
    end
    @printf(io, "chosen θcut = %.1f deg\n", θcut)
    adir = joinpath(outdir, "analysis"); mkpath(adir)
    write(joinpath(adir, "v2_$(mapper)_thetacut.txt"), String(take!(io)))

    rows = String[]
    for lay in sort(readdir(root))
        evs = load_layout(joinpath(root, lay))
        any(Objective.kind(s) == :signal for s in keys(evs)) || continue
        nbg = count(s -> Objective.kind(s) == :background, keys(evs))
        for sp in (:train, :select, :test)
            m = Objective.metrics(evs, sp, θcut)
            push!(rows, @sprintf("| %s | %s | %d/%d | %.4g | %.4g%s | %.3g | %.3f | %.2f |", lay, sp,
                  count(s -> Objective.kind(s) == :signal, keys(evs)), nbg, m.S, m.Beff,
                  m.floor_used ? " (floor)" : "", m.Bfe, m.Z, m.median_dpsi))
        end
    end
    open(joinpath(adir, "v2_$(mapper)_leaderboard.md"), "w") do f
        @printf(f, "θcut = %.1f° (chosen on select split for %s). Z = Asimov significance, 10 yr.\n\n", θcut, ref)
        println(f, "| layout | split | signal/bkg sets | S /yr | B /yr | Fe-56 B /yr (x-check) | Z (10 yr) | median Δψ (°) |")
        println(f, "|---|---|---|---|---|---|---|---|")
        foreach(l -> println(f, l), rows)
    end
    return θcut, rows
end
θcut, rows = main(root, ref)
println(read(joinpath(outdir, "analysis", "v2_$(mapper)_thetacut.txt"), String))
foreach(println, rows)
