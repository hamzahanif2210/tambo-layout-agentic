# Local search over detector positions, seeded from an existing layout
# ((1+λ) evolution strategy, search tier, same reward/jobs as 12_de_driver.jl).
#
#   julia --project=<env> 16_local_search.jl <outdir> <name> <init> [λ=20] [gens=30] [seed=1]
#
# init: official125 | sub:<n>:<seed> (random n-subset of official125) | file:<units.jls>
# Each generation writes λ mutants (Vector{OBB} .jls) of the current layout, scores them
# with 09_eval_batch.jl (train split) + 10_reward.jl (forced channel, raw), and moves to
# the best mutant if it beats the current layout. Rewards are deterministic for a fixed
# layout and event set, so comparisons are paired. A mutant moves m detectors (m drawn
# from MOVES): each move is a local nudge (20–100 m, p = 0.6) or a jump to a random
# allowed site (area-weighted, slope ≤ TAMBO_MAX_SLOPE), keeping ≥ 40 m horizontal spacing.
# Env: TASKFILE, NPLANPART, TAMBO_MAX_SLOPE (inherited by the eval/reward jobs).
# Resumable from search/ls_<name>/state.jls.
include(joinpath(@__DIR__, "LayoutScore.jl"))
include(joinpath(@__DIR__, "Placer.jl"))
using .LayoutScore, .Placer, TamboSim, AstroParticleUnits, Random, Serialization, TOML, Printf

const W = ARGS[1]; const NAME = ARGS[2]; const INIT = ARGS[3]
const LAMBDA = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 20
const GENS = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 30
const SEED = length(ARGS) >= 6 ? parse(Int, ARGS[6]) : 1
const MOVES = (5, 10, 20, 40)
const MINSP = 40.0
const QCUT = ("90", "0.1", "1e9")
const JENV = "/n/home04/hhanif/tambo-org/LlhamaFit.jl/examples"
const TASKFILE = get(ENV, "TASKFILE", "tasks_eval_sets_split43.txt")
const NPLANPART = parse(Int, get(ENV, "NPLANPART", "2"))
const ROOT = joinpath(W, "search", "ls_$NAME"); mkpath(ROOT)

geo = LayoutScore.load_geometry(joinpath(LayoutScore.OFFICIAL, "datasets/v1/volcano_valley/nu_tau/00010"))
const R = Placer.Region(geo.d_frame, geo.cs)
const OK = findall(!isempty, R.zone)
const CUMA = cumsum(R.area[OK])
const TRIS = R.bvh.triangles
xy(u) = (c = ustrip.(u"m", u.center.point); (c[1], c[2]))

function initial(rng)
    off = copy(geo.d_frame["detector_units"])
    INIT == "official125" && return off
    if startswith(INIT, "sub:")
        _, n, s = split(INIT, ":"); r2 = MersenneTwister(parse(Int, s))
        return off[sort(randperm(r2, length(off))[1:parse(Int, n)])]
    end
    startswith(INIT, "file:") && return open(deserialize, INIT[6:end])
    error("unknown init $INIT")
end

"Random point inside an allowed triangle (area-weighted)."
function random_site(rng)
    t = TRIS[OK[searchsortedfirst(CUMA, rand(rng) * CUMA[end])]]
    a, b = rand(rng), rand(rng); a + b > 1 && ((a, b) = (1 - a, 1 - b))
    p(v, i) = ustrip(u"m", v.point[i])
    return (p(t.v1, 1) + a * (p(t.v2, 1) - p(t.v1, 1)) + b * (p(t.v3, 1) - p(t.v1, 1)),
            p(t.v1, 2) + a * (p(t.v2, 2) - p(t.v1, 2)) + b * (p(t.v3, 2) - p(t.v1, 2)))
end

function mutate(units, rng)
    u = copy(units); P = [xy(x) for x in u]; m = rand(rng, MOVES)
    for k in randperm(rng, length(u))[1:min(m, length(u))]
        for _ in 1:100
            x, y = if rand(rng) < 0.6
                ρ, φ = 20 + 80rand(rng), 2π * rand(rng); (P[k][1] + ρ * cos(φ), P[k][2] + ρ * sin(φ))
            else
                random_site(rng)
            end
            h = Placer.project(R, x, y)
            (isnothing(h) || isempty(R.zone[h[1]])) && continue
            any(j != k && (P[j][1] - x)^2 + (P[j][2] - y)^2 < MINSP^2 for j in eachindex(P)) && continue
            u[k] = Placer.make_obb(R, h[2], h[3]); P[k] = (x, y); break
        end
    end
    return u, m
end

sh(cmd) = strip(read(cmd, String))
function wait_job(id)
    while true
        st = split(sh(`sacct -j $id -X -n -o State`))
        !isempty(st) && all(s -> !(s in ("PENDING", "RUNNING", "REQUEUED", "CONFIGURING", "COMPLETING")), st) && return st
        sleep(20)
    end
end

"Score the .jls layouts already in dir; returns Dict name => rewards.toml entry."
function evaluate(dir)
    tasks = joinpath(W, "slurm", TASKFILE); ntask = countlines(tasks)
    ids = [sh(`sbatch --parsable --array=1-$ntask --export=ALL,TASKS=$tasks,BATCH=$dir,SPLIT=train,PLANPART=$pp,NPLANPART=$NPLANPART $(joinpath(W, "slurm", "09_eval_batch.sh"))`) for pp in 1:NPLANPART]
    foreach(wait_job, ids)
    rid = sh(`sbatch --parsable -J ls_reward -p serial_requeue --requeue -c 1 --mem=4g -t 0:30:00
              -o $(joinpath(W, "logs", "%x-%j.out"))
              --wrap "julia --project=$JENV $(joinpath(W, "scripts", "10_reward.jl")) $W $dir $(QCUT[1]) $(QCUT[2]) $(QCUT[3]) forced raw"`)
    wait_job(rid)
    rf = joinpath(dir, "rewards.toml")
    isfile(rf) || return Dict{String,Any}()
    return Dict(p["name"] => p for p in TOML.parsefile(rf)["plans"])
end
gendir(g) = (d = joinpath(ROOT, @sprintf("gen_%03d", g)); isdir(d) && !isfile(joinpath(d, "rewards.toml")) && (d *= "_r"); mkpath(d); d)

function main()
    sf = joinpath(ROOT, "state.jls"); logf = joinpath(ROOT, "log.csv")
    if isfile(sf)
        st = open(deserialize, sf); cur, best, g0, rng = st.cur, st.best, st.gen + 1, st.rng
        println("resuming at generation $g0, current reward $best")
    else
        rng = MersenneTwister(SEED); cur = initial(rng)
        d = gendir(0); open(io -> serialize(io, cur), joinpath(d, "ind_00.jls"), "w")
        r = evaluate(d); best = get(get(r, "ind_00", Dict()), "reward", -Inf)
        open(logf, "w") do io
            println(io, "gen,ind,reward,S_raw,proton_leak,median_dpsi,moved,accepted_best")
            println(io, "0,0,$best,$(get(get(r, "ind_00", Dict()), "S_raw_per_yr", NaN)),$(get(get(r, "ind_00", Dict()), "proton_leak", -1)),$(get(get(r, "ind_00", Dict()), "median_dpsi_deg", NaN)),0,$best")
        end
        open(io -> serialize(io, (; cur, best, gen = 0, rng = copy(rng))), sf, "w")
        open(io -> serialize(io, cur), joinpath(ROOT, "best_units.jls"), "w")
        println("generation 0: $(length(cur)) units, reward $best"); g0 = 1
    end
    for g in g0:GENS
        d = gendir(g); muts = [mutate(cur, rng) for _ in 1:LAMBDA]
        for (k, (u, _)) in enumerate(muts); open(io -> serialize(io, u), joinpath(d, @sprintf("ind_%02d.jls", k)), "w"); end
        r = evaluate(d)
        rw = [get(get(r, @sprintf("ind_%02d", k), Dict()), "reward", -Inf) for k in 1:LAMBDA]
        b = argmax(rw); acc = rw[b] > best
        acc && (cur = muts[b][1]; best = rw[b]; open(io -> serialize(io, cur), joinpath(ROOT, "best_units.jls"), "w"))
        open(logf, "a") do io
            for k in 1:LAMBDA
                e = get(r, @sprintf("ind_%02d", k), Dict())
                println(io, "$g,$k,$(rw[k]),$(get(e, "S_raw_per_yr", NaN)),$(get(e, "proton_leak", -1)),$(get(e, "median_dpsi_deg", NaN)),$(muts[k][2]),$best")
            end
        end
        @printf("gen %3d: best mutant %.2f (moved %d) -> %s, current %.2f\n", g, rw[b], muts[b][2], acc ? "ACCEPT" : "keep", best)
        open(io -> serialize(io, (; cur, best, gen = g, rng = copy(rng))), sf, "w")
    end
    println("DONE; best reward $best -> $(joinpath(ROOT, "best_units.jls"))")
end
if get(ENV, "DRYRUN", "0") == "1"           # mutate the initial layout a few times, check spacing/slope, no jobs
    rng = MersenneTwister(SEED); u0 = initial(rng); println("initial units: ", length(u0))
    for k in 1:3
        u, m = mutate(u0, rng); P = [xy(x) for x in u]
        moved = count(i -> xy(u[i]) != xy(u0[i]), eachindex(u))
        dmin = minimum(sqrt((P[i][1] - P[j][1])^2 + (P[i][2] - P[j][2])^2) for i in eachindex(P) for j in i+1:length(P))
        sl = [R.slope[Placer.project(R, p...)[1]] for p in P]
        @printf("mutant %d: m=%d moved=%d n=%d min_spacing=%.1f m max_slope=%.1f\n", k, m, moved, length(u), dmin, maximum(sl))
    end
else
    main()
end
