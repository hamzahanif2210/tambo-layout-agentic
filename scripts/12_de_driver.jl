# Differential-evolution baseline over plan parameters (search tier).
#
#   julia --project=<env> 12_de_driver.jl <outdir> <budget> [NP=20] [gens=25] [seed=1]
#
# Runs as a long single-core SLURM job; each generation it writes NP plan TOMLs,
# submits 09_eval_batch.jl as an array over the 10 train sets (6 ν_τ + 4 proton),
# waits, runs 10_reward.jl, and does DE/rand/1/bin selection. Resumable: state is
# saved after every generation in <outdir>/search/de_b<budget>_s<seed>/state.jls.
#
# Encoding (all genes in [0,1]); per zone (11 zones from site/zones.toml):
#   share, pattern (tri_grid|contour|cluster|random), spacing (40–250 m), orient (0–60°)
# Decoding keeps plans feasible: shares < 0.1 → zone unused; counts by largest
# remainder; each zone capped at its capacity at 40 m packing (excess handed to
# zones with spare room); grid/contour spacing capped so the count fits the zone.
using TOML, Random, Serialization, Printf

const W = ARGS[1]; const BUDGET = parse(Int, ARGS[2])
const NP = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 20
const GENS = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 25
const SEED = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 1
const F_DE, CR = 0.6, 0.7
const QCUT = ("90", "0.1", "1e9")                       # θcut, u3min, sψmax (frozen 2026-11-18)
const JENV = "/n/home04/hhanif/tambo-org/LlhamaFit.jl/examples"
# OBJ env: "" (v2: all channels, binned) or "forced_raw" (v3: forced channel, raw sum)
const OBJ = get(ENV, "OBJ", "")
const RARGS = OBJ == "forced_raw" ? ("forced", "raw") : ("all", "binned")
const TASKFILE = get(ENV, "TASKFILE", "tasks_eval_sets_split.txt")
const NPLANPART = parse(Int, get(ENV, "NPLANPART", "1"))   # arrays per generation, each scoring 1/NPLANPART of the plans
# SEEDPLAN env: "east" (default: hand-written east grid as individual 1) or "uniform"
# (individual 1 = area-proportional tri_grid in every zone, the closest encodable plan to
# the official uniform arrays). TAG env: extra suffix for the run directory.
const SEEDPLAN = get(ENV, "SEEDPLAN", "east")
const ROOT = joinpath(W, "search", "de_b$(BUDGET)_s$(SEED)" * (OBJ == "" ? "" : "_$OBJ") * get(ENV, "TAG", "")); mkpath(ROOT)
# SITE env selects the zone set (site = slope ≤35°, site_s50 = ≤50° with TAMBO_MAX_SLOPE=50, inherited by all jobs)
const ZONES = TOML.parsefile(joinpath(W, get(ENV, "SITE", "site"), "zones.toml"))["zones"]
const NZ = length(ZONES); const D = 4NZ
const PATTERNS = ("tri_grid", "contour", "cluster", "random")
const MINSP = 40.0

capacity(z) = floor(Int, 0.4 * z["area_km2"] * 1e6 / (sqrt(3) / 2 * MINSP^2))

function decode(x)
    sh = [x[4(i-1)+1] for i in 1:NZ]; sh = [s < 0.1 ? 0.0 : s for s in sh]
    sum(sh) == 0 && (sh .= 1.0)
    raw = BUDGET .* sh ./ sum(sh); n = floor.(Int, raw)
    for i in sortperm(raw .- n; rev = true)[1:BUDGET-sum(n)]; n[i] += 1; end
    cap = capacity.(ZONES)
    for _ in 1:NZ                                         # hand excess to zones with room
        ex = sum(max.(n .- cap, 0)); ex == 0 && break
        n = min.(n, cap); room = findall(i -> sh[i] > 0 && n[i] < cap[i], 1:NZ)
        isempty(room) && (room = findall(i -> n[i] < cap[i], 1:NZ))
        for k in 1:ex; i = room[mod1(k, length(room))]; n[i] < cap[i] && (n[i] += 1); end
    end
    zs = Dict{String,Any}[]
    for (i, z) in enumerate(ZONES)
        n[i] == 0 && continue
        g = x[4(i-1)+1:4i]
        pat = PATTERNS[min(4, floor(Int, g[2] * 4) + 1)]
        sp = 40 + 210 * g[3]
        fit = sqrt(2 * z["area_km2"] * 1e6 / (sqrt(3) * n[i])) * 0.95    # tri-lattice spacing that fits n
        pat in ("tri_grid", "contour") && (sp = max(MINSP, min(sp, fit)))
        e = Dict{String,Any}("name" => z["name"], "n" => n[i], "pattern" => pat,
                             "spacing_m" => round(sp; digits = 1), "orient_deg" => round(60 * g[4]; digits = 1))
        pat == "cluster" && (e["radius_m"] = round(max(2sp, MINSP * sqrt(n[i])); digits = 1))
        pat == "contour" && (e["row_step_m"] = round(max(20.0, sp / 2); digits = 1))
        push!(zs, e)
    end
    return Dict{String,Any}("budget" => BUDGET, "rules" => Dict("min_spacing_m" => MINSP, "seed" => 1), "zones" => zs)
end

sh(cmd) = strip(read(cmd, String))
function wait_job(id)
    while true
        st = split(sh(`sacct -j $id -X -n -o State`))
        !isempty(st) && all(s -> !(s in ("PENDING", "RUNNING", "REQUEUED", "CONFIGURING", "COMPLETING")), st) && return st
        sleep(60)
    end
end

"Evaluate a set of genomes in one batch dir; returns rewards (−Inf if missing)."
function evaluate(dir, X)
    mkpath(dir)
    # a generation interrupted under the old 10-task scheme has whole-set CSVs; use a fresh dir
    isdir(joinpath(dir, "eval")) && !isfile(joinpath(dir, "rewards.toml")) && (dir = dir * "_r"; mkpath(dir))
    for (k, x) in enumerate(X)
        open(io -> TOML.print(io, decode(x)), joinpath(dir, @sprintf("ind_%02d.toml", k)), "w")
    end
    tasks = joinpath(W, "slurm", TASKFILE)
    ntask = countlines(tasks)
    ids = [sh(`sbatch --parsable --array=1-$ntask --export=ALL,TASKS=$tasks,BATCH=$dir,SPLIT=train,PLANPART=$pp,NPLANPART=$NPLANPART $(joinpath(W, "slurm", "09_eval_batch.sh"))`) for pp in 1:NPLANPART]
    for id in ids
        st = wait_job(id); println("  eval array $id: ", join(unique(st), ","))
    end
    rid = sh(`sbatch --parsable -J reward -p serial_requeue -c 1 --mem=8g -t 0:30:00
              -o $(joinpath(W, "logs", "%x-%j.out"))
              --wrap "julia --project=$JENV $(joinpath(W, "scripts", "10_reward.jl")) $W $dir $(QCUT[1]) $(QCUT[2]) $(QCUT[3]) $(RARGS[1]) $(RARGS[2])"`)
    wait_job(rid)
    rf = joinpath(dir, "rewards.toml")
    isfile(rf) || (@warn "no rewards.toml in $dir"; return fill(-Inf, length(X)), [Dict{String,Any}() for _ in X])
    r = Dict(p["name"] => p for p in TOML.parsefile(rf)["plans"])
    return [get(get(r, @sprintf("ind_%02d", k), Dict()), "reward", -Inf) for k in 1:length(X)],
           [get(r, @sprintf("ind_%02d", k), Dict()) for k in 1:length(X)]
end

function seed_uniform_east()                              # the hand-written east-face grid, encoded
    x = fill(0.0, D)
    for (i, z) in enumerate(ZONES)
        startswith(z["name"], "east_facing") && (x[4(i-1)+1] = 1.0)
        x[4(i-1)+2] = 0.1; x[4(i-1)+3] = (110 - 40) / 210; x[4(i-1)+4] = 0.0
    end
    return x
end

function seed_uniform_zones()                           # every zone, share ∝ area, tri_grid at the spacing that fits
    amax = maximum(z["area_km2"] for z in ZONES); x = fill(0.0, D)
    for (i, z) in enumerate(ZONES)
        x[4(i-1)+1] = max(0.1, z["area_km2"] / amax)
        x[4(i-1)+2] = 0.1; x[4(i-1)+3] = 1.0; x[4(i-1)+4] = 0.0   # spacing gene 1 → decoder caps to the fitting spacing
    end
    return x
end

function main()
    sf = joinpath(ROOT, "state.jls"); logf = joinpath(ROOT, "log.csv")
    rng = MersenneTwister(SEED)
    if isfile(sf)
        st = open(deserialize, sf); pop, fit, info, g0 = st.pop, st.fit, st.info, st.gen + 1
        rng = st.rngstate
        println("resuming at generation $g0")
    else
        pop = [rand(rng, D) for _ in 1:NP]; pop[1] = SEEDPLAN == "uniform" ? seed_uniform_zones() : seed_uniform_east()
        println("generation 0: evaluating $NP initial plans (budget $BUDGET)")
        fit, info = evaluate(joinpath(ROOT, "gen_000"), pop); g0 = 1
        open(io -> serialize(io, (; pop, fit, info, gen = 0, rngstate = copy(rng))), sf, "w")
        open(logf, "w") do io
            println(io, "gen,ind,reward,S_good,proton_leak,median_dpsi,best_so_far")
            for k in 1:NP
                println(io, "0,$k,$(fit[k]),$(get(info[k], "S_good_per_yr", NaN)),$(get(info[k], "proton_leak", -1)),$(get(info[k], "median_dpsi_deg", NaN)),$(maximum(fit))")
            end
        end
    end
    for g in g0:GENS
        trials = map(1:NP) do i
            a, b, c = shuffle(rng, setdiff(1:NP, i))[1:3]
            v = clamp.(pop[a] .+ F_DE .* (pop[b] .- pop[c]), 0.0, 1.0)
            jr = rand(rng, 1:D)
            [rand(rng) < CR || j == jr ? v[j] : pop[i][j] for j in 1:D]
        end
        tf, ti = evaluate(joinpath(ROOT, @sprintf("gen_%03d", g)), trials)
        for i in 1:NP
            if tf[i] >= fit[i]; pop[i] = trials[i]; fit[i] = tf[i]; info[i] = ti[i]; end
        end
        b = argmax(fit)
        @printf("gen %3d: best=%.2f (S=%.4g leak=%s) median=%.2f\n", g, fit[b],
                get(info[b], "S_good_per_yr", NaN), get(info[b], "proton_leak", "?"), sort(fit)[NP ÷ 2])
        open(logf, "a") do io
            for k in 1:NP
                println(io, "$g,$k,$(tf[k]),$(get(ti[k], "S_good_per_yr", NaN)),$(get(ti[k], "proton_leak", -1)),$(get(ti[k], "median_dpsi_deg", NaN)),$(fit[b])")
            end
        end
        open(io -> TOML.print(io, decode(pop[b])), joinpath(ROOT, "best_plan.toml"), "w")
        open(io -> serialize(io, (; pop, fit, info, gen = g, rngstate = copy(rng))), sf, "w")
    end
    println("DONE; best reward $(maximum(fit)) -> $(joinpath(ROOT, "best_plan.toml"))")
end
if get(ENV, "DRYRUN", "0") == "1"                        # write 40 decoded random plans, no jobs
    d = joinpath(ROOT, "dryrun"); mkpath(d); r = MersenneTwister(99)
    for k in 1:40
        x = k == 1 ? (SEEDPLAN == "uniform" ? seed_uniform_zones() : seed_uniform_east()) : rand(r, D)
        open(io -> TOML.print(io, decode(x)), joinpath(d, @sprintf("dry_%02d.toml", k)), "w")
    end
    println("wrote 40 plans to $d")
else
    main()
end
