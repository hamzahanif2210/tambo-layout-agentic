# Inventory: events per run, available files, and what weighting info a Q/M frame carries.
using TamboSim, TOML
const D = "/n/holylfs05/LABS/arguelles_delgado_lab/Lab/TAMBO/official_simulations/datasets/v1/volcano_valley/nu_tau"

for s in ("00019", "00010", "00011", "00012", "00013", "00014")
    sd = joinpath(D, s)
    runs = sort([f[1:5] for f in readdir(joinpath(sd, "injection")) if endswith(f, ".jld2")])
    nh = isdir(joinpath(sd, "hits_125m_large")) ? length(readdir(joinpath(sd, "hits_125m_large"))) : 0
    np = count(f -> startswith(f, "planned_showers_"), readdir(joinpath(sd, "showers")))
    println("set $s: runs=$(length(runs)) hits_files=$nh planned_showers=$np")
end

# one planned_showers file: event count + keys
ps = joinpath(D, "00011", "showers", "planned_showers_00000.jld2")
tf = load_frames(ps)
println("planned_showers 00011/00000: Q=", length(tf.q_frames))
q = tf.q_frames[1]
println("Q keys: ", sort(collect(keys(q.data))))
m = tf.m_frames[1]
println("M keys: ", sort(collect(keys(m.data))))
inj = m["injection"]
println("M[injection] type: ", typeof(inj))
inj isa AbstractDict && println("M[injection] keys: ", sort(collect(keys(inj))))
println("injection_initial_state: ", q["injection_initial_state"])
haskey(q.data, "phase_space_point") && println("phase_space_point: ", q["phase_space_point"])
