# Benchmarks vs optimized layouts (train vs held-out select), colour = S_binned(select)
using CairoMakie
W = ARGS[1]
function load(sp)
    d = Dict{String,Tuple{Float64,Float64}}()
    for l in eachline(joinpath(W, "analysis", "bootstrap_$(sp)_forced_raw.md"))
        c = strip.(split(l, "|"))
        length(c) > 6 && any(startswith(c[2], p) for p in ("plan_", "official", "hex_", "random")) &&
            (d[c[2]] = (parse(Float64, c[3]), parse(Float64, c[6])))
    end; d
end
tr, se = load("train"), load("select")
keep(n) = haskey(se, n) && (occursin("917", n) || n in ("official125","official150","hex_100","hex_125","hex_150","hex_200","random_0.5_1","plan_ls_off125_s50"))
names = sort(filter(keep, collect(keys(tr))))
kind(n) = startswith(n, "plan_llm") ? :llm : startswith(n, "plan_") ? :opt : :bench
mk = Dict(:llm=>:star5, :opt=>:rect, :bench=>:circle)
col(n) = se[n][1]*1e3
cr = (minimum(col.(names)), maximum(col.(names)))
f = Figure(size=(620,540), fontsize=13)
ax = Axis(f[1,1], xlabel="S_raw, train split [10⁻³/yr]", ylabel="S_raw, select split (held out) [10⁻³/yr]",
          title="917-unit layouts: benchmarks vs optimized", titlesize=14)
lines!(ax, [0,40], [0,40], color=:gray80, linewidth=1)
for k in (:bench, :opt, :llm)
    ns = filter(n -> kind(n)==k, names)
    scatter!(ax, [tr[n][2]*1e3 for n in ns], [se[n][2]*1e3 for n in ns], color=col.(ns), colormap=:turbo, colorrange=cr,
             marker=mk[k], markersize=k==:llm ? 22 : 14, strokewidth=0.7, strokecolor=:black, label=string(k))
end
for (n,lab) in [("official125","official125"),("official150","official150"),("hex_100","hex_100"),("plan_llm917_r20_ind_20","LLM (best)"),("plan_de5s50_917_g015_i27","DE (≤50°)")]
    haskey(tr,n) && text!(ax, tr[n][2]*1e3, se[n][2]*1e3, text=lab, offset=(8,-14), fontsize=12)
end
xlims!(ax, 12, 34); ylims!(ax, 12, 34)
Colorbar(f[1,2], colormap=:turbo, limits=cr, label="S_binned, select [10⁻³/yr]")
save(joinpath(W,"talk","figs","scatter_benchmarks.pdf"), f); save(joinpath(W,"talk","figs","scatter_benchmarks.png"), f, px_per_unit=2)
println(length(names), " layouts")
