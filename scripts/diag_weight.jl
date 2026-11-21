include("/n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout/scripts/LayoutScore.jl")
using .LayoutScore, TamboSim, AstroParticleUnits, Statistics
dir = "/n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout/library/full/nu_tau/00010"
for f in filter(f -> startswith(f, "00001_") && endswith(f, ".jls"), readdir(dir))
    lib = LayoutScore.load_library(joinpath(dir, f))
    ow = LayoutScore.oneweights(lib)
    println(f, " ", typeof(lib.psp[1]), " fields=", fieldnames(typeof(lib.psp[1])))
    cds = [ustrip(u"g/cm^2", p.cd) for p in lib.psp if !isnothing(p)]
    println("  cd [g/cm2] median=", median(cds), " max=", maximum(cds))
    for (k, ev) in enumerate(first(lib.evrange):last(lib.evrange))
        ev in (2198, 2895) || continue
        p = lib.psp[k]
        println("  event $ev OW=", ow[k], " E=", lib.truthE[k], " θ=", rad2deg(lib.truthθ[k]))
        for fn in fieldnames(typeof(p)); println("    $fn = ", getfield(p, fn)); end
        σ = p.sigma; m = speedoflight^(-2) * (938.27208816u"MeV" + 939.5654133u"MeV") / 2
        println("    N = σ·cd/m_N = ", uconvert(NoUnits, σ * p.cd / m))
    end
end
