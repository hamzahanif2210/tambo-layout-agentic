import re,sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
def load(sp):
    d={}
    for l in open(f"analysis/bootstrap_{sp}_forced_raw.md"):
        c=[x.strip() for x in l.split("|")]
        if len(c)>6 and (c[1].startswith(("plan_","official","hex_","random"))):
            d[c[1]]=(float(c[2]),float(c[5]))
    return d
tr,se=load("train"),load("select")
names=[n for n in tr if n in se and ("917" in n or "125" in n or "150" in n or n.startswith("hex_1") or n.startswith("hex_2") or n.startswith("random_0.5"))]
def kind(n):
    if n.startswith("plan_llm"): return "LLM"
    if n.startswith(("plan_de","plan_ls")): return "opt"
    return "bench"
fig,ax=plt.subplots(figsize=(6.2,5.4))
vm,vM=min(se[n][0] for n in names),max(se[n][0] for n in names)
for n in names:
    k=kind(n); x=tr[n][1]*1e3; y=se[n][1]*1e3
    sc=ax.scatter(x,y,c=[se[n][0]*1e3],cmap="turbo",vmin=vm*1e3,vmax=vM*1e3,s=70 if k!="bench" else 55,
                  marker={"LLM":"*","opt":"s","bench":"o"}[k],edgecolors="k",linewidths=.5)
for n,lab in [("official125","official125"),("official150","official150"),("hex_100","hex_100"),("plan_llm917_r20_ind_20","LLM (best)"),("plan_de5s50_917_g015_i27","DE (≤50°)")]:
    if n in tr: ax.annotate(lab,(tr[n][1]*1e3,se[n][1]*1e3),xytext=(6,-12),textcoords="offset points",fontsize=9)
ax.set_xlabel("S_raw, train split  [10$^{-3}$ /yr]"); ax.set_ylabel("S_raw, select split (held out)  [10$^{-3}$ /yr]")
ax.plot([0,40],[0,40],color="#bbb",lw=.8,zorder=0); ax.set_xlim(left=10); ax.set_ylim(bottom=10)
ax.set_title("917-unit layouts: benchmarks vs optimized (★ LLM, ■ DE/LS, ● official/hex)",fontsize=9)
cb=fig.colorbar(sc); cb.set_label("S_binned, select  [10$^{-3}$ /yr]")
fig.tight_layout(); fig.savefig("talk/figs/scatter_benchmarks.pdf"); fig.savefig("talk/figs/scatter_benchmarks.png",dpi=160)
print(len(names),"layouts")
