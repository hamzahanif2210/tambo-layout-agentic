# LLM plan-writing loop (PCB-LLM style) over the same plan space, placer, eval jobs and
# reward as the DE baseline (12_de_driver.jl), so results are directly comparable.
#
#   python3 17_llm_driver.py <outdir> <budget> [rounds=20] [K=40] [model=glm-5]
#
# Each round: ask the LLM (NRP managed endpoint, OpenAI-compatible) for K plans in
# batches of PER_CALL, given the site zones, the objective, and the scored history;
# repair plans to feasibility (zone names, budget, capacity, spacing ranges); write
# ind_XX.toml; score them with 09_eval_batch.sh (train split) + 10_reward.jl (forced, raw);
# append to history. Round 0 also contains the uniform-zones seed plan that DE starts from.
# Env: TASKFILE, NPLANPART, SITE, TAMBO_MAX_SLOPE (inherited by the eval jobs), TAG.
# Resumable: finished rounds are read back from their rewards.toml.
# Prompts and raw responses are saved per round (llm_XX.json); the API key is never logged.
import json, math, os, random, re, subprocess, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor

W = sys.argv[1]; BUDGET = int(sys.argv[2])
ROUNDS = int(sys.argv[3]) if len(sys.argv) > 3 else 20
K = int(sys.argv[4]) if len(sys.argv) > 4 else 40
MODEL = sys.argv[5] if len(sys.argv) > 5 else "glm-5"
PER_CALL = 10
MINSP = 40.0
PATTERNS = ("tri_grid", "contour", "cluster", "random")
QCUT = ("90", "0.1", "1e9")
JENV = "/n/home04/hhanif/tambo-org/LlhamaFit.jl/examples"
TASKFILE = os.environ.get("TASKFILE", "tasks_eval_sets_split43.txt")
NPLANPART = int(os.environ.get("NPLANPART", "4"))
SITE = os.environ.get("SITE", "site_s50")
ROOT = os.path.join(W, "search", "llm_b%d_%s%s" % (BUDGET, MODEL, os.environ.get("TAG", "")))
os.makedirs(ROOT, exist_ok=True)
URL = "https://ellm.nrp-nautilus.io/v1/chat/completions"
KEY = open("/n/home04/hhanif/.nrp_api").read().strip()
rng = random.Random(1)

# ---------------------------------------------------------------- site zones (minimal TOML read)
def read_zones(path):
    zones, cur = [], None
    for line in open(path):
        line = line.strip()
        if line == "[[zones]]":
            cur = {}; zones.append(cur)
        elif "=" in line and cur is not None:
            k, v = [s.strip() for s in line.split("=", 1)]
            cur[k] = v.strip('"') if v.startswith('"') else float(v)
    return zones

ZONES = read_zones(os.path.join(W, SITE, "zones.toml"))
ZN = {z["name"]: z for z in ZONES}
def capacity(z): return int(math.floor(0.4 * z["area_km2"] * 1e6 / (math.sqrt(3) / 2 * MINSP ** 2)))

# ---------------------------------------------------------------- plans
def repair(p):
    """Feasible plan dict from LLM output (or None)."""
    zs = []
    for e in p.get("zones", []):
        if not isinstance(e, dict) or e.get("name") not in ZN: continue
        try: n = float(e.get("n", 0))
        except (TypeError, ValueError): continue
        if n > 0: zs.append(dict(e, n=n))
    if not zs: return None
    merged = {}                                                   # one entry per zone
    for e in zs:
        merged.setdefault(e["name"], dict(e, n=0.0))["n"] += e["n"]
    zs = list(merged.values()); tot = sum(e["n"] for e in zs)
    raw = [BUDGET * e["n"] / tot for e in zs]; n = [int(math.floor(r)) for r in raw]
    for i in sorted(range(len(zs)), key=lambda i: raw[i] - n[i], reverse=True)[:BUDGET - sum(n)]: n[i] += 1
    cap = [capacity(ZN[e["name"]]) for e in zs]
    for _ in range(len(zs)):                                      # hand excess to zones with room
        ex = sum(max(a - c, 0) for a, c in zip(n, cap))
        if ex == 0: break
        n = [min(a, c) for a, c in zip(n, cap)]; room = [i for i in range(len(zs)) if n[i] < cap[i]]
        if not room: break
        for k in range(ex):
            i = room[k % len(room)]
            if n[i] < cap[i]: n[i] += 1
    out = []
    for e, ni in zip(zs, n):
        if ni == 0: continue
        z = ZN[e["name"]]
        pat = e.get("pattern") if e.get("pattern") in PATTERNS else "tri_grid"
        def num(key, d, lo, hi):
            try: return min(hi, max(lo, float(e.get(key, d))))
            except (TypeError, ValueError): return d
        sp = num("spacing_m", 100.0, MINSP, 250.0)
        fit = math.sqrt(2 * z["area_km2"] * 1e6 / (math.sqrt(3) * ni)) * 0.95
        if pat in ("tri_grid", "contour"): sp = max(MINSP, min(sp, fit))
        q = {"name": e["name"], "n": ni, "pattern": pat, "spacing_m": round(sp, 1),
             "orient_deg": round(num("orient_deg", 0.0, 0.0, 60.0), 1)}
        if pat == "cluster": q["radius_m"] = round(num("radius_m", max(2 * sp, MINSP * math.sqrt(ni)), MINSP, 2000.0), 1)
        if pat == "contour": q["row_step_m"] = round(num("row_step_m", max(20.0, sp / 2), 20.0, 250.0), 1)
        out.append(q)
    return {"budget": BUDGET, "zones": out}

def uniform_plan():                                               # = DE's seed_uniform_zones, decoded
    amax = max(z["area_km2"] for z in ZONES)
    return repair({"zones": [{"name": z["name"], "n": max(0.1, z["area_km2"] / amax), "pattern": "tri_grid",
                              "spacing_m": 250.0, "orient_deg": 0.0} for z in ZONES]})

def key(p): return json.dumps(sorted((e["name"], e["n"], e["pattern"], e["spacing_m"], e["orient_deg"]) for e in p["zones"]))

def write_toml(p, path):
    with open(path, "w") as f:
        f.write("budget = %d\n" % p["budget"])
        for e in p["zones"]:
            f.write("\n[[zones]]\n")
            for k, v in e.items():
                f.write('%s = "%s"\n' % (k, v) if isinstance(v, str) else "%s = %s\n" % (k, v))
        f.write("\n[rules]\nmin_spacing_m = 40.0\nseed = 1\n")

def read_plan_toml(path):
    zs, cur = [], None
    for line in open(path):
        line = line.strip()
        if line == "[[zones]]": cur = {}; zs.append(cur)
        elif line.startswith("[rules]"): cur = None
        elif "=" in line and cur is not None:
            k, v = [s.strip() for s in line.split("=", 1)]
            cur[k] = v.strip('"') if v.startswith('"') else (int(v) if k == "n" else float(v))
    return {"budget": BUDGET, "zones": zs}

# ---------------------------------------------------------------- SLURM
def sh(cmd): return subprocess.run(cmd, stdout=subprocess.PIPE, universal_newlines=True, check=True).stdout.strip()
def wait_job(jid):
    while True:
        st = sh(["sacct", "-j", jid, "-X", "-n", "-o", "State"]).split()
        if st and all(s not in ("PENDING", "RUNNING", "REQUEUED", "CONFIGURING", "COMPLETING") for s in st): return st
        time.sleep(30)

def read_rewards(path):
    plans, cur = {}, None
    for line in open(path):
        line = line.strip()
        if line == "[[plans]]": cur = {}
        elif "=" in line and cur is not None:
            k, v = [s.strip() for s in line.split("=", 1)]
            if v.startswith('"'): val = v.strip('"')
            elif v in ("true", "false"): val = (v == "true")
            else:
                try: val = float(v)
                except ValueError: val = v
            cur[k] = val
            if k == "name": plans[val] = cur
    return plans

def evaluate(d, plans):
    os.makedirs(d, exist_ok=True)
    for k, p in enumerate(plans, 1): write_toml(p, os.path.join(d, "ind_%02d.toml" % k))
    tasks = os.path.join(W, "slurm", TASKFILE); ntask = sum(1 for _ in open(tasks))
    ids = [sh(["sbatch", "--parsable", "--array=1-%d" % ntask,
               "--export=ALL,TASKS=%s,BATCH=%s,SPLIT=train,PLANPART=%d,NPLANPART=%d" % (tasks, d, pp, NPLANPART),
               os.path.join(W, "slurm", "09_eval_batch.sh")]) for pp in range(1, NPLANPART + 1)]
    for i in ids: print("  eval array %s: %s" % (i, ",".join(sorted(set(wait_job(i))))), flush=True)
    rid = sh(["sbatch", "--parsable", "-J", "llm_reward", "-p", "serial_requeue", "--requeue", "-c", "1", "--mem=8g",
              "-t", "0:30:00", "-o", os.path.join(W, "logs", "%x-%j.out"), "--wrap",
              "julia --project=%s %s %s %s %s %s %s forced raw" % (JENV, os.path.join(W, "scripts", "10_reward.jl"), W, d, *QCUT)])
    wait_job(rid)
    rf = os.path.join(d, "rewards.toml")
    return read_rewards(rf) if os.path.isfile(rf) else {}

def feedback_sections(d):
    f = os.path.join(d, "feedback.md")
    if not os.path.isfile(f): return {}
    secs = re.split(r"\n## ", open(f).read())[1:]
    return {s.split(" ", 1)[0]: "\n".join(s.split("\n")[1:]).strip() for s in secs}

# ---------------------------------------------------------------- LLM
SYSTEM = """You design detector layouts for TAMBO, a planned tau-neutrino observatory in a deep valley in the Peruvian Andes (Colca valley). Earth-skimming tau neutrinos interact in the mountain, the tau lepton exits and decays in the air, and the resulting air shower crosses the valley and is recorded by an array of small particle-detector units on the valley slopes. Showers are reconstructed from detector hit times (plane fit); a layout is good if many tau-neutrino showers trigger, reconstruct with an up-going-compatible direction (zenith <= 90 deg), have a good fit quality, and if downward cosmic-ray proton showers do not pass the same cut.

You write layout PLANS, not coordinates. A deterministic placer turns a plan into detector positions on the real terrain. You must answer with JSON only."""

def site_text():
    rows = ["| zone | area km2 | mean slope deg | elevation m (min-max) | centroid x,y m (+x east, +y north) | capacity at 40 m |", "|---|---|---|---|---|---|"]
    for z in ZONES:
        rows.append("| %s | %.2f | %.1f | %.0f-%.0f | %.0f, %.0f | %d |" % (z["name"], z["area_km2"], z["mean_slope"], z["zmin"], z["zmax"], z["cx"], z["cy"], capacity(z)))
    return "\n".join(rows)

def plan_line(p):
    return json.dumps(p["zones"], separators=(",", ":"))

def build_prompt(hist, last_dir, nplans, call_idx, rnd):
    scored = sorted([h for h in hist if h["reward"] > -1e9], key=lambda h: -h["reward"])
    top = scored[:12]; fb = feedback_sections(last_dir) if last_dir else {}
    s = []
    s.append("## Task\nWrite %d new layout plans with a budget of exactly %d detector units (round %d of %d)." % (nplans, BUDGET, rnd, ROUNDS))
    s.append("## Site zones (terrain with slope <= 50 deg; the valley floor runs roughly north-south; most zones face east across the valley)\n" + site_text())
    s.append("""## Plan format
A plan is a list of zone entries: {"name": <zone>, "n": <units>, "pattern": <pattern>, "spacing_m": 40-250, "orient_deg": 0-60, optional "row_step_m" (contour, 20-250), optional "radius_m" (cluster)}.
Patterns: tri_grid = triangular lattice with the given spacing, rotated by orient_deg; contour = rows along elevation contours, row_step_m vertical between rows, spacing_m along a row; cluster = units packed around the zone centre within radius_m; random = uniform random positions. Units are kept >= 40 m apart. If a pattern cannot fit n units in the zone, the placer fills the remainder randomly in that zone; "placed" < budget means units were lost. The n values of all entries must sum to %d; zones not listed get no units. Each zone holds at most its capacity.""" % BUDGET)
    s.append("""## Objective (higher is better)
reward = 1000 * S_good - 10 * proton_leak - 2 * max(0, median_dpsi - 5 deg)
S_good = expected tau-neutrino events per year that trigger, are fitted, have reconstructed zenith <= 90 deg and fit quality u3 > 0.1 (Monte Carlo weighted). proton_leak = number of simulated proton showers passing the same cut (should be 0). median_dpsi = median angular error of passing signal events. Scores come from a fixed simulated event sample; small differences (< ~5%) are within noise.""")
    if top:
        s.append("## Best plans so far (reward, S_good/yr, proton_leak, median_dpsi, placed)")
        for h in top:
            s.append("- R%d %s: reward %.2f, S %.4f, leak %d, dpsi %.2f, placed %d\n  plan: %s" % (
                h["round"], h["name"], h["reward"], h["S"], h["leak"], h["dpsi"], h["placed"], plan_line(h["plan"])))
        if fb:
            last = sorted([h for h in hist if h["round"] == rnd - 1 and h["reward"] > -1e9], key=lambda h: -h["reward"])[:3]
            s.append("## Per-energy-set breakdown for the 3 best plans of the last round (nu_tau sets 00010..00014 rise in energy 0.1-10 PeV, 00019 is the lowest-energy set; proton = background)")
            for h in last:
                s.append("### %s (reward %.2f)\n%s" % (h["name"], h["reward"], fb.get(h["name"], "")))
        prev = [h for h in hist if h["round"] == rnd - 1]
        if prev:
            rs = sorted(h["reward"] for h in prev)
            s.append("## Last round summary\n%d plans, reward min/median/max = %.2f / %.2f / %.2f; %d had proton leaks; %d lost units (placed < budget)." % (
                len(prev), rs[0], rs[len(rs) // 2], rs[-1], sum(h["leak"] > 0 for h in prev), sum(h["placed"] < BUDGET for h in prev)))
            worst = sorted(prev, key=lambda h: h["reward"])[:3]
            s.append("Weakest last-round plans: " + "; ".join("%.2f %s" % (h["reward"], plan_line(h["plan"])) for h in worst))
    mix = ["Mostly refine the best plans (move units between zones, change spacing/pattern/orientation), plus 2-3 bold new ideas.",
           "Mostly explore: try clearly different allocations and patterns from the best plans, plus a few refinements.",
           "Make careful single-change edits of the top 3 plans so the effect of each change can be learned, plus 1-2 new ideas.",
           "Combine ideas from different high-scoring plans; include at least 2 plans that concentrate units where the per-set breakdown suggests most signal comes from."]
    s.append("## Instructions\n%s\nThink about the physics briefly (which zones see showers from tau leptons exiting the opposite slope, density vs coverage trade-off), then output ONLY a JSON object:\n{\"plans\": [{\"idea\": \"<=15 words\", \"zones\": [...]}, ...]} with exactly %d plans. No duplicates of plans already listed." % (mix[call_idx % len(mix)], nplans))
    return "\n\n".join(s)

REASONING = os.environ.get("REASONING", "low")                   # reasoning_effort; high effort hits the ~15 min server cap
def call_llm(prompt, temperature):
    """Streamed chat completion (keeps the connection alive); returns (content, usage-like dict)."""
    body = {"model": MODEL, "temperature": temperature, "max_tokens": 32000, "stream": True,
            "reasoning_effort": REASONING,
            "messages": [{"role": "system", "content": SYSTEM}, {"role": "user", "content": prompt}]}
    partial = ""
    for attempt in range(4):
        try:
            req = urllib.request.Request(URL, data=json.dumps(body).encode(),
                                         headers={"Authorization": "Bearer " + KEY, "Content-Type": "application/json"})
            txt, nreason, finish = "", 0, None
            with urllib.request.urlopen(req, timeout=180) as r:
                for line in r:
                    line = line.decode().strip()
                    if not line.startswith("data:") or line == "data: [DONE]": continue
                    ch = json.loads(line[5:]).get("choices") or []
                    if not ch: continue
                    finish = ch[0].get("finish_reason") or finish
                    d = ch[0].get("delta", {})
                    txt += d.get("content") or ""; nreason += len(d.get("reasoning_content") or d.get("reasoning") or "")
            if txt.strip() and finish == "stop": return txt, {"content_chars": len(txt), "reasoning_chars": nreason, "attempt": attempt + 1}
            if len(parse_plans(txt)) > len(parse_plans(partial)): partial = txt
            print("  LLM stream ended early (finish=%s, content %d chars, reasoning %d), retry %d" % (finish, len(txt), nreason, attempt + 1), flush=True)
        except Exception as e:
            print("  LLM call failed (%s), retry %d" % (str(e)[:120], attempt + 1), flush=True); time.sleep(30 * (attempt + 1))
    return partial, {"truncated": True}

def parse_plans(txt):
    """Plans from the reply; if the JSON is truncated (stream cut), keep every complete plan object."""
    txt = re.sub(r"<think>.*?</think>", "", txt, flags=re.S)
    m = re.search(r"\{.*\}", txt, flags=re.S)
    if m:
        try: return [p for p in json.loads(m.group(0)).get("plans", []) if isinstance(p, dict)]
        except (ValueError, AttributeError): pass
    i = txt.find("[")                                              # salvage: balanced {...} objects inside the plans array
    out, depth, start, instr, esc = [], 0, None, False, False
    for j in range(i + 1 if i >= 0 else 0, len(txt)):
        ch = txt[j]
        if instr:
            if esc: esc = False
            elif ch == "\\": esc = True
            elif ch == '"': instr = False
            continue
        if ch == '"': instr = True
        elif ch == "{":
            if depth == 0: start = j
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0 and start is not None:
                try:
                    o = json.loads(txt[start:j + 1])
                    if isinstance(o, dict) and "zones" in o: out.append(o)
                except ValueError: pass
                start = None
    return out

def perturb(p):                                                   # fallback filler: small random edit of a good plan
    q = json.loads(json.dumps(p)); zs = q["zones"]
    e = rng.choice(zs); e["spacing_m"] = e["spacing_m"] * rng.uniform(0.8, 1.25); e["orient_deg"] = rng.uniform(0, 60)
    if len(zs) > 1:
        a, b = rng.sample(range(len(zs)), 2); m = max(1, int(0.1 * zs[a]["n"])); zs[a]["n"] -= m; zs[b]["n"] += m
    return repair(q)

def propose(hist, last_dir, rnd, d, nplans):
    seen = {key(h["plan"]) for h in hist}; out, ideas, log = [], [], []
    calls = int(math.ceil(nplans / PER_CALL))
    def one(c):
        prompt = build_prompt(hist, last_dir, PER_CALL, c, rnd)
        t = time.time(); txt, usage = call_llm(prompt, 0.7 + 0.1 * (c % 3))
        return c, prompt, txt, usage, time.time() - t
    def absorb(res):
        c, prompt, txt, usage, sec = res
        got = parse_plans(txt); nvalid = 0
        for p in got:
            q = repair(p)
            if q and key(q) not in seen and len(out) < nplans:
                seen.add(key(q)); out.append(q); ideas.append(str(p.get("idea", ""))[:120]); nvalid += 1
        log.append({"call": c, "seconds": round(sec, 1), "usage": usage, "n_parsed": len(got), "n_valid_new": nvalid,
                    "prompt": prompt, "response": txt})
        print("  call %d: %.0fs, %d parsed, %d valid new" % (c, sec, len(got), nvalid), flush=True)
    with ThreadPoolExecutor(max_workers=calls) as ex:              # calls of one round run in parallel
        for res in ex.map(one, range(calls)): absorb(res)
    for c in range(calls, calls + 2):                              # up to 2 extra calls if plans were invalid/duplicate
        if len(out) >= nplans: break
        absorb(one(c))
    json.dump(log, open(os.path.join(d, "llm_calls.json"), "w"), indent=1)
    src = ["llm"] * len(out)
    good = sorted([h for h in hist if h["reward"] > -1e9], key=lambda h: -h["reward"])[:5] or [{"plan": uniform_plan()}]
    while len(out) < nplans:
        q = perturb(rng.choice(good)["plan"])
        if q and key(q) not in seen: seen.add(key(q)); out.append(q); ideas.append("filler perturbation"); src.append("filler")
    return out, ideas, src

# ---------------------------------------------------------------- main loop
def rounddir(r): return os.path.join(ROOT, "round_%03d" % r)

def load_round(d, ideas_src=None):
    r = read_rewards(os.path.join(d, "rewards.toml")); meta = json.load(open(os.path.join(d, "meta.json")))
    hist = []
    for k, (idea, src) in enumerate(zip(meta["ideas"], meta["src"]), 1):
        nm = "ind_%02d" % k; e = r.get(nm, {})
        hist.append({"round": meta["round"], "name": nm, "plan": read_plan_toml(os.path.join(d, nm + ".toml")),
                     "reward": e.get("reward", -1e18) if isinstance(e.get("reward"), float) else -1e18,
                     "S": e.get("S_raw_per_yr", float("nan")), "leak": int(e.get("proton_leak", -1)),
                     "dpsi": e.get("median_dpsi_deg", float("nan")), "placed": int(e.get("placed", 0)), "idea": idea, "src": src})
    return hist

def main():
    hist = []; logf = os.path.join(ROOT, "log.csv")
    for r in range(0, ROUNDS + 1):
        d = rounddir(r)
        if os.path.isfile(os.path.join(d, "rewards.toml")) and os.path.isfile(os.path.join(d, "meta.json")):
            hist += load_round(d); continue
        if os.path.isfile(os.path.join(d, "meta.json")):           # plans written, scoring interrupted: rescore the same plans
            meta = json.load(open(os.path.join(d, "meta.json")))
            plans = [read_plan_toml(os.path.join(d, "ind_%02d.toml" % k)) for k in range(1, len(meta["ideas"]) + 1)]
            d2 = d + "_r%d" % int(time.time()); os.makedirs(d2)
            json.dump(meta, open(os.path.join(d2, "meta.json"), "w"))
            evaluate(d2, plans); os.rename(d, d + "_interrupted_%d" % int(time.time())); os.rename(d2, d)
        else:
            os.makedirs(d, exist_ok=True)
            print("round %d: asking %s for plans" % (r, MODEL), flush=True)
            if r == 0:
                plans, ideas, src = propose(hist, None, r, d, K - 1)
                plans = [uniform_plan()] + plans; ideas = ["uniform seed (same as DE)"] + ideas; src = ["seed"] + src
            else:
                plans, ideas, src = propose(hist, rounddir(r - 1), r, d, K)
            json.dump({"round": r, "ideas": ideas, "src": src}, open(os.path.join(d, "meta.json"), "w"), indent=1)
            evaluate(d, plans)
        new = load_round(d); hist += new
        best = max(h["reward"] for h in hist)
        nb = max(new, key=lambda h: h["reward"])
        print("round %d: best this round %.2f (%s: %s), best so far %.2f, fillers %d" % (
            r, nb["reward"], nb["name"], nb["idea"], best, sum(h["src"] == "filler" for h in new)), flush=True)
        newfile = not os.path.isfile(logf)
        with open(logf, "a") as f:
            if newfile: f.write("round,ind,reward,S_raw,proton_leak,median_dpsi,placed,src,best_so_far,idea\n")
            for h in new:
                f.write('%d,%s,%s,%s,%d,%s,%d,%s,%s,"%s"\n' % (r, h["name"], h["reward"], h["S"], h["leak"], h["dpsi"], h["placed"], h["src"], best, h["idea"].replace('"', "'")))
        bh = max(hist, key=lambda h: h["reward"]); write_toml(bh["plan"], os.path.join(ROOT, "best_plan.toml"))
    print("DONE; best reward %.2f -> %s" % (max(h["reward"] for h in hist), os.path.join(ROOT, "best_plan.toml")))

if __name__ == "__main__":
    if os.environ.get("DRYRUN") == "1":                            # one LLM call, repair, write plans; no SLURM
        d = os.path.join(ROOT, "dryrun"); os.makedirs(d, exist_ok=True)
        plans, ideas, src = propose([], None, 0, d, PER_CALL)
        for k, (p, i, s) in enumerate(zip(plans, ideas, src), 1):
            write_toml(p, os.path.join(d, "dry_%02d.toml" % k)); print(s, sum(e["n"] for e in p["zones"]), i, plan_line(p)[:200])
        u = uniform_plan(); print("uniform seed:", sum(e["n"] for e in u["zones"]), plan_line(u)[:300])
    else:
        main()
