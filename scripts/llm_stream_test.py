# Test: streamed chat completion with the real round-0 prompt; reports time, tokens and parse result per model.
import json, sys, time, urllib.request
sys.argv = ["x", "/n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout", "917", "20", "40", sys.argv[1]]
src = open("/n/netscratch/arguelles_delgado_lab/Everyone/hhanif/tambo_layout/scripts/17_llm_driver.py").read().split("if __name__")[0]
exec(src)
prompt = build_prompt([], None, PER_CALL, 0, 0)
body = {"model": MODEL, "temperature": 0.7, "max_tokens": 32000, "stream": True, **json.loads(os.environ.get("EXTRA", "{}")),
        "messages": [{"role": "system", "content": SYSTEM}, {"role": "user", "content": prompt}]}
t = time.time(); txt, nr = "", 0
req = urllib.request.Request(URL, data=json.dumps(body).encode(), headers={"Authorization": "Bearer " + KEY, "Content-Type": "application/json"})
with urllib.request.urlopen(req, timeout=120) as r:
    for line in r:
        line = line.decode().strip()
        if not line.startswith("data:") or line == "data: [DONE]": continue
        ch = json.loads(line[5:]).get("choices") or []
        if not ch: continue
        d = ch[0].get("delta", {})
        txt += d.get("content") or ""; nr += len(d.get("reasoning_content") or d.get("reasoning") or "")
print(os.environ.get("EXTRA",""), "%s: %.0fs reasoning_chars=%d content_chars=%d parsed=%d" % (MODEL, time.time() - t, nr, len(txt), len(parse_plans(txt))))
