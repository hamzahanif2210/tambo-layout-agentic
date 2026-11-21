# Connectivity / model check for the NRP managed LLM endpoint (key read from ~/.nrp_api, never printed).
import json, sys, time, urllib.request
KEY = open("/n/home04/hhanif/.nrp_api").read().strip()
URL = "https://ellm.nrp-nautilus.io/v1/chat/completions"
for model in sys.argv[1:]:
    body = {"model": model, "temperature": 0.7, "max_tokens": 4000,
            "messages": [{"role": "user", "content": 'Reply with only this JSON, filled in: {"plans": [{"zone": "a", "n": 3}]} using two plans.'}]}
    t = time.time()
    try:
        req = urllib.request.Request(URL, data=json.dumps(body).encode(), headers={"Authorization": f"Bearer {KEY}", "Content-Type": "application/json"})
        r = json.load(urllib.request.urlopen(req, timeout=300))
        txt = r["choices"][0]["message"]["content"] or ""
        print(f"{model}: {time.time()-t:.1f}s tokens={r.get('usage',{}).get('total_tokens')} -> {txt[:150]!r}")
    except Exception as e:
        print(f"{model}: FAIL {time.time()-t:.1f}s {str(e)[:200]}")
