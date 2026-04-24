import requests

print("=== CTD Graves ophthalmopathy (MESH:D049970) ===")
# Use MeSH disease ID for Thyroid Eye Disease / Graves Ophthalmopathy
for term, meshid in [("Graves' Ophthalmopathy", "D049970"), ("Thyroid Eye Disease", "D049970")]:
    try:
        url = f"https://ctdbase.org/tools/batchQuery.go?report=genes_by_disease&inputType=disease&inputTerms={meshid}&format=tsv"
        r = requests.get(url, timeout=20)
        lines = [l for l in r.text.split('\n') if l and not l.startswith('#')]
        if len(lines) > 1:
            print(f"CTD ({term}): {len(lines)-1} gene associations found")
            for l in lines[1:6]: print(' ', l.split('\t')[0] if '\t' in l else l)
        else:
            print(f"CTD ({term}): 0 associations (empty response)")
    except Exception as e:
        print(f"CTD error: {e}")

print("\n=== DisGeNET Thyroid Eye Disease ===")
try:
    for cui in ["C0342021", "C0011847"]:
        r2 = requests.get(f"https://www.disgenet.org/api/gda/disease/{cui}?source=ALL&format=json&limit=5",
                          headers={"accept": "application/json"}, timeout=20)
        if r2.status_code == 200 and r2.text.strip() not in ["[]", ""]:
            data = r2.json()
            print(f"DisGeNET CUI {cui}: {len(data)} associations found")
            for d in data[:5]: print(' ', d.get('gene_symbol','?'), d.get('score',''))
        else:
            print(f"DisGeNET CUI {cui}: 0 associations. Status={r2.status_code}")
except Exception as e:
    print(f"DisGeNET error: {e}")
