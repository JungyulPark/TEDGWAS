import requests
import pandas as pd

token = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJvcGhqeXBAbmF2ZXIuY29tIiwiaWF0IjoxNzc2MTU2MzExLCJleHAiOjE3NzczNjU5MTF9.j2juecbi9UMa6CHv5zb9C7Rh_r2DLZR9SIAkVD2lO40ydAGowTyQLSL2mT4zLgphD1rlDrkh482uRqn2-PnsNLublr0WiowkcepycDElYNRc_L1sdu0IhfvthQqoibRhQdeuibgJeFxjkalE9c1Mmrf8E12eISuFK0Mq9N9juThjAfeJ-lgakOVl0QcF0f9n2Q-2qEssc5vTgUMq6mg9CVrn5IpW-Yj5G4FbNEaXiPCOIOZV0fXRknoXYAt6lGQBu4tMrD9lv8P7ETavrtW-kQBJWy7yyT8UCXfLB6jWCXu7H1OYKNfxD4BWedJVCAdfHACRtaNFnh6BlWKBcJqLVA"
headers = {"Authorization": f"Bearer {token}"}

print("=== ieu-a-1098 ===")
r = requests.get("https://gwas-api.mrcieu.ac.uk/gwasinfo/ieu-a-1098", headers=headers)
if r.ok:
    d = r.json()
    for k, v in d.items():
        print(f"{v.get('id', '')} | {v.get('trait', '')} | N={v.get('sample_size', 'NA')} (cases={v.get('ncase', 'NA')}) | SNPs={v.get('nsnp', 'NA')} | Year={v.get('year', 'NA')}")
else:
    print(r.status_code, r.text)

print("\n=== Search API ===")
res = requests.get("https://gwas-api.mrcieu.ac.uk/gwasinfo", headers=headers)
if res.ok:
    df = pd.DataFrame(res.json()).T
    if 'trait' in df.columns:
        graves = df[df['trait'].str.contains('Graves', case=False, na=False)]
        print("\nGraves datasets:")
        print(graves[['id', 'trait', 'sample_size', 'ncase', 'ncontrol', 'nsnp', 'year']].to_string())
        
        hyper = df[df['trait'].str.contains('hyperthyroidism', case=False, na=False)]
        print("\nHyperthyroidism datasets:")
        print(hyper[['id', 'trait', 'sample_size', 'ncase', 'ncontrol', 'nsnp', 'year']].to_string())

        thyro = df[df['trait'].str.contains('thyrotoxicosis', case=False, na=False)]
        print("\nThyrotoxicosis datasets:")
        print(thyro[['id', 'trait', 'sample_size', 'ncase', 'ncontrol', 'nsnp', 'year']].to_string())
