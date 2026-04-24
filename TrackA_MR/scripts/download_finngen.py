import requests, sys, os
url = 'https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_E4_GRAVES.gz'
print("Checking if FinnGen Graves summary stats exists...")
if os.path.exists('c:/ProjectTEDGWAS/TrackA_MR/data/finngen_R10_E4_GRAVES.gz'):
    print("Already exists.")
    sys.exit(0)

print(f"Downloading {url}...")
with requests.get(url, stream=True) as r:
    r.raise_for_status()
    with open('c:/ProjectTEDGWAS/TrackA_MR/data/finngen_R10_E4_GRAVES.gz', 'wb') as f:
        for chunk in r.iter_content(chunk_size=1024*1024): 
            if chunk:
                f.write(chunk)
print("FinnGen download complete.")
