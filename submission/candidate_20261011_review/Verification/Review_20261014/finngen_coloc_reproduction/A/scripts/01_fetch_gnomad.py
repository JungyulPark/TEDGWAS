#!/usr/bin/env python3
"""Regional extract of gnomAD v2.1.1 genome sites (GRCh37) at eQTLGen positions.

Usage: python -I 01_fetch_gnomad.py GENE OUTDIR

Streams the gnomAD genome sites VCF over HTTPS with tabix (pysam) across the
eQTL region in 50-kb windows (with retries), keeps only records whose POS is one
of the eQTLGen SNP positions for GENE, and writes a small TSV with the fields
needed for the reference frequency (AF, AF_nfe, AF_fin and allele numbers).
The full gnomAD file is never saved.
"""
import os
import sys
import time
import json
import hashlib

import pandas as pd
import pysam

EQTL_DIR = '/home/user/TEDGWAS/TrackA_MR/data'
URL = ('https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/'
       'vcf/genomes/gnomad.genomes.r2.1.1.sites.{chrom}.vcf.bgz')
KEEP_INFO = ['AC', 'AN', 'AF', 'AC_nfe', 'AN_nfe', 'AF_nfe', 'AC_fin', 'AN_fin', 'AF_fin']
WINDOW = 50_000


def parse_info(info):
    out = {}
    for item in info.split(';'):
        if '=' not in item:
            continue
        k, v = item.split('=', 1)
        if k in KEEP_INFO:
            out[k] = v
    return out


def main():
    gene, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    os.chdir(outdir)  # pysam caches the .tbi in cwd
    e = pd.read_csv(f'{EQTL_DIR}/eqtl_{gene}_all.csv')
    chrom = str(int(e.SNPChr.iloc[0]))
    positions = set(int(p) for p in e.SNPPos)
    lo, hi = int(e.SNPPos.min()), int(e.SNPPos.max())
    url = URL.format(chrom=chrom)
    rows = []
    seen = set()
    n_streamed = 0
    t0 = time.time()
    tb = None
    start = lo - 1  # 0-based
    while start < hi:
        end = min(start + WINDOW, hi)
        for attempt in range(6):
            try:
                if tb is None:
                    tb = pysam.TabixFile(url, index=url + '.tbi')
                chunk = []
                k = 0
                for line in tb.fetch(chrom, start, end):
                    k += 1
                    f = line.split('\t', 8)
                    pos = int(f[1])
                    if pos not in positions or not (start < pos <= end):
                        continue
                    info = parse_info(f[7])
                    chunk.append(dict(chrom=f[0], pos=pos, id=f[2], ref=f[3], alt=f[4], filter=f[6], **info))
                break
            except Exception as exc:  # network hiccup: reopen and retry the window
                print(f'{gene} window {start}-{end} attempt {attempt} failed: {exc!r}', flush=True)
                tb = None
                time.sleep(5 * (attempt + 1))
        else:
            raise RuntimeError(f'{gene}: window {start}-{end} failed after retries')
        n_streamed += k
        for r in chunk:
            key = (r['pos'], r['ref'], r['alt'])
            if key not in seen:
                seen.add(key)
                rows.append(r)
        print(f'{gene} {chrom}:{start + 1}-{end} streamed={n_streamed} kept={len(rows)} '
              f'elapsed={time.time() - t0:.0f}s', flush=True)
        start = end
    out = pd.DataFrame(rows)
    fn = f'gnomad_r2.1.1_genomes_{gene}_eqtl_positions.tsv'
    out.to_csv(fn, sep='\t', index=False)
    meta = dict(gene=gene, url=url, chrom=chrom, region_GRCh37=[lo, hi], n_eqtl_positions=len(positions),
                n_records_streamed=n_streamed, n_records_kept=len(out),
                n_positions_with_record=int(out.pos.nunique()) if len(out) else 0,
                elapsed_s=time.time() - t0,
                sha256=hashlib.sha256(open(fn, 'rb').read()).hexdigest())
    with open(f'gnomad_{gene}_manifest.json', 'w') as fh:
        json.dump(meta, fh, indent=2)
    print(json.dumps(meta), flush=True)


if __name__ == '__main__':
    main()
