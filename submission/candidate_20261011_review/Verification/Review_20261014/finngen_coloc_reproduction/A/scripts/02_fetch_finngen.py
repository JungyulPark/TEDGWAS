#!/usr/bin/env python3
"""Regional extract of FinnGen R12 summary statistics (GRCh38) via remote tabix.

Usage: python -I 02_fetch_finngen.py OUTDIR

For each gene, the eQTLGen region (GRCh37) is shifted by an approximate
GRCh37->GRCh38 offset and padded by 1.5 Mb on each side; every FinnGen row in
that window is saved (regional extract only). Matching to eQTLs is by rsID in
the analysis step, never by position.

Note: the R12 manifest phenocode for the original GO endpoint (858 cases /
499,490 controls) is GRAVES_OPHT, not E4_GRAVES_OPHT (that URL returns 404).
"""
import os
import sys
import time
import json
import hashlib

import pandas as pd
import pysam

EQTL_DIR = '/home/user/TEDGWAS/TrackA_MR/data'
URL = 'https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_{ep}.gz'
ENDPOINTS = ['GRAVES_OPHT', 'E4_GRAVES_STRICT', 'E4_GRAVES_OPHT_STRICT']
OFFSET = {'TSHR': -466_000, 'IGF1R': -543_000, 'CTLA4': -865_000}
PAD = 1_500_000
WINDOW = 250_000


def main():
    outdir = sys.argv[1]
    os.makedirs(outdir, exist_ok=True)
    os.chdir(outdir)
    manifest = []
    for gene in ['TSHR', 'IGF1R', 'CTLA4']:
        e = pd.read_csv(f'{EQTL_DIR}/eqtl_{gene}_all.csv')
        chrom = str(int(e.SNPChr.iloc[0]))
        lo = max(1, int(e.SNPPos.min()) + OFFSET[gene] - PAD)
        hi = int(e.SNPPos.max()) + OFFSET[gene] + PAD
        for ep in ENDPOINTS:
            url = URL.format(ep=ep)
            t0 = time.time()
            tb = None
            header = None
            lines = []
            start = lo - 1
            while start < hi:
                end = min(start + WINDOW, hi)
                for attempt in range(6):
                    try:
                        if tb is None:
                            tb = pysam.TabixFile(url, index=url + '.tbi')
                            header = tb.header[-1].lstrip('#').split('\t')
                        chunk = [ln for ln in tb.fetch(chrom, start, end)
                                 if start < int(ln.split('\t', 2)[1]) <= end]
                        break
                    except Exception as exc:
                        print(f'{gene}/{ep} {start}-{end} attempt {attempt} failed: {exc!r}', flush=True)
                        tb = None
                        time.sleep(5 * (attempt + 1))
                else:
                    raise RuntimeError(f'{gene}/{ep}: window {start}-{end} failed')
                lines.extend(chunk)
                start = end
            fn = f'finngen_R12_{ep}_{gene}_chr{chrom}_{lo}_{hi}.tsv'
            with open(fn, 'w') as fh:
                fh.write('\t'.join(header) + '\n')
                for ln in lines:
                    fh.write(ln + '\n')
            meta = dict(gene=gene, endpoint=ep, url=url, chrom=chrom, window_GRCh38=[lo, hi],
                        n_rows=len(lines), elapsed_s=round(time.time() - t0, 1), file=fn,
                        sha256=hashlib.sha256(open(fn, 'rb').read()).hexdigest())
            manifest.append(meta)
            print(json.dumps(meta), flush=True)
    with open('finngen_manifest.json', 'w') as fh:
        json.dump(manifest, fh, indent=2)


if __name__ == '__main__':
    main()
