#!/usr/bin/env python3
"""Variant B, step 1: fetch regional extracts from primary public sources.

Written from scratch for the independent reproduction (no code shared with the
Codex pipeline).  Run with `python3 -I` from the idx/ directory so that the
.tbi index files pysam downloads land there and nothing in cwd is importable.

Outputs (all regional, nothing genome-wide is saved):
  data/finngen/<ENDPOINT>__<GENE>.tsv.gz   every FinnGen R12 row in a GRCh38 window
  data/finngen/chrom_scan__<ENDPOINT>__<GENE>.json  whole-chromosome rsID scan (counts only)
  data/gnomad/<GENE>.tsv.gz                gnomAD v2.1.1 genome records at eQTL GRCh37 positions
  logs/fetch_manifest.json
"""
import gzip, json, os, sys, time, hashlib, datetime
import pysam

os.environ.setdefault("CURL_CA_BUNDLE", "<home>")
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EQTL = "/home/user/TEDGWAS/TrackA_MR/data/eqtl_{g}_all.csv"
FG_URL = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_{e}.gz"
GN_URL = "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.{c}.vcf.bgz"
# NOTE: the task text names the original endpoint "E4_GRAVES_OPHT"; that object does not
# exist in the R12 bucket (HTTP 404).  The R12 manifest lists GRAVES_OPHT (858/499,490),
# which is also the exact file (798,476,004 bytes) recorded by the September provenance.
ENDPOINTS = ["GRAVES_OPHT", "E4_GRAVES_STRICT", "E4_GRAVES_OPHT_STRICT"]
GENES = {"TSHR": -466_000, "IGF1R": -543_000, "CTLA4": -865_000}  # approx GRCh37->38 offsets
MARGIN = 1_000_000
GN_KEYS = ["AF", "AC", "AN", "AF_nfe", "AC_nfe", "AN_nfe", "AF_fin", "AC_fin", "AN_fin",
           "allele_type", "n_alt_alleles"]


def log(msg):
    line = datetime.datetime.now(datetime.timezone.utc).isoformat() + " " + msg
    print(line, flush=True)
    with open(os.path.join(BASE, "logs", "fetch.log"), "a") as fh:
        fh.write(line + "\n")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def read_eqtl(g):
    """Minimal CSV reader (no pandas) returning rsID set, positions, chromosome."""
    rs, pos, chrom = set(), set(), None
    with open(EQTL.format(g=g)) as fh:
        hdr = fh.readline().rstrip("\n").split(",")
        i_snp, i_chr, i_pos = hdr.index("SNP"), hdr.index("SNPChr"), hdr.index("SNPPos")
        for line in fh:
            f = line.rstrip("\n").split(",")
            rs.add(f[i_snp]); pos.add(int(f[i_pos]))
            chrom = f[i_chr] if chrom is None else chrom
            assert chrom == f[i_chr]
    return rs, pos, chrom


def fetch_with_retry(tb, chrom, start, end, tries=4):
    for k in range(tries):
        try:
            return list(tb.fetch(chrom, start, end))
        except Exception as exc:  # network hiccup: retry the whole region
            log(f"retry {k+1} {chrom}:{start}-{end}: {exc!r}")
            time.sleep(5 * (k + 1))
    raise RuntimeError(f"fetch failed {chrom}:{start}-{end}")


def main():
    manifest = {"created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                "pysam": pysam.__version__, "finngen": [], "gnomad": []}
    eq = {g: read_eqtl(g) for g in GENES}

    # ---------------- FinnGen R12 regional windows (GRCh38) ----------------
    for e in ENDPOINTS:
        url = FG_URL.format(e=e)
        tb = pysam.TabixFile(url)
        header = tb.header[-1] if tb.header else None
        assert header.startswith("#chrom\tpos\tref\talt\trsids"), header
        for g, off in GENES.items():
            rs, pos37, chrom = eq[g]
            lo = max(1, min(pos37) + off - MARGIN); hi = max(pos37) + off + MARGIN
            out = os.path.join(BASE, "data", "finngen", f"{e}__{g}.tsv.gz")
            t0 = time.time()
            rows = fetch_with_retry(tb, chrom, lo - 1, hi)
            with gzip.open(out, "wt") as fh:
                fh.write(header + "\n")
                for r in rows:
                    fh.write(r + "\n")
            hits = [r.split("\t") for r in rows]
            single = [f for f in hits if f[4] in rs]
            alias = [f for f in hits if "," in f[4] and set(f[4].split(",")) & rs]
            mpos = [int(f[1]) for f in single]
            rec = dict(endpoint=e, gene=g, url=url, chrom=chrom, grch38_start=lo, grch38_end=hi,
                       rows_in_window=len(rows), rows_exact_rsid=len(single),
                       unique_exact_rsid=len({f[4] for f in single}), alias_rows=len(alias),
                       matched_pos_min=min(mpos), matched_pos_max=max(mpos),
                       seconds=round(time.time() - t0, 1), file=out, sha256=sha256(out))
            manifest["finngen"].append(rec); log(json.dumps(rec))

        # whole-chromosome rsID scan (streamed, counts only; nothing saved) to make
        # sure no eQTL rsID sits outside the window (rsID join, not position join)
        for g, off in GENES.items():
            rs, pos37, chrom = eq[g]
            lo = max(1, min(pos37) + off - MARGIN); hi = max(pos37) + off + MARGIN
            t0 = time.time(); n = 0; inside = set(); outside = []
            for k in range(4):
                try:
                    n = 0; inside = set(); outside = []
                    for r in tb.fetch(chrom):
                        n += 1
                        f = r.split("\t", 6)
                        if f[4] in rs:
                            p = int(f[1])
                            if lo <= p <= hi:
                                inside.add(f[4])
                            else:
                                outside.append((f[4], p, f[2], f[3]))
                    break
                except Exception as exc:
                    log(f"chrom scan retry {k+1} {e} {chrom}: {exc!r}"); time.sleep(5 * (k + 1))
            else:
                raise RuntimeError("chromosome scan failed")
            rec = dict(endpoint=e, gene=g, chrom=chrom, rows_on_chrom=n,
                       exact_rsids_inside_window=len(inside), rows_outside_window=len(outside),
                       outside_examples=outside[:20], seconds=round(time.time() - t0, 1))
            with open(os.path.join(BASE, "data", "finngen", f"chrom_scan__{e}__{g}.json"), "w") as fh:
                json.dump(rec, fh, indent=1)
            log(json.dumps({k: v for k, v in rec.items() if k != "outside_examples"}))
        tb.close()

    # ---------------- gnomAD v2.1.1 genomes (GRCh37) ----------------
    for g in GENES:
        rs, pos37, chrom = eq[g]
        url = GN_URL.format(c=chrom)
        tb = pysam.TabixFile(url)
        out = os.path.join(BASE, "data", "gnomad", f"{g}.tsv.gz")
        t0 = time.time(); n_rec = 0; kept = 0
        lo, hi = min(pos37), max(pos37)
        for k in range(4):
            try:
                n_rec = 0; kept_rows = []
                for r in tb.fetch(chrom, lo - 1, hi):
                    n_rec += 1
                    f = r.split("\t", 8)
                    p = int(f[1])
                    if p not in pos37:
                        continue
                    info = {}
                    for item in f[7].split(";"):
                        if "=" in item:
                            key, val = item.split("=", 1)
                            if key in GN_KEYS:
                                info[key] = val
                    kept_rows.append([f[0], f[1], f[2], f[3], f[4], f[6]] + [info.get(x, "") for x in GN_KEYS])
                break
            except Exception as exc:
                log(f"gnomAD retry {k+1} {g}: {exc!r}"); time.sleep(5 * (k + 1))
        else:
            raise RuntimeError("gnomAD fetch failed")
        with gzip.open(out, "wt") as fh:
            fh.write("\t".join(["chrom", "pos", "id", "ref", "alt", "filter"] + GN_KEYS) + "\n")
            for row in kept_rows:
                fh.write("\t".join(row) + "\n")
        rec = dict(gene=g, url=url, chrom=chrom, grch37_start=lo, grch37_end=hi, records_scanned=n_rec,
                   records_at_eqtl_positions=len(kept_rows),
                   eqtl_positions_with_record=len({int(r[1]) for r in kept_rows}),
                   eqtl_positions=len(pos37), seconds=round(time.time() - t0, 1), file=out, sha256=sha256(out))
        manifest["gnomad"].append(rec); log(json.dumps(rec))
        tb.close()

    with open(os.path.join(BASE, "logs", "fetch_manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=1)
    log("FETCH DONE")


if __name__ == "__main__":
    main()
