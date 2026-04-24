"""
TED-TRAP Upgrade — Literature cross-check
Verify INSR novelty claim against Kim 2024 JCI Insight bioRxiv preprint
=============================================================================
Purpose:
    Supplement to R script 13. Uses the publicly available bioRxiv preprint
    of Kim 2024 (which is free and has open full text) to programmatically
    search for INSR mentions in the text before the user spends time
    downloading the journal supplementary tables.

Reference:
    Kim DW et al. 2024 JCI Insight 9(24):e182352
    bioRxiv preprint: doi:10.1101/2024.04.19.590238

Output:
    TrackA_MR/results/13b_kim2024_insr_scan.md
"""

import urllib.request
import re
import os
from datetime import datetime

PROJECT_ROOT = os.environ.get("TEDTRAP_ROOT", r"c:\ProjectTEDGWAS")
OUT_DIR = os.path.join(PROJECT_ROOT, "TrackA_MR", "results")
os.makedirs(OUT_DIR, exist_ok=True)

# bioRxiv full-text URLs
# Note: These URLs may need to be updated; always check the DOI landing page
PAPER_URLS = {
    "Kim2024_biorxiv_abstract": "https://www.biorxiv.org/content/10.1101/2024.04.19.590238v1.abstract",
    "Kim2024_biorxiv_full":     "https://www.biorxiv.org/content/10.1101/2024.04.19.590238v1.full",
    "Kim2021_IOVS_abstract":    "https://iovs.arvojournals.org/article.aspx?articleid=2777099",
}


def fetch_page(url: str) -> str:
    """Fetch URL and return text."""
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            return resp.read().decode("utf-8", errors="ignore")
    except Exception as e:
        print(f"  [ERROR] Fetch failed: {e}")
        return ""


def scan_for_gene(text: str, gene: str) -> dict:
    """Find gene mentions and extract surrounding context."""
    # Case-sensitive word-boundary match (INSR is uppercase gene symbol)
    pattern = r"\b" + re.escape(gene) + r"\b"
    matches = list(re.finditer(pattern, text))
    return {
        "n_mentions": len(matches),
        "first_contexts": [text[max(0, m.start() - 100):m.end() + 100]
                           for m in matches[:5]]
    }


def main():
    print("=" * 70)
    print("Literature cross-check: INSR in Kim 2021/2024")
    print(f"Date: {datetime.now()}")
    print("=" * 70 + "\n")

    report_lines = [
        "# INSR novelty — Literature scan (bioRxiv + abstracts)",
        f"Date: {datetime.now()}",
        ""
    ]

    for label, url in PAPER_URLS.items():
        print(f"Fetching {label}...")
        text = fetch_page(url)
        if not text:
            report_lines.append(f"## {label}\n- Could not fetch ({url})\n")
            continue

        # Remove HTML tags for cleaner search (simple regex — not perfect)
        plain = re.sub(r"<[^>]+>", " ", text)
        plain = re.sub(r"\s+", " ", plain)

        insr_result = scan_for_gene(plain, "INSR")
        print(f"  INSR mentions: {insr_result['n_mentions']}")

        report_lines.append(f"## {label}")
        report_lines.append(f"- URL: {url}")
        report_lines.append(f"- INSR mentions: **{insr_result['n_mentions']}**")

        if insr_result["n_mentions"] > 0:
            report_lines.append("- Contexts (first 5):")
            for i, ctx in enumerate(insr_result["first_contexts"]):
                # Strip and truncate
                ctx_clean = " ".join(ctx.split())[:250]
                report_lines.append(f"  {i+1}. `...{ctx_clean}...`")
        else:
            report_lines.append("- INSR not mentioned in scanned text.")

        # Also check for related terms
        for related in ["insulin receptor", "INS-R", "insulin signaling"]:
            if related.lower() in plain.lower():
                report_lines.append(f"- '{related}' mentioned: YES")
            else:
                report_lines.append(f"- '{related}' mentioned: NO")

        report_lines.append("")

    # --- Verdict ---
    report_lines.append("## Verdict")
    report_lines.append("")
    report_lines.append("Based on the literature scan above, the v3 manuscript should")
    report_lines.append("frame INSR as follows:")
    report_lines.append("")
    report_lines.append("- If INSR is mentioned ≥1 time in Kim 2021 OR Kim 2024: **soften to**")
    report_lines.append('  "*We provide direct transcriptomic evidence of INSR upregulation in')
    report_lines.append('  TED orbital adipose tissue, complementing prior reports of insulin-')
    report_lines.append('  signaling pathway enrichment in TED (Kim 2021; Kim 2024).*"')
    report_lines.append("")
    report_lines.append("- If INSR is NOT mentioned in either: acceptable to write")
    report_lines.append('  "*To our knowledge, direct transcriptomic upregulation of INSR in')
    report_lines.append('  TED orbital adipose tissue has not been previously reported.*"')
    report_lines.append("")
    report_lines.append("**Important**: Scan above is of abstracts/full-text where available")
    report_lines.append("and does NOT include supplementary DEG tables. Final confirmation")
    report_lines.append("requires downloading supplementary files and running script 13.")

    out = os.path.join(OUT_DIR, "13b_kim2024_insr_scan.md")
    with open(out, "w", encoding="utf-8") as f:
        f.write("\n".join(report_lines))
    print(f"\n✅ Saved: {out}")


if __name__ == "__main__":
    main()
