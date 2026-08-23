#!/usr/bin/env python3
"""
Build a single self-contained HTML review page for the TED-TRAP v5 manuscript.

Converts the locked master markdown to HTML (pandoc), embeds the five figures as
data URIs at their legend positions, and wraps everything in a print-like scholarly
layout. Output is content-only (no <html>/<head>/<body>) for the Artifact wrapper.

Usage: python3 scripts/24_build_review_artifact.py
Out  : /tmp/artifact/manuscript_review.html
"""
import base64, os, re, subprocess, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MASTER = f"{ROOT}/MANUSCRIPT_TED_TRAP_v5_MASTER.md"
FIGDIR = f"{ROOT}/submission/figures"
OUTDIR = "/tmp/artifact"
os.makedirs(OUTDIR, exist_ok=True)

FIGS = {  # legend key -> (file, alt)
    "Figure 1.":  ("Figure1.png",  "Study design and analytical workflow"),
    "Figure 2.":  ("Figure2.png",  "Druggable-gene-wide BBJ discovery MR volcano plot"),
    "Figure 3.":  ("Figure3.png",  "Multi-layer evidence integration for backbone genes"),
    "Figure S1.": ("FigureS1.png", "TSHR fine-mapping and LD rationale"),
    "Figure S2.": ("FigureS2.png", "Outcome-dependent colocalization of candidate hits"),
}

def md_to_html(path):
    r = subprocess.run(["pandoc", path, "-t", "html", "--no-highlight"],
                       capture_output=True, text=True, check=True)
    return r.stdout

def data_uri(p):
    return "data:image/png;base64," + base64.b64encode(open(p, "rb").read()).decode()

def insert_figures(html):
    """Put each figure immediately BEFORE its legend paragraph."""
    for key, (fname, alt) in FIGS.items():
        fp = os.path.join(FIGDIR, fname)
        if not os.path.exists(fp):
            print(f"  ! missing {fp}", file=sys.stderr); continue
        # legend paragraphs look like: <p><strong>Figure 1. ...</strong> ...</p>
        pat = re.compile(r'(<p><strong>' + re.escape(key) + r')')
        img = (f'<figure class="fig"><img src="{data_uri(fp)}" alt="{alt}" loading="lazy">'
               f'</figure>\n')
        html, n = pat.subn(img + r"\1", html, count=1)
        print(f"  {key:11} embedded: {'yes' if n else 'NO'} ({os.path.getsize(fp)//1024} KB)")
    return html

def build_toc(html):
    """Anchor h2s and return a jump nav."""
    items = []
    def slug(t):
        return re.sub(r"[^a-z0-9]+", "-", re.sub("<[^>]+>", "", t).lower()).strip("-")
    def repl(m):
        text = m.group(2); s = slug(text)
        items.append((s, re.sub("<[^>]+>", "", text)))
        return f'<h2 id="{s}"{m.group(1)}>{text}</h2>'
    html = re.sub(r"<h2([^>]*)>(.*?)</h2>", repl, html, flags=re.S)
    nav = "\n".join(f'<a href="#{s}">{t}</a>' for s, t in items)
    return html, nav

CSS = """
<style>
/* ---- tokens: light (complete palette lives here) ---- */
:root{
  --paper:#fcfcfd; --panel:#f4f6f9; --ink:#14181f; --ink-soft:#4d5766;
  --rule:#dbe1e9; --rule-soft:#e9edf2;
  --accent:#1f4e79; --accent-soft:#e7eef6; --accent-ink:#1f4e79;
  --anchor:#a4262c;   /* TSHR red, from the figures */
  --effector:#2e75b6; /* IGF1R blue, from the figures */
  --shadow:0 1px 2px rgba(20,24,31,.05), 0 8px 24px rgba(20,24,31,.05);
}
@media (prefers-color-scheme: dark){
  :root:not([data-theme="light"]){
    --paper:#12151a; --panel:#1a1f27; --ink:#e4e8ee; --ink-soft:#9aa5b4;
    --rule:#2b323c; --rule-soft:#232932;
    --accent:#7fb0dd; --accent-soft:#18232f; --accent-ink:#9cc4e8;
    --anchor:#e5787d; --effector:#7fb0dd;
    --shadow:0 1px 2px rgba(0,0,0,.4), 0 8px 24px rgba(0,0,0,.35);
  }
}
:root[data-theme="dark"]{
  --paper:#12151a; --panel:#1a1f27; --ink:#e4e8ee; --ink-soft:#9aa5b4;
  --rule:#2b323c; --rule-soft:#232932;
  --accent:#7fb0dd; --accent-soft:#18232f; --accent-ink:#9cc4e8;
  --anchor:#e5787d; --effector:#7fb0dd;
  --shadow:0 1px 2px rgba(0,0,0,.4), 0 8px 24px rgba(0,0,0,.35);
}

*{box-sizing:border-box}
body{
  margin:0; background:var(--paper); color:var(--ink);
  font-family:Georgia,"Iowan Old Style","Palatino Linotype",Palatino,"Source Serif 4",serif;
  font-size:17px; line-height:1.68; -webkit-text-size-adjust:100%;
}
.sans{font-family:system-ui,-apple-system,"Segoe UI",Roboto,"Helvetica Neue",sans-serif}

/* ---- masthead ---- */
.masthead{border-bottom:1px solid var(--rule); background:var(--panel)}
.masthead .inner{max-width:74rem;margin:0 auto;padding:2.6rem 1.5rem 2rem;
  display:flex;flex-direction:column;gap:1.1rem}
.eyebrow{font-family:system-ui,sans-serif;font-size:.7rem;letter-spacing:.14em;
  text-transform:uppercase;color:var(--ink-soft);display:flex;gap:.75rem;flex-wrap:wrap;align-items:center}
.eyebrow .dot{width:4px;height:4px;border-radius:50%;background:var(--accent);display:inline-block}
h1.title{margin:0;font-size:clamp(1.5rem,3.6vw,2.3rem);line-height:1.22;font-weight:600;
  letter-spacing:-.012em;text-wrap:balance;max-width:34ch}
.byline{font-family:system-ui,sans-serif;font-size:.9rem;color:var(--ink-soft);line-height:1.6}
.byline b{color:var(--ink);font-weight:600}

/* at-a-glance strip: the locked headline numbers reviewers check first */
.glance{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:1px;
  background:var(--rule);border:1px solid var(--rule);border-radius:3px;overflow:hidden;margin-top:.4rem}
.glance div{background:var(--paper);padding:.8rem .95rem}
.glance dt{font-family:system-ui,sans-serif;font-size:.68rem;letter-spacing:.09em;
  text-transform:uppercase;color:var(--ink-soft);margin:0 0 .3rem}
.glance dd{margin:0;font-size:1.18rem;font-weight:600;font-variant-numeric:tabular-nums;letter-spacing:-.01em}
.glance dd small{display:block;font-size:.74rem;font-weight:400;color:var(--ink-soft);
  font-family:system-ui,sans-serif;letter-spacing:0;margin-top:.15rem;line-height:1.4}
.g-anchor dd{color:var(--anchor)} .g-effector dd{color:var(--effector)}

/* ---- jump nav ---- */
nav.jump{position:sticky;top:0;z-index:10;background:var(--paper);
  border-bottom:1px solid var(--rule);backdrop-filter:saturate(1.6) blur(6px)}
nav.jump .inner{max-width:74rem;margin:0 auto;padding:.55rem 1.5rem;
  display:flex;gap:.15rem;overflow-x:auto;scrollbar-width:thin}
nav.jump a{font-family:system-ui,sans-serif;font-size:.78rem;color:var(--ink-soft);
  text-decoration:none;padding:.34rem .6rem;border-radius:3px;white-space:nowrap}
nav.jump a:hover{background:var(--accent-soft);color:var(--accent-ink)}
nav.jump a:focus-visible{outline:2px solid var(--accent);outline-offset:1px}

/* ---- article ---- */
main{max-width:74rem;margin:0 auto;padding:2.2rem 1.5rem 5rem}
article{max-width:68ch}
article h2{font-size:1.32rem;font-weight:600;letter-spacing:-.008em;margin:2.9rem 0 .2rem;
  padding-top:1.1rem;border-top:1px solid var(--rule);scroll-margin-top:3.4rem;text-wrap:balance}
article h2:first-of-type{margin-top:1rem}
article h3{font-family:system-ui,sans-serif;font-size:.79rem;font-weight:650;letter-spacing:.1em;
  text-transform:uppercase;color:var(--accent-ink);margin:2rem 0 .1rem}
article p{margin:.95rem 0}
article strong{font-weight:650}
article em{font-style:italic}
article a{color:var(--accent-ink)}
article ul{padding-left:1.15rem;margin:.9rem 0}
article li{margin:.45rem 0}
article li>p{margin:.4rem 0}
hr{border:0;border-top:1px solid var(--rule);margin:2.6rem 0}
code{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:.87em;
  background:var(--panel);padding:.1em .35em;border-radius:3px}

/* tables break out of the measure and scroll on their own */
.tablewrap{overflow-x:auto;margin:1.1rem 0 1.4rem;border:1px solid var(--rule);border-radius:3px;
  max-width:min(100%,74rem);box-shadow:var(--shadow)}
table{border-collapse:collapse;width:100%;font-family:system-ui,sans-serif;font-size:.795rem;
  font-variant-numeric:tabular-nums;background:var(--paper)}
thead th{background:var(--panel);text-align:left;font-weight:650;color:var(--ink);
  padding:.6rem .7rem;border-bottom:1px solid var(--rule);white-space:nowrap;vertical-align:bottom}
tbody td{padding:.5rem .7rem;border-bottom:1px solid var(--rule-soft);vertical-align:top}
tbody tr:last-child td{border-bottom:0}
tbody tr:hover td{background:var(--accent-soft)}
colgroup{display:none}
table em{font-style:italic}

/* figures */
figure.fig{margin:2rem 0 .6rem;max-width:min(100%,74rem)}
figure.fig img{width:100%;height:auto;display:block;border:1px solid var(--rule);
  border-radius:3px;background:#fff;box-shadow:var(--shadow)}
/* the legend paragraph that follows a figure */
figure.fig + p{font-family:system-ui,sans-serif;font-size:.83rem;line-height:1.6;
  color:var(--ink-soft);max-width:68ch;margin-top:.7rem}
figure.fig + p strong{color:var(--ink)}

footer.note{max-width:74rem;margin:0 auto;padding:1.6rem 1.5rem 3rem;border-top:1px solid var(--rule);
  font-family:system-ui,sans-serif;font-size:.78rem;color:var(--ink-soft);line-height:1.65}

@media (max-width:640px){
  body{font-size:16px}
  .masthead .inner{padding:1.9rem 1.1rem 1.5rem}
  main{padding:1.6rem 1.1rem 3.5rem}
  nav.jump .inner{padding:.5rem 1.1rem}
}
@media print{
  nav.jump{display:none} .masthead{background:none}
  body{font-size:10.5pt} .tablewrap{box-shadow:none} figure.fig img{box-shadow:none}
  article h2{break-after:avoid} figure.fig{break-inside:avoid}
}
@media (prefers-reduced-motion:reduce){*{animation:none!important;transition:none!important}}
</style>
"""

def main():
    print("building review artifact…")
    html = md_to_html(MASTER)

    # strip pandoc's own title block (h1 + author paragraphs) — the masthead carries it
    html = re.sub(r"<h1[^>]*>.*?</h1>", "", html, count=1, flags=re.S)

    html = insert_figures(html)
    html, nav = build_toc(html)

    # wrap tables for horizontal scroll
    html = re.sub(r"(<table[^>]*>.*?</table>)", r'<div class="tablewrap">\1</div>', html, flags=re.S)

    page = f"""<title>TSHR versus IGF1R</title>
{CSS}
<header class="masthead">
  <div class="inner">
    <p class="eyebrow"><span>Manuscript · v5</span><span class="dot"></span><span>Review copy — the journal file is the .docx, not this page</span><span class="dot"></span><span>Endocrine Connections</span></p>
    <h1 class="title">Genetic susceptibility and therapeutic target biology diverge at <em>TSHR</em> and <em>IGF1R</em> in Graves disease and thyroid eye disease</h1>
    <p class="byline"><b>Jungyul Park</b><sup>1</sup> · <b>Min-Seon Kim</b><sup>2</sup> · <b>Kyung-Hwa Shin</b><sup>3</sup>* · <b>Suk-Woo Yang</b><sup>1</sup>*<br>
    <sup>1</sup>Seoul St. Mary's Hospital, The Catholic University of Korea · <sup>2</sup>College of Medicine, The Catholic University of Korea · <sup>3</sup>Pusan National University Hospital · *corresponding</p>
    <dl class="glance">
      <div class="g-anchor"><dt>TSHR — anchor</dt><dd>0.951 / 0.986<small>coloc PP.H4, BBJ / FinnGen · shared variant rs179252</small></dd></div>
      <div class="g-effector"><dt>IGF1R — effector-compatible profile</dt><dd>0.79 / 0.68<small>coloc PP.H2 — no shared-variant expression colocalization</small></dd></div>
      <div><dt>Druggable genes screened</dt><dd>2,545<small>6,136 instruments · 13 Bonferroni hits</small></dd></div>
      <div><dt>Robust novel targets</dt><dd>0<small>after cross-outcome colocalization filtering</small></dd></div>
      <div><dt>Detectable effect</dt><dd>OR 2.55<small>median, 80% power at the discovery threshold (Table S8)</small></dd></div>
    </dl>
  </div>
</header>

<nav class="jump" aria-label="Section navigation"><div class="inner">{nav}</div></nav>

<main><article>
{html}
</article></main>

<footer class="note">
  Rendered from the locked master <code>MANUSCRIPT_TED_TRAP_v5_MASTER.md</code>. Tables 1–3 and
  Supplementary Tables S1–S8 appear in full; Figures 1–3 and S1–S2 are embedded above their legends.
  Unpublished manuscript prepared for peer review — please do not redistribute.
</footer>
"""
    out = f"{OUTDIR}/manuscript_review.html"
    open(out, "w", encoding="utf-8").write(page)
    print(f"  wrote {out}  ({os.path.getsize(out)/1024/1024:.2f} MB)")

if __name__ == "__main__":
    main()
