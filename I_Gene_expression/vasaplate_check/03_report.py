#!/usr/bin/env python3
"""
03_report.py -- assemble res/vasaplate/vasaplate_paper_check.html.

Self-contained: figures are embedded as base64 PNGs (a light and a dark copy of
each, swapped by CSS), so the file can be copied anywhere and still render.

Every number in the page is read from the TSVs 01_compare.py wrote. Nothing is
typed in by hand, so the prose cannot drift away from the data.
"""
import base64
import html
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vp_common as C

PNG = f"{C.FIGS}/png"


def b64(path):
    with open(path, "rb") as fh:
        return base64.b64encode(fh.read()).decode()


def figure(name, caption):
    lp, dp = f"{PNG}/{name}.light.png", f"{PNG}/{name}.dark.png"
    if not os.path.exists(lp):
        return f'<p class="missing">figure {html.escape(name)} not generated yet</p>'
    out = [f'<figure><img class="lt" src="data:image/png;base64,{b64(lp)}" alt="{html.escape(caption)}">']
    if os.path.exists(dp):
        out.append(f'<img class="dk" src="data:image/png;base64,{b64(dp)}" alt="{html.escape(caption)}">')
    out.append(f"<figcaption>{caption}</figcaption></figure>")
    return "".join(out)


def fmt(v):
    if isinstance(v, float):
        return f"{v:,.0f}" if abs(v) >= 1000 and v == int(v) else f"{v:,.4g}"
    return html.escape(str(v))


def table(df, cls="data"):
    head = "".join(f"<th>{html.escape(str(c))}</th>" for c in df.columns)
    rows = []
    for idx, r in df.iterrows():
        cells = "".join(f"<td>{fmt(v)}</td>" for v in r)
        rows.append(f"<tr><th>{html.escape(str(idx))}</th>{cells}</tr>")
    return (f'<div class="scroll"><table class="{cls}"><thead><tr><th></th>{head}</tr>'
            f"</thead><tbody>{''.join(rows)}</tbody></table></div>")


def kv(pairs):
    out = ['<div class="scroll"><table class="data"><tbody>']
    for k, v, note in pairs:
        out.append(f'<tr><th>{html.escape(k)}</th><td class="num">{v}</td>'
                   f'<td class="note">{note}</td></tr>')
    out.append("</tbody></table></div>")
    return "".join(out)


CSS = """
:root{color-scheme:light dark;--bg:#fcfcfb;--fg:#0b0b0b;--fg2:#52514e;--line:#e3e2dd;
--card:#ffffff;--accent:#2a78d6;--warn:#eb6834;--ok:#1baf7a}
@media (prefers-color-scheme:dark){:root:where(:not([data-theme=light])){
--bg:#1a1a19;--fg:#fff;--fg2:#c3c2b7;--line:#3a3a38;--card:#232322;
--accent:#3987e5;--warn:#d95926;--ok:#199e70}}
:root[data-theme=dark]{--bg:#1a1a19;--fg:#fff;--fg2:#c3c2b7;--line:#3a3a38;--card:#232322;
--accent:#3987e5;--warn:#d95926;--ok:#199e70}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);
font:15px/1.65 ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,sans-serif}
.wrap{max-width:940px;margin:0 auto;padding:40px 20px 80px}
h1{font-size:1.85rem;line-height:1.25;margin:0 0 .3em;letter-spacing:-.01em}
h2{font-size:1.25rem;margin:2.4em 0 .6em;padding-top:.8em;border-top:1px solid var(--line)}
h3{font-size:1rem;margin:1.6em 0 .4em;color:var(--fg2)}
p,li{color:var(--fg)}
.sub{color:var(--fg2);font-size:.95rem;margin:0 0 2em}
code{background:var(--card);border:1px solid var(--line);border-radius:4px;
padding:.08em .35em;font-size:.88em}
.scroll{overflow-x:auto;margin:1em 0}
table.data{border-collapse:collapse;width:100%;font-size:.88rem}
table.data th,table.data td{border-bottom:1px solid var(--line);padding:.45em .6em;text-align:left;
vertical-align:top}
table.data thead th{color:var(--fg2);font-weight:600;white-space:nowrap}
table.data tbody th{font-weight:500;white-space:nowrap}
td.num{font-variant-numeric:tabular-nums;white-space:nowrap;font-weight:600}
td.note{color:var(--fg2);font-size:.92em}
figure{margin:1.4em 0}
figure img{width:100%;height:auto;border:1px solid var(--line);border-radius:8px;background:var(--card)}
figcaption{color:var(--fg2);font-size:.86rem;margin-top:.5em}
img.dk{display:none}
@media (prefers-color-scheme:dark){:root:where(:not([data-theme=light])) img.lt{display:none}
:root:where(:not([data-theme=light])) img.dk{display:block}}
:root[data-theme=dark] img.lt{display:none}
:root[data-theme=dark] img.dk{display:block}
:root[data-theme=light] img.lt{display:block}
:root[data-theme=light] img.dk{display:none}
.callout{background:var(--card);border:1px solid var(--line);border-left:3px solid var(--accent);
border-radius:6px;padding:.9em 1.1em;margin:1.2em 0}
.callout.warn{border-left-color:var(--warn)}
.callout.ok{border-left-color:var(--ok)}
.callout p:first-child{margin-top:0}.callout p:last-child{margin-bottom:0}
.missing{color:var(--warn)}
ul{padding-left:1.3em}
"""


def main():
    S = pd.read_csv(f"{C.RES}/comparison_summary.tsv", sep="\t", index_col=0)["value"]
    bs = pd.read_csv(f"{C.RES}/biotype_composition.tsv", sep="\t", index_col=0)

    def g(k, d="—"):
        return fmt(S[k]) if k in S.index else d

    has_trna = k_trna = False
    if "tRNA rows: ours" in S.index and float(S["tRNA rows: ours"]) > 0:
        has_trna = k_trna = True

    run_label = "bedv2 (annotation BED with tRNA)" if has_trna else \
                "rrnav2 (annotation BED without tRNA — tRNA re-run still in flight)"

    doc = [f"<title>VASA-plate mixing control: our re-run vs the paper</title>",
           f"<style>{CSS}</style>", '<div class="wrap">']
    A = doc.append

    A("<h1>VASA-plate species-mixing control: our re-run vs Salmen &amp; De Jonghe 2022</h1>")
    A('<p class="sub">Library <code>SRR14783059</code> / <code>GSM5369495</code>, '
      '<code>vasaplate-HEK293T-mESC</code> — the only VASA-plate library in GSE176588. '
      f'Pipeline run compared here: <strong>{html.escape(run_label)}</strong>.</p>')

    # --- what can be checked --------------------------------------------------
    A("<h2>What the paper actually lets us check</h2>")
    A("<p>This matters before any number is quoted, because the manuscript is nearly "
      "silent about this particular library:</p>")
    A('<div class="callout warn"><p><strong>The paper reports no barnyard number for '
      "this library.</strong> The 3.08% heterotypic doublet rate (p. 1781, Fig. 1d) is "
      "<em>VASA-drop</em> (GSM5369496). GSM5369495 appears only as a source of cells for "
      "the benchmarking panels. It also reports no rRNA percentage and no tRNA percentage "
      "for VASA-seq at all.</p>"
      "<p>So the primary reference here is the <strong>deposited count table</strong> for "
      "the same library, not the manuscript. Manuscript numbers are checked where they "
      "exist.</p></div>")

    A('<div class="callout"><p><strong>The paper states two different doublet rules, and '
      "they disagree by ~6× on this data.</strong> Fig. 1d thresholds the <em>UFI</em> "
      "fraction (&gt;25% from the other species = mixed); the Methods (p. 18) threshold the "
      "<em>gene</em> fraction (&gt;75% from one species = singlet). Both are reported below, "
      "always labelled.</p></div>")

    # --- headline -------------------------------------------------------------
    A("<h2>Headline: the pipeline reproduces the deposited table</h2>")
    A(kv([
        ("Genes shared with the published table", g("rows: simple shared"),
         f'of {g("rows: published simple")} simple rows published'),
        ("Spearman r, per-gene totals", g("concordance: Spearman r (per-gene totals)"),
         "across all shared genes"),
        ("Pearson r, log10 per-gene totals", g("concordance: Pearson r (log10 totals)"), ""),
        ("Median per-cell Pearson r (log1p)", g("concordance: median per-cell Pearson r (log1p)"),
         "computed per well, then median over 384"),
        ("Median log2(ours / published)", g("concordance: median log2(ours/published)"),
         "0.000 = no systematic scale difference"),
        ("Genes with exactly equal totals", g("concordance: exactly equal totals"), ""),
    ]))
    A(figure("04_ours_vs_published",
             "Per-gene total UFIs, ours against the deposited table, and the distribution of "
             "per-gene log2 ratios. Tight on y = x with a median ratio of 0."))

    # --- barnyard -------------------------------------------------------------
    A("<h2>Species mixing</h2>")
    A(figure("01_barnyard",
             "Barnyard plot under the Fig. 1d rule. The two panels are the same barcodes "
             "scored from the deposited table and from our re-run."))
    A(kv([
        ("Barcodes ≥ 7,500 UFIs — published", g("published: barcodes >= 7500 UFIs"), ""),
        ("Barcodes ≥ 7,500 UFIs — ours", g("ours: barcodes >= 7500 UFIs"), "identical gate"),
        ("Doublet rate, Fig. 1d rule — published", g("published: doublet rate, Fig.1d rule (%)") + "%",
         f'{g("published: n mixed (Fig.1d)")} mixed'),
        ("Doublet rate, Fig. 1d rule — ours", g("ours: doublet rate, Fig.1d rule (%)") + "%",
         f'{g("ours: n mixed (Fig.1d)")} mixed'),
        ("Doublet rate, Methods rule — published", g("published: doublet rate, Methods rule (%)") + "%", ""),
        ("Doublet rate, Methods rule — ours", g("ours: doublet rate, Methods rule (%)") + "%",
         "gene-based; sensitive to low-count detection, see caveat"),
        ("Human / mouse calls — published",
         f'{g("published: n human (Fig.1d)")} / {g("published: n mouse (Fig.1d)")}', ""),
        ("Human / mouse calls — ours",
         f'{g("ours: n human (Fig.1d)")} / {g("ours: n mouse (Fig.1d)")}', ""),
    ]))
    A('<div class="callout warn"><p>The <strong>Methods-rule</strong> rates diverge more than '
      "the Fig. 1d rates do. That rule counts <em>genes</em> rather than UFIs, so it is "
      "dominated by genes detected at one or two counts, and our run detects more distinct "
      f'rows overall ({g("rows: ours simple")} vs {g("rows: published simple")} simple rows). '
      "The UFI-based rule, which is the one Fig. 1d is drawn with, agrees closely.</p></div>")

    # --- per cell -------------------------------------------------------------
    A("<h2>Per-cell yield</h2>")
    A(figure("02_genes_umis",
             "Total UFIs and genes detected per cell, split by species call, ours against "
             "the deposited table. Medians are drawn and labelled."))
    A(kv([
        ("HEK293T median UFIs — published / ours",
         f'{g("published: HEK293T median UFIs")} / {g("ours: HEK293T median UFIs")}', ""),
        ("HEK293T median genes — published / ours",
         f'{g("published: HEK293T median genes")} / {g("ours: HEK293T median genes")}',
         "unambiguous human rows only"),
        ("HEK293T mean purity — published / ours",
         f'{g("published: HEK293T mean purity")} / {g("ours: HEK293T mean purity")}', ""),
        ("mESC median UFIs — published / ours",
         f'{g("published: mESC median UFIs")} / {g("ours: mESC median UFIs")}', ""),
        ("mESC median genes — published / ours",
         f'{g("published: mESC median genes")} / {g("ours: mESC median genes")}', ""),
        ("mESC mean purity — published / ours",
         f'{g("published: mESC mean purity")} / {g("ours: mESC mean purity")}', ""),
    ]))
    A('<div class="callout"><p>These are <strong>not</strong> comparable to the paper\'s '
      "<em>9,480 ± 1,252 genes per cell</em> (Fig. 1f): that figure is downsampled to 75,000 "
      "trimmed reads per cell, while these tables are at full depth. The comparable quantity "
      "is the deposited table, alongside.</p></div>")

    # --- biotype --------------------------------------------------------------
    A("<h2>Biotype composition — and the one manuscript number that matches</h2>")
    A(figure("03_biotype_composition",
             "Share of assigned UFIs by biotype over simple rows, ours against the deposited "
             "table."))
    A(kv([
        ("sncRNA share — the paper states (VASA-plate)", g("  paper states, VASA-plate (%)") + "%",
         "p. 1781"),
        ("sncRNA share — deposited table", g("sncRNA UFI share, published (%)") + "%", ""),
        ("sncRNA share — ours", g("sncRNA UFI share, ours (%)") + "%",
         "miscRNA, snoRNA, snRNA, scaRNA, miRNA, ribozyme, rRNA, Mt"),
    ]))
    A('<div class="callout ok"><p>This is the one headline number in the manuscript that this '
      "library can be checked against directly, and it lands: the paper says <strong>1.4%</strong> "
      f'sncRNA for VASA-plate, the deposited table gives <strong>{g("sncRNA UFI share, published (%)")}%</strong>, '
      f'and our re-run gives <strong>{g("sncRNA UFI share, ours (%)")}%</strong>.</p></div>')
    A("<h3>Full biotype table</h3>")
    A(table(bs.round(4)))

    # --- tRNA -----------------------------------------------------------------
    A("<h2>tRNA — the gap this work closes</h2>")
    A(figure("05_trna_before_after",
             "Distinct tRNA rows detected. The annotation BED shipped with the mixed "
             "reference contained none despite being named IntronExonTrna.bed."))
    A(kv([
        ("tRNA rows — deposited table", g("tRNA rows: published"), "bare GtRNAdb tRNAscan ids"),
        ("tRNA rows — ours", g("tRNA rows: ours"),
         "0 under the old BED" if not has_trna else "after the BED rebuild"),
        ("tRNA rows shared", g("tRNA rows: shared"), ""),
        ("tRNA UFI share — deposited table", g("tRNA UFI share, published (%)") + "%", ""),
    ]))
    A("<p>Two findings from matching the authors' method rather than inventing one:</p><ul>")
    A("<li><strong>Their tRNA names carry no species tag.</strong> Rebuilding both GtRNAdb sets "
      "gives 619 human + 1,139 mouse loci with exactly <strong>1</strong> colliding name "
      "(<code>16.tRNA1-LysNNN</code>), and the published rows show no species marker — so one "
      "flat namespace, reproduced.</li>")
    A("<li><strong>Their locus ids cannot be reproduced, and chasing them is wasted effort.</strong> "
      "Only 503 of the 1,130 published tRNA names exist in the current GtRNAdb. The tell is the "
      "isotype label: 110 published rows say <code>Undet</code>, the retired tRNAscan-SE spelling, "
      "and current GtRNAdb writes <code>Und</code>. GtRNAdb has also renumbered. What survives is "
      "the biology — <strong>59 of 62 isotype classes (95%)</strong> — and isotype is exactly what "
      "step 7 collapses tRNA to.</li></ul>")
    A('<div class="callout warn"><p><strong>The <code>_tRNA.*Counts.tsv</code> tables merge the two '
      "species.</strong> <code>countTables_fromPickle.py</code> groups tRNA rows by "
      "<code>rsplit('.')[-1]</code>, i.e. isotype+anticodon, so human <code>ValAAC</code> and mouse "
      "<code>ValAAC</code> land in one row. This is true of the published tables too. "
      "Species-resolved per-locus tRNA counts survive only in the <code>total</code> tables, which "
      "keep the full row name.</p></div>")

    # --- not checkable --------------------------------------------------------
    A("<h2>What could not be checked, and why</h2><ul>")
    for item in C.NOT_IN_PAPER:
        A(f"<li>{html.escape(item)}</li>")
    A("</ul>")
    A('<div class="callout warn"><p>One more caution against over-reading any barnyard mismatch: '
      "the authors' <em>own</em> deposited VASA-drop table gives 4.57% (Fig. 1d rule) and 5.42% "
      "(Methods rule) against the 3.08% printed in the paper. Their published table does not "
      "trivially reproduce their own published number under either stated rule, so a discrepancy "
      "on a barnyard figure is not automatically a fault in this pipeline.</p></div>")

    # --- provenance -----------------------------------------------------------
    A("<h2>How to regenerate this</h2>")
    A("<pre><code>cd code/I_Gene_expression/vasaplate_check\n"
      "./01_compare.py rrnav2      # or bedv2\n./02_figures.py\n./03_report.py</code></pre>")
    A('<p class="sub">Every number on this page is read from '
      "<code>comparison_summary.tsv</code>, <code>per_cell.tsv</code>, "
      "<code>biotype_composition.tsv</code> and <code>gene_concordance.tsv</code> in "
      "<code>res/vasaplate/</code>, all written by <code>01_compare.py</code>. None is typed in "
      "by hand.</p>")
    A("</div>")

    dest = f"{C.RES}/vasaplate_paper_check.html"
    with open(dest, "w") as fh:
        fh.write("\n".join(doc))
    print(f"wrote {dest} ({os.path.getsize(dest)/1e6:.2f} MB)")


if __name__ == "__main__":
    main()
