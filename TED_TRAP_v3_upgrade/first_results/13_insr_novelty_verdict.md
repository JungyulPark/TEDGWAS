# INSR Novelty Verification — Final Verdict

**Date**: 2026-04-23
**Method**: Full-text scan of Kim 2024 JCI Insight (via web_fetch) + Kim 2021 IOVS reference analysis
**Checked by**: Literature cross-reference scan

---

## Summary of Findings

### Kim 2024 JCI Insight (Kim DW et al., 2024;9(24):e182352)
**Full-text scanned.** Key findings related to INSR:

1. **INSR is NOT reported as a differentially expressed gene** in TED fibroblasts or adipocytes in main text or figures.
2. **Insulin receptor is mentioned ONLY in reference to linsitinib's off-target binding**: "Considering linsitinib's known affinity for the insulin receptor ([ref 31])..." — this is pharmacological context, not expression data.
3. "Insulin response" pathway appears in GO enrichment (Figure 5D) for TED adipocytes — pathway-level, not individual INSR expression.
4. **Main IGF-related DEGs reported**: IGF1, IGF1R (in fibroblasts only), IGFBP5 (up), IGFBP6 (down), IGFBP2 (down), IGFL4 (down).
5. Supplementary tables referenced in Kim 2024 do **not** mention INSR in the main text narrative.

### Kim 2021 IOVS (Kim DW et al., 2021;62(9):24) — the prior bulk RNA-seq paper
Cited as reference #18 by Kim 2024. Based on Kim 2024's own narrative citations:
- Focus was on **adipogenic transition** in TED orbital fibroblasts
- Primary findings cited: IGF1/IGF1R upregulation during adipogenesis
- No INSR-specific claim was propagated into the Kim 2024 citation trail
- **Requires direct supplementary Table S2 verification** for full confirmation

---

## Verdict: **UPHOLD with caveat**

Based on thorough scan of Kim 2024 full text (the newer, more comprehensive snRNA-Seq study):

> The v3 manuscript CAN report INSR upregulation as an observation that has not been individually highlighted in the dominant TED transcriptomic literature. However, the framing should avoid absolute "first report" language because:
> 1. Kim 2021 IOVS supplementary Table 2 has not been directly verified
> 2. Pathway-level "insulin response" enrichment IS reported by Kim 2024
> 3. Linsitinib's known INSR binding has been discussed at a mechanistic level

---

## Recommended Manuscript Language

### ✅ APPROVED wording (use this)

> *"We observed upregulation of the insulin receptor (INSR) transcript in TED orbital adipose tissue. While prior TED transcriptomic studies have documented enrichment of insulin-response pathways at the level of Gene Ontology (Kim et al., 2024, JCI Insight) and have emphasized the pharmacological cross-reactivity of IGF-1R inhibitors with INSR (Mulvihill et al., 2009), direct transcript-level INSR elevation in orbital adipocytes has not been a primary reported finding in prior published analyses to our knowledge. This observation provides one additional data point potentially relevant to the insulin-signaling off-target hypothesis for IGF-1R antagonism."*

### ❌ AVOID these phrasings

- "First report of INSR in TED"
- "Novel discovery of INSR upregulation"
- "Previously undescribed INSR finding"

### 🟡 CAVEAT to keep in Discussion limitations

> *"A systematic supplementary-table-level comparison with all prior TED transcriptomic datasets was not performed; this framing is provisional pending such verification."*

---

## Key references verified against

| Reference | DOI | Verified? | INSR mentioned? |
|-----------|-----|-----------|------------------|
| Kim DW et al. 2024 JCI Insight | 10.1172/jci.insight.182352 | ✅ Full text | Only in pharmacology context (linsitinib off-target) |
| Kim DW et al. 2021 IOVS | 10.1167/iovs.62.9.24 | ⚠️ Only via Kim 2024 citation trail | Unknown in supp Table 2 |
| Mulvihill et al. 2009 Future Med Chem | 10.4155/fmc.09.89 | Not fetched | Yes — this is the paper on INSR/IGF1R dual inhibition by OSI-906 |

---

## Action Items

- [x] Scanned Kim 2024 full text — no INSR DEG finding
- [ ] User should download Kim 2021 IOVS Supp Table 2 (https://iovs.arvojournals.org/article.aspx?articleid=2777099) and grep for "INSR"
- [ ] Update manuscript Discussion to use the APPROVED wording above
- [ ] Add Kim 2024 reference to support insulin-response pathway enrichment as prior finding
