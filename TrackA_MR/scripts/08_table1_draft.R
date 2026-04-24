# =============================================================================
# TED-TRAP — Table 1 Draft Generator
# Study design: data sources, cohorts, candidate genes
# =============================================================================

suppressPackageStartupMessages({
    library(dplyr)
})

have_officer <- requireNamespace("officer", quietly = TRUE) &&
    requireNamespace("flextable", quietly = TRUE)

if (!dir.exists("TrackA_MR/manuscript")) dir.create("TrackA_MR/manuscript", recursive = TRUE)

# =============================================================================
# SECTION A — GWAS Summary Statistics
# =============================================================================

section_A <- data.frame(
    Role = c(
        "Exposure (cis-eQTL)",
        "Primary outcome",
        "Replication outcome",
        "Sensitivity outcome",
        "MVMR Exposure 2"
    ),
    Dataset = c(
        "eQTLGen (blood cis-eQTL)",
        "Graves disease (BBJ)",
        "Hyperthyroidism (UKB)",
        "Graves ophthalmopathy (FinnGen R12)",
        "Graves disease, strict (FinnGen R12)"
    ),
    Ancestry = c(
        "European (87.5%)", "East Asian", "European",
        "European (Finnish)", "European (Finnish)"
    ),
    Sample_Size = c(
        "n = 31,684",
        "n = 212,453 (2,809 cases)",
        "n = 484,598 (3,731 cases)",
        "n = 500,348 (858 cases)",
        "n = 500,348 (3,962 cases)"
    ),
    Access = c(
        "OpenGWAS: eqtl-a-*",
        "OpenGWAS: ebi-a-GCST90018627",
        "OpenGWAS: ebi-a-GCST90038636",
        "FinnGen R12 public download",
        "FinnGen R12 public download"
    ),
    Reference = c(
        "Võsa et al., 2021 Nat Genet",
        "Sakaue et al., 2021 Nat Genet",
        "Dönertaş et al., 2021 Nat Aging",
        "Kurki et al., 2023 Nature",
        "Kurki et al., 2023 Nature"
    ),
    stringsAsFactors = FALSE
)

cat("=== Section A: GWAS Data Sources ===\n")
print(section_A, row.names = FALSE)

# =============================================================================
# SECTION B — Orbital Transcriptomic Cohort
# =============================================================================

section_B <- data.frame(
    Variable = c(
        "Total samples",
        "TED samples",
        "Control sample",
        "Sample type",
        "Sequencing",
        "Technical replication",
        "Ethics approval"
    ),
    Description = c(
        "n = 5 patients (10 libraries)",
        "n = 4 (inactive TED, orbital adipose tissue during rehabilitative surgery)",
        "n = 1 (nasal fat pad, non-TED blepharoplasty)",
        "Bulk RNA-seq from orbital adipose tissue",
        "Illumina platform; paired-end 100-150 bp",
        "2 technical replicates per patient (averaged for analysis)",
        "Catholic University of Korea IRB; written informed consent"
    ),
    stringsAsFactors = FALSE
)

cat("\n=== Section B: Orbital Transcriptomic Cohort ===\n")
print(section_B, row.names = FALSE)

# =============================================================================
# SECTION C — Pre-specified Candidate Gene Set
# =============================================================================

section_C <- data.frame(
    Group = c(
        "MR Primary (3 genes)",
        "MR Secondary (6 genes)",
        "Insulin Signaling Cassette (5 genes)",
        "TED-Specific Biology (6 genes)"
    ),
    Genes = c(
        "TSHR, IGF1R, IGF1",
        "ARRB1, PPARG, IRS1, AKT1, TNF, CTLA4",
        "INSR, IRS2, FOXO1, PIK3R1, PDPK1",
        "HAS1, HAS2, HAS3, ADIPOQ, FABP4, CEBPA"
    ),
    Rationale = c(
        "Primary drug targets (TSHR; IGF1R = teprotumumab target) + ligand",
        "Secondary MR candidates from prior TED literature and autoimmunity pathways",
        "Receptor + downstream adapters (IR/IGF-1R hybrid receptor + signaling)",
        "TED orbital pathology hallmarks: HA synthesis + adipogenesis"
    ),
    stringsAsFactors = FALSE
)

cat("\n=== Section C: Pre-specified Candidate Genes (20 total) ===\n")
print(section_C, row.names = FALSE)

# =============================================================================
# SECTION D — Analysis Framework
# =============================================================================

section_D <- data.frame(
    Endpoint = c(
        "Primary EP",
        "SE1 — Mediation",
        "SE2 — Off-target",
        "SE3 — Tissue validation",
        "SE4 — Colocalization",
        "SE5 — Triangulation"
    ),
    Description = c(
        "Cis-MR causal effect: TSHR vs IGF1R on Graves disease",
        "MVMR of TSHR on TED, adjusted for GD liability",
        "IGF1R → glucose/HbA1c/T2D/hearing",
        "Pre-specified 20-gene DESeq2 in orbital tissue",
        "Shared causal variant at TSHR locus (eQTL × GD)",
        "Integration of genetic + locus + tissue evidence"
    ),
    Method = c(
        "Inverse variance weighted + 4 sensitivity",
        "MVMR-IVW (Sanderson 2019)",
        "IVW + MR Egger + Weighted median + mode",
        "DESeq2 raw + TPM; binomial direction test",
        "coloc.abf (Giambartolomei 2014)",
        "Concordance matrix"
    ),
    Correction = c(
        "Bonferroni: α = 0.05 / 24 = 2.08×10⁻³",
        "—",
        "—",
        "BH within 20-gene set",
        "PP.H4 ≥ 0.80 threshold",
        "—"
    ),
    stringsAsFactors = FALSE
)

cat("\n=== Section D: Analysis Framework ===\n")
print(section_D, row.names = FALSE)

# =============================================================================
# SAVE CSVs
# =============================================================================

write.csv(section_A, "TrackA_MR/manuscript/Table1_SectionA_GWAS.csv", row.names = FALSE)
write.csv(section_B, "TrackA_MR/manuscript/Table1_SectionB_Orbital.csv", row.names = FALSE)
write.csv(section_C, "TrackA_MR/manuscript/Table1_SectionC_Candidates.csv", row.names = FALSE)
write.csv(section_D, "TrackA_MR/manuscript/Table1_SectionD_Framework.csv", row.names = FALSE)

cat("\n💾 Saved 4 CSV files to TrackA_MR/manuscript/\n")

# =============================================================================
# DOCX Export
# =============================================================================

if (have_officer) {
    library(officer)
    library(flextable)

    doc <- read_docx()
    doc <- body_add_par(doc, "Table 1. Study design, data sources, and pre-specified analysis framework.",
        style = "heading 1"
    )
    doc <- body_add_par(doc, "")

    # Section A
    doc <- body_add_par(doc, "A. GWAS summary statistics", style = "heading 2")
    doc <- body_add_flextable(
        doc,
        flextable(section_A) %>% autofit() %>% theme_booktabs()
    )
    doc <- body_add_par(doc, "")

    # Section B
    doc <- body_add_par(doc, "B. Orbital transcriptomic cohort", style = "heading 2")
    doc <- body_add_flextable(
        doc,
        flextable(section_B) %>% autofit() %>% theme_booktabs()
    )
    doc <- body_add_par(doc, "")

    # Section C
    doc <- body_add_par(doc, "C. Pre-specified candidate gene set (n = 20)", style = "heading 2")
    doc <- body_add_flextable(
        doc,
        flextable(section_C) %>% autofit() %>% theme_booktabs()
    )
    doc <- body_add_par(doc, "")

    # Section D
    doc <- body_add_par(doc, "D. Analysis framework", style = "heading 2")
    doc <- body_add_flextable(
        doc,
        flextable(section_D) %>% autofit() %>% theme_booktabs()
    )

    print(doc, target = "TrackA_MR/manuscript/Table1_draft.docx")
    cat("💾 Saved DOCX: TrackA_MR/manuscript/Table1_draft.docx\n")
}

cat("\n✅ Table 1 complete.\n")
