# =============================================================================
# Task 2 — Quantify the effect of the wrong FinnGen N on coloc.abf posteriors
# =============================================================================
# Question : taskE_01_coloc.R used N = 520,387 (an unreplaced placeholder with
#            no traceable source) instead of the correct N = 500,348.
#            How much does that move PP.H4?
#
# Answer   : essentially not at all (|dPP.H4| ~ 1e-08), because coloc's
#            case-control variance approximation depends on N and s only
#            through the product N*s*(1-s), and both parameterizations encode
#            the same case count (858).
#
# Usage    : Rscript internal/scripts/verify_coloc_sN_impact.R
#
# Result recorded in internal/FINNGEN_endpoint_provenance.md section 6.
# =============================================================================

suppressMessages(library(coloc))

cat("coloc version:", as.character(packageVersion("coloc")), "\n\n")

# ---- 1. The variance approximation actually used -----------------------------
cat("=== Var.data.cc (from installed package source) ===\n")
print(coloc:::Var.data.cc)
cat("\nN and s enter only via the product N*s*(1-s).\n\n")

# ---- 2. Parameterizations ----------------------------------------------------
CASES <- 858
used_N <- 520387; used_s <- CASES / used_N   # taskE_01_coloc.R:76-77 (placeholder)
tru_N  <- 500348; tru_s  <- CASES / tru_N    # FinnGen R12 GRAVES_OPHT (confirmed)

a <- used_N * used_s * (1 - used_s)
b <- tru_N  * tru_s  * (1 - tru_s)

cat("=== N*s*(1-s) ===\n")
cat(sprintf("  used    : N=%d  s=%.8f  -> %.4f\n", used_N, used_s, a))
cat(sprintf("  correct : N=%d  s=%.8f  -> %.4f\n", tru_N,  tru_s,  b))
cat(sprintf("  ratio   : %.8f   (V differs by %.4f%%)\n\n", a / b, (a / b - 1) * 100))

# ---- 3. Empirical check on a synthetic shared-causal locus --------------------
cat("=== Empirical coloc.abf comparison ===\n")
set.seed(42)
n   <- 1500
maf <- runif(n, 0.05, 0.5)
z_e <- rnorm(n)
z_g <- rnorm(n)

causal <- 400
z_e[causal] <- 9.5
z_g[causal] <- 7.2
for (j in (causal - 5):(causal + 5)) {      # crude LD block around the causal SNP
  z_e[j] <- z_e[causal] * 0.9
  z_g[j] <- z_g[causal] * 0.9
}

p_e <- 2 * pnorm(-abs(z_e))
p_g <- 2 * pnorm(-abs(z_g))
snp <- paste0("rs", seq_len(n))

d_eqtl <- list(snp = snp, pvalues = p_e, type = "quant", N = 31684, MAF = maf)

run <- function(N, s) {
  d_gwas <- list(snp = snp, pvalues = p_g, type = "cc", s = s, N = N, MAF = maf)
  invisible(capture.output(r <- coloc.abf(d_gwas, d_eqtl)))
  r$summary
}

r_used <- run(used_N, used_s)
r_true <- run(tru_N,  tru_s)

cmp <- data.frame(
  hypothesis     = names(r_used)[-1],
  used_520387    = as.numeric(r_used)[-1],
  correct_500348 = as.numeric(r_true)[-1]
)
cmp$abs_diff <- cmp$correct_500348 - cmp$used_520387
print(cmp, digits = 10, row.names = FALSE)

cat(sprintf("\nPP.H4 absolute difference = %.3e\n",
            cmp$abs_diff[cmp$hypothesis == "PP.H4.abf"]))
cat("\nConclusion: provenance defect, not a numerical one. The reported\n")
cat("PP.H4 values are unchanged at the precision reported (3 decimals).\n")
cat("This says NOTHING about the other Task 4 defects (p-value-based ABF,\n")
cat("EUR MAF applied to the East Asian BBJ outcome, missing UKB coloc).\n")
