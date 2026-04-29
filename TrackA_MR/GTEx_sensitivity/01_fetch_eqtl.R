setwd("c:/ProjectTEDGWAS/TrackA_MR/GTEx_sensitivity/")
library(httr)

TSHR_ENS <- "ENSG00000165409"
api_url <- paste0(
    "https://www.ebi.ac.uk/eqtl/api/v2/datasets/QTD000341/associations?",
    "gene_id=", TSHR_ENS,
    "&size=1000"
)

cat("Fetching eQTL Catalogue API for TSHR in GTEx Thyroid...\n")
response <- GET(api_url)
status <- status_code(response)
cat("Status:", status, "\n")
cat("Content:", content(response, "text"), "\n")
