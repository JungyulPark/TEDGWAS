# 1. Replication scale 확인
library(ieugwasr)
token <- Sys.getenv("OPENGWAS_JWT")
if (nchar(token) < 100) stop("Token not loaded.")
rep_info <- gwasinfo(id = "ebi-a-GCST90038636", opengwas_jwt = token)
pri_info <- gwasinfo(id = "ebi-a-GCST90018627", opengwas_jwt = token)
print(rbind(rep_info, pri_info))

# 2. Sensitivity 결과 요약
sens <- read.csv("TrackA_MR/results/02_mr_sensitivity.csv")
primary_sens <- sens[sens$outcome_role == "Primary", ]
print(primary_sens, row.names = FALSE)
