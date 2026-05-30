library(data.table)
dt <- fread('c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/01_pqtl_check/results/TaskA_01_pqtl_availability_v1.csv')
dt[gene_symbol == 'IGF1' & platform == 'deCODE', c('availability_status', 'notes') := list('cis_pqtl_available', 'web search 기반 잠정; deCODE 패널 목록 원본 대조 필요')]
dt[gene_symbol == 'IGF1' & platform == 'UKBPPP', c('availability_status', 'notes') := list('not_measured', 'Olink Explore 3072 패널 부재 (web search 확인)')]
dt[gene_symbol %in% c('TSHR', 'IGF1R', 'INSR') & availability_status == 'pending', notes := '패널 목록 원본 확인 후 확정 요망']
fwrite(dt, 'c:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/01_pqtl_check/results/TaskA_01_pqtl_availability_v1.csv')
print(dt[, c('gene_symbol', 'platform', 'availability_status', 'notes')])
