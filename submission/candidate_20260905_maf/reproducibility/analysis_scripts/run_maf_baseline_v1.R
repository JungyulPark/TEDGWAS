#!/usr/bin/env Rscript
# Reproduce the verified 2544-gene primary MR with unchanged reference frequency.
# Original files are read only. Re-harmonization uses TwoSampleMR action=2.
suppressPackageStartupMessages({library(data.table); library(TwoSampleMR)})
setDTthreads(2)
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
root <- normalizePath(args[1], winslash="/", mustWork=TRUE)
inp <- file.path(root,"inputs")
out <- file.path(root,"baseline_v1")
if (dir.exists(out)) stop("Versioned output already exists; do not overwrite.")
dir.create(out)
logfile <- file.path(out,"run.log")
logmsg <- function(...) {s <- paste(format(Sys.time(),tz="UTC",usetz=TRUE),...);cat(s,"\n");cat(s,"\n",file=logfile,append=TRUE)}
options(error=function(){cat("FAILED; see console error.\n",file=logfile,append=TRUE);q(status=1)})
logmsg("Baseline reproduction started.")
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
iv <- fread(file.path(inp,"instruments_verified.csv"))
fr <- fread(file.path(inp,"eur_freq_selected.tsv"))
stopifnot(nrow(iv)==6135,uniqueN(iv$gene_symbol)==2544,!anyDuplicated(fr$SNP))
ex <- merge(iv,fr[,.(snp=SNP,frq_A1=A1,frq_A2=A2,MAF)],by="snp",all.x=TRUE)
ex[, eaf:=fifelse(effect_allele==frq_A1,MAF,fifelse(effect_allele==frq_A2,1-MAF,NA_real_))]
ex[, se:=1/sqrt(2*eaf*(1-eaf)*(n_samples+zscore^2))]
ex[, beta:=zscore*se]
ex[, genome_build:="GRCh37"]
fwrite(ex,file.path(out,"exposures_reference.csv"))
outnames <- c(BBJ="BBJ_Graves",UKB="UKB_hyperthyroid",FinnGen="FinnGen_GO")
outcomes <- lapply(names(outnames),function(nm){
 d<-fread(file.path(inp,paste0(nm,"_selected.tsv")))
 d<-d[!is.na(snp)&snp!=""]
 if(nm=="UKB") {mu<-3731/484598;d[,c("beta","se"):=.(beta/(mu*(1-mu)),se/(mu*(1-mu)))]}
 d[,.(SNP=snp,effect_allele.outcome=toupper(ea),other_allele.outcome=toupper(oa),
       beta.outcome=beta,se.outcome=se,pval.outcome=pvalue,eaf.outcome=eaf)]
})
names(outcomes)<-unname(outnames)
results<-list(); harmonized<-list(); qc<-list(); warns<-character()
genes<-unique(iv$gene_symbol)
for(i in seq_along(genes)) {
 g<-genes[i]; eg<-ex[gene_symbol==g]
 ef<-data.frame(SNP=eg$snp,beta.exposure=eg$beta,se.exposure=eg$se,
   effect_allele.exposure=eg$effect_allele,other_allele.exposure=eg$other_allele,
   eaf.exposure=eg$eaf,pval.exposure=eg$pvalue,exposure=g,id.exposure=g)
 for(oc in names(outcomes)) {
   of<-as.data.frame(outcomes[[oc]][SNP %in% eg$snp]);nfound<-nrow(of)
   if(nfound==0) {qc[[length(qc)+1]]<-data.table(gene_symbol=g,outcome=oc,n_pre=nrow(eg),n_matched=0,n_kept=0,n_ambiguous=0,n_remove=0);next}
   of$outcome<-oc;of$id.outcome<-oc
   h<-withCallingHandlers(suppressMessages(harmonise_data(ef,of,action=2)),warning=function(w){warns<<-c(warns,paste(g,oc,conditionMessage(w)));invokeRestart("muffleWarning")})
   harmonized[[length(harmonized)+1]]<-as.data.table(h)
   nk<-sum(h$mr_keep)
   qc[[length(qc)+1]]<-data.table(gene_symbol=g,outcome=oc,n_pre=nrow(eg),n_matched=nfound,n_kept=nk,n_ambiguous=sum(h$ambiguous),n_remove=sum(h$remove))
   if(nk==0) next
   method<-if(nk==1) "mr_wald_ratio" else "mr_ivw"
   r<-suppressMessages(mr(h,method_list=method))
   stopifnot(nrow(r)==1,is.finite(r$b),is.finite(r$se),r$se>0)
   results[[length(results)+1]]<-data.table(gene_symbol=g,outcome=oc,method=r$method,n_iv=r$nsnp,beta=r$b,se=r$se,pvalue=r$pval)
 }
 if(i%%250==0 || i==length(genes)) logmsg("Genes",i,"/",length(genes),"primary estimates",length(results))
}
res<-rbindlist(results);hh<-rbindlist(harmonized,fill=TRUE)
fwrite(res,file.path(out,"primary_MR_reproduced.csv"))
saveRDS(hh,file.path(out,"harmonized_all.rds"))
fwrite(hh,file.path(out,"harmonized_all.csv"))
fwrite(rbindlist(qc),file.path(out,"harmonization_qc.csv"))
writeLines(unique(warns),file.path(out,"warnings.txt"))
canonical<-fread(file.path(dirname(root),"revision/submission/candidate_20260905/provenance/MR_primary_canonical.csv"))
cmp<-merge(res,canonical[,.(gene_symbol,outcome,method,n_iv,beta,se,pvalue)],by=c("gene_symbol","outcome"),all=TRUE,suffixes=c("_rerun","_canonical"))
cmp[, `:=`(delta_beta=beta_rerun-beta_canonical,delta_se=se_rerun-se_canonical,delta_p=pvalue_rerun-pvalue_canonical)]
cmp[, pass:=!is.na(method_rerun)&!is.na(method_canonical)&method_rerun==method_canonical&n_iv_rerun==n_iv_canonical&abs(delta_beta)<1e-8&abs(delta_se)<1e-8&abs(delta_p)<1e-8]
fwrite(cmp,file.path(out,"canonical_comparison.csv"))
summary<-data.table(test=c("primary_rows","paired_rows","all_estimates_match","max_abs_beta_diff","max_abs_se_diff","max_abs_p_diff"),
 value=c(nrow(res),nrow(cmp),all(cmp$pass),max(abs(cmp$delta_beta),na.rm=TRUE),max(abs(cmp$delta_se),na.rm=TRUE),max(abs(cmp$delta_p),na.rm=TRUE)))
fwrite(summary,file.path(out,"verification_summary.csv"));print(summary)
stopifnot(nrow(res)==7219,nrow(cmp)==7219,all(cmp$pass))
logmsg("PASS: all 7219 primary estimates and instrument counts reproduced.")
