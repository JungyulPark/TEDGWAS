#!/usr/bin/env Rscript
# Frequency-only (paired variants) and full action=2 re-harmonization analyses.
# Requires a normalized, provenance-verified eQTLGen AF extract; no imputation.
suppressPackageStartupMessages({library(data.table);library(TwoSampleMR)})
setDTthreads(2)
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==2)
root<-normalizePath(args[1],winslash="/",mustWork=TRUE)
afpath<-normalizePath(args[2],winslash="/",mustWork=TRUE)
out<-file.path(root,"sensitivity_v1")
if(dir.exists(out)) stop("Versioned sensitivity output exists; do not overwrite.")
af<-fread(afpath)
stopifnot(all(c("snp","allele_a","allele_b","af_b","genome_build") %in% names(af)),!anyDuplicated(af$snp),all(af$genome_build=="GRCh37"))
dir.create(out)
logfile<-file.path(out,"run.log")
logmsg<-function(...){s<-paste(format(Sys.time(),tz="UTC",usetz=TRUE),...);cat(s,"\n");cat(s,"\n",file=logfile,append=TRUE)}
options(error=function(){cat("FAILED; inspect console error.\n",file=logfile,append=TRUE);q(status=1)})
logmsg("eQTLGen AF sensitivity started; input:",afpath)
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
ex<-fread(file.path(root,"baseline_v1/exposures_reference.csv"))
stopifnot(nrow(ex)==6135,uniqueN(ex$gene_symbol)==2544)
setnames(ex,c("eaf","se","beta"),c("eaf_ref","se_ref","beta_ref"))
ex<-merge(ex,af,by="snp",all.x=TRUE,suffixes=c("","_af"))
ex[, pair_ok:=(effect_allele==allele_b&other_allele==allele_a)|(effect_allele==allele_a&other_allele==allele_b)]
ex[, eaf_new:=fifelse(pair_ok,fifelse(effect_allele==allele_b,af_b,1-af_b),NA_real_)]
ex[, eligible_af:=!is.na(eaf_new)&eaf_new>0&eaf_new<1]
ex[, se_new:=fifelse(eligible_af,1/sqrt(2*eaf_new*(1-eaf_new)*(n_samples+zscore^2)),NA_real_)]
ex[, beta_new:=zscore*se_new]
ex[, scale_new_ref:=se_new/se_ref]
ex[, delta_z:=beta_new/se_new-beta_ref/se_ref]
stopifnot(max(abs(ex$delta_z),na.rm=TRUE)<1e-8)
fwrite(ex,file.path(out,"exposure_frequency_comparison.csv"))
logmsg("Usable AF rows",sum(ex$eligible_af),"/",nrow(ex),"max |delta Z|",max(abs(ex$delta_z),na.rm=TRUE))
baseline_h<-readRDS(file.path(root,"baseline_v1/harmonized_all.rds"))
ocs<-c(BBJ="BBJ_Graves",UKB="UKB_hyperthyroid",FinnGen="FinnGen_GO")
outcomes<-lapply(names(ocs),function(nm){
 d<-fread(file.path(root,"inputs",paste0(nm,"_selected.tsv")))
 d<-d[!is.na(snp)&snp!=""]
 if(nm=="UKB"){mu<-3731/484598;d[,c("beta","se"):=.(beta/(mu*(1-mu)),se/(mu*(1-mu)))]}
 d[,.(SNP=snp,effect_allele.outcome=toupper(ea),other_allele.outcome=toupper(oa),beta.outcome=beta,se.outcome=se,pval.outcome=pvalue,eaf.outcome=eaf)]
});names(outcomes)<-unname(ocs)
fit<-function(h,g,oc,scenario){
 n<-sum(h$mr_keep)
 if(n==0) return(data.table(gene_symbol=g,outcome=oc,scenario=scenario,method=NA_character_,n_iv=0L,beta=NA_real_,se=NA_real_,pvalue=NA_real_))
 r<-suppressMessages(mr(as.data.frame(h),method_list=if(n==1) "mr_wald_ratio" else "mr_ivw"))
 stopifnot(nrow(r)==1,is.finite(r$b),is.finite(r$se),r$se>0)
 data.table(gene_symbol=g,outcome=oc,scenario=scenario,method=r$method,n_iv=r$nsnp,beta=r$b,se=r$se,pvalue=r$pval)
}
results<-list();harmonized<-list();qc<-list();warnings_seen<-character();genes<-unique(ex$gene_symbol)
for(i in seq_along(genes)){
 g<-genes[i];eg<-ex[gene_symbol==g];eg_valid<-eg[eligible_af==TRUE]
 for(oc in names(outcomes)){
   bh<-copy(baseline_h[id.exposure==g&id.outcome==oc])
   common<-bh[mr_keep&SNP %in% eg_valid$snp]
   results[[length(results)+1]]<-fit(common,g,oc,"paired_reference")
   nh<-copy(common)
   scales<-eg_valid$scale_new_ref[match(nh$SNP,eg_valid$snp)]
   nh[,beta.exposure:=beta.exposure*scales];nh[,se.exposure:=se.exposure*scales]
   results[[length(results)+1]]<-fit(nh,g,oc,"paired_eqtlgen")
   ef<-data.frame(SNP=eg_valid$snp,beta.exposure=eg_valid$beta_new,se.exposure=eg_valid$se_new,effect_allele.exposure=eg_valid$effect_allele,other_allele.exposure=eg_valid$other_allele,eaf.exposure=eg_valid$eaf_new,pval.exposure=eg_valid$pvalue,exposure=rep(g,nrow(eg_valid)),id.exposure=rep(g,nrow(eg_valid)))
   of<-as.data.frame(outcomes[[oc]][SNP %in% eg_valid$snp])
   fh<-bh[0]
   if(nrow(ef)>0&nrow(of)>0){
     of$outcome<-oc;of$id.outcome<-oc
     fh<-as.data.table(withCallingHandlers(suppressMessages(harmonise_data(ef,of,action=2)),warning=function(w){warnings_seen<<-c(warnings_seen,paste(g,oc,conditionMessage(w)));invokeRestart("muffleWarning")}))
   }
   results[[length(results)+1]]<-fit(fh,g,oc,"reharmonized_eqtlgen")
   harmonized[[length(harmonized)+1]]<-fh
   comparison<-merge(bh[,.(SNP,keep_ref=mr_keep,palindromic_ref=palindromic,ambiguous_ref=ambiguous,beta_out_ref=beta.outcome)],fh[,.(SNP,keep_new=mr_keep,palindromic_new=palindromic,ambiguous_new=ambiguous,beta_out_new=beta.outcome)],by="SNP",all=TRUE)
   comparison[,`:=`(gene_symbol=g,outcome=oc)]
   qc[[length(qc)+1]]<-comparison
 }
 if(i%%250==0|i==length(genes))logmsg("Genes",i,"/",length(genes))
}
res<-rbindlist(results);fwrite(res,file.path(out,"primary_MR_all_scenarios.csv"))
fwrite(res[gene_symbol %in% c("TSHR","IGF1R","CTLA4")],file.path(out,"anchor_MR_comparison.csv"))
saveRDS(rbindlist(harmonized,fill=TRUE),file.path(out,"harmonized_eqtlgen.rds"))
fwrite(rbindlist(qc,fill=TRUE),file.path(out,"harmonization_variant_comparison.csv"))
writeLines(unique(warnings_seen),file.path(out,"warnings.txt"))
paired<-merge(res[scenario=="paired_reference"],res[scenario=="paired_eqtlgen"],by=c("gene_symbol","outcome"),suffixes=c("_ref","_new"))
stopifnot(all(paired$n_iv_ref==paired$n_iv_new))
single<-paired[n_iv_ref==1]
stopifnot(all(sign(single$beta_ref)==sign(single$beta_new)),max(abs(single$pvalue_ref-single$pvalue_new))<1e-8)
fwrite(paired,file.path(out,"paired_MR_comparison.csv"))
baseline<-fread(file.path(root,"baseline_v1/primary_MR_reproduced.csv"));baseline[,scenario:="original_reference"]
power_input<-rbindlist(list(baseline,res),fill=TRUE)[n_iv>0]
power<-power_input[,{
 alpha<-if(outcome[1]=="BBJ_Graves")0.05/2544 else .05
 cutoff<-(qnorm(1-alpha/2)+qnorm(.8))*se
 qu<-quantile(cutoff,c(.25,.5,.75))
 list(n_genes=.N,alpha=alpha,or_q1=exp(qu[1]),or_median=exp(qu[2]),or_q3=exp(qu[3]),frac_OR1_5=mean(cutoff<=log(1.5)),frac_OR2=mean(cutoff<=log(2)),frac_OR3=mean(cutoff<=log(3)))
},by=.(scenario,outcome)]
fwrite(power,file.path(out,"power_all_scenarios.csv"))
logmsg("Completed MR sensitivity; same-variant single-IV direction/P invariance passed.")
