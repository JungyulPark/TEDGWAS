#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(TwoSampleMR)})
setDTthreads(2)
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1)
t<-normalizePath(args[1],winslash="/",mustWork=TRUE);w<-dirname(t)
stopifnot(file.exists(file.path(t,"inputs_v2/extraction_manifest.json")))
out<-file.path(t,"mr_v1");if(dir.exists(out))stop("Version exists; do not overwrite")
dir.create(out)
logmsg<-function(...){s<-paste(format(Sys.time(),tz="UTC",usetz=TRUE),...);cat(s,"\n");cat(s,"\n",file=file.path(out,"run.log"),append=TRUE)}
options(error=function(){logmsg("FAILED stage; diagnose before continuing");q(status=1)})
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
old<-file.path(w,"maf_sensitivity_20260905")
iv<-fread(file.path(old,"inputs/instruments_verified.csv"))
fr<-fread(file.path(old,"inputs/eur_freq_selected.tsv"))
af<-fread(file.path(old,"inputs/eqtlgen_af_selected_v1.csv"))
stopifnot(nrow(iv)==6135,uniqueN(iv$gene_symbol)==2544,!anyDuplicated(fr$SNP),!anyDuplicated(af$snp))
ex<-merge(iv,fr[,.(snp=SNP,frq_A1=A1,frq_A2=A2,MAF)],by="snp",all.x=TRUE)
ex[,eaf_ref:=fifelse(effect_allele==frq_A1,MAF,fifelse(effect_allele==frq_A2,1-MAF,NA_real_))]
ex<-merge(ex,af[,.(snp,allele_a,allele_b,af_b)],by="snp",all.x=TRUE)
ex[,pair_ok:=(effect_allele==allele_b&other_allele==allele_a)|(effect_allele==allele_a&other_allele==allele_b)]
ex[,eaf_new:=fifelse(pair_ok,fifelse(effect_allele==allele_b,af_b,1-af_b),NA_real_)]
ex[,eligible_af:=!is.na(eaf_new)&eaf_new>0&eaf_new<1]
ex[,se_ref:=1/sqrt(2*eaf_ref*(1-eaf_ref)*(n_samples+zscore^2))]
ex[,se_new:=fifelse(eligible_af,1/sqrt(2*eaf_new*(1-eaf_new)*(n_samples+zscore^2)),NA_real_)]
ex[,`:=`(beta_ref=zscore*se_ref,beta_new=zscore*se_new,scale_new_ref=se_new/se_ref)]
ex[,delta_z:=beta_new/se_new-beta_ref/se_ref]
stopifnot(max(abs(ex$delta_z),na.rm=TRUE)<1e-8)
fwrite(ex,file.path(out,"exposure_comparison.csv"))
genes<-unique(iv$gene_symbol)
hits<-fread(file.path(w,"expanded_analysis_20260930/results/Stage1_13_discovery_genes.csv"))[outcome=="BBJ_Graves",gene_symbol]
focus<-c(hits,"IGF1R");stopifnot(length(focus)==14)
warns<-character();verify_delta<-0;results<-list();harm<-list();qc<-list();loo<-list()
fit<-function(h,g,oc,sc,removed=NA_character_){
  h<-h[mr_keep==TRUE];n<-nrow(h);stopifnot(!anyDuplicated(h$SNP))
  if(n==0)return(data.table(gene_symbol=g,outcome=oc,scenario=sc,n_iv=0L,method=NA_character_,beta=NA_real_,se=NA_real_,pvalue=NA_real_,OR=NA_real_,CI_lower=NA_real_,CI_upper=NA_real_,Q=NA_real_,Q_df=NA_integer_,Q_pvalue=NA_real_,removed_snp=removed))
  r<-suppressMessages(mr(as.data.frame(h),method_list=if(n==1)"mr_wald_ratio" else "mr_ivw"))
  # Independent closed-form calculation from the retained SNP-level estimates.
  bx<-h$beta.exposure;by<-h$beta.outcome;sy<-h$se.outcome
  b<-sum(bx*by/sy^2)/sum(bx^2/sy^2)
  q<-sum((by-b*bx)^2/sy^2)
  se<-sqrt(if(n==1)1/sum(bx^2/sy^2) else max(1,q/(n-1))/sum(bx^2/sy^2))
  p<-2*pnorm(abs(b/se),lower.tail=FALSE)
  delta<-max(abs(c(r$b-b,r$se-se,r$pval-p)))
  verify_delta<<-max(verify_delta,delta)
  stopifnot(nrow(r)==1,r$nsnp==n,delta<1e-8,is.finite(r$b),is.finite(r$se),r$se>0)
  data.table(gene_symbol=g,outcome=oc,scenario=sc,n_iv=n,method=r$method,beta=r$b,se=r$se,pvalue=r$pval,OR=exp(r$b),CI_lower=exp(r$b-qnorm(.975)*r$se),CI_upper=exp(r$b+qnorm(.975)*r$se),Q=if(n>1)q else NA_real_,Q_df=if(n>1)n-1L else NA_integer_,Q_pvalue=if(n>1)pchisq(q,n-1,lower.tail=FALSE) else NA_real_,removed_snp=removed)
}
harmonize<-function(eg,d,g,oc,kind){
  if(nrow(eg)==0 || nrow(d)==0)return(data.table(SNP=character(),mr_keep=logical()))
  ef<-data.frame(SNP=eg$snp,beta.exposure=eg[[paste0("beta_",kind)]],se.exposure=eg[[paste0("se_",kind)]],effect_allele.exposure=eg$effect_allele,other_allele.exposure=eg$other_allele,eaf.exposure=eg[[paste0("eaf_",kind)]],pval.exposure=eg$pvalue,exposure=g,id.exposure=g)
  of<-data.frame(SNP=d$snp,effect_allele.outcome=toupper(d$ea),other_allele.outcome=toupper(d$oa),beta.outcome=d$beta,se.outcome=d$se,pval.outcome=d$pvalue,eaf.outcome=d$eaf,outcome=oc,id.outcome=oc)
  as.data.table(withCallingHandlers(suppressMessages(harmonise_data(ef,of,action=2)),warning=function(x){warns<<-c(warns,paste(g,oc,conditionMessage(x)));invokeRestart("muffleWarning")}))
}
for(oc in c("FinnGen_GO","E4_GRAVES_STRICT","E4_GRAVES_OPHT_STRICT")){
  d<-fread(if(oc=="FinnGen_GO")file.path(old,"inputs/FinnGen_selected.tsv") else file.path(t,"inputs_v2",paste0(oc,"_selected.tsv")))
  d<-d[snp %in% iv$snp]
  # Do not compare raw positions: expression GRCh37, disease GRCh38.
  for(i in seq_along(genes)){
    g<-genes[i];eg<-ex[gene_symbol==g];dg<-d[snp %in% eg$snp]
    stopifnot(all(dg$chr %in% eg$chr))
    h<-harmonize(eg,dg,g,oc,"ref")
    results[[length(results)+1]]<-fit(h,g,oc,"original_reference")
    qc[[length(qc)+1]]<-data.table(gene_symbol=g,outcome=oc,n_selected=nrow(eg),n_matched=nrow(dg),n_retained=sum(h$mr_keep),n_ambiguous=sum(h$ambiguous),n_incompatible=sum(h$remove))
    if(nrow(h)>0){hh<-copy(h);hh[,scenario:="original_reference"];harm[[length(harm)+1]]<-hh}
    if(oc!="FinnGen_GO"){
      ev<-eg[eligible_af==TRUE];common<-copy(h[mr_keep==TRUE & SNP %in% ev$snp])
      results[[length(results)+1]]<-fit(common,g,oc,"paired_reference")
      nh<-copy(common)
      if(nrow(nh)>0){scales<-ev$scale_new_ref[match(nh$SNP,ev$snp)];nh[,beta.exposure:=beta.exposure*scales];nh[,se.exposure:=se.exposure*scales]}
      results[[length(results)+1]]<-fit(nh,g,oc,"paired_eqtlgen")
      newh<-harmonize(ev,dg[snp %in% ev$snp],g,oc,"new")
      results[[length(results)+1]]<-fit(newh,g,oc,"reharmonized_eqtlgen")
      if(nrow(newh)>0){hh<-copy(newh);hh[,scenario:="reharmonized_eqtlgen"];harm[[length(harm)+1]]<-hh}
      if(g %in% focus && sum(h$mr_keep)>=2){
        for(s in h[mr_keep==TRUE,SNP])loo[[length(loo)+1]]<-fit(h[SNP!=s],g,oc,"leave_one_out",s)
      }
    }
    if(i%%500==0 || i==length(genes))logmsg(oc,"genes",i,"/",length(genes))
  }
  if(oc=="FinnGen_GO"){
    base<-rbindlist(results)[outcome==oc & n_iv>0]
    canonical<-fread(file.path(old,"baseline_v1/primary_MR_reproduced.csv"))[outcome==oc]
    cmp<-merge(base,canonical,by=c("gene_symbol","outcome"),all=TRUE,suffixes=c("_new","_old"))
    cmp[,pass:=n_iv_new==n_iv_old & abs(beta_new-beta_old)<1e-8 & abs(se_new-se_old)<1e-8 & abs(pvalue_new-pvalue_old)<1e-8]
    fwrite(cmp,file.path(out,"identity_control_original_FinnGen.csv"))
    stopifnot(nrow(cmp)==2480,all(cmp$pass))
    logmsg("IDENTITY CONTROL PASS: all 2480 original FinnGen estimates")
  }
}
res<-rbindlist(results);res[,`:=`(BBJ_discovery=gene_symbol %in% hits,nominal=pvalue<.05,full_screen_bonferroni=pvalue<.05/2544)]
bb<-fread(file.path(w,"expanded_analysis_20260930/results/Stage1_13_discovery_genes.csv"))[outcome=="BBJ_Graves",.(gene_symbol,BBJ_beta=beta)]
res<-merge(res,bb,by="gene_symbol",all.x=TRUE)
res[,focused_GD_corroboration:=outcome=="E4_GRAVES_STRICT" & BBJ_discovery & sign(beta)==sign(BBJ_beta) & pvalue<.05/13]
fwrite(res,file.path(out,"MR_all_genes_all_scenarios.csv"))
fwrite(res[gene_symbol %in% focus],file.path(out,"MR_14genes_all_scenarios.csv"))
fwrite(rbindlist(harm,fill=TRUE),file.path(out,"harmonized_snps.csv"))
fwrite(rbindlist(qc),file.path(out,"harmonization_counts.csv"))
fwrite(rbindlist(loo),file.path(out,"leave_one_out_14genes.csv"))
paired<-merge(res[scenario=="paired_reference"],res[scenario=="paired_eqtlgen"],by=c("gene_symbol","outcome"),suffixes=c("_ref","_new"))
single<-paired[n_iv_ref==1]
stopifnot(all(paired$n_iv_ref==paired$n_iv_new),all(sign(single$beta_ref)==sign(single$beta_new)),max(abs(single$pvalue_ref-single$pvalue_new))<1e-8)
writeLines(unique(warns),file.path(out,"warnings.txt"))
fwrite(data.table(test=c("reference_identity_rows","independent_MR_max_delta","max_delta_exposure_Z","same_SNP_single_IV_max_delta_P","new_primary_gene_endpoint_rows"),value=c(2480,verify_delta,max(abs(ex$delta_z),na.rm=TRUE),max(abs(single$pvalue_ref-single$pvalue_new)),nrow(res[outcome!="FinnGen_GO" & scenario=="original_reference"]))),file.path(out,"verification.csv"))
logmsg("MR PASS; independent maximum difference",verify_delta)
