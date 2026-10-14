.libPaths(c('${LOCAL_USER_HOME}/AppData/Local/R/win-library/4.3',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(TwoSampleMR)})
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1)
t<-normalizePath(args[1],winslash='/',mustWork=TRUE);w<-dirname(t);o<-file.path(t,'instruments_v1')
stopifnot(file.exists(file.path(t,'candidates_v1/R_verification.csv')),!dir.exists(o));dir.create(o)
options(error=function(){cat('FAILED: diagnose before continuing\n',file=file.path(o,'run.log'),append=TRUE);q(status=1)})
writeLines(capture.output(sessionInfo()),file.path(o,'R_sessionInfo.txt'))
h<-fread(file.path(w,'reanalysis_20261002/mr_v1/harmonized_snps.csv'))
h<-h[exposure %in% c('C4B','HLA-G') & outcome %in% c('E4_GRAVES_STRICT','E4_GRAVES_OPHT_STRICT') & mr_keep==TRUE]
fwrite(h,file.path(o,'retained_SNP_inputs.csv'))
old<-fread(file.path(w,'reanalysis_20261002/mr_v1/MR_all_genes_all_scenarios.csv'))
old<-old[gene_symbol %in% c('C4B','HLA-G') & outcome %in% c('E4_GRAVES_STRICT','E4_GRAVES_OPHT_STRICT')]
ex<-fread(file.path(w,'reanalysis_20261002/mr_v1/exposure_comparison.csv'))[gene_symbol %in% c('C4B','HLA-G')]
fwrite(ex,file.path(o,'exposure_frequency_scaling.csv'))
full<-list();loo<-list();snps<-list();maxdelta<-0
fit<-function(d,removed='none'){
 n<-nrow(d);stopifnot(n>=1,!anyDuplicated(d$SNP))
 r<-suppressMessages(mr(as.data.frame(d),method_list=if(n==1)'mr_wald_ratio' else 'mr_ivw'))
 b<-r$b;s<-r$se;bx<-d$beta.exposure;by<-d$beta.outcome;sy<-d$se.outcome
 q<-sum((by-b*bx)^2/sy^2)
 b2<-sum(bx*by/sy^2)/sum(bx^2/sy^2);s2<-sqrt(if(n==1)1/sum(bx^2/sy^2) else max(1,q/(n-1))/sum(bx^2/sy^2))
 delta<-max(abs(c(b-b2,s-s2,r$pval-2*pnorm(abs(b2/s2),lower.tail=FALSE))));maxdelta<<-max(maxdelta,delta);stopifnot(delta<1e-8)
 data.table(gene=d$exposure[1],outcome=d$outcome[1],scenario=d$scenario[1],removed_snp=removed,n_iv=n,beta=b,se=s,pvalue=r$pval,OR=exp(b),CI_lower=exp(b-qnorm(.975)*s),CI_upper=exp(b+qnorm(.975)*s),Q=if(n>1)q else NA_real_,Q_df=if(n>1)n-1 else NA_integer_,Q_pvalue=if(n>1)pchisq(q,n-1,lower.tail=FALSE) else NA_real_)
}
for(g in c('C4B','HLA-G'))for(oc in c('E4_GRAVES_STRICT','E4_GRAVES_OPHT_STRICT')){
 a<-h[exposure==g & outcome==oc & scenario=='original_reference'];b<-h[exposure==g & outcome==oc & scenario=='reharmonized_eqtlgen']
 stopifnot(setequal(a$SNP,b$SNP))
 pair<-merge(a,b,by='SNP',suffixes=c('_ref','_new'))
 stopifnot(max(abs(pair$beta.outcome_ref-pair$beta.outcome_new))<1e-12,max(abs(pair$se.outcome_ref-pair$se.outcome_new))<1e-12)
 stopifnot(max(abs(pair$beta.exposure_ref/pair$se.exposure_ref-pair$beta.exposure_new/pair$se.exposure_new))<1e-8)
 # Check the same-SNP frequency-only estimator equals the separately re-harmonized result here.
 paired<-old[gene_symbol==g & outcome==oc & scenario=='paired_eqtlgen'];rh<-old[gene_symbol==g & outcome==oc & scenario=='reharmonized_eqtlgen']
 stopifnot(nrow(paired)==1,nrow(rh)==1,paired$n_iv==rh$n_iv,max(abs(c(paired$beta-rh$beta,paired$se-rh$se,paired$pvalue-rh$pvalue)))<1e-8)
 for(sc in c('original_reference','reharmonized_eqtlgen')){
  d<-h[exposure==g & outcome==oc & scenario==sc];r<-fit(d);full[[length(full)+1]]<-r
  prior<-old[gene_symbol==g & outcome==oc & scenario==sc]
  stopifnot(max(abs(c(r$beta-prior$beta,r$se-prior$se,r$pvalue-prior$pvalue)))<1e-8)
  bx<-d$beta.exposure;by<-d$beta.outcome;sy<-d$se.outcome;bw<-by/bx;sw<-sy/abs(bx);ww<-bx^2/sy^2
  ss<-data.table(gene=g,outcome=oc,scenario=sc,SNP=d$SNP,beta_exposure=bx,se_exposure=d$se.exposure,eaf_exposure=d$eaf.exposure,beta_outcome=by,se_outcome=sy,beta_wald=bw,se_wald=sw,OR_wald=exp(bw),CI_lower=exp(bw-qnorm(.975)*sw),CI_upper=exp(bw+qnorm(.975)*sw),P_wald=2*pnorm(abs(bw/sw),lower.tail=FALSE),weight_fraction=ww/sum(ww),Q_contribution=(by-r$beta*bx)^2/sy^2,F_exposure=(bx/d$se.exposure)^2)
  snps[[length(snps)+1]]<-ss
  for(s in d$SNP)loo[[length(loo)+1]]<-fit(d[SNP!=s],s)
 }
}
fwrite(rbindlist(full),file.path(o,'full_MR.csv'));fwrite(rbindlist(loo),file.path(o,'leave_one_out.csv'));fwrite(rbindlist(snps),file.path(o,'SNP_estimates_weights.csv'))
fwrite(data.table(test=c('full_estimates','leave_one_out_estimates','SNP_estimates','R_vs_formula_max_delta','same_SNP_sets_and_outcome_estimates','paired_frequency_equals_reharmonized','PASS'),value=c(length(full),length(loo),nrow(rbindlist(snps)),maxdelta,1,1,1)),file.path(o,'R_verification.csv'))
cat('PASS: full',length(full),'; leave-one-out',length(loo),'; maximum formula delta',maxdelta,'\n')
