args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1)
o<-normalizePath(args[1],winslash='/',mustWork=TRUE)
d<-read.csv(file.path(o,'Tables/Harmonized_instruments_all_outcomes.csv'),stringsAsFactors=FALSE)
r<-read.csv(file.path(o,'Tables/Tissue_MR_all_estimable.csv'),stringsAsFactors=FALSE)
check<-list();computed<-list()
for(i in seq_len(nrow(r))){
 x<-d[d$gene==r$gene[i]&d$tissue==r$tissue[i]&d$outcome==r$outcome[i],]
 stopifnot(nrow(x)==1)
 b<-x$beta_outcome/x$beta_exposure
 se<-x$se_outcome/abs(x$beta_exposure)
 p<-2*pnorm(abs(b/se),lower.tail=FALSE)
 v<-c(beta=b,se=se,OR=exp(b),CI_lower=exp(b-qnorm(.975)*se),CI_upper=exp(b+qnorm(.975)*se),P=p,P_adjusted=min(1,p*r$planned_family_tests[i]))
 want<-as.numeric(r[i,names(v)])
 delta<-max(abs(v-want)/pmax(abs(v),1e-300))
 stopifnot(delta<1e-10)
 check[[i]]<-data.frame(gene=r$gene[i],tissue=r$tissue[i],outcome=r$outcome[i],max_relative_difference=delta,passed=TRUE)
 computed[[i]]<-cbind(r[i,c('gene','tissue','outcome')],as.data.frame(as.list(v)))
}
write.csv(do.call(rbind,check),file.path(o,'R_independent_verification.csv'),row.names=FALSE)
write.csv(do.call(rbind,computed),file.path(o,'Tables/R_recomputed_estimates.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(o,'R_sessionInfo.txt'))
cat('PASS: independent R checks of',nrow(r),'estimates. Max relative difference:',max(do.call(rbind,check)$max_relative_difference),'\n')
