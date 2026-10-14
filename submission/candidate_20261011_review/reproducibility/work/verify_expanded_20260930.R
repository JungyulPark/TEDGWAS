args <- commandArgs(trailingOnly=TRUE)
root <- args[1]
ladd <- function(a,b) {m<-pmax(a,b); ifelse(is.infinite(m),m,m+log(exp(a-m)+exp(b-m)))}
lsum <- function(x) {m<-max(x);m+log(sum(exp(x-m)))}
rows<-list()
for (f in list.files(file.path(root,'inputs/r_validation'),pattern='csv$',full.names=TRUE)) {
 d<-read.csv(f);keys<-strsplit(sub('\\.csv$','',basename(f)),'__',fixed=TRUE)[[1]]
 z1<-d$beta1/sqrt(d$varbeta1);z2<-d$beta2/sqrt(d$varbeta2)
 r1<-.04/(.04+d$varbeta1);r2<-.0225/(.0225+d$varbeta2)
 a<-.5*(log1p(-r1)+r1*z1*z1);b<-.5*(log1p(-r2)+r2*z2*z2)
 # Sum all distinct-variant pairs with prefix/suffix sums, avoiding subtraction
 # of almost equal total and shared-variant Bayes factors.
 n<-length(a);pre<-rep(-Inf,n+1);suf<-rep(-Inf,n+1)
 for(i in seq_len(n))pre[i+1]<-ladd(pre[i],b[i])
 for(i in n:1)suf[i]<-ladd(suf[i+1],b[i])
 distinct<-lsum(a+ladd(pre[seq_len(n)],suf[seq_len(n)+1]))
 for(p12 in c(1e-6,5e-6,1e-5)) {
  h<-c(0,log(1e-4)+lsum(a),log(1e-4)+lsum(b),2*log(1e-4)+distinct,log(p12)+lsum(a+b))
  pp<-exp(h-lsum(h))
  rows[[length(rows)+1]]<-data.frame(gene=keys[1],outcome=keys[2],scenario=keys[3],p12=p12,H0=pp[1],H1=pp[2],H2=pp[3],H3=pp[4],H4=pp[5])
 }
}
out<-do.call(rbind,rows);write.csv(out,file.path(root,'results/independent_R_posteriors.csv'),row.names=FALSE)
canonical<-read.csv(file.path(root,'results/Stage2_all_504_colocalization_settings.csv'),check.names=FALSE)
cmp<-merge(out,canonical,by=c('gene','outcome','scenario','p12'))
delta<-max(sapply(0:4,function(i)max(abs(cmp[[paste0('H',i)]]-cmp[[paste0('PP.H',i)]]))))
stopifnot(nrow(cmp)==504,delta<1e-9)
cat('Independent base-R pair-enumeration equations: PASS; 504 rows; maximum absolute difference:',format(delta,digits=16),'\n')
writeLines(c(paste('R',getRversion()),'Implementation: independent base-R ABF and prefix/suffix sum of distinct-variant pairs; not a coloc package run.',paste('Rows',nrow(cmp)),paste('Max absolute posterior difference',format(delta,digits=16)),'PASS'),file.path(root,'results/independent_R_verification.txt'))
