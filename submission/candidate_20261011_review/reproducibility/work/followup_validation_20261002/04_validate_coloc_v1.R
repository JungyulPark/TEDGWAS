.libPaths(c('${LOCAL_USER_HOME}/AppData/Local/R/win-library/4.3',.libPaths()))
suppressPackageStartupMessages({library(data.table);library(coloc)})
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==3)
o<-normalizePath(args[1],winslash='/',mustWork=TRUE);expected<-fread(file.path(o,args[2]));target<-as.integer(args[3])
stopifnot(!file.exists(file.path(o,'R_verification.csv')))
writeLines(capture.output(sessionInfo()),file.path(o,'R_sessionInfo.txt'))
ans<-list();warns<-character();output<-character()
for(f in list.files(file.path(o,'r_validation'),pattern='csv$',full.names=TRUE)){
 d<-fread(f);parts<-strsplit(sub('.csv$','',basename(f)),'__',fixed=TRUE)[[1]]
 for(p in c(1e-5,5e-6,1e-6)){
  output<-c(output,capture.output(r<-withCallingHandlers(suppressMessages(coloc.abf(
   dataset1=list(beta=d$beta1,varbeta=d$varbeta1,snp=d$snp,type='cc'),
   dataset2=list(beta=d$beta2,varbeta=d$varbeta2,snp=d$snp,type='quant',sdY=1),
   p1=1e-4,p2=1e-4,p12=p)),warning=function(x){warns<<-c(warns,paste(basename(f),conditionMessage(x)));invokeRestart('muffleWarning')})))
  s<-r$summary;row<-data.table(gene=parts[1],outcome=parts[2],window=parts[3],scenario=parts[4],p12=p,n_overlap=unname(s['nsnps']))
  for(k in 0:4)row[[paste0('PP.H',k)]]<-unname(s[paste0('PP.H',k,'.abf')])
  ans[[length(ans)+1]]<-row
 }
}
got<-rbindlist(ans);fwrite(got,file.path(o,'R_posteriors.csv'));writeLines(unique(warns),file.path(o,'R_warnings.txt'));writeLines(output,file.path(o,'R_console.txt'))
m<-merge(got,expected,by=c('gene','outcome','window','scenario','p12'),suffixes=c('_R','_Python'))
delta<-max(sapply(0:4,function(k)max(abs(m[[paste0('PP.H',k,'_R')]]-m[[paste0('PP.H',k,'_Python')]]))))
stopifnot(nrow(m)==target,all(m$n_overlap_R==m$n_overlap_Python),delta<1e-8)
fwrite(data.table(test=c('vectors_validated','max_abs_delta','PASS'),value=c(nrow(m),delta,1)),file.path(o,'R_verification.csv'))
cat('PASS:',nrow(m),'vectors; maximum absolute posterior difference',delta,'\n')
