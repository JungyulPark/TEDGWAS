suppressPackageStartupMessages({library(data.table);library(coloc)})
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1)
t<-normalizePath(args[1],winslash="/",mustWork=TRUE);o<-file.path(t,"coloc_v1")
if(file.exists(file.path(o,"R_coloc_verification.csv")))stop("Validation already completed")
writeLines(capture.output(sessionInfo()),file.path(o,"R_sessionInfo.txt"))
expected<-fread(file.path(o,"Coloc_all_settings_with_identity.csv"));ans<-list();warnings<-character()
for(f in list.files(file.path(o,"r_validation"),pattern="csv$",full.names=TRUE)){
 d<-fread(f);parts<-strsplit(sub(".csv$","",basename(f)),"__",fixed=TRUE)[[1]]
 for(p in c(1e-5,5e-6,1e-6)){
  r<-withCallingHandlers(suppressMessages(coloc.abf(dataset1=list(beta=d$beta1,varbeta=d$varbeta1,snp=d$snp,type="cc"),dataset2=list(beta=d$beta2,varbeta=d$varbeta2,snp=d$snp,type="quant",sdY=1),p1=1e-4,p2=1e-4,p12=p)),warning=function(x){warnings<<-c(warnings,paste(basename(f),conditionMessage(x)));invokeRestart("muffleWarning")})
  s<-r$summary;row<-data.table(gene=parts[1],outcome=parts[2],scenario=parts[3],p12=p,n_overlap=unname(s['nsnps']))
  for(k in 0:4)row[[paste0('PP.H',k)]]<-unname(s[paste0('PP.H',k,'.abf')])
  ans[[length(ans)+1]]<-row
 }
}
got<-rbindlist(ans);fwrite(got,file.path(o,"R_coloc_posteriors.csv"))
m<-merge(got,expected,by=c("gene","outcome","scenario","p12"),suffixes=c("_R","_Python"))
delta<-max(sapply(0:4,function(k)max(abs(m[[paste0("PP.H",k,"_R")]]-m[[paste0("PP.H",k,"_Python")]]))))
stopifnot(nrow(m)==504,all(m$n_overlap_R==m$n_overlap_Python),delta<1e-8)
fwrite(data.table(test=c("posterior_vectors_validated","maximum_absolute_posterior_difference","PASS"),value=c(nrow(m),delta,1)),file.path(o,"R_coloc_verification.csv"))
writeLines(unique(warnings),file.path(o,"R_warnings.txt"));cat("PASS R coloc:",nrow(m),"vectors; max delta",delta,"\n")
