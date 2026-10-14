# Independent base-R validation of the post hoc SNP-exclusion results.
# Arguments: local harmonized CSV, output CSV. No source variant rows exported.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2, !file.exists(args[2]))
h <- read.csv(args[1], check.names=FALSE, stringsAsFactors=FALSE)
h <- h[h$mr_keep & h$exposure %in% c('TSHR','IGF1R','CTLA4'), ]
estimate <- function(d) {
  n <- nrow(d)
  x <- d$beta.exposure; y <- d$beta.outcome; sy <- d$se.outcome
  stopifnot(n>=1, all(is.finite(c(x,y,sy))), all(sy>0), all(x!=0))
  if (n==1) {
    b <- y/x; se <- sy/abs(x); Q <- NA_real_
  } else {
    m <- lm(y ~ x - 1, weights=1/sy^2)
    s <- summary(m)
    b <- unname(coef(m)[1])
    se <- s$coef[1,2]/min(1,s$sigma)
    Q <- sum((residuals(m)/sy)^2)
  }
  data.frame(n_iv=n, beta=b, se=se, pvalue=2*pnorm(abs(b/se), lower.tail=FALSE),
             OR=exp(b), CI_lower=exp(b-qnorm(.975)*se), CI_upper=exp(b+qnorm(.975)*se), Q=Q)
}
out <- list()
for (g in c('TSHR','IGF1R','CTLA4')) {
  for (o in c('BBJ_Graves','UKB_hyperthyroid','FinnGen_GO')) {
    d <- h[h$exposure==g & h$outcome==o, ]
    stopifnot(nrow(d)>0, !anyDuplicated(d$SNP))
    excludes <- c('None (all instruments)',if (nrow(d)>1) d$SNP else character())
    for (snp in excludes) {
      selected <- if (snp=='None (all instruments)') d else d[d$SNP!=snp, ]
      out[[length(out)+1]] <- cbind(data.frame(gene=g,outcome=o,excluded_SNP=snp),estimate(selected))
    }
  }
}
out <- do.call(rbind,out)
stopifnot(nrow(out)==24, sum(out$excluded_SNP!='None (all instruments)')==15)
write.csv(out,args[2],row.names=FALSE,na='NA')
cat('Validated by independent base-R regression: 9 full-set and 15 omission estimates.\n')
cat(R.version.string,'\n')
