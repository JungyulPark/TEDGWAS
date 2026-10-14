args<-commandArgs(trailingOnly=TRUE)
src<-normalizePath(args[1],winslash='/',mustWork=TRUE)
out<-args[2];dir.create(out,recursive=TRUE,showWarnings=FALSE)
d<-read.csv(file.path(src,'Tables/Tissue_MR_all_estimable.csv'),stringsAsFactors=FALSE)
av<-read.csv(file.path(src,'Tables/Tissue_availability_10_combinations.csv'),stringsAsFactors=FALSE)
ord<-c('TSHR|Adipose_Subcutaneous','TSHR|Adipose_Visceral_Omentum','IGF1R|Thyroid','IGF1R|Cells_Cultured_fibroblasts')
p<-d[d$outcome=='E4_GRAVES_STRICT',];p<-p[match(ord,paste(p$gene,p$tissue,sep='|')),]
fmt<-function(x){if(x<.001)formatC(x,format='e',digits=2) else if(x<.1)sprintf('%.4f',x) else sprintf('%.3f',x)}
draw<-function(){
 layout(matrix(c(1,2),1,2),widths=c(1,2.05))
 par(mar=c(2.8,.5,3.8,.4),oma=c(3.8,.3,1.5,.3),family='sans',fg='black',col.axis='black',col.lab='black',col.main='black')
 plot.new();plot.window(xlim=c(0,1),ylim=c(0,7),xaxs='i',yaxs='i')
 text(.01,6.85,'A  Instrument availability',adj=0,font=2,cex=1.15)
 text(.01,6.20,'GTEx v8 tissue',adj=0,font=2,cex=.88)
 text(.65,6.20,'TSHR',font=4,cex=.95);text(.88,6.20,'IGF1R',font=4,cex=.95)
 tissue<-c('Thyroid','Cells_Cultured_fibroblasts','Adipose_Subcutaneous','Adipose_Visceral_Omentum','Whole_Blood')
 labs<-c('Thyroid','Cultured fibroblasts','Subcutaneous adipose','Visceral adipose','Whole blood')
 for(i in 1:5){
  yy<-5.4-(i-1)*.9;text(.01,yy,labs[i],adj=0,cex=.88)
  for(j in 1:2){
   gene<-c('TSHR','IGF1R')[j];a<-av[av$gene==gene & av$tissue==tissue[i],];xx<-c(.65,.88)[j]
   if(a$n_selected_instruments==1){rect(xx-.062,yy-.19,xx+.062,yy+.19,col='#C7DDE5',border='#4E7181');text(xx,yy,'1',font=2,cex=.95)}
   else text(xx,yy,if(a$tested)'x' else '-',cex=1.15)
  }
 }
 text(.01,1.07,'1  One qualifying instrument',adj=0,cex=.79)
 text(.01,.68,'x  No eQTL meeting P < 5e-8',adj=0,cex=.79)
 text(.01,.29,'-  Gene absent from the tested-gene file',adj=0,cex=.79)
 par(mar=c(2.8,.3,3.8,.4))
 plot.new();plot.window(xlim=c(0,1),ylim=c(0,7),xaxs='i',yaxs='i')
 text(.015,6.85,'B  Primary comparison within the tissue extension',adj=0,font=2,cex=1.15)
 text(.015,6.40,'Disease dataset: FinnGen Graves disease',adj=0,cex=.89)
 text(.015,5.85,'Expression instrument',adj=0,font=2,cex=.88)
 text(.64,5.85,'OR (95% CI)',font=2,cex=.88);text(.81,5.85,'P value',font=2,cex=.88);text(.949,5.85,'Adjusted P',font=2,cex=.88)
 sx<-function(v).345+(log(v)-log(.25))/(log(2)-log(.25))*.21
 yy<-c(5.04,4.02,2.72,1.70)
 segments(sx(1),1.15,sx(1),5.46,lty=2,col='#888888',lwd=1)
 lab<-c('TSHR\nSubcutaneous adipose','TSHR\nVisceral adipose','IGF1R\nThyroid','IGF1R\nCultured fibroblasts')
 for(i in 1:4){
  color<-if(p$gene[i]=='TSHR')'#236787' else '#AC5034';hit<-p$P_adjusted[i]<.05
  text(.015,yy[i],lab[i],adj=0,cex=.88)
  segments(sx(p$CI_lower[i]),yy[i],sx(p$CI_upper[i]),yy[i],col=color,lwd=2)
  points(sx(p$OR[i]),yy[i],pch=if(hit)15 else 0,col=color,cex=1.05,lwd=1.5)
  text(.64,yy[i],sprintf('%.2f (%.2f-%.2f)',p$OR[i],p$CI_lower[i],p$CI_upper[i]),cex=.84)
  text(.81,yy[i],fmt(p$P[i]),font=if(hit)2 else 1,cex=.84)
  text(.949,yy[i],paste0(fmt(p$P_adjusted[i]),if(hit)'*' else ''),font=if(hit)2 else 1,cex=.84)
 }
 axis(1,at=sx(c(.25,.5,1,2)),labels=c('0.25','0.5','1','2'),pos=.99,tck=-.015,cex.axis=.8)
 text(.45,.34,'Odds ratio',cex=.84)
 mtext('TSHR and IGF1R tissue extension',outer=TRUE,side=3,line=.2,font=2,cex=1.32)
 mtext('* Corrected support: raw P < 0.005; adjusted P = min(1, 10 x P). Filled squares mark corrected support.',outer=TRUE,side=1,line=.5,adj=0,cex=.83)
 mtext('Ten primary combinations were planned; six were unestimable. All 40 secondary comparisons are reported separately.',outer=TRUE,side=1,line=1.6,adj=0,cex=.83)
 mtext('GTEx blood differs from eQTLGen. ORs are per normalized expression unit. Non-orbital tissues; no tissue-specific colocalization.',outer=TRUE,side=1,line=2.7,adj=0,cex=.83)
}
pdf(file.path(out,'Tissue_primary_summary.pdf'),width=15,height=7.5,useDingbats=FALSE);draw();dev.off()
png(file.path(out,'Tissue_primary_summary.png'),width=3600,height=1800,res=240,type='cairo');draw();dev.off()
