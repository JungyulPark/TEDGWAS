args<-commandArgs(trailingOnly=TRUE);o<-normalizePath(args[1],winslash='/',mustWork=TRUE)
d<-read.csv(file.path(o,'Tables/Tissue_MR_all_estimable.csv'),stringsAsFactors=FALSE)
dir.create(file.path(o,'Figures'),showWarnings=FALSE)
oc<-c('E4_GRAVES_STRICT','BBJ_Graves','UKB_hyperthyroid','FinnGen_GO','E4_GRAVES_OPHT_STRICT')
labs<-c('FinnGen GD','BBJ GD','UKB hyperthyroidism','FinnGen ophthalmopathy','FinnGen strict ophthalmopathy')
fmt<-function(p){if(p<.001)formatC(p,format='e',digits=2) else sprintf('%.4f',p)}
draw<-function(){
 par(mfrow=c(1,2),mar=c(4.8,0.5,3.5,0.5),oma=c(4.2,0.5,1.4,0.5),family='sans',fg='black',col.axis='black',col.lab='black',col.main='black')
 for(g in c('TSHR','IGF1R')){
  tissues<-if(g=='TSHR')c('Adipose_Subcutaneous','Adipose_Visceral_Omentum') else c('Thyroid','Cells_Cultured_fibroblasts')
  tissue_labels<-if(g=='TSHR')c('Subcutaneous adipose','Visceral adipose') else c('Thyroid','Cultured fibroblasts')
  color<-if(g=='TSHR')'#236787' else '#AC5034'
  plot.new();plot.window(xlim=c(0,1),ylim=c(-.3,14),xaxs='i',yaxs='i')
  text(.02,13.7,if(g=='TSHR')'A  TSHR' else 'B  IGF1R',adj=0,font=2,cex=1.25)
  text(.02,12.8,'Disease dataset',adj=0,font=2,cex=.83);text(.73,12.8,'OR (95% CI)',font=2,cex=.83);text(.945,12.8,'P value',font=2,cex=.83)
  sx<-function(or).395+(log(or)-log(.2))/(log(3)-log(.2))*.22
  segments(sx(1),.55,sx(1),12.1,lty=2,col='#888888',lwd=1)
  for(j in 1:2){
   s<-d[d$gene==g & d$tissue==tissues[j],];s<-s[match(oc,s$outcome),]
   yy<-if(j==1)c(10.8,9.85,8.9,7.95,7) else c(4.75,3.8,2.85,1.9,.95)
   text(.02,max(yy)+.77,tissue_labels[j],adj=0,font=2,cex=.98)
   for(i in 1:5){
    text(.02,yy[i],labs[i],adj=0,cex=.76,font=if(i==1)2 else 1)
    segments(sx(s$CI_lower[i]),yy[i],sx(s$CI_upper[i]),yy[i],col=color,lwd=1.8)
    points(sx(s$OR[i]),yy[i],pch=if(s$statistical_category[i]=='family_corrected')15 else 0,col=color,cex=.85,lwd=1.4)
    text(.73,yy[i],sprintf('%.2f (%.2f-%.2f)',s$OR[i],s$CI_lower[i],s$CI_upper[i]),cex=.78)
    text(.945,yy[i],paste0(fmt(s$P[i]),if(s$statistical_category[i]=='family_corrected')'*' else ''),cex=.76,font=if(s$statistical_category[i]=='family_corrected')2 else 1)
   }
  }
  axis(1,at=sx(c(.25,.5,1,2)),labels=c('0.25','0.5','1','2'),cex.axis=.82,tck=-.013)
  mtext('Odds ratio',side=1,line=2.1,at=.50,cex=.82)
 }
 mtext('Tissue-specific expression associations across five disease datasets',outer=TRUE,side=3,line=0,font=2,cex=1.2)
 mtext('Each estimate uses one cis-eQTL. ORs are per unit of normalized, inverse-normal-transformed expression.',outer=TRUE,side=1,line=1.0,cex=.8,adj=0)
 mtext('* Family-corrected: FinnGen GD P < 0.005; other datasets P < 0.00125. Unmarked estimates do not pass these thresholds.',outer=TRUE,side=1,line=2.05,cex=.8,adj=0)
 mtext('FinnGen endpoints overlap. Adipose and cultured fibroblasts are not orbital tissue. Thirty unestimable planned comparisons are tabulated separately.',outer=TRUE,side=1,line=3.1,cex=.8,adj=0)
}
pdf(file.path(o,'Figures/Tissue_extension_forest.pdf'),width=16,height=8.8,useDingbats=FALSE);draw();dev.off()
png(file.path(o,'Figures/Tissue_extension_forest.png'),width=3200,height=1760,res=200,type='cairo');draw();dev.off()
