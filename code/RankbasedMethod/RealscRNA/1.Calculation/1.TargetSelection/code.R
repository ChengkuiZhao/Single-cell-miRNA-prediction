load(file = 'TNBC_miRNA.Rdata')
load(file = 'TNBC_RNA.Rdata')
hsa=read.csv('hsa_RNA_miRNA_score')
hsa=as.matrix(hsa[,-1])
miRNAname_Targetscan=unique(hsa[,2])
miRNAname_select=intersect(rownames(TNBC_miRNA),miRNAname_Targetscan)
TNBC_miRNA_select=TNBC_miRNA[miRNAname_select,]
TNBC_RNA_rsum=apply(TNBC_RNA,1,sum)
C=list()
P=list()
X=list()
Y=list()
for (i in 1:length(rownames(TNBC_miRNA_select))){
  hsa_RNA=hsa[which(hsa[,2]==rownames(TNBC_miRNA_select)[i]),1]
  RNA_intersect=intersect(rownames(TNBC_RNA),hsa_RNA)
  if(length(RNA_intersect)>0){
    c=rep(0,length(RNA_intersect))
    p=rep(0,length(RNA_intersect))
    for (j in 1:length(RNA_intersect)){
      c[j]=cor.test(TNBC_miRNA_select[i,],TNBC_RNA[RNA_intersect[j],],'l','pearson')$estimate
      p[j]=cor.test(TNBC_miRNA_select[i,],TNBC_RNA[RNA_intersect[j],],'l','pearson')$p.value
    }
    RNA_intersect_expression=TNBC_RNA_rsum[RNA_intersect]
    C[i]=list(c)
    P[i]=list(p)
    if (length(which(c<(-0.2)))>0){
      X[i]=list(RNA_intersect[which(c<(-0.2))])
    }
    Y[i]=list(RNA_intersect_expression)
  }
  print(i)
}
##Plot
TNBC_miRNA_select_rsum=apply(TNBC_miRNA_select,1,sum)
o=order(TNBC_miRNA_select_rsum,decreasing = T)
index=rep(0,348)
for (i in o[1:348]){
  RNA_intersect=unlist(X[i])
  if(length(RNA_intersect)>0){
    c=rep(0,length(RNA_intersect))
    for (j in 1:length(RNA_intersect)){
      c[j]=cor.test(TNBC_miRNA_select[i,],TNBC_RNA[RNA_intersect[j],],'l','pearson')$estimate
    }
    RNA_intersect_expression=TNBC_RNA_rsum[RNA_intersect]
  }
  index[which(o==i)]=mean(c)
}
library(ggplot2)
dat=data.frame(xvar=c(1:348),yvar=index)
ggplot(dat, aes(x=xvar, y=yvar)) +
  geom_point(shape=1) +
  geom_smooth(method=lm)+xlab('miRNA ordered by expression')+ylab('Mean of the miRNA-RNA Pearson correlations')+ggtitle('Overall Correlation')+ylim(-0.5,0)
save(X,file='Targetgene_selected.Rdata')
save(TNBC_miRNA_select,file = 'TNBC_miRNA_select.Rdata')
save(TNBC_miRNAname_select,file = 'TNBC_miRNAname_select.Rdata')


