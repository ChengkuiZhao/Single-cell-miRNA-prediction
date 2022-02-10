load(file = 'TNBC_miRNA.Rdata')
miRNA=read.csv('miRNA')
index=rep(0,832)
for (i in c(1:832)){
  index[i]=unlist(strsplit(colnames(miRNA)[i],split = '[.]'))[4]
}
Normal_miRNA=miRNA[,which(index=='11')]

library(ggplot2)
Phenotype=read.csv('Phenotype.csv')
Stage=Phenotype[,2]
names(Stage)=Phenotype[,1]
###Check miRNA between Normal&TNBC 
status=c(rep('Normal',76),rep('TNBC',67))
#mir='hsa-miR-98-3p'
#mir='hsa-let-7d-5p'
#mir='hsa-miR-18a-3p'
mir='hsa-let-7f-1-3p'
miR=c(Normal_miRNA[mir,],TNBC_miRNA[mir,])
miRNA=data.frame(status,miR)
t.test(Normal_miRNA[mir,],TNBC_miRNA[mir,],alternative = 'less')
p=ggplot(miRNA,aes(x=status,y=miR,fill=status))+geom_boxplot(outlier.shape=NA)
p=p+ggtitle(paste(mir,'in Normal and TNBC\nt-test P value = 5.978e-10',round(median(Normal_miRNA[mir,]),3),round(median(TNBC_miRNA[mir,]),3)))
p=p+theme(plot.title = element_text(hjust = 0.5))+scale_x_discrete(limits=c('Normal','TNBC'))
p















