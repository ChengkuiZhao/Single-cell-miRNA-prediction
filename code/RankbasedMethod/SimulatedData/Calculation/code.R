/N/u/cz18/Carbonate/software/R-3.5.3/bin/R
####Extract 3000 genes per cell and get score
load(file='data_var0.2.RData')
load(file='genes_selected(-0.2).Rdata')
load(file='miRNAname_select.Rdata')
data=data_var0.2
rm(data_var0.2)
data_names=rownames(data)
miRNA_name=miRNAname_select
s=c(3000:1)
scores_m=matrix(0,nrow=344,ncol=ncol(data))
n=0
for (i in 1:344){
miRNA_targetnames=unlist(X[i])
    for (j in 1:ncol(data)){
               o=order(data[,j],decreasing=T)
               g_n=data_names[o[1:3000]]
               names(s)=g_n
               x=intersect(g_n,miRNA_targetnames)
               score=sum(s[x])
               scores_m[i,j]=score
                                    }
               n=n+1
               print(n)

                    }
#######filter and get miRNA score
rownames(scores_m)=miRNA_name[1:344]
colnames(scores_m)=c(1:10000)
save(scores_m,file='scores_m.Rdata')
reverse=function(x)((max(x)-x)/(max(x)-min(x)))
miRNA_score=apply(scores_m,1,reverse)
miRNA_score=t(miRNA_score)
miRNA_score[is.na(miRNA_score)]=0
m=apply(scores_m,1,mean)
a=which(m>0)
scores_m_select=scores_m[a,]
miRNA_score_select=miRNA_score[rownames(scores_m_select),]
miRNA_score_select=round(miRNA_score_select,3)
save(miRNA_score_select,file='miRNA.Rdata')


