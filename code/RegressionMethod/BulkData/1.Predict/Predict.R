load(file='Correlation_abs.Rdata')
load(file = 'label_1500.Rdata')
load(file='miRNA_test.Rdata')
load(file = 'miRNA_train_mean.Rdata')
load(file='miRNA_train_sd.Rdata')
load(file='miRNA_train_sd.Rdata')
load(file='RNA_test.Rdata')
load(file = 'RNA_train_mean.Rdata')
load(file='RNA_train_sd.Rdata')
load(file='TrainedWeights.Rdata')
RNA=RNA_test
miRNA=miRNA_test
remove(RNA_test)
remove(miRNA_test)
RNA=as.matrix(RNA)
miRNA=as.matrix(miRNA)

RNA_standardize=matrix(0,nrow(RNA),ncol(RNA))
for (i in c(1:25696)){
  RNA_standardize[i,]=((RNA[i,]-RNA_train_mean[i])/RNA_train_sd[i])
}
RNA_standardize[is.na(RNA_standardize)]=0
rownames(RNA_standardize)=rownames(RNA)
colnames(RNA_standardize)=colnames(RNA)
##Kmeans cluster k=1500
RNA_st_s1=matrix(0,1500,ncol(RNA))
for (j in c(1:1500)) {
  if(length(which(label_1500==j))==1){
    RNA_st_s1[j,]=RNA_standardize[which(label_1500==j),]
  }
  if(length(which(label_1500==j))>1){
    RNA_st_s1[j,]=apply(RNA_standardize[which(label_1500==j),],2,mean)
  }
}
##Use Correlation_abs,select top30 groups
TopN=function(x,n){o=order(x,decreasing = T) 
                   top=o[1:n]}
topn=apply(Correlation_abs,2,TopN,n=30)
##Predict with selected groups
pred=function(x,w){sum(c(1,x)*w)}
Y_p=matrix(0,nrow(miRNA),ncol(miRNA))
for (k in c(1:2238)){
  x=t(RNA_st_s1[topn[,k],])
  y_p=apply(x,1,pred,w=W[k,])
  Y_p[k,]=y_p
}
##Recover prediction
miRNA_pred=matrix(0,nrow(miRNA),ncol(miRNA))
for(q in c(1:nrow(Y_p))){
  miRNA_pred[q,]=Y_p[q,]*miRNA_train_sd[q]+miRNA_train_mean[q]
}
##Record the result
Cor_row=rep(0,nrow(miRNA))
Cor_col=rep(0,ncol(miRNA))
Cor_train_mean=rep(0,ncol(miRNA))
for(r1 in c(1:nrow(miRNA))){
  Cor_row[r1]=cor.test(miRNA[r1,],miRNA_pred[r1,])$estimate
}
for (r2 in c(1:ncol(miRNA))) {
  Cor_col[r2]=cor.test(miRNA[,r2],miRNA_pred[,r2])$estimate
}
for (r3 in c(1:ncol(miRNA))) {
  Cor_train_mean[r3]=cor.test(miRNA[,r3],miRNA_train_mean)$estimate
}
result1=Cor_row
result2=Cor_col
result3=Cor_train_mean

