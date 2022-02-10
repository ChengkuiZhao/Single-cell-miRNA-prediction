load(file='data_final.RData')
load(file='RNA_train_mean.RData')
load(file='RNA_train_sd.RData')
load(file='miRNA_train_mean.RData')
load(file='miRNA_train_sd.RData')
load(file='Correlation_abs.RData')
load(file='label_2000.RData')
load(file='TrainedWeights.RData')
miRNA_10=read.csv('miRNA_10')

RNA_standardize=matrix(0,nrow(RNA),ncol(RNA))
for (i in c(1:20530)){
  RNA_standardize[i,]=(RNA[i,]-RNA_train_mean[i])/RNA_train_sd[i]
}
rownames(RNA_standardize)=rownames(RNA)
colnames(RNA_standardize)=colnames(RNA)
##Kmeans cluster k=2000
RNA_st_s1=matrix(0,2000,ncol(RNA))
for (j in c(1:2000)) {
  if(length(which(label_2000==j))==1){
    RNA_st_s1[j,]=RNA_standardize[which(label_2000==j),]
  }
  if(length(which(label_2000==j))>1){
    RNA_st_s1[j,]=apply(RNA_standardize[which(label_2000==j),],2,mean)
  }
}
##Use Correlation_abs,select top15 groups
TopN=function(x,n){o=order(x,decreasing = T) 
                   top=o[1:n]}
topn=apply(Correlation_abs,2,TopN,n=15)
##Predict with selected groups
pred=function(x,w){sum(c(1,x)*w)}
Y_p=matrix(0,2238,10000)
for (k in c(1:2238)){
  x=t(RNA_st_s1[topn[,k],])
  y_p=apply(x,1,pred,w=W[k,])
  Y_p[k,]=y_p
  print(k)
}
##Recover prediction
miRNA_pred=matrix(0,2238,10000)
for(q in c(1:nrow(Y_p))){
  miRNA_pred[q,]=Y_p[q,]*miRNA_train_sd[q]+miRNA_train_mean[q]
}
miRNA_pred_10=matrix(0,nrow(miRNA_10),ncol(miRNA_10))
for (l in c(0:9)) {
  miRNA_pred_10[,(l+1)]=apply(miRNA_pred[,(l*1000+1):((l+1)*1000)],1,mean)
}

##Record the result
Cor_row=rep(0,nrow(miRNA_10))
Cor_col=rep(0,ncol(miRNA_10))
Cor_train_mean=rep(0,ncol(miRNA_10))
for(r1 in c(1:nrow(miRNA_10))){
  Cor_row[r1]=cor.test(miRNA_10[r1,],miRNA_pred_10[r1,])$estimate
}
for (r2 in c(1:ncol(miRNA_10))) {
  Cor_col[r2]=cor.test(miRNA_10[,r2],miRNA_pred_10[,r2])$estimate
}
for (r3 in c(1:ncol(miRNA_10))) {
  Cor_train_mean[r3]=cor.test(miRNA_10[,r3],miRNA_train_mean)$estimate
}
result1=Cor_row
result2=Cor_col
result3=Cor_train_mean

