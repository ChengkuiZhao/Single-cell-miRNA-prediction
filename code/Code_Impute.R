library(Seurat)
library(SAVER)
pbmc1.saver <- saver(pbmc1.data, ncores = 64, estimates.only = TRUE)
save(pbmc1.saver,file='1_fresh_imputed.RData')


data=matrix(0,20530,10000)
for (i in c(1:10)){
      d=data_var0.2[,((i-1)*1000+1):(i*1000)]
      data[,((i-1)*1000+1):(i*1000)]=saver(d,ncores=64,estimates.only=TRUE)
      print(i)
}