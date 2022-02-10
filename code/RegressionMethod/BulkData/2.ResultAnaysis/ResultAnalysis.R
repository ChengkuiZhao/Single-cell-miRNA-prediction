library(ggplot2)
load(file='miRNA_train.Rdata')
load(file='C_mi_Cor.Rdata')
load(file='C_P_Cor.Rdata')
load(file='Train_Mean_Cor.Rdata')
miRNA_train_mean=apply(miRNA_train,1,mean)
o=order(miRNA_train_mean,decreasing = T)
boxplot(Train_mean_Cor,C_mi_Cor,names = c('Mean','Predict'),col=c('lightblue','green'),ylab='Cross miRNA Correlation',main='Cross miRNAs Correlation(128 Patients, 2238 miRNAs)')
text(c(0.7,1.7),c(0.992,0.992),c(0.978,0.983))
boxplot(0,C_P_Cor,names =c('Mean','Predict'),col=c('lightblue','green'),ylab='Cross Patients Correlation',main='Cross Patients Correlation(128 Patients, 2238 miRNAs)')
median(Cor_row,na.rm = T)
text(c(0.8,1.8),c(0.1,0.45),c(0,0.134))
Noise=C_P_Cor[apply(miRNA_train,1,function(x)(all(x<1)))]
Nonnoise=C_P_Cor[apply(miRNA_train,1,function(x)(any(x>1)))]
boxplot(Noise,Nonnoise,names =c('Noise','NonNoise'),col = c('red','green'),ylab='Cross Patients Correlation')
median(Noise,na.rm = T)
median(Nonnoise,na.rm = T)
boxplot(Noise,Nonnoise,names =c('Noise','NonNoise'),col = c('red','green'),ylab='Cross Patients Correlation',main='Filter Noise')
text(c(0.8,1.8),c(0.1,0.5),c(-0.008,0.222))
NonNoise_index=apply(miRNA_train,1,function(x)(any(x>1)))
miRNA_train_NonNoise=miRNA_train[NonNoise_index,]
miRNA_train_NonNoise_mean=apply(miRNA_train_NonNoise,1,mean)
C_P_Cor_NonNoise=C_P_Cor[rownames(miRNA_train_NonNoise)]
Expression_high=C_P_Cor_NonNoise[which(miRNA_train_NonNoise_mean>1)]
Expression_low=C_P_Cor_NonNoise[which(miRNA_train_NonNoise_mean<1)]
boxplot(Expression_low,Expression_high,names =c('Low','High'),col = c('red','green'),ylab='Cross Patients Correlation',main='Influenced by Expression')
median(Expression_high,na.rm = T)
median(Expression_low,na.rm = T)
text(c(0.8,1.8),c(0.3,0.6),c(0.123,0.485))
spread=apply(miRNA_train_NonNoise,1,function(x)(diff(range(x))))
qplot(spread)
Spread_low=C_P_Cor_NonNoise[which(spread<summary(spread)[2])]
Spread_high=C_P_Cor_NonNoise[which(spread>=summary(spread)[2])]
boxplot(Spread_low,Spread_high,names =c('Low','High'),col =c('red','green'),ylab='Cross Patients Correlation',main='Influenced by Spread')
median(Spread_low,na.rm = T)
median(Spread_high,na.rm = T)
text(c(0.8,1.8),c(0.15,0.55),c(0.029,0.318))
CV=apply(miRNA_train_NonNoise,1,function(x)(sd(x)/mean(x)))
CV_high=C_P_Cor_NonNoise[which(CV>=summary(CV)[2])]
CV_low=C_P_Cor_NonNoise[which(CV<summary(CV)[2])]
boxplot(CV_low,CV_high,names =c('Low','High'),col =c('green','red'),ylab='Cross Patients Correlation',main='Influenced by CV')
median(CV_low,na.rm = T)
median(CV_high,na.rm = T)
text(c(0.8,1.8),c(0.6,0.4),c(0.49,0.139))
boxplot(Noise,Expression_low,Expression_high,Spread_low,Spread_high,CV_low,CV_high,names = c('Noise','Expression_low','Expression_high','Spread_low','Spread_high','CV_low','CV_high'),col = c('red','lightblue','green','lightblue','green','green','lightblue'),ylab='Cross Patients Correlation',main='Influence Factor')
text(c(0.8:6.8),c(0.1,0.3,0.6,0.15,0.55,0.6,0.35),labels = c(-0.008,0.123,0.485,0.029,0.318,0.49,0.139))
miRNA_selected=C_P_Cor_NonNoise[Reduce(intersect,list(names(Expression_high),names(Spread_high),names(CV_low)))]
miRNA_filtered=C_P_Cor[setdiff(names(C_P_Cor),names(miRNA_selected))]
boxplot(C_P_Cor,miRNA_filtered,miRNA_selected,names = c('All miRNAs','Filtered miRNA','Selected miRNA'),col=c('brown','cornflowerblue','green'),ylab='Cross Patients Correlation',main='Comparison')
round(median(C_P_Cor,na.rm = T),3)
round(median(miRNA_filtered,na.rm = T),3)
round(median(miRNA_selected,na.rm = T),3)
text(c(0.8,1.8,2.8),c(0.4,0.4,0.6),labels = c(0.111,0.086,0.248))
miRNAname_selected=names(miRNA_selected)
write.csv(miRNAname_selected,'miRNAname_selected')












