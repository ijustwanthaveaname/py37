rm(list=ls())
options(stringsAsFactors = F)
counts=read.table('Count.txt',header=T)
rownames(counts)=counts[,1]
counts=counts[,-1]
meta=counts[,1:5]
exprSet=counts[,6:ncol(counts)]
library(corrplot)#相关性图
png('cor.png')
corrplot(cor(log2(exprSet+1)))
dev.off()
library(pheatmap)#热图
png('heatmap.png')
m=cor(log2(exprSet+1))#变为相关性矩阵
pheatmap(scale(cor(log2(exprSet+1))))#归一化后画热图
dev.off()
###以上为矩阵的基本探索