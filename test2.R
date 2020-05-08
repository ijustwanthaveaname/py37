rm(list=ls())
options(stringsAsFactors = F)
library(airway)
data(airway)
exprSet=assay(airway)
group_list=colData(airway)[,3]
library(DESeq2)
colData=data.frame(row.names=colnames(exprSet),
                   group_list=group_list)
dds=DESeqDataSetFromMatrix(countData = exprSet,
                           colData=colData,
                           design = ~group_list)
dds=DESeq(dds)
res=results(dds,
            contrast=c('group_list','trt','untrt'))
resOrdered=res[order(res$padj),]
head(resOrdered)
DEG=as.data.frame(resOrdered)
