rm(list=ls())
options(stringsAsFactors = F)
library(airway)
data(airway)
exprSet=assay(airway)#得到表达矩阵
group_list=colData(airway)[,3]#分组信息
###用DESeq2来寻找差异表达的基因
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
DEG=na.omit(DEG)
DEG=DEG[substring(rownames(DEG),1,6)!='LOC107',]#处理壁虎预测基因
###画热图
library(pheatmap)
choose_gene=head(rownames(DEG),100)##50 maybe better
choose_matrix=exprSet[choose_gene,]
choose_matrix=t(scale(t(choose_matrix)))
pheatmap(choose_matrix,filename='DEG_top100_heatmap.png',fontsize=4)
###绘制火山图
library(ggplot2)
logFC_cutoff=with(DEG,mean(abs(log2FoldChange))+2*sd(abs(log2FoldChange)))
#logFC_cutoff=1
DEG$Change=as.factor(ifelse(DEG$pvalue<0.05&abs(DEG$log2FoldChange)>logFC_cutoff,
                            ifelse(DEG$log2FoldChange>logFC_cutoff,'UP','DOWN'),'NOT')
)
this_tile=paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                 '\nThe number of up gene is ',nrow(DEG[DEG$Change=='UP',]),
                 '\nThe number of down gene is ',nrow(DEG[DEG$Change=='DOWN',])
)

g=ggplot(data=DEG,
         aes(x=log2FoldChange,y=-log10(pvalue),
             color=Change))+
  geom_point(alpha=0.4,size=1.75)+
  theme_set(theme_set(theme_bw(base_size = 20)))+
  xlab('log2 fold change')+ylab('-log10 p-value')+
  ggtitle(this_tile)+theme(plot.title = element_text(size=15,hjust=0.5))+
  scale_colour_manual(values=c('blue','black','red'))##corresponding to the level
print(g)
ggsave(g,filename = 'volcano.png')
###用edgeR来寻找差异表达的基因
library(edgeR)
dge=DGEList(counts=exprSet,group = factor(group_list))
dge$samples$lib.size=colSums(dge$counts)#样本的基因总计数
dge=calcNormFactors(dge)#计算标准化的因子
design=model.matrix(~0+factor(group_list))
rownames(design)=colnames(dge)
colnames(design)=levels(factor(group_list))
dge=estimateGLMCommonDisp(dge,design)#计算普通的离散度
dge=estimateGLMTrendedDisp(dge,design)
dge=estimateGLMTagwiseDisp(dge,design)#计算基因间范围内的离散度
fit=glmFit(dge,design)#拟合广义线性模型
lrt=glmLRT(fit,contrast = c(1,0))#利用似然比检验来检验差异表达基因
nrDEG=topTags(lrt,n=nrow(exprSet)) #输出靠前的差异表达基因
nrDEG=as.data.frame(nrDEG)
###经典edgeR的使用方法
x<-read.delim("fileofcounts.txt",row.names="Symbol")# 读取reads count文件

group<-factor(c(1,1,2,2))#分组变量 前两个为一组， 后一个为一组， 每个有两个重复

y<-DGEList(counts=x,group=group)# 构建基因表达列表

y<-calcNormFactors(y)# 计算样本内标准化因子

y<-estimateCommonDisp(y)#计算普通的离散度

y<-estimateTagwiseDisp(y)#计算基因间范围内的离散度

et<-exactTest(y)# 进行精确检验

topTags(et)# 输出排名靠前的差异表达基因信息

###广义线性模型edgeR(可用于无重复样本)
rawdata<-read.delim("TableS1.txt",check.names=FALSE,stringsAsFactors=FALSE);  #读取原始文件

y<-DGEList(counts=rawdata[,4:9],genes=rawdata[,1:3]);  # 建立基因表达列表

y$samples$lib.size<-colSums(y$counts);  #设置各个样本的库大小

y<-calcNormFactors(y);  #利用库的大小来进行标准化TMM

Patient<-factor(c(8,8,33,33,51,51));  #因素1

Tissue<-factor(c("N","T","N","T","N","T"));#因素2

data.frame(Sample=colnames(y),Patient,Tissue);

design<-model.matrix(~Patient+Tissue);  #建立分组变量

rownames(design)<-colnames(y);

y<-estimateGLMCommonDisp(y,design,verbose=TRUE);  #估计离散度

y<-estimateGLMTrendedDisp(y,design);  

y<-estimateGLMTagwiseDisp(y,design);

fit<-glmFit(y,design);  #拟合广义线性模型

lrt<-glmLRT(fit) #利用似然比检验来检验差异表达基因

topTags(lrt)    #输出靠前的差异表达基因

