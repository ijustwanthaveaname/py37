rm(list=ls())
options(stringsAsFactors = F)
library(edgeR)
#计算TMM标准化因子
f=calcNormFactors(exprSet,method = "TMM")
#将表达矩阵进行TMM标准化
exprSet=t(t(exprSet)/f)
library(WGCNA)
gsg=goodSamplesGenes(exprSet,verbose = 1)
if(!gsg$allOK)
{
  # optionally,print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste('Removing genes:',paste(names(exprSet)[!gsg$goodGenes],collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste('Removing samples:',paste(rownames(exprSet)[!gsg$goodSamples],collapse=", ")))
  # Remove the offending genes and samples from the data:
  exprSet=exprSet[gsg$goodSamples,gsg$goodGenes]
}
##filter可不进行
meanFPKM=0.5 ##过滤标准，可以修改
n=nrow(exprSet)
exprSet=as.data.frame(exprSet)
exprSet[n+1,]=apply(exprSet[c(1:nrow(exprSet)),],2,mean)
exprSet=exprSet[1:n,exprSet[n+1,]>meanFPKM]

filtered_fpkm=t(exprSet)
filtered_fpkm=data.frame(rownames(filtered_fpkm),filtered_fpkm)
names(filtered_fpkm)[1]='sample'
write.table(filtered_fpkm,file='FPKM_filter.xls',row.names=F,col.names=T,quote=FALSE,sep='\t')
##############Sample_cluster#######样品聚类，可不进行#############研究样品的相关性或去除异常样品
sampleTree=hclust(dist(filtered_fpkm),method='average')
sizeGrWindow(12,9)
#pdf('hclust.pdf',width=12,height=9)
par(cex=0.6)
par(mar=c(0,4,2,0))
plot(sampleTree,main="Sample clustering to detect outliers",sub = "",xlab="",cex.lab=1.5,
     cex.axis=1.5,cex.main=2)
#dev.off()
#####plot a line to show the cut
abline(h=8000,col='red')###设置剪切高度红线，保留红线以下
####Determine cluster under the line 
clust =cutreeStatic(sampleTree,cutHeight = 8000,minSize = 10)
table(clust)
keepsamples=(clust==1)
filtered_fpkm=filtered_fpkm[keepsamples]
###Loading clinical trait data
traitData=read.table('traits.txt',row.names=1,header=T,comment.char = '',check.names = F)
dim(traitData)
names(traitData)
#remove column that hold information we do not need
allTraits=traitData
dim(allTraits)
names(allTraits)
#Form a data frame analogous to expression data that will hod the clinical traits
fpkmSamples = rownames(filtered_fpkm)
traitSamples = rownames(allTraits)
traitRows = match(fpkmSamples,traitSamples)
datTraits=allTraits[traitRows,]
datTraits=as.data.frame(datTraits)
dim(datTraits)
rownames(datTraits)=as.data.frame(rownames(allTraits))[traitRows,]
collectGarbage()

##cluster samples
sampleTree2=hclust(dist(filtered_fpkm),method='average')
##covert traits to a color representation:white means low,red means high,grey means missing entry
traitColors=numbers2colors(datTraits,signed=FALSE)
# plot the sample dendrogram and the colors underneath
sizeGrWindow(12,12)
#pdf(file='Sample dendgrogram and trait heatmap.pdf',width=12,height=12)
plotDendroAndColors(sampleTree2,traitColors,
                    groupLabels = names(datTraits),
                    main='Sample dendrogram and trait heatmap')
#dev.off()
##允许多线程
enableWGCNAThreads()
##choose a set of soft-thresholding powers
powers=c(1:30)
#call the network topology analysis function
filtered_fpkm=filtered_fpkm[,-1]
sft=pickSoftThreshold(filtered_fpkm,powerVector = powers,verbose = 5)
#plot the results
sizeGrWindow(9,5)
#pdf(file='Scale independence.pdf',width=9,height=5)
par(mfrow=c(1,2))
cex1=0.9
#scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft threshold (power)",ylab = "Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
#this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col='red')
#Mean connectivity as a funtion of the soft-threshoding power
plot(sft$fitIndices[,1],sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab = "Mean Connectivity",type = "n",
     main=paste("Mean connectivity"))
text(sft$fitIndices[,1],sft$fitIndices[,5],labels=powers,cex = cex1,col='red')
#dev.off()
#####chose the softpower
softPower=sft$powerEstimate
adjacency=adjacency(filtered_fpkm,power = softPower)
### Turn adjacency into topological overlap
TOM=TOMsimilarity(adjacency)
dissTOM=1-TOM
#call the hierarchical cllustering function
geneTree=hclust(as.dist(dissTOM),method='average')
#plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
#pdf(file='Gene clustering on TOM-based dissimilarity.pdf',width=12,height=9)
plot(geneTree,xlab="",sub = "",main="Gene clustering of TOM-based dissimilarity",
     labels=FALSE,hang=0.04)
#dev.off()

#we like large modules ,so we set the minimum module size relatively high:
minmodulesize=30
#Module identification using dynamic tree cut:
dynamicMods=cutreeDynamic(dendro = geneTree,distM = dissTOM,
                          deepSplit = 2,pamRespectsDendro = FALSE,
                          minClusterSize = minmodulesize)
table(dynamicMods)

#Convert numeric lables into colors
dynamicColors=labels2colors(dynamicMods)
table(dynamicColors)
#plot the dendrogram and colors underneath
sizeGrWindow(8,6)
#pdf(file='Dynamic Tree Cut.pdf',width=8,height=6)
plotDendroAndColors(geneTree,dynamicColors,"Dynamic Tree Cut",
                    dendroLabels = FALSE,hang = 0.03,
                    addGuide = TRUE,guideHang = 0.05,
                    main="Gene dendrogram and module colors")
#dev.off()

#Calculate eigengene
MEList=moduleEigengenes(filtered_fpkm,colors=dynamicColors)
MEs=MEList$eigengenes
#Calculate dissimilarity of module eigengenes
MEDiss=1-cor(MEs)
#cluster module eigengenes
METree=hclust(as.dist(MEDiss),method = 'average')
# plot the result
sizeGrWindow(7,6)
#pdf(file='clustering of module eigengenes.pdf',width=7,height=6)
plot(METree,main='clustering of module eigengenes',
     xlab='',sub='')
MEDissThres=0.25#####剪切高度可修改
# plot the cut line into the dendrogram
abline(h = MEDissThres,col='red')
#dev.off()
# Call an automatic merging function
merge=mergeCloseModules(filtered_fpkm,dynamicColors,cutHeight = MEDissThres,verbose=3)
#The merged modue colors
mergedColors=merge$colors
#Eigengenes of the new merged modules:
mergedMEs=merge$newMEs
sizeGrWindow(12,9)
#pdf(file='merged dynamic.pdf',width=9,height=6)
plotDendroAndColors(geneTree,cbind(dynamicColors,mergedColors),
                    c('Dynamic Tree Cut','Merged dynamic'),
                    dendroLabels = FALSE,hang=0.03,
                    addGuide = TRUE,guideHang = 0.05)
#dev.off()

#Rename to moduleColors
moduleColors=mergedColors
table(moduleColors)
#Construct numerical labels corresponding to the colors
colorOrder=c('grey',standardColors(50))
moduleLabels=match(moduleColors,colorOrder)-1
MEs=mergedMEs
table(moduleLabels)


#####################relate modules to external clinical triats########################
# Define numbers of genes and samples
nGenes = ncol(filtered_fpkm)
nSamples=nrow(filtered_fpkm)

moduleTraitCor=cor(MEs,datTraits,use = 'p')
moduleTraitPvalue=corPvalueStudent(moduleTraitCor,nSamples)

sizeGrWindow(10,6)
#pdf(file='Module-trait relationships.pdf',width=10,height=6)
#will display correlations and their p-values
textMatrix=paste(signif(moduleTraitCor,2),'\n(',
                 signif(moduleTraitPvalue,1),')',sep='')

dim(textMatrix)=dim(moduleTraitCor)
par(mar=c(6,8.5,3,3))

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix=moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors=greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim=c(-1,1),
               main=paste('Module-trait relationships'))
#dev.off()


####### Define variable weight containing all column of datTraits

###MM and GS

# names (colors) of the modules
modNames=substring(names(MEs),3)

geneModuleMembership=as.data.frame(cor(filtered_fpkm,MEs,use='p'))
MMPvalue=as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),nSamples))

names(geneModuleMembership)=paste('MM',modNames,sep='')
names(MMPvalue)=paste("p.MM",modNames,sep='')

#names of those trait
traitNames=names(filtered_fpkm)

geneTraitSignificance=as.data.frame(cor(filtered_fpkm,datTraits,use='p'))
GSPvalue=as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples))

names(geneTraitSignificance)=paste('GS.',traitNames,sep='')
names(GSPvalue)=paste('p.GS',traitNames,sep='')

###plot MM vs GS for each trait vs each module

#####example:royalblue and CK
module='royalblue'
column=match(module,modNames)
moduleGenes= moduleColors==module

trait='CK'
traitColumn=match(trait,traitNames)

sizeGrWindow(7,7)
par(mfrow=c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
                   abs(geneTraitSignificance[moduleGenes,traitColumn]),
                   xlab=paste('Module membership in',module,'module'),
                   ylab=paste('Gene significance for ',trait),
                   main=paste('Module membership vs. gene significance\n'),
                   cex.main = 1.2,cex.lab = 1.2,cex.axis = 1.2,col=module)
#####
names(filtered_fpkm)
probes=names(filtered_fpkm)

#######################export GS and MM########################
geneInfo0=data.frame(probes=probes,
                     moduleColor=moduleColors)

for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames=names(geneInfo0)
  geneInfo0=data.frame(geneInfo0,geneTraitSignificance[,Tra],
                       GSPvalue[,Tra])
  names(geneInfo0)=c(oldNames,names(geneTraitSignificance)[Tra],
                     names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames=names(geneInfo0)
  geneInfo0=data.frame(geneInfo0,geneModuleMembership[,mod],
                       MMPvalue[,mod])
  names(geneInfo0)=c(oldNames,names(geneModuleMembership)[mod],
                     names(MMPvalue)[mod])
}
geneOrder=order(geneInfo0$moduleColor)
geneInfo=geneInfo0[geneOrder,]

write.table(geneInfo,file='GS and MM.xls',sep='\t',rownames=F)

####################visualizing the gene network############################

nGenes=ncol(filtered_fpkm)
nSamples=nrow(filtered_fpkm)

# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM=dissTOM^7
#Set diagonal to NA for a nicer plot
diag(plotTOM)=NA

# Call the plot function

sizeGrWindow(9,9)
#pdf(file='Network hatmap plot all gene.pdf',width=9,height=9)
TOMplot(plotTOM,geneTree,moduleColors,main='Netwrok heatmap plot,all genes')
#dev.off()

nSelect=400#选400个基因，同上步，只是选择少部分基因。
#FOr reproducibility,we set the random seed
set.seed(10)
select=sample(nGenes,size=nSelect)
selectTOM=dissTOM[select,select]
# There's no simple way of restricting a clustering tree to a subset of genes,so we must re-cluster
selectTree=hclust(as.dist(selectTOM),method='average')
selectColors=moduleColors(select)

# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissmilarity to a power,say 10,makes the plot more informative by effectively changing
# the color palette:setting the diagonal to NA also improvves the clarity of the plot
plotDiss =selectTOM^7
diag(plotDiss)=NA
#pdf(file='Network heatmap plot_selected genes.pdf',width=9,height=9)
TOMplot(plotDiss,selectTree,selectColors,main='Network heatmap plot,selected genes')
#dev.off()

###############visualizing the gene network of eigengenes############################

sizeGrWindow(5,7.5)
#pdf(file='Eigengene dendrogram and Eigengene adjacency heatmap.pdf',width=5,heigth=7.5)
par(cex=0.9)
plotEigengeneNetworks(MEs,'',marDendro = c(0,4,1,2),marHeatmap = c(3,4,1,2),cex.lab=0.8,xLabelsAngle=90)
#dev.off()

#将以上结果分别绘制
#or devide into two parts
sizeGrWindow(6,6)
#pdf(file='Eigengene dendrogram_2.pdf',width=6,height=6)
par(cex=1.0)
plotEigengeneNetworks(MEs,'Eigengene dendrogram',marDendro = c(0,4,2,0),plotHeatmaps = FALSE)
#dev.off()

#pdf(file='Eigengene adjacency heatmap_2.pdf',width=6,height=6)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex=1.0)
plotEigengeneNetworks(MEs,'Eigengene adjacency heatmap',marHeatmap = c(3,4,2,2),plotDendrograms = FALSE,xLabelsAngle=90)
#dev.off()



####################Exporting to Cytoscape all one by one ########################
#Select each module
for (mod in 1:nrow(table(moduleCOlors)))
{
  modules=names(table(moduleColors))[mod]
  # Select module probes
  probes = names(filtered_fpkm)
  inModule=(moduleColors==modules)
  modProbes=probes[inModule]
  modGenes=modProbes
  # Select the corresponding Topological Overlap
  modTOM=TOM[inModule,inModule]
  
  dimnames(modTOM)=list(modProbes,modProbes)
  # Export the network into edge and node list files Cytoscape can read
  cyt=exportNetworkToCytoscape(modTOM,
                               edgeFile = paste('CytoscapeInput-edges-',modules,'.txt',sep = ''),
                               nodeFile = paste('CytoscapeInput-nodes-',modules,'.txt',sep=''),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr=moduleColors[inModule])
}
