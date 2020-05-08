library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(pathview)
library(AnnotationHub)
library(biomaRt)
library(org.Hs.eg.db)
hub=AnnotationHub::AnnotationHub()
query(hub,"Gekko")
Gekko.orgdb=hub[["AH76465"]]
upgene=rownames(DEG[DEG$Change=='UP',])
downgene=rownames(DEG[DEG$Change=='DOWN',])
#ÊäÈësymbolnameµÄdataframe
entrezid=bitr(upgene,fromType="SYMBOL",
              toType ="ENTREZID",
              OrgDb = Gekko.orgdb)
geneid=na.omit(entrezid)#É¾³ýÈ±Ê§Öµ
goenrich=enrichGO(gene=geneid,
                  OrgDb=Gekko.orgdb,
                  keyType = "ENTREZID",
                  ont="BP",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2)
#dotplot(goenrich)
#barplot(goenrich)
#keggenrich=enrichKEGG(gene=geneid,
#                      organism = gja,
#                      keyType = "kegg",
#                      pvalueCutoff = 0.05,
#                      qvalueCutoff = 0.2)
#dotplot(keggenrich)
#barplot(keggenrich)