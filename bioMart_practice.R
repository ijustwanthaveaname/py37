rm(list=ls())
options(stringsAsFactors = FALSE)
library(biomaRt)
#We need to choose a BioMart database
listMarts()
#Choices include ensembl, vega, unimart, or many others
ensembl<-useMart('ensembl')
listDatasets(ensembl)
#We can browse the datasets and select human
ensembl<-useDataset('hsapiens_gene_ensembl',mart=ensembl)
filters<-listFilters(ensembl)
#Look at the first seven rows of filters
#then at the last few rows with the tail function
filters[1:7,]
tail(filters)
attributes<-listAttributes(ensembl)
attributes[1:5,]
myentrez<-c('3043','3045','3046','3047','3051')
mydata<-getBM(attributes = c('entrezgene','hgnc_symbol','percentage_gc_content'),
filters='entrezgene',values=myentrez,mart=ensem)
