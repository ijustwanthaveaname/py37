#假设原来的表达矩阵fpkm_expr.txt中行为基因，列为样本，中间数值是FPKM计算得到的值
#先读取自己的表达矩阵
expMatrix<-read.table("fpkm_expr.txt",header = T,row.names = 1)
#原始counts转化为Tpm
countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}
#原始counts转化为FPKM
countToFpkm <- function(counts, effLen)
{
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}
#fpkm转为tpm
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
#counts转为effcounts
countToEffCounts <- function(counts, len, effLen)
{
  counts * (len / effLen)
}
#如果要计算TPM值，只需要用一下apply函数
tpms <- apply(expMatrix,2,fpkmToTpm)
tpms[1:3,]
#最后可以根据TPM的特征进行检查，看每列加和是否一致
colSums(tpms)