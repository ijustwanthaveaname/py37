options(digits = 1)
library(VIM)
sleep[!complete.cases(sleep),]#显示一个或多个缺失值的行

sum(is.na(sleep$Dream))
mean(is.na(sleep$Dream))
mean(!complete.cases(sleep))
library(mice)
md.pattern(sleep)#生成缺失值模式表格

aggr(sleep,prop=FALSE,numbers=TRUE)
aggr(sleep,prop=TRUE,numbers=TRUE)

matrixplot(sleep)

marginplot(sleep[c("Gest","Dream")],pch=C(20),
           col=c("darkgray","red","blue"))

#用相关性探索缺失值
x<-as.data.frame(abs(is.na(sleep)))
head(sleep,n=5)
head(x,n=5)
y<-x[which(apply(x,2,sum)>0)]
cor(y)

cor(sleep,y,use="pairwise.complete.obs")

#多重插补(MI)
library(mice)
imp<-mice(sleep,seed=1234)
fit<-with(imp,lm(Dream~Span+Gest))
pooled<-pool(fit)
summary(pooled)
imp
imp$imp$Dream #展示Dream变量上的插补信息
dataset3<-complete(imp,action=3)#选择多重插补中的第三个完整数据集
