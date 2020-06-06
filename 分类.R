loc<-"http://archive.ics.uci.edu/ml/machine-learning-databases/"
ds<-"breast-cancer-wisconsin/breast-cancer-wisconsin.data"
url<-paste(loc,ds,sep="")

breast<-read.table(url,sep=',',header=FALSE,na.strings='?')
names(breast)<-c('ID','clumpThickness','sizeUniformity',
                 'shapeUniformity','maginalAdhesion',
                 'singleEpithelialCellSize','bareNuclei',
                 'blandChromatin','normalNucleoli','mitosis','class')

df<-breast[-1]
df$class<-factor(df$class,levels=c(2,4),
                 labels=c("benign","malignant"))

set.seed(1234)
train<-sample(nrow(df),0.7*nrow(df))
df.train<-df[train,]
df.validate<-df[-train,]
table(df.train$class)
table(df.validate$class)

####logistic regression
fit.logit<-glm(class~.,data=df.train,family=binomial())
summary(fit.logit)
prob<-predict(fit.logit,df.validate,type='response')
logit.pred<-factor(prob > .5,levels=c(FALSE,TRUE),
                   labels = c("benign","malignant"))
logit.pref<-table(df.validate$class,logit.pred,
                  dnn=c("Actual","Predicted"))
logit.pref

logit.fit.reduced<-step(fit.logit)

####Decision Tree
library(rpart)
set.seed(1234)
dtree<-rpart(class~.,data=df.train,method="class",
             parms=list(split="information"))
dtree$cptable
plotcp(dtree)
dtree.pruned<-prune(dtree,cp=.0125)
library(rpart.plot)
prp(dtree.pruned,type=2,extra=104,
    fallen.leaves=TRUE,main="Deciosion Tree")
dtree.pred<-predict(dtree.pruned,df.validate,type="class")
dtree.perf<-table(df.validate$class,dtree.pred,
                  dnn=c("Actual","Predicted"))
dtree.perf

####conditional inference tree
library(party)
fit.ctree<-ctree(class~.,data=df.train)
plot(fit.ctree,main="Conditional Inference Tree")
ctree.pred<-predict(fit.ctree,df.validate,type="response")
ctree.perf<-table(df.validate$class,ctree.pred,dnn=c("Actual","Predicted"))
ctree.perf

####random forest
library(randomForest)
set.seed(1234)
fit.forest<-randomForest(class~.,data=df.train,na.action=na.roughfix,importance=TRUE)
fit.forest
importance(fit.forest,type=2)
forest.pred<-predict(fit.forest,df.validate)
forest.pref<-table(df.validate$class,forest.pred,dnn=c("Actual","Predicted"))
forest.perf

####SVM
library(e1071)
set.seed(1234)
fit.svm<-svm(class~.,data=df.train)
fit.svm
svm.pred<-predict(fit.svm,na.omit(df.validate))
svm.pref<-table(na.omit(df.validate)$class,
                svm.pred,dnn=c("Actual","Predicted"))
svm.pref

####RBF-SVM
library(e1071)
set.seed(1234)
tuned<-tune.svm(class~.,data=df.train,
                gamma=10^(-6:1),
                cost=10^(-10:10))
tuned
fit.svm<-svm(class~.,data=df.train,gamma=.01,cost=1)
svm.pred<-predict(fit.svm,na.omit(df.validate))
svm.perf<-table(na.omit(df.validate)$class,
                svm.pred,dnn=c("Actual","Predicted"))
svm.perf

####Accuracy evaluation
performance<-function(table,n=2){
  if(!all(dim(table)==c(2,2)))
    stop("Must be a 2x2 table")
  tn = table[1,1]
  fp = table[1,2]
  fn = table[2,1]
  tp = table[2,2]
  sensitivity = tp/(tp+fn)#正类样本被成功预测的概率
  specificity = tn/(tn+fp)#负类样本被成功预测的概率
  ppp = tp/(tp+fp)#被预测为正类的样本单元中，预测正确的样本单元占比
  npp = tn/(tn+fn)#被预测为负累的样本单元中，预测正确的样本单元占比
  hitrate=(tp+tn)/(tp+tn+fp+fn)#准确率
  result<-paste("Sensitivity = ",round(sensitivity,n),
                "\nSpecificity = ", round(specificity,n),
                "\nPositive Predictive Value = ",round(ppp,n),
                "\nNegative Predictive Value = ",round(npp,n),
                "\nAccuracy = ",round(hitrate,n),"\n",sep="")
  cat(result)
}
####Accuracy of each model
performance(logit.perf)
performance(dtree.perf)
performance(ctree.perf)
performance(forest.perf)
performance(svm.perf)