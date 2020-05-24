wssplot<-function(data,nc=15,seed=1234){
  wss<-(nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i]<-sum(kmeans(data,centers=i)$withinss)}
    plot(1:nc,wss,type="b",xlab = "Number of Clusters",
         ylab="Within groups sum of squares")
}