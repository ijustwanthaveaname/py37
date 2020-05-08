shrinkage<-function(fit,k=10){
  require(bootstrap)
  
  theta.fit<-function(x,y){lsfit(x,y)}
  theta.predict<-function(fit,x){cbind(1,x)%*%fit$coef}
  
  x<-fit$model[,2:ncol(fit$model)]
  y<-fit$model[,1]
  
  results<-crossval(x,y,theta.fit,theta.predict,ngroup=k)
  r2<-cor(y,fit$fitted.values)^2
  r2cv<-cor(y,results$cv.fit)^2
  cat('Original R-square =',r2,'\n')
  cat(k,'Fold Cross-Validated R-square =',r2cv,'\n')
  cat('Change =',r2-r2cv,'\n')
}

states<-as.data.frame(state.x77[,c('Murder','Population',"Illiteracy",'Income','Frost')])
fit<-lm(Murder~Population+Income+Illiteracy+Frost,data=states)
shrinkage(fit)

fit2<-lm(Murder~Population+Illiteracy,data=states)
shrinkage(fit2)
