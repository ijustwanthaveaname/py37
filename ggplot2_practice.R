##fig1
library(ggplot2)
ggplot(data=mtcars,aes(x=wt,y=mpg))+
  geom_point()+
  labs(title="Automobile Data",x="Weight",y="Miles Per Gallon")
##fig2
library(ggplot2)
ggplot(data=mtcars,aes(x=wt,y=mpg))+
  geom_point(pch=17,color="blue",size=2)+
  geom_smooth(method="lm",color="red",linetype=2)+
  labs(ttitle = "Automobile Data",x="Weight",y="Miles Per Gallon")
