if(!require(polynom))install.packages("polynom");library(polynom)
if(!require(ggplot2))install.packages("ggplot2");library(ggplot2)

#LEAST SQUARES APPROXIMATION
#least squares approximation with polynomial of degree m
least_sq<-function(x,y,m){
  dat<-data.frame(cbind(x,y))
  n<-length(y)
  X<-rep(1,n)
  for (i in 1:m){
    X<-cbind(X,x^i)
  }
  beta<-solve(t(X)%*%X)%*%t(X)%*%y
  beta<-round(beta,digits=5)
  lsq<-polynomial(beta)
  label<-paste0("p",paste(m),"(x)==",paste(lsq))
  print(as.character(lsq))
  lsq<-as.function(lsq)
  plot<-ggplot(dat, aes(x=x, y=y)) + ggtitle("Least Square Approximation Polynomial")+
    geom_point(size=5, col='green') + stat_function(fun = lsq, size=1.25, alpha=0.4)+
    annotate("label",x=-Inf,y=-Inf,hjust=0,vjust=-2,label = label,parse=TRUE)
  print(plot)
  ggsave(filename=paste(n,"least_sq.pdf"), plot=plot,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")
  return(lsq)
}

x<-c(1,2,3,4,5)
y<-c(5,4,6,3,2)
m<-4
least_sq(x,y,m)
