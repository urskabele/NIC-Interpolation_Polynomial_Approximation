if(!require(interp))install.packages("interp");library(interp)

#PIECEWISE LINEAR INTERPOLATION
piece_lin<-function(x,y){
  dat<-data.frame(cbind(x,y)) #create data frame
  af <- approxfun(x,y)
  plot<-ggplot(dat, aes(x=x, y=y)) + ggtitle("Piecewise Linear Interpolation")+
    geom_point(size=5, col='red') + stat_function(fun = af, size=1.25, alpha=0.4)
  print(plot)
  ggsave(filename=paste("linear.pdf"), plot=plot,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")
  assign("af",af, envir = .GlobalEnv)
  return(invisible())
}

x<-c(1:4)
y<-c(1,3,-1,2)

piece_lin(x,y)
