if(!require(ggplot2))install.packages("ggplot2");library(ggplot2)

#CUBIC SPLINE INTERPOLATION
cub_spline<-function(x,y){
  dat<-data.frame(cbind(x,y)) #create data frame
  spl <- splinefun(x, y, method="fmm")
  plot<-ggplot(dat, aes(x=x, y=y)) + ggtitle("Cubic Spline Interpolation")+
    geom_point(size=5, col='purple') + stat_function(fun = spl, size=1.25, alpha=0.4)
  print(plot)
  ggsave(filename=paste("spline.pdf"), plot=plot,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")
  assign("spl",spl, envir = .GlobalEnv)
  return(invisible())
}


x <- c(1:8)
y<- c(1,3,4,2,-2,-3,0,1)
cub_spline(x,y)
