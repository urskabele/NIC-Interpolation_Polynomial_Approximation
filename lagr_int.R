if(!require(ggplot2))install.packages("ggplot2");library(ggplot2)
if(!require(rSymPy))install.packages("rSymPy");library(rSymPy)
if(!require(DeducerExtras))install.packages("DeducerExtras");library(DeducerExtras)
if(!require(polynom))install.packages("polynom");library(polynom)


#FUNCTION THAT RETURNS LAGRANGE INTERPOLATION POLYNOMIAL AND PLOT
lagrange<-function(x,y){
  n<-length(x)-1
  dat<-data.frame(cbind(x,y)) #create data frame
  int_pol<-poly.calc(x, y) #Lagrange interpolation polynomial
  coefi<-round(coef(int_pol),digits=5)
  int_pol<-polynomial(coefi)
  print(int_pol)  
  label<-paste0("L",paste(n),"(x)==",paste(int_pol))
  int_pol<-as.function(int_pol)
  plot<-ggplot(dat, aes(x=x, y=y)) + ggtitle("Interpolation Polynomial in Lagrange Form")+
    geom_point(size=5, col='blue') + stat_function(fun = int_pol, size=1.25, alpha=0.4)+
      annotate("label",x=-Inf,y=Inf,hjust=0,vjust=2,label = label,parse = TRUE)
  print(plot)
  ggsave(filename=paste(n,"lagrange.pdf"), plot=plot,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")
  assign("int_pol",int_pol, envir = .GlobalEnv)
  return(invisible())
  }

x <- c(0, 2, 3, 4)
y <- c(7,76,-3,32)
# lagrange(x,y)

lagrange(x,y)
