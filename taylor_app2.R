if(!require(pracma))install.packages("pracma");library(pracma)
if(!require(ggplot2))install.packages("ggplot2");library(ggplot2)
if(!require(polynom))install.packages("polynom");library(polynom)

#TAYLOR APPROXIMATION
#Local polynomial approximation through Taylor series, returns taylor polynomial 
#ADD LEGEND AND EXPRESSION OF TAYL POLY
tayl_app<-function(f,a,k){
  x<-a
  y<-f(a)
  dat<-data.frame(cbind(x,y))
  coefi<-rev(taylor(f,a,n=k))
  coefi<-round(coefi,digits=5)
  g<-polynomial(coefi)
  label<-paste0("T",paste(k),"(f,",paste(a),")(x)==",paste(g))
  print(as.character(g))
  g<-as.function(g)
  plot<-ggplot(dat, aes(x=x, y=y)) + ggtitle("Taylor Polynomial Approximation")+
    geom_point(size=4, col='blue') +
    stat_function(fun = f, size=1.25,aes(colour = "f(x)")) + 
    stat_function(fun = g, size=1.25,alpha=0.4,aes(colour = "T(x)"))+
    #scale_colour_manual("Lgend title", values = c("purple", "blue"))+
    xlim(a-5,a+5)+
    annotate("label",x=-Inf,y=Inf,hjust=0,vjust=2,label = label,parse = TRUE)+
    scale_colour_manual("", values = c("#767676","red"))
  print(plot)
  ggsave(filename=paste("taylor.pdf"), plot=plot,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")
  assign("taylor_poly",g, envir = .GlobalEnv)
  return(invisible())
}

tayl_app(sin,0,5)

#tayl_app(exp,0,6)
