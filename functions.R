if(!require(ggplot2))install.packages("ggplot2");library(ggplot2)
if(!require(rSymPy))install.packages("rSymPy");library(rSymPy)
if(!require(DeducerExtras))install.packages("DeducerExtras");library(DeducerExtras)
if(!require(polynom))install.packages("polynom");library(polynom)
if(!require(pracma))install.packages("pracma");library(pracma)
if(!require(interp))install.packages("interp");library(interp)

#LAGRANGE INTERPOLATION POLYNOMIAL
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
    xlim(a-7,a+7)+
    annotate("label",x=-Inf,y=Inf,hjust=0,vjust=2,label = label,parse = TRUE)+
    scale_colour_manual("", values = c("#767676","red"))
  print(plot)
  ggsave(filename=paste("taylor.pdf"), plot=plot,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")
  assign("taylor_poly",g, envir = .GlobalEnv)
  return(invisible())
}

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
  assign("lsq",lsq, envir = .GlobalEnv)
  return(invisible())
}

################################################################################################################
#TRY WITH DATA
f<-function(x){sin(2*x)}
x<-c(-7,-5,-2.3,-0.89,4,6.34)
y<-f(x)
dat<-data.frame(cbind(x,y))

#CALCULATE ALL APPROXIMATIONS
lagrange(x,y)
piece_lin(x,y)
cub_spline(x,y)
tayl_app(f,0,5)
least_sq(x,y,5)

#PLOT
plot_tog<-ggplot(dat, aes(x=x, y=y)) + ggtitle("Interpolation and Approximation")+geom_point(size=5, col='red')+
  stat_function(fun = int_pol, size=0.8, aes(colour="Lagr. int."))+
  stat_function(fun = af, size=0.8, aes(colour="Piec. lin."))+
  stat_function(fun = spl, size=0.8, aes(colour="Cub. spl."))+
  stat_function(fun = lsq, size=0.8, aes(colour="Least Sq."))+
  stat_function(fun = f, size=0.8, aes(colour="f(x)"))+
  scale_colour_manual("", values = c("blue","red","#CC00CC","#009999","#00CC00"))
print(plot_tog)

ggsave(filename="together.pdf", plot=plot_tog,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")

############################################################################################################### 
#TAYLOR DIFFERENT 
f<-function(x){sin(2*x)}

tayl_app(f,0,3)
tay_poly_3<-taylor_poly
tayl_app(f,0,5)
tay_poly_5<-taylor_poly
tayl_app(f,0,7)
tay_poly_7<-taylor_poly
x<-0
y<-f(0)
dat<-data.frame(cbind(x,y))

plot_tay<-ggplot(dat, aes(x=x, y=y)) + ggtitle("Taylor Polynomial Approximation")+geom_point(size=5, col='red')+
  stat_function(fun = tay_poly_3, size=0.8, aes(colour="T3(x)"))+
  stat_function(fun = tay_poly_5, size=0.8, aes(colour="T5(x)"))+
  stat_function(fun = tay_poly_7, size=0.8, aes(colour="T7(x)"))+
  stat_function(fun = f, size=0.8, aes(colour="f(x)"))+
  xlim(-3,3)+
  scale_colour_manual("", values = c("red","blue","#CC00CC","#009999"))
print(plot_tay)

ggsave(filename="taylor_diff.pdf", plot=plot_tay,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")

##############################################################################################################
#LEAST SQUARES DIFFERENT
x<-c(1,2,3,4,5,6,7)
y<-c(5,4,6,3,2,0,-1)
dat<-data.frame(cbind(x,y))

least_sq(x,y,1)
least_1<-lsq
least_sq(x,y,2)
least_2<-lsq
least_sq(x,y,3)
least_3<-lsq

plot_lsq<-ggplot(dat, aes(x=x, y=y)) + ggtitle("Least Square Approximation")+geom_point(size=5, col='green')+
  stat_function(fun = least_1, size=0.8, aes(colour="p1(x)"))+
  stat_function(fun = least_2, size=0.8, aes(colour="p2(x)"))+
  stat_function(fun = least_3, size=0.8, aes(colour="p3(x)"))+
  #stat_function(fun = f, size=0.8, aes(colour="f(x)"))+
  scale_colour_manual("", values = c("blue","#CC00CC","#009999"))
print(plot_lsq)

ggsave(filename="lsq_diff.pdf", plot=plot_lsq,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")





























