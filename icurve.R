#I-CURVE

#data from 3/7/2018
maturities<-c(1/12,0.25,0.5,1,2,3,5,7,10,20,30)
yields<-c(1.91,1.98,2.12,2.33,2.53,2.63,2.72,2.79,2.83,2.89,2.96)
dat<-data.frame(cbind(maturities,yields))

#adapted function for graphs
lagrange_fix<-function(x,y){
  n<-length(x)-1
  dat<-data.frame(cbind(x,y)) #create data frame
  int_pol<-poly.calc(x, y) #Lagrange interpolation polynomial
  coefi<-coef(int_pol)
  int_pol<-polynomial(coefi)
  print(int_pol)  
  label<-paste0("L",paste(n),"(x)==",paste(int_pol))
  int_pol<-as.function(int_pol)
  plot<-ggplot(dat, aes(x=x, y=y)) + ggtitle("I-Curve: Lagrange Interpolation")+
    geom_point(size=5, col='blue') + stat_function(fun = int_pol, size=1.25, alpha=0.4)+
    labs(x = "Maturities [years]",y="Yields [%]")
    #annotate("label",x=-Inf,y=Inf,hjust=0,vjust=2,label = label,parse = TRUE)
  print(plot)
  ggsave(filename=paste("lagrange_ic.pdf"), plot=plot,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")
  assign("int_pol",int_pol, envir = .GlobalEnv)
  return(invisible())
}
piece_lin_fix<-function(x,y){
  dat<-data.frame(cbind(x,y)) #create data frame
  af <- approxfun(x,y)
  plot<-ggplot(dat, aes(x=x, y=y)) + ggtitle("I-Curve: Piecewise Linear Interpolation")+
    geom_point(size=5, col='red') + stat_function(fun = af, size=1.25, alpha=0.4)+
    labs(x = "Maturities [years]",y="Yields [%]")
  print(plot)
  ggsave(filename=paste("linear_ic.pdf"), plot=plot,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")
  assign("af",af, envir = .GlobalEnv)
  return(invisible())
}
cub_spline_fix<-function(x,y){
  dat<-data.frame(cbind(x,y)) #create data frame
  spl <- splinefun(x, y, method="fmm")
  plot<-ggplot(dat, aes(x=x, y=y)) + ggtitle("I-Curve: Cubic Spline Interpolation")+
    geom_point(size=5, col='purple') + stat_function(fun = spl, size=1.25, alpha=0.4)+
    labs(x = "Maturities [years]",y="Yields [%]")
  print(plot)
  ggsave(filename=paste("spline_ic.pdf"), plot=plot,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")
  assign("spl",spl, envir = .GlobalEnv)
  return(invisible())
}

lagrange_fix(maturities,yields)
piece_lin_fix(maturities,yields)
cub_spline_fix(maturities,yields)

#TOGETHER PLOTS
plot_tog<-ggplot(dat, aes(x=maturities, y=yields)) + ggtitle("Interpolated Yield Curve")+geom_point(size=3, col='red')+
  labs(x = "Maturities [years]",y="Yields [%]")+
  stat_function(fun = int_pol, size=0.8, aes(colour="Lagr. int."))+
  stat_function(fun = af, size=0.8, aes(colour="Piec. lin."))+
  stat_function(fun = spl, size=0.8, aes(colour="Cub. spl."))+
  scale_colour_manual("", values = c("#CC00CC","blue","red"))
print(plot_tog)

ggsave(filename="together_ic.pdf", plot=plot_tog,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")

#CALCULATE ERRORS ACCORDING TO CUBIC SPLINE FUNCTION
err_lagr_ic<-absfun(diff(spl,int_pol))
err_lin_ic<-absfun(diff(spl,af))

#PLOT ERRORS
plot_err<-ggplot(dat, aes(x=maturities, y=yields)) + ggtitle("Error Functions for I-Curve")+
  labs(x = "Maturities [years]",y="Absolute error function")+
  stat_function(fun = err_lagr_ic, size=0.8, aes(colour="Lagr. int."))+
  stat_function(fun = err_lin_ic, size=0.8, aes(colour="Piec. lin."))+
  scale_colour_manual("", values = c("blue","red"))
print(plot_err)

ggsave(filename="errors_ic.pdf", plot=plot_err,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")




#LEAST SQUARE APPROXIMATIONS POLYNOMIALS
least_sq(maturities,yields,1)
least_1<-lsq
least_sq(maturities,yields,2)
least_2<-lsq
least_sq(maturities,yields,3)
least_3<-lsq

plot_lsq<-ggplot(dat, aes(x=maturities, y=yields)) + ggtitle("Least Square Approximation of I-Curve")+geom_point(size=3, col='green')+
  labs(x = "Maturities [years]",y="Yields [%]")+
  stat_function(fun = least_1, size=0.8, aes(colour="p1(x)"))+
  stat_function(fun = least_2, size=0.8, aes(colour="p2(x)"))+
  stat_function(fun = least_3, size=0.8, aes(colour="p3(x)"))+
  scale_colour_manual("", values = c("blue","#CC00CC","#009999"))
print(plot_lsq)

ggsave(filename="lsq_icu.pdf", plot=plot_lsq,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")


