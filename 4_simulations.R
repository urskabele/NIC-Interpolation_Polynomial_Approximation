#SIMULATIONS

################################################################################################################
#TRY WITH DATA, 7 NODES
f<-function(x){sin(2*x)}
x<-c(-5,-3,-1,0.1,1,3,5)
y<-f(x)
dat<-data.frame(cbind(x,y))

#CALCULATE ALL APPROXIMATIONS
lagrange(x,y)
piece_lin(x,y)
cub_spline(x,y)
tayl_app(f,0,5)
least_sq(x,y,5)

#PLOT
plot_tog<-ggplot(dat, aes(x=x, y=y)) + ggtitle("Interpolation and Approximation")+ylim(-7,7)+xlim(-7,7)+
  stat_function(fun = int_pol, size=0.8, aes(colour="Lagr. int."))+
  stat_function(fun = af, size=0.8, aes(colour="Piec. lin."))+
  stat_function(fun = spl, size=0.8, aes(colour="Cub. spl."))+
  stat_function(fun = lsq, size=0.8, aes(colour="Least Sq."))+
  stat_function(fun = f, size=0.8, aes(colour="f(x)"))+
  geom_point(size=3, col='red')+
  scale_colour_manual("", values = c("blue","red","#CC00CC","#009999","#00CC00"))
print(plot_tog)

ggsave(filename="together.pdf", plot=plot_tog,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")

#CALCULATE ERRORS
err_lagr<-absfun(diff(f,int_pol))
err_lin<-absfun(diff(f,af))
err_spl<-absfun(diff(f,spl))
err_lsq<-absfun(diff(f,lsq))

#PLOT ERRORS
plot_err<-ggplot(dat, aes(x=x, y=y)) + ggtitle("Error Functions")+ylim(-0,3)+
  stat_function(fun = err_lagr, size=0.8, aes(colour="Lagr. int."))+
  stat_function(fun = err_lin, size=0.8, aes(colour="Piec. lin."))+
  stat_function(fun = err_spl, size=0.8, aes(colour="Cub. spl."))+
  stat_function(fun = err_lsq, size=0.8, aes(colour="Least Sq."))+
  scale_colour_manual("", values = c("blue","#CC00CC","#009999","#00CC00"))
print(plot_err)

ggsave(filename="errors.pdf", plot=plot_err,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")

#######################################################################################################3
#TRY WITH DATA, 11 NODES
f<-function(x){sin(2*x)}
x<-c(-5:5)
y<-f(x)
dat<-data.frame(cbind(x,y))

#CALCULATE ALL APPROXIMATIONS
lagrange(x,y)
piece_lin(x,y)
cub_spline(x,y)
tayl_app(f,0,5)
least_sq(x,y,5)

#PLOT
plot_tog<-ggplot(dat, aes(x=x, y=y)) + ggtitle("Interpolation and Approximation")+geom_point(size=3, col='red')+ylim(-7,7)+xlim(-7,7)+
  stat_function(fun = int_pol, size=0.8, aes(colour="Lagr. int."))+
  stat_function(fun = af, size=0.8, aes(colour="Piec. lin."))+
  stat_function(fun = spl, size=0.8, aes(colour="Cub. spl."))+
  stat_function(fun = lsq, size=0.8, aes(colour="Least Sq."))+
  stat_function(fun = f, size=0.8, aes(colour="f(x)"))+
  scale_colour_manual("", values = c("blue","red","#CC00CC","#009999","#00CC00"))
print(plot_tog)

ggsave(filename="together2.pdf", plot=plot_tog,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")

#CALCULATE ERRORS
err_lagr2<-absfun(diff(f,int_pol))
err_lin2<-absfun(diff(f,af))
err_spl2<-absfun(diff(f,spl))
err_lsq2<-absfun(diff(f,lsq))

#PLOT ERRORS
plot_err<-ggplot(dat, aes(x=x, y=y)) + ggtitle("Error Functions")+ylim(-0,3)+
  stat_function(fun = err_lagr2, size=0.8, aes(colour="Lagr. int."))+
  stat_function(fun = err_lin2, size=0.8, aes(colour="Piec. lin."))+
  stat_function(fun = err_spl2, size=0.8, aes(colour="Cub. spl."))+
  stat_function(fun = err_lsq2, size=0.8, aes(colour="Least Sq."))+
  scale_colour_manual("", values = c("blue","#CC00CC","#009999","#00CC00"))
print(plot_err)

ggsave(filename="errors2.pdf", plot=plot_err,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")

#######################################################################################################3
#TRY WITH DATA, RANDOM 11 NODES
f<-function(x){sin(2*x)}
x<-sort(runif(11, -5,5))
y<-f(x)
dat<-data.frame(cbind(x,y))

#CALCULATE ALL APPROXIMATIONS
lagrange(x,y)
piece_lin(x,y)
cub_spline(x,y)
tayl_app(f,0,5)
least_sq(x,y,5)

#PLOT
plot_tog<-ggplot(dat, aes(x=x, y=y)) + ggtitle("Interpolation and Approximation")+geom_point(size=3, col='red')+ylim(-7,7)+xlim(-7,7)+
  stat_function(fun = int_pol, size=0.8, aes(colour="Lagr. int."))+
  stat_function(fun = af, size=0.8, aes(colour="Piec. lin."))+
  stat_function(fun = spl, size=0.8, aes(colour="Cub. spl."))+
  stat_function(fun = lsq, size=0.8, aes(colour="Least Sq."))+
  stat_function(fun = f, size=0.8, aes(colour="f(x)"))+
  scale_colour_manual("", values = c("blue","red","#CC00CC","#009999","#00CC00"))
print(plot_tog)

ggsave(filename="together_r.pdf", plot=plot_tog,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")

#CALCULATE ERRORS
err_lagr_r<-absfun(diff(f,int_pol))
err_lin_r<-absfun(diff(f,af))
err_spl_r<-absfun(diff(f,spl))
err_lsq_r<-absfun(diff(f,lsq))

#PLOT ERRORS
plot_err<-ggplot(dat, aes(x=x, y=y)) + ggtitle("Error Functions")+ylim(-0,3)+
  stat_function(fun = err_lagr_r, size=0.8, aes(colour="Lagr. int."))+
  stat_function(fun = err_lin_r, size=0.8, aes(colour="Piec. lin."))+
  stat_function(fun = err_spl_r, size=0.8, aes(colour="Cub. spl."))+
  stat_function(fun = err_lsq_r, size=0.8, aes(colour="Least Sq."))+
  scale_colour_manual("", values = c("blue","#CC00CC","#009999","#00CC00"))
print(plot_err)

ggsave(filename="errors_r.pdf", plot=plot_err,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")


#############################################################################################################

#PLOTS OF ERRORS WITH DIFFERENT NUMBER OF NODES
#LAGRANGE
plot_err_lag<-ggplot(dat, aes(x=x, y=y)) + ggtitle("Error Functions for Lagrange Interpolation")+ylim(0,3)+
  stat_function(fun = err_lagr, size=0.8, aes(colour="7 nodes"))+
  stat_function(fun = err_lagr2, size=0.8, aes(colour="11 nodes"))+
  stat_function(fun = err_lagr_r, size=0.8, aes(colour="11 random nodes"))+
  scale_colour_manual("", values = c("green","red","blue"))
print(plot_err_lag)

ggsave(filename="errors_lag.pdf", plot=plot_err_lag,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")

#PIECEWISE LINEAR
plot_err_lin<-ggplot(dat, aes(x=x, y=y)) + ggtitle("Error Functions for Piecewise Linear Interpolation")+ylim(0,3)+
  stat_function(fun = err_lin, size=0.8, aes(colour="7 nodes"))+
  stat_function(fun = err_lin2, size=0.8, aes(colour="11 nodes"))+
  stat_function(fun = err_lin_r, size=0.8, aes(colour="11 random nodes"))+
  scale_colour_manual("", values = c("green","red","blue"))
print(plot_err_lin)

ggsave(filename="errors_lin.pdf", plot=plot_err_lin,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")

#CUBIC SPLINE
plot_err_spl<-ggplot(dat, aes(x=x, y=y)) + ggtitle("Error Functions for Cubic Spline Interpolation")+ylim(0,3)+
  stat_function(fun = err_spl, size=0.8, aes(colour="7 nodes"))+
  stat_function(fun = err_spl2, size=0.8, aes(colour="11 nodes"))+
  stat_function(fun = err_spl_r, size=0.8, aes(colour="11 random nodes"))+
  scale_colour_manual("", values = c("green","red","blue"))
print(plot_err_spl)

ggsave(filename="errors_spl.pdf", plot=plot_err_spl,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")


#LEAST SQUARE
plot_err_lsq<-ggplot(dat, aes(x=x, y=y)) + ggtitle("Error Functions for Least Square Approximation")+ylim(0,3)+
  stat_function(fun = err_lsq, size=0.8, aes(colour="7 nodes"))+
  stat_function(fun = err_lsq2, size=0.8, aes(colour="11 nodes"))+
  stat_function(fun = err_lsq_r, size=0.8, aes(colour="11 random nodes"))+
  scale_colour_manual("", values = c("green","red","blue"))
print(plot_err_lsq)

ggsave(filename="errors_lsq.pdf", plot=plot_err_lsq,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")

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

#CALCULATE AND PLOT ERRORS
err_tay3<-diff(f,tay_poly_3)
err_tay5<-diff(f,tay_poly_5)
err_tay7<-diff(f,tay_poly_7)

plot_err_tay<-ggplot(dat, aes(x=x, y=y)) + ggtitle("Error Functions for Taylor Polynomials")+
  stat_function(fun = err_tay3, size=0.8, aes(colour="T3"))+
  stat_function(fun = err_tay5, size=0.8, aes(colour="T5"))+
  stat_function(fun = err_tay7, size=0.8, aes(colour="T7"))+
  xlim(-3,3)+
  scale_colour_manual("", values = c("blue","#CC00CC","#009999"))
print(plot_err_tay)

ggsave(filename="errors_tay.pdf", plot=plot_err_tay,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")

##############################################################################################################
#LEAST SQUARES DIFFERENT
x<-sort(runif(30,0,10))
y<-sort(runif(30,-30,50))
dat<-data.frame(cbind(x,y))

least_sq(x,y,1)
least_1<-lsq
least_sq(x,y,2)
least_2<-lsq
least_sq(x,y,3)
least_3<-lsq

plot_lsq<-ggplot(dat, aes(x=x, y=y)) + ggtitle("Least Square Approximation")+geom_point(size=3, col='green')+
  stat_function(fun = least_1, size=0.8, aes(colour="p1(x)"))+
  stat_function(fun = least_2, size=0.8, aes(colour="p2(x)"))+
  stat_function(fun = least_3, size=0.8, aes(colour="p3(x)"))+
  #stat_function(fun = f, size=0.8, aes(colour="f(x)"))+
  scale_colour_manual("", values = c("blue","#CC00CC","#009999"))
print(plot_lsq)

ggsave(filename="lsq_diff.pdf", plot=plot_lsq,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")
