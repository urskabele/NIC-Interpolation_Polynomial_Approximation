if(!require(ggplot2))install.packages("ggplot2");library(ggplot2)
if(!require(rSymPy))install.packages("rSymPy");library(rSymPy)
if(!require(DeducerExtras))install.packages("DeducerExtras");library(DeducerExtras)
if(!require(polynom))install.packages("polynom");library(polynom)
if(!require(pracma))install.packages("pracma");library(pracma)
if(!require(interp))install.packages("interp");library(interp)

###################################################################################################################################
#ALL DEFINED FUNCTIONS

#multiplies two functions
multiply<-function(a,b){
  force(a)
  force(b)
  function(x){a(x)*b(x)}
}

#sums functions
sumf<-function(a,b){
  force(a)
  force(b)
  function(x){a(x)+b(x)}
}

#subtracts functions
diff<-function(a,b){
  force(a)
  force(b)
  function(x){a(x)-b(x)}
}

#function for vector with alternating sign
alt<-function(len){
  vec<-c()
  for (i in 1:len){
    vec[i]<-(-1)^(i-1)
  }
  return(vec)}

#function for absolute value
absfun<-function(g){
  force(g)
  function(x){abs(g(x))}
}

#RETURNS THE LIST OF LAGRANGE BASIS POLYNOMIALS
lagrange.basis<-function(x,y){
  l<-list() #List to store lagr. basis polynomials
  k<-1
  for (i in x) {
    # Set the numerator and denominator of the Lagrange polynomials to 1 and build them up
    num <- 1
    denom <- 1
    # Remove the current x value from the iterated list
    p <- x[! x %in% i]
    # For the remaining points, construct the Lagrangian polynomial by successively appending each x value
    for (j in p) {
      num <- num*(Var("z")-j)
      denom <-denom*(i-j)
    }
    # Set each Lagrangian polynomial in rSymPy to simplify later.
    l[k] <- num/denom
    k <- k + 1
  }
  args<-"z"
  for (i in 1:length(l)){
    l[[i]]<-eval(parse(text = paste('f <- function(', args, ') {l[[i]]}', sep='')))
  }
  return(l)
}

#LAGRANGE INTERPOLATION POLYNOMIAL
lagrange<-function(x,y){
  n<-length(x)-1
  dat<-data.frame(cbind(x,y)) #create data frame
  int_pol<-poly.calc(x, y) #Lagrange interpolation polynomial
  coefi<-coef(int_pol)
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
    xlim(a-5,a+5)+
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

#LEAST MAX APPROXIMATION
#REMEZ ALGORITHM
#we are looking for the approximating polynomial of degree n
#we are approximating function f
#E is the initial set of n+2 nodes

remez<-function(f,n,E){
  E2<-rep(0,length(E))
  # x<-c(1)
  # y<-c(1)
  # dat<-data.frame(cbind(x,y))
  
  while (all(E==E2)==FALSE){
    E2<-E
    #solve linear system
    fun<-f(E)
    mones<-alt(length(E))
    ones<-rep(1,length(E))
    A<-cbind(mones,ones)
    
    for (i in 1:n){
      A<-cbind(A,E^i)
    }
    
    unknowns<-solve(A,fun)
    coeff<-unknowns[-1]
    
    #derive polynomial
    pol<-polynomial(coeff)
    label<-paste0("p",paste(n),"(x)==",paste(pol))
    pol<-as.function(pol)
    
    #define residual function
    r<-diff(f,pol)
    rabs<-absfun(r)
    
    #find maximum
    maxim<-optimize(rabs, interval=c(E[1], E[length(E)]), maximum=TRUE)$maximum
    
    #find neighbours
    closest<-which.min(abs(E-maxim))
    
    if (maxim>E[closest]){
      second<-closest+1
    } else{
      second<-closest-1
    }
    
    #compute values
    rmaxim<-r(maxim)
    rclosest<-r(E[closest])
    rsecond<-r(E[second])
    
    #change E
    if (sign(rmaxim)==sign(rclosest)){
      E[closest]<-maxim
    } else if (sign(rmaxim)==sign(rsecond)){
      E[second]<-maxim}
    print(E)
  }
  plot(pol) #je ok
  #final polynomial is pol
  #plot
  x<-c(1)
  y<-c(1)
  dat<-data.frame(cbind(x,y))
  #label<-paste0("p",paste(n),"(x)==",paste(pol))
  plot<-ggplot(dat,aes(x=x, y=y)) + ggtitle("Least Maximum Approximation Polynomial")+
    stat_function(fun = f, size=0.8, aes(colour="f(x)"))+
    stat_function(fun = pol, size=0.8, aes(colour="p(x)"))+
    scale_colour_manual("", values = c("red","blue"))+xlim(min(E)-4,max(E)+4)+
    annotate("label",x=-Inf,y=Inf,hjust=0,vjust=2,label = label,parse = TRUE)
  print(plot)
  ggsave(filename=paste("least_max.pdf"), plot=plot,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")
  assign("lmax",pol, envir = .GlobalEnv)
  return(invisible())
}

#################################################################################################################
#YIELDS (PLOT)
mat1<-c(1:11)
yields1<-c(0.1,0.15,0.25,0.5,1,1.5,2.5,3,3.4,4.1,4.2)
dat1<-data.frame(cbind(mat1,yields1))
piece_lin(mat1,yields1)
af1<-af
ticks<-c("1 mo","3 mo", "6 mo", "1 yr", "2 yr", "3 yr", "5 yr", "7 yr", "10 yr", "20 yr", "30 yr")

mat2<-c(1,2,3,3.5,4,5,6,6.5,7,8,9,9.5,10,10.5,11)
yields2<-c(0.17,0.36,0.65,0.95,1.35,1.85,2.8,3.4,3.9,4.2,4.9,5,5.2,5.7,5.9)
dat2<-data.frame(cbind(mat2,yields2))
piece_lin(mat2,yields2)
af2<-af

missing<-c(3.5,6.5,9.5,10.5)
missingy<-af1(missing)
dat_m<-data.frame(cbind(missing,missingy))

yplot<-ggplot(dat2,aes(x=mat2, y=yields2)) + ggtitle("Yield Curve")+
  geom_point(size=3, col='red')+
  stat_function(fun = af1, size=0.8, aes(colour="Government Treasury"))+
  geom_point(data=dat1,aes(x=mat1, y=yields1),size=3, col='blue')+
  geom_point(data=dat_m,aes(x=missing, y=missingy),size=3, col='black',shape=17)+
  ylim(0,6.5)+
  stat_function(fun = af2, size=0.8, aes(colour="Corporate bond"))+
  xlab("Time to maturity")+ylab("Yield [%]")+
  scale_x_continuous(breaks=c(1:11),labels=ticks)+
  geom_segment(aes(x=mat2,y=af1(mat2),xend=mat2,yend=af2(mat2)),linetype="dotted",size=1.3)+
  geom_segment(aes(x=3.5,y=af1(3.5),xend=3.5,yend=af2(3.5)),size=1.3)+
  geom_segment(aes(x=6.5,y=af1(6.5),xend=6.5,yend=af2(6.5)),size=1.3)+
  geom_segment(aes(x=9.5,y=af1(9.5),xend=9.5,yend=af2(9.5)),size=1.3)+
  geom_segment(aes(x=10.5,y=af1(10.5),xend=10.5,yend=af2(10.5)),size=1.3)+
  scale_colour_manual("", values = c("red","blue"))+ theme(legend.position="top")

print(yplot)

ggsave(filename="yields.pdf", plot=yplot,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")


##########################################################################################################################
#PLOT SOME BASIC METHODS

#LAGRANGE INTERPOLATION
x <- c(1:8)
y<- c(1,3,4,2,-2,-3,0,1)
lagrange(x,y)

#PIECEWISE LINEAR EXAMPLE 
x <- c(1:8)
y<- c(1,3,4,2,-2,-3,0,1)
piece_lin(x,y)

#CUBIC SPLINE EXAMPLE
x <- c(1:8)
y<- c(1,3,4,2,-2,-3,0,1)
cub_spline(x,y)

#TAYLOR APPROXIMATION
tayl_app(sin,0,5)

#LEAST MAX OF X^2 APPROXIMATION WITH REMEZ ALGORITHM
E<-c(0,0.25,1)
f<-function(x){x^2}
n<-1
remez(f,n,E)

#LEAST SQUARE APPROXIMATION
x <- c(1:8)
y<- c(1,3,4,2,-2,-3,0,1)
least_sq(x,y,3)


#######################################################################################################################
#SIMULATION 1

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

##########################################################
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

#########################################################
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


########################################################

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

########################################################
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

##########################################################
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

##########################################################################################################################################
#SIMULATION 2
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











