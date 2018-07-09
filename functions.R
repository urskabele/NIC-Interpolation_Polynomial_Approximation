if(!require(ggplot2))install.packages("ggplot2");library(ggplot2)
if(!require(rSymPy))install.packages("rSymPy");library(rSymPy)
if(!require(DeducerExtras))install.packages("DeducerExtras");library(DeducerExtras)
if(!require(polynom))install.packages("polynom");library(polynom)
if(!require(pracma))install.packages("pracma");library(pracma)
if(!require(interp))install.packages("interp");library(interp)

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
























