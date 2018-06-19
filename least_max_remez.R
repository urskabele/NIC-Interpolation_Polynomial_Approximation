#LEAST MAX APPROXIMATION
#REMEZ ALGORITHM
#we are looking for the approximating polynomial of degree n
#we are approximating function f
#E is the initial set of n+2 nodes

remez<-function(f,n,E){
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
    pol<-as.function(polynomial(coeff))
    
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
  # plot<-ggplot(dat,aes(x=x, y=y)) + ggtitle("Interpolation and Approximation")+
  #   stat_function(fun = f, size=0.8, aes(colour="f(x)"))+
  #   stat_function(fun = pol, size=0.8, aes(colour="p",paste(n),"(x)"))+
  #   scale_colour_manual("", values = c("red","blue"))
  # print(plot)
  # ggsave(filename=paste("least_max.pdf"), plot=plot,path="C:\\Users\\urska\\Desktop\\Numerical Intr Course\\R codes")
  assign("lmax",pol, envir = .GlobalEnv)
  return(invisible())
  
}

#######
E<-c(0,0.25,1)
f<-function(x){x^2}
n<-1

remez(f,n,E)























