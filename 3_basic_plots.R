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