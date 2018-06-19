#HERMITE INTERPOLATION
if(!require(Deriv))install.packages("Deriv");library(Deriv)

f<-sin
df<-Deriv(f)
x<-c(-2:2)
y<-f(x)
dy<-df(x)

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

basis<-lagrange.basis(x,y)


#create Hs
hs<-list()
for (i in 1:length(x)){
  #odv<-odvajas basis[[i]]!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #odv<-odv(x[i])
  r1<-function(t){1-2*(t-x[i])}
  r2<-multiply(basis[[i]],basis[[i]])
  hs[[i]]<-multiply(r1,r2)
}

#create H_hats
hats<-list()
for (i in 1:length(x)){
  s1<-function(t){(t-x[i])}
  s2<-multiply(basis[[i]],basis[[i]])
  hats[[i]]<-multiply(s1,s2)
}

#Sum up!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ni dobr!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
for (i in 1:length(x)){
sumf(multiply(y[i],hs[[i]]),multiply(dy[i],hats[[i]]))
}
