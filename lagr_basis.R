#returns the list of lagrange basis polynomials
if(!require(rSymPy))install.packages("rSymPy");library(rSymPy)

lagrange.basis<-function(x,y){
  l<-list() #List to store lagr. basis polynomials
  k<-1
  for (i in x) {
    # Set the numerator and denominator of the Lagrangian polynomials to 1 and build them up
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
    l[[i]]<-eval(parse(text = paste('f <- function(', args, ') { return(' , l[[i]] , ')}', sep='')))
  }
  return(l)
}

x <- c(1, 0, 3, 4)
y <- c(7, 11, 28, 63)
basis<-lagrange.basis(x,y)
basis[[1]](7)

