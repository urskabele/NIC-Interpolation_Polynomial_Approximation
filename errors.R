#CALCULATE ERRORS

#subtracts functions
diff<-function(a,b){
  force(a)
  force(b)
  function(x){a(x)-b(x)}
}

#CALCULATE ERRORS
err_lagr<-diff(f,int_pol)
err_lin<-diff(f,af)
err_spl<-diff(f,spl)
err_tayl<-diff(f,taylor_poly)


