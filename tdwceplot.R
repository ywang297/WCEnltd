tdwceplot<-function(x,m,p,coef,knottime,var_order){
  tdfct<-0
  for (k in 1:(m+p+1)){
    tdfct<-tdfct+coef[k,var_order]*spli(x,k,p,knottime)
  }
  return(tdfct)
}