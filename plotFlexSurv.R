plot.FlexSurv<-function(model.FlexSurv,variable,TD,NL,TimePoint=-999,ref.value.NL=0,...){
  
  
  variableNEW<-match(variable,model.FlexSurv$variables)
  
  tdestim1<-function(x){
    tdfct<-0
    for (k in 1:(model.FlexSurv$Number_knots+model.FlexSurv$Degree_of_splines+1)){
      tdfct<-tdfct+model.FlexSurv$coefficients_splines_TD[k,variableNEW]*spli(x,k,model.FlexSurv$Degree_of_splines,model.FlexSurv$knots_time)
    }
    return(tdfct)
  }
  
  
  nlestim1<-function(y){
    nlfct<-0
    for (k in 1:(model.FlexSurv$Number_knots+model.FlexSurv$Degree_of_splines+1)){
      nlfct<-nlfct+model.FlexSurv$coefficients_splines_NL[k,variableNEW]*spli(y,k,model.FlexSurv$Degree_of_splines,model.FlexSurv$knots_covariates[variableNEW,])
    }
    return(nlfct)
  }
  
  if(TD==1){
    axist<-seq(0,model.FlexSurv$knots_time[model.FlexSurv$Degree_of_splines+2+model.FlexSurv$Number_knots],0.1)
  }
  if(NL==1){
    axisx<-seq(model.FlexSurv$knots_covariates[variableNEW,1],model.FlexSurv$knots_covariates[variableNEW,model.FlexSurv$Degree_of_splines+2+model.FlexSurv$Number_knots],0.1)
  }
  
  if (TD==1 & NL==0){
    plot(axist,tdestim1(axist),...)
  }
  
  
  if (TD==0 & NL==1){
    plot(axisx,nlestim1(axisx)-nlestim1(ref.value.NL),...)
  }
  
  
  if (TD==1 & NL==1 & TimePoint==-999){
    op<-par(mfrow=c(2,1))
    plot(axist,tdestim1(axist),...)
    plot(axisx,nlestim1(axisx)-nlestim1(ref.value.NL),...)
    par(op)
  }
  
  if (TD==1 & NL==1 & TimePoint!=-999){
    plot(axisx,(nlestim1(axisx)-nlestim1(ref.value.NL))*tdestim1(TimePoint),...)
  }
}
