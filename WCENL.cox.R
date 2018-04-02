library(plyr)
library(survival)
library(splines)

load("/Users/yishu/Dropbox/TEL/FlexSurv_30MAR2017.RData")


## run last_prog, last_prog3, .knots.equi, .augm.knots, DvlpMatrix, .my.bic 

spli <- function(x, j, p, knots) { # regression B-spline base
  if (p == 0) {
    b <- ifelse(x >= knots[j] & x < knots[j + 1], 1, 
                0)
    return(b)
  }
  else {
    a1 <- ifelse(rep(knots[j] != knots[j + p], length(x)), 
                 (x - knots[j])/(knots[j + p] - knots[j]), 0)
    a2 <- ifelse(rep(knots[j + p + 1] != knots[j + 1], 
                     length(x)), (knots[j + p + 1] - x)/(knots[j + 
                                                                 p + 1] - knots[j + 1]), 0)
    return(a1 * spli(x, j, p - 1, knots) + a2 * spli(x, 
                                                     j + 1, p - 1, knots))
  }
}

#add p compared to the original function from WCE package
.augm.knots <- function(inner, f.up, p){ # augments the set of interior knots for spline basis
  ret <- c(seq(-p,0,1), inner, seq(f.up,f.up+p,1))
  names(ret) <- NULL
  ret
}

.knots.equi <- function(n.knots, m){ # use quantiles for knots placement
  if (n.knots==1){
    f <- round(quantile(seq(1,m),seq(0,1, by=1/(n.knots+1))),0)[2]} else {
      f <- round(quantile(seq(1,m),seq(0,1, by=1/(n.knots+1))),0)[-1]}
  return(f[1:(length(f)-1)])}  

.my.bic <- function(PL, n.events, npara , aic=FALSE){ # estimate AIC OR BIC for the model 
      

                    if (aic == TRUE) 	{	
                      bic <- -2*PL + npara * 2} else {
                        bic <- -2*PL + npara * log(n.events)}
  return(bic)}

.wcecalc <- function(ev, dose, stop, Bbasis, cutoff){#  calculated the D_j(u) for a given individual for all relevant risk sets
  fup <- length(dose) ##length of fup for each individual 
  myev <- ev[ev<=stop[fup]]  ## stop[fup]=last fup of this individual, myev stores all the event time before the end of fup of this patient
  if (length(myev)>0)    { ## if the patients are in risk sets before the end of his follow up 
    linesfor1 <- matrix(NA, ncol=dim(Bbasis)[2], nrow=length(myev)) ## store the calculated D_j(u) of this patients at each event time               
    for (i in 1:length(myev)){
      vec <- dose[stop <= myev[i]] 
      pos <- length(vec)
      if (pos<cutoff) { # if the previous time of the patient at some event time less then the predefined time window 
        vec <- c(rep(0, (cutoff-length(vec))), vec)} else {
          pos <- length(vec)
          vec <- vec[(pos-cutoff+1):pos]} 
      linesfor1[i,] <- rev(vec)%*%Bbasis}
  }   else {linesfor1 <- rep(NA, dim(Bbasis)[2])} # if the patients never in risk set, then his D_j(u) is NA
  linesfor1}

DvlpMatrix <- function(data, listeT, ncol, TypeNEW) {
  data <- matrix(data, ncol = ncol) # data is a list of full data on each individual, first transform this list to the matrix 
  if (max(data[, TypeNEW[2]]) < min(listeT[listeT != 0])) {  
    XX <- data  ## if max(stop time) is less then min(event time), then data remain the same
  }
  else { ## if max(stop time of certain individual) > min(event time)
    aindex <- rep(0, (sum(listeT <= max(data[, TypeNEW[2]])) - 
                        1))
    # apb <- rep(0, (sum(listeT <= max(data[, TypeNEW[2]])) - 
    #                  1))
    for (i in 1:(sum(listeT <= max(data[, TypeNEW[2]])) -  ## for the length of unique event time 
                 1)) {
      for (j in 1:(dim(data)[1])) { #for each row of certain individual, if start time < ith order of event time<=stop time of row j, stores
        # the largest row satisfy this condition in aindex[i]
        if (as.numeric(data[j, TypeNEW[1]]) < as.numeric(listeT[1 + i]) & as.numeric(data[j, TypeNEW[2]]) >= as.numeric(listeT[1 + i]))
          aindex[i] <- j
      }
    }   ## xx first coloum repeat ID, coloum for start day assigns c(0,unique event time-last one), colum for stop day assigns unique event time 
    XX <- matrix(nrow = sum(listeT <= max(data[, TypeNEW[2]])) - 
                   1, ncol = ncol)
    XX[1:(sum(listeT <= max(data[, TypeNEW[2]])) - 1), 1] <- rep(data[1, 
                                                                      1], sum(listeT <= max(data[, TypeNEW[2]])) - 1)
    XX[1:(sum(listeT <= max(data[, TypeNEW[2]])) - 1), TypeNEW[1]] <- listeT[1:(sum(listeT <= 
                                                                                      max(data[, TypeNEW[2]])) - 1)]
    XX[1:(sum(listeT <= max(data[, TypeNEW[2]])) - 1), TypeNEW[2]] <- listeT[2:(sum(listeT <= 
                                                                                      max(data[, TypeNEW[2]])))]
    XX[1:(sum(listeT <= max(data[, TypeNEW[2]])) - 1), TypeNEW[3]] <- c(rep(0,(sum(listeT <= max(data[, TypeNEW[2]])) - 2)), data[dim(data)[1],TypeNEW[3]])
    XX[1:(sum(listeT <= max(data[, TypeNEW[2]])) - 1), -c(1, 
                                                          TypeNEW[1], TypeNEW[2], TypeNEW[3])] <- as.matrix(data[aindex, 
                                                                                                                 -c(1, TypeNEW[1], TypeNEW[2], TypeNEW[3])])
    
  }
  X <- XX
  list(X)
}

###### The program only work for unconstrained estiamte of WCE for now 

flexplot<-function(x, pTDNL, nknots,coef,knot){
  tdfct<-0
  for (k in 1:(nknots+pTDNL+1)){  
    tdfct<-tdfct+coef[k]*spli(x,k,pTDNL,knot)
  }
  return(tdfct)
}


# only work for default knots and only allow estimate of WCE with 1 choice of knot for now


# WCE <- function(data, analysis, nknots, cutoff, constrained = FALSE, int.knots = NULL, aic = FALSE, MatchedSet = NULL, id='Id', event = 'Event',  start='Start', stop='Stop', expos ='dose', covariates = NULL, controls = NULL,...) UseMethod("WCE")
# 
# WCE.default <- function(data, analysis, nknots, cutoff, constrained = FALSE, int.knots = NULL, aic = FALSE, id='Id', event = 'Event',  start='Start', stop='Stop', expos ='dose', covariates = NULL, MatchedSet = NULL, controls = NULL, ...) {
#   if (constrained == 'right') constrained = 'Right'
#   if (constrained == 'left') constrained = 'Left'
#   if (is.data.frame(data) == F)  stop("ERROR: data must be a data frame")
#   if (is.null(covariates) == F & sum(covariates %in% names(data))!= length(covariates)) stop("ERROR: At least one covariate does not belong to the data set supplied") 
#   if (analysis == 'Cox' | analysis == 'cox') {WCE.cox(data, nknots, cutoff, constrained, int.knots = NULL, aic, id, event, start, stop, expos, covariates, controls)} else 
#     if (analysis == 'NCC' | analysis == 'ncc') {stop("Methods for nested case control designs are not implemented yet.")} else 
#       if  (analysis == 'CC' | analysis == 'cc') {stop("Methods for case control designs are not implemented yet")} else stop("ERROR: Requested analysis is invalid")
# }






WCEnltd.cox<-function (data, Type, variables, WCE, TD, NL, nknotsWCE, pWCE, cutoff, nknotsTDNL,pTDNL,constrained, 
          aic = FALSE){ 
  V <- length(variables) # number of total covariates
  
  variablesNEW<-match(variables,names(data))  # stores the postition of each variable in the data set 
  TypeNEW<-match(Type,names(data)) # stored the position of the Type varialbes in the data set
  
  if(sum(WCE)==0){ # if no WCE 
    
   
    last_prog3(data, Type, variables,WCE, TD,NL, nknotsTDNL,pTDNL, knots = -999)
        #else if(length(int.knots) != length(nknotsWCE)+2) {## if user defined knots, then it must be a list of length equal to number of different knots considered for WCE and TD NL
      ## Here we assume all the variables that have WCE effects will only be considered with the same knots at one time, ie. can not estimate two WCEs where one has one knot and the other
      ### have two knots at the same time; also for estimate of TD/NL effects, only one choice of knots will be considered at a time, 
      # }
    
    
  }else{ # if any variable have WCE effects 
    
    i1 <- sum((NL+TD+WCE) == 0)    ##number of non-nl non-td non-WCE variables
    i2 <- sum(((NL == 1) & (TD == 0) & (WCE==0))) ##number of NL non-td non-WCE variables
    i3 <- sum(((NL == 0) & (TD == 1) & (WCE==0)))  ##number of non-NL td non-WCE variables
    i4 <- sum((NL + TD) == 2 & (WCE==0)) ###number of NL TD non-WCE variables
    i5 <- sum(((NL == 0) & (TD == 0) & (WCE==1))) ##number of non-NL non-td WCE variables
    i6 <- sum(((NL == 1) & (TD == 0) & (WCE==1))) ##number of NL non-td WCE variables
    i7 <- sum(((NL == 0) & (TD == 1) & (WCE==1))) ##number of non-NL td WCE variables
    i8 <- sum(((NL == 1) & (TD == 1) & (WCE==1))) ##number of NL td WCE variables
    
    nNLb<-nknotsTDNL+pTDNL #number of NL basis
    nWCEb<-nknotsWCE+pWCE+1 #number of WCE basis
    
    if(length(cutoff)!=sum(WCE))
      stop("ERROR:length of cutoff must equal to the number of specified WCE effects")
    if(length(constrained)!=sum(WCE))
      stop("ERROR:length of constrained must equal to the number of specified WCE effects")
    maxTime <- max(data[,TypeNEW[2]])# maximum FUP time 
    if (max(cutoff) > maxTime) ## cutoff can be a vector with length equal to sum(WCE) 
      stop("ERROR: cutoff must be smaller than the longest follow-up time")
    
    if ( sum(constrained %in% c(FALSE, "Right","R", "right", "Left", "L", "left"))!=length(constrained)) 
      stop("ERROR: constrained has to be one of : FALSE, 'Right', or 'Left'.")
    if ( sum(constrained %in% c("R", "right"))>0 ) {
      constrained[constrained %in% c("R", "right")] = rep("Right", sum(constrained %in% c("R", "right")))
    }
    if ( sum(constrained %in% c("L", "left"))>0 ) {
      constrained[constrained %in% c("L", "left")] = rep("Left",  sum(constrained %in% c("L", "left")))
    }
     #### do not allow user specified knots for now, only allow default knots ###
    
    listeprobaquantile <- seq(1, nknotsTDNL)/(nknotsTDNL + 1) ## place the interior knots according to quantiles  
    knotsNEW <- list() ## store both the interior and exterior knots for each model 
    
   ## first, store the knots for WCE estimate 
    for(i in 1:sum(WCE)){ #different variable for WCE may have different cutoff windows thus result in different knots
      knotsNEW[[paste("WCE",nknotsWCE, "knot(s) for",variables[WCE==1][i], sep = " ")]] <- .augm.knots(.knots.equi(nknotsWCE,  
                                                                       cutoff[i]), cutoff[i], pWCE)
    }
   ###2nd, store the knots for all the variables for possible NL estimate, the first coefficient for 
      for (i in 1:V) {
        knotsNEW[[paste(variables[i],nknotsTDNL, "knot(s)", sep = " ")]] <- c(rep(min(data[, variables[i]]), pTDNL + 1), ## place the first p+1 exterior knots at the min(variables[i])
                           quantile(data[, variables[i]], probs = listeprobaquantile), ## place the interior knots at quantile 
                           seq(max(data[, variables[i]]), max(data[, variables[i]])+pTDNL, 1)) ### place the last p+1 exterior knots equally spaced between (max(variables[i]), ...+p)
      }
   ##3rd, store knots for time, TD estimate 
 
     knotsNEW[[paste("TD",nknotsTDNL, "knot(s)", sep = " ")]] <- c(rep(0, pTDNL + 1), ## place the exterior knots for time at 0
                               quantile(data[data[, TypeNEW[3]] == 1, TypeNEW[2]], probs = listeprobaquantile), ### place the interior knots at quantile of event times 
                               seq(max(data[data[,TypeNEW[3]] == 1, TypeNEW[2]]), max(data[data[,TypeNEW[3]] == 1, TypeNEW[2]]) + pTDNL, 1)) 
        # place exterior knots equally spaced between max(event time) and max(event time)+p 
   
     
    
    data <- as.matrix(data)
    listeT <- c(0, sort(unique(data[data[, TypeNEW[3]] == 1, TypeNEW[2]]))) # sort the unique event time and stored in listeT
    ncol <- dim(data)[2] # number of coloum 
    
    
    
    ### remove values not at risk will cause problem in the WCE estimate ????!!!! ####
    #### the calculation of WCE is based on the full data then remove the rows from the risk set
    ##so get rid of values not belonging to risk set is OK here ####
    
    X <- split(data, data[, 1]) ##lists store full data for each individual respectively
    matX <- sapply(X, DvlpMatrix, listeT = listeT, ncol = ncol, TypeNEW=TypeNEW) ## transform the data structure
    QWR <- do.call(rbind, matX) #combine the lists to a full data matrix, get rid of observation not in any risk set
    
  ##Always estimate WCE first, create bspline base for each specified WCE effects###
  ## bspline basis are different when the variable only have WCE and have both WCE and TD
  ### the data still have to be in 1 line per day format since to estiamte the TD effect, the data need to be transformed to this format first
    
    nWCE<-sum(WCE)# number of total variables having WCE effect, use to index the knotsNEW since each variable were generated a knot
    Id <- unique(data[,1])
     #i5=number of varialbes only have WCE effect, i7=no. of varialbes have both WCE+TD
    Bbasis<-list() # (i5+i7)*(nknotsWCE+pWCE+1) spline basis 
    if(i5!=0){
      pos1<-match(variablesNEW[NL == 0 & TD == 0 & WCE==1], variablesNEW[WCE == 1])
      #pos1:mark position of variable only with WCE among all varialbes having WCE
      for(i in 1:i5){
      Bbasis[[i]]<-splineDesign(knots = knotsNEW[[pos1[i]]], x = 1:cutoff[pos1[i]], 
                                ord =pWCE+1)
      kal <- do.call("rbind", lapply(1:length(Id),  # combine the calcluated D_j(u) of each paitents at risk set together
                                     function(j) .wcecalc(listeT[-1], data[data[,1]==Id[j], variablesNEW[NL == 0 &TD==0 & WCE==1][i]], 
                                                          data[data[,1]==Id[j],TypeNEW[2]], Bbasis[[i]], cutoff[pos1[i]])))
      kal <- kal[is.na(kal[, 1]) == FALSE, ] # get rid of all the patients that do not belong to any risk set 
      QWR<-cbind(QWR, kal)
      }
    }
    if(i7!=0){
     pos2<-match(variablesNEW[NL == 0 & TD == 1 & WCE==1], variablesNEW[WCE == 1])
     #pos2:mark position of variable with WCE+TD among all varialbes having WCE
     for(i in 1:i7){
       Bbasis[[(i5+i)]]<-splineDesign(knots = knotsNEW[[pos2[i]]], x = 1:cutoff[pos2[i]], 
                                 ord =pWCE+1)
       kal <- do.call("rbind", lapply(1:length(Id),  # combine the calcluated D_j(u) of each paitents at risk set together
                                      function(j) .wcecalc(listeT[-1], data[data[,1]==Id[j], variablesNEW[NL == 0 & TD==1& WCE==1][i]], 
                                                           data[data[,1]==Id[j],TypeNEW[2]], Bbasis[[(i5+i)]], cutoff[pos2[i]])))
       kal <- kal[is.na(kal[, 1]) == FALSE, ] # get rid of all the patients that do not belong to any risk set 
       QWR<-cbind(QWR, kal)
     }
    }

    # generate the spline values for NL effect not having WCE, since the cumulative sum of the spline basis are not needed, the spline
    ## base can be generated based on QWR where the observation not from risk set are removed 
   if(i2!=0){
     for (i in 1:i2){
       QWR <- cbind(QWR, splineDesign(knotsNEW[[seq(nWCE+1, nWCE+V, 1)[NL == 1&WCE==0&TD==0][i]]], x=QWR[, variablesNEW[NL == 1&WCE==0&TD==0][i]],ord=pTDNL+1)[,-1] )    
     }
   }
    if(i4!=0){
      for (i in 1:i4){
        QWR <- cbind(QWR, splineDesign(knotsNEW[[seq(nWCE+1, nWCE+V, 1)[NL == 1&WCE==0&TD==1][i]]], x=QWR[, variablesNEW[NL == 1&WCE==0&TD==1][i]],ord=pTDNL+1)[,-1] )    
      }
    }

  # spline value for time, pTDNL+nknotTDNL+1 
    if((i3+i4+i7+i8)!=0){
    QWR <- cbind(QWR, splineDesign(knotsNEW[[nWCE+V+1]], x=QWR[, TypeNEW[2]],ord=pTDNL+1))    
 # check the basis for time again, it is generated after the observation not in the risk set are removed?
    }
    
    #After this step, the order of variables in QWR: variables, (i5+i7)*(nknotsWCE+pWCE+1)+(i2+i4)*(nNLb)+nNLb+1
    
    
   # for the variables having both WCE and NL effect, need to generate NL basis based on the original data instead of the reduced data
    Dbasis<-list()
    Nbasis<-list() #include NL basis for NL+WCE and NL+WCE+TD
    if(i6!=0){
      pos8<-match(variablesNEW[NL == 1 & TD == 0 & WCE==1], variablesNEW[WCE == 1])
      #pos8: positions of variables having WCE+NL effects among those having WCE to locate the position for cutoff
      for(i in 1:i6){
        Dbasis[[i]]<-splineDesign(knots = knotsNEW[[pos8[i]]], x = 1:cutoff[pos8[i]], ord =pWCE+1)
        Nbasis[[i]]<-splineDesign(knotsNEW[[seq(nWCE+1, nWCE+V, 1)[NL == 1&WCE==1&TD==0][i]]], x=data[, variablesNEW[NL == 1&WCE==1&TD==0][i]],ord=pTDNL+1)[,-1]
      }
    }
    if(i8!=0){
      pos5<-match(variablesNEW[NL == 1 & TD == 1 & WCE==1], variablesNEW[WCE == 1])
      #pos5: positions of variables having 3 effects among those having WCE to locate the position for cutoff
      for(i in 1:i8){
        Dbasis[[(i6+i)]]<-splineDesign(knots = knotsNEW[[pos5[i]]], x = 1:cutoff[pos5[i]], ord =pWCE+1)
        Nbasis[[(i6+i)]]<-splineDesign(knotsNEW[[seq(nWCE+1, nWCE+V, 1)[NL == 1&WCE==1&TD==1][i]]], x=data[, variablesNEW[NL == 1&WCE==1&TD==1][i]],ord=pTDNL+1)[,-1]
      }
    }  
    
    
    
    tt <- paste("modX<-coxph(Surv(QWR[,", TypeNEW[1], "],QWR[,", 
                TypeNEW[2], "],QWR[,", TypeNEW[3], "])~", sep = "")
    
    
    
    
    ### res stores the final rescale coefficient####### 
    res<-rep(NA,V)
    
    
   ############## if there is none or only 1 special effect #####
    Nn <- 0 
    if (i1 != 0) { # if there are variables have neither NL nor TD nor WCE effects 
      for (k in 1:i1) {
        tt <- paste(tt, "QWR[,", variablesNEW[NL == 0 & TD == 
                                                0 & WCE==0][k], "]+", sep = "")
        Nn <- Nn + 1
      }
    }##Nn=number (NL=TD=WCE=0)*i1
    
    indW<-c(0,rep(NA,i5))
    if (i5 != 0) { #only WCE
      pos3<-match(variablesNEW[NL == 0 & TD == 0 & WCE==1], variablesNEW[WCE == 1&NL==0])
    ## pos3: find the position of varialbe have only WCE among variables have WCE and WCE+TD
      #pos1:mark position of variable only with WCE among all varialbes having WCE
        for (k in 1:i5) {
          if(constrained[pos1[k]]=="FALSE" ){
            indW[k+1]<-nWCEb
            covp<-paste( "QWR[,", (dim(data)[2] + (pos3[k] -1) * nWCEb + 1):(dim(data)[2] + pos3[k] * nWCEb)  , "]", sep = "")
            tt<-paste(tt, paste(c(covp,""), collapse= "+") )
            Nn<-Nn+nWCEb
          }else if(constrained[pos1[k]]=="Right"){
            indW[k+1]<-nWCEb-2
            covp<-paste( "QWR[,", (dim(data)[2] + (pos3[k] -1) * nWCEb + 1):(dim(data)[2] + pos3[k] * nWCEb-2)  , "]", sep = "")
            tt<-paste(tt, paste(c(covp,""), collapse= "+") )
            Nn<-Nn+nWCEb-2
          }else{
            indW[k+1]<-nWCEb-2
            covp<-paste( "QWR[,", (dim(data)[2] + (pos3[k] -1) * nWCEb + 3):(dim(data)[2] + pos3[k] * nWCEb)  , "]", sep = "")
            tt<-paste(tt, paste(c(covp,""), collapse= "+") )
            Nn<-Nn+nWCEb-2
          }
        }
      
      }          
    cumindW<-cumsum(indW)  

    ##Nn=number (NL=TD=WCE=0)*i1 +(TD=NL=0,WCE=1)*i5
    
    
    if (i2 != 0) { #only NL
      for (k in 1:i2) { # all the WCE bases are created, no need to use cumind to index
        covp<-paste( "QWR[,", (dim(data)[2] + (i5+i7)*nWCEb + (k -1) * nNLb + 1):(dim(data)[2] + (i5+i7)*nWCEb + k * nNLb)  , "]", sep = "")
        tt<-paste(tt, paste(c(covp,""), collapse= "+") )
        Nn<-Nn+nNLb
      }
    }##Nn=number (NL=TD=WCE=0)*i1 +(TD=NL=0,WCE=1)*i5+(TD=WCE=0, NL=1)*i2
    
    if (i3 != 0) {# only TD effect
      for (k in 1:i3) {
        flag<-dim(QWR)[2]+1
        QWR <- cbind(QWR, QWR[, variablesNEW[NL == 0 & TD ==1& WCE==0][k]] * QWR[, (dim(data)[2] + (i5+i7)*nWCEb+(i2+i4)*nNLb 
                                      +1): (dim(data)[2] + (i5+i7)*nWCEb+(i2+i4+1)*nNLb+1)]) 
        covp<-paste("QWR[,", flag: dim(QWR)[2], "]", sep = "" )
        tt<-paste(tt, paste(c(covp,""), collapse= "+"))  #### only the TD spline base combined with covariate values 
        # will enter the model 
        Nn <- Nn + nNLb+1
      }##Nn=number (NL=TD=WCE=0)*i1 +(TD=NL=0,WCE=1)*i5+(TD=WCE=0, NL=1)*i2 +(NL=WCE=0,TD=1)*i3
    } #Nn mark the number of parameters having no or only 1 special effect
    tt2 <- tt ## tt mark the non or only 1 effect estimates
    
    
    ### having more than one special effects###
    
    vrais <- c()
    
    ######### previously tt take care of all the variables having only 1 or no special effect, in this case, no ACE if needed, so 
    # they can be directly estimated, now starting to considere the case when ACE is needed, 2 different situations, if there are 
    #variables having 3 effects, then 3 fold ACE is needed, if only exists variables having 2 effects, then only 2 fold of ACE is needed
    # so look at these two situations seperately if(3 effects) else(any of the 2 effects) 
   ## starting with estimating WCE effect first  #####  
    indWNT<-c(0,rep(NA,i8))
    indWN<-c(0,rep(NA,i6))
    indWT<-c(0,rep(NA,i7))
    cumindWNT<-cumsum(indWNT)
    cumindWN<-cumsum(indWN)
    cumindWT<-cumsum(indWT)
    if(i8!=0){ # having 3 effects
       # 1st, the WCE effect of the variables with 3 effects 
      
          #pos5: positions of variables having 3 effects among those having WCE to locate the position for cutoff
         pos6<-match(variablesNEW[NL == 1 & TD == 1 & WCE==1], variablesNEW[WCE == 1 & NL==1])
         #pos6: positions of variables having 3 effect among those having 3 effects and WCE+NL to locate the position 
         #of the spline basis for Nbasis and Dbasis
        
        for(k in 1:i8){ # for the covariates having 3 effects, need to combine the WCE basis and the NL TVC first 
          # in the first step, assume the TVC have Linear effect 
          # calculate the Dj(u) with WCE basis and TVC
          
           kal<- do.call("rbind", lapply(1:length(Id),  # combine the calcluated D_j(u) of each paitents at risk set together
                                    function(j) .wcecalc(listeT[-1], data[data[,1]==Id[j], variablesNEW[NL == 1 & TD == 
                                    1 & WCE==1][k]], data[data[,1]==Id[j],TypeNEW[2]], Dbasis[[pos6[k]]], cutoff[pos5[k]])))
           kal <- kal[is.na(kal[, 1]) == FALSE, ] # get rid of all the patients that do not belong to any risk set 
          
           flag<-dim(QWR)[2]+1
          if(constrained[pos5[k]]=="FALSE"){
            indWNT[k+1]<-nWCEb
            QWR<-cbind(QWR, kal)
            covp<-paste("QWR[,", flag: dim(QWR)[2], "]", sep = "" )
            tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))
          }else if(constrained[pos5[k]]=="Right"){
            indWNT[k+1]<-nWCEb-2
            QWR<-cbind(QWR, kal[,1:(dim(kal)[2]-2)])
            covp<-paste("QWR[,", flag: dim(QWR)[2], "]", sep = "" )
            tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))
          }else{
            indWNT[k+1]<-nWCEb-2
            QWR<-cbind(QWR, kal[,3:(dim(kal)[2])])
            covp<-paste("QWR[,", flag: dim(QWR)[2], "]", sep = "" )
            tt2<-paste(tt2, paste(c(covp,""), collapse= "+")) 
          }
           cumindWNT<-cumsum(indWNT)
        }
        
         
        #2nd, the WCE effects of those only having 2 effects
         if(i6!=0){ #if exist variables with only WCE and NL effect 
           pos7<-match(variablesNEW[NL == 1 & TD == 0 & WCE==1], variablesNEW[WCE == 1 & NL==1])
           #pos7: positions of variables having only WCE+NL effect among those having 3 effects and WCE+NL to locate the position 
           #of the spline basis for Nbasis
           #pos8: positions of variables having WCE+NL effects among those having WCE to locate the position for cutoff
          
           for(k in 1:i6){ # for the covariates having WCE+NL effects, need to combine the WCE basis and the TVC first 
            #in 1st step, assume TVC have linear effects
             
             kal<- do.call("rbind", lapply(1:length(Id),  # combine the calcluated D_j(u) of each paitents at risk set together
                              function(j) .wcecalc(listeT[-1], data[data[,1]==Id[j], variablesNEW[NL == 1 & TD == 
                                            0 & WCE==1][k]], data[data[,1]==Id[j],TypeNEW[2]], Dbasis[[pos7[k]]], cutoff[pos8[k]])))
             kal <- kal[is.na(kal[, 1]) == FALSE, ] # get rid of all the patients that do not belong to any risk set 
             
             flag<-dim(QWR)[2]+1
             if(constrained[pos8[k]]=="FALSE"){
               indWN[k+1]<-nWCEb
               QWR<-cbind(QWR, kal)
             }else if(constrained[pos8[k]]=="Right"){
               indWN[k+1]<-nWCEb-2
               QWR<-cbind(QWR, kal[,1:(dim(kal)[2]-2)])
             }else{
               indWN[k+1]<-nWCEb-2
               QWR<-cbind(QWR, kal[,3:(dim(kal)[2])])
             }
             covp<-paste("QWR[,", flag: dim(QWR)[2], "]", sep = "" )
             tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))
           }
           cumindWN<-cumsum(indWN)
         }
        
         
         if(i7!=0){ #if exist variable having WCE+TD
           pos9<-match(variablesNEW[NL == 0 & TD == 1 & WCE==1], variablesNEW[WCE == 1 & NL==0])
           #pos9: positions of variables having only WCE+TD effect among those having WCE no NL to locate the position 
           #of the spline basis for QWR
           #pos2:mark position of variable with WCE+TD among all varialbes having WCE
           
           for(k in 1:i7){
             if(constrained[pos2[k]]=="FALSE"){
               indWT[k+1]<-nWCEb
             covp<-paste( "QWR[,", (dim(data)[2] + (pos9[k] -1) * nWCEb + 1):(dim(data)[2] + pos9[k] * nWCEb)  , "]", sep = "")
             }else if(constrained[pos2[k]]=="Right"){
               indWT[k+1]<-nWCEb-2
               covp<-paste( "QWR[,", (dim(data)[2] + (pos9[k] -1) * nWCEb + 1):(dim(data)[2] + pos9[k] * nWCEb-2)  , "]", sep = "")
             }else{
               indWT[k+1]<-nWCEb-2
               covp<-paste( "QWR[,", (dim(data)[2] + (pos9[k] -1) * nWCEb + 3):(dim(data)[2] + pos9[k] * nWCEb)  , "]", sep = "")
             }
             tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
           }
           cumindWT<-cumsum(indWT)
         }
         
         
         modX <- coxph(eval(parse(text = substr(tt2, 13, nchar(tt2) - 
                                                  1))), method = "efron")
         vrais <- c(vrais, modX$loglik[2]) # estimate all the WCE effects first
         
         
         #rescale the estimated weights and nonlinear effects when they both exist
         #new revision, always rescale the weights when there is weight+additonal effects
         # when there is only weight function, the rescale should be done after the model is fitted???
         ## yes, when there is only weights, the rescale should be done at the last step, and store the 
         ## rescale factor at the coefficients part
         
         #after each estimation, the estimated spline coefficients for weights (only for weight, no NL coefficients
         # are rescaled ????) are 
         #rescaled by divided by the sum of the corresponding estimated functions
         #for weights, divided by sum of the resulting weight at each time point
        
         for (j in 1:i8){## rescale for the weight when there is both TD, NL and weight, the spline basis is different for both weight and NL effect and for weight or weight and TD  
             if(constrained[pos5[j]]=="FALSE"){
               rescale<-as.vector(Dbasis[[(i6+i)]]%*%modX$coef[(i1 + cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) + cumindWNT[j]
                                                                + 1):(i1 + cumindW[i5+1]+ i2 * nNLb +i3* (nNLb+1)+cumindWNT[j]+indWNT[j+1])])
               modX$coef[(i1 + cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) + cumindWNT[j]
                          + 1):(i1 + cumindW[i5+1]+ i2 * nNLb +i3* (nNLb+1)+cumindWNT[j]+indWNT[j+1])]<-modX$coef[(i1 + cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) + cumindWNT[j]
                                                                                                                   + 1):(i1 + cumindW[i5+1]+ i2 * nNLb +i3* (nNLb+1)+cumindWNT[j]+indWNT[j+1])]/sum(rescale)
             }else if(constrained[pos5[j]]=="Right"){
               rescale<-as.vector(Dbasis[[(i6+i)]]%*%c(modX$coef[(i1 + cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) + cumindWNT[j]
                                                                  + 1):(i1 + cumindW[i5+1]+ i2 * nNLb +i3* (nNLb+1)+cumindWNT[j]+indWNT[j+1])],0,0))
               modX$coef[(i1 + cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) + cumindWNT[j]
                          + 1):(i1 + cumindW[i5+1]+ i2 * nNLb +i3* (nNLb+1)+cumindWNT[j]+indWNT[j+1])]<- modX$coef[(i1 + cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) + cumindWNT[j]
                                                                                                                    + 1):(i1 + cumindW[i5+1]+ i2 * nNLb +i3* (nNLb+1)+cumindWNT[j]+indWNT[j+1])]/sum(rescale)
             }else{
               rescale<-as.vector(Dbasis[[(i6+i)]]%*%c(0,0,modX$coef[(i1 + cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) + cumindWNT[j]
                                                                      + 1):(i1 + cumindW[i5+1]+ i2 * nNLb +i3* (nNLb+1)+cumindWNT[j]+indWNT[j+1])]))
               modX$coef[(i1 + cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) + cumindWNT[j]
                          + 1):(i1 + cumindW[i5+1]+ i2 * nNLb +i3* (nNLb+1)+cumindWNT[j]+indWNT[j+1])] <-modX$coef[(i1 + cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) + cumindWNT[j]
                                                                                                                    + 1):(i1 + cumindW[i5+1]+ i2 * nNLb +i3* (nNLb+1)+cumindWNT[j]+indWNT[j+1])]/sum(rescale)
             }
           res[i5+j]<-sum(rescale)
         }
         
          if(i6!=0){ # NLWCE
           for (j in 1:i6){ ## rescale for the weight when there is both NL and TD effects
             if(constrained[pos8[j]]=="FALSE" ){
               rescale<-as.vector(Dbasis[[(i)]]%*%modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                                                             + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])])
               modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                          + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])] <- modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                                                                                                                                      + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])]/sum(rescale)
             }else if(constrained[pos8[j]]=="Right"){
               rescale<-as.vector(Dbasis[[(i)]]%*%c(modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                                                               + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])],0,0))
               modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                          + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])] <- modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                                                                                                                                      + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])]/sum(rescale)
             }else{
               rescale<-as.vector(Dbasis[[(i)]]%*%c(0,0, modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                                                                    + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])]))
               modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                          + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])]<-modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                                                                                                                                    + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])]/sum(rescale)
             }
             res[i5+i8+j]<-sum(rescale)
           }
            
          }
         
         if(i7!=0){ # rescale for the weights when there is both weight and TD effect 
           for(j in 1:i7){
             if(constrained[pos9[j]]=="FALSE" ){
               rescale<-as.vector(Bbasis[[(i5+i)]]%*%modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2*nNLb+i3*(nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                                                                + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2*nNLb+i3* (nNLb+1) +cumindWN[i6+1]+cumindWT[j]+indWT[j+1])])
               modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                          + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[i6+1]+cumindWT[j]+indWT[j+1])] <- modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                                                                                                                                                     + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[i6+1]+cumindWT[j]+indWT[j+1])]/sum(rescale)
             }else if(constrained[pos9[j]]=="Right"){
               rescale<-as.vector( Bbasis[[(i5+i)]]%*%c(modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                                                                   + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[i6+1]+cumindWT[j]+indWT[j+1])],0,0))
               modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                          + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[i6+1]+cumindWT[j]+indWT[j+1])] <- modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                                                                                                                                                     + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[i6+1]+cumindWT[j]+indWT[j+1])]/sum(rescale)
             }else{
               rescale<-as.vector( Bbasis[[(i5+i)]]%*%c(0,0, modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                                                                        + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]+indWT[j+1])]))
               modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                          + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[i6+1]+cumindWT[j]+indWT[j+1])]<-modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                                                                                                                                                   + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[i6+1]+cumindWT[j]+indWT[j+1])]/sum(rescale)
             }
             res[i5+i8+i6+j]<-sum(rescale)
           }
           
         }
         
         # print(modX$coefficients)
         # print(modX$coef)
         #############################################
         # Then estimate NL effects (TDNLWCE+TDNL+NLWCE) (i8+i6+i4)
         ######################################
         modWCE<-modX
         tt2<-tt
         VV<-matrix(ncol=(i4+i6+i8)*(nknotsTDNL+pTDNL),nrow=dim(QWR)[1]) #number of colum is the number of variables having NL
         # effect multiply by the number of spline basis for the NL effect 
         
         sumWCEbase3<-list() #sum WCE base for variables having 3 effects
         
         for (k in 1:i8){ #estimate NL effect of variables having 3 effects 
           #pos5: positions of variables having 3 effects among those having WCE to locate the position for cutoff
           #pos6: positions of variables having 3 effect among those having 3 effects and WCE+NL to locate the position 
           #of the spline basis for Nbasis and Dbasis
           if(constrained[pos5[k]]=="FALSE"){
           sumWCEbase3[[k]]<- matrix(unlist(Dbasis[[pos6[k]]]), ncol = (nknotsWCE + pWCE+1), byrow = FALSE)%*%modWCE$coef[(Nn+cumindWNT[k]+1):(Nn+cumindWNT[k]+indWNT[k+1])]
           #dimension of the matrix in sumWCEbase is cutoff*1 
           } else if(constrained[pos5[k]]=="Right"){
             sumWCEbase3[[k]]<- matrix(unlist(Dbasis[[pos6[k]]]), ncol = nWCEb, byrow = FALSE)[,1:(nWCEb-2)]%*%modWCE$coef[(Nn+cumindWNT[k]+1):(Nn+cumindWNT[k]+indWNT[k+1])]
           }else{
             sumWCEbase3[[k]]<- matrix(unlist(Dbasis[[pos6[k]]]), ncol = nWCEb, byrow = FALSE)[,3:(nWCEb)]%*%modWCE$coef[(Nn+cumindWNT[k]+1):(Nn+cumindWNT[k]+indWNT[k+1])]
           }
           NLbase<-matrix(unlist(Nbasis[[pos6[k]]]), ncol = (nknotsTDNL+pTDNL), byrow = FALSE)
           
           for (k2 in 1:(nknotsTDNL+pTDNL)){ 
             
            
           kal<- do.call("rbind", lapply(1:length(Id),  # combine the calcluated D_j(u) of each paitents at risk set together
                                         function(j) .wcecalc(listeT[-1], NLbase[data[,1]==Id[j],k2] , 
                                                              data[data[,1]==Id[j],TypeNEW[2]], sumWCEbase3[[k]], cutoff[pos5[k]])))
           kal <- kal[is.na(kal[, 1]) == FALSE, ] # get rid of all the patients that do not belong to any risk set 
           
           
           VV[,(nknotsTDNL+pTDNL)*(k-1)+k2]<-kal #calculat E(u) for each of the NL spline coefficient
             
             
            tt2<-paste(tt2,"VV[,",(nknotsTDNL+pTDNL)*(k-1)+k2,"]+",sep="")
           }
         }
       
         if(i6!=0){
           sumWCEbaseWN<-list()  # sum WCE base for variables having only 2 effects with WCE
         for(k in 1:i6){#estimate NL effect of variables having NL+WCE effects
           #pos7: positions of variables having only WCE+NL effect among those having 3 effects and WCE+NL to locate the position 
           #of Nbasis and Dbasis
           #pos8: positions of variables having WCE+NL effects among those having WCE to locate the position for cutoff
           if(constrained[pos8[k]]=="FALSE"){
           sumWCEbaseWN[[k]]<- matrix(unlist(Dbasis[[pos7[k]]]), ncol = nWCEb, byrow = FALSE)%*%modWCE$coef[(Nn+cumindWNT[i8+1]+cumindWN[k]+1):(Nn+ cumindWNT[i8+1]+cumindWN[k]+indWN[k+1])]
           }else if(constrained[[pos8[k]]]=="Right"){
             sumWCEbaseWN[[k]]<- matrix(unlist(Dbasis[[pos7[k]]]), ncol =nWCEb, byrow = FALSE)[,1:(nWCEb-2)]%*%modWCE$coef[(Nn+cumindWNT[i8+1]+cumindWN[k]+1):(Nn+cumindWNT[i8+1]+cumindWN[k]+indWN[k+1])]
           }else{
             sumWCEbaseWN[[k]]<- matrix(unlist(Dbasis[[pos7[k]]]), ncol =nWCEb, byrow = FALSE)[,3:nWCEb]%*%modWCE$coef[(Nn+cumindWNT[i8+1]+cumindWN[k]+1):(Nn+cumindWNT[i8+1]+cumindWN[k]+indWN[k+1])]
           }
           NLbase<-matrix(unlist(Nbasis[[pos7[k]]]), ncol = (nknotsTDNL+pTDNL), byrow = FALSE)
           
            for (k2 in 1:(nknotsTDNL+pTDNL)){ 
             
             
            kal<- do.call("rbind", lapply(1:length(Id),  # combine the calcluated D_j(u) of each paitents at risk set together
                                           function(j) .wcecalc(listeT[-1], NLbase[data[,1]==Id[j],k2] , 
                                                                data[data[,1]==Id[j],TypeNEW[2]], sumWCEbaseWN[[k]], cutoff[pos8[k]])))
             kal <- kal[is.na(kal[, 1]) == FALSE, ] # get rid of all the patients that do not belong to any risk set 
             
             
             VV[,(nknotsTDNL+pTDNL)*(i8+k-1)+k2]<-kal #calculat E(u) for each of the NL spline coefficient
             
             
             tt2<-paste(tt2,"VV[,",(nknotsTDNL+pTDNL)*(i8+k-1)+k2,"]+",sep="")
           }
          }
         }   
         
         if(i4!=0){ # estimate NL effect for variables having NLTD effects
           pos10<-match(variablesNEW[NL == 1 & TD == 1 & WCE==0], variablesNEW[WCE == 0&NL==1])
           #pos10 marks the position of variables having NL TD among those having NL no WCE to locate the NL spline basis in QWR
         for(k in 1:i4){ 
           VV[,(nNLb*(i8+i6+k-1)+1):(nNLb*(i8+i6+k))]<-QWR[,(dim(data)[2]+(i5+i7)*nWCEb +(pos10[k]-1) *nNLb +1):(dim(data)[2]+(i5+i7)*nWCEb +pos10[k]*nNLb)]
           covp<-paste( "VV[,", (nNLb*(i8+i6+k-1)+1):(nNLb*(i8+i6+k)), "]", sep = "")
           tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
          }
        } 
         

      
         modX<-coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
         vrais<-c(vrais,modX$loglik[2]) # estimate the NL effects conditional on the WCE effects
         
         ############################################
         # estimate TD effects conditional on NL WCE(TDNLWCE+TDWCE+TDNL) (i8+i4+i7)
         modNL<-modX
         tt2<-tt ## VV stores the TD spline
         VV<-matrix(ncol=(i8+i4+i7)*(nNLb+1),nrow=dim(QWR)[1]) 
         
         sumNLbaseW<-matrix(NA,ncol=i8,nrow = length(data[,1]))
         sumNLbaseNW<-matrix(NA,ncol=i4,nrow = length(QWR[,1]))
         
         ##########################################################
         ## double check whether sumWCEbase need to be redefined here
         #########################################################
         
         for (k in 1:i8){ #estimate TD effect of variables with 3 effects
           #pos6: positions of variables having 3 effect among those having 3 effects and WCE+NL to locate the position 
           #of the spline basis for Nbasis and Dbasis
           #pos5: positions of variables having 3 effects among those having WCE to locate the position for cutoff
           sumNLbaseW[,k]<- matrix(unlist(Nbasis[[pos6[k]]]), ncol = (nknotsTDNL+pTDNL),
                                  byrow = FALSE)%*%modNL$coef[(Nn+(k-1)*(nknotsTDNL+pTDNL)+1):(Nn+k*(nknotsTDNL+pTDNL))]
  
           
           kal<- do.call("rbind", lapply(1:length(Id),  #calculate the WCE incorperationg both estimated coefficients for WCE and NL spline
                                         function(j) .wcecalc(listeT[-1], sumNLbaseW[data[,1]==Id[j],k] , 
                                                     data[data[,1]==Id[j],TypeNEW[2]], sumWCEbase3[[k]], cutoff[pos5[k]])))
           kal <- kal[is.na(kal[, 1]) == FALSE, ] # get rid of all the patients that do not belong to any risk set 
           
           VV[,((1+nNLb)*(k-1)+1):((1+nNLb)*k)]<-kal*QWR[,(dim(data)[2]+(i5+i7)*nWCEb+(i2+i4)*nNLb+1):(dim(data)[2]+(i5+i7)*nWCEb+(i2+i4+1)*nNLb+1)]
                                                   
           covp<-paste( "VV[,", ((1+nNLb)*(k-1)+1):((1+nNLb)*k), "]", sep = "")
           tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
         }
        
         if(i7!=0){ #estimate TD effect of variables having TD+WCE
           #pos9: positions of variables having only WCE+TD effect among those having WCE no NL to locate the position 
           #of the spline basis for QWR 
           sumWCEbaseWT<-list()
           for(k in 1:i7){
             if(constrained[pos2[k]]=="FALSE"){
             sumWCEbaseWT[[k]]<- QWR[,c((dim(data)[2]+(pos9[k]-1)*nWCEb+1):(dim(data)[2]+pos9[k]*nWCEb))]%*%modWCE$coef[(Nn+cumindWNT[i8+1]
                                          +cumindWN[i6+1]+cumindWT[k]+1):(Nn+cumindWNT[i8+1]+cumindWN[i6+1]+cumindWT[k]+indWT[k+1])]
             }else if(constrained[pos2[k]]=="Right"){
             sumWCEbaseWT[[k]]<- QWR[,c((dim(data)[2]+(pos9[k]-1)*nWCEb+1):(dim(data)[2]+pos9[k]*nWCEb))][,1:(nWCEb-2)]%*%modWCE$coef[(Nn+cumindWNT[i8+1]
                                                               +cumindWN[i6+1]+cumindWT[k]+1):(Nn+cumindWNT[i8+1]+cumindWN[i6+1]+cumindWT[k]+indWT[k+1])]
             }else{
               sumWCEbaseWT[[k]]<- QWR[,c((dim(data)[2]+(pos9[k]-1)*nWCEb+1):(dim(data)[2]+pos9[k]*nWCEb))][,3:nWCEb]%*%modWCE$coef[(Nn+cumindWNT[i8+1]
                                                              +cumindWN[i6+1]+cumindWT[k]+1):(Nn+cumindWNT[i8+1]+cumindWN[i6+1]+cumindWT[k]+indWT[k+1])]
             }
             VV[,((1+nNLb)*(i8+k-1)+1):((1+nNLb)*(i8+k))]<-QWR[,(dim(data)[2]+(i5+i7)*nWCEb+(i2+i4)*nNLb+1):(dim(data)[2]+(i5+i7)*nWCEb+(i2+i4+1)*nNLb+1)
                                                               ]*as.vector(matrix(unlist(sumWCEbaseWT[[k]]), ncol =1, byrow = FALSE))
             covp<-paste( "VV[,", ((1+nNLb)*(i8+k-1)+1):((1+nNLb)*(i8+k)), "]", sep = "")
             tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
           }
         }
         if(i4!=0){#estimate TD effect of variables having TD+NL
           #pos10 marks the position of variables having NL TD among those having NL no WCE to locate the NL spline basis in QWR
           for(k in 1:i4){
             sumNLbaseNW[,k]<- QWR[,c((dim(data)[2]+(i5+i7)*nWCEb+(pos10[k]-1)*nNLb+1):(dim(data)[2]
                                +(i5+i7)*nWCEb+pos10[k]*nNLb))]%*%modNL$coef[(Nn+(i8+i6+k-1)*nNLb+1):(Nn+(i8+i6+k)*nNLb)]
             VV[,((1+nNLb)*(i8+i7+k-1)+1):((1+nNLb)*(i8+i7+k))]<-sumNLbaseNW[,k]*QWR[,(dim(data)[2]+(i5+i7)*nWCEb+(i2+i4)*nNLb+1):(dim(data)[2]+(i5+i7)*nWCEb+(i2+i4+1)*nNLb+1)]
            
             covp<-paste( "VV[,", ((1+nNLb)*(i8+i7+k-1)+1):((1+nNLb)*(i8+i7+k)), "]", sep = "")
             tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
           }
         }  
         

         
         modX<-coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
         vrais<-c(vrais,modX$loglik[2]) # estimate the TEL effects conditional on the TD and NL effects
         
         diff<-1
         #################################
         #### double checked upto here
         ################################   
         while(diff>0.01){ ## ACE algorithm until converge 
           
           ##########################################################################
           # estimate WCE effect first conditional on previouly estimated TD, NL effect
           modTD<-modX
           tt2<-tt
           VV<-matrix(ncol=cumindWNT[i8+1]+cumindWN[i6+1]+cumindWT[i7+1],nrow=dim(QWR)[1]) 
           
           sumNLbase<-matrix(NA,ncol=(i8+i6),nrow = length(data[,1]))
           sumTDbase<-matrix(NA,ncol=(i8+i7),nrow = length(QWR[,1]))
           for (k in 1:i8){
             #pos6: positions of variables having 3 effect among those having 3 effects and WCE+NL to locate the position 
             #of the spline basis for Nbasis and Dbasis
             sumNLbase[,k]<- matrix(unlist(Nbasis[[pos6[k]]]), ncol = (nknotsTDNL+pTDNL)
                                    ,byrow = FALSE)%*%modNL$coef[(Nn+(k-1)*(nknotsTDNL+pTDNL)+1):(Nn+k*(nknotsTDNL+pTDNL))]
             
             sumTDbase[,k]<- QWR[,c((dim(data)[2]+(i5+i7)*(nknotsWCE + pWCE+1)+(i2+i4)*(nknotsTDNL+pTDNL)
                              +1):(dim(data)[2]+(i5+i7)*(nknotsWCE + pWCE+1)+(i2+i4+1)*(nknotsTDNL+pTDNL)+1))]%*%modTD$coef[(Nn+(k-1)*(nknotsTDNL+pTDNL+1)+1):(Nn+k*(nknotsTDNL+pTDNL+1))]
             
             kal<- do.call("rbind", lapply(1:length(Id), #combine the WCE spline basis with estimated NL effects
                                           function(j) .wcecalc(listeT[-1],  sumNLbase[data[,1]==Id[j],k] , 
                                                                data[data[,1]==Id[j],TypeNEW[2]], Dbasis[[pos6[k]]], cutoff[pos5[k]])))
             kal <- kal[is.na(kal[, 1]) == FALSE, ] # get rid of all the patients that do not belong to any risk set 
            
             
             if(constrained[pos5[k]]=="FALSE"){
               VV[,(cumindWNT[k]+1):(cumindWNT[k]+indWNT[k+1])]<-sumTDbase[,k]*kal
               covp<-paste( "VV[,",(cumindWNT[k]+1):(cumindWNT[k]+indWNT[k+1]), "]", sep = "")
             }else if(constrained[pos5[k]]=="Right"){
               VV[,(cumindWNT[k]+1):(cumindWNT[k]+indWNT[k+1])]<-sumTDbase[,k]*kal[,1:(nWCEb-2)]
               covp<-paste( "VV[,",(cumindWNT[k]+1):(cumindWNT[k]+indWNT[k+1]), "]", sep = "")
             }else{
               VV[,(cumindWNT[k]+1):(cumindWNT[k]+indWNT[k+1])]<-sumTDbase[,k]*kal[,3:nWCEb]
               covp<-paste( "VV[,",(cumindWNT[k]+1):(cumindWNT[k]+indWNT[k+1]), "]", sep = "")
             }
             tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
           }
           
           
           if(i6!=0){ #if exist variables with only WCE and NL effect 
             #pos7: positions of variables having only WCE+NL effect among those having 3 effects and WCE+NL to locate the position 
             #of the spline basis for Nbasis
             #pos8: positions of variables having WCE+NL effects among those having WCE to locate the position for cutoff
             for(k in 1:i6){ # for the covariates having WCE+NL effects, need to combine the WCE basis and the NL basis first 
               # pNLTDWCE indicate the position of the varialbes having WCE+NL+TD among the variables with NL 
               # calculate the Dj(u) with WCE basis and NL basis
               ## set the initial coefficents for NL basis to 1 
               sumNLbase[,(i8+k)]<- matrix(unlist(Nbasis[[pos7[k]]]), ncol = (nknotsTDNL+pTDNL)
                                      , byrow = FALSE)%*%modNL$coef[(Nn+(i8+k-1)*(nknotsTDNL+pTDNL)+1):(Nn+(i8+k)*(nknotsTDNL+pTDNL))]
               
               kal<- do.call("rbind", lapply(1:length(Id),  # combine the calcluated D_j(u) of each paitents at risk set together
                                             function(j) .wcecalc(listeT[-1],  sumNLbase[data[,1]==Id[j],(i8+k)] , 
                                                                  data[data[,1]==Id[j],TypeNEW[2]], Dbasis[[pos7[k]]], cutoff[pos8[k]])))
               kal <- kal[is.na(kal[, 1]) == FALSE, ] # get rid of all the patients that do not belong to any risk set 
               
               if(constrained[pos8[k]]=="FALSE"){
                 VV[,(cumindWNT[i8+1]+cumindWN[k]+1):(cumindWNT[i8+1]+cumindWN[k]+indWNT[k+1])]<-kal
                 covp<-paste( "VV[,",(cumindWNT[i8+1]+cumindWN[k]+1):(cumindWNT[i8+1]+cumindWN[k]+indWNT[k+1]), "]", sep = "")
               }else if(constrained[pos8[k]]=="Right"){
                 VV[,(cumindWNT[i8+1]+cumindWN[k]+1):(cumindWNT[i8+1]+cumindWN[k]+indWNT[k+1])]<-kal[,1:(nWCEb-2)]
                 covp<-paste( "VV[,",(cumindWNT[i8+1]+cumindWN[k]+1):(cumindWNT[i8+1]+cumindWN[k]+indWNT[k+1]), "]", sep = "")
               }else{
                 VV[,(cumindWNT[i8+1]+cumindWN[k]+1):(cumindWNT[i8+1]+cumindWN[k]+indWNT[k+1])]<-sumTDbase[,k]*kal[,3:nWCEb]
                 covp<-paste( "VV[,",(cumindWNT[i8+1]+cumindWN[k]+1):(cumindWNT[i8+1]+cumindWN[k]+indWNT[k+1]), "]", sep = "")
               }
               tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
             }
           }
             
           if(i7!=0){ #if exist variable having WCE+TD
             #pos9: positions of variables having only WCE+TD effect among those having WCE no NL to locate the position 
             #of the spline basis for QWR
             for(k in 1:i7){
               sumTDbase[,(i8+k)]<- QWR[,c((dim(data)[2]+(i5+i7)*(nknotsWCE + pWCE+1)+(i2+i4)*(nknotsTDNL+pTDNL)
                                       +1):(dim(data)[2]+(i5+i7)*(nknotsWCE + pWCE+1)+(i2+i4+1)*(nknotsTDNL+pTDNL)+1))]%*%modTD$coef[(Nn+(i8+k-1)*(nknotsTDNL+pTDNL+1)+1):(Nn+(i8+k)*(nknotsTDNL+pTDNL+1))]
              
              if(constrained[pos2[k]]=="FALSE"){
                 VV[,(cumindWNT[i8+1]+cumindWN[i6+1]+cumindWT[k]+1):(cumindWNT[i8+1]+cumindWN[i6+1]+cumindWT[k]+indWNT[k+1])]<-sumTDbase[
                                                                                          ,(i8+k)]*QWR[, (dim(data)[2] + (pos9[k] -1)*nWCEb +1):(dim(data)[2] +pos9[k]*nWCEb)]
                 covp<-paste( "VV[,",(cumindWNT[i8+1]+cumindWN[i6+1]+cumindWT[k]+1):(cumindWNT[i8+1]+cumindWN[i6+1]+cumindWT[k]+indWNT[k+1]), "]", sep = "")
               
               }else if(constrained[pos2[k]]=="Right"){
                 VV[,(cumindWNT[i8+1]+cumindWN[i6+1]+cumindWT[k]+1):(cumindWNT[i8+1]+cumindWN[i6+1]+cumindWT[k]+indWNT[k+1])]<-sumTDbase[
                   ,(i8+k)]*QWR[,(dim(data)[2] + (pos9[k] -1)*nWCEb +1):(dim(data)[2] +pos9[k]*nWCEb-2)]
                 covp<-paste( "VV[,",(cumindWNT[i8+1]+cumindWN[i6+1]+cumindWT[k]+1):(cumindWNT[i8+1]+cumindWN[i6+1]+cumindWT[k]+indWNT[k+1]), "]", sep = "")
          
               }else{
                 VV[,(cumindWNT[i8+1]+cumindWN[i6+1]+cumindWT[k]+1):(cumindWNT[i8+1]+cumindWN[i6+1]+cumindWT[k]+indWNT[k+1])]<-sumTDbase[
                   ,(i8+k)]*QWR[,(dim(data)[2] + (pos9[k] -1)*nWCEb +3):(dim(data)[2] +pos9[k]*nWCEb)]
                 covp<-paste( "VV[,",(cumindWNT[i8+1]+cumindWN[i6+1]+cumindWT[k]+1):(cumindWNT[i8+1]+cumindWN[i6+1]+cumindWT[k]+indWNT[k+1]), "]", sep = "")
             
               }
               tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
             } #WCE basis multiply the estimated value of TD effects
           } 
             
             
     
           modX<-coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
           vrais<-c(vrais,modX$loglik[2]) 
           
            #NL+TD+WCE
             for (j in 1:i8) {
               if(constrained[pos5[j]]=="FALSE"){
                 rescale<-as.vector(Dbasis[[(i6+i)]]%*%modX$coef[(i1 + cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) + cumindWNT[j]
                                                                  + 1):(i1 + cumindW[i5+1]+ i2 * nNLb +i3* (nNLb+1)+cumindWNT[j]+indWNT[j+1])])
                 modX$coef[(i1 + cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) + cumindWNT[j]
                            + 1):(i1 + cumindW[i5+1]+ i2 * nNLb +i3* (nNLb+1)+cumindWNT[j]+indWNT[j+1])]<-modX$coef[(i1 + cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) + cumindWNT[j]
                                                                                                                     + 1):(i1 + cumindW[i5+1]+ i2 * nNLb +i3* (nNLb+1)+cumindWNT[j]+indWNT[j+1])]/sum(rescale)
               }else if(constrained[pos5[j]]=="Right"){
                 rescale<-as.vector(Dbasis[[(i6+i)]]%*%c(modX$coef[(i1 + cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) + cumindWNT[j]
                                                                    + 1):(i1 + cumindW[i5+1]+ i2 * nNLb +i3* (nNLb+1)+cumindWNT[j]+indWNT[j+1])],0,0))
                 modX$coef[(i1 + cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) + cumindWNT[j]
                            + 1):(i1 + cumindW[i5+1]+ i2 * nNLb +i3* (nNLb+1)+cumindWNT[j]+indWNT[j+1])]<- modX$coef[(i1 + cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) + cumindWNT[j]
                                                                                                                      + 1):(i1 + cumindW[i5+1]+ i2 * nNLb +i3* (nNLb+1)+cumindWNT[j]+indWNT[j+1])]/sum(rescale)
               }else{
                 rescale<-as.vector(Dbasis[[(i6+i)]]%*%c(0,0,modX$coef[(i1 + cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) + cumindWNT[j]
                                                                        + 1):(i1 + cumindW[i5+1]+ i2 * nNLb +i3* (nNLb+1)+cumindWNT[j]+indWNT[j+1])]))
                 modX$coef[(i1 + cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) + cumindWNT[j]
                            + 1):(i1 + cumindW[i5+1]+ i2 * nNLb +i3* (nNLb+1)+cumindWNT[j]+indWNT[j+1])] <-modX$coef[(i1 + cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) + cumindWNT[j]
                                                                                                                      + 1):(i1 + cumindW[i5+1]+ i2 * nNLb +i3* (nNLb+1)+cumindWNT[j]+indWNT[j+1])]/sum(rescale)
               }
               res[i5+j]<-sum(rescale)
             }
           
           if(i6!=0){ # NLWCE
             for (j in 1:i6){
               if(constrained[pos8[j]]=="FALSE" ){
                 rescale<-as.vector(Dbasis[[(i)]]%*%modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                                                               + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])])
                 modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                            + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])] <- modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                                                                                                                                        + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])]/sum(rescale)
               }else if(constrained[pos8[j]]=="Right"){
                 rescale<-as.vector(Dbasis[[(i)]]%*%c(modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                                                                 + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])],0,0))
                 modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                            + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])] <- modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                                                                                                                                        + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])]/sum(rescale)
               }else{
                 rescale<-as.vector(Dbasis[[(i)]]%*%c(0,0, modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                                                                      + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])]))
                 modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                            + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])]<-modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                                                                                                                                      + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])]/sum(rescale)
               }
               res[i5+i8+j]<-sum(rescale)
             }                     
           }
           
           
           if(i7!=0){
             for(j in 1:i7){
               if(constrained[pos9[j]]=="FALSE" ){
                 rescale<-as.vector(Bbasis[[(i5+i)]]%*%modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2*nNLb+i3*(nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                                                                  + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2*nNLb+i3* (nNLb+1) +cumindWN[i6+1]+cumindWT[j]+indWT[j+1])])
                 modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                            + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[i6+1]+cumindWT[j]+indWT[j+1])] <- modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                                                                                                                                                       + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[i6+1]+cumindWT[j]+indWT[j+1])]/sum(rescale)
               }else if(constrained[pos9[j]]=="Right"){
                 rescale<-as.vector( Bbasis[[(i5+i)]]%*%c(modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                                                                     + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[i6+1]+cumindWT[j]+indWT[j+1])],0,0))
                 modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                            + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[i6+1]+cumindWT[j]+indWT[j+1])] <- modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                                                                                                                                                       + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[i6+1]+cumindWT[j]+indWT[j+1])]/sum(rescale)
               }else{
                 rescale<-as.vector( Bbasis[[(i5+i)]]%*%c(0,0, modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                                                                          + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]+indWT[j+1])]))
                 modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                            + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[i6+1]+cumindWT[j]+indWT[j+1])]<-modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                                                                                                                                                     + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[i6+1]+cumindWT[j]+indWT[j+1])]/sum(rescale)
               }
               res[i5+i8+i6+j]<-sum(rescale)
             }
           }
           
           
           
           
           # print(modX$coefficients)
           # print(modX$coef)
           ########################################################
           # estimate NL effect conditional on the WCE and TD effects
           
           modWCE<-modX
           tt2<-tt
           VV<-matrix(ncol=(i4+i6+i8)*(nknotsTDNL+pTDNL),nrow=dim(QWR)[1])
           sumWCEbase3<-list()
           sumTDbase<-matrix(NA,ncol=(i8+i4),nrow = length(QWR[,1]))
           for (k in 1:i8){ #estimate NL effect of variables having 3 effects 
             #pos5: positions of variables having 3 effects among those having WCE to locate the position for cutoff
             #pos6: positions of variables having 3 effect among those having 3 effects and WCE+NL to locate the position 
             #of the spline basis for Nbasis and Dbasis
             
             if(constrained[pos5[k]]=="FALSE"){
               sumWCEbase3[[k]]<- matrix(unlist(Dbasis[[pos6[k]]]), ncol = nWCEb, byrow = FALSE)%*%modWCE$coef[(Nn+cumindWNT[k]+1):(Nn+cumindWNT[k]+indWNT[k+1])]
               #dimension of the matrix in sumWCEbase is cutoff*1 
             } else if(constrained[pos5[k]]=="Right"){
               sumWCEbase3[[k]]<- matrix(unlist(Dbasis[[pos6[k]]]), ncol = nWCEb, byrow = FALSE)[,1:(nWCEb-2)]%*%modWCE$coef[(Nn+cumindWNT[k]+1):(Nn+cumindWNT[k]+indWNT[k+1])]
             }else{
               sumWCEbase3[[k]]<- matrix(unlist(Dbasis[[pos6[k]]]), ncol = nWCEb, byrow = FALSE)[,3:(nWCEb)]%*%modWCE$coef[(Nn+cumindWNT[k]+1):(Nn+cumindWNT[k]+indWNT[k+1])]
             }
             #dimension of the matrix in sumWCEbase3 is cutoff*1 
             
             
             NLbase<-matrix(unlist(Nbasis[[pos6[k]]]), ncol = (nknotsTDNL+pTDNL), byrow = FALSE)
             
             sumTDbase[,k]<-QWR[,c( (dim(data)[2]+(i5+i7)*(nknotsWCE + pWCE+1)+(i2+i4)*(nknotsTDNL+pTDNL)
                        +1): (dim(data)[2]+(i5+i7)*(nknotsWCE + pWCE+1)+(i2+i4+1)*(nknotsTDNL+pTDNL)+1))]%*%modTD$coef[(Nn+(k-1)*(nknotsTDNL+pTDNL+1)+1):(Nn+k*(nknotsTDNL+pTDNL+1))]

             for (k2 in 1:(nknotsTDNL+pTDNL)){ 
               
               kal<- do.call("rbind", lapply(1:length(Id),  # combine the calcluated D_j(u) of each paitents at risk set together
                                             function(j) .wcecalc(listeT[-1], NLbase[data[,1]==Id[j],k2] , 
                                                      data[data[,1]==Id[j],TypeNEW[2]], sumWCEbase3[[k]], cutoff[pos5[k]])))
               kal <- kal[is.na(kal[, 1]) == FALSE, ] # get rid of all the patients that do not belong to any risk set 
               
               
               VV[,(nknotsTDNL+pTDNL)*(k-1)+k2]<-sumTDbase[,k]*kal #calculat E(u) for each of the NL spline coefficient
               
               
               tt2<-paste(tt2,"VV[,",(nknotsTDNL+pTDNL)*(k-1)+k2,"]+",sep="")
             }
           }
          
           if(i6!=0){
             sumWCEbaseWN<-list()
             for(k in 1:i6){#estimate NL effect of variables having NL+WCE effects
               #pos7: positions of variables having only WCE+NL effect among those having 3 effects and WCE+NL to locate the position 
               #of Nbasis and Dbasis
               #pos8: positions of variables having WCE+NL effects among those having WCE to locate the position for cutoff
               if(constrained[pos8[k]]=="FALSE"){
                 sumWCEbaseWN[[k]]<- matrix(unlist(Dbasis[[pos7[k]]]), ncol = nWCEb, byrow = FALSE)%*%modWCE$coef[(Nn+cumindWNT[i8+1]+cumindWN[k]+1):(Nn+ cumindWNT[i8+1]+cumindWN[k]+indWN[k+1])]
               }else if(constrained[[pos8[k]]]=="Right"){
                 sumWCEbaseWN[[k]]<- matrix(unlist(Dbasis[[pos7[k]]]), ncol =nWCEb, byrow = FALSE)[,1:(nWCEb-2)]%*%modWCE$coef[(Nn+cumindWNT[i8+1]+cumindWN[k]+1):(Nn+cumindWNT[i8+1]+cumindWN[k]+indWN[k+1])]
               }else{
                 sumWCEbaseWN[[k]]<- matrix(unlist(Dbasis[[pos7[k]]]), ncol =nWCEb, byrow = FALSE)[,3:nWCEb]%*%modWCE$coef[(Nn+cumindWNT[i8+1]+cumindWN[k]+1):(Nn+cumindWNT[i8+1]+cumindWN[k]+indWN[k+1])]
               }
               
               NLbase<-matrix(unlist(Nbasis[[pos7[k]]]), ncol = (nknotsTDNL+pTDNL), byrow = FALSE)
               
               for (k2 in 1:(nknotsTDNL+pTDNL)){ 
                 
                 
                 kal<- do.call("rbind", lapply(1:length(Id),  # combine the calcluated D_j(u) of each paitents at risk set together
                                               function(j) .wcecalc(listeT[-1], NLbase[data[,1]==Id[j],k2] , 
                                                                    data[data[,1]==Id[j],TypeNEW[2]], sumWCEbaseWN[[k]], cutoff[pos8[k]])))
                 kal <- kal[is.na(kal[, 1]) == FALSE, ] # get rid of all the patients that do not belong to any risk set 
                 
                 
                 VV[,(nknotsTDNL+pTDNL)*(i8+k-1)+k2]<-kal #calculat E(u) for each of the NL spline coefficient
                 
                 
                 tt2<-paste(tt2,"VV[,",(nknotsTDNL+pTDNL)*(i8+k-1)+k2,"]+",sep="")
               }
             }
           }   
           
           if(i4!=0){ # estimate NL effect for variables having NLTD effects
             #pos10 marks the position of variables having NL TD among those having NL no WCE to locate the NL spline basis in QWR
             for(k in 1:i4){ 
               sumTDbase[,(i8+k)]<-QWR[,c( (dim(data)[2]+(i5+i7)*(nknotsWCE + pWCE+1)+(i2+i4)*(nknotsTDNL+pTDNL)
                                       +1): (dim(data)[2]+(i5+i7)*(nknotsWCE + pWCE+1)+(i2+i4+1)*(nknotsTDNL+pTDNL)+1))]%*%modTD$coef[(Nn+(i8+i7+k-1)*(nknotsTDNL+pTDNL+1)+1):(Nn+(i8+i7+k)*(nknotsTDNL+pTDNL+1))]
               
               VV[,(nNLb*(i8+i6+k-1)+1):(nNLb*(i8+i6+k))]<-sumTDbase[,(i8+k)]*QWR[,(dim(data)[2]+(i5+i7)*nWCEb +(pos10[k]-1) *nNLb +1):(dim(data)[2]+(i5+i7)*nWCEb +pos10[k]*nNLb)]
               covp<-paste( "VV[,", (nNLb*(i8+i6+k-1)+1):(nNLb*(i8+i6+k)), "]", sep = "")
               tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))
             }
           } 
           
        
           
           modX<-coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
           vrais<-c(vrais,modX$loglik[2])
           
           ###############################################################################
           # estimate TD effects conditional on the NL and WCE effects
           
           
           modNL<-modX
           tt2<-tt   
           VV<-matrix(ncol=(i8+i4+i7)*(nknotsTDNL+pTDNL+1),nrow=dim(QWR)[1]) 
           
           sumNLbaseW<-matrix(NA,ncol=i8,nrow = length(data[,1]))
           sumNLbaseNW<-matrix(NA,ncol=i4,nrow = length(QWR[,1]))
           
           for (k in 1:i8){ #estimate TD effect of variables with 3 effects
             #pos6: positions of variables having 3 effect among those having 3 effects and WCE+NL to locate the position 
             #of the spline basis for Nbasis and Dbasis
             #pos5: positions of variables having 3 effects among those having WCE to locate the position for cutoff
             sumNLbaseW[,k]<- matrix(unlist(Nbasis[[pos6[k]]]), ncol = (nknotsTDNL+pTDNL),
                                    byrow = FALSE)%*%modNL$coef[(Nn+(k-1)*(nknotsTDNL+pTDNL)+1):(Nn+k*(nknotsTDNL+pTDNL))]
             
              
             kal<- do.call("rbind", lapply(1:length(Id),  #calculate the WCE incorperationg both estimated coefficients for WCE and NL spline
                                           function(j) .wcecalc(listeT[-1], sumNLbaseW[data[,1]==Id[j],k] , 
                                                                data[data[,1]==Id[j],TypeNEW[2]], sumWCEbase3[[k]], cutoff[pos5[k]])))
             kal <- kal[is.na(kal[, 1]) == FALSE, ] # get rid of all the patients that do not belong to any risk set 
             
             
             VV[,((1+nNLb)*(k-1)+1):((1+nNLb)*k)]<-kal*QWR[,(dim(data)[2]+(i5+i7)*nWCEb+(i2+i4)*nNLb+1):(dim(data)[2]+(i5+i7)*nWCEb+(i2+i4+1)*nNLb+1)]
             
             covp<-paste( "VV[,", ((1+nNLb)*(k-1)+1):((1+nNLb)*k), "]", sep = "")
             tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
           }
           
           if(i7!=0){ #estimate TD effect of variables having TD+WCE
             #pos9: positions of variables having only WCE+TD effect among those having WCE no NL to locate the position 
             #of the spline basis for QWR 
             sumWCEbaseWT<-list()
             for(k in 1:i7){
               if(constrained[pos2[k]]=="FALSE"){
                 sumWCEbaseWT[[k]]<- QWR[,c((dim(data)[2]+(pos9[k]-1)*nWCEb+1):(dim(data)[2]+pos9[k]*nWCEb))]%*%modWCE$coef[(Nn+cumindWNT[i8+1]
                                                                                                                             +cumindWN[i6+1]+cumindWT[k]+1):(Nn+cumindWNT[i8+1]+cumindWN[i6+1]+cumindWT[k]+indWT[k+1])]
               }else if(constrained[pos2[k]]=="Right"){
                 sumWCEbaseWT[[k]]<- QWR[,c((dim(data)[2]+(pos9[k]-1)*nWCEb+1):(dim(data)[2]+pos9[k]*nWCEb))][,1:(nWCEb-2)]%*%modWCE$coef[(Nn+cumindWNT[i8+1]
                                                                                                                                           +cumindWN[i6+1]+cumindWT[k]+1):(Nn+cumindWNT[i8+1]+cumindWN[i6+1]+cumindWT[k]+indWT[k+1])]
               }else{
                 sumWCEbaseWT[[k]]<- QWR[,c((dim(data)[2]+(pos9[k]-1)*nWCEb+1):(dim(data)[2]+pos9[k]*nWCEb))][,3:nWCEb]%*%modWCE$coef[(Nn+cumindWNT[i8+1]
                                                                                                                                       +cumindWN[i6+1]+cumindWT[k]+1):(Nn+cumindWNT[i8+1]+cumindWN[i6+1]+cumindWT[k]+indWT[k+1])]
               }
               VV[,((1+nNLb)*(i8+k-1)+1):((1+nNLb)*(i8+k))]<-QWR[,(dim(data)[2]+(i5+i7)*nWCEb+(i2+i4)*nNLb+1):(dim(data)[2]+(i5+i7)*nWCEb+(i2+i4+1)*nNLb+1)
                                                                 ]*as.vector(matrix(unlist(sumWCEbaseWT[[k]]), ncol =1, byrow = FALSE)) #Td spline basis *estimated WCE(u)
               covp<-paste( "VV[,", ((1+nNLb)*(i8+k-1)+1):((1+nNLb)*(i8+k)), "]", sep = "")
               tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
                       
             }
           }
           if(i4!=0){#estimate TD effect of variables having TD+NL
             #pos10 marks the position of variables having NL TD among those having NL no WCE to locate the NL spline basis in QWR
             for(k in 1:i4){
               
               sumNLbaseNW[,k]<- QWR[,c((dim(data)[2]+(i5+i7)*(nknotsWCE + pWCE+1)+(pos10[k]-1)*(nknotsTDNL+pTDNL)+1):(dim(data)[2]
                        +(i5+i7)*(nknotsWCE +pWCE+1)+pos10[k]*(nknotsTDNL+pTDNL)))]%*%modNL$coef[(Nn+(i8+i6+k-1)*(nknotsTDNL+pTDNL)+1):(Nn+(i8+i6+k)*(nknotsTDNL+pTDNL))]
               #TD spline basis * estimated NL effect
               VV[,((1+nNLb)*(i8+i7+k-1)+1):((1+nNLb)*(i8+i7+k))]<-sumNLbaseNW[,k]*QWR[,(dim(data)[2]+(i5+i7)*nWCEb+(i2+i4)*nNLb+1):(dim(data)[2]+(i5+i7)*nWCEb+(i2+i4+1)*nNLb+1)]
               
               covp<-paste( "VV[,", ((1+nNLb)*(i8+i7+k-1)+1):((1+nNLb)*(i8+i7+k)), "]", sep = "")
               tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
             }
           }  
      
          modX<-coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
           vrais<-c(vrais,modX$loglik[2]) # estimate the TD effects conditional on the WCE and NL effects
           
           
           diff<-vrais[length(vrais)]-vrais[length(vrais)-3]
           
            print(diff)
         }
  
    } else if( (i4+i6+i7) !=0){ # having only two of the three effects
      ##################################### 
      ## estimate only the WCE effect first if there is any#
      
      
      if( (i6+i7) !=0){ #if there is WCE effect
        
        if ( i6!= 0) { # if there is both WCE+NL
          pos7<-match(variablesNEW[NL == 1 & TD == 0 & WCE==1], variablesNEW[WCE == 1 & NL==1])
          #pos7: positions of variables having only WCE+NL effect among those having 3 effects and WCE+NL to locate the position 
          #of the spline basis for Nbasis
          #pos8: positions of variables having WCE+NL effects among those having WCE to locate the position for cutoff
          for(k in 1:i6){ # for the covariates having WCE+NL effects, need to combine the WCE basis and the TVC first 
            #in 1st step, assume TVC have linear effects
            
            kal<- do.call("rbind", lapply(1:length(Id),  # combine the calcluated D_j(u) of each paitents at risk set together
                                          function(j) .wcecalc(listeT[-1], data[data[,1]==Id[j], variablesNEW[NL == 1 & TD == 
                                                                   0 & WCE==1][k]], data[data[,1]==Id[j],TypeNEW[2]], Dbasis[[pos7[k]]], cutoff[pos8[k]])))
            kal <- kal[is.na(kal[, 1]) == FALSE, ] # get rid of all the patients that do not belong to any risk set 
            flag<-dim(QWR)[2]+1
            if(constrained[pos8[k]]=="FALSE"){
              indWN[k+1]<-nWCEb
              QWR<-cbind(QWR, kal)
             }else if(constrained[pos8[k]]=="Right"){
              indWN[k+1]<-nWCEb-2
              QWR<-cbind(QWR, kal[,1:(dim(kal)[2]-2)])
            }else{
              indWN[k+1]<-nWCEb-2
              QWR<-cbind(QWR, kal[,3:(dim(kal)[2])])
            }
            covp<-paste("QWR[,", flag: dim(QWR)[2], "]", sep = "" )
            tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))
          }
          cumindWN<-cumsum(indWN)
        }
        
        if(i7!=0){ # if there is WCE+TD
          pos9<-match(variablesNEW[NL == 0 & TD == 1 & WCE==1], variablesNEW[WCE == 1 & NL==0])
           #pos9: positions of variables having only WCE+TD effect among those having WCE no NL to locate the position 
          #of the spline basis for QWR
          #pos2:mark position of variable with WCE+TD among all varialbes having WCE
          for(k in 1:i7){
            if(constrained[pos2[k]]=="FALSE"){
              indWT[k+1]<-nWCEb
              covp<-paste( "QWR[,", (dim(data)[2] + (pos9[k] -1) * nWCEb + 1):(dim(data)[2] + pos9[k] * nWCEb)  , "]", sep = "")
           }else if(constrained[pos2[k]]=="Right"){
              indWT[k+1]<-nWCEb-2
              covp<-paste( "QWR[,", (dim(data)[2] + (pos9[k] -1) * nWCEb + 1):(dim(data)[2] + pos9[k] * nWCEb-2)  , "]", sep = "")
            }else{
              indWT[k+1]<-nWCEb-2
              covp<-paste( "QWR[,", (dim(data)[2] + (pos9[k] -1) * nWCEb + 3):(dim(data)[2] + pos9[k] * nWCEb)  , "]", sep = "")
            }
            tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
          }
          cumindWT<-cumsum(indWT) 
        }
        
        modX <- coxph(eval(parse(text = substr(tt2, 13, nchar(tt2) - 
                                                 1))), method = "efron")
        vrais <- c(vrais, modX$loglik[2])
        
        if(i6!=0){ # NLWCE
          for (j in 1:i6){
            if(constrained[pos8[j]]=="FALSE" ){
              rescale<-as.vector(Dbasis[[(i)]]%*%modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                                                            + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])])
              modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                         + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])] <- modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                                                                                                                                     + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])]/sum(rescale)
            }else if(constrained[pos8[j]]=="Right"){
              rescale<-as.vector(Dbasis[[(i)]]%*%c(modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                                                              + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])],0,0))
              modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                         + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])] <- modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                                                                                                                                     + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])]/sum(rescale)
            }else{
              rescale<-as.vector(Dbasis[[(i)]]%*%c(0,0, modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                                                                   + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])]))
              modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                         + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])]<-modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                                                                                                                                   + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])]/sum(rescale)
            }
            res[i5+j]<-sum(rescale)
          }                     
        }
        
        if(i7!=0){
          for(j in 1:i7){
            if(constrained[pos9[j]]=="FALSE" ){
              rescale<-as.vector(Bbasis[[(i5+i)]]%*%modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2*nNLb+i3*(nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                                                                + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2*nNLb+i3* (nNLb+1) +cumindWN[i6+1]+cumindWT[j]+indWT[j+1])])
              modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                         + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[i6+1]+cumindWT[j]+indWT[j+1])] <- modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                                                                                                                                                    + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[i6+1]+cumindWT[j]+indWT[j+1])]/sum(rescale)
            }else if(constrained[pos9[j]]=="Right"){
              rescale<-as.vector( Bbasis[[(i5+i)]]%*%c(modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                                                                  + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[i6+1]+cumindWT[j]+indWT[j+1])],0,0))
              modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                         + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[i6+1]+cumindWT[j]+indWT[j+1])] <- modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                                                                                                                                     + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[i6+1]+cumindWT[j]+indWT[j+1])]/sum(rescale)
            }else{
              rescale<-as.vector( Bbasis[[(i5+i)]]%*%c(0,0, modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                                                                       + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]+indWT[j+1])]))
              modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                         + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[i6+1]+cumindWT[j]+indWT[j+1])]<-modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                                                                                                                                                  + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[i6+1]+cumindWT[j]+indWT[j+1])]/sum(rescale)
            }
            res[i5+i6+j]<-sum(rescale)
          }
        }
        
        # print(modX$coefficients)
        # print(modX$coef)
        
        modWCE<-modX
        tt2<-tt
      } 
      
      
      ##########################################
      ## estimate NL effect conditional on the TD and WCE 
      if((i6+i4)!=0){ #if there is a NL effect
        VV<-matrix(ncol=(i4+i6)*(nknotsTDNL+pTDNL),nrow=dim(QWR)[1]) #number of colum is the number of variables having NL
        
        if(i6!=0){
          sumWCEbaseWN<-list()
         for(k in 1:i6){#estimate NL effect of variables having NL+WCE effects
          #pos7: positions of variables having only WCE+NL effect among those having 3 effects and WCE+NL to locate the position 
          #of Nbasis and Dbasis
          #pos8: positions of variables having WCE+NL effects among those having WCE to locate the position for cutoff
           if(constrained[pos8[k]]=="FALSE"){
             sumWCEbaseWN[[k]]<- matrix(unlist(Dbasis[[pos7[k]]]), ncol = nWCEb, byrow = FALSE)%*%modWCE$coef[(Nn+cumindWN[k]+1):(Nn+cumindWN[k]+indWN[k+1])]
           }else if(constrained[[pos8[k]]]=="Right"){
             sumWCEbaseWN[[k]]<- matrix(unlist(Dbasis[[pos7[k]]]), ncol =nWCEb, byrow = FALSE)[,1:(nWCEb-2)]%*%modWCE$coef[(Nn+cumindWN[k]+1):(Nn+cumindWN[k]+indWN[k+1])]
           }else{
             sumWCEbaseWN[[k]]<- matrix(unlist(Dbasis[[pos7[k]]]), ncol =nWCEb, byrow = FALSE)[,3:nWCEb]%*%modWCE$coef[(Nn+cumindWN[k]+1):(Nn+cumindWN[k]+indWN[k+1])]
           }
             
          NLbase<-matrix(unlist(Nbasis[[pos7[k]]]), ncol = (nknotsTDNL+pTDNL), byrow = FALSE)
          
          for (k2 in 1:(nknotsTDNL+pTDNL)){ 
            
            
            kal<- do.call("rbind", lapply(1:length(Id),  # combine the calcluated D_j(u) of each paitents at risk set together
                                          function(j) .wcecalc(listeT[-1], NLbase[data[,1]==Id[j],k2] , 
                                                               data[data[,1]==Id[j],TypeNEW[2]], sumWCEbaseWN[[k]], cutoff[pos8[k]])))
            kal <- kal[is.na(kal[, 1]) == FALSE, ] # get rid of all the patients that do not belong to any risk set 
            
            
            VV[,(nknotsTDNL+pTDNL)*(k-1)+k2]<-kal #calculat E(u) for each of the NL spline coefficient
            
            
            tt2<-paste(tt2,"VV[,",(nknotsTDNL+pTDNL)*(k-1)+k2,"]+",sep="")
          }
        }
        }
        
        if(i4!=0){
          pos10<-match(variablesNEW[NL == 1 & TD == 1 & WCE==0], variablesNEW[WCE == 0&NL==1])
          #pos10 marks the position of variables having NL TD among those having NL no WCE to locate the NL spline basis in QWR
          for(k in 1:i4){ 
            VV[,(nNLb*(i6+k-1)+1):(nNLb*(i6+k))]<-QWR[,(dim(data)[2]+(i5+i7)*nWCEb +(pos10[k]-1) *nNLb +1):(dim(data)[2]+(i5+i7)*nWCEb +pos10[k]*nNLb)]
            covp<-paste( "VV[,", (nNLb*(i6+k-1)+1):(nNLb*(i6+k)), "]", sep = "")
            tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
          }
        }
        
        
        modX<-coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
        vrais<-c(vrais,modX$loglik[2]) 
        
        modNL<-modX
        tt2<-tt
        }
        
      ############################
      #estimate TD conditional on WCE NL
      
      if( (i4+i7) !=0){ #if there is a TD effect
        VV<-matrix(ncol=(i4+i7)*(nknotsTDNL+pTDNL+1),nrow=dim(QWR)[1]) 
        
        sumNLbase<-matrix(NA,ncol=i4,nrow = length(QWR[,1]))
        
        if(i7!=0){
          sumWCEbaseWT<-list()
          #pos9: positions of variables having only WCE+TD effect among those having WCE no NL to locate the position 
          #of the spline basis for QWR 
          for(k in 1:i7){
            if(constrained[pos2[k]]=="FALSE"){
              sumWCEbaseWT[[k]]<- QWR[,c((dim(data)[2]+(pos9[k]-1)*nWCEb+1):(dim(data)[2]+pos9[k]*nWCEb))]%*%modWCE$coef[(Nn
                                                                        +cumindWN[i6+1]+cumindWT[k]+1):(Nn+cumindWN[i6+1]+cumindWT[k]+indWT[k+1])]
            }else if(constrained[pos2[k]]=="Right"){
              sumWCEbaseWT[[k]]<- QWR[,c((dim(data)[2]+(pos9[k]-1)*nWCEb+1):(dim(data)[2]+pos9[k]*nWCEb))][,1:(nWCEb-2)]%*%modWCE$coef[(Nn
                                                                  +cumindWN[i6+1]+cumindWT[k]+1):(Nn+cumindWN[i6+1]+cumindWT[k]+indWT[k+1])]
            }else{
              sumWCEbaseWT[[k]]<- QWR[,c((dim(data)[2]+(pos9[k]-1)*nWCEb+1):(dim(data)[2]+pos9[k]*nWCEb))][,3:nWCEb]%*%modWCE$coef[(Nn
                                                                  +cumindWN[i6+1]+cumindWT[k]+1):(Nn+cumindWN[i6+1]+cumindWT[k]+indWT[k+1])]
            }
            VV[,((1+nNLb)*(k-1)+1):((1+nNLb)*k)]<-QWR[,(dim(data)[2]+(i5+i7)*nWCEb+(i2+i4)*nNLb+1):(dim(data)[2]+(i5+i7)*nWCEb+(i2+i4+1)*nNLb+1)
                                                              ]*as.vector(matrix(unlist(sumWCEbaseWT[[k]]), ncol =1, byrow = FALSE))
            covp<-paste( "VV[,", ((1+nNLb)*(k-1)+1):((1+nNLb)*k), "]", sep = "")
            tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))
          }
        } 
        
        if(i4!=0){
          #pos10 marks the position of variables having NL TD among those having NL no WCE to locate the NL spline basis in QWR
          for(k in 1:i4){
            sumNLbase[,k]<- QWR[,c((dim(data)[2]+(i5+i7)*nWCEb+(pos10[k]-1)*nNLb+1):(dim(data)[2]
                                              +(i5+i7)*nWCEb+pos10[k]*nNLb))]%*%modNL$coef[(Nn+(i6+k-1)*nNLb+1):(Nn+(i6+k)*nNLb)]
            VV[,((1+nNLb)*(i7+k-1)+1):((1+nNLb)*(i7+k))]<-sumNLbase[,k]*QWR[,(dim(data)[2]+(i5+i7)*nWCEb+(i2+i4)*nNLb+1):(dim(data)[2]+(i5+i7)*nWCEb+(i2+i4+1)*nNLb+1)]
            
            covp<-paste( "VV[,", ((1+nNLb)*(i7+k-1)+1):((1+nNLb)*(i7+k)), "]", sep = "")
            tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
          } 
        }
        
        modX<-coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
        vrais<-c(vrais,modX$loglik[2]) # estimate the TEL effects conditional on the TD and NL effects
        
        modTD<-modX
        tt2<-tt
      }
      
      diff<-1
      #w<-1
      #w<-0
      while(diff>0.01){ ## ACE algorithm until converge 
        #############################
        #estimate WCE effect conditional on TD, NL effects
       
        if( (i6+i7) !=0){ #if there is a WCE effect 
          VV<-matrix(ncol=cumindWN[i6+1]+cumindWT[i7+1],nrow=dim(QWR)[1]) 
        
         if(i6!=0){
            sumNLbase<-matrix(NA,ncol=i6,nrow = length(data[,1]))
          for(k in 1:i6){ # for the covariates having WCE+NL effects, need to combine the WCE basis and the NL basis first 
            # pNLTDWCE indicate the position of the varialbes having WCE+NL+TD among the variables with NL 
            # calculate the Dj(u) with WCE basis and NL basis
            ## set the initial coefficents for NL basis to 1 
            sumNLbase[,k]<- matrix(unlist(Nbasis[[pos7[k]]]), ncol = (nknotsTDNL+pTDNL)
                                        , byrow = FALSE)%*%modNL$coef[(Nn+(k-1)*(nknotsTDNL+pTDNL)+1):(Nn+k*(nknotsTDNL+pTDNL))]
            
            kal<- do.call("rbind", lapply(1:length(Id),  # combine the calcluated D_j(u) of each paitents at risk set together
                                          function(j) .wcecalc(listeT[-1],  sumNLbase[data[,1]==Id[j],k] , 
                                                               data[data[,1]==Id[j],TypeNEW[2]], Dbasis[[pos7[k]]], cutoff[pos8[k]])))
            kal <- kal[is.na(kal[, 1]) == FALSE, ] # get rid of all the patients that do not belong to any risk set 
            
            
            if(constrained[pos8[k]]=="FALSE"){
              VV[,(cumindWN[k]+1):(cumindWN[k]+indWN[k+1])]<-kal
              covp<-paste( "VV[,",(cumindWN[k]+1):(cumindWN[k]+indWN[k+1]), "]", sep = "")
            }else if(constrained[pos8[k]]=="Right"){
              VV[,(cumindWN[k]+1):(cumindWN[k]+indWN[k+1])]<-kal[,1:(nWCEb-2)]
              covp<-paste( "VV[,",(cumindWN[k]+1):(cumindWN[k]+indWN[k+1]), "]", sep = "")
             }else{
              VV[,(cumindWN[k]+1):(cumindWN[k]+indWN[k+1])]<-sumTDbase[,k]*kal[,3:nWCEb]
              covp<-paste( "VV[,",(cumindWN[k]+1):(cumindWN[k]+indWN[k+1]), "]", sep = "")
              }
            tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
          }
         }  
          
          if(i7!=0){ 
            sumTDbase<-matrix(NA,ncol=i7,nrow = length(QWR[,1]))
          #pos9: positions of variables having only WCE+TD effect among those having WCE no NL to locate the position 
          #of the spline basis for QWR
          for(k in 1:i7){
            sumTDbase[,k]<- QWR[,c((dim(data)[2]+(i5+i7)*(nknotsWCE + pWCE+1)+(i2+i4)*(nknotsTDNL+pTDNL)
                                         +1):(dim(data)[2]+(i5+i7)*(nknotsWCE + pWCE+1)+(i2+i4+1)*(nknotsTDNL+pTDNL)+1))]%*%modTD$coef[(Nn+(k-1)*(nknotsTDNL+pTDNL+1)+1):(Nn+k*(nknotsTDNL+pTDNL+1))]
            
          if(constrained[pos2[k]]=="FALSE"){
              VV[,(cumindWN[i6+1]+cumindWT[k]+1):(cumindWN[i6+1]+cumindWT[k]+indWT[k+1])]<-sumTDbase[
                ,k]*QWR[, (dim(data)[2] + (pos9[k] -1)*nWCEb +1):(dim(data)[2] +pos9[k]*nWCEb)]
               }else if(constrained[pos2[k]]=="Right"){
              VV[,(cumindWN[i6+1]+cumindWT[k]+1):(cumindWN[i6+1]+cumindWT[k]+indWT[k+1])]<-sumTDbase[
                ,k]*QWR[,(dim(data)[2] + (pos9[k] -1)*nWCEb +1):(dim(data)[2] +pos9[k]*nWCEb-2)]
            }else{
              VV[,(cumindWN[i6+1]+cumindWT[k]+1):(cumindWN[i6+1]+cumindWT[k]+indWT[k+1])]<-sumTDbase[
                ,k]*QWR[,(dim(data)[2] + (pos9[k] -1)*nWCEb +3):(dim(data)[2] +pos9[k]*nWCEb)]
           }
            covp<-paste( "VV[,",(cumindWN[i6+1]+cumindWT[k]+1):(cumindWN[i6+1]+cumindWT[k]+indWT[k+1]), "]", sep = "")
            tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
           } #WCE basis multiply the estimated value of TD effects
          }

          modX<-coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
          vrais<-c(vrais,modX$loglik[2])
          
          if(i6!=0){ # NLWCE
            for (j in 1:i6){
              if(constrained[pos8[j]]=="FALSE" ){
                rescale<-as.vector(Dbasis[[(i)]]%*%modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                                                              + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])])
                modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                           + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])] <- modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                                                                                                                                       + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])]/sum(rescale)
              }else if(constrained[pos8[j]]=="Right"){
                rescale<-as.vector(Dbasis[[(i)]]%*%c(modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                                                                + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])],0,0))
                modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                           + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])] <- modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                                                                                                                                       + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])]/sum(rescale)
              }else{
                rescale<-as.vector(Dbasis[[(i)]]%*%c(0,0, modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                                                                     + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])]))
                modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                           + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])]<-modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                                                                                                                                     + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])]/sum(rescale)
              }
              res[i5+j]<-sum(rescale)
            }                     
          }
          
          
          if(i7!=0){
            for(j in 1:i7){
              if(constrained[pos9[j]]=="FALSE" ){
                rescale<-as.vector(Bbasis[[(i5+i)]]%*%modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2*nNLb+i3*(nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                                                                 + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2*nNLb+i3* (nNLb+1) +cumindWN[i6+1]+cumindWT[j]+indWT[j+1])])
                modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                           + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[i6+1]+cumindWT[j]+indWT[j+1])] <- modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                                                                                                                                                      + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[i6+1]+cumindWT[j]+indWT[j+1])]/sum(rescale)
              }else if(constrained[pos9[j]]=="Right"){
                rescale<-as.vector( Bbasis[[(i5+i)]]%*%c(modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                                                                    + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[i6+1]+cumindWT[j]+indWT[j+1])],0,0))
                modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                           + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[i6+1]+cumindWT[j]+indWT[j+1])] <- modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                                                                                                                                                      + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[i6+1]+cumindWT[j]+indWT[j+1])]/sum(rescale)
              }else{
                rescale<-as.vector( Bbasis[[(i5+i)]]%*%c(0,0, modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                                                                         + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]+indWT[j+1])]))
                modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                           + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[i6+1]+cumindWT[j]+indWT[j+1])]<-modX$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[i6+1]+cumindWT[j]
                                                                                                                                                    + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[i6+1]+cumindWT[j]+indWT[j+1])]/sum(rescale)
              }
              res[i5+i6+j]<-sum(rescale)
            }
          }
          
          # print(modX$coefficients)
          # print(modX$coef)
          modWCE<-modX
          tt2<-tt
        }
        ##########################################################
        # estimate WCE conditional on NL,TD
        
        if( (i6+i4) !=0){ #if there is a NL effect
          
          VV<-matrix(ncol=(i4+i6)*(nknotsTDNL+pTDNL),nrow=dim(QWR)[1])
          
          
          if(i6!=0){
            sumWCEbaseWN<-list()
            for(k in 1:i6){#estimate NL effect of variables having NL+WCE effects
              #pos7: positions of variables having only WCE+NL effect among those having 3 effects and WCE+NL to locate the position 
              #of Nbasis and Dbasis
              #pos8: positions of variables having WCE+NL effects among those having WCE to locate the position for cutoff
              if(constrained[pos8[k]]=="FALSE"){
                sumWCEbaseWN[[k]]<- matrix(unlist(Dbasis[[pos7[k]]]), ncol = nWCEb, byrow = FALSE)%*%modWCE$coef[(Nn+cumindWN[k]+1):(Nn+ cumindWN[k]+indWN[k+1])]
              }else if(constrained[[pos8[k]]]=="Right"){
                sumWCEbaseWN[[k]]<- matrix(unlist(Dbasis[[pos7[k]]]), ncol =nWCEb, byrow = FALSE)[,1:(nWCEb-2)]%*%modWCE$coef[(Nn+cumindWN[k]+1):(Nn+cumindWN[k]+indWN[k+1])]
              }else{
                sumWCEbaseWN[[k]]<- matrix(unlist(Dbasis[[pos7[k]]]), ncol =nWCEb, byrow = FALSE)[,3:nWCEb]%*%modWCE$coef[(Nn+cumindWN[k]+1):(Nn+cumindWN[k]+indWN[k+1])]
              }
            
              NLbase<-matrix(unlist(Nbasis[[pos7[k]]]), ncol = (nknotsTDNL+pTDNL), byrow = FALSE)
              
              for (k2 in 1:(nknotsTDNL+pTDNL)){ 
                
                
                kal<- do.call("rbind", lapply(1:length(Id),  # combine the calcluated D_j(u) of each paitents at risk set together
                                              function(j) .wcecalc(listeT[-1], NLbase[data[,1]==Id[j],k2] , 
                                                                   data[data[,1]==Id[j],TypeNEW[2]], sumWCEbaseWN[[k]], cutoff[pos8[k]])))
                kal <- kal[is.na(kal[, 1]) == FALSE, ] # get rid of all the patients that do not belong to any risk set 
                
                
                VV[,(nknotsTDNL+pTDNL)*(k-1)+k2]<-kal #calculat E(u) for each of the NL spline coefficient
                
                
                tt2<-paste(tt2,"VV[,",(nknotsTDNL+pTDNL)*(k-1)+k2,"]+",sep="")
              }
            }
          }  
          
          if(i4!=0){ # estimate NL effect for variables having NLTD effects
            #pos10 marks the position of variables having NL TD among those having NL no WCE to locate the NL spline basis in QWR
            sumTDbase<-matrix(NA,ncol=i4,nrow = length(QWR[,1]))
            for(k in 1:i4){ 
              sumTDbase[,k]<-QWR[,c( (dim(data)[2]+(i5+i7)*(nknotsWCE + pWCE+1)+(i2+i4)*(nknotsTDNL+pTDNL)
                                           +1): (dim(data)[2]+(i5+i7)*(nknotsWCE + pWCE+1)+(i2+i4+1)*(nknotsTDNL+pTDNL)+1))]%*%modTD$coef[(Nn+(i7+k-1)*(nknotsTDNL+pTDNL+1)+1):(Nn+(i7+k)*(nknotsTDNL+pTDNL+1))]
              
              VV[,(nNLb*(i6+k-1)+1):(nNLb*(i6+k))]<-sumTDbase[,k]*QWR[,(dim(data)[2]+(i5+i7)*nWCEb +(pos10[k]-1) *nNLb +1):(dim(data)[2]+(i5+i7)*nWCEb +pos10[k]*nNLb)]
              covp<-paste( "VV[,", (nNLb*(i6+k-1)+1):(nNLb*(i6+k)), "]", sep = "")
              tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))
            }
          } 
          modX<-coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
          vrais<-c(vrais,modX$loglik[2])
          
          modNL<-modX
          tt2<-tt 
          
        }
        
        ###############################################################################
        # estimate TD effects conditional on the NL and WCE effects (NLTEL+TDTEL)
        
        if( (i7+i4) !=0){ #if there is a TD effect
          
          VV<-matrix(ncol=(i4+i7)*(nknotsTDNL+pTDNL+1),nrow=dim(QWR)[1]) 
         
          if(i7!=0){ #estimate TD effect of variables having TD+WCE
            #pos9: positions of variables having only WCE+TD effect among those having WCE no NL to locate the position 
            #of the spline basis for QWR 
            sumWCEbaseWT<-list()
            for(k in 1:i7){
              if(constrained[pos2[k]]=="FALSE"){
                sumWCEbaseWT[[k]]<- QWR[,c((dim(data)[2]+(pos9[k]-1)*nWCEb+1):(dim(data)[2]+pos9[k]*nWCEb))]%*%modWCE$coef[(Nn
                                                                                        +cumindWN[i6+1]+cumindWT[k]+1):(Nn+cumindWN[i6+1]+cumindWT[k]+indWT[k+1])]
              }else if(constrained[pos2[k]]=="Right"){
                sumWCEbaseWT[[k]]<- QWR[,c((dim(data)[2]+(pos9[k]-1)*nWCEb+1):(dim(data)[2]+pos9[k]*nWCEb))][,1:(nWCEb-2)]%*%modWCE$coef[(Nn
                                                                                +cumindWN[i6+1]+cumindWT[k]+1):(Nn+cumindWN[i6+1]+cumindWT[k]+indWT[k+1])]
              }else{
                sumWCEbaseWT[[k]]<- QWR[,c((dim(data)[2]+(pos9[k]-1)*nWCEb+1):(dim(data)[2]+pos9[k]*nWCEb))][,3:nWCEb]%*%modWCE$coef[(Nn
                                                                            +cumindWN[i6+1]+cumindWT[k]+1):(Nn+cumindWN[i6+1]+cumindWT[k]+indWT[k+1])]
              }
              VV[,((1+nNLb)*(k-1)+1):((1+nNLb)*k)]<-QWR[,(dim(data)[2]+(i5+i7)*nWCEb+(i2+i4)*nNLb+1):(dim(data)[2]+(i5+i7)*nWCEb+(i2+i4+1)*nNLb+1)
                                                                ]*as.vector(matrix(unlist(sumWCEbaseWT[[k]]), ncol =1, byrow = FALSE)) #Td spline basis *estimated WCE(u)
              covp<-paste( "VV[,", ((1+nNLb)*(k-1)+1):((1+nNLb)*k), "]", sep = "")
              tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
           }
          }
          
          if(i4!=0){#estimate TD effect of variables having TD+NL
            #pos10 marks the position of variables having NL TD among those having NL no WCE to locate the NL spline basis in QWR
            sumNLbase<-matrix(NA,ncol=i4,nrow = length(QWR[,1]))
             for(k in 1:i4){
              sumNLbase[,k]<- QWR[,c((dim(data)[2]+(i5+i7)*(nknotsWCE + pWCE+1)+(pos10[k]-1)*(nknotsTDNL+pTDNL)+1):(dim(data)[2]
                                        +(i5+i7)*(nknotsWCE +pWCE+1)+pos10[k]*(nknotsTDNL+pTDNL)))]%*%modNL$coef[(Nn+(i6+k-1)*(nknotsTDNL+pTDNL)+1):(Nn+(i6+k)*(nknotsTDNL+pTDNL))]
            
              VV[,((1+nNLb)*(i7+k-1)+1):((1+nNLb)*(i7+k))]<-sumNLbase[,k]*QWR[,(dim(data)[2]+(i5+i7)*nWCEb+(i2+i4)*nNLb+1):(dim(data)[2]+(i5+i7)*nWCEb+(i2+i4+1)*nNLb+1)]
              
              covp<-paste( "VV[,", ((1+nNLb)*(i7+k-1)+1):((1+nNLb)*(i7+k)), "]", sep = "")
              tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
            }
          } 
          
          
          modX<-coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
          vrais<-c(vrais,modX$loglik[2]) 
          
          modTD<-modX #### check if this line is needed 
          tt2<-tt 
        }
        ### since there will only be two of the effects, so compare every two models
        ### there could be cases where one variable has WCE+TD and the other variable has NL+TD
        ### in this case, there were three effects estimated, ie, three models 
        if((i4!=0&i6!=0)|(i4!=0&i7!=0)|(i6!=0&i7!=0)){
          diff<-vrais[length(vrais)]-vrais[length(vrais)-3]}
        else{ 
          diff<-vrais[length(vrais)]-vrais[length(vrais)-2]
        }
         print(diff)
        
      }
      
      
      
    }    else {
      modX <- coxph(eval(parse(text = substr(tt2, 13, nchar(tt2) - 
                                               1))), method = "efron")
      vrais <- c(modX$loglik[2])
    } 
    rm(QWR)
    
    
    MAT <- matrix(ncol=V,nrow = 2*(1 + nNLb) +nWCEb) # store the estimated parameter 
    if (i1 != 0) { #if no NL and TD no WCE
      for (j in 1:i1) { # first row store estimated coefficients for those having no NL TD WCE effects
        MAT[1, j] <- modX$coef[j]
      }
    }
    if (i5 != 0) {# only WCE no TD no NL
      for (j in 1:i5) { # 2nd to nknotsWCE+pWCE+2 row stores the nknotsWCE+pWCE+1 coefficients for thoes only having WCE
        # effects, each column represents one variable
       
        if(constrained[pos1[j]]=="FALSE"){
          
          rescale<-as.vector(Bbasis[[i]]%*%modX$coef[(i1 + cumindW[j]+1):(i1+cumindW[j]+indW[j+1])])
          MAT[1, j] <- sum(rescale)
          MAT[2:(nWCEb+1), i1 + j] <- modX$coef[(i1 + cumindW[j]+1):(i1+cumindW[j]+indW[j+1])]/sum(rescale)
          
        }else if(constrained[pos1[j]]=="Right"){
          
          rescale<-as.vector(Bbasis[[i]]%*%c(modX$coef[(i1 + cumindW[j]+1):(i1+cumindW[j]+indW[j+1])],0,0))
          MAT[1, j] <-sum(rescale)
          MAT[2:(nWCEb+1), i1 + j] <- c(modX$coef[(i1 + cumindW[j]+1):(i1+cumindW[j]+indW[j+1])],0,0)/sum(rescale)
        
        }else{
          rescale<-as.vector(Bbasis[[i]]%*%c(0,0, modX$coef[(i1 + cumindW[j]+1):(i1+cumindW[j]+indW[j+1])]))
          MAT[1, j] <-sum(rescale)
          MAT[2:(nWCEb+1), i1 + j] <- c(0,0, modX$coef[(i1 + cumindW[j]+1):(i1+cumindW[j]+indW[j+1])])/sum(rescale)
        }
       
      }
    }
    if (i2 != 0) { # only NL no TD no WCE
      for (j in 1:i2) { # row nWCEb+1+1 to nWCEb+1+nNLb stores TD effects
        MAT[(nWCEb+1+1):(nWCEb+1+nNLb), i1 + i5 + j] <- modX$coef[(i1 + 
                    cumindW[i5+1] + (j - 1) * nNLb + 1):(i1 +cumindW[i5+1] + j * nNLb )]
        
      }
    }
    
    if (i3 != 0) { # only TD no NL no WCE
      for (j in 1:i3) { # row nWCEb+nNLb+1+1 to nWCEb+2*nNLb+1+1 stores TD effects
        MAT[(nWCEb+nNLb+2):(nWCEb+2*nNLb+ 2), i1 + i5 + i2+j] <- modX$coef[(i1 + 
                 cumindW[i5+1] + i2 * nNLb +(j - 1) * (nNLb+1)+ 1):(i1 + cumindW[i5+1] + i2 * nNLb+ j * (nNLb+ 1) )]
      }
    }
    
    if(i8!=0){ #NL+TD+WCE
      for (j in 1:i8) {
        MAT[1, i1 +i5 +i2+i3 + j]<-res[i5+j] 
        if(constrained[pos5[j]]=="FALSE"){
          
        MAT[2:(nWCEb+1), i1 +i5 +i2+i3 + j] <- modWCE$coef[(i1 + cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) + cumindWNT[j]
                                  + 1):(i1 + cumindW[i5+1]+ i2 * nNLb +i3* (nNLb+1)+cumindWNT[j]+indWNT[j+1])]
        }else if(constrained[pos5[j]]=="Right"){
        MAT[2:(nWCEb+1), i1 +i5 +i2+i3 + j] <- c(modWCE$coef[(i1 + cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) + cumindWNT[j]
                                           + 1):(i1 + cumindW[i5+1]+ i2 * nNLb +i3* (nNLb+1)+cumindWNT[j]+indWNT[j+1])],0,0)
        }else{
        MAT[2:(nWCEb+1), i1 +i5 +i2+i3 + j] <- c(0,0,modWCE$coef[(i1 + cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) + cumindWNT[j]
                                                              + 1):(i1 + cumindW[i5+1]+ i2 * nNLb +i3* (nNLb+1)+cumindWNT[j]+indWNT[j+1])]) 
        }
        MAT[(nWCEb+2):(nWCEb+1+nNLb), i1 + i2 + i3 + i5+j] <- modNL$coef[(i1 + cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) +(j-1)*nNLb 
                           + 1):(i1 + cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) + j*nNLb )]
        MAT[(nWCEb+nNLb+2):(nWCEb+2*(nNLb)+2), i1 + i2 + i3 + i5+
               j] <- modTD$coef[(i1 +cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) + (j-1)*(nNLb+1) 
                      + 1):(i1 +cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) +j*(nNLb+1))]
        
      }
    }
    if(i6!=0){ # NLWCE
      for (j in 1:i6){
        MAT[1, i1 + i2 + i3 + i5 + i8+j]<-res[i5+i8+j]
        if(constrained[pos8[j]]=="FALSE" ){
        MAT[2:(nWCEb+ 1), i1 + i2 + i3 + i5 + i8+j] <- modWCE$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                                  + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])]
        }else if(constrained[pos8[j]]=="Right"){
        MAT[2:(nWCEb+ 1), i1 + i2 + i3 + i5 + i8+j] <- c(modWCE$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                                 + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])],0,0)
          
        }else{
          MAT[2:(nWCEb+ 1), i1 + i2 + i3 + i5 + i8+j] <-c(0,0, modWCE$coef[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                                       + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])])
          
        }
        MAT[(nWCEb+2):(nWCEb+1+nNLb), i1 + i2 + i3 + i5+ i8+ j] <- modNL$coef[(i1 + cumindW[i5+1] + (i2+i8) * nNLb+i3* (nNLb+1)
                                 + (j - 1) * nNLb + 1):(i1 + cumindW[i5+1] + (i2+i8) * nNLb+i3* (nNLb+1) + j*nNLb )]
        
      }                     
    }
    if(i7!=0){ # WCETD
      for(j in 1:i7){
        MAT[1, i1 + i2 + i3 + i5 + i8+i6+j]<-res[i5+i8+i6+j]
        if(constrained[pos2[j]]=="FALSE" ){
        MAT[2:(nWCEb + 1), i1 + i2 + i3 + i5 + i8+i6+j] <- modWCE$coef[(i1 +cumindW[i5+1]+cumindWNT[i8+1]+cumindWN[i6+1]+i2 * nNLb+i3* (nNLb+1)
                                    + cumindWT[j] + 1):(i1 +cumindW[i5+1]+cumindWNT[i8+1]+cumindWN[i6+1] + i2 * nNLb+ i3* (nNLb+1) + cumindWT[j]+indWT[j+1])]
        }else if(constrained[pos2[j]]=="Right"){
         MAT[2:(nWCEb + 1), i1 + i2 + i3 + i5 + i8+i6+j] <- c(modWCE$coef[(i1 +cumindW[i5+1]+cumindWNT[i8+1]+cumindWN[i6+1]+i2 * nNLb+i3* (nNLb+1)
                                + cumindWT[j] + 1):(i1 +cumindW[i5+1]+cumindWNT[i8+1]+cumindWN[i6+1] + i2 * nNLb+ i3* (nNLb+1) + cumindWT[j]+indWT[j+1])],0,0)
        }else{
          MAT[2:(nWCEb + 1), i1 + i2 + i3 + i5 + i8+i6+j] <- c(NA,NA,modWCE$coef[(i1 +cumindW[i5+1]+cumindWNT[i8+1]+cumindWN[i6+1]+i2 * nNLb+i3* (nNLb+1)
                                                                    + cumindWT[j] + 1):(i1 +cumindW[i5+1]+cumindWNT[i8+1]+cumindWN[i6+1] + i2 * nNLb+ i3* (nNLb+1) + cumindWT[j]+indWT[j+1])])
       }
        MAT[(nWCEb+nNLb+2):(nWCEb+2*nNLb+2), i1 + i2 + i3 + i5+ i8+i6+j] <-modTD$coef[(i1 + cumindW[i5+1]+ i2 * nNLb+(i3+i8)* (nNLb+1)+ (j - 1) * (nNLb+ 1)
                                + 1):(i1 +cumindW[i5+1] + i2 * nNLb+(i3+i8)* (nNLb+1) + j * (nNLb + 1))]
      }
    }
    if(i4!=0){ #TDNL
      for(j in 1:i4){
        MAT[(nWCEb+ 2):(nWCEb+1+nNLb), i1 + i2 + i3 + i5 + i8+i6+i7+j] <- modNL$coef[(i1 + cumindW[i5+1] + (i2+i8+i6) * nNLb+i3* (nNLb+1)
                           + (j - 1) * nNLb + 1):(i1 + cumindW[i5+1] + (i2+i8+i6) * nNLb+i3* (nNLb+1) + j * nNLb )]
        MAT[(nWCEb+nNLb + 2):(nWCEb+2*nNLb+2), i1 + i2 + i3 + i5+ i8+i6+i7+j] <- modTD$coef[(i1 +cumindW[i5+1] + i2 * nNLb+(i3+i8+i7)* (nNLb+1)
                          + (j - 1) * (nNLb + 1) + 1):(i1 +cumindW[i5+1] + i2 * nNLb+(i3+i8+i7)* (nNLb+1) + j * (nNLb + 1))]
        
      }
    }
    
   ### modified stops here, run the model first ##### 
    
    MATse <- matrix(ncol = V, nrow = 2*(1 + nNLb) +nWCEb ) # store the estimated SE
    if (i1 != 0) {
      for (j in 1:i1) {
        MATse[1, j] <- sqrt(diag(modX$var)[j])
      }
    }
    if (i5 != 0) {
      for (k in 1:i5) {
        if(constrained[pos1[k]]=="FALSE"){
          MATse[2:(nWCEb+1), i1 + k] <- sqrt(diag(modX$var)[(i1 + cumindW[k]+1):(i1+cumindW[k]+indW[k+1])])
        }else if(constrained[pos1[k]]=="Right"){
          MATse[2:(nWCEb+1), i1 + k] <- c(sqrt(diag(modX$var)[(i1 + cumindW[k]+1):(i1+cumindW[k]+indW[k+1])]),NA,NA)
        }else{
          MATse[2:(nWCEb+1), i1 + k] <- c(NA,NA, sqrt(diag(modX$var)[(i1 + cumindW[k]+1):(i1+cumindW[k]+indW[k+1])]))
        }
      }
    }
    if (i2 != 0) {
      for (j in 1:i2) {
        MATse[(nWCEb+ 2):(nWCEb+nNLb+ 1), i1 + i5 + j] <- sqrt(diag(modX$var)[(i1 + 
                       cumindW[i5+1] + (j - 1) * nNLb + 1):(i1 + cumindW[i5+1] + j * nNLb )])
      }
    }
    
    if (i3 != 0) { # only TD no NL no WCE
      for (j in 1:i3) { # row m+p+2 to 2*(m+p+1) stores TD effects
        MATse[(nWCEb+nNLb+2):(nWCEb+2*nNLb+ 2), i1 + i2 + i5+j] <-sqrt(diag(modX$var)[(i1 + 
               cumindW[i5+1] + i2 * nNLb +(j - 1) * (nNLb+1)+ 1):(i1 +  cumindW[i5+1] + i2 * nNLb+ j * (nNLb+ 1) )])
      }
    }
      
    if(i8!=0){ #NL+TD+WCE
      for (j in 1:i8) {
        
        if(constrained[pos5[j]]=="FALSE"){
          MATse[2:(nWCEb+1), i1 +i5 +i2+i3 + j] <- sqrt(diag(modWCE$var)[(i1 + cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) + cumindWNT[j]
                                                              + 1):(i1 + cumindW[i5+1]+ i2 * nNLb +i3* (nNLb+1)+cumindWNT[j]+indWNT[j+1])])
        }else if(constrained[pos5[j]]=="Right"){
          MATse[2:(nWCEb+1), i1 +i5 +i2+i3 + j] <- c(sqrt(diag(modWCE$var)[(i1 + cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) + cumindWNT[j]
                                                                + 1):(i1 + cumindW[i5+1]+ i2 * nNLb +i3* (nNLb+1)+cumindWNT[j]+indWNT[j+1])]),NA,NA)
        }else{
          MATse[2:(nWCEb+1), i1 +i5 +i2+i3 + j] <- c(NA,NA,sqrt(diag(modWCE$var)[(i1 + cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) + cumindWNT[j]
                                                                    + 1):(i1 + cumindW[i5+1]+ i2 * nNLb +i3* (nNLb+1)+cumindWNT[j]+indWNT[j+1])])) 
        }
        
        MATse[(nWCEb+2):(nWCEb+1+nNLb), i1 + i2 + i3 + i5+j] <-sqrt(diag(modNL$var)[(i1 + cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) +(j-1)*nNLb 
                                                                          + 1):(i1 + cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) + j*nNLb )])
        MATse[(nWCEb+nNLb+2):(nWCEb+2*(nNLb)+2), i1 + i2 + i3 + i5+
              j] <- sqrt(diag(modX$var)[(i1 +cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) + (j-1)*(nNLb+1) 
                               + 1):(i1 +cumindW[i5+1] + i2 * nNLb +i3* (nNLb+1) +j*(nNLb+1))])
        }
    }
    if(i6!=0){ # NLWCE
      for (j in 1:i6){
        if(constrained[pos8[j]]=="FALSE" ){
          MATse[2:(nWCEb+ 1), i1 + i2 + i3 + i5 + i8+j] <-sqrt(diag(modWCE$var)[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                                                                      + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])])
        }else if(constrained[pos8[j]]=="Right"){
          MATse[2:(nWCEb+ 1), i1 + i2 + i3 + i5 + i8+j] <- c(sqrt(diag(modWCE$var)[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                                                                        + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])]),NA,NA)
          
        }else{
          MATse[2:(nWCEb+ 1), i1 + i2 + i3 + i5 + i8+j] <-c(NA,NA,sqrt(diag(modWCE$var)[(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) + cumindWN[j] 
                                                                            + 1):(i1 + cumindW[i5+1]+cumindWNT[i8+1] + i2 * nNLb+i3* (nNLb+1) +  cumindWN[j]+indWN[j+1])]))
          
        }
        MATse[(nWCEb+2):(nWCEb+1+nNLb), i1 + i2 + i3 + i5+ i8+ j] <-sqrt(diag(modNL$var)[(i1 + cumindW[i5+1] + (i2+i8) * nNLb+i3* (nNLb+1)
                                                                               + (j - 1) * nNLb + 1):(i1 + cumindW[i5+1] + (i2+i8) * nNLb+i3* (nNLb+1) + j*nNLb )])
        
       }                     
    }
    
   if(i7!=0){ # WCETD
      for(j in 1:i7){
        if(constrained[pos2[j]]=="FALSE" ){
          MATse[2:(nWCEb + 1), i1 + i2 + i3 + i5 + i8+i6+j] <- sqrt(diag(modWCE$var)[(i1 +cumindW[i5+1]+cumindWNT[i8+1]+cumindWN[i6+1]+i2 * nNLb+i3* (nNLb+1)
                                                                          + cumindWT[j] + 1):(i1 +cumindW[i5+1]+cumindWNT[i8+1]+cumindWN[i6+1] + i2 * nNLb+ i3* (nNLb+1) + cumindWT[j]+indWT[j+1])])
        }else if(constrained[pos2[j]]=="Right"){
          MATse[2:(nWCEb + 1), i1 + i2 + i3 + i5 + i8+i6+j] <- c(sqrt(diag(modWCE$var)[(i1 +cumindW[i5+1]+cumindWNT[i8+1]+cumindWN[i6+1]+i2 * nNLb+i3* (nNLb+1)
                                                                            + cumindWT[j] + 1):(i1 +cumindW[i5+1]+cumindWNT[i8+1]+cumindWN[i6+1] + i2 * nNLb+ i3* (nNLb+1) + cumindWT[j]+indWT[j+1])]),NA,NA)
          
        }else{
          MATse[2:(nWCEb + 1), i1 + i2 + i3 + i5 + i8+i6+j] <- c(NA,NA,sqrt(diag(modWCE$var)[(i1 +cumindW[i5+1]+cumindWNT[i8+1]+cumindWN[i6+1]+i2 * nNLb+i3* (nNLb+1)
                                                                                        + cumindWT[j] + 1):(i1 +cumindW[i5+1]+cumindWNT[i8+1]+cumindWN[i6+1] + i2 * nNLb+ i3* (nNLb+1) + cumindWT[j]+indWT[j+1])]))
       }
        
       MATse[(nWCEb+nNLb + 2):(nWCEb+2*nNLb+2), i1 + i2 + i3 + i5+ i8+i6+ j] <-sqrt(diag(modX$var)[(i1 +cumindW[i5+1] + i2 * nNLb+ (i3+i8) * (nNLb+ 1)
                          + (j - 1) * (nNLb + 1) + 1):(i1 +cumindW[i5+1] + i2 * nNLb+ (i3+i8) * (nNLb+ 1) + j * (nNLb + 1) )])
      }
    }
    if(i4!=0){ #TDNL
      for(j in 1:i4){
        MATse[(nWCEb + 2):(nWCEb+nNLb+ 1), i1 + i2 + i3 + i5 + i8+i7+i6+j] <-sqrt(diag(modNL$var)[(i1 + 
                     cumindW[i5+1] + (i2+i8+i6) * nNLb+ i3 * (nNLb+ 1)+ (j - 1) * nNLb + 1):(i1 + cumindW[i5+1] + (i2+i8+i6) * nNLb+ i3 * (nNLb+ 1) 
                                                                             + j * nNLb )])
        MATse[(nWCEb+nNLb+ 2):(nWCEb+2*nNLb+2), i1 + i2 + i3 + i5+ i8+i6+i7+
                j] <-sqrt(diag(modX$var)[(i1 + cumindW[i5+1]+ i2 * nNLb+ (i3+i8+i7) * (nNLb+ 1)+ (j - 1) * (nNLb + 1) + 1):(i1  
                                               + cumindW[i5+1] + i2 * nNLb+ (i3+i8+i7) * (nNLb+ 1) + j * (nNLb + 1) )])
        
      }
    }
    
    
    var_order <- c(variables[NL == 0 & TD == 0 &WCE==0], variables[NL == 
                      0 & TD == 0 & WCE==1], variables[NL == 1 & TD == 0&WCE==0],variables[NL == 0 & TD == 1&WCE==0],
                   variables[NL == 1 & TD == 1&WCE==1],variables[NL ==1 & TD == 0 & WCE==1],
                   variables[NL == 0 & TD == 1&WCE==1], variables[NL == 1 & TD == 1&WCE==0])
    
    
    coefficients<-MAT[1,] 
    names(coefficients)<-variables[match(var_order, variables)]
    
    se_coef<-MATse[1,] 
    names(se_coef)<-variables[match(var_order, variables)]
    
    coefficients_splines_WCE<-as.matrix(MAT[2:(nWCEb+1),])

    colnames(coefficients_splines_WCE)<-variables[match(var_order, variables)]
    
    coefficients_splines_NL<-as.matrix(MAT[(nWCEb+2):(nWCEb+nNLb+1),])
    coefficients_splines_NL<-rbind(rep(0,V),coefficients_splines_NL)
    # coefficients_splines_NL[1,(NL==0)]<-NA
     colnames(coefficients_splines_NL)<-variables[match(var_order, variables)]
    
    coefficients_splines_TD<-as.matrix(MAT[(nWCEb+nNLb+2):(nWCEb+2*nNLb+2),])
    colnames(coefficients_splines_TD)<-variables[match(var_order, variables)]
    
    knots_WCE<-knotsNEW[c(1:sum(WCE))]
    knots_covariates<-knotsNEW[c((sum(WCE)+1):(sum(WCE)+V))[match(var_order, variables)]]
 

    if (V>1) {knots_covariates[ match(variables,var_order)[NL==0]]<-NA} else { if (NL==0) {knots_covariates<-NA} }
    knots_time<-knotsNEW[[sum(WCE)+V+1]]
   

    gnu.c <-list() ### store the estimated weights at each time point over the past cutoff window for each variable
    if(i1!=0){
      for(i in 1:i1){
        gnu.c[[paste("WCE weights for variable",variables[c(NL==0&TD==0&WCE==0)][i], sep = " ")]]<-NA
      }
    }
    if(i5!=0){
      for(i in 1:i5){
      gnu.c[[paste("WCE weights for variable",variables[c(NL==0&TD==0&WCE==1)][i], sep = " ")]]<-as.vector(Bbasis[[i]]%*%coefficients_splines_WCE[,i1+i])
      }
    }
    if(i2!=0){
      for(i in 1:i2){
        gnu.c[[paste("WCE weights for variable",variables[c(NL==1&TD==0&WCE==0)][i], sep = " ")]]<-NA
      }
    }
    if(i3!=0){
      for(i in 1:i3){
        gnu.c[[paste("WCE weights for variable",variables[c(NL==0&TD==1&WCE==0)][i], sep = " ")]]<-NA
      }
    }
    if(i8!=0){
      for(i in 1:i8){
        gnu.c[[paste("WCE weights for variable",variables[c(NL==1&TD==1&WCE==1)][i], sep = " ")]]<-as.vector(Dbasis[[(i6+i)]]%*%coefficients_splines_WCE[,i1+i5+i2+i3+i])
      }
    }
    if(i6!=0){
    for(i in 1:i6){
      gnu.c[[paste("WCE weights for variable",variables[c(NL==1&TD==0&WCE==1)][i], sep = " ")]]<-as.vector(Dbasis[[(i)]]%*%coefficients_splines_WCE[,(i1+i5+i2+i3+i8+i)])
    }
    }
    if(i7!=0){
      for(i in 1:i7){
        gnu.c[[paste("WCE weights for variable",variables[c(NL==0&TD==1&WCE==1)][i], sep = " ")]]<-as.vector(Bbasis[[i5+i]]%*%coefficients_splines_WCE[,(i1+i5+i2+i3+i8+i6+i)])
      }
    }
    if(i4!=0){
      for(i in 1:i4){
        gnu.c[[paste("WCE weights for variable",variables[c(NL==1&TD==1&WCE==0)][i], sep = " ")]]<-NA
      }
    }
  

    n.events<-sum(data[,TypeNEW[3]]==1)
    PL.c<-vrais[length(vrais)]
    npara<-i1 + cumindWNT[i8+1]+cumindW[i5+1]+cumindWN[i6+1]+cumindWT[i7+1]+ (i2+i4+i6+i8) * nNLb + (i3+i4+i7+i8) * (nNLb + 1)
    
    BIC.c <- .my.bic(PL.c, n.events , npara , aic)
   
    rm(data)
    
    list(Partial_Log_Likelihood =PL.c,
          aic=aic,info.criterion=BIC.c,Number_of_parameters =npara, Number_events=n.events, 
         WCEconstrained=constrained,
         variables=variables[match(var_order, variables)],WCEmat = gnu.c,
         knots_WCE = knots_WCE,
         knots_covariates = knots_covariates, 
         knots_time = knots_time, 
         WCE=WCE[match(var_order, variables)],
         TD=TD[match(var_order, variables)],
         NL=NL[match(var_order, variables)],
         coefficients = coefficients, Standard_Error=se_coef,
         coefficients_splines_WCE = coefficients_splines_WCE,
         coefficients_splines_NL = coefficients_splines_NL,
         coefficients_splines_TD = coefficients_splines_TD)
    
         
  }
  
}

      
      
      
      
      
    





