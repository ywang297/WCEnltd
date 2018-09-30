sparseWCE<-function (data, Type,variables,WCE,TD,NL,weights,cutoff, nknotsTDNL,pTDNL,
                       aic = FALSE){ 
  
  

  V <- length(variables) # number of total covariates
  
  variablesNEW<-match(variables,names(data))  # stores the postition of each variable in the data set 
  TypeNEW<-match(Type,names(data)) # stored the position of the Type varialbes in the data set
  
  # Only need this program when there are both WCE effects and NL effects. If the NL effect is rejected for certain 
  # variable, then only need to used the estimated weight combined with the LOCF to get the wce, then use Cox PH 
  # to estimate this metric, the estimated coefficient for linear effect with WCE indicate the magnitude of this effect
    

  
    i1 <- sum((NL+TD+WCE) == 0)    ##number of non-nl non-td non-WCE variables
    i2 <- sum(((NL == 1) & (TD == 0) & (WCE==0))) ##number of NL non-td non-WCE variables
    i3 <- sum(((NL == 0) & (TD == 1) & (WCE==0)))  ##number of non-NL td non-WCE variables
    i4 <- sum((NL + TD) == 2 & (WCE==0)) ###number of NL TD non-WCE variables
    i5 <- sum(((NL == 0) & (TD == 0) & (WCE==1))) ##number of non-NL non-td WCE variables
    i6 <- sum(((NL == 1) & (TD == 0) & (WCE==1))) ##number of NL non-td WCE variables
    i7 <- sum(((NL == 0) & (TD == 1) & (WCE==1))) ##number of non-NL td WCE variables
    i8 <- sum(((NL == 1) & (TD == 1) & (WCE==1))) ##number of NL td WCE variables
    
    nNLb<-nknotsTDNL+pTDNL #number of NL basis
   
    if(length(cutoff)!=sum(WCE))
      stop("ERROR:length of cutoff must equal to the number of specified WCE effects")
    maxTime <- max(data[,TypeNEW[2]])# maximum FUP time 
    if (max(cutoff) > maxTime) ## cutoff can be a vector with length equal to sum(WCE) 
      stop("ERROR: cutoff must be smaller than the longest follow-up time")
     #### do not allow user specified knots for now, only allow default knots ###
    
    listeprobaquantile <- seq(1, nknotsTDNL)/(nknotsTDNL + 1) ## place the interior knots according to quantiles  
    knotsNEW <- list() ## store both the interior and exterior knots for each model 
    
    ## No need to store the knots for WCE, will use weights directly ##
 
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
    
    X <- split(data, data[, 1]) ##lists store full data for each individual respectively
    matX <- sapply(X, DvlpMatrix, listeT = listeT, ncol = ncol, TypeNEW=TypeNEW) ## transform the data structure
    QWR <- do.call(rbind, matX) #combine the lists to a full data matrix, get rid of observation not in any risk set
    
    ##Always estimate WCE first, create bspline base for each specified WCE effects###
    ## bspline basis are different when the variable only have WCE and have both WCE and TD
    ### the data still have to be in 1 line per day format since to estiamte the TD effect, the data need to be transformed to this format first
    
     Id <- unique(data[,1])
    #i5=number of varialbes only have WCE effect, i7=no. of varialbes have both WCE+TD
    #.wcecalc() calculated the D_j(u) when applied to the spline basis for weights and TVC
    ## calculated the wce when applied to the weights and the TVC directly
    ## No need to use spline base for weights in this function 
    # weights is a list store the estimated weights only for those having WCE effect,so similar useage as cutoff
    
     
     var_order <- c(variables[NL == 0 & TD == 0 &WCE==0], variables[NL == 
                   0 & TD == 0 & WCE==1], variables[NL == 1 & TD == 0&WCE==0],variables[NL == 0 & TD == 1&WCE==0],
                    variables[NL == 1 & TD == 1&WCE==1],variables[NL ==1 & TD == 0 & WCE==1],
                    variables[NL == 0 & TD == 1&WCE==1], variables[NL == 1 & TD == 1&WCE==0])
     
     
    if(i5!=0){#calculate the weighted cumulative TVC for variables having only WCE effect
     pos12<-match(variablesNEW[NL == 0 & TD == 0 & WCE==1], variablesNEW[WCE == 1])
      pos1<-match(variables[NL == 0 & TD == 0 & WCE==1],var_order)
      #pos1:mark position of variable only with WCE among all varialbes having WCE
      for(i in 1:i5){
        kal <- do.call("rbind", lapply(1:length(Id),  # combine the calcluated D_j(u) of each paitents at risk set together
                                       function(j) .wcecalc(listeT[-1], data[data[,1]==Id[j], variablesNEW[NL == 0 &TD==0 & WCE==1][i]], 
                                                            data[data[,1]==Id[j],TypeNEW[2]], as.matrix(weights[[pos1[i]]]), cutoff[pos12[i]])))
        kal <- kal[is.na(kal[, 1]) == FALSE, ] # get rid of all the patients that do not belong to any risk set 
        QWR<-cbind(QWR, kal)
      }
    }
    if(i7!=0){## calculate the weighted cumulative TVC for variables having WCE+TD effect
      pos22<-match(variablesNEW[NL == 0 & TD == 1 & WCE==1], variablesNEW[WCE == 1])
      pos2<-match(variables[NL == 0 & TD == 1 & WCE==1],var_order)
     #pos2:mark position of variable with WCE+TD among all varialbes having WCE
      for(i in 1:i7){
          kal <- do.call("rbind", lapply(1:length(Id),  # combine the calcluated D_j(u) of each paitents at risk set together
                                       function(j) .wcecalc(listeT[-1], data[data[,1]==Id[j], variablesNEW[NL == 0 & TD==1& WCE==1][i]], 
                                                            data[data[,1]==Id[j],TypeNEW[2]], as.matrix(weights[[pos2[i]]]), cutoff[pos22[i]])))
        kal <- kal[is.na(kal[, 1]) == FALSE, ] # get rid of all the patients that do not belong to any risk set 
        QWR<-cbind(QWR, kal)
      }
    }
    
    # generate the spline values for NL effect not having WCE, since the cumulative sum of the spline basis are not needed, the spline
    ## base can be generated based on QWR where the observation not from risk set are removed 
    if(i2!=0){
      for (i in 1:i2){
        QWR <- cbind(QWR, splineDesign(knotsNEW[[seq(1, V, 1)[NL == 1&WCE==0&TD==0][i]]], x=QWR[, variablesNEW[NL == 1&WCE==0&TD==0][i]],ord=pTDNL+1)[,-1] )    
      }
    }
    if(i4!=0){
      for (i in 1:i4){
        QWR <- cbind(QWR, splineDesign(knotsNEW[[seq(1, V, 1)[NL == 1&WCE==0&TD==1][i]]], x=QWR[, variablesNEW[NL == 1&WCE==0&TD==1][i]],ord=pTDNL+1)[,-1] )    
      }
    }
    
    # spline value for time, pTDNL+nknotTDNL+1 
    if((i3+i4+i7+i8)!=0){
      QWR <- cbind(QWR, splineDesign(knotsNEW[[V+1]], x=QWR[, TypeNEW[2]],ord=pTDNL+1))    
      # check the basis for time again, it is generated after the observation not in the risk set are removed?
    }
    
    #After this step, the order of variables in QWR: variables, (i5+i7)+(i2+i4)*(nNLb)+nNLb+1
    
    
    # for the variables having both WCE and NL effect, need to generate NL basis based on the original data instead of the reduced data
    Nbasis<-list() #include NL basis for NL+WCE and NL+WCE+TD
    if(i6!=0){#spline basis for NL combined with the estimated weights, i6:WCE+NL
      pos82<-match(variablesNEW[NL == 1 & TD == 0 & WCE==1], variablesNEW[WCE == 1])
      pos8<-match(variables[NL == 1 & TD == 0 & WCE==1],var_order)
      
      #pos8: positions of variables having WCE+NL effects among those having WCE to locate the position for cutoff
      for(i in 1:i6){
         Nbasis[[i]]<-splineDesign(knotsNEW[[seq(1, V, 1)[NL == 1&WCE==1&TD==0][i]]], x=data[, variablesNEW[NL == 1&WCE==1&TD==0][i]],ord=pTDNL+1)[,-1]
         
         NLbase<-matrix(unlist(Nbasis[[i]]), ncol = (nknotsTDNL+pTDNL), byrow = FALSE)
         for (k2 in 1:(nknotsTDNL+pTDNL)){ 
           
          kal<- do.call("rbind", lapply(1:length(Id),  # combine each of the NL basis with the weights
                                         function(j) .wcecalc(listeT[-1], NLbase[data[,1]==Id[j],k2] , 
                                                              data[data[,1]==Id[j],TypeNEW[2]], as.matrix(weights[[pos8[i]]]), cutoff[pos82[i]])))
           kal <- kal[is.na(kal[, 1]) == FALSE, ] # get rid of all the patients that do not belong to any risk set 
          
           QWR<-cbind(QWR, kal)
         
          }
       }
    }
    if(i8!=0){## spline basis for NL combined with weights for variable having WCE+NL+TD
      pos52<-match(variablesNEW[NL == 1 & TD == 1 & WCE==1], variablesNEW[WCE == 1])
      pos5<-match(variables[NL == 1 & TD == 1 & WCE==1],var_order)
      
       #pos5: positions of variables having 3 effects among those having WCE to locate the position for cutoff
      for(i in 1:i8){
         Nbasis[[(i6+i)]]<-splineDesign(knotsNEW[[seq(1,V, 1)[NL == 1&WCE==1&TD==1][i]]], x=data[, variablesNEW[NL == 1&WCE==1&TD==1][i]],ord=pTDNL+1)[,-1]
         NLbase<-matrix(unlist(Nbasis[[i6+i]]), ncol = (nknotsTDNL+pTDNL), byrow = FALSE)
         for (k2 in 1:(nknotsTDNL+pTDNL)){ 
           
           kal<- do.call("rbind", lapply(1:length(Id),  # combine each of the NL basis with the weights
                                         function(j) .wcecalc(listeT[-1], NLbase[data[,1]==Id[j],k2] , 
                                                              data[data[,1]==Id[j],TypeNEW[2]], as.matrix(weights[[pos5[i]]]), cutoff[pos52[i]])))
           kal <- kal[is.na(kal[, 1]) == FALSE, ] # get rid of all the patients that do not belong to any risk set 
           
           QWR<-cbind(QWR, kal)
           
         }
      }
    }  
    
    
    
    tt <- paste("modX<-coxph(Surv(QWR[,", TypeNEW[1], "],QWR[,", 
                TypeNEW[2], "],QWR[,", TypeNEW[3], "])~", sep = "")
    
    
    ############## if there is none or only 1 special effect #####
    Nn <- 0 
    if (i1 != 0) { # if there are variables have neither NL nor TD nor WCE effects 
      for (k in 1:i1) {
        tt <- paste(tt, "QWR[,", variablesNEW[NL == 0 & TD == 
                                                0 & WCE==0][k], "]+", sep = "")
        Nn <- Nn + 1
      }
    }##Nn=number (NL=TD=WCE=0)*i1
    

    if (i5 != 0) { #only WCE
      for (k in 1:i5) {
        tt <- paste(tt, "QWR[,", dim(data)[2]+k, "]+", sep = "")
        Nn <- Nn + 1
      }
    }  
  
    
    if (i2 != 0) { #only NL
      for (k in 1:i2) { # all the WCE bases are created, no need to use cumind to index
        covp<-paste( "QWR[,", (dim(data)[2] + (i5+i7) + (k-1) * nNLb + 1):(dim(data)[2] + (i5+i7) + k * nNLb)  , "]", sep = "")
        tt<-paste(tt, paste(c(covp,""), collapse= "+") )
        Nn<-Nn+nNLb
      }
    }##Nn=number (NL=TD=WCE=0)*i1 +(TD=NL=0,WCE=1)*i5+(TD=WCE=0, NL=1)*i2
    
    if (i3 != 0) {# only TD effect
      for (k in 1:i3) {
        flag<-dim(QWR)[2]+1
        QWR <- cbind(QWR, QWR[, variablesNEW[NL == 0 & TD ==1& WCE==0][k]] * QWR[, (dim(data)[2] + (i5+i7)+(i2+i4)*nNLb 
                                                                                    +1): (dim(data)[2] + (i5+i7)+(i2+i4+1)*nNLb+1)]) 
        covp<-paste("QWR[,", flag: dim(QWR)[2], "]", sep = "" )
        tt<-paste(tt, paste(c(covp,""), collapse= "+"))  #### only the TD spline base combined with covariate values 
        # will enter the model 
        Nn <- Nn + nNLb+1
      }##Nn=number (NL=TD=WCE=0)*i1 +(TD=NL=0,WCE=1)*i5+(TD=WCE=0, NL=1)*i2 +(NL=WCE=0,TD=1)*i3
    } #Nn mark the number of parameters having no or only 1 special effect
    
    
    ### having more than one special effects, in this case, for variables having WCE+NL, and WCE+TD, 
    ## no ACE is needed either, just need to estimate the spline coefficients directly
    ## the varialbes need ACE is the one having WCE+TD+NL, and NL+TD 
    
    if(i6!=0){
      for(k in 1:i6){ #NL basis for WCE+NL starts at dim(data)[2]+(i5+i7)+(i2+i4)*nNLb+i3*(nNLb+1)
        covp<-paste( "QWR[,", (dim(data)[2] + (i5+i7) + (i2+i4+i3)*nNLb + i3+(k-1)*nNLb+1 ):(dim(data)[2] + (i5+i7) +(i2+i4+i3+k)* nNLb+i3)  , "]", sep = "")
        tt<-paste(tt, paste(c(covp,""), collapse= "+") )
        Nn<-Nn+nNLb
      }
    }
    
   if(i7!=0){##variables having WCE+TD effect, multiply the weigted TVC with TD spline basis
     for(k in 1:i7){
       flag<-dim(QWR)[2]+1
       QWR <- cbind(QWR, QWR[,dim(data)[2]+i5+k] * QWR[, (dim(data)[2] + (i5+i7)+(i2+i4)*nNLb 
                                                                                   +1): (dim(data)[2] + (i5+i7)+(i2+i4+1)*nNLb+1)]) 
       covp<-paste("QWR[,", flag: dim(QWR)[2], "]", sep = "" )
       tt<-paste(tt, paste(c(covp,""), collapse= "+"))  #### only the TD spline base combined with covariate values 
       # will enter the model 
       Nn <- Nn + nNLb+1
     }
   }
    
   tt2 <- tt ## tt mark the effects that do not need ACE
    
    
    ### for variables having WCE+NL+TD and NL+TD effects, only 2 steps ACE iterating TD+NL effect is needed in both cases###
    
    vrais <- c()
    

    ## starting with estimating NL effect first  #####  
 
    if((i4+i8)!=0){ # having either WCE+NL+TD, or TD+NL

            ##########################################
      ## estimate NL effect conditional on the TD 
      VV<-matrix(ncol=(i4+i8)*(nknotsTDNL+pTDNL),nrow=dim(QWR)[1]) #number of colum is the number of variables having NL
        
     if(i4!=0){ #estimate NL splines for NL+TD effect
       ##QWR order: dim(data)[2]+(i5+i7)+(i2+i4)nNLb+(nNLb+1)+(i6+i8)nNLb+i7*nNLb   
       for(k in 1:i4){ 
            VV[,(nNLb*(k-1)+1):(nNLb*(k))]<-QWR[,(dim(data)[2]+(i5+i7) +(i2+k-1) *nNLb +1):(dim(data)[2]+(i5+i7) +(i2+k)*nNLb)]
            covp<-paste( "VV[,", (nNLb*(k-1)+1):(nNLb*(k)), "]", sep = "")
            tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
          }
     }
      if(i8!=0){
      for (k in 1:i8){ #estimate NL effect of variables having 3 effects 
        VV[,(nNLb*(i4+k-1)+1):(nNLb*(i4+k))]<-QWR[,(dim(data)[2]+(i5+i7)+(i2+i4+1+i6+k-1) *nNLb +1+1):(dim(data)[2]+(i5+i7)+(i2+i4+1+i6+k)*nNLb+1)]
        covp<-paste( "VV[,", (nNLb*(i4+k-1)+1):(nNLb*(i4+k)), "]", sep = "")
        tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
           }
      }  
        
        modX<-coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
        vrais<-c(vrais,modX$loglik[2]) 
        
        modNL<-modX
        tt2<-tt
      
      
      ############################
      #estimate TD conditional on  NL
      
      VV<-matrix(ncol=(i4+i8)*(nknotsTDNL+pTDNL+1),nrow=dim(QWR)[1]) # store TD basis for each variable
        
        sumNLbase<-matrix(NA,ncol=(i4+i8),nrow = length(QWR[,1]))
    
     if(i4!=0){
      
          for(k in 1:i4){
            sumNLbase[,k]<- QWR[,c((dim(data)[2]+(i5+i7)+(i2+k-1)*nNLb+1):(dim(data)[2]
                                                                              +(i5+i7)+(i2+k)*nNLb))]%*%modNL$coef[(Nn+(k-1)*nNLb+1):(Nn+(k)*nNLb)]
            VV[,((1+nNLb)*(k-1)+1):((1+nNLb)*(k))]<-sumNLbase[,k]*QWR[,(dim(data)[2]+(i5+i7)+(i2+i4)*nNLb+1):(dim(data)[2]+(i5+i7)+(i2+i4+1)*nNLb+1)]
            
            covp<-paste( "VV[,", ((1+nNLb)*(k-1)+1):((1+nNLb)*(k)), "]", sep = "")
            tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
          } 
     }
        
    if(i8!=0){
         
        for(k in 1:i8){
            sumNLbase[,(i4+k)]<- QWR[,c((dim(data)[2]+(i5+i7)+(i2+i4+1+i6+k-1)*nNLb+1+1):(dim(data)[2]
                                        +(i5+i7)+(i2+i4+1+i6+k)*nNLb+1))]%*%modNL$coef[(Nn+(i4+k-1)*nNLb+1):(Nn+(i4+k)*nNLb)]
            VV[,((1+nNLb)*(i4+k-1)+1):((1+nNLb)*(i4+k))]<-sumNLbase[,(i4+k)]*QWR[,(dim(data)[2]+(i5+i7)+(i2+i4)*nNLb+1):(dim(data)[2]+(i5+i7)+(i2+i4+1)*nNLb+1)]
            
            covp<-paste( "VV[,", ((1+nNLb)*(i4+k-1)+1):((1+nNLb)*(i4+k)), "]", sep = "")
            tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
          } 
        }
        
        modX<-coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
        vrais<-c(vrais,modX$loglik[2]) # estimate the TEL effects conditional on the TD and NL effects
        
        modTD<-modX
        tt2<-tt
    
      
      diff<-1
      #w<-1
      #w<-0
      while(diff>0.01){ ## ACE algorithm until converge 
         
          VV<-matrix(ncol=(i4+i8)*(nknotsTDNL+pTDNL),nrow=dim(QWR)[1])
          sumTDbase<-matrix(NA,ncol=(i4+i8),nrow = length(QWR[,1]))## combine TD spline basis with estimated TD coefficient
          
          
          if(i4!=0){ # estimate NL effect for variables having NLTD effects
       
           for(k in 1:i4){ 
              sumTDbase[,k]<-QWR[,c( (dim(data)[2]+(i5+i7)+(i2+i4)*nNLb
                                      +1): (dim(data)[2]+(i5+i7)+(i2+i4+1)*nNLb+1))]%*%modTD$coef[(Nn+(k-1)*(nNLb+1)+1):(Nn+(k)*(nNLb+1))]
              
              VV[,(nNLb*(k-1)+1):(nNLb*(k))]<-sumTDbase[,k]*QWR[,(dim(data)[2]+(i5+i7)+(i2+k-1) *nNLb +1):(dim(data)[2]+(i5+i7) +(i2+k)*nNLb)]
              covp<-paste( "VV[,", (nNLb*(k-1)+1):(nNLb*(k)), "]", sep = "")
              tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))
            }
          } 
          
          
          if(i8!=0){ # estimate NL effect for variables having WCENLTD effects
            for(k in 1:i8){ 
              sumTDbase[,(i4+k)]<-QWR[,c( (dim(data)[2]+(i5+i7)+(i2+i4)*nNLb
                                      +1): (dim(data)[2]+(i5+i7)+(i2+i4+1)*nNLb+1))]%*%modTD$coef[(Nn+(i4+k-1)*(nNLb+1)+1):(Nn+(i4+k)*(nNLb+1))]
              
              VV[,(nNLb*(i4+k-1)+1):(nNLb*(i4+k))]<-sumTDbase[,(i4+k)]*QWR[,(dim(data)[2]+(i5+i7)+(i2+i4+1+i6+k-1) *nNLb +1+1):(dim(data)[2]+(i5+i7)+(i2+i4+1+i6+k)*nNLb+1)]
              covp<-paste( "VV[,", (nNLb*(i4+k-1)+1):(nNLb*(i4+k)), "]", sep = "")
              tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))
            }
          }
          
          modX<-coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
          vrais<-c(vrais,modX$loglik[2])
          
          modNL<-modX
          tt2<-tt 
          
  
        
        ###############################################################################
        # estimate TD effects conditional on the NL and WCE effects (NLTEL+TDTEL)
        

       VV<-matrix(ncol=(i4+i8)*(nknotsTDNL+pTDNL+1),nrow=dim(QWR)[1]) # store TD basis for each variable
          
          sumNLbase<-matrix(NA,ncol=(i4+i8),nrow = length(QWR[,1]))
          
          if(i4!=0){
            
            for(k in 1:i4){
              sumNLbase[,k]<- QWR[,c((dim(data)[2]+(i5+i7)+(i2+k-1)*nNLb+1):(dim(data)[2]
                                                                             +(i5+i7)+(i2+k)*nNLb))]%*%modNL$coef[(Nn+(k-1)*nNLb+1):(Nn+(k)*nNLb)]
              VV[,((1+nNLb)*(k-1)+1):((1+nNLb)*(k))]<-sumNLbase[,k]*QWR[,(dim(data)[2]+(i5+i7)+(i2+i4)*nNLb+1):(dim(data)[2]+(i5+i7)+(i2+i4+1)*nNLb+1)]
              
              covp<-paste( "VV[,", ((1+nNLb)*(k-1)+1):((1+nNLb)*(k)), "]", sep = "")
              tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
            } 
          }
          
          if(i8!=0){
            
            for(k in 1:i8){
              sumNLbase[,(i4+k)]<- QWR[,c((dim(data)[2]+(i5+i7)+(i2+i4+1+i6+k-1)*nNLb+1+1):(dim(data)[2]
                                                                                            +(i5+i7)+(i2+i4+1+i6+k)*nNLb+1))]%*%modNL$coef[(Nn+(i4+k-1)*nNLb+1):(Nn+(i4+k)*nNLb)]
              VV[,((1+nNLb)*(i4+k-1)+1):((1+nNLb)*(i4+k))]<-sumNLbase[,(i4+k)]*QWR[,(dim(data)[2]+(i5+i7)+(i2+i4)*nNLb+1):(dim(data)[2]+(i5+i7)+(i2+i4+1)*nNLb+1)]
              
              covp<-paste( "VV[,", ((1+nNLb)*(i4+k-1)+1):((1+nNLb)*(i4+k)), "]", sep = "")
              tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
            } 
          }
          
          modX<-coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
          vrais<-c(vrais,modX$loglik[2]) 
          
          modTD<-modX
          tt2<-tt        
  
          diff<-abs(vrais[length(vrais)]-vrais[length(vrais)-2])
  
        print(diff)
      }
      
      
      
    }    else {
      modX <- coxph(eval(parse(text = substr(tt2, 13, nchar(tt2) - 
                                               1))), method = "efron")
      vrais <- c(modX$loglik[2])
    } 
    rm(QWR)
    
    ###modX order: i1+i5+i2*nNLb+i3*(nNLb+1)+i6*(nNLb)+i7*(nNLb+1)+ace coefficients
    
    MAT <- matrix(ncol=V,nrow = 2*(1 + nNLb)) # store the estimated parameter 
    if (i1 != 0) { #if no NL and TD no WCE
      for (j in 1:i1) { # first row store estimated coefficients for those having no NL TD WCE effects 
        MAT[1, j] <- modX$coef[j]
      }
    }
    if (i5 != 0) {# first row also stores for estimates of only WCE effect
      for (j in 1:i5) {
        MAT[1,i1+j] <- modX$coef[i1+j]
              }
    }
    if (i2 != 0) { # only NL no TD no WCE
      for (j in 1:i2) { # row 1+1 to 1+nNLb stores TD effects
        MAT[(1+1):(1+nNLb), i1 + i5 + j] <- modX$coef[(i1 +i5 + (j - 1) * nNLb + 1):(i1 +i5 + j * nNLb )]
        
      }
    }
    
    if (i3 != 0) { # only TD no NL no WCE
      for (j in 1:i3) { # row nWCEb+nNLb+1+1 to nWCEb+2*nNLb+1+1 stores TD effects
        MAT[(nNLb+2):(2*nNLb+ 2), i1 + i5 + i2+j] <- modX$coef[(i1 + i5 + i2 * nNLb +(j - 1)*(nNLb+1)+ 1):(i1+i5+ i2*nNLb+ j*(nNLb+1))]
      }
    }
    
    if(i6!=0){ # NLWCE
      for (j in 1:i6){
       MAT[(2):(1+nNLb), i1+i2+i3+i5+j]<- modX$coef[(i1+i5+(i2+i3)* nNLb+i3 +(j-1)* nNLb + 1):(i1+i5+(i2+i3)*nNLb+i3 +j*nNLb )]
       }                     
    }
    
    if(i7!=0){ # WCETD
      for(j in 1:i7){
      MAT[(nNLb+2):(2*nNLb+2),i1 +i2+i3+i5+i6+j]<-modX$coef[(i1+i5+(i2+i6+i3)*nNLb+i3+(j - 1)*(nNLb+ 1)
                                  + 1):(i1+i5+(i2+i6+i3)*nNLb+i3+j*(nNLb + 1))]
      }
    }
    if(i4!=0){ #TDNL
      for(j in 1:i4){
        MAT[(2):(1+nNLb), i1+i2+i3+i5+i6+i7+j]<-modNL$coef[(i1+i5+(i2+i3+i6+i7)*nNLb+i3+i7
                            + (j-1)*nNLb+1):(i1+i5+(i2+i3+i6+i7)*nNLb+i3+i7 +j*nNLb)]
        MAT[(nNLb+2):(2*nNLb+2), i1+i2+i3+i5+i6+i7+j]<-modTD$coef[(i1+i5+(i2+i3+i6+i7)*nNLb+(i3+i7)
                      + (j - 1) * (nNLb + 1) + 1):(i1 +i5 + (i2+i3+i6+i7)*nNLb+(i3+i7)+j*(nNLb + 1))]
        
      }
    }
    
    if(i8!=0){ #NL+TD+WCE
      for (j in 1:i8) {
      MAT[(2):(1+nNLb), i1+i2+i3+i4+i5+i6+i7+j] <- modNL$coef[(i1 + i5 + (i2+i3+i6+i7+i4)*nNLb +i3+i7+(j-1)*nNLb 
                                  + 1):(i1 +i5 + (i2+i3+i6+i7+i4+j)*nNLb +i3+i7 )]
      MAT[(nNLb+2):(2*(nNLb)+2), i1 + i2 + i3+i4+i5+i6+i7+j]<-modTD$coef[(i1 +i5 + (i2+i3+i6+i7+i4)*nNLb +i3+i7+i4+(j-1)*(nNLb+1) 
                               + 1):(i1 +i5 + (i2+i3+i6+i7+i4)*nNLb+i3+i7+i4 +j*(nNLb+1))]
        
      }
    }
    
    ### modified stops here, run the model first ##### 
    
    MATse <- matrix(ncol = V, nrow = 2*(1 + nNLb)) # store the estimated SE
    if (i1 != 0) {
      for (j in 1:i1) {
        MATse[1, j] <- sqrt(diag(modX$var)[j])
      }
    }
    if (i5 != 0) {
      for (j in 1:i5) {
        MATse[1, i1+j]<-sqrt(diag(modX$var)[i1+j])
      }
    }
    if (i2 != 0) {
      for (j in 1:i2) {
        MATse[(2):(nNLb+ 1), i1 + i5 + j] <- sqrt(diag(modX$var)[(i1 + i5 + (j - 1) * nNLb + 1):(i1 +i5 + j * nNLb )])
      }
    }
    
    if (i3 != 0) { # only TD no NL no WCE
      for (j in 1:i3) { # row m+p+2 to 2*(m+p+1) stores TD effects
        MATse[(nNLb+2):(2*nNLb+ 2), i1 + i2 + i5+j] <-sqrt(diag(modX$var)[(i1 + 
                                          i5 + i2 * nNLb +(j - 1) * (nNLb+1)+ 1):(i1 +i5 +i2 * nNLb+ j * (nNLb+ 1) )])
      }
    }
    if(i6!=0){ # NLWCE
      for (j in 1:i6){
        MATse[(2):(1+nNLb), i1+i2+i3+i5+j]<- sqrt(diag(modX$var)[(i1+i5+(i2+i3)* nNLb+i3 +(j-1)* nNLb + 1):(i1+i5+(i2+i3)*nNLb+i3 +j*nNLb )])
      }                     
    }
    if(i7!=0){ # WCETD
      for(j in 1:i7){
        MATse[(nNLb+2):(2*nNLb+2),i1 +i2+i3+i5+i6+j]<-sqrt(diag(modX$var)[(i1+i5+(i2+i6+i3)*nNLb+i3+(j - 1)*(nNLb+ 1)
                                                               + 1):(i1+i5+(i2+i6+i3)*nNLb+i3+j*(nNLb + 1))])
      }
    }
    
    
    if(i4!=0){ #TDNL
      for(j in 1:i4){
        MATse[(2):(1+nNLb), i1+i2+i3+i5+i6+i7+j]<-sqrt(diag(modNL$var)[(i1+i5+(i2+i3+i6+i7)*nNLb+i3+i7
                                                            + (j-1)*nNLb+1):(i1+i5+(i2+i3+i6+i7)*nNLb+i3+i7 +j*nNLb)])
        MATse[(nNLb+2):(2*nNLb+2), i1+i2+i3+i5+i6+i7+j]<-sqrt(diag(modTD$var)[(i1+i5+(i2+i3+i6+i7)*nNLb+(i3+i7)
                                                           + (j - 1) * (nNLb + 1) + 1):(i1 +i5 + (i2+i3+i6+i7)*nNLb+(i3+i7)+j*(nNLb + 1))])
        
      }
    }
    
    if(i8!=0){ #NL+TD+WCE
      for (j in 1:i8) {
        MATse[(2):(1+nNLb), i1+i2+i3+i4+i5+i6+i7+j] <-sqrt(diag(modNL$var)[(i1 + i5 + (i2+i3+i6+i7+i4)*nNLb +i3+i7+(j-1)*nNLb 
                                                                 + 1):(i1 +i5 + (i2+i3+i6+i7+i4+j)*nNLb +i3+i7 )])
        MATse[(nNLb+2):(2*(nNLb)+2), i1 + i2 + i3+i4+i5+i6+i7+j]<-sqrt(diag(modTD$var)[(i1 +i5 + (i2+i3+i6+i7+i4)*nNLb +i3+i7+i4+(j-1)*(nNLb+1) 
                                                                            + 1):(i1 +i5 + (i2+i3+i6+i7+i4)*nNLb+i3+i7+i4 +j*(nNLb+1))])
        
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
    
    coefficients_splines_NL<-as.matrix(MAT[(2):(nNLb+1),])
    coefficients_splines_NL<-rbind(rep(0,V),coefficients_splines_NL)
    # coefficients_splines_NL[1,(NL==0)]<-NA
    colnames(coefficients_splines_NL)<-variables[match(var_order, variables)]
    
    coefficients_splines_TD<-as.matrix(MAT[(nNLb+2):(2*nNLb+2),])
    colnames(coefficients_splines_TD)<-variables[match(var_order, variables)]
    
    knots_covariates<-knotsNEW[c(1:V)[match(var_order, variables)]]
    
    
    if (V>1) {knots_covariates[ match(variables,var_order)[NL==0]]<-NA} else { if (NL==0) {knots_covariates<-NA} }
    knots_time<-knotsNEW[[V+1]]
    
    

    n.events<-sum(data[,TypeNEW[3]]==1)
    PL.c<-vrais[length(vrais)]
    npara<-i1 +i5+(i2+i4+i6+i8)*nNLb +(i3+i4+i7+i8)*(nNLb + 1)
    
    BIC.c <- .my.bic(PL.c, n.events , npara , aic)
    
    rm(data)
    
    list(Partial_Log_Likelihood =PL.c,
         aic=aic,info.criterion=BIC.c,Number_of_parameters =npara, Number_events=n.events, 
         variables=variables[match(var_order, variables)],
         knots_covariates = knots_covariates, 
         knots_time = knots_time, 
         WCEmat=weights,
         WCE=WCE[match(var_order, variables)],
         TD=TD[match(var_order, variables)],
         NL=NL[match(var_order, variables)],
         pTDNL=pTDNL,
         coefficients = coefficients, Standard_Error=se_coef,
         coefficients_splines_NL = coefficients_splines_NL,
         coefficients_splines_TD = coefficients_splines_TD)
    
}
