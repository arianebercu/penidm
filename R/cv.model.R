### Code:
##' @title idm cv model update of beta in penalised regression
##' @param beta  parameters on explanatory variables
##' @param nva01 number of variables for transition 0 -->1 
##' @param nva02 number of variables for transition 0 -->2
##' @param nva12 number of variables for transition 1 -->2
##' @param fix indicators of fixed and unfixed parameters
##' @param penalty base::which penalty to consider
##' @param penalty.factor base::which variable should be penalised
##' @param v variance covariance matrix 
##' @param fu -loglikelihood
##' @param lambda lambda penalised parameter
##' @param alpha alpha penalised parameter
#' @useDynLib SmoothHazardoptim9
#' @author R: Ariane Bercu <ariane.bercu@@u-bordeaux.fr>  
cv.model<-function(beta,
                   nva01,
                   nva02,
                   nva12,
                   fix,
                   penalty.factor,
                   penalty,
                   v,
                   fu,
                   lambda,
                   alpha){
  
# add to do base::which for CRAN check 
  BETA<-beta[fix==0]
  NEW.BETA.all<-beta
  penalty.factor<-penalty.factor[fix==0]
  
  num<-sapply(c(1:dim(v)[1]),FUN=function(x){
    fu[x]-sum(v[x,-x]*BETA[-x])+sum(BETA*v[x,])
  }) 

  
  sign<-ifelse(num<0,-1,
               ifelse(num>0,1,0))
  denum<-diag(v)
  num<-abs(num)
  
  num01<-NULL
  denum01<-NULL
  sign01<-NULL
  
  num02<-NULL
  denum02<-NULL
  sign02<-NULL
  
  num12<-NULL
  denum12<-NULL
  sign12<-NULL
  
  if(nva01>0){
    num01<-num[1:nva01]
    denum01<-denum[1:nva01]
    sign01<-sign[1:nva01]
  }
  
  if(nva02>0){
    num02<-num[(nva01+1):(nva01+nva02)]
    denum02<-denum[(nva01+1):(nva01+nva02)]
    sign02<-sign[(nva01+1):(nva01+nva02)]
    
  }
  
  if(nva12>0){
    num12<-num[(nva01+nva02+1):length(num)]
    denum12<-denum[(nva01+nva02+1):length(denum)]
    sign12<-sign[(nva01+nva02+1):length(num)]
  }
  
  NEWBETA<-rep(NA,length(num))
  idbeta<-NULL
  idbeta<-base::which(penalty.factor==0)
  # if no penalty on parameter, beta_k=A_k/-x_kk
  NEWBETA[idbeta]<-sign[idbeta]*num[idbeta]/denum[idbeta]

  # if penalty update beta all at once 
  if(penalty%in%c("lasso","ridge","elasticnet")){
    # 0 -> 1
    if(nva01>0){
    idbeta<-base::which(num01>(lambda[,1]*alpha))
    NEWBETA[idbeta]<-sign01[idbeta]*(num01[idbeta]-lambda[,1]*alpha)/(denum01[idbeta]+2*lambda[,1]*(1-alpha))
    idbeta<-base::which(num01<=(lambda[,1]*alpha))
    NEWBETA[idbeta]<-0}

    
    # 0 ->2
    if(nva02>0){
    idbeta<-base::which(num02>(lambda[,2]*alpha))
    NEWBETA[idbeta+nva01]<-sign02[idbeta]*(num02[idbeta]-lambda[,2]*alpha)/(denum02[idbeta]+2*lambda[,2]*(1-alpha))
    idbeta<-base::which(num02<=(lambda[,2]*alpha))
    NEWBETA[idbeta+nva01]<-0}
    

    # 1 ->2
    if(nva12>0){
    idbeta<-base::which(num12>(lambda[,3]*alpha))
    NEWBETA[idbeta+nva01+nva02]<-sign12[idbeta]*(num12[idbeta]-lambda[,3]*alpha)/(denum12[idbeta]+2*lambda[,3]*(1-alpha))
    idbeta<-base::which(num12<=(lambda[,3]*alpha))
    NEWBETA[idbeta+nva01+nva02]<-0}
  }
    
  
  
  
  if(penalty=="mcp"){
    
    
    # 0 -> 1, 
    idbeta<-base::which((num01>=(alpha*lambda[,1]*denum01)) & (num01<lambda[,1]) & (denum01<(1/alpha)))
    NEWBETA[idbeta]<-sign01[idbeta]*alpha*lambda[,1]
    # no definition put 0 ? 
    idbeta<-base::which(((num01<(alpha*lambda[,1]*denum01)) | (num01>=lambda[,1])) & (denum01<(1/alpha)))
    NEWBETA[idbeta]<-0
    
    idbeta<-base::which((num01<=(alpha*lambda[,1]*denum01)) & num01>lambda[,1] & (denum01>=(1/alpha)))
    NEWBETA[idbeta]<-sign01[idbeta]*(num01[idbeta]-lambda[,1])/(denum01[idbeta]-(1/alpha))
    idbeta<-base::which(((num01>(alpha*lambda[,1]*denum01)) | num01<=lambda[,1]) & (denum01>=(1/alpha)))
    NEWBETA[idbeta]<-0
    
    # 0 -> 2
    idbeta<-base::which((num02>=(alpha*lambda[,2]*denum02)) & (num02<lambda[,2]) & (denum02<(1/alpha)))
    NEWBETA[idbeta+nva01]<-sign02[idbeta]*alpha*lambda[,2]
    idbeta<-base::which(((num02<(alpha*lambda[,2]*denum02)) | (num02>=lambda[,2])) & (denum02<(1/alpha)))
    NEWBETA[idbeta+nva01]<-0
    
    idbeta<-base::which((num02<=(alpha*lambda[,2]*denum02)) & num02>lambda[,2] & (denum02>=(1/alpha)))
    NEWBETA[idbeta+nva01]<-sign02[idbeta]*(num02[idbeta]-lambda[,2])/(denum02[idbeta]-(1/alpha))
    idbeta<-base::which(((num02>(alpha*lambda[,2]*denum02)) | num02<=lambda[,2]) & (denum02>=(1/alpha)))
    NEWBETA[idbeta+nva01]<-0
    
    
    # 1 -> 2
    idbeta<-base::which((num12>=(alpha*lambda[,3]*denum12)) & (num12<lambda[,3]) & (denum12<(1/alpha)))
    NEWBETA[idbeta+nva01+nva02]<-sign12[idbeta]*alpha*lambda[,3]
    idbeta<-base::which(((num12<(alpha*lambda[,3]*denum12)) | (num12>=lambda[,3])) & (denum12<(1/alpha)))
    NEWBETA[idbeta+nva01+nva02]<-0
    
    idbeta<-base::which((num12<=(alpha*lambda[,3]*denum12)) & num12>lambda[,3] & (denum12>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02]<-sign12[idbeta]*(num12[idbeta]-lambda[,3])/(denum12[idbeta]-(1/alpha))
    idbeta<-base::which(((num12>(alpha*lambda[,3]*denum12)) | num12<=lambda[,3]) & (denum12>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02]<-0
    
    }
  
  
  if(penalty=="scad"){
    
    
    
    # 0 -> 1
    idbeta<-base::which((num01<=(lambda[,1]*(1+denum01))) & (num01 >lambda[,1]) & (denum01>=(1/(alpha-1))))
    NEWBETA[idbeta]<-sign01[idbeta]*(num01[idbeta]-lambda[,1])/denum01[idbeta]
    idbeta<-base::which(((num01>(lambda[,1]*(1+denum01))) | (num01 <= lambda[,1]*denum01)) & (denum01>=(1/(alpha-1))))
    NEWBETA[idbeta]<-0
    
    idbeta<-base::which((num01<=(alpha*lambda[,1]*denum01)) & (num01 > lambda[,1]) & (denum01<(1/(alpha-1))) & (denum01>=(1/alpha)))
    NEWBETA[idbeta]<-sign01[idbeta]*(num01[idbeta]-lambda[,1])/denum01[idbeta]
    idbeta<-base::which(((num01>(alpha*lambda[,1]*denum01)) | (num01 <= lambda[,1])) & (denum01<(1/(alpha-1))) & (denum01>=(1/alpha)))
    NEWBETA[idbeta]<-0
    idbeta<-base::which( (denum01<(1/(alpha-1))) & (denum01<(1/alpha)))
    NEWBETA[idbeta]<-0
    
    # 0 ->2
    idbeta<-base::which((num02<=(lambda[,2]*(1+denum02))) & (num02 >lambda[,2]) & (denum02>=(1/(alpha-1))))
    NEWBETA[idbeta+nva01]<-sign02[idbeta]*(num02[idbeta]-lambda[,2])/denum02[idbeta]
    idbeta<-base::which(((num02>(lambda[,2]*(1+denum02))) | (num02 <= lambda[,2]*denum02)) & (denum02>=(1/(alpha-1))))
    NEWBETA[idbeta+nva01]<-0
    
    idbeta<-base::which((num02<=(alpha*lambda[,2]*denum02)) & (num02 > lambda[,2]) & (denum02<(1/(alpha-1))) & (denum02>=(1/alpha)))
    NEWBETA[idbeta+nva01]<-sign02[idbeta]*(num02[idbeta]-lambda[,2])/denum02[idbeta]
    idbeta<-base::which(((num02>(alpha*lambda[,2]*denum02)) | (num02 <= lambda[,2])) & (denum02<(1/(alpha-1))) & (denum02>=(1/alpha)))
    NEWBETA[idbeta+nva01]<-0
    idbeta<-base::which( (denum02<(1/(alpha-1))) & (denum02<(1/alpha)))
    NEWBETA[idbeta+nva01]<-0
    
    # 1 ->2
    
    idbeta<-base::which((num12<=(lambda[,3]*(1+denum12))) & (num12 >lambda[,3]) & (denum12>=(1/(alpha-1))))
    NEWBETA[idbeta+nva01+nva02]<-sign12[idbeta]*(num12[idbeta]-lambda[,3])/denum12[idbeta]
    idbeta<-base::which(((num12>(lambda[,3]*(1+denum12))) | (num12 <= lambda[,3]*denum12)) & (denum12>=(1/(alpha-1))))
    NEWBETA[idbeta+nva01+nva02]<-0
    
    idbeta<-base::which((num12<=(alpha*lambda[,3]*denum12)) & (num12 > lambda[,3]) & (denum12<(1/(alpha-1))) & (denum12>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02]<-sign12[idbeta]*(num12[idbeta]-lambda[,3])/denum12[idbeta]
    idbeta<-base::which(((num12>(alpha*lambda[,3]*denum12)) | (num12 <= lambda[,3])) & (denum12<(1/(alpha-1))) & (denum12>=(1/alpha)))
    NEWBETA[idbeta+nva01+nva02]<-0
    idbeta<-base::which( (denum12<(1/(alpha-1))) & (denum12<(1/alpha)))
    NEWBETA[idbeta+nva01+nva02]<-0
    }

  
  NEW.BETA.all[fix==0]<-NEWBETA
  
  return(list(b=NEW.BETA.all))
}

  