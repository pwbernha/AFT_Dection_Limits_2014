#####################################################################################
#NOTE: PLEASE DO NOT USE CODE FOR ANY PUBLICATION PURPOSES WITHOUT FIRST CONTACTING #
#      PAUL BERNHARDT AT PAUL.BERNHARDT@VILLANOVA.EDU                               #
#####################################################################################

#####################Explanation of Main Functions in Program Below##############################
# SNPMI, PARMI: main program functions: obtain estimates for the multiple imputatio 		#
# 		    algorithm and iterated multiple imputation algorithm (depending on  	#
#		    whether "M" is defined in the "Simulation Parameters" section below as  	#
# 		    a number or vector); SNPMI is used when SNP distribution is assumed for	#
#		    error in AFT model, PARMI is used when other error distributions are 	#
#		    used; see below in "Main Computational Functions" section for more 	     	#
#		    details								        #                                                        					
# SNPchoice:    used to obtain the SNP kernel (exponential or normal), the degree of the 	#
#		    SNP polynomial, K (=0,1,2), and the initial value and variance estimates	#
# 		    for the parameters in the AFT-SNP; this function must be called before	#
#		    the SNPMI since we need the initial values to obtain imputations; note      #
#		    that no similar function exists for PARMI since the finding of initial      #
#		    values is easily conducted (and is thus included in the PARMI function      #
#		    below using the survreg function in the survival package)		        #                  	
# MHMI, MHMI2:  called within SNPMI and PARMI, respectively, and are the functions which      	#
#		    actually draw the multiple imputations				        #                                	
# cendata:      takes in the generated survival data, covariate data, and censoring data      	#
#		    and orders in chunks corresponding to which surivival/covariate 	        #
#		    variables are censored; given a set of data, this function must be applied  #
#		    first in order to properly format the data for input in other functions     #
# Datagen:	generates the survival, covariate, and censoring variables for all 		#
#		    datasets in the simulation							#                                                
# CensBN:       likelihood function for censored f(x|z) distribution called within an		#
#                   optim function (built-in R optimization function) and used to get maximum 	#
#		    likelihood estimates for the parameters in f(x|z)			        #                        	
#################################################################################################


###################Necessary Packages########################
library(mvtnorm) 	#generate/evaluate points of a multivariate normal
library(tmvtnorm) 	#generate from a truncated normal
library(msm) 		#generated from a truncated normal
library(corpcor)	#checking positive definiteness of matrices
library(Matrix) 	#create diagonal matrix of two or more matrices
library(evd) 		#generate/evaluate points from an extreme value distribution
library(survival) 	#run AFT analysis with popular error distributions via survreg function

##################Simulation Parameters######################
N <- 1000		#number of simulated datasets
n <- 500		#number of obsevations per dataset
M <- c(15,15,30,50)	#number of imputations per iteration
divs <- 31		#number of grid point evaluations for each dimension of K (higher number should increase accuracy of ML)
Iters <- length(M)	#number of iterations per dataset

#####################Model Parameters#########################
MU1 <- c(1.4,2.4)		 #mean vector for covariates X1 and X2 when Z=0 (X1,X2 representing cytokine biomarkers and Z representing age)
MU2 <- c(1/115,-1/110)		 #effect of Z (representing age) on X1 and X2 respectively as Z increases by 1
VARI <- matrix(c(2.5,1,1,5),2,2) #covariance matrix for X1,X2
d <- c(log(4),log(5))		 #detection limits for X1,X2
Beta <- c(5,-0.2,0.2,-0.1)	 #Parameter vector for intercept and effect of X1,X2 and Z on the log-survival


#########################Simulation###########################
set.seed(8089305)
Data <- Datagen(N,n,MU1,MU2,VARI,Beta,d,distribution="Mixture")  #Generates data for all N datasets with chosen distribution for T
Xgen <- Data[[1]]
Tgen <- Data[[2]]
Cgen <- Data[[3]]

#Initialization of vectors and matrices for storing simulation data
InitBN <- InitSNP <- FinalSNP1 <- FinalSNP2 <- FinalSNP3 <- FinalSNP4 <-  FinalBN <- FinalBN2 <- matrix(0,7,N)
InitPar <- FinalPar1 <- FinalPar2 <- FinalPar3 <- FinalPar4 <- Initdist <- Finaldist <- matrix(0,5,N)
SEsnp1_a <- SEsnp2_a <- SEsnp3_a <- SEsnp4_a <- SEsnp1_c <- SEsnp2_c <- SEsnp3_c <- SEsnp4_c <- SEPar1 <- SEPar2 <- SEPar3 <- SEPar4 <-  matrix(0,5,N)
K <- EXPd <- dists <- rep(0,N)
Aexp <- matrix(c(1,1,2,1,2,6,2,6,24),3,3)	#A matrix in SNP fitting for exponential (see Zhang and Davidian 2001)
Anorm <- matrix(c(1,0,1,0,1,0,1,0,3),3,3) #A matrix in SNP fitting for normal 

#Loops through all N datasets and calculates Beta and SE estimates for each method
for(j in 1:N){ 
  
  #Gives an idea on progress
  print(j)
  
  #Dividing data into observed and censored data
  CD <- cendata(Xgen[,((4*j-3):(4*j))],Tgen[,j],Cgen[,j])
  T <- CD[[1]]
  C <- CD[[2]]
  X <- CD[[3]]
  cc <- CD[[4]][1]	#no censoring
  mt <- CD[[4]][2]	#censored t's
  m1 <- CD[[4]][3]	#censored X1's
  mt1 <- CD[[4]][4]	#censored t's, X1's
  m2 <- CD[[4]][5] 	#censored X2's
  mt2 <- CD[[4]][6]	#censored t's, X2's
  m12 <- CD[[4]][7]	#censored X1's, X2's
  mt12 <- CD[[4]][8]	#all censored
  
  #Maximum Likelihood Estimates for BVN only (starting values being correct for simulation purposes only
  MAX <- optim(c(MU1[1],MU2[1],MU1[2],MU2[2],VARI[1,1],VARI[1,2],VARI[2,2]),CensBN(X,(cc+mt),(m1+mt1),(m2+mt2),(m12+mt12)),method="BFGS",hessian=TRUE)
  InitBN[,j] <- as.vector(MAX$par)
  InitVarBN <- solve(MAX$hessian)
  
  ########### Obtaining Estimates using SNP distribution ###########
  SNPch <- SNPchoice(T,C,X,divs,cc,mt) #Find complete case parameter estimates for SNP for multiple imputation
  
  InitSNP[1:length(SNPch[[1]]),j] <- SNPch[[1]]
  InitVar <- as.matrix(bdiag(SNPch[[2]], InitVarBN))
  K[j] <- SNPch[[3]]
  EXPd[j] <- SNPch[[4]]
  
  #Obtains final IMI (or MI if Iters=1)SNP AFT estimates using complete case estimates
  SNPmi <- SNPMI(T,C,X,cc,mt,m1,mt1,m2,mt2,m12,mt12,InitSNP[1:length(SNPch[[1]]),j],InitBN[,j],InitVar,K[j],EXPd[j],M)
  FinalSNP1[1:length(SNPch[[1]]),j] <- SNPmi[[1]][1,]
  FinalSNP2[1:length(SNPch[[1]]),j] <- SNPmi[[1]][2,]
  FinalSNP3[1:length(SNPch[[1]]),j] <- SNPmi[[1]][3,]
  FinalSNP4[1:length(SNPch[[1]]),j] <- SNPmi[[1]][4,]
  FinalBN[,j] <- rowMeans(SNPmi[[4]])
  SEsnp1_a[,j] <- SNPmi[[2]][1,]
  SEsnp2_a[,j] <- SNPmi[[2]][2,]
  SEsnp3_a[,j] <- SNPmi[[2]][3,]
  SEsnp4_a[,j] <- SNPmi[[2]][4,]
  #SEsnp1_c[,j] <- SNPmi[[3]][1,]
  #SEsnp2_c[,j] <- SNPmi[[3]][2,]
  #SEsnp3_c[,j] <- SNPmi[[3]][3,]
  #SEsnp4_c[,j] <- SNPmi[[3]][4,]
  
  #Obtains final IMI (or MI if Iters=1) AFT estimates using parametric assumptions
  Parmi <- PARMI(T,C,X,cc,mt,m1,mt1,m2,mt2,m12,mt12,InitBN[,j],InitVarBN,M,Weib=FALSE)
  FinalPar1[,j] <- Parmi[[1]][1,]
  FinalPar2[,j] <- Parmi[[1]][2,]
  FinalPar3[,j] <- Parmi[[1]][3,]
  FinalPar4[,j] <- Parmi[[1]][4,]
  FinalBN2[,j] <- rowMeans(Parmi[[3]])
  SEPar1[,j] <- Parmi[[2]][1,]
  SEPar2[,j] <- Parmi[[2]][2,]
  SEPar3[,j] <- Parmi[[2]][3,]
  SEPar4[,j] <- Parmi[[2]][4,]
}



#####################Data Generation Function########################
#Generates Weibull, exponential, lognormal, loglogistic, and mixture survival variables
Datagen <- function(N,n,MU1,MU2,VARI,Beta,d,distribution,Scale=1,Scale2=1){
  Tgen <- matrix(0,n,N)
  Cgen <- matrix(0,n,N)
  Xgen <- matrix(0,n,4*N)
  
  for(j in 1:N){
    Xgen[,(4*j)] <- rbeta(n,3,2)*84+18  #Generates "age" covariate
    Mean <- (matrix(MU1,n,2,byrow=TRUE)+matrix(c(MU2[1]*Xgen[,(4*j)],MU2[2]*Xgen[,(4*j)]),n,2)) #Mean for "cytokine" covariates
    Vals <- rmvnorm(n,c(0,0),VARI)	#Generates "cytokine" covariates
    Xgen[,((4*j-3):(4*j-1))] <- cbind(1,(Mean+Vals))
    
    #T ~ Weibull
    if(distribution=="Weibull"){
      Tgen[,j] <- Xgen[,((4*j-3):(4*j))]%*%Beta - rgumbel(n,scale=Scale) 
      Cgen[,j] <- runif(n,-10,7)  #can be played with/altered to get desired censoring
    }
    
    #T ~ Exponential
    if(distribution=="Exponential"){
      Tgen[,j] <- Xgen[,((4*j-3):(4*j))]%*%Beta - rgumbel(n,scale=1)
      Cgen[,j] <- runif(n,-9,6)  #can be played with/altered to get desired censoring
    }
    
    #T ~ Lognormal
    if(distribution=="Lognormal"){
      Tgen[,j] <- rnorm(n,Xgen[,((4*j-3):(4*j))]%*%Beta,Scale)
      Cgen[,j] <- runif(n,-4,8) #can be played with/altered to get desired censoring
    }
    
    #T ~ Loglogistic
    if(distribution=="Loglogistic"){
      Tgen[,j] <- rlogis(n,Xgen[,((4*j-3):(4*j))]%*%Beta,Scale)
      Cgen[,j] <- runif(n,-10,12) #can be played with/altered to get desired censoring
    }
    
    #T ~ Mixture of log-normal and Weibull
    if(distribution=="Mixture"){
      p <- rbinom(n,1,0.6)
      Tgen[,j] <- Xgen[,((4*j-3):(4*j))]%*%Beta+ p*rnorm(n,mean=0,Scale)+(1-p)*(5+rgumbel(n,scale=Scale2)) 
      Cgen[,j] <- runif(n,-5,22) #can be played with/altered to get desired censoring
    }
  }
  return(list(Xgen,Tgen,Cgen))
}


##################Main Computational Functions#####################

#Obtains IMI (or MI if M is of length 1) estimate for parameters using SNP model; inputs survival (T) and covariate (X) 
#data, censoring values (C), initial parameter estimates for f(x|z) - InitBN, initial parameter estimates
#for f(t|x,z) - InitSNP, K, and EXP (kernel, == TRUE -> exponential, == FALSE -> normal) obtained via the 
#SNPchoice function (see below);
#outputs parameter and standard error estimates for all model parameters (when M is a vector, indicating the
#iterated multiple imputation method, a matrix of estimates is outputed corresponding to the update at
#each iteration)
#Note: cc, mt, m1, mt1, m2, mt2, m12, and mt12 are found via the censdata function (see below) and correspond 
#to the amount of observed data points with censoring on no variables (cc), the survival variable (t) and the
#first and second covariates subject to censoring (1 and 2)
SNPMI <- function(T,C,X,cc,mt,m1,mt1,m2,mt2,m12,mt12,InitSNP,InitBN,InitVar,K,EXP,M){
  Iters <- length(M)
  EstUpdates <- matrix(0,Iters,length(InitSNP))
  SEUpdates_a  <- matrix(0,Iters,5)
  SEUpdates_c  <- matrix(0,Iters,5)
  
  for(i in 1:Iters){
    if(i==1){
      Initsnp <- InitSNP
      Initbn <- InitBN
    }
    
    if(i>1){
      Initsnp <- rowMeans(SNPmis)
      Initbn <- rowMeans(BNmis)
      #if(EXP==1) InitVar[c(1:5,(1+length(InitSNP)):(7+length(InitSNP))),c(1:5,(1+length(InitSNP)):(7+length(InitSNP)))] <- VarExp(Timped,Cimped,Ximped,(cc+m1+m2+m12),(mt+mt1+mt2+mt12),rbind(SNPmis[1:5,],BNmis),matrix(SNPmis[6:(6+K-1),],ncol=M[i]),K,M[i-1],InitVar[c(1:5,(1+length(InitSNP)):(7+length(InitSNP))),c(1:5,(1+length(InitSNP)):(7+length(InitSNP)))],FinalVar[c(1:5,(1+length(InitSNP)):(7+length(InitSNP))),c(1:5,(1+length(InitSNP)):(7+length(InitSNP)))]) else InitVar[c(1:5,(1+length(InitSNP)):(7+length(InitSNP))),c(1:5,(1+length(InitSNP)):(7+length(InitSNP)))] <- VarNorm(Timped,Cimped,Ximped,(cc+m1+m2+m12),(mt+mt1+mt2+mt12),rbind(SNPmis[1:5,],BNmis),matrix(SNPmis[6:(6+K-1),],ncol=M[i]),K,M[i-1],InitVar[c(1:5,(1+length(InitSNP)):(7+length(InitSNP))),c(1:5,(1+length(InitSNP)):(7+length(InitSNP)))],FinalVar[c(1:5,(1+length(InitSNP)):(7+length(InitSNP))),c(1:5,(1+length(InitSNP)):(7+length(InitSNP)))])
      InitVar <- ApproxVar(SNPmis,BNmis,InitVar,FinalVar,M[i-1])
      
    }
    BNmis <- matrix(0,7,M[i])
    FinalVar <-matrix(0,(length(InitSNP)+7),(length(InitSNP)+7))
    Int <- rep(0,M[i])
    SNPmis <- matrix(0,length(InitSNP), M[i])
    Full <- MHMI(T,C,X,K,Initsnp,Initbn,cc,mt,m1,mt1,m2,mt2,m12,mt12,EXP,M[i])
    Timped <- Full[[3]]
    Cimped <- Full[[4]]
    Ximped <- Full[[5]]
    for(d in 1:M[i]){
      Tnew <- Timped[,d]
      Cnew <- Cimped[,d]
      Xnew <- Ximped[,(d*4-3):(d*4)]
      FINAL <- Optimize(InitVals(Tnew,Cnew,Xnew,(cc+m1+m2+m12),(mt+mt1+mt2+mt12),InitSNP,PhiVals,K,EXP,FINAL=TRUE),Tnew,Cnew,Xnew,K,(cc+m1+m2+m12),(mt+mt1+mt2+mt12),EXP)
      SNPmis[,d] <- FINAL$par
      if(EXP==1) VAR <- nlm(CensSNPexp(Tnew,Cnew,Xnew,K,(cc+m1+m2+m12),(mt+mt1+mt2+mt12)),c(FINAL$par),hessian=TRUE,iterlim=1)$hessian else VAR <- nlm(CensSNPnorm(Tnew,Cnew,Xnew,K,(cc+m1+m2+m12),(mt+mt1+mt2+mt12)),c(FINAL$par),hessian=TRUE,iterlim=1)$hessian
      MAX2 <- optim(InitBN,CensBN(Xnew,n,0,0,0), method="BFGS", hessian=TRUE)
      FinalVar <- FinalVar + solve(as.matrix(bdiag(VAR,MAX2$hessian)))/M[i]
      BNmis[,d] <- as.vector(MAX2$par)
      if(K>0 & EXP==1) Int[d] <- SNPmis[1,d] + exp(SNPmis[5,d])*log(integrate(SNPexpX,lower=0,upper=Inf,a=c(acoef(SNPmis[6:length(InitSNP),d],K,EXP=EXP),rep(0,2-K)))$value) #Update Intercept
      if(K>0 & EXP==0) Int[d] <- SNPmis[1,d] + SNPmis[5,d]*integrate(SNPnormX,lower=-Inf,upper=Inf,a=c(acoef(SNPmis[6:length(InitSNP),d],K,EXP=EXP),rep(0,2-K)))$value
    }
    EstUpdates[i,] <- rowMeans(SNPmis)
    EstUpdates[i,1] <- mean(Int)
    SEUpdates_a[i,] <- sqrt(diag(ApproxVar(SNPmis[1:5,],BNmis,InitVar[c(1:5,(1+length(InitSNP)):(7+length(InitSNP))),c(1:5,(1+length(InitSNP)):(7+length(InitSNP)))],FinalVar[c(1:5,(1+length(InitSNP)):(7+length(InitSNP))),c(1:5,(1+length(InitSNP)):(7+length(InitSNP)))],M[i])[1:5,1:5]))
    #if(EXP==1) SEUpdates_c[i,] <- sqrt(diag(VarExp(Timped,Cimped,Ximped,(cc+m1+m2+m12),(mt+mt1+mt2+mt12),rbind(SNPmis[1:5,],BNmis),matrix(SNPmis[6:(6+K-1),],ncol=M[i]),K,M[i],InitVar[c(1:5,(1+length(InitSNP)):(7+length(InitSNP))),c(1:5,(1+length(InitSNP)):(7+length(InitSNP)))],FinalVar[c(1:5,(1+length(InitSNP)):(7+length(InitSNP))),c(1:5,(1+length(InitSNP)):(7+length(InitSNP)))])))[1:5] else SEUpdates[i,] <- sqrt(diag(VarNorm(Timped,Cimped,Ximped,(cc+m1+m2+m12),(mt+mt1+mt2+mt12),rbind(SNPmis[1:5,],BNmis),matrix(SNPmis[6:(6+K-1),],ncol=M[i]),K,M[i],InitVar[c(1:5,(1+length(InitSNP)):(7+length(InitSNP))),c(1:5,(1+length(InitSNP)):(7+length(InitSNP)))],FinalVar[c(1:5,(1+length(InitSNP)):(7+length(InitSNP))),c(1:5,(1+length(InitSNP)):(7+length(InitSNP)))])))[1:5]
  }
  
  return(list(EstUpdates,SEUpdates_a,SEUpdates_c,BNmis))
}




#Obtains IMI (or MI) estimates using one of the typical parametric models; inputs survival (T) and 
#covariate (X) data, censoring values (C), and initial parameter estimates for f(x|z) - InitBN, InitVarBN
#outputs parameter and standard error estimates for all model parameters (when M is a vector, indicating the
#iterated multiple imputation method, a matrix of estimates is outputed corresponding to the update at
#each iteration)
PARMI <- function(T,C,X,cc,mt,m1,mt1,m2,mt2,m12,mt12,InitBN,InitVarBN,M,Weib=FALSE){
  Iters <- length(M)
  InitPar <- rep(0,5)
  survreg.control(maxiter=500)
  ###Getting Complete Case Estimates###
  AFT1<-survreg(Surv(exp(c(T[1:cc],C[(cc+1):(cc+mt)])),c(rep(1,cc),rep(0,mt))) ~ X[1:(cc+mt),2] + X[1:(cc+mt),3] + X[1:(cc+mt),4],dist="exponential")
  AFT2<-survreg(Surv(exp(c(T[1:cc],C[(cc+1):(cc+mt)])),c(rep(1,cc),rep(0,mt))) ~ X[1:(cc+mt),2] + X[1:(cc+mt),3] + X[1:(cc+mt),4],dist="weibull")
  AFT3<-survreg(Surv(exp(c(T[1:cc],C[(cc+1):(cc+mt)])),c(rep(1,cc),rep(0,mt))) ~ X[1:(cc+mt),2] + X[1:(cc+mt),3] + X[1:(cc+mt),4],dist="lognorm")
  AFT4<-survreg(Surv(exp(c(T[1:cc],C[(cc+1):(cc+mt)])),c(rep(1,cc),rep(0,mt))) ~ X[1:(cc+mt),2] + X[1:(cc+mt),3] + X[1:(cc+mt),4],dist="loglogistic")
  Models <- list(AFT1,AFT2,AFT3,AFT4)
  Params <- c(-1,0,0,0)
  aic <- Inf
  test <- which(c(is.nan(Models[[1]]$loglik[2])==FALSE,is.nan(Models[[2]]$loglik[2])==FALSE,is.nan(Models[[3]]$loglik[2])==FALSE,is.nan(Models[[4]]$loglik[2])==FALSE)==TRUE)
  for(i in test){
    if(AIC(Models[[i]]$loglik[2],Params[i])<aic){
      aic <- AIC(Models[[i]]$loglik[2],Params[i])
      InitPar[1:4] <- Models[[i]]$coefficients
      InitPar[5] <- Models[[i]]$scale
      dists <- i
      InitVar <- as.matrix(Models[[i]]$var)
    }
  }
  
  if(Weib==TRUE){
    InitPar[1:4] <- Models[[2]]$coefficients
    InitPar[5] <- Models[[2]]$scale
    dists <- 2
    InitVar <- as.matrix(Models[[2]]$var)
  }
  
  InitVar <- as.matrix(bdiag(InitVar,InitVarBN))
  EstUpdates <- matrix(0,Iters,5)
  SEUpdates  <- matrix(0,Iters,5)
  
  for(i in 1:Iters){
    if(i==1){
      Initpar <- InitPar
      Initbn <- InitBN
    }
    if(i>1){
      Initpar <- rowMeans(Parmis)
      Initbn <- rowMeans(BNmis)
      if(dists==1) InitVar <- ApproxVar(Parmis[1:4,],BNmis,InitVar,FinalVar,M[(i-1)]) else InitVar <- ApproxVar(Parmis,BNmis,InitVar,FinalVar,M[(i-1)])
    }
    
    BNmis <- matrix(0,7,M[i])
    if(dists==1) FinalVar <-matrix(0,11,11) else FinalVar <-matrix(0,12,12)
    Parmis <- matrix(0,5, M[i])
    Full <- MHMI2(T,C,X,Initpar,Initbn,cc,mt,m1,mt1,m2,mt2,m12,mt12,dists,M[i])
    Timped <- Full[[3]]
    Cimped <- Full[[4]]
    Ximped <- Full[[5]]
    
    q <- 1
    while(q<=M[i]){
      Tnew <- Timped[,q]
      Cnew <- Cimped[,q]
      Xnew <- Ximped[,(q*4-3):(q*4)]
      if(dists==1)	FINAL <- survreg(Surv(exp(c(Tnew[1:(cc+m1+m2+m12)],Cnew[(cc+m1+m2+m12+1):n])),c(rep(1,(cc+m1+m2+m12)),rep(0,(mt+mt1+mt2+mt12)))) ~ Xnew[,2]+Xnew[,3]+ Xnew[,4],dist="exponential")
      if(dists==2)	FINAL <- survreg(Surv(exp(c(Tnew[1:(cc+m1+m2+m12)],Cnew[(cc+m1+m2+m12+1):n])),c(rep(1,(cc+m1+m2+m12)),rep(0,(mt+mt1+mt2+mt12)))) ~ Xnew[,2]+Xnew[,3]+ Xnew[,4],dist="weibull")
      if(dists==3)	FINAL <- survreg(Surv(exp(c(Tnew[1:(cc+m1+m2+m12)],Cnew[(cc+m1+m2+m12+1):n])),c(rep(1,(cc+m1+m2+m12)),rep(0,(mt+mt1+mt2+mt12))))~ Xnew[,2]+Xnew[,3]+ Xnew[,4],dist="lognorm")
      if(dists==4)	FINAL <- survreg(Surv(exp(c(Tnew[1:(cc+m1+m2+m12)],Cnew[(cc+m1+m2+m12+1):n])),c(rep(1,(cc+m1+m2+m12)),rep(0,(mt+mt1+mt2+mt12)))) ~ Xnew[,2]+Xnew[,3]+ Xnew[,4],dist="loglogistic")
      Parmis[,q] <- c(FINAL$coefficients,FINAL$scale)
      MAX2 <- optim(Initbn,CensBN(Xnew,n,0,0,0), method="BFGS", hessian=TRUE)
      BNmis[,q] <-as.vector(MAX2$par)
      #protecting against various issues which arise when maximizing survival data
      if(sum(is.nan(FINAL$var))==0 & sum(is.na(FINAL$var))==0) FinalVar <- FinalVar + as.matrix(bdiag(FINAL$var/M[i],solve(MAX2$hessian)/M[i]))
      if(sum(is.nan(FINAL$var))>0 | sum(is.na(FINAL$var))>0) {NEW <- MHMI2(T,C,X,Initpar,Initbn,cc,mt,m1,mt1,m2,mt2,m12,mt12,dists,2)
      Timped[,q] <- NEW[[3]][,1]
      Cimped[,q] <- NEW[[4]][,1]
      Ximped[,(q*4-3):(q*4)] <- NEW[[5]][,1:4]
      q <- q-1}
      if(Parmis[1,q]<0.5) {NEW <- MHMI2(T,C,X,Initpar,Initbn,cc,mt,m1,mt1,m2,mt2,m12,mt12,dists,2)
      Timped[,q] <- NEW[[3]][,1]
      Cimped[,q] <- NEW[[4]][,1]
      Ximped[,(q*4-3):(q*4)] <- NEW[[5]][,1:4]
      q <- q-1}
      q <- q + 1
    }
    
    EstUpdates[i,] <- rowMeans(Parmis)
    if(dists==1) SEUpdates[i,] <-  c(sqrt(diag(ApproxVar(Parmis[1:4,],BNmis,InitVar,FinalVar,M[i])[1:4,1:4])),1) else SEUpdates[i,] <- sqrt(diag(ApproxVar(Parmis,BNmis,InitVar,FinalVar,M[i])[1:5,1:5]))
  }
  return(list(EstUpdates,SEUpdates,BNmis))
}




#Obtains multiple imputation when using SNP AFT model - called within SNPMI function
MHMI <- function(T,C,X,k,InitSNP,ests,cc,mt,m1,mt1,m2,mt2,m12,mt12,EXP=TRUE,M) {
  Beta <- InitSNP[1:4]
  Sigma <- InitSNP[5]
  Phi <- InitSNP[6:length(InitSNP)]
  a <- acoef(Phi,k,EXP=EXP)
  s <- sd(T)
  
  X1imp <- matrix(0,M,(m1+mt1))
  if((m1+mt1)==0) X1imp<- NULL
  if((m1+mt1)>0){
    for(i in (cc+mt+1):(cc+mt+m1+mt1)){
      Xi <- rep(100+20*M)
      Xi[1] <- rtnorm(1,ests[1]+ests[2]*X[i,4],sqrt(ests[5])/2,upper=d[1])
      rden <- GenTsnp(InitSNP,k,a,EXP,T[i],C[i],c(1,Xi[1],X[i,3:4]))*dnorm(Xi[1],((ests[1]+ests[2]*X[i,4])+ests[6]/ests[7]*(X[i,3]-(ests[3]+ests[4]*X[i,4]))),sqrt(ests[5]-ests[6]^2/ests[7]))
      for(t in 2:(100+20*M)) {
        Xstar <- rnorm(1,mean=Xi[(t-1)],sqrt(ests[5])/2)
        rnum <- GenTsnp(InitSNP,k,a,EXP,T[i],C[i],c(1,Xstar,X[i,3:4]))*dnorm(Xstar,((ests[1]+ests[2]*X[i,4])+ests[6]/ests[7]*(X[i,3]-(ests[3]+ests[4]*X[i,4]))),sqrt(ests[5]-ests[6]^2/ests[7]))
        r <- rnum/max(rden,1e-300)
        u <- runif(1)
        if(r>u & Xstar<d[1]) {rden <- rnum
        Xi[t] <- Xstar}   else Xi[t] <- Xi[(t-1)]
        
      }
      X1imp[,(i-(cc+mt))] <- sample(Xi[seq(110,(100+20*M),20)])
    }}
  
  X2imp <- matrix(0,M,(m2+mt2))
  if((m2+mt2)==0) X2imp<- NULL
  if((m2+mt2)>0){
    for(i in (cc+mt+m1+mt1+1):(cc+mt+m1+mt1+m2+mt2)){
      Xi <- rep(100+20*M)
      Xi[1] <- rtnorm(1,ests[3]+ests[4]*X[i,4],sqrt(ests[7])/2,upper=d[2])
      rden <- GenTsnp(InitSNP,k,a,EXP,T[i],C[i],c(1,X[i,2],Xi[1],X[i,4]))*dnorm(Xi[1],((ests[3]+ests[4]*X[i,4])+ests[6]/ests[5]*(X[i,2]-(ests[1]+ests[2]*X[i,4]))),sqrt(ests[7]-ests[6]^2/ests[5]))
      for(t in 2:(100+20*M)) {
        Xstar <- rnorm(1,mean=Xi[(t-1)],sqrt(ests[7])/2)
        rnum <- GenTsnp(InitSNP,k,a,EXP,T[i],C[i],c(1,X[i,2],Xstar,X[i,4]))*dnorm(Xstar,((ests[3]+ests[4]*X[i,4])+ests[6]/ests[5]*(X[i,2]-(ests[1]+ests[2]*X[i,4]))),sqrt(ests[7]-ests[6]^2/ests[5]))
        r <- rnum/max(rden,1e-300)
        u <- runif(1)
        if(r>u & Xstar<d[2]) {rden <- rnum
        Xi[t] <- Xstar} else Xi[t] <- Xi[(t-1)]
      }
      X2imp[,(i-(cc+mt+m1+mt1))] <- sample(Xi[seq(120,(100+20*M),20)])
    }}
  
  Genvals <- rmvnorm((50*M+200)*(m12+mt12),c(0,0),make.positive.definite(matrix(c(ests[5]/4,ests[6]/4,ests[6]/4,ests[7]/4),2,2),tol=0.001))
  l <- 1
  X12imp <- matrix(0,M,(m12+mt12)*2)
  if((m12+mt12)==0) X12imp<- NULL
  if((m12+mt12)>0){
    for(i in (cc+mt+m1+mt1+m2+mt2+1):(cc+mt+m1+mt1+m2+mt2+m12+mt12)){
      Xi <- matrix(0,(50*M+200),2)
      Xi[1,] <- rtmvnorm(1,mean=c(ests[1]+ests[2]*(X[i,4]),ests[3]+ests[4]*(X[i,4])),matrix(c(ests[5],ests[6],ests[6],ests[7]),2,2)/2, upper=c(d[1],d[2]),algorithm="gibbs")
      rden <- GenTsnp(InitSNP,k,a,EXP,T[i],C[i],c(1,Xi[1,],X[i,4]))*dnorm(Xi[1,1],(ests[1]+ests[2]*X[i,4]+ests[6]/ests[7]*(Xi[1,2]-(ests[3]+ests[4]*X[i,4]))),sqrt(ests[5]-ests[6]^2/ests[7]))*dnorm(Xi[1,2],(ests[3]+ests[4]*X[i,4]),sqrt(ests[7]))
      for(t in 2:(50*M+200)) {
        Xstar <- Xi[(t-1),]+Genvals[l,]
        l <- l+1
        rnum <- GenTsnp(InitSNP,k,a,EXP,T[i],C[i],c(1,Xstar,X[i,4]))*dnorm(Xstar[1],(ests[1]+ests[2]*X[i,4]+ests[6]/ests[7]*(Xstar[2]-(ests[3]+ests[4]*X[i,4]))),sqrt(ests[5]-ests[6]^2/ests[7]))*dnorm(Xstar[2],(ests[3]+ests[4]*X[i,4]),sqrt(ests[7]))
        r <- rnum/max(rden,1e-300)
        u <- runif(1)
        if(r>u & Xstar[1]<d[1] & Xstar[2]<d[2]) {rden <- rnum
        Xi[t,] <- Xstar}  else Xi[t,] <- Xi[(t-1),]		 
      } 
      if(M>1) X12imp[,((i-(cc+mt+m1+mt1+m2+mt2))*2-1):((i-(cc+mt+m1+mt1+m2+mt2))*2)] <- Xi[sample(seq(250,(200+50*M),50)),]
      if(M==1) X12imp[,((i-(cc+mt+m1+mt1+m2+mt2))*2-1):((i-(cc+mt+m1+mt1+m2+mt2))*2)] <- Xi[seq(250,(200+50*M),50),]
    }}
  
  #Still need fixed for an EM program
  #TEM <- c(T[1:cc],as.vector(Timp),rep(T[(cc+mt+1):(cc+mt+m1)],each=M),as.vector(TX1imp[,seq(1,2*mt1,2)]),rep(T[(cc+mt+m1+mt1+1):(cc+mt+m1+mt1+m2)],each=M),as.vector(TX2imp[,seq(1,2*mt2,2)]),rep(T[(cc+mt+m1+mt1+m2+mt2+1):(cc+mt+m1+mt1+m2+mt2+m12)],each=M),as.vector(TX12imp[,seq(1,3*mt12,3)]))
  #XEM <- rbind(X[1:cc,],cbind(1,rbind(cbind(rep(X[(cc+1):(cc+mt),2],each=M),rep(X[(cc+1):(cc+mt),3],each=M)),cbind(as.vector(X1imp),rep(X[(cc+mt+1):(cc+mt+m1),3],each=M)),cbind(as.vector(TX1imp[,seq(2,(mt1*2),2)]),rep(X[(cc+mt+m1+1):(cc+mt+m1+mt1),3],each=M)),cbind(rep(X[(cc+mt+m1+mt1+1):(cc+mt+m1+mt1+m2),2],each=M),as.vector(X2imp)),cbind(rep(X[(cc+mt+m1+mt1+m2+1):(cc+mt+m1+mt1+m2+mt2),2],each=M),as.vector(TX2imp[,seq(2,(mt2*2),2)])),cbind(as.vector(X12imp[,seq(1,(m12*2),2)]),as.vector(X12imp[,seq(2,(m12*2),2)])),cbind(as.vector(TX12imp[,seq(2,(mt12*3),3)]),as.vector(TX12imp[,seq(3,(mt12*3),3)]))),rep(X[(cc+1):n,4],each=M)))
  XEM <- TEM <- NULL
  
  Timputed <- NULL
  Cimputed <- NULL
  Ximputed <- NULL
  for(i in 1:M){
    Timputed <- cbind(Timputed,T)
    Cimputed <- cbind(Cimputed,C)
    Ximputed <- cbind(Ximputed,rbind(X[1:(cc+mt),],cbind(1,X1imp[i,],X[(cc+mt+1):(cc+mt+m1+mt1),3:4]),cbind(X[(cc+mt+m1+mt1+1):(cc+mt+m1+mt1+m2+mt2),1:2],X2imp[i,],X[(cc+mt+m1+mt1+1):(cc+mt+m1+mt1+m2+mt2),4]),cbind(1,X12imp[i,seq(1,(2*(m12+mt12)),2)],X12imp[i,seq(2,(2*(m12+mt12)),2)],X[(cc+mt+m1+mt1+m2+mt2+1):n,4])))
  } 
  
  if(M>1) {Timputed <- rbind(Timputed[1:cc,],Timputed[(cc+mt+1):(cc+mt+m1),],Timputed[(cc+mt+m1+mt1+1):(cc+mt+m1+mt1+m2),],Timputed[(cc+mt+m1+mt1+m2+mt2+1):(cc+mt+m1+mt1+m2+mt2+m12),],Timputed[(cc+1):(cc+mt),],Timputed[(cc+mt+m1+1):(cc+mt+m1+mt1),],Timputed[(cc+mt+m1+mt1+m2+1):(cc+mt+m1+mt1+m2+mt2),],Timputed[(cc+mt+m1+mt1+m2+mt2+m12+1):n,])
  Cimputed <- rbind(Cimputed[1:cc,],Cimputed[(cc+mt+1):(cc+mt+m1),],Cimputed[(cc+mt+m1+mt1+1):(cc+mt+m1+mt1+m2),],Cimputed[(cc+mt+m1+mt1+m2+mt2+1):(cc+mt+m1+mt1+m2+mt2+m12),],Cimputed[(cc+1):(cc+mt),],Cimputed[(cc+mt+m1+1):(cc+mt+m1+mt1),],Cimputed[(cc+mt+m1+mt1+m2+1):(cc+mt+m1+mt1+m2+mt2),],Cimputed[(cc+mt+m1+mt1+m2+mt2+m12+1):n,])}
  if(M==1) {Timputed <- c(T[1:cc],T[(cc+mt+1):(cc+mt+m1)],T[(cc+mt+m1+mt1+1):(cc+mt+m1+mt1+m2)],T[(cc+mt+m1+mt1+m2+mt2+1):(cc+mt+m1+mt1+m2+mt2+m12)],T[(cc+1):(cc+mt)],T[(cc+mt+m1+1):(cc+mt+m1+mt1)],T[(cc+mt+m1+mt1+m2+1):(cc+mt+m1+mt1+m2+mt2)],T[(cc+mt+m1+mt1+m2+mt2+m12+1):n])
  Cimputed <- c(C[1:cc],C[(cc+mt+1):(cc+mt+m1)],C[(cc+mt+m1+mt1+1):(cc+mt+m1+mt1+m2)],C[(cc+mt+m1+mt1+m2+mt2+1):(cc+mt+m1+mt1+m2+mt2+m12)],C[(cc+1):(cc+mt)],C[(cc+mt+m1+1):(cc+mt+m1+mt1)],C[(cc+mt+m1+mt1+m2+1):(cc+mt+m1+mt1+m2+mt2)],C[(cc+mt+m1+mt1+m2+mt2+m12+1):n])}
  Ximputed <- rbind(Ximputed[1:cc,],Ximputed[(cc+mt+1):(cc+mt+m1),],Ximputed[(cc+mt+m1+mt1+1):(cc+mt+m1+mt1+m2),],Ximputed[(cc+mt+m1+mt1+m2+mt2+1):(cc+mt+m1+mt1+m2+mt2+m12),],Ximputed[(cc+1):(cc+mt),],Ximputed[(cc+mt+m1+1):(cc+mt+m1+mt1),],Ximputed[(cc+mt+m1+mt1+m2+1):(cc+mt+m1+mt1+m2+mt2),],Ximputed[(cc+mt+m1+mt1+m2+mt2+m12+1):n,])
  
  list(TEM, XEM,Timputed,Cimputed, Ximputed)
}




#Obtains multiple imputations when assuming Exponential, Weibull, lognormal, or loglogistic AFT model - called
#within PARMI function
MHMI2 <- function(T,C,X,InitPar,ests,cc,mt,m1,mt1,m2,mt2,m12,mt12,dist,M) {
  Beta <- InitPar[1:4]
  s <- sd(T)
  
  X1imp <- matrix(0,M,(m1+mt1))
  if((m1+mt1)==0) X1imp<- NULL
  if((m1+mt1)>0){
    for(i in (cc+mt+1):(cc+mt+m1+mt1)){
      Xi <- rep(100+20*M)
      Xi[1] <- rtnorm(1,ests[1]+ests[2]*X[i,4],sqrt(ests[5])/2,upper=d[1])
      rden <- GenT(dist,InitPar,T[i],C[i],c(X[i,1],Xi[1],X[i,3:4]))*dnorm(Xi[1],((ests[1]+ests[2]*X[i,4])+ests[6]/ests[7]*(X[i,3]-(ests[3]+ests[4]*X[i,4]))),sqrt(ests[5]-ests[6]^2/ests[7]))
      for(t in 2:(100+20*M)) {
        Xstar <- rnorm(1,mean=Xi[(t-1)],sqrt(ests[5])/2)
        rnum <- GenT(dist,InitPar,T[i],C[i],c(X[i,1],Xstar,X[i,3:4]))*dnorm(Xstar,((ests[1]+ests[2]*X[i,4])+ests[6]/ests[7]*(X[i,3]-(ests[3]+ests[4]*X[i,4]))),sqrt(ests[5]-ests[6]^2/ests[7]))
        r <- rnum/max(rden,1e-300)
        u <- runif(1)
        if(r>u & Xstar<d[1]) {rden <- rnum
        Xi[t] <- Xstar}   else Xi[t] <- Xi[(t-1)]
        
      }
      X1imp[,(i-(cc+mt))] <- sample(Xi[seq(120,(100+20*M),20)])
    }}
  
  X2imp <- matrix(0,M,(m2+mt2))
  if((m2+mt2)==0) X2imp<- NULL
  if((m2+mt2)>0){
    for(i in (cc+mt+m1+mt1+1):(cc+mt+m1+mt1+m2+mt2)){
      Xi <- rep(100+20*M)
      Xi[1] <- rtnorm(1,ests[3]+ests[4]*X[i,4],sqrt(ests[7])/2,upper=d[2])
      rden <- GenT(dist,InitPar,T[i],C[i],c(X[i,1:2],Xi[1],X[i,4]))*dnorm(Xi[1],((ests[3]+ests[4]*X[i,4])+ests[6]/ests[5]*(X[i,2]-(ests[1]+ests[2]*X[i,4]))),sqrt(ests[7]-ests[6]^2/ests[5]))
      for(t in 2:(100+20*M)) {
        Xstar <- rnorm(1,mean=Xi[(t-1)],sqrt(ests[7])/2)
        rnum <- GenT(dist,InitPar,T[i],C[i],c(X[i,1:2],Xstar,X[i,4]))*dnorm(Xstar,((ests[3]+ests[4]*X[i,4])+ests[6]/ests[5]*(X[i,2]-(ests[1]+ests[2]*X[i,4]))),sqrt(ests[7]-ests[6]^2/ests[5]))
        r <- rnum/max(rden,1e-300)
        u <- runif(1)
        if(r>u & Xstar<d[2]) {rden <- rnum
        Xi[t] <- Xstar} else Xi[t] <- Xi[(t-1)]
      }
      X2imp[,(i-(cc+mt+m1+mt1))] <- sample(Xi[seq(120,(100+20*M),20)])
    }}
  
  Genvals <- rmvnorm((50*M+200)*(m12+mt12),c(0,0),make.positive.definite(matrix(c(ests[5]/4,ests[6]/4,ests[6]/4,ests[7]/4),2,2),tol=0.001))
  l <- 1
  X12imp <- matrix(0,M,(m12+mt12)*2)
  if((m12+mt12)==0) X12imp<- NULL
  if((m12+mt12)>0){
    for(i in (cc+mt+m1+mt1+m2+mt2+1):(cc+mt+m1+mt1+m2+mt2+m12+mt12)){
      Xi <- matrix(0,(50*M+200),2)
      Xi[1,] <- rtmvnorm(1,mean=c(ests[1]+ests[2]*(X[i,4]),ests[3]+ests[4]*(X[i,4])),matrix(c(ests[5],ests[6],ests[6],ests[7]),2,2)/2, upper=c(d[1],d[2]),algorithm="gibbs")
      rden <- GenT(dist,InitPar,T[i],C[i],c(X[i,1],Xi[1,],X[i,4]))*dnorm(Xi[1,1],(ests[1]+ests[2]*X[i,4]+ests[6]/ests[7]*(Xi[1,2]-(ests[3]+ests[4]*X[i,4]))),sqrt(ests[5]-ests[6]^2/ests[7]))*dnorm(Xi[1,2],(ests[3]+ests[4]*X[i,4]),sqrt(ests[7]))
      for(t in 2:(50*M+200)) {
        Xstar <- Xi[(t-1),]+Genvals[l,]
        l <- l+1
        rnum <- GenT(dist,InitPar,T[i],C[i],c(X[i,1],Xstar,X[i,4]))*dnorm(Xstar[1],(ests[1]+ests[2]*X[i,4]+ests[6]/ests[7]*(Xstar[2]-(ests[3]+ests[4]*X[i,4]))),sqrt(ests[5]-ests[6]^2/ests[7]))*dnorm(Xstar[2],(ests[3]+ests[4]*X[i,4]),sqrt(ests[7]))
        r <- rnum/max(rden,1e-300)
        u <- runif(1)
        if(r>u & Xstar[1]<d[1] & Xstar[2]<d[2]) {rden <- rnum
        Xi[t,] <- Xstar}  else Xi[t,] <- Xi[(t-1),]		 
      } 
      X12imp[,((i-(cc+mt+m1+mt1+m2+mt2))*2-1):((i-(cc+mt+m1+mt1+m2+mt2))*2)] <- Xi[sample(seq(250,(200+50*M),50)),]
    }}
  
  #Still need fixed for an EM program
  #TEM <- c(T[1:cc],as.vector(Timp),rep(T[(cc+mt+1):(cc+mt+m1)],each=M),as.vector(TX1imp[,seq(1,2*mt1,2)]),rep(T[(cc+mt+m1+mt1+1):(cc+mt+m1+mt1+m2)],each=M),as.vector(TX2imp[,seq(1,2*mt2,2)]),rep(T[(cc+mt+m1+mt1+m2+mt2+1):(cc+mt+m1+mt1+m2+mt2+m12)],each=M),as.vector(TX12imp[,seq(1,3*mt12,3)]))
  #XEM <- rbind(X[1:cc,],cbind(1,rbind(cbind(rep(X[(cc+1):(cc+mt),2],each=M),rep(X[(cc+1):(cc+mt),3],each=M)),cbind(as.vector(X1imp),rep(X[(cc+mt+1):(cc+mt+m1),3],each=M)),cbind(as.vector(TX1imp[,seq(2,(mt1*2),2)]),rep(X[(cc+mt+m1+1):(cc+mt+m1+mt1),3],each=M)),cbind(rep(X[(cc+mt+m1+mt1+1):(cc+mt+m1+mt1+m2),2],each=M),as.vector(X2imp)),cbind(rep(X[(cc+mt+m1+mt1+m2+1):(cc+mt+m1+mt1+m2+mt2),2],each=M),as.vector(TX2imp[,seq(2,(mt2*2),2)])),cbind(as.vector(X12imp[,seq(1,(m12*2),2)]),as.vector(X12imp[,seq(2,(m12*2),2)])),cbind(as.vector(TX12imp[,seq(2,(mt12*3),3)]),as.vector(TX12imp[,seq(3,(mt12*3),3)]))),rep(X[(cc+1):n,4],each=M)))
  XEM <- TEM <- NULL
  
  Timputed <- NULL
  Cimputed <- NULL
  Ximputed <- NULL
  for(i in 1:M){
    Timputed <- cbind(Timputed,T)
    Cimputed <- cbind(Cimputed,C)
    Ximputed <- cbind(Ximputed,rbind(X[1:(cc+mt),],cbind(1,X1imp[i,],X[(cc+mt+1):(cc+mt+m1+mt1),3:4]),cbind(X[(cc+mt+m1+mt1+1):(cc+mt+m1+mt1+m2+mt2),1:2],X2imp[i,],X[(cc+mt+m1+mt1+1):(cc+mt+m1+mt1+m2+mt2),4]),cbind(1,X12imp[i,seq(1,(2*(m12+mt12)),2)],X12imp[i,seq(2,(2*(m12+mt12)),2)],X[(cc+mt+m1+mt1+m2+mt2+1):n,4])))
  } 
  
  Timputed <- rbind(Timputed[1:cc,],Timputed[(cc+mt+1):(cc+mt+m1),],Timputed[(cc+mt+m1+mt1+1):(cc+mt+m1+mt1+m2),],Timputed[(cc+mt+m1+mt1+m2+mt2+1):(cc+mt+m1+mt1+m2+mt2+m12),],Timputed[(cc+1):(cc+mt),],Timputed[(cc+mt+m1+1):(cc+mt+m1+mt1),],Timputed[(cc+mt+m1+mt1+m2+1):(cc+mt+m1+mt1+m2+mt2),],Timputed[(cc+mt+m1+mt1+m2+mt2+m12+1):n,])
  Cimputed <- rbind(Cimputed[1:cc,],Cimputed[(cc+mt+1):(cc+mt+m1),],Cimputed[(cc+mt+m1+mt1+1):(cc+mt+m1+mt1+m2),],Cimputed[(cc+mt+m1+mt1+m2+mt2+1):(cc+mt+m1+mt1+m2+mt2+m12),],Cimputed[(cc+1):(cc+mt),],Cimputed[(cc+mt+m1+1):(cc+mt+m1+mt1),],Cimputed[(cc+mt+m1+mt1+m2+1):(cc+mt+m1+mt1+m2+mt2),],Cimputed[(cc+mt+m1+mt1+m2+mt2+m12+1):n,])
  Ximputed <- rbind(Ximputed[1:cc,],Ximputed[(cc+mt+1):(cc+mt+m1),],Ximputed[(cc+mt+m1+mt1+1):(cc+mt+m1+mt1+m2),],Ximputed[(cc+mt+m1+mt1+m2+mt2+1):(cc+mt+m1+mt1+m2+mt2+m12),],Ximputed[(cc+1):(cc+mt),],Ximputed[(cc+mt+m1+1):(cc+mt+m1+mt1),],Ximputed[(cc+mt+m1+mt1+m2+1):(cc+mt+m1+mt1+m2+mt2),],Ximputed[(cc+mt+m1+mt1+m2+mt2+m12+1):n,])
  
  list(TEM, XEM,Timputed,Cimputed, Ximputed)
}




#Finds complete case parameter estimates for SNP distribution
SNPchoice <- function(T,C,X,divs,cc,mt){
  
  PhiVals <- seq(-1.50,1.50,3/(divs-1)) #grid points to be evaluated over range of phi
  Ins <- optim(c(Beta,1),CensNorm(T,C,X,cc,mt),method="BFGS")$par[1:5] #Initial-Initial Values SNP-EXP (still need updating)
  
  #k=0,EXP (Ins[1:4] is initial estimate for Beta, and Ins[5] is initial estimate for variation for normal SNP)
  K0exp <- optim(c(Ins[1:4],0.01),CensSNPexp(T,C,X,k=0,cc,mt),method="BFGS")
  
  #k=0,NORM
  K0norm <-optim(Ins,CensSNPnorm(T,C,X,k=0,cc,mt),method="BFGS")
  
  #k=1, EXP
  Initial <- InitVals(T,C,X,cc,mt,Ins,PhiVals,k=1,EXP=1)
  K1exp <- Optimize(Initial,T,C,X,k=1,cc,mt,EXP=1)
  
  #k=1, NORM
  Initial <- InitVals(T,C,X,cc,mt,Ins,PhiVals,k=1,EXP=0)
  K1norm <- Optimize(Initial,T,C,X,k=1,cc,mt,EXP=0)
  
  #k=2, EXP
  Initial <- InitVals(T,C,X,cc,mt,Ins,PhiVals,k=2,EXP=1)
  K2exp <- Optimize(Initial,T,C,X,k=2,cc,mt,EXP=1)
  
  #k=2, NORM
  Initial <- InitVals(T,C,X,cc,mt,Ins,PhiVals,k=2,EXP=0)
  K2norm <- Optimize(Initial,T,C,X,k=2,cc,mt,EXP=0)
  
  minim <-min(BIC(-K0exp$value,k=0),BIC(-K0norm$value,k=0),BIC(-K1exp$value,k=1),BIC(-K1norm$value,k=1),BIC(-K2exp$value,k=2),BIC(-K2norm$value,k=2))
  if(minim==BIC(-K0exp$value,0)){ Est <- K0exp
  EXP=TRUE
  K <- 0}
  if(minim==BIC(-K0norm$value,0)) {Est <- K0norm
  EXP=FALSE
  K <- 0}
  if(minim==BIC(-K1exp$value,1)) {Est <- K1exp
  EXP=TRUE
  K <- 1}
  if(minim==BIC(-K1norm$value,1)) {Est <- K1norm
  EXP=FALSE
  K <- 1}
  if(minim==BIC(-K2exp$value,2)) {Est <- K2exp
  EXP=TRUE
  K <- 2}
  if(minim==BIC(-K2norm$value,2)) {Est <- K2norm
  EXP=FALSE
  K <- 2}
  
  if(EXP==1) InitVar <- solve(optim(Est$par,CensSNPexp(T,C,X,K,cc,mt),method="BFGS",hessian=TRUE,control=list(maxit=1))$hessian) else InitVar <- solve(optim(Est$par,CensSNPnorm(T,C,X,K,cc,mt),method="BFGS",hessian=TRUE,control=list(maxit=1))$hessian)
  
  return(list(Est$par, InitVar, K, EXP))
}




#Puts data in order by CCs, censored T2, censored X1s, censored X2s
cendata <-function(X1,T1,C1){
  Xcc <- X1[(T1<C1 & X1[,2]>d[1] & X1[,3]>d[2]),]
  Xmt <- X1[(T1>C1 & X1[,2]>d[1] & X1[,3]>d[2]),]
  Xm1 <- X1[(T1<C1 & X1[,2]<d[1] & X1[,3]>d[2]),]
  Xm2 <- X1[(T1<C1 & X1[,2]>d[1] & X1[,3]<d[2]),]
  Xmt1 <- X1[(T1>C1 & X1[,2]<d[1] & X1[,3]>d[2]),]
  Xmt2 <- X1[(T1>C1 & X1[,2]>d[1] & X1[,3]<d[2]),]
  Xm12<- X1[(T1<C1 & X1[,2]<d[1] & X1[,3]<d[2]),]
  Xmt12 <- X1[(T1>C1 & X1[,2]<d[1] & X1[,3]<d[2]),]
  X <- rbind(Xcc,Xmt,Xm1,Xmt1,Xm2,Xmt2,Xm12,Xmt12)
  Tcc <- T1[(T1<C1 & X1[,2]>d[1] & X1[,3]>d[2])]
  Tmt <- T1[(T1>C1 & X1[,2]>d[1] & X1[,3]>d[2])]
  Tm1 <- T1[(T1<C1 & X1[,2]<d[1] & X1[,3]>d[2])]
  Tm2 <- T1[(T1<C1 & X1[,2]>d[1] & X1[,3]<d[2])]
  Tmt1 <- T1[(T1>C1 & X1[,2]<d[1] & X1[,3]>d[2])]
  Tmt2 <- T1[(T1>C1 & X1[,2]>d[1] & X1[,3]<d[2])]
  Tm12<- T1[(T1<C1 & X1[,2]<d[1] & X1[,3]<d[2])]
  Tmt12 <- T1[(T1>C1 & X1[,2]<d[1] & X1[,3]<d[2])]
  T <- c(Tcc,Tmt,Tm1,Tmt1,Tm2,Tmt2,Tm12,Tmt12)
  Ccc <- C1[(T1<C1 & X1[,2]>d[1] & X1[,3]>d[2])]
  Cmt <- C1[(T1>C1 & X1[,2]>d[1] & X1[,3]>d[2])]
  Cm1 <- C1[(T1<C1 & X1[,2]<d[1] & X1[,3]>d[2])]
  Cm2 <- C1[(T1<C1 & X1[,2]>d[1] & X1[,3]<d[2])]
  Cmt1 <- C1[(T1>C1 & X1[,2]<d[1] & X1[,3]>d[2])]
  Cmt2 <- C1[(T1>C1 & X1[,2]>d[1] & X1[,3]<d[2])]
  Cm12<- C1[(T1<C1 & X1[,2]<d[1] & X1[,3]<d[2])]
  Cmt12 <- C1[(T1>C1 & X1[,2]<d[1] & X1[,3]<d[2])]
  C <- c(Ccc,Cmt,Cm1,Cmt1,Cm2,Cmt2,Cm12,Cmt12)
  m <- c(length(Tcc),length(Tmt),length(Tm1),length(Tmt1),length(Tm2),length(Tmt2),length(Tm12),length(Tmt12))
  list(T,C,X,m)
}




#negative Log-likelihood function for censored bivariation normal (for maximizing covariate distribution parameters)
CensBN <- function(X,mn,m1,m2,mb) {
  like<-function(t) {
    fc<-rep(0,n)
    if(mn>0) fc[1:mn] <- -log(dnorm(X[1:mn,3],(t[3]+t[4]*X[1:mn,4]),sqrt(max(t[7],0.001)))*dnorm(X[1:mn,2],((t[1]+t[2]*X[1:mn,4])+t[6]/t[7]*(X[1:mn,3]-(t[3]+t[4]*X[1:mn,4]))),sqrt(max((t[5]-t[6]^2/t[7]),0.001))))
    if(m1>0) fc[(mn+1):(mn+m1)] <- -log(dnorm(X[(mn+1):(mn+m1),3],((t[3]+t[4]*X[(mn+1):(mn+m1),4])),sqrt(max(t[7],0.001)))*pnorm(d[1],((t[1]+t[2]*X[(mn+1):(mn+m1),4])+t[6]/t[7]*(X[(mn+1):(mn+m1),3]-(t[3]+t[4]*X[(mn+1):(mn+m1),4]))),sqrt(max((t[5]-t[6]^2/t[7]),0.001))))
    if(m2>0) fc[(mn+m1+1):(mn+m1+m2)] <- -log(dnorm(X[(mn+m1+1):(mn+m1+m2),2],((t[1]+t[2]*X[(mn+m1+1):(mn+m1+m2),4])),sqrt(max(t[5],0.001)))*pnorm(d[2],((t[3]+t[4]*X[(mn+m1+1):(mn+m1+m2),4])+t[6]/t[5]*(X[(mn+m1+1):(mn+m1+m2),2]-(t[1]+t[2]*X[(mn+m1+1):(mn+m1+m2),4]))),sqrt(max((t[7]-t[6]^2/t[5]),0.001))))	
    if(mb>0) for(i in (mn+m1+m2+1):n) fc[i] <- -log(pmvnorm(lower=c(-Inf,-Inf),upper=c(d[1],d[2]),mean=c((t[1]+t[2]*X[i,4]),(t[3]+t[4]*X[i,4])),sigma=matrix(c(max(t[5],0.001),t[6],t[6],max(t[7],0.001)),2,2)))
    like <- sum(fc)
  }
  like
}




#Log Likelihood function for censored bivariate normal (used in maximization for initial values) - called within
#SNPchoice function
CensNorm <- function(T,C,X,cnone,cone) {
  like<-function(t) {
    fc<-rep(0,cnone+cone)
    if(cnone>0) fc[1:cnone] <- -log(1/(t[5])*dnorm((T[1:cnone]-as.real(X[1:cnone,]%*%c(t[1],t[2],t[3],t[4])))/t[5]))
    if(cone>0) fc[(cnone+1):(cnone+cone)] <- -log(pnorm((C[(cnone+1):(cnone+cone)]-as.real(X[(cnone+1):(cnone+cone),]%*%c(t[1],t[2],t[3],t[4])))/t[5],lower.tail=FALSE))
    like <- sum(fc)
  }
  like
}




#Log Likelihood function for censored SNP-exponential (used in maximization) - called within the SNPchoice 
#and Optimize functions
CensSNPexp <- function(T,C,X,k,cnone,cone) {
  like<-function(t) {
    fc<-rep(0,(cnone+cone))
    a <- acoef(t[6:length(t)],k,EXP=TRUE)
    P_k <- Pk(T[1:cnone],X[1:cnone,],t,k,a,cnone,EXP=TRUE)
    a <- c(a, rep(0, 3-length(a))) 
    fc[1:cnone] <- t[5]-(T[1:cnone]-as.vector(X[1:cnone,]%*%c(t[1],t[2],t[3],t[4])))/exp(t[5])-log((P_k[1:cnone])^2)+exp((T[1:cnone]-as.vector(X[1:cnone,]%*%c(t[1],t[2],t[3],t[4])))/exp(t[5]))
    if(cone>0)  fc[(cnone+1):(cnone+cone)] <- -log(pSNPexp(exp((C[(cnone+1):(cnone+cone)]-as.vector(X[(cnone+1):(cnone+cone),]%*%c(t[1],t[2],t[3],t[4])))/exp(t[5])),a))			    
    like <- sum(fc)
    if(is.nan(like)==TRUE | is.na(like)==TRUE ) like <- 10000 else if(abs(like)>1e100) like <- 10000 
    return(like)
  }
  like
}




#Log Likelihood function for censored SNP-normal (used in maximization) - called within the SNPchoice 
#and Optimize functions
CensSNPnorm <- function(T,C,X,k,cnone,cone) {
  like<-function(t) {
    fc<-rep(0,(cnone+cone))
    a <- acoef(t[6:length(t)],k,EXP=FALSE)
    P_k <- Pk(T[1:cnone],X[1:cnone,],t,k,a,cnone,EXP=FALSE)
    a <- c(a, rep(0, 3-length(a))) 
    if(cnone>0) fc[1:cnone] <-  -log(1/(t[5]*exp(T[1:cnone]))*(P_k[1:cnone])^2*dnorm(((T[1:cnone]-as.vector(X[1:cnone,]%*%c(t[1],t[2],t[3],t[4])))/t[5])))
    if(cone>0)  fc[(cnone+1):(cnone+cone)] <- -log(pSNPnorm(((C[(cnone+1):(cnone+cone)]-as.vector(X[(cnone+1):(cnone+cone),]%*%c(t[1],t[2],t[3],t[4])))/t[5]),a))			    
    like <- sum(fc)
    if(is.nan(like)==TRUE | is.na(like)==TRUE ) like <- 10000 else if(abs(like)>1e100) like <- 10000 
    return(like)
  }
  like
}




#Calculates 3 to 5 sets of initial values for K=1,2 which are then used as starting points in maximization - called
#within the SNPMI and SNPchoice functions 
InitVals <- function(T,C,X,cc,mt,Ins,PhiVals,k,EXP,FINAL=FALSE) {
  PhiValsM <- PhiI <- SigmaI <- Beta0I <- NULL
  divs <- length(PhiVals)
  if(k>0) {for(i in 1:k) PhiValsM <- cbind(rep(rep(PhiVals,each=length(PhiVals)^(i-1)),length(PhiVals)^(k-i)),PhiValsM)
  Grid <- array(0,rep(length(PhiVals),k))
  Grid <- array(apply(PhiValsM,1,IV,T1=T,C1=C,X1=X,Ins,k,EXP,cc,mt,FINAL),rep(length(PhiVals),k))}
  if(k==1) for(i in 1:3){ y <- as.vector(which(Grid==min(Grid), arr.in=TRUE))
  Grid[(max(1,y[1]-round(divs/6,0)):min(length(PhiVals),y[1]+round(divs/6,0)))] <- max(Grid) 
  PhiI <- c(PhiI,PhiVals[y[1]])
  Beta0I <- c(Beta0I,IV2(PhiI[i],T,C,X,Ins,k,EXP,cc,mt,FINAL)[[1]])
  SigmaI <- c(SigmaI,IV2(PhiI[i],T,C,X,Ins,k,EXP,cc,mt,FINAL)[[2]])
  }
  if(k==2) for(i in 1:5){ y <- as.vector(which(Grid==min(Grid), arr.in=TRUE))
  Grid[(max(1,y[1]-round(divs/6,0)):min(length(PhiVals),y[1]+round(divs/6,0))),(max(1,y[2]-round(divs/6,0)):min(length(PhiVals),y[2]+round(divs/6,0)))] <- max(Grid)
  PhiI <- rbind(PhiI,c(PhiVals[y[2]],PhiVals[y[1]]))
  Beta0I <- c(Beta0I,IV2(PhiI[i,],T,C,X,Ins,k,EXP,cc,mt,FINAL)[[1]])
  SigmaI <- c(SigmaI,IV2(PhiI[i,],T,C,X,Ins,k,EXP,cc,mt,FINAL)[[2]])
  }
  if(k==0){Beta0I <- Ins[1]
  SigmaI <- Ins[5]}
  return(cbind(Beta0I,matrix(Ins[2:4],(1+k*2),3,byrow=TRUE),SigmaI,PhiI))
}




#Optimization function for SNP given data and initial values - called within the SNPMI and SNPchoice functions
Optimize <- function(InitSNP,T,C,X,k,cc,mt,EXP){
  if(EXP==TRUE){
    if(k==0) opt <- optim(InitSNP,CensSNPexp(T,C,X,k,cc,mt),method="BFGS")
    if(k==1) {
      for(i in 1:dim(InitSNP)[[1]]){
        possopt <- optim(InitSNP[i,],CensSNPexp(T,C,X,k,cc,mt),method="BFGS", control=list(maxit=1000))
        if(i==1) opt <- possopt else{ if(possopt$value<opt$value & possopt$par[6]<pi/2 & possopt$par[6]>-pi/2) opt <- possopt}		
      }
    }	
    if(k==2) {
      for(i in 1:dim(InitSNP)[[1]]){
        possopt <- optim(InitSNP[i,],CensSNPexp(T,C,X,k,cc,mt),method="BFGS", control=list(maxit=1000))
        if(i==1) opt <- possopt else{ if(possopt$value<opt$value & possopt$par[6]<pi/2 & possopt$par[6]>-pi/2& possopt$par[7]< pi/2 & possopt$par[7]> -pi/2) opt <- possopt}		
      }
    }
    if(k==1 & (opt$par[6]>pi/2 | opt$par[6]< -pi/2)) opt <- optim(InitSNP[1,],CensSNPexp(T,C,X,k,cc,mt),lower = c(-Inf,-Inf,-Inf,-Inf,-Inf,-pi/2), upper = c(Inf,Inf,Inf,Inf,Inf,pi/2),method="L-BFGS-B", control=list(maxit=1000))
    if(k==2 & (opt$par[6]>pi/2 | opt$par[6]< -pi/2 | opt$par[7]>pi/2 | opt$par[7]< -pi/2)) opt <- optim(InitSNP[1,],CensSNPexp(T,C,X,k,cc,mt),lower = c(-Inf,-Inf,-Inf,-Inf,-Inf,-pi/2,-pi/2), upper = c(Inf,Inf,Inf,Inf,Inf,pi/2,pi/2),method="L-BFGS-B",control=list(maxit=1000))
  }
  if(EXP==FALSE){
    if(k==0) opt <- optim(InitSNP,CensSNPnorm(T,C,X,k,cc,mt),method="BFGS")
    if(k==1) {
      for(i in 1:dim(InitSNP)[[1]]){
        possopt <- optim(InitSNP[i,],CensSNPnorm(T,C,X,k,cc,mt),method="BFGS", control=list(maxit=1000))
        if(i==1) opt <- possopt else{ if(possopt$value<opt$value  & possopt$par[6]<pi/2 & possopt$par[6]>-pi/2) opt <- possopt}		
      }
    }	
    if(k==2) {
      for(i in 1:dim(InitSNP)[[1]]){
        possopt <- optim(InitSNP[i,],CensSNPnorm(T,C,X,k,cc,mt),method="BFGS", control=list(maxit=1000))
        if(i==1) opt <- possopt else{ if(possopt$value<opt$value & possopt$par[6]<pi/2 & possopt$par[6]>-pi/2 & possopt$par[7]<pi/2 & possopt$par[7]>-pi/2) opt <- possopt}		
      }
    }
    if(k==1 & (opt$par[6]>pi/2 | opt$par[6]< -pi/2)) opt <- optim(InitSNP[1,],CensSNPnorm(T,C,X,k,cc,mt),lower = c(-Inf,-Inf,-Inf,-Inf,-Inf,-pi/2), upper = c(Inf,Inf,Inf,Inf,Inf,pi/2),method="L-BFGS-B", control=list(maxit=1000))
    if(k==2 & (opt$par[6]>pi/2 | opt$par[6]< -pi/2 | opt$par[7]>pi/2 | opt$par[7]< -pi/2)) opt <- optim(InitSNP[1,],CensSNPnorm(T,C,X,k,cc,mt),lower = c(-Inf,-Inf,-Inf,-Inf,-Inf,-pi/2,-pi/2), upper = c(Inf,Inf,Inf,Inf,Inf,pi/2,pi/2),method="L-BFGS-B",control=list(maxit=1000))
  }
  return(opt)
}




#Calculates approximate Standard Error for MI - called within the SNPMI and PARMI function
ApproxVar <- function(beta,ests,I,F,M) {
  par <-rbind(beta,ests)
  BW <- matrix(0,nrow(par),nrow(par))
  for(i in 1:M) BW <- BW + (par[,i]-rowMeans(par))%*%t(par[,i]-rowMeans(par))/(M-1)
  (F + (1+1/M)*BW+BW%*%solve(F)%*%I%*%solve(F)%*%BW)
}




#Calculates consistent standard error for MI with normal kernel SNP - can be called within the SNPMI function
VarNorm <- function(Timp,Cimp,Ximp,cnone,cone,par,phi,k,M,I,F) {
  
  ScoreMatrix <- matrix(0,dim(par)[[1]],M*n)
  
  for(d in 1:M){
    X <- Ximp[,((d-1)*4+1):(d*4)]
    T <- as.vector(Timp[,d])
    C <- as.vector(Cimp[,d])
    beta <- as.vector(par[1:4,d])
    sigma <- as.numeric(par[5,d])
    if(k>0) a <- acoef(phi[,d],k,EXP=FALSE) else a <- c(1,0,0)
    ests <- as.vector(par[6:dim(par)[[1]],d])
    
    S <- (T[1:cnone]-as.vector(X[1:cnone,]%*%beta))/sigma
    PK <- (a[1]+a[2]*S+a[3]*S^2)
    V <- (X[,2]-ests[1]-ests[2]*X[,4]-ests[6]/ests[7]*(X[,3]-ests[3]-ests[4]*X[,4]))
    Sig <- (ests[5]-ests[6]^2/ests[7])
    Num <- (X[,3]-ests[3]-ests[4]*X[,4])
    
    for(j in 1:4) ScoreMatrix[j,(d*n+1-n):(d*n-cone)] <- -2/PK*(a[2]/sigma*X[1:cnone,j]+2*a[3]*S*X[1:cnone,j]/sigma)+S*X[1:cnone,j]/sigma
    ScoreMatrix[5,(d*n+1-n):(d*n-cone)] <- -1/sigma - 2/PK*(a[2]*S/sigma+2*a[3]*S^2/sigma)+S^2/sigma
    ScoreMatrix[6,(d*n+1-n):(d*n)] <- V/Sig
    ScoreMatrix[7,(d*n+1-n):(d*n)] <- X[,4]*V/Sig
    ScoreMatrix[8,(d*n+1-n):(d*n)] <- -V/Sig*ests[6]^2/ests[7]+Num/ests[7]
    ScoreMatrix[9,(d*n+1-n):(d*n)] <- -V/Sig*X[,4]*ests[6]^2/ests[7]+X[,4]*Num/ests[7]
    ScoreMatrix[10,(d*n+1-n):(d*n)] <- -1/(2*Sig)+V^2/(2*Sig^2)
    ScoreMatrix[11,(d*n+1-n):(d*n)] <- -ests[6]^2/ests[7]^2/(2*Sig)-1/(2*ests[7])-ests[6]/ests[7]^2*V/Sig+V^2*ests[6]^2/ests[7]^2/(2*Sig^2)+Num^2/(2*ests[7]^2)
    ScoreMatrix[12,(d*n+1-n):(d*n)] <- ests[6]/ests[7]/Sig + V*1/ests[7]^2/Sig-V^2*ests[6]/ests[7]/Sig^2
    
    S2 <- (C[(cnone+1):(cnone+cone)]-as.vector(X[(cnone+1):(cnone+cone),]%*%beta))/sigma
    for(j in 1:4) ScoreMatrix[j,(d*n+1+cnone-n):(d*n)] <- -SNPnorm(S2,a)*X[(cnone+1):(cnone+cone),j]/sigma
    ScoreMatrix[5,(d*n+1+cnone-n):(d*n)] <- -SNPnorm(S2,a)*S2/sigma
  }
  
  IFmIN <- matrix(0,12,12)
  for(i in 1:n){
    Means <- rowMeans(ScoreMatrix[,seq(i,M*n,n)])
    for(d in 1:M){
      IFmIN <- IFmIN + (ScoreMatrix[,(n*(d-1)+i)]-Means)%*%t((ScoreMatrix[,(n*(d-1)+i)]-Means)) 
    }
  }
  
  IFmIN <- IFmIN/(M-1)
  F + (1+1/M)*F%*%IFmIN%*%F + F%*%IFmIN%*%I%*%IFmIN%*%F
}




#Calculates consistent standard error for MI with exponential kernel SNP - can be called within the SNPMI function
VarExp <- function(Timp,Cimp,Ximp,cnone,cone,par,phi,k,M,I,F) {
  
  ScoreMatrix <- matrix(0,dim(par)[[1]],M*n)
  
  for(d in 1:M){
    X <- Ximp[,((d-1)*4+1):(d*4)]
    T <- as.vector(Timp[,d])
    C <- as.vector(Cimp[,d])
    beta <- as.vector(par[1:4,d])
    sigma <- as.numeric(par[5,d])
    if(k>0) a <-acoef(phi[,d],k,EXP=TRUE) else a <- c(1,0,0)
    ests <- as.vector(par[6:dim(par)[[1]],d])
    
    S <- exp((T[1:cnone]-as.vector(X[1:cnone,]%*%beta))/exp(sigma))
    SB <- ((T[1:cnone]-as.vector(X[1:cnone,]%*%beta))/exp(2*sigma))
    PK <- (a[1]+a[2]*S+a[3]*S^2)
    V <- (X[,2]-ests[1]-ests[2]*X[,4]-ests[6]/ests[7]*(X[,3]-ests[3]-ests[4]*X[,4]))
    Sig <- (ests[5]-ests[6]^2/ests[7])
    Num <- (X[,3]-ests[3]-ests[4]*X[,4])
    
    for(j in 1:4) ScoreMatrix[j,(d*n+1-n):(d*n-cone)] <- -X[1:cnone,j]/exp(sigma)+S*X[1:cnone,j]/exp(sigma)-2/PK*(a[2]*S*X[1:cnone,j]/exp(sigma)+a[3]*S^2*2*X[1:cnone,j]/exp(sigma))
    ScoreMatrix[5,(d*n+1-n):(d*n-cone)] <- -1 - SB+S*SB - 2/PK*(a[2]*S*SB+2*a[3]*S^2*SB)
    ScoreMatrix[6,(d*n+1-n):(d*n)] <- V/Sig
    ScoreMatrix[7,(d*n+1-n):(d*n)] <- X[,4]*V/Sig
    ScoreMatrix[8,(d*n+1-n):(d*n)] <- -V/Sig*ests[6]^2/ests[7]+Num/ests[7]
    ScoreMatrix[9,(d*n+1-n):(d*n)] <- -V/Sig*X[,4]*ests[6]^2/ests[7]+X[,4]*Num/ests[7]
    ScoreMatrix[10,(d*n+1-n):(d*n)] <- -1/(2*Sig)+V^2/(2*Sig^2)
    ScoreMatrix[11,(d*n+1-n):(d*n)] <- -ests[6]^2/ests[7]^2/(2*Sig)-1/(2*ests[7])-ests[6]/ests[7]^2*V/Sig+V^2*ests[6]^2/ests[7]^2/(2*Sig^2)+Num^2/(2*ests[7]^2)
    ScoreMatrix[12,(d*n+1-n):(d*n)] <- ests[6]/ests[7]/Sig + V*1/ests[7]^2/Sig-V^2*ests[6]/ests[7]/Sig^2
    
    S2 <- exp((C[(cnone+1):(cnone+cone)]-as.vector(X[(cnone+1):(cnone+cone),]%*%beta))/exp(sigma))
    S2B <- ((T[(cnone+1):(cnone+cone)]-as.vector(X[(cnone+1):(cnone+cone),]%*%beta))/exp(2*sigma))
    for(j in 1:4) ScoreMatrix[j,(d*n+1+cnone-n):(d*n)] <- -SNPexp(S2,a)*S2*X[(cnone+1):(cnone+cone),j]/exp(sigma)
    ScoreMatrix[5,(d*n+1+cnone-n):(d*n)] <- -SNPnorm(S2,a)*S2*S2B
  }
  
  IFmIN <- matrix(0,12,12)
  for(i in 1:n){
    Means <- rowMeans(ScoreMatrix[,seq(i,M*n,n)])
    for(d in 1:M){
      IFmIN <- IFmIN + (ScoreMatrix[,(n*(d-1)+i)]-Means)%*%t((ScoreMatrix[,(n*(d-1)+i)]-Means)) 
    }
  }
  
  IFmIN <- IFmIN/(M-1)
  F + (1+1/M)*F%*%IFmIN%*%F + F%*%IFmIN%*%I%*%IFmIN%*%F
}



#############Additional Computational Functions###################
##(minor functions called within other functions for computational details)##

#Gets c coefficients, then a coefficients (See Zhang and Davidian 2001)
acoef <- function(Phi,k,EXP=TRUE) {
  if(EXP==TRUE) {B <- chol(Aexp[1:(k+1),1:(k+1)])
  Binv <- solve(B)}
  if(EXP==FALSE) {B <- chol(Anorm[1:(k+1),1:(k+1)])
  Binv <- solve(B)}
  c <- rep(1,k+1)
  if(k>0) { for(i in 2:(k+1)) c[i] <- c[i-1]*cos(Phi[i-1])
  for(i in 1:(k)) c[i] <- c[i]*sin(Phi[i])
  }   		   
  as.real(Binv%*%c)
}

#Finds P_K values
Pk <- function(T,X,t,k,a,cc,EXP=TRUE) {
  s <- matrix(0,cc,(k+1))
  if(EXP==TRUE) for(i in 1:(k+1)) s[,i] <- (exp((T-as.real(X%*%c(t[1],t[2],t[3],t[4])))/exp(t[5])))^(i-1)
  if(EXP==FALSE) for(i in 1:(k+1)) s[,i] <- ((T-as.real(X%*%c(t[1],t[2],t[3],t[4])))/t[5])^(i-1)
  s%*%a
}

#Calculates AIC or BIC 
AIC <- function(logLike,k) -2*logLike+2*(k+5)+ 2*(k+5)*(k+6)/(n-k-5-1)
BIC <- function(logLike,k) -2*logLike+(k+5)*log(n)

#pdf functions used in GenT below (needed for the MI algorithm)
EV <- function(x,loc=0,scale=1) {z=(x-loc)/scale
1/scale*exp(z-exp(z))}
logistic <- function(x,loc=0,scale=1) {z=(x-loc)/scale
1/scale*exp(z)/(1+exp(z))^2}

#Used in MH-algorithm for imputing values when assuming an exponential, weibull, lognormal or log-logistic model for T
GenT <- function(dist,InitPar,T,C,X){
  Beta <- InitPar[1:4]
  Rest <- InitPar[5:length(InitPar)]
  if(dist==1){if(T<C) gen <- EV(T,loc=X%*%Beta) #standard extreme value (exponential model for T)
  if(T>C) gen <- pgumbel(-C,loc=-X%*%Beta)}
  if(dist==2){if(T<C)  gen <- EV(T,loc=X%*%Beta,scale=Rest) #Gumbel/log-weibull  (weibull model for T)
  if(T>C) gen <- pgumbel(-C,loc=-X%*%Beta,scale=Rest)}
  if(dist==3){if(T<C) gen <- dnorm(T,mean=X%*%Beta,Rest)   #normal (log-normal model for T)
  if(T>C) gen <-pnorm(C,mean=X%*%Beta,Rest) }
  if(dist==4){if(T<C)  gen <- logistic(T,loc=X%*%Beta,scale=Rest) #logistic (log-logisitic model for T)
  if(T>C) gen <-dlogis(C,location=X%*%Beta,Rest) }
  
  gen
}

#CDF fucionts used in GenTsnp below (needed for MI algorithm)
pSNPexp <- function(x,a) a[1]^2*pgamma(x,1,1,lower.tail=FALSE)+2*a[1]*a[2]*pgamma(x,2,1,lower.tail=FALSE)+2*(a[2]^2+2*a[1]*a[3])*pgamma(x,3,1,lower.tail=FALSE)+12*a[2]*a[3]*pgamma(x,4,1,lower.tail=FALSE)+24*a[3]^2*pgamma(x,5,1,lower.tail=FALSE)
pSNPnorm <-function(x,a) a[1]^2*pnorm(x,lower.tail=FALSE)+2*a[1]*a[2]*dnorm(x)+(a[2]^2+2*a[1]*a[3])*(x*dnorm(x)+pnorm(x,lower.tail=FALSE))+2*a[2]*a[3]*(x^2+2)*dnorm(x)+a[3]^2*(x^3*dnorm(x)+3*x*dnorm(x)+3*pnorm(x,lower.tail=FALSE))

#Used in MH-algorithm for imputing values when assuming an SNP model
GenTsnp <- function(InitSNP,k,a,EXP=TRUE,T,C,X){
  Beta <- InitSNP[1:4]
  Sigma <- InitSNP[5]
  if(k==0) a1 <- c(a,0,0)
  if(k==1) a1 <- c(a,0)
  if(k==2) a1 <- a
  if(EXP==TRUE){if(T<C) gen <- 1/(exp(Sigma))*exp((T-as.real(X%*%Beta))/exp(Sigma))*(Pk(T,X,InitSNP,k,a,1))^2*dexp(exp((T-as.real(X%*%Beta))/exp(Sigma)))
  if(T>C) gen <- pSNPexp(min(exp((C-as.real(X%*%Beta))/exp(Sigma)),1e50),a1)
  }
  if(EXP==FALSE){if(T<C) gen <- 1/Sigma*(Pk(T,X,InitSNP,k,a,1,EXP=FALSE))^2*dnorm((T-as.real(X%*%Beta))/Sigma)
  if(T>C) gen <- pSNPnorm(min(((C-as.real(X%*%Beta))/Sigma),1e50),a1)						    
  }
  gen
}

#Functions for finding moments of various quantities
SNPexpX <- function(x,a) x*dexp(x)*(a[1]+a[2]*x+a[3]*x^2)^2
SNPnormX <- function(x,a) x*dnorm(x)*(a[1]+a[2]*x+a[3]*x^2)^2
SNPlog <- function(x,a) log(x)*dexp(x)*(a[1]+a[2]*x+a[3]*x^2)^2
SNPlog2 <- function(x,a) (log(x))^2*dexp(x)*(a[1]+a[2]*x+a[3]*x^2)^2
SNPx <- function(x,a) x*dnorm(x)*(a[1]+a[2]*x+a[3]*x^2)^2
SNPx2 <- function(x,a) x^2*dnorm(x)*(a[1]+a[2]*x+a[3]*x^2)^2

#Pdf of SNP distributions 
SNPexp <- function(x,a) dexp(x)*(a[1]+a[2]*x+a[3]*x^2)^2
SNPnorm <- function(x,a) dnorm(x)*(a[1]+a[2]*x+a[3]*x^2)^2


#Finds initial Values for Intecepts
Beta0If <- function(Phi,Sigma,Beta,a,EXP=TRUE,FINAL=FALSE,a2=NULL) {
  if(FINAL==FALSE){
    if(EXP==TRUE) Beta0I <- Beta[1]-exp(Sigma)*integrate(SNPlog,lower=0,upper=Inf,stop.on.error=FALSE,a)$value
    if(EXP==FALSE) Beta0I <- Beta[1]-Sigma*integrate(SNPx,lower=-Inf,upper=Inf,stop.on.error=FALSE,a)$value
  }
  if(FINAL==TRUE){
    if(EXP==TRUE) Beta0I <- Beta[1]+exp(Sigma)*integrate(SNPlog,lower=0,upper=Inf,stop.on.error=FALSE,a2)$value-exp(Sigma)*integrate(SNPlog,lower=0,upper=Inf,stop.on.error=FALSE,a)$value
    if(EXP==FALSE) Beta0I <- Beta[1]+Sigma*integrate(SNPx,lower=-Inf,upper=Inf,stop.on.error=FALSE,a2)$value-Sigma*integrate(SNPx,lower=-Inf,upper=Inf,stop.on.error=FALSE,a)$value
  }
  
  Beta0I
}

#Calculates likelihood Values at Initial Values
IV <- function(Phi,T1,C1,X1,Ins,k,EXP=TRUE,cc,mt,FINAL=FALSE) {
  BetaI <- Ins[1:4]
  TauI <- Ins[5]
  a <- c(acoef(Phi,k,EXP), rep(0, 2-k)) 
  if(FINAL==FALSE){
    if(EXP==TRUE) SigmaI <- log(TauI/sqrt(integrate(SNPlog2,lower=0,upper=Inf,stop.on.error=FALSE,a)$value-(integrate(SNPlog,lower=0,upper=Inf,stop.on.error=FALSE,a)$value)^2))
    if(EXP==FALSE) SigmaI <- TauI/sqrt(integrate(SNPx2,lower=-Inf,upper=Inf,stop.on.error=FALSE,a)$value-(integrate(SNPx,lower=-Inf,upper=Inf,stop.on.error=FALSE,a)$value)^2)
    Beta0I <- Beta0If(Phi,SigmaI,BetaI,a,EXP,FINAL)
  }
  if(FINAL==TRUE){EstPhi <- Ins[6:length(Ins)]
  a2 <- c(acoef(EstPhi,k,EXP), rep(0, 2-k)) 
  if(EXP==TRUE) SigmaI <- log(exp(TauI)*sqrt((integrate(SNPlog2,lower=0,upper=Inf,stop.on.error=FALSE,a2)$value-(integrate(SNPlog,lower=0,upper=Inf,stop.on.error=FALSE,a2)$value)^2)/(integrate(SNPlog2,lower=0,upper=Inf,stop.on.error=FALSE,a)$value-(integrate(SNPlog,lower=0,upper=Inf,stop.on.error=FALSE,a)$value)^2)))
  if(EXP==FALSE) SigmaI <- (TauI*sqrt((integrate(SNPx2,lower=-Inf,upper=Inf,stop.on.error=FALSE,a2)$value-(integrate(SNPx,lower=-Inf,upper=Inf,stop.on.error=FALSE,a2)$value)^2)/(integrate(SNPx2,lower=-Inf,upper=Inf,stop.on.error=FALSE,a)$value-(integrate(SNPx,lower=-Inf,upper=Inf,stop.on.error=FALSE,a)$value)^2)))
  Beta0I <- Beta0If(Phi,SigmaI,BetaI,a,EXP,FINAL,a2)
  }
  Init(c(Beta0I,BetaI[2:4],SigmaI,Phi),T1,C1,X1,k,cc,mt,EXP)
}

#Calculates likelihood Values at Initial Values and outputs only updated initial values
IV2 <- function(Phi,T1,C1,X1,Ins,k,EXP=TRUE,cc,mt,FINAL=FALSE) {
  BetaI <- Ins[1:4]
  TauI <- Ins[5]
  a <- c(acoef(Phi,k,EXP), rep(0, 2-k)) 
  if(FINAL==FALSE){
    if(EXP==TRUE) SigmaI <- log(TauI/sqrt(integrate(SNPlog2,lower=0,upper=Inf,stop.on.error=FALSE,a)$value-(integrate(SNPlog,lower=0,upper=Inf,stop.on.error=FALSE,a)$value)^2))
    if(EXP==FALSE) SigmaI <- TauI/sqrt(integrate(SNPx2,lower=-Inf,upper=Inf,stop.on.error=FALSE,a)$value-(integrate(SNPx,lower=-Inf,upper=Inf,stop.on.error=FALSE,a)$value)^2)
    Beta0I <- Beta0If(Phi,SigmaI,BetaI,a,EXP,FINAL)
  }
  if(FINAL==TRUE){EstPhi <- Ins[6:length(Ins)]
  a2 <- c(acoef(EstPhi,k,EXP), rep(0, 2-k)) 
  if(EXP==TRUE) SigmaI <- log(exp(TauI)*sqrt((integrate(SNPlog2,lower=0,upper=Inf,stop.on.error=FALSE,a2)$value-(integrate(SNPlog,lower=0,upper=Inf,stop.on.error=FALSE,a2)$value)^2)/(integrate(SNPlog2,lower=0,upper=Inf,stop.on.error=FALSE,a)$value-(integrate(SNPlog,lower=0,upper=Inf,stop.on.error=FALSE,a)$value)^2)))
  if(EXP==FALSE) SigmaI <- (TauI*sqrt((integrate(SNPx2,lower=-Inf,upper=Inf,stop.on.error=FALSE,a2)$value-(integrate(SNPx,lower=-Inf,upper=Inf,stop.on.error=FALSE,a2)$value)^2)/(integrate(SNPx2,lower=-Inf,upper=Inf,stop.on.error=FALSE,a)$value-(integrate(SNPx,lower=-Inf,upper=Inf,stop.on.error=FALSE,a)$value)^2)))
  Beta0I <- Beta0If(Phi,SigmaI,BetaI,a,EXP,FINAL,a2)
  }
  list(Beta0I,SigmaI)
}


#Outputs likelihood value for each grid point of phi which is being evaluated using function InitVals
Init <- function(t,T,C,X,k,cnone,cone,EXP=TRUE) {
  fc<-rep(0,(cnone+cone))
  if(EXP==TRUE){a <- acoef(t[6:length(t)],k,EXP=TRUE)
  P_k <- Pk(T[1:cnone],X[1:cnone,],t,k,a,cnone,EXP=TRUE)
  a <- c(a, rep(0, 3-length(a))) 
  if(cnone>0) fc[1:cnone] <- t[5]-(T[1:cnone]-as.vector(X[1:cnone,]%*%c(t[1],t[2],t[3],t[4])))/exp(t[5])-log((P_k[1:cnone])^2)+exp((T[1:cnone]-as.vector(X[1:cnone,]%*%c(t[1],t[2],t[3],t[4])))/exp(t[5]))
  if(cone>0) for(i in (cnone+1):(cnone+cone)) fc[i] <- -log(tryCatch(integrate(SNPexp,lower=min(exp((C[i]-as.real(X[i,]%*%c(t[1],t[2],t[3],t[4])))/exp(t[5])),1e50),upper=Inf,stop.on.error=FALSE,rel.tol = .Machine$double.eps^0.5,a)$value, error=function(...) 1))}		    
  if(EXP==FALSE){a <- acoef(t[6:length(t)],k,EXP=FALSE)
  P_k <- Pk(T[1:cnone],X[1:cnone,],t,k,a,cnone,EXP=FALSE)
  a <- c(a, rep(0, 3-length(a))) 
  if(cnone>0) fc[1:cnone] <-  -log(1/t[5]*(P_k[1:cnone])^2*dnorm(((T[1:cnone]-as.vector(X[1:cnone,]%*%c(t[1],t[2],t[3],t[4])))/t[5])))
  if(cone>0) for(i in (cnone+1):(cnone+cone)) fc[i] <- -log(tryCatch(integrate(SNPnorm,lower=min(((C[i]-as.real(X[i,]%*%c(t[1],t[2],t[3],t[4])))/t[5]),1e50),upper=Inf,stop.on.error=FALSE,rel.tol = .Machine$double.eps^0.5,a)$value, error=function(...) 1))}					    
  sum(fc)
}
