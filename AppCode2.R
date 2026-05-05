###################Necessary Packages########################
library(mnormt)
library(tmvtnorm)
library(msm)
library(corpcor)
library(Matrix)
library(evd)
library(survival)

set.seed(8089205)
####################Data Read and Clean-Up####################Reading in  GenIMS longitudinal data
#Reading in  GenIMS longitudinal data from laptop
GenIMS <- read.csv("C://Users/Paul/Documents/Data.csv")

#Reading in GenIMS single time indepdent data from laptop
GenIMSInd <- read.csv("C://Users/Paul/Documents/Data2.csv")

#Attaching sex,race and age covariates from time independent data and restricting to first day measurements
GenIMS1 <- cbind(GenIMS[GenIMS[,2]==1,],GenIMSInd[,c(4,5,69)])

#Removing patients who did not actually have CAP 
GenIMS1 <- GenIMS1[GenIMSInd[,68]==1,]

#Removing patients who were immediately discharged
GenIMS1 <- GenIMS1[GenIMS1[,66]=="INPATIENT",]

#Defining Response (survival time and survival indicator)
TnC <- as.numeric(GenIMS1[,60])
Dead <- as.numeric(GenIMS1[,56])-1    #>90 day death status GenIMS1[,59] 

#Defining Covariate matrix: covariate 1 is TNF, covariate 2 is IL-6, covariate 3 is IL-10, covariate 4 is sex (1=male), covariate 5 is race, covariate 6 is age
X <- cbind(1,as.numeric(GenIMS1[,22]),as.numeric(GenIMS1[,23]),as.numeric(GenIMS1[,24]),as.numeric(GenIMS1[,82]-1),as.numeric(GenIMS1[,83]),as.numeric(GenIMS1[,84]))

#Eliminating those without first day observations (should not affect analysis greatly since 
#lack of first day observation occured only when patient arrived on a weekend 
ResCov <- cbind(TnC,Dead,X)
ResCov <- na.omit(ResCov)
TnC <- ResCov[,1]
Dead <- ResCov[,2]
Dead[Dead==1][which(TnC[Dead==1]>90)]<-0
X <- ResCov[,3:9]

#Defining n=number of individuals (1418) and censoring values after log-transformation; note that for 16 individuals, the censoring value for the 
#IL-6 covariate was 2 instead of 5;  this is taken into account in the methods below
n <- length(ResCov[,1])
d <- c(log(4), log(5), log(5))


#Taking log tranformation of covariates to make them approximately normal; defining race as a binary variabl
for(i in 1:n){
  if(X[i,2]>exp(d[1])) X[i,2] <- log(X[i,2])
  if(X[i,3]>exp(d[2])) X[i,3] <- log(X[i,3])
  if(X[i,4]>exp(d[3])) X[i,4] <- log(X[i,4])
  if(X[i,6]>1) X[i,6] <- 0
}

#creating an n x 3 matrix of the detection limits for convenience below (all rows are the same except for those individuals with censoring at 2 for IL-6
twocens <- which(X[,3]==-2)  #Finds which observations for IL-6 are censored at 2 rather than 5 (5 is much more common); none occur with just IL-6 censored
dmatrix <- matrix(d,n,3,byrow=TRUE)
dmatrix[twocens,2] <- log(2)

###########Parameters and Matrix Initializations##############
M <- c(15,15,30,50)
Iters <- length(M)
divs <- 31
Aexp <- matrix(c(1,1,2,1,2,6,2,6,24),3,3)
Anorm <- matrix(c(1,0,1,0,1,0,1,0,3),3,3)

######################Analysis of Data########################
CD <- Censdata(X,TnC,Dead)  #Ordering censored data according to censoring
X <- CD[[1]]
T <- CD[[2]]
Cind <- CD[[3]]
cc <- CD[[4]][1] #no censoring
mt <- CD[[4]][2] #no censoring
m1 <- CD[[4]][3] #TNF censored
mt1 <- CD[[4]][4] #TNF censored
m2 <- CD[[4]][5] #IL-6 censored
mt2 <- CD[[4]][6] #IL-6 censored
m3 <- CD[[4]][7] #IL-10 censored
mt3 <- CD[[4]][8] #IL-10 censored
m12 <- CD[[4]][9] #TNF and IL-6 censored
mt12 <- CD[[4]][10] #TNF and IL-6 censored
m13 <- CD[[4]][11] #TNF and IL-10 censored
mt13 <- CD[[4]][12] #TNF and IL-10 censored
m23 <- CD[[4]][13] #IL-6 and IL-10 censored
mt23 <- CD[[4]][14] #IL-6 and IL-10 censored
m123 <- CD[[4]][15] #TNF, IL-6, and IL-10 censored
mt123 <- CD[[4]][16] #TNF, IL-6, and IL-10 censored
C <- (1-Cind)*log(90) #(1-Cind)*log(T) for C>90
T <- log(T)

#maximizes censored trivariate normal
Ests <- optim(c(1.7,0,0,0,3.6,0,0,0,1.7,0,0,0, 1.09, 0.96, 0.56, 5.28, 1.80, 2.88),Maxing(X,cc+mt,m1+mt1,m2+mt2,m3+mt3,m12+mt12,m13+mt13,m23+mt23,m123+mt123,dmatrix),hessian=TRUE,method="BFGS",control=list(trace=3,reltol=1e-10,maxit=1000))
InitBN <- Ests$par  #Initial Parameter Estimates in trivariate normal
InitVarBN <- solve(Ests$hessian) #Initial variance estimates

##################################################
############## Analysis using SNP ################
##################################################
SNPch <- SNPchoice(T,C,X,divs,cc,mt) #Find complete case parameter estimates for SNP for multiple imputation
InitSNP <- SNPch[[1]]
InitVar <- as.matrix(bdiag(SNPch[[2]][1:8,1:8], InitVarBN))
K <- SNPch[[3]]
EXPd <- SNPch[[4]]
SNPmi <- SNPMI(T,C,X,cc,mt,m1,mt1,m2,mt2,m3,mt3,m12,mt12,m13,mt13,m23,mt23,m123,mt123,InitSNP,InitBN,InitVar,K,EXPd,M)
FinalSNP1 <- SNPmi[[1]][1,]
FinalSNP2 <- SNPmi[[1]][2,]
FinalSNP3 <- SNPmi[[1]][3,]
FinalSNP4 <- SNPmi[[1]][4,]
FinalBN <- rowMeans(SNPmi[[3]])
SEsnp1_a <- SNPmi[[2]][1,]
SEsnp2_a <- SNPmi[[2]][2,]
SEsnp3_a <- SNPmi[[2]][3,]
SEsnp4_a <- SNPmi[[2]][4,]



##################################################
########### Analysis using Parametric ############
##################################################
Parmi <- PARMI(T,C,X,cc,mt,m1,mt1,m2,mt2,m3,mt3,m12,mt12,m13,mt13,m23,mt23,m123,mt123,InitBN,InitVarBN,M,Weib=FALSE)
FinalPar1 <- Parmi[[1]][1,]
FinalPar2 <- Parmi[[1]][2,]
FinalPar3 <- Parmi[[1]][3,]
FinalPar4 <- Parmi[[1]][4,]
FinalBN2 <- rowMeans(Parmi[[3]])
SEPar1 <- Parmi[[2]][1,]
SEPar2 <- Parmi[[2]][2,]
SEPar3 <- Parmi[[2]][3,]
SEPar4 <- Parmi[[2]][4,]


###################Computational Functions#####################

#Creating ordered data matrix (no censoring - all three covariates censored); uses 0 as cut-off since 
#censored values are coded as -4,-2/-5,-5 and logs of all covariates here are positive
Censdata <-function(X1,T1,C1){
  Xcc <- X1[(C1==1 & X1[,2]>0 & X1[,3]>0 & X1[,4]>0),]
  Xmt <- X1[(C1==0 & X1[,2]>0 & X1[,3]>0 & X1[,4]>0),]
  Xm1 <- X1[(C1==1 & X1[,2]<0 & X1[,3]>0 & X1[,4]>0),]
  Xmt1 <- X1[(C1==0 & X1[,2]<0 & X1[,3]>0 & X1[,4]>0),]
  Xm2 <- X1[(C1==1 & X1[,2]>0 & X1[,3]<0 & X1[,4]>0),]
  Xmt2 <- X1[(C1==0 & X1[,2]>0 & X1[,3]<0 & X1[,4]>0),]
  Xm3 <- X1[(C1==1 & X1[,2]>0 & X1[,3]>0 & X1[,4]<0),]
  Xmt3 <- X1[(C1==0 & X1[,2]>0 & X1[,3]>0 & X1[,4]<0),]
  Xm12 <- X1[(C1==1 & X1[,2]<0 & X1[,3]<0 & X1[,4]>0),]
  Xmt12 <- X1[(C1==0 & X1[,2]<0 & X1[,3]<0 & X1[,4]>0),]
  Xm13 <- X1[(C1==1 & X1[,2]<0 & X1[,3]>0 & X1[,4]<0),]
  Xmt13 <- X1[(C1==0 & X1[,2]<0 & X1[,3]>0 & X1[,4]<0),]
  Xm23 <- X1[(C1==1 & X1[,2]>0 & X1[,3]<0 & X1[,4]<0),]
  Xmt23 <- X1[(C1==0 & X1[,2]>0 & X1[,3]<0 & X1[,4]<0),]
  Xm123 <- X1[(C1==1 & X1[,2]<0 & X1[,3]<0 & X1[,4]<0),]
  Xmt123 <- X1[(C1==0 & X1[,2]<0 & X1[,3]<0 & X1[,4]<0),]
  X <- rbind(Xcc,Xmt,Xm1,Xmt1,Xm2,Xmt2,Xm3,Xmt3,Xm12,Xmt12,Xm13,Xmt13,Xm23,Xmt23,Xm123,Xmt123)
  Tcc <- T1[(C1==1 & X1[,2]>0 & X1[,3]>0 & X1[,4]>0)]
  Tmt <- T1[(C1==0 & X1[,2]>0 & X1[,3]>0 & X1[,4]>0)]
  Tm1 <- T1[(C1==1 & X1[,2]<0 & X1[,3]>0 & X1[,4]>0)]
  Tmt1 <- T1[(C1==0 & X1[,2]<0 & X1[,3]>0 & X1[,4]>0)]
  Tm2 <- T1[(C1==1 & X1[,2]>0 & X1[,3]<0 & X1[,4]>0)]
  Tmt2 <- T1[(C1==0 & X1[,2]>0 & X1[,3]<0 & X1[,4]>0)]
  Tm3 <- T1[(C1==1 & X1[,2]>0 & X1[,3]>0 & X1[,4]<0)]
  Tmt3 <- T1[(C1==0 & X1[,2]>0 & X1[,3]>0 & X1[,4]<0)]
  Tm12 <- T1[(C1==1 & X1[,2]<0 & X1[,3]<0 & X1[,4]>0)]
  Tmt12 <- T1[(C1==0 & X1[,2]<0 & X1[,3]<0 & X1[,4]>0)]
  Tm13 <- T1[(C1==1 & X1[,2]<0 & X1[,3]>0 & X1[,4]<0)]
  Tmt13 <- T1[(C1==0 & X1[,2]<0 & X1[,3]>0 & X1[,4]<0)]
  Tm23 <- T1[(C1==1 & X1[,2]>0 & X1[,3]<0 & X1[,4]<0)]
  Tmt23 <- T1[(C1==0 & X1[,2]>0 & X1[,3]<0 & X1[,4]<0)]
  Tm123 <- T1[(C1==1 & X1[,2]<0 & X1[,3]<0 & X1[,4]<0)]
  Tmt123 <- T1[(C1==0 & X1[,2]<0 & X1[,3]<0 & X1[,4]<0)]
  T <- c(Tcc,Tmt,Tm1,Tmt1,Tm2,Tmt2,Tm3,Tmt3,Tm12,Tmt12,Tm13,Tmt13,Tm23,Tmt23,Tm123,Tmt123)
  Ccc <- C1[(C1==1 & X1[,2]>0 & X1[,3]>0 & X1[,4]>0)]
  Cmt <- C1[(C1==0 & X1[,2]>0 & X1[,3]>0 & X1[,4]>0)]
  Cm1 <- C1[(C1==1 & X1[,2]<0 & X1[,3]>0 & X1[,4]>0)]
  Cmt1 <- C1[(C1==0 & X1[,2]<0 & X1[,3]>0 & X1[,4]>0)]
  Cm2 <- C1[(C1==1 & X1[,2]>0 & X1[,3]<0 & X1[,4]>0)]
  Cmt2 <- C1[(C1==0 & X1[,2]>0 & X1[,3]<0 & X1[,4]>0)]
  Cm3 <- C1[(C1==1 & X1[,2]>0 & X1[,3]>0 & X1[,4]<0)]
  Cmt3 <- C1[(C1==0 & X1[,2]>0 & X1[,3]>0 & X1[,4]<0)]
  Cm12 <- C1[(C1==1 & X1[,2]<0 & X1[,3]<0 & X1[,4]>0)]
  Cmt12 <- C1[(C1==0 & X1[,2]<0 & X1[,3]<0 & X1[,4]>0)]
  Cm13 <- C1[(C1==1 & X1[,2]<0 & X1[,3]>0 & X1[,4]<0)]
  Cmt13 <- C1[(C1==0 & X1[,2]<0 & X1[,3]>0 & X1[,4]<0)]
  Cm23 <- C1[(C1==1 & X1[,2]>0 & X1[,3]<0 & X1[,4]<0)]
  Cmt23 <- C1[(C1==0 & X1[,2]>0 & X1[,3]<0 & X1[,4]<0)]
  Cm123 <- C1[(C1==1 & X1[,2]<0 & X1[,3]<0 & X1[,4]<0)]
  Cmt123 <- C1[(C1==0 & X1[,2]<0 & X1[,3]<0 & X1[,4]<0)]
  C <- c(Ccc,Cmt,Cm1,Cmt1,Cm2,Cmt2,Cm3,Cmt3,Cm12,Cmt12,Cm13,Cmt13,Cm23,Cmt23,Cm123,Cmt123)
  m <- c(length(Tcc),length(Tmt),length(Tm1),length(Tmt1),length(Tm2),length(Tmt2),length(Tm3),length(Tmt3),length(Tm12),length(Tmt12),length(Tm13),length(Tmt13),length(Tm23),length(Tmt23),length(Tm123),length(Tmt123))
  list(X,T,C,m)
}

#Log Likelihood function for censored trivariate normal with mean function Az (where A is 3 x 4 matrix of coefficients on intercept,sex,race,age)
Maxing <- function(X,mn,m1,m2,m3,m12,m13,m23,ma,dmatrix) {
  like<-function(t) {
    fc<-rep(0,n)
    sig123 <- matrix(c(t[13],t[14],t[15],t[14],t[16],t[17],t[15],t[17],t[18]),3,3)
    sig12 <- matrix(c(t[13],t[14],t[14],t[16]),2,2)
    sig13 <- matrix(c(t[13],t[15],t[15],t[18]),2,2)
    sig23 <- matrix(c(t[16],t[17],t[17],t[18]),2,2)
    if(mn>0) fc[1:mn] <- -log(dnorm(X[1:mn,2], mean=(t[1]+t[2]*X[1:mn,5]+t[3]*X[1:mn,6]+t[4]*X[1:mn,7]+as.real(matrix(c(t[14],t[15]),1,2)%*%solve(sig23)%*%t((X[1:mn,3:4]-matrix(c(t[5]+t[6]*X[1:mn,5]+t[7]*X[1:mn,6]+t[8]*X[1:mn,7],t[9]+t[10]*X[1:mn,5]+t[11]*X[1:mn,6]+t[12]*X[1:mn,7]),mn,2))))),sd=sqrt((t[13]-as.real(matrix(c(t[14],t[15]),1,2)%*%solve(sig23)%*%matrix(c(t[14],t[15]),2,1)))))*dnorm(X[1:mn,3],mean=(t[5]+t[6]*X[1:mn,5]+t[7]*X[1:mn,6]+t[8]*X[1:mn,7]+t[17]/t[18]*(X[1:mn,4]-t[9]-t[10]*X[1:mn,5]-t[11]*X[1:mn,6]-t[12]*X[1:mn,7])),sd=sqrt(t[16]-t[17]^2/t[18]))*dnorm(X[1:mn,4],mean=(t[9]+t[10]*X[1:mn,5]+t[11]*X[1:mn,6]+t[12]*X[1:mn,7]),sd=sqrt(t[18])))	
    if(m1>0) fc[(mn+1):(mn+m1)] <- -log(dnorm(X[(mn+1):(mn+m1),3],mean=(t[5]+t[6]*X[(mn+1):(mn+m1),5]+t[7]*X[(mn+1):(mn+m1),6]+t[8]*X[(mn+1):(mn+m1),7]+t[17]/t[18]*(X[(mn+1):(mn+m1),4]-t[9]-t[10]*X[(mn+1):(mn+m1),5]-t[11]*X[(mn+1):(mn+m1),6]-t[12]*X[(mn+1):(mn+m1),7])),sd=sqrt(t[16]-t[17]^2/t[18]))*dnorm(X[(mn+1):(mn+m1),4],mean=(t[9]+t[10]*X[(mn+1):(mn+m1),5]+t[11]*X[(mn+1):(mn+m1),6]+t[12]*X[(mn+1):(mn+m1),7]),sd=sqrt(t[18]))*pnorm(d[1],mean=(t[1]+t[2]*X[(mn+1):(mn+m1),5]+t[3]*X[(mn+1):(mn+m1),6]+t[4]*X[(mn+1):(mn+m1),7]+as.real(matrix(c(t[14],t[15]),1,2)%*%solve(sig23)%*%t((X[(mn+1):(mn+m1),3:4]-matrix(c(t[5]+t[6]*X[(mn+1):(mn+m1),5]+t[7]*X[(mn+1):(mn+m1),6]+t[8]*X[(mn+1):(mn+m1),7],t[9]+t[10]*X[(mn+1):(mn+m1),5]+t[11]*X[(mn+1):(mn+m1),6]+t[12]*X[(mn+1):(mn+m1),7]),(m1),2))))),sd=sqrt((t[13]-as.real(matrix(c(t[14],t[15]),1,2)%*%solve(sig23)%*%matrix(c(t[14],t[15]),2,1))))))
    if(m2>0) fc[(mn+m1+1):(mn+m1+m2)] <- -log(dnorm(X[(mn+m1+1):(mn+m1+m2),2],mean=(t[1]+t[2]*X[(mn+m1+1):(mn+m1+m2),5]+t[3]*X[(mn+m1+1):(mn+m1+m2),6]+t[4]*X[(mn+m1+1):(mn+m1+m2),7]+t[15]/t[18]*(X[(mn+m1+1):(mn+m1+m2),4]-t[9]-t[10]*X[(mn+m1+1):(mn+m1+m2),5]-t[11]*X[(mn+m1+1):(mn+m1+m2),6]-t[12]*X[(mn+m1+1):(mn+m1+m2),7])),sd=sqrt(t[13]-t[15]^2/t[18]))*dnorm(X[(mn+m1+1):(mn+m1+m2),4],mean=(t[9]+t[10]*X[(mn+m1+1):(mn+m1+m2),5]+t[11]*X[(mn+m1+1):(mn+m1+m2),6]+t[12]*X[(mn+m1+1):(mn+m1+m2),7]),sd=sqrt(t[18]))*pnorm(d[2],mean=(t[5]+t[6]*X[(mn+m1+1):(mn+m1+m2),5]+t[7]*X[(mn+m1+1):(mn+m1+m2),6]+t[8]*X[(mn+m1+1):(mn+m1+m2),7]+as.real(matrix(c(t[14],t[17]),1,2)%*%solve(sig13)%*%t((X[(mn+m1+1):(mn+m1+m2),c(2,4)]-matrix(c(t[1]+t[2]*X[(mn+m1+1):(mn+m1+m2),5]+t[3]*X[(mn+m1+1):(mn+m1+m2),6]+t[4]*X[(mn+m1+1):(mn+m1+m2),7],t[9]+t[10]*X[(mn+m1+1):(mn+m1+m2),5]+t[11]*X[(mn+m1+1):(mn+m1+m2),6]+t[12]*X[(mn+m1+1):(mn+m1+m2),7]),m2,2))))),sd=sqrt((t[16]-as.real(matrix(c(t[14],t[17]),1,2)%*%solve(sig13)%*%matrix(c(t[14],t[17]),2,1))))))
    if(m3>0) fc[(mn+m1+m2+1):(mn+m1+m2+m3)] <- -log(dnorm(X[(mn+m1+m2+1):(mn+m1+m2+m3),2],mean=(t[1]+t[2]*X[(mn+m1+m2+1):(mn+m1+m2+m3),5]+t[3]*X[(mn+m1+m2+1):(mn+m1+m2+m3),6]+t[4]*X[(mn+m1+m2+1):(mn+m1+m2+m3),7]+t[14]/t[16]*(X[(mn+m1+m2+1):(mn+m1+m2+m3),3]-t[5]-t[6]*X[(mn+m1+m2+1):(mn+m1+m2+m3),5]-t[7]*X[(mn+m1+m2+1):(mn+m1+m2+m3),6]-t[8]*X[(mn+m1+m2+1):(mn+m1+m2+m3),7])),sd=sqrt(t[13]-t[14]^2/t[16]))*dnorm(X[(mn+m1+m2+1):(mn+m1+m2+m3),3],mean=(t[5]+t[6]*X[(mn+m1+m2+1):(mn+m1+m2+m3),5]+t[7]*X[(mn+m1+m2+1):(mn+m1+m2+m3),6]+t[8]*X[(mn+m1+m2+1):(mn+m1+m2+m3),7]),sd=sqrt(t[16]))*pnorm(d[3],mean=(t[9]+t[10]*X[(mn+m1+m2+1):(mn+m1+m2+m3),5]+t[11]*X[(mn+m1+m2+1):(mn+m1+m2+m3),6]+t[12]*X[(mn+m1+m2+1):(mn+m1+m2+m3),7]+as.real(matrix(c(t[15],t[17]),1,2)%*%solve(sig12)%*%t((X[(mn+m1+m2+1):(mn+m1+m2+m3),2:3]-matrix(c(t[1]+t[2]*X[(mn+m1+m2+1):(mn+m1+m2+m3),5]+t[3]*X[(mn+m1+m2+1):(mn+m1+m2+m3),6]+t[4]*X[(mn+m1+m2+1):(mn+m1+m2+m3),7],t[5]+t[6]*X[(mn+m1+m2+1):(mn+m1+m2+m3),5]+t[7]*X[(mn+m1+m2+1):(mn+m1+m2+m3),6]+t[8]*X[(mn+m1+m2+1):(mn+m1+m2+m3),7]),m3,2))))),sd=sqrt((t[18]-as.real(matrix(c(t[15],t[17]),1,2)%*%solve(sig12)%*%matrix(c(t[15],t[17]),2,1))))))
    if(m12>0) fc[(mn+m1+m2+m3+1):(mn+m1+m2+m3+m12)] <- -log(dnorm(X[(mn+m1+m2+m3+1):(mn+m1+m2+m3+m12),4],mean=(t[9]+t[10]*X[(mn+m1+m2+m3+1):(mn+m1+m2+m3+m12),5]+t[11]*X[(mn+m1+m2+m3+1):(mn+m1+m2+m3+m12),6]+t[12]*X[(mn+m1+m2+m3+1):(mn+m1+m2+m3+m12),7]),sd=sqrt(t[18]))*apply(cbind(X[(mn+m1+m2+m3+1):(mn+m1+m2+m3+m12),],dmatrix[(mn+m1+m2+m3+1):(mn+m1+m2+m3+m12),1:2]),1,Int1,Lower=c(-Inf,-Inf),t=t,sig12=sig12))
    if(m13>0) fc[(mn+m1+m2+m3+m12+1):(mn+m1+m2+m3+m12+m13)] <- -log(dnorm(X[(mn+m1+m2+m3+m12+1):(mn+m1+m2+m3+m12+m13),3],mean=(t[5]+t[6]*X[(mn+m1+m2+m3+m12+1):(mn+m1+m2+m3+m12+m13),5]+t[7]*X[(mn+m1+m2+m3+m12+1):(mn+m1+m2+m3+m12+m13),6]+t[8]*X[(mn+m1+m2+m3+m12+1):(mn+m1+m2+m3+m12+m13),7]),sd=sqrt(20))*apply(X[(mn+m1+m2+m3+m12+1):(mn+m1+m2+m3+m12+m13),],1,Int2,Lower=c(-Inf,-Inf),Upper=c(d[1],d[3]),t=t,sig13=sig13))
    if(m23>0) fc[(mn+m1+m2+m3+m12+m13+1):(mn+m1+m2+m3+m12+m13+m23)] <- -log(dnorm(X[(mn+m1+m2+m3+m12+m13+1):(mn+m1+m2+m3+m12+m13+m23),2],mean=(t[1]+t[2]*X[(mn+m1+m2+m3+m12+m13+1):(mn+m1+m2+m3+m12+m13+m23),5]+t[3]*X[(mn+m1+m2+m3+m12+m13+1):(mn+m1+m2+m3+m12+m13+m23),6]+t[4]*X[(mn+m1+m2+m3+m12+m13+1):(mn+m1+m2+m3+m12+m13+m23),7]),sd=sqrt(t[13]))*apply(cbind(X[(mn+m1+m2+m3+m12+m13+1):(mn+m1+m2+m3+m12+m13+m23),],dmatrix[(mn+m1+m2+m3+m12+m13+1):(mn+m1+m2+m3+m12+m13+m23),2:3]),1,Int3,Lower=c(-Inf,-Inf),t=t,sig23=sig23))
    if(ma>0) fc[(mn+m1+m2+m3+m12+m13+m23+1):n] <- -log(apply(cbind(X[(mn+m1+m2+m3+m12+m13+m23+1):n,],dmatrix[(mn+m1+m2+m3+m12+m13+m23+1):n,]),1,Int4,t=t,sig123=sig123,Lower=c(-Inf,-Inf,-Inf)))
    like <- sum(fc)
  }
  like
}

SNPchoice <- function(T,C,X,divs,cc,mt){
  
  PhiVals <- seq(-1.50,1.50,3/(divs-1)) #grid points to be evaluated over range of phi
  Ins <- optim(c(8,-.1,-.05,0.002,-.25,.04,-.025,1),CensNorm(T,C,X,cc,mt),method="BFGS")$par[1:8]
  
  #k=0,EXP (Ins[1:7] is initial estimate for Beta, and Ins[8] is initial estimate for variation for normal SNP)
  K0exp <- optim(c(Ins[1:7],0.01),CensSNPexp(T,C,X,k=0,cc,mt),method="BFGS")
  
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
  EXPd=TRUE
  K <- 0}
  if(minim==BIC(-K0norm$value,0)) {Est <- K0norm
  EXPd=FALSE
  K <- 0}
  if(minim==BIC(-K1exp$value,1)) {Est <- K1exp
  EXPd=TRUE
  K <- 1}
  if(minim==BIC(-K1norm$value,1)) {Est <- K1norm
  EXPd=FALSE
  K <- 1}
  if(minim==BIC(-K2exp$value,2)) {Est <- K2exp
  EXPd=TRUE
  K <- 2}
  if(minim==BIC(-K2norm$value,2)) {Est <- K2norm
  EXPd=FALSE
  K <- 2}
  
  if(EXPd==1) InitVar<- solve(optim(Est$par,CensSNPexp(T,C,X,K,cc,mt),method="BFGS",hessian=TRUE,control=list(maxit=1))$hessian) else InitVar <- solve(optim(Est$par,CensSNPnorm(T,C,X,K,cc,mt),method="BFGS",hessian=TRUE,control=list(maxit=1))$hessian)	
  return(list(Est$par, InitVar, K, EXP))
}

#Obtains IMI (or simple MI) estimate for parameters using SNP model
SNPMI <- function(T,C,X,cc,mt,m1,mt1,m2,mt2,m3,mt3,m12,mt12,m13,mt13,m23,mt23,m123,mt123,InitSNP,InitBN,InitVar,K,EXP,M){
  Iters <- length(M)
  EstUpdates <- matrix(0,Iters,length(InitSNP))
  SEUpdates_a  <- matrix(0,Iters,8)
  
  for(i in 1:Iters){
    if(i==1){
      Initsnp <- InitSNP
      Initbn <- InitBN
    }
    
    if(i>1){
      Initsnp <- rowMeans(SNPmis)
      Initbn <- rowMeans(BNmis)
      InitVar <- VarNorm(Timped,Cimped,Ximped,(cc+m1+m2+m3+m12+m13+m23+m123),(mt+mt1+mt2+mt3+mt12+mt13+mt23+mt123),rbind(SNPmis[1:8,],BNmis),matrix(SNPmis2[9,],ncol=M[i-1]),K,M[i-1],InitVar,FinalVar)
    }
    BNmis <- matrix(0,18,M[i])
    FinalVar <-matrix(0,26,26)
    Int <- rep(0,M[i])
    SNPmis <- matrix(0,length(InitSNP), M[i])
    Full <- XimpsSNP(T,C,X,K,Initsnp,Initbn,cc,mt,m1,mt1,m2,mt2,m3,mt3,m12,mt12,m13,mt13,m23,mt23,m123,mt123,EXPd,M[i],d,twocens)
    Timped <- Full[[1]]
    Cimped <- Full[[2]]
    Ximped <- Full[[3]]
    Int <- rep(0,M[i])
    for(q in 1:M[i]){
      Tnew <- Timped[,q]
      Cnew <- Cimped[,q]
      Xnew <- Ximped[,(q*7-6):(q*7)]
      FINAL <- Optimize(InitVals(Tnew,Cnew,Xnew,(cc+m1+m2+m3+m12+m13+m23+m123),(mt+mt1+mt2+mt3+mt12+mt13+mt23+mt123),Initsnp,PhiVals,K,EXPd,FINAL=TRUE),Tnew,Cnew,Xnew,K,(cc+m1+m2+m3+m12+m13+m23+m123),(mt+mt1+mt2+mt3+mt12+mt13+mt23+mt123),EXPd)
      SNPmis[,q] <- FINAL$par
      if(EXPd==1) VAR <- nlm(CensSNPexp(Tnew,Cnew,Xnew,K,(cc+m1+m2+m3+m12+m13+m23+m123),(mt+mt1+mt2+mt3+mt12+mt13+mt23+mt123)),c(FINAL$par),hessian=TRUE,iterlim=1)$hessian else VAR <- nlm(CensSNPnorm(Tnew,Cnew,Xnew,K,(cc+m1+m2+m3+m12+m13+m23+m123),(mt+mt1+mt2+mt3+mt12+mt13+mt23+mt123)),c(FINAL$par),hessian=TRUE,iterlim=1)$hessian
      MAX2 <- optim(Initbn,Maxing(Xnew,n,0,0,0,0,0,0,0,dmatrix), method="BFGS", hessian=TRUE)
      BNmis[,q] <- MAX2$par
      FinalVar <- FinalVar + solve(as.matrix(bdiag(VAR[1:8,1:8],MAX2$hessian)))/M[i]
      Int[q] <- SNPmis2[1,q] + SNPmis2[8,q]*integrate(SNPnormX,lower=-Inf,upper=Inf,a=c(acoef(SNPmis[9:length(Est$par),q],K,EXP=EXPd),rep(0,2-K)))$value
    }
    
    EstUpdates[i,] <- rowMeans(SNPmis)
    EstUpdates[i,1] <- mean(Int)
    SEUpdates_a[i,] <- sqrt(diag(ApproxVar(SNPmis[1:8,],BNmis,InitVar,FinalVar,M[i])))[1:8]
  }
  
  return(list(EstUpdates,SEUpdates_a,BNmis))
}

#Obtains IMI (or MI) estimates using one of the typical parametric models
PARMI <- function(T,C,X,cc,mt,m1,mt1,m2,mt2,m3,mt3,m12,mt12,m13,mt13,m23,mt23,m123,mt123,InitBN,InitVarBN,M,Weib=FALSE){
  Iters <- length(M)
  InitPar <- rep(0,8)
  survreg.control(maxiter=500)
  ###Getting Complete Case Estimates###
  AFT1<-survreg(Surv(exp(c(T[1:cc],C[(cc+1):(cc+mt)])),c(rep(1,cc),rep(0,mt))) ~ X[1:(cc+mt),2] + X[1:(cc+mt),3] + X[1:(cc+mt),4]+ X[1:(cc+mt),5]+ X[1:(cc+mt),6]+ X[1:(cc+mt),7],dist="exponential")
  AFT2<-survreg(Surv(exp(c(T[1:cc],C[(cc+1):(cc+mt)])),c(rep(1,cc),rep(0,mt))) ~ X[1:(cc+mt),2] + X[1:(cc+mt),3] + X[1:(cc+mt),4]+ X[1:(cc+mt),5]+ X[1:(cc+mt),6]+ X[1:(cc+mt),7],dist="weibull")
  AFT3<-survreg(Surv(exp(c(T[1:cc],C[(cc+1):(cc+mt)])),c(rep(1,cc),rep(0,mt))) ~ X[1:(cc+mt),2] + X[1:(cc+mt),3] + X[1:(cc+mt),4]+ X[1:(cc+mt),5]+ X[1:(cc+mt),6]+ X[1:(cc+mt),7],dist="lognorm")
  AFT4<-survreg(Surv(exp(c(T[1:cc],C[(cc+1):(cc+mt)])),c(rep(1,cc),rep(0,mt))) ~ X[1:(cc+mt),2] + X[1:(cc+mt),3] + X[1:(cc+mt),4]+ X[1:(cc+mt),5]+ X[1:(cc+mt),6]+ X[1:(cc+mt),7],dist="loglogistic")
  Models <- list(AFT1,AFT2,AFT3,AFT4)
  Params <- c(-1,0,0,0)
  aic <- Inf
  test <- which(c(is.nan(Models[[1]]$loglik[2])==FALSE,is.nan(Models[[2]]$loglik[2])==FALSE,is.nan(Models[[3]]$loglik[2])==FALSE,is.nan(Models[[4]]$loglik[2])==FALSE)==TRUE)
  for(i in test){
    if(AIC(Models[[i]]$loglik[2],Params[i])<aic){
      aic <- AIC(Models[[i]]$loglik[2],Params[i])
      InitPar[1:7] <- Models[[i]]$coefficients
      InitPar[8] <- Models[[i]]$scale
      dists <- i
      InitVar <- as.matrix(Models[[i]]$var)
    }
  }
  
  if(Weib==TRUE){
    InitPar[1:7] <- Models[[2]]$coefficients
    InitPar[8] <- Models[[2]]$scale
    dists <- 2
    InitVar <- as.matrix(Models[[2]]$var)
  }
  
  InitVar <- as.matrix(bdiag(InitVar,InitVarBN))
  EstUpdates <- matrix(0,Iters,8)
  SEUpdates  <- matrix(0,Iters,8)
  
  for(i in 1:Iters){
    if(i==1){
      Initpar <- InitPar
      Initbn <- InitBN
    }
    if(i>1){
      Initpar <- rowMeans(Parmis)
      Initbn <- rowMeans(BNmis)
      InitVar <- ApproxVar(Parmis,BNmis,InitVar,FinalVar,M[(i-1)])
    }
    
    BNmis <- matrix(0,18,M[i])
    FinalVar <-matrix(0,26,26)
    Parmis <- matrix(0,8, M[i])
    Full <- XimpsPar(T,C,X,Initpar,Initbn,cc,mt,m1,mt1,m2,mt2,m3,mt3,m12,mt12,m13,mt13,m23,mt23,m123,mt123,dists,M[i],d,twocens)
    Timped <- Full[[1]]
    Cimped <- Full[[2]]
    Ximped <- Full[[3]]
    
    q <- 1
    while(q<=M[i]){
      Tnew <- Timped[,q]
      Cnew <- Cimped[,q]
      Xnew <- Ximped[,(q*7-6):(q*7)]
      if(dists==1)	FINAL <- survreg(Surv(exp(c(Tnew[1:(cc+m1+m2+m3+m12+m13+m23+m123)],Cnew[(cc+m1+m2+m3+m12+m13+m23+m123+1):n])),c(rep(1,(cc+m1+m2+m3+m12+m13+m23+m123)),rep(0,(mt+mt1+mt2+mt3+mt12+mt13+mt23+mt123)))) ~ Xnew[,2]+Xnew[,3]+ Xnew[,4] + Xnew[,5] + Xnew[,6] + Xnew[,7],dist="exponential")
      if(dists==2)	FINAL <- survreg(Surv(exp(c(Tnew[1:(cc+m1+m2+m3+m12+m13+m23+m123)],Cnew[(cc+m1+m2+m3+m12+m13+m23+m123+1):n])),c(rep(1,(cc+m1+m2+m3+m12+m13+m23+m123)),rep(0,(mt+mt1+mt2+mt3+mt12+mt13+mt23+mt123)))) ~ Xnew[,2]+Xnew[,3]+ Xnew[,4] + Xnew[,5] + Xnew[,6] + Xnew[,7],dist="weibull")
      if(dists==3)	FINAL <- survreg(Surv(exp(c(Tnew[1:(cc+m1+m2+m3+m12+m13+m23+m123)],Cnew[(cc+m1+m2+m3+m12+m13+m23+m123+1):n])),c(rep(1,(cc+m1+m2+m3+m12+m13+m23+m123)),rep(0,(mt+mt1+mt2+mt3+mt12+mt13+mt23+mt123))))~ Xnew[,2]+Xnew[,3]+ Xnew[,4] + Xnew[,5] + Xnew[,6] + Xnew[,7],dist="lognorm")
      if(dists==4)	FINAL <- survreg(Surv(exp(c(Tnew[1:(cc+m1+m2+m3+m12+m13+m23+m123)],Cnew[(cc+m1+m2+m3+m12+m13+m23+m123+1):n])),c(rep(1,(cc+m1+m2+m3+m12+m13+m23+m123)),rep(0,(mt+mt1+mt2+mt3+mt12+mt13+mt23+mt123)))) ~ Xnew[,2]+Xnew[,3]+ Xnew[,4] + Xnew[,5] + Xnew[,6] + Xnew[,7],dist="loglogistic")
      Parmis[,q] <- c(FINAL$coefficients,FINAL$scale)
      MAX2 <- optim(InitBN,Maxing(Xnew,n,0,0,0,0,0,0,0,dmatrix),method="BFGS",hessian=TRUE)
      BNmis[,q] <-as.vector(MAX2$par)
      if(sum(is.nan(FINAL$var))==0 & sum(is.na(FINAL$var))==0) FinalVar <- FinalVar + as.matrix(bdiag(FINAL$var/M[i],solve(MAX2$hessian)/M[i]))
      if(sum(is.nan(FINAL$var))>0 | sum(is.na(FINAL$var))>0) {NEW <- XimpsPar(T,C,X,Initpar,Initbn,cc,mt,m1,mt1,m2,mt2,m3,mt3,m12,mt12,m13,mt13,m23,mt23,m123,mt123,dists,2,d,twocens)
      Timped[,q] <- NEW[[1]][,1]
      Cimped[,q] <- NEW[[2]][,1]
      Ximped[,(q*7-6):(q*7)] <- NEW[[3]][,1:7]
      q <- q-1}
      if(Parmis[1,q]<0.5) {NEW <- XimpsPar(T,C,X,Initpar,Initbn,cc,mt,m1,mt1,m2,mt2,m3,mt3,m12,mt12,m13,mt13,m23,mt23,m123,mt123,dists,2,d,twocens)
      Timped[,q] <- NEW[[1]][,1]
      Cimped[,q] <- NEW[[2]][,1]
      Ximped[,(q*7-6):(q*7)] <- NEW[[3]][,1:7]
      q <- q-1}
      q <- q + 1
    }
    
    EstUpdates[i,] <- rowMeans(Parmis)
    SEUpdates[i,] <- sqrt(diag(ApproxVar(Parmis,BNmis,InitVar,FinalVar,M[i])[1:8,1:8]))
    
  }
  return(list(EstUpdates,SEUpdates,BNmis))
  
}

#Log Likelihood function for censored bivariate normal
CensNorm <- function(T,C,X,cnone,cone) {
  like<-function(t) {
    fc<-rep(0,cnone+cone)
    if(cnone>0) fc[1:cnone] <- -log(1/(t[8])*dnorm((T[1:cnone]-as.real(X[1:cnone,]%*%c(t[1],t[2],t[3],t[4],t[5],t[6],t[7])))/t[8]))
    if(cone>0) fc[(cnone+1):(cnone+cone)] <- -log(pnorm((C[(cnone+1):(cnone+cone)]-as.real(X[(cnone+1):(cnone+cone),]%*%c(t[1],t[2],t[3],t[4],t[5],t[6],t[7])))/t[8],lower.tail=FALSE))
    like <- sum(fc)
  }
  like
}

#Log Likelihood function for censored SNP-exponential
CensSNPexp <- function(T,C,X,k,cnone,cone) {
  like<-function(t) {
    fc<-rep(0,(cnone+cone))
    a <- acoef(t[9:length(t)],k,EXP=TRUE)
    P_k <- Pk(T[1:cnone],X[1:cnone,],t,k,a,cnone,EXP=TRUE)
    a <- c(a, rep(0, 3-length(a))) 
    fc[1:cnone] <- t[8]-(T[1:cnone]-as.vector(X[1:cnone,]%*%c(t[1],t[2],t[3],t[4],t[5],t[6],t[7])))/exp(t[8])-log((P_k[1:cnone])^2)+exp((T[1:cnone]-as.vector(X[1:cnone,]%*%c(t[1],t[2],t[3],t[4],t[5],t[6],t[7])))/exp(t[8]))
    if(cone>0)  fc[(cnone+1):(cnone+cone)] <- -log(pSNPexp(exp((C[(cnone+1):(cnone+cone)]-as.vector(X[(cnone+1):(cnone+cone),]%*%c(t[1],t[2],t[3],t[4],t[5],t[6],t[7])))/exp(t[8])),a))			    
    like <- sum(fc)
    if(is.nan(like)==TRUE | is.na(like)==TRUE ) like <- 10000 else if(abs(like)>1e100) like <- 10000 
    return(like)
  }
  like
}

#Log Likelihood function for censored SNP-normal
CensSNPnorm <- function(T,C,X,k,cnone,cone) {
  like<-function(t) {
    fc<-rep(0,(cnone+cone))
    a <- acoef(t[9:length(t)],k,EXP=FALSE)
    P_k <- Pk(T[1:cnone],X[1:cnone,],t,k,a,cnone,EXP=FALSE)
    a <- c(a, rep(0, 3-length(a))) 
    if(cnone>0) fc[1:cnone] <-  -log(1/t[8]*(P_k[1:cnone])^2*dnorm(((T[1:cnone]-as.vector(X[1:cnone,]%*%c(t[1],t[2],t[3],t[4],t[5],t[6],t[7])))/t[8])))
    if(cone>0)  fc[(cnone+1):(cnone+cone)] <- -log(pSNPnorm(((C[(cnone+1):(cnone+cone)]-as.vector(X[(cnone+1):(cnone+cone),]%*%c(t[1],t[2],t[3],t[4],t[5],t[6],t[7])))/t[8]),a))			    
    like <- sum(fc)
    if(is.nan(like)==TRUE | is.na(like)==TRUE ) like <- 10000 else if(abs(like)>1e100) like <- 10000 
    return(like)
  }
  like
}



#Function to find imputations for censored values using SNP
XimpsSNP <- function(T,C,X,k,InitSNP,ests,cc,mt,m1,mt1,m2,mt2,m3,mt3,m12,mt12,m13,mt13,m23,mt23,m123,mt123,EXP=TRUE,M,d,twocens) {
  Phi <- InitSNP[9:length(InitSNP)]
  a <- acoef(Phi,k,EXP=EXP)
  sig123 <- matrix(c(ests[13],ests[14],ests[15],ests[14],ests[16],ests[17],ests[15],ests[17],ests[18]),3,3)
  sig12 <- matrix(c(ests[13],ests[14],ests[14],ests[16]),2,2)
  sig13 <- matrix(c(ests[13],ests[15],ests[15],ests[18]),2,2)
  sig23 <- matrix(c(ests[16],ests[17],ests[17],ests[18]),2,2)
  
  X1imp <- matrix(0,M,(m1+mt1))
  if((m1+mt1)==0) X1imp<- NULL
  if((m1+mt1)>0){
    for(i in (cc+mt+1):(cc+mt+m1+mt1)){
      Xi <- rep(100+20*M)
      Xi[1] <- rtnorm(1,ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],sqrt(ests[13])/2,upper=d[1])
      rden <- GenTsnp(InitSNP,k,a,EXP,T[i],C[i],c(1,Xi[1],X[i,3:7]))*dnorm(Xi[1],mean=(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7]+as.real(matrix(c(ests[14],ests[15]),1,2)%*%solve(sig23)%*%matrix((X[i,3:4]-c(ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7])),2,1))),sd=sqrt((ests[13]-as.real(matrix(c(ests[14],ests[15]),1,2)%*%solve(sig23)%*%matrix(c(ests[14],ests[15]),2,1)))))
      for(t in 2:(100+20*M)) {
        Xstar <- rnorm(1,mean=Xi[(t-1)],sqrt(ests[13])/2)
        rnum <- GenTsnp(InitSNP,k,a,EXP,T[i],C[i],c(1,Xstar,X[i,3:7]))*dnorm(Xstar,mean=(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7]+as.real(matrix(c(ests[14],ests[15]),1,2)%*%solve(sig23)%*%matrix((X[i,3:4]-c(ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7])),2,1))),sd=sqrt((ests[13]-as.real(matrix(c(ests[14],ests[15]),1,2)%*%solve(sig23)%*%matrix(c(ests[14],ests[15]),2,1)))))
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
      Xi[1] <- rtnorm(1,ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7],sqrt(ests[16])/2,upper=d[2])
      rden <- GenTsnp(InitSNP,k,a,EXP,T[i],C[i],c(X[i,1:2],Xi[1],X[i,4:7]))*dnorm(Xi[1],mean=(ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7]+as.real(matrix(c(ests[14],ests[17]),1,2)%*%solve(sig13)%*%matrix((X[i,c(2,4)]-c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7])),2,1))),sd=sqrt((ests[16]-as.real(matrix(c(ests[14],ests[17]),1,2)%*%solve(sig13)%*%matrix(c(ests[14],ests[17]),2,1)))))
      for(t in 2:(100+20*M)) {
        Xstar <- rnorm(1,mean=Xi[(t-1)],sqrt(ests[16])/2)
        rnum <- GenTsnp(InitSNP,k,a,EXP,T[i],C[i],c(X[i,1:2],Xstar,X[i,4:7]))*dnorm(Xstar,mean=(ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7]+as.real(matrix(c(ests[14],ests[17]),1,2)%*%solve(sig13)%*%matrix((X[i,c(2,4)]-c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7])),2,1))),sd=sqrt((ests[16]-as.real(matrix(c(ests[14],ests[17]),1,2)%*%solve(sig13)%*%matrix(c(ests[14],ests[17]),2,1)))))
        r <- rnum/max(rden,1e-300)
        u <- runif(1)
        if(r>u & Xstar<d[2]) {rden <- rnum
        Xi[t] <- Xstar}   else Xi[t] <- Xi[(t-1)]	 
      }
      X2imp[,(i-(cc+mt+m1+mt1))] <- sample(Xi[seq(110,(100+20*M),20)])
    }}
  
  X3imp <- matrix(0,M,(m3+mt3))
  if((m3+mt3)==0) X3imp <- NULL
  if((m3+mt3)>0){
    for(i in (cc+mt+m1+mt1+m2+mt2+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3)){
      Xi <- rep(100+20*M)
      Xi[1] <- rtnorm(1,ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7],sqrt(ests[18])/2,upper=d[3])
      rden <- GenTsnp(InitSNP,k,a,EXP,T[i],C[i],c(X[i,1:3],Xi[1],X[i,5:7]))*dnorm(Xi[1],mean=(ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7]+as.real(matrix(c(ests[15],ests[17]),1,2)%*%solve(sig12)%*%matrix((X[i,2:3]-c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7])),2,1))),sd=sqrt((ests[18]-as.real(matrix(c(ests[15],ests[17]),1,2)%*%solve(sig12)%*%matrix(c(ests[15],ests[17]),2,1)))))
      for(t in 2:(100+20*M)) {
        Xstar <- rnorm(1,mean=Xi[(t-1)],sqrt(ests[18])/2)
        rnum <- GenTsnp(InitSNP,k,a,EXP,T[i],C[i],c(X[i,1:3],Xstar,X[i,5:7]))*dnorm(Xstar,mean=(ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7]+as.real(matrix(c(ests[15],ests[17]),1,2)%*%solve(sig12)%*%matrix((X[i,2:3]-c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7])),2,1))),sd=sqrt((ests[18]-as.real(matrix(c(ests[15],ests[17]),1,2)%*%solve(sig12)%*%matrix(c(ests[15],ests[17]),2,1)))))
        r <- rnum/max(rden,1e-300)
        u <- runif(1)
        if(r>u & Xstar<d[3]) {rden <- rnum
        Xi[t] <- Xstar}   else Xi[t] <- Xi[(t-1)]	 
      }
      X3imp[,(i-(cc+mt+m1+mt1+m2+mt2))] <- sample(Xi[seq(110,(100+20*M),20)])
    }}
  
  Genvals <- rmvnorm((50*M+200)*(m12+mt12),c(0,0),make.positive.definite(sig12/4,tol=0.001))
  l <- 1
  X12imp <- matrix(0,M,(m12+mt12)*2)
  if((m12+mt12)==0) X12imp<- NULL
  if((m12+mt12)>0){
    for(i in (cc+mt+m1+mt1+m2+mt2+1):(cc+mt+m1+mt1+m2+mt2+m12+mt12)){
      if(any(twocens==i)) Upper=c(d[1],log(2)) else Upper=c(d[1],d[2])
      Xi <- matrix(0,(50*M+200),2)
      Xi[1,] <- rtmvnorm(1,mean=c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7]),sig12/2, upper=Upper,algorithm="gibbs")
      rden <- GenTsnp(InitSNP,k,a,EXP,T[i],C[i],c(X[i,1],Xi[1,],X[i,4:7]))*dmvnorm(Xi[1,],mean=(c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7])+as.vector(matrix(c(ests[15],ests[17]),2,1)*1/ests[18]*(X[i,4]-ests[9]-ests[10]*X[i,5]-ests[11]*X[i,6]-ests[12]*X[i,7]))),sigma=sig12-1/ests[18]*matrix(c(ests[15],ests[17]),2,1)%*%t(matrix(c(ests[15],ests[17]),2,1)))
      for(t in 2:(50*M+200)) {
        Xstar <- Xi[(t-1),]+Genvals[l,]
        l <- l+1
        rnum <- GenTsnp(InitSNP,k,a,EXP,T[i],C[i],c(X[i,1],Xstar,X[i,4:7]))*dmvnorm(Xstar,mean=(c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7])+as.vector(matrix(c(ests[15],ests[17]),2,1)*1/ests[18]*(X[i,4]-ests[9]-ests[10]*X[i,5]-ests[11]*X[i,6]-ests[12]*X[i,7]))),sigma=sig12-1/ests[18]*matrix(c(ests[15],ests[17]),2,1)%*%t(matrix(c(ests[15],ests[17]),2,1)))
        r <- rnum/max(rden,1e-300)
        u <- runif(1)
        if(r>u & Xstar[1]<Upper[1] & Xstar[2]<Upper[2]) {rden <- rnum
        Xi[t,] <- Xstar}  else Xi[t,] <- Xi[(t-1),]		 
      } 
      if(M>1) X12imp[,((i-(cc+mt+m1+mt1+m2+mt2))*2-1):((i-(cc+mt+m1+mt1+m2+mt2))*2)] <- Xi[sample(seq(250,(200+50*M),50)),]
      if(M==1) X12imp[,((i-(cc+mt+m1+mt1+m2+mt2))*2-1):((i-(cc+mt+m1+mt1+m2+mt2))*2)] <- Xi[seq(250,(200+50*M),50),]
    }}
  
  Genvals <- rmvnorm((50*M+200)*(m13+mt13),c(0,0),make.positive.definite(sig13/4,tol=0.001))
  l <- 1
  X13imp <- matrix(0,M,(m13+mt13)*2)
  if((m13+mt13)==0) X13imp<- NULL
  if((m13+mt13)>0){
    for(i in (cc+mt+m1+mt1+m2+mt2+m12+mt12+1):(cc+mt+m1+mt1+m2+mt2+m12+mt12+m13+mt13)){
      Xi <- matrix(0,(50*M+200),2)
      Xi[1,] <- rtmvnorm(1,mean=c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7]),sig13/2, upper=c(d[1],d[3]),algorithm="gibbs")
      rden <- GenTsnp(InitSNP,k,a,EXP,T[i],C[i],c(X[i,1],Xi[1,1],X[i,3],Xi[1,2],X[i,5:7]))*dmvnorm(Xi[1,],mean=(c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7])+as.vector(matrix(c(ests[14],ests[17]),2,1)*1/ests[16]*(X[i,3]-ests[5]-ests[6]*X[i,5]-ests[7]*X[i,6]-ests[8]*X[i,7]))),sigma=sig13-1/ests[16]*matrix(c(ests[14],ests[17]),2,1)%*%t(matrix(c(ests[14],ests[17]),2,1)))
      for(t in 2:(50*M+200)) {
        Xstar <- Xi[(t-1),]+Genvals[l,]
        l <- l+1
        rnum <- GenTsnp(InitSNP,k,a,EXP,T[i],C[i],c(X[i,1],Xstar[1],X[i,3],Xstar[2],X[i,5:7]))*dmvnorm(Xstar,mean=(c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7])+as.vector(matrix(c(ests[14],ests[17]),2,1)*1/ests[16]*(X[i,3]-ests[5]-ests[6]*X[i,5]-ests[7]*X[i,6]-ests[8]*X[i,7]))),sigma=sig13-1/ests[16]*matrix(c(ests[14],ests[17]),2,1)%*%t(matrix(c(ests[14],ests[17]),2,1)))
        r <- rnum/max(rden,1e-300)
        u <- runif(1)
        if(r>u & Xstar[1]<d[1] & Xstar[2]<d[3]) {rden <- rnum
        Xi[t,] <- Xstar}  else Xi[t,] <- Xi[(t-1),]		 
      } 
      if(M>1) X13imp[,((i-(cc+mt+m1+mt1+m2+mt2+m12+mt12))*2-1):((i-(cc+mt+m1+mt1+m2+mt2+m12+mt12))*2)] <- Xi[sample(seq(250,(200+50*M),50)),]
      if(M==1) X13imp[,((i-(cc+mt+m1+mt1+m2+mt+m12+mt12))*2-1):((i-(cc+mt+m1+mt1+m2+mt2+m12+mt12))*2)] <- Xi[seq(250,(200+50*M),50),]
    }}
  
  Genvals <- rmvnorm((50*M+200)*(m23+mt23),c(0,0),make.positive.definite(sig23/4,tol=0.001))
  l <- 1
  X23imp <- matrix(0,M,(m23+mt23)*2)
  if((m23+mt23)==0) X23imp<- NULL
  if((m23+mt23)>0){
    for(i in (cc+mt+m1+mt1+m2+mt2+m12+mt12+m13+mt13+1):(cc+mt+m1+mt1+m2+mt2+m12+mt12+m13+mt13+m23+mt23)){
      if(any(twocens==i)) Upper=c(log(2),d[3]) else Upper=c(d[2],d[3])
      Xi <- matrix(0,(50*M+200),2)
      Xi[1,] <- rtmvnorm(1,mean=c(ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7]),sig23/2, upper=Upper,algorithm="gibbs")
      rden <- GenTsnp(InitSNP,k,a,EXP,T[i],C[i],c(X[i,1:2],Xi[1,],X[i,5:7]))*dmvnorm(Xi[1,],mean=(c(ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7])+as.vector(matrix(c(ests[14],ests[15]),2,1)*1/ests[13]*(X[i,2]-ests[1]-ests[2]*X[i,5]-ests[3]*X[i,6]-ests[4]*X[i,7]))),sigma=(sig23-1/ests[13]*matrix(c(ests[14],ests[15]),2,1)%*%t(matrix(c(ests[14],ests[15]),2,1))))
      for(t in 2:(50*M+200)) {
        Xstar <- Xi[(t-1),]+Genvals[l,]
        l <- l+1
        rnum <- GenTsnp(InitSNP,k,a,EXP,T[i],C[i],c(X[i,1:2],Xstar,X[i,5:7]))*dmvnorm(Xstar,mean=(c(ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7])+as.vector(matrix(c(ests[14],ests[15]),2,1)*1/ests[13]*(X[i,2]-ests[1]-ests[2]*X[i,5]-ests[3]*X[i,6]-ests[4]*X[i,7]))),sigma=(sig23-1/ests[13]*matrix(c(ests[14],ests[15]),2,1)%*%t(matrix(c(ests[14],ests[15]),2,1))))
        r <- rnum/max(rden,1e-300)
        u <- runif(1)
        if(r>u & Xstar[1]<Upper[1] & Xstar[2]<Upper[2]) {rden <- rnum
        Xi[t,] <- Xstar}  else Xi[t,] <- Xi[(t-1),]		 
      } 
      if(M>1) X23imp[,((i-(cc+mt+m1+mt1+m2+mt2+m12+mt12+m13+mt13))*2-1):((i-(cc+mt+m1+mt1+m2+mt2+m12+mt12+m13+mt13))*2)] <- Xi[sample(seq(250,(200+50*M),50)),]
      if(M==1) X23imp[,((i-(cc+mt+m1+mt1+m2+mt+m12+mt12+m13+mt13))*2-1):((i-(cc+mt+m1+mt1+m2+mt2+m12+mt12+m13+mt13))*2)] <- Xi[seq(250,(200+50*M),50),]
    }}
  
  
  Genvals <- rmvnorm((50*M+200)*(m123+mt123),c(0,0,0),make.positive.definite(sig123/4,tol=0.001))
  l <- 1
  X123imp <- matrix(0,M,(m123+mt123)*3)
  if((m123+mt123)==0) X123imp<- NULL
  if((m123+mt123)>0){
    for(i in (cc+mt+m1+mt1+m2+mt2+m12+mt12+m13+mt13+m23+mt23+1):(cc+mt+m1+mt1+m2+mt2+m12+mt12+m13+mt13+m23+mt23+m123+mt123)){
      if(any(twocens==i)) Upper=c(d[1],log(2),d[3]) else Upper=d
      Xi <- matrix(0,(50*M+200),3)
      Xi[1,] <- rtmvnorm(1,mean=c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7]),sig123/2, upper=Upper,algorithm="gibbs")
      rden <- GenTsnp(InitSNP,k,a,EXP,T[i],C[i],c(X[i,1],Xi[1,],X[i,5:7]))*dmvnorm(Xi[1,],mean=c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7]),sigma=matrix(c(ests[13],ests[14],ests[15],ests[14],ests[16],ests[17],ests[15],ests[17],ests[18]),3,3))
      for(t in 2:(50*M+200)) {
        Xstar <- Xi[(t-1),]+Genvals[l,]
        l <- l+1
        rnum <- GenTsnp(InitSNP,k,a,EXP,T[i],C[i],c(X[i,1],Xstar,X[i,5:7]))*dmvnorm(Xstar,mean=c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7]),sigma=matrix(c(ests[13],ests[14],ests[15],ests[14],ests[16],ests[17],ests[15],ests[17],ests[18]),3,3))
        r <- rnum/max(rden,1e-300)
        u <- runif(1)
        if(r>u & Xstar[1]<Upper[1] & Xstar[2]<Upper[2] & Xstar[3]<Upper[3]) {rden <- rnum
        Xi[t,] <- Xstar}  else Xi[t,] <- Xi[(t-1),]		 
      } 
      if(M>1) X123imp[,((i-(cc+mt+m1+mt1+m2+mt2+m12+mt12+m13+mt13+m23+mt23))*3-2):((i-(cc+mt+m1+mt1+m2+mt2+m12+mt12+m13+mt13+m23+mt23))*3)] <- Xi[sample(seq(250,(200+50*M),50)),]
    }}
  
  Timputed <- NULL
  Cimputed <- NULL
  Ximputed <- NULL
  for(i in 1:M){
    Timputed <- cbind(Timputed,T)
    Cimputed <- cbind(Cimputed,C)
    Ximputed <- cbind(Ximputed,rbind(X[1:(cc+mt),],cbind(1,X1imp[i,],X[(cc+mt+1):(cc+mt+m1+mt1),3:7]),cbind(X[(cc+mt+m1+mt1+1):(cc+mt+m1+mt1+m2+mt2),1:2],X2imp[i,],X[(cc+mt+m1+mt1+1):(cc+mt+m1+mt1+m2+mt2),4:7]),cbind(X[(cc+mt+m1+mt1+m2+mt2+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3),1:3],X3imp[i,],X[(cc+mt+m1+mt1+m2+mt2+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3),5:7]),cbind(1,X12imp[i,seq(1,(2*(m12+mt12)),2)],X12imp[i,seq(2,(2*(m12+mt12)),2)],X[(cc+mt+m1+mt1+m2+mt2+m3+mt3+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12),4:7]),cbind(1,X13imp[i,seq(1,(2*(m13+mt13)),2)],X[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13),3],X13imp[i,seq(2,(2*(m13+mt13)),2)],X[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13),5:7]),cbind(1,X[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23),2],X23imp[i,seq(1,(2*(m23+mt23)),2)],X23imp[i,seq(2,(2*(m23+mt23)),2)],X[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23),5:7]),cbind(1,X123imp[i,seq(1,(3*(m123+mt123)),3)],X123imp[i,seq(2,(3*(m123+mt123)),3)],X123imp[i,seq(3,(3*(m123+mt123)),3)],X[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23+m123+mt123),5:7])))
  } 
  
  Timputed <- rbind(Timputed[1:cc,],Timputed[(cc+mt+1):(cc+mt+m1),],Timputed[(cc+mt+m1+mt1+1):(cc+mt+m1+mt1+m2),],Timputed[(cc+mt+m1+mt1+m2+mt2+1):(cc+mt+m1+mt1+m2+mt2+m3),],Timputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12),],Timputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13),],Timputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23),],Timputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23+m123),],Timputed[(cc+1):(cc+mt),],Timputed[(cc+mt+m1+1):(cc+mt+m1+mt1),],Timputed[(cc+mt+m1+mt1+m2+1):(cc+mt+m1+mt1+m2+mt2),],Timputed[(cc+mt+m1+mt1+m2+mt2+m3+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3),],Timputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12),],Timputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13),],Timputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23),],Timputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23+m123+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23+m123+mt123),])
  Cimputed <- rbind(Cimputed[1:cc,],Cimputed[(cc+mt+1):(cc+mt+m1),],Cimputed[(cc+mt+m1+mt1+1):(cc+mt+m1+mt1+m2),],Cimputed[(cc+mt+m1+mt1+m2+mt2+1):(cc+mt+m1+mt1+m2+mt2+m3),],Cimputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12),],Cimputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13),],Cimputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23),],Cimputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23+m123),],Cimputed[(cc+1):(cc+mt),],Cimputed[(cc+mt+m1+1):(cc+mt+m1+mt1),],Cimputed[(cc+mt+m1+mt1+m2+1):(cc+mt+m1+mt1+m2+mt2),],Cimputed[(cc+mt+m1+mt1+m2+mt2+m3+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3),],Cimputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12),],Cimputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13),],Cimputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23),],Cimputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23+m123+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23+m123+mt123),])
  Ximputed <- rbind(Ximputed[1:cc,],Ximputed[(cc+mt+1):(cc+mt+m1),],Ximputed[(cc+mt+m1+mt1+1):(cc+mt+m1+mt1+m2),],Ximputed[(cc+mt+m1+mt1+m2+mt2+1):(cc+mt+m1+mt1+m2+mt2+m3),],Ximputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12),],Ximputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13),],Ximputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23),],Ximputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23+m123),],Ximputed[(cc+1):(cc+mt),],Ximputed[(cc+mt+m1+1):(cc+mt+m1+mt1),],Ximputed[(cc+mt+m1+mt1+m2+1):(cc+mt+m1+mt1+m2+mt2),],Ximputed[(cc+mt+m1+mt1+m2+mt2+m3+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3),],Ximputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12),],Ximputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13),],Ximputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23),],Ximputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23+m123+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23+m123+mt123),])
  
  list(Timputed,Cimputed,Ximputed)
}  


#Function to find imputations for censored values using parametric distribution
XimpsPar <- function(T,C,X,InitPar,ests,cc,mt,m1,mt1,m2,mt2,m3,mt3,m12,mt12,m13,mt13,m23,mt23,m123,mt123,dist,M,d,twocens) {
  
  sig123 <- matrix(c(ests[13],ests[14],ests[15],ests[14],ests[16],ests[17],ests[15],ests[17],ests[18]),3,3)
  sig12 <- matrix(c(ests[13],ests[14],ests[14],ests[16]),2,2)
  sig13 <- matrix(c(ests[13],ests[15],ests[15],ests[18]),2,2)
  sig23 <- matrix(c(ests[16],ests[17],ests[17],ests[18]),2,2)
  
  X1imp <- matrix(0,M,(m1+mt1))
  if((m1+mt1)==0) X1imp<- NULL
  if((m1+mt1)>0){
    for(i in (cc+mt+1):(cc+mt+m1+mt1)){
      Xi <- rep(100+20*M)
      Xi[1] <- rtnorm(1,ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],sqrt(ests[13])/2,upper=d[1])
      rden <- GenT(dist,InitPar,T[i],C[i],c(1,Xi[1],X[i,3:7]))*dnorm(Xi[1],mean=(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7]+as.real(matrix(c(ests[14],ests[15]),1,2)%*%solve(sig23)%*%matrix((X[i,3:4]-c(ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7])),2,1))),sd=sqrt((ests[13]-as.real(matrix(c(ests[14],ests[15]),1,2)%*%solve(sig23)%*%matrix(c(ests[14],ests[15]),2,1)))))
      for(t in 2:(100+20*M)) {
        Xstar <- rnorm(1,mean=Xi[(t-1)],sqrt(ests[13])/2)
        rnum <- GenT(dist,InitPar,T[i],C[i],c(1,Xstar,X[i,3:7]))*dnorm(Xstar,mean=(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7]+as.real(matrix(c(ests[14],ests[15]),1,2)%*%solve(sig23)%*%matrix((X[i,3:4]-c(ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7])),2,1))),sd=sqrt((ests[13]-as.real(matrix(c(ests[14],ests[15]),1,2)%*%solve(sig23)%*%matrix(c(ests[14],ests[15]),2,1)))))
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
      Xi[1] <- rtnorm(1,ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7],sqrt(ests[16])/2,upper=d[2])
      rden <- GenT(dist,InitPar,T[i],C[i],c(X[i,1:2],Xi[1],X[i,4:7]))*dnorm(Xi[1],mean=(ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7]+as.real(matrix(c(ests[14],ests[17]),1,2)%*%solve(sig13)%*%matrix((X[i,c(2,4)]-c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7])),2,1))),sd=sqrt((ests[16]-as.real(matrix(c(ests[14],ests[17]),1,2)%*%solve(sig13)%*%matrix(c(ests[14],ests[17]),2,1)))))
      for(t in 2:(100+20*M)) {
        Xstar <- rnorm(1,mean=Xi[(t-1)],sqrt(ests[16])/2)
        rnum <- GenT(dist,InitPar,T[i],C[i],c(X[i,1:2],Xstar,X[i,4:7]))*dnorm(Xstar,mean=(ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7]+as.real(matrix(c(ests[14],ests[17]),1,2)%*%solve(sig13)%*%matrix((X[i,c(2,4)]-c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7])),2,1))),sd=sqrt((ests[16]-as.real(matrix(c(ests[14],ests[17]),1,2)%*%solve(sig13)%*%matrix(c(ests[14],ests[17]),2,1)))))
        r <- rnum/max(rden,1e-300)
        u <- runif(1)
        if(r>u & Xstar<d[2]) {rden <- rnum
        Xi[t] <- Xstar}   else Xi[t] <- Xi[(t-1)]	 
      }
      X2imp[,(i-(cc+mt+m1+mt1))] <- sample(Xi[seq(110,(100+20*M),20)])
    }}
  
  X3imp <- matrix(0,M,(m3+mt3))
  if((m3+mt3)==0) X3imp <- NULL
  if((m3+mt3)>0){
    for(i in (cc+mt+m1+mt1+m2+mt2+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3)){
      Xi <- rep(100+20*M)
      Xi[1] <- rtnorm(1,ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7],sqrt(ests[18])/2,upper=d[3])
      rden <- GenT(dist,InitPar,T[i],C[i],c(X[i,1:3],Xi[1],X[i,5:7]))*dnorm(Xi[1],mean=(ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7]+as.real(matrix(c(ests[15],ests[17]),1,2)%*%solve(sig12)%*%matrix((X[i,2:3]-c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7])),2,1))),sd=sqrt((ests[18]-as.real(matrix(c(ests[15],ests[17]),1,2)%*%solve(sig12)%*%matrix(c(ests[15],ests[17]),2,1)))))
      for(t in 2:(100+20*M)) {
        Xstar <- rnorm(1,mean=Xi[(t-1)],sqrt(ests[18])/2)
        rnum <- GenT(dist,InitPar,T[i],C[i],c(X[i,1:3],Xstar,X[i,5:7]))*dnorm(Xstar,mean=(ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7]+as.real(matrix(c(ests[15],ests[17]),1,2)%*%solve(sig12)%*%matrix((X[i,2:3]-c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7])),2,1))),sd=sqrt((ests[18]-as.real(matrix(c(ests[15],ests[17]),1,2)%*%solve(sig12)%*%matrix(c(ests[15],ests[17]),2,1)))))
        r <- rnum/max(rden,1e-300)
        u <- runif(1)
        if(r>u & Xstar<d[3]) {rden <- rnum
        Xi[t] <- Xstar}   else Xi[t] <- Xi[(t-1)]	 
      }
      X3imp[,(i-(cc+mt+m1+mt1+m2+mt2))] <- sample(Xi[seq(110,(100+20*M),20)])
    }}
  
  Genvals <- rmvnorm((50*M+200)*(m12+mt12),c(0,0),make.positive.definite(sig12/4,tol=0.001))
  l <- 1
  X12imp <- matrix(0,M,(m12+mt12)*2)
  if((m12+mt12)==0) X12imp<- NULL
  if((m12+mt12)>0){
    for(i in (cc+mt+m1+mt1+m2+mt2+1):(cc+mt+m1+mt1+m2+mt2+m12+mt12)){
      if(any(twocens==i)) Upper=c(d[1],log(2)) else Upper=c(d[1],d[2])
      Xi <- matrix(0,(50*M+200),2)
      Xi[1,] <- rtmvnorm(1,mean=c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7]),sig12/2, upper=Upper,algorithm="gibbs")
      rden <- GenT(dist,InitPar,T[i],C[i],c(X[i,1],Xi[1,],X[i,4:7]))*dmvnorm(Xi[1,],mean=(c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7])+as.vector(matrix(c(ests[15],ests[17]),2,1)*1/ests[18]*(X[i,4]-ests[9]-ests[10]*X[i,5]-ests[11]*X[i,6]-ests[12]*X[i,7]))),sigma=sig12-1/ests[18]*matrix(c(ests[15],ests[17]),2,1)%*%t(matrix(c(ests[15],ests[17]),2,1)))
      for(t in 2:(50*M+200)) {
        Xstar <- Xi[(t-1),]+Genvals[l,]
        l <- l+1
        rnum <- GenT(dist,InitPar,T[i],C[i],c(X[i,1],Xstar,X[i,4:7]))*dmvnorm(Xstar,mean=(c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7])+as.vector(matrix(c(ests[15],ests[17]),2,1)*1/ests[18]*(X[i,4]-ests[9]-ests[10]*X[i,5]-ests[11]*X[i,6]-ests[12]*X[i,7]))),sigma=sig12-1/ests[18]*matrix(c(ests[15],ests[17]),2,1)%*%t(matrix(c(ests[15],ests[17]),2,1)))
        r <- rnum/max(rden,1e-300)
        u <- runif(1)
        if(r>u & Xstar[1]<Upper[1] & Xstar[2]<Upper[2]) {rden <- rnum
        Xi[t,] <- Xstar}  else Xi[t,] <- Xi[(t-1),]		 
      } 
      if(M>1) X12imp[,((i-(cc+mt+m1+mt1+m2+mt2))*2-1):((i-(cc+mt+m1+mt1+m2+mt2))*2)] <- Xi[sample(seq(250,(200+50*M),50)),]
      if(M==1) X12imp[,((i-(cc+mt+m1+mt1+m2+mt2))*2-1):((i-(cc+mt+m1+mt1+m2+mt2))*2)] <- Xi[seq(250,(200+50*M),50),]
    }}
  
  Genvals <- rmvnorm((50*M+200)*(m13+mt13),c(0,0),make.positive.definite(sig13/4,tol=0.001))
  l <- 1
  X13imp <- matrix(0,M,(m13+mt13)*2)
  if((m13+mt13)==0) X13imp<- NULL
  if((m13+mt13)>0){
    for(i in (cc+mt+m1+mt1+m2+mt2+m12+mt12+1):(cc+mt+m1+mt1+m2+mt2+m12+mt12+m13+mt13)){
      Xi <- matrix(0,(50*M+200),2)
      Xi[1,] <- rtmvnorm(1,mean=c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7]),sig13/2, upper=c(d[1],d[3]),algorithm="gibbs")
      rden <- GenT(dist,InitPar,T[i],C[i],c(X[i,1],Xi[1,1],X[i,3],Xi[1,2],X[i,5:7]))*dmvnorm(Xi[1,],mean=(c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7])+as.vector(matrix(c(ests[14],ests[17]),2,1)*1/ests[16]*(X[i,3]-ests[5]-ests[6]*X[i,5]-ests[7]*X[i,6]-ests[8]*X[i,7]))),sigma=sig13-1/ests[16]*matrix(c(ests[14],ests[17]),2,1)%*%t(matrix(c(ests[14],ests[17]),2,1)))
      for(t in 2:(50*M+200)) {
        Xstar <- Xi[(t-1),]+Genvals[l,]
        l <- l+1
        rnum <- GenT(dist,InitPar,T[i],C[i],c(X[i,1],Xstar[1],X[i,3],Xstar[2],X[i,5:7]))*dmvnorm(Xstar,mean=(c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7])+as.vector(matrix(c(ests[14],ests[17]),2,1)*1/ests[16]*(X[i,3]-ests[5]-ests[6]*X[i,5]-ests[7]*X[i,6]-ests[8]*X[i,7]))),sigma=sig13-1/ests[16]*matrix(c(ests[14],ests[17]),2,1)%*%t(matrix(c(ests[14],ests[17]),2,1)))
        r <- rnum/max(rden,1e-300)
        u <- runif(1)
        if(r>u & Xstar[1]<d[1] & Xstar[2]<d[3]) {rden <- rnum
        Xi[t,] <- Xstar}  else Xi[t,] <- Xi[(t-1),]		 
      } 
      if(M>1) X13imp[,((i-(cc+mt+m1+mt1+m2+mt2+m12+mt12))*2-1):((i-(cc+mt+m1+mt1+m2+mt2+m12+mt12))*2)] <- Xi[sample(seq(250,(200+50*M),50)),]
      if(M==1) X13imp[,((i-(cc+mt+m1+mt1+m2+mt+m12+mt12))*2-1):((i-(cc+mt+m1+mt1+m2+mt2+m12+mt12))*2)] <- Xi[seq(250,(200+50*M),50),]
    }}
  
  Genvals <- rmvnorm((50*M+200)*(m23+mt23),c(0,0),make.positive.definite(sig23/4,tol=0.001))
  l <- 1
  X23imp <- matrix(0,M,(m23+mt23)*2)
  if((m23+mt23)==0) X23imp<- NULL
  if((m23+mt23)>0){
    for(i in (cc+mt+m1+mt1+m2+mt2+m12+mt12+m13+mt13+1):(cc+mt+m1+mt1+m2+mt2+m12+mt12+m13+mt13+m23+mt23)){
      if(any(twocens==i)) Upper=c(log(2),d[3]) else Upper=c(d[2],d[3])
      Xi <- matrix(0,(50*M+200),2)
      Xi[1,] <- rtmvnorm(1,mean=c(ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7]),sig23/2, upper=Upper,algorithm="gibbs")
      rden <- GenT(dist,InitPar,T[i],C[i],c(X[i,1:2],Xi[1,],X[i,5:7]))*dmvnorm(Xi[1,],mean=(c(ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7])+as.vector(matrix(c(ests[14],ests[15]),2,1)*1/ests[13]*(X[i,2]-ests[1]-ests[2]*X[i,5]-ests[3]*X[i,6]-ests[4]*X[i,7]))),sigma=(sig23-1/ests[13]*matrix(c(ests[14],ests[15]),2,1)%*%t(matrix(c(ests[14],ests[15]),2,1))))
      for(t in 2:(50*M+200)) {
        Xstar <- Xi[(t-1),]+Genvals[l,]
        l <- l+1
        rnum <- GenT(dist,InitPar,T[i],C[i],c(X[i,1:2],Xstar,X[i,5:7]))*dmvnorm(Xstar,mean=(c(ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7])+as.vector(matrix(c(ests[14],ests[15]),2,1)*1/ests[13]*(X[i,2]-ests[1]-ests[2]*X[i,5]-ests[3]*X[i,6]-ests[4]*X[i,7]))),sigma=(sig23-1/ests[13]*matrix(c(ests[14],ests[15]),2,1)%*%t(matrix(c(ests[14],ests[15]),2,1))))
        r <- rnum/max(rden,1e-300)
        u <- runif(1)
        if(r>u & Xstar[1]<Upper[1] & Xstar[2]<Upper[2]) {rden <- rnum
        Xi[t,] <- Xstar}  else Xi[t,] <- Xi[(t-1),]		 
      } 
      if(M>1) X23imp[,((i-(cc+mt+m1+mt1+m2+mt2+m12+mt12+m13+mt13))*2-1):((i-(cc+mt+m1+mt1+m2+mt2+m12+mt12+m13+mt13))*2)] <- Xi[sample(seq(250,(200+50*M),50)),]
      if(M==1) X23imp[,((i-(cc+mt+m1+mt1+m2+mt+m12+mt12+m13+mt13))*2-1):((i-(cc+mt+m1+mt1+m2+mt2+m12+mt12+m13+mt13))*2)] <- Xi[seq(250,(200+50*M),50),]
    }}
  
  
  Genvals <- rmvnorm((50*M+200)*(m123+mt123),c(0,0,0),make.positive.definite(sig123/4,tol=0.001))
  l <- 1
  X123imp <- matrix(0,M,(m123+mt123)*3)
  if((m123+mt123)==0) X123imp<- NULL
  if((m123+mt123)>0){
    for(i in (cc+mt+m1+mt1+m2+mt2+m12+mt12+m13+mt13+m23+mt23+1):(cc+mt+m1+mt1+m2+mt2+m12+mt12+m13+mt13+m23+mt23+m123+mt123)){
      if(any(twocens==i)) Upper=c(d[1],log(2),d[3]) else Upper=d
      Xi <- matrix(0,(50*M+200),3)
      Xi[1,] <- rtmvnorm(1,mean=c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7]),sig123/2, upper=Upper,algorithm="gibbs")
      rden <- GenT(dist,InitPar,T[i],C[i],c(X[i,1],Xi[1,],X[i,5:7]))*dmvnorm(Xi[1,],mean=c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7]),sigma=matrix(c(ests[13],ests[14],ests[15],ests[14],ests[16],ests[17],ests[15],ests[17],ests[18]),3,3))
      for(t in 2:(50*M+200)) {
        Xstar <- Xi[(t-1),]+Genvals[l,]
        l <- l+1
        rnum <- GenT(dist,InitPar,T[i],C[i],c(X[i,1],Xstar,X[i,5:7]))*dmvnorm(Xstar,mean=c(ests[1]+ests[2]*X[i,5]+ests[3]*X[i,6]+ests[4]*X[i,7],ests[5]+ests[6]*X[i,5]+ests[7]*X[i,6]+ests[8]*X[i,7],ests[9]+ests[10]*X[i,5]+ests[11]*X[i,6]+ests[12]*X[i,7]),sigma=matrix(c(ests[13],ests[14],ests[15],ests[14],ests[16],ests[17],ests[15],ests[17],ests[18]),3,3))
        r <- rnum/max(rden,1e-300)
        u <- runif(1)
        if(r>u & Xstar[1]<Upper[1] & Xstar[2]<Upper[2] & Xstar[3]<Upper[3]) {rden <- rnum
        Xi[t,] <- Xstar}  else Xi[t,] <- Xi[(t-1),]		 
      } 
      if(M>1) X123imp[,((i-(cc+mt+m1+mt1+m2+mt2+m12+mt12+m13+mt13+m23+mt23))*3-2):((i-(cc+mt+m1+mt1+m2+mt2+m12+mt12+m13+mt13+m23+mt23))*3)] <- Xi[sample(seq(250,(200+50*M),50)),]
      if(M==1) X123imp[,((i-(cc+mt+m1+mt1+m2+mt+m12+mt12+m13+mt13+m23+mt23))*3-2):((i-(cc+mt+m1+mt1+m2+mt2+m12+mt12+m13+mt13+m23+mt23))*3)] <- Xi[seq(250,(200+50*M),50),]
    }}
  
  Timputed <- NULL
  Cimputed <- NULL
  Ximputed <- NULL
  for(i in 1:M){
    Timputed <- cbind(Timputed,T)
    Cimputed <- cbind(Cimputed,C)
    Ximputed <- cbind(Ximputed,rbind(X[1:(cc+mt),],cbind(1,X1imp[i,],X[(cc+mt+1):(cc+mt+m1+mt1),3:7]),cbind(X[(cc+mt+m1+mt1+1):(cc+mt+m1+mt1+m2+mt2),1:2],X2imp[i,],X[(cc+mt+m1+mt1+1):(cc+mt+m1+mt1+m2+mt2),4:7]),cbind(X[(cc+mt+m1+mt1+m2+mt2+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3),1:3],X3imp[i,],X[(cc+mt+m1+mt1+m2+mt2+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3),5:7]),cbind(1,X12imp[i,seq(1,(2*(m12+mt12)),2)],X12imp[i,seq(2,(2*(m12+mt12)),2)],X[(cc+mt+m1+mt1+m2+mt2+m3+mt3+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12),4:7]),cbind(1,X13imp[i,seq(1,(2*(m13+mt13)),2)],X[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13),3],X13imp[i,seq(2,(2*(m13+mt13)),2)],X[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13),5:7]),cbind(1,X[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23),2],X23imp[i,seq(1,(2*(m23+mt23)),2)],X23imp[i,seq(2,(2*(m23+mt23)),2)],X[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23),5:7]),cbind(1,X123imp[i,seq(1,(3*(m123+mt123)),3)],X123imp[i,seq(2,(3*(m123+mt123)),3)],X123imp[i,seq(3,(3*(m123+mt123)),3)],X[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23+m123+mt123),5:7])))
  } 
  
  Timputed <- rbind(Timputed[1:cc,],Timputed[(cc+mt+1):(cc+mt+m1),],Timputed[(cc+mt+m1+mt1+1):(cc+mt+m1+mt1+m2),],Timputed[(cc+mt+m1+mt1+m2+mt2+1):(cc+mt+m1+mt1+m2+mt2+m3),],Timputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12),],Timputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13),],Timputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23),],Timputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23+m123),],Timputed[(cc+1):(cc+mt),],Timputed[(cc+mt+m1+1):(cc+mt+m1+mt1),],Timputed[(cc+mt+m1+mt1+m2+1):(cc+mt+m1+mt1+m2+mt2),],Timputed[(cc+mt+m1+mt1+m2+mt2+m3+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3),],Timputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12),],Timputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13),],Timputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23),],Timputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23+m123+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23+m123+mt123),])
  Cimputed <- rbind(Cimputed[1:cc,],Cimputed[(cc+mt+1):(cc+mt+m1),],Cimputed[(cc+mt+m1+mt1+1):(cc+mt+m1+mt1+m2),],Cimputed[(cc+mt+m1+mt1+m2+mt2+1):(cc+mt+m1+mt1+m2+mt2+m3),],Cimputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12),],Cimputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13),],Cimputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23),],Cimputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23+m123),],Cimputed[(cc+1):(cc+mt),],Cimputed[(cc+mt+m1+1):(cc+mt+m1+mt1),],Cimputed[(cc+mt+m1+mt1+m2+1):(cc+mt+m1+mt1+m2+mt2),],Cimputed[(cc+mt+m1+mt1+m2+mt2+m3+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3),],Cimputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12),],Cimputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13),],Cimputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23),],Cimputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23+m123+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23+m123+mt123),])
  Ximputed <- rbind(Ximputed[1:cc,],Ximputed[(cc+mt+1):(cc+mt+m1),],Ximputed[(cc+mt+m1+mt1+1):(cc+mt+m1+mt1+m2),],Ximputed[(cc+mt+m1+mt1+m2+mt2+1):(cc+mt+m1+mt1+m2+mt2+m3),],Ximputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12),],Ximputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13),],Ximputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23),],Ximputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23+m123),],Ximputed[(cc+1):(cc+mt),],Ximputed[(cc+mt+m1+1):(cc+mt+m1+mt1),],Ximputed[(cc+mt+m1+mt1+m2+1):(cc+mt+m1+mt1+m2+mt2),],Ximputed[(cc+mt+m1+mt1+m2+mt2+m3+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3),],Ximputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12),],Ximputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13),],Ximputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23),],Ximputed[(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23+m123+1):(cc+mt+m1+mt1+m2+mt2+m3+mt3+m12+mt12+m13+mt13+m23+mt23+m123+mt123),])
  
  list(Timputed,Cimputed,Ximputed)
}  

#Calculates 3 to 5 sets of initial values for K=1,2 which are then used as starting points in maximization
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
  SigmaI <- Ins[8]}
  return(cbind(Beta0I,matrix(Ins[2:7],(1+k*2),6,byrow=TRUE),SigmaI,PhiI))
}

#Optimization function for SNP given data and initial values
Optimize <- function(InitSNP,T,C,X,k,cc,mt,EXP){
  if(EXP==TRUE){
    if(k==0) opt <- optim(InitSNP,CensSNPexp(T,C,X,k,cc,mt),method="BFGS")
    if(k==1) {
      for(i in 1:dim(InitSNP)[[1]]){
        possopt <- optim(InitSNP[i,],CensSNPexp(T,C,X,k,cc,mt),method="BFGS", control=list(maxit=1000))
        if(i==1) opt <- possopt else{ if(possopt$value<opt$value & possopt$par[9]<pi/2 & possopt$par[9]>-pi/2) opt <- possopt}		
      }
    }	
    if(k==2) {
      for(i in 1:dim(InitSNP)[[1]]){
        possopt <- optim(InitSNP[i,],CensSNPexp(T,C,X,k,cc,mt),method="BFGS", control=list(maxit=1000))
        if(i==1) opt <- possopt else{ if(possopt$value<opt$value & possopt$par[9]<pi/2 & possopt$par[9]>-pi/2& possopt$par[10]< pi/2 & possopt$par[10]> -pi/2) opt <- possopt}		
      }
    }
    if(k==1 & (opt$par[9]>pi/2 | opt$par[9]< -pi/2)) opt <- optim(InitSNP[1,],CensSNPexp(T,C,X,k,cc,mt),lower = c(-Inf,-Inf,-Inf,-Inf,-Inf,-pi/2), upper = c(Inf,Inf,Inf,Inf,Inf,pi/2),method="L-BFGS-B", control=list(maxit=1000))
    if(k==2 & (opt$par[9]>pi/2 | opt$par[9]< -pi/2 | opt$par[10]>pi/2 | opt$par[10]< -pi/2)) opt <- optim(InitSNP[1,],CensSNPexp(T,C,X,k,cc,mt),lower = c(-Inf,-Inf,-Inf,-Inf,-Inf,-pi/2,-pi/2), upper = c(Inf,Inf,Inf,Inf,Inf,pi/2,pi/2),method="L-BFGS-B",control=list(maxit=1000))
  }
  if(EXP==FALSE){
    if(k==0) opt <- optim(InitSNP,CensSNPnorm(T,C,X,k,cc,mt),method="BFGS")
    if(k==1) {
      for(i in 1:dim(InitSNP)[[1]]){
        possopt <- optim(InitSNP[i,],CensSNPnorm(T,C,X,k,cc,mt),method="BFGS", control=list(maxit=1000))
        if(i==1) opt <- possopt else{ if(possopt$value<opt$value  & possopt$par[9]<pi/2 & possopt$par[9]>-pi/2) opt <- possopt}		
      }
    }	
    if(k==2) {
      for(i in 1:dim(InitSNP)[[1]]){
        possopt <- optim(InitSNP[i,],CensSNPnorm(T,C,X,k,cc,mt),method="BFGS", control=list(maxit=1000))
        if(i==1) opt <- possopt else{ if(possopt$value<opt$value & possopt$par[9]<pi/2 & possopt$par[9]>-pi/2 & possopt$par[10]<pi/2 & possopt$par[10]>-pi/2) opt <- possopt}		
      }
    }
    if(k==1 & (opt$par[9]>pi/2 | opt$par[9]< -pi/2)) opt <- optim(InitSNP[1,],CensSNPnorm(T,C,X,k,cc,mt),lower = c(-Inf,-Inf,-Inf,-Inf,-Inf,-pi/2), upper = c(Inf,Inf,Inf,Inf,Inf,pi/2),method="L-BFGS-B", control=list(maxit=1000))
    if(k==2 & (opt$par[9]>pi/2 | opt$par[9]< -pi/2 | opt$par[10]>pi/2 | opt$par[10]< -pi/2)) opt <- optim(InitSNP[1,],CensSNPnorm(T,C,X,k,cc,mt),lower = c(-Inf,-Inf,-Inf,-Inf,-Inf,-pi/2,-pi/2), upper = c(Inf,Inf,Inf,Inf,Inf,pi/2,pi/2),method="L-BFGS-B",control=list(maxit=1000))
  }
  return(opt)
}

#Function to calculate correct Standard Error for MI
ApproxVar <- function(beta,ests,I,F,M) {
  par <-rbind(beta,ests)
  BW <- matrix(0,nrow(par),nrow(par))
  for(i in 1:M) BW <- BW + (par[,i]-rowMeans(par))%*%t(par[,i]-rowMeans(par))/(M-1)
  (F + (1+1/M)*BW+BW%*%solve(F)%*%I%*%solve(F)%*%BW)
}


#Calculates consistent standard error for MI with normal kernel SNP
VarNorm <- function(Timp,Cimp,Ximp,cnone,cone,par,phi,k,M,I,F) {
  
  ScoreMatrix <- matrix(0,dim(par)[[1]],M*n)
  
  for(d in 1:M){
    X <- Ximp[,((d-1)*7+1):(d*7)]
    T <- as.vector(Timp[,d])
    C <- as.vector(Cimp[,d])
    beta <- as.vector(par[1:7,d])
    sigma <- as.numeric(par[8,d])
    if(k>0) a <- acoef(phi[,d],k,EXP=FALSE) else a <- c(1,0,0)
    if(k>0) ests <- as.vector(par[(8+dim(phi)[[1]]):dim(par)[[1]],d]) else ests <- as.vector(par[8:dim(par)[[1]],d])
    a <- c(a,rep(0,2-k))
    
    S <- (T[1:cnone]-as.vector(X[1:cnone,]%*%beta))/sigma
    PK <- (a[1]+a[2]*S+a[3]*S^2)
    for(j in 1:7) ScoreMatrix[j,(d*n+1-n):(d*n-cone)] <- -2/PK*(a[2]/sigma*X[1:cnone,j]+2*a[3]*S*X[1:cnone,j]/sigma)+S*X[1:cnone,j]/sigma
    
    S2 <- (C[(cnone+1):(cnone+cone)]-as.vector(X[(cnone+1):(cnone+cone),]%*%beta))/sigma
    for(j in 1:7) ScoreMatrix[j,(d*n+1+cnone-n):(d*n)] <- -SNPnorm(S2,a)*X[(cnone+1):(cnone+cone),j]/sigma
    ScoreMatrix[8,(d*n+1+cnone-n):(d*n)] <- -SNPnorm(S2,a)*S2/sigma
    
    B1 <- ests[1]
    B2 <- ests[2]
    B3 <- ests[3]
    B4 <- ests[4]
    B5 <- ests[5]
    B6 <- ests[6]
    B7 <- ests[7]
    B8 <- ests[8]
    B9 <- ests[9]
    B10 <- ests[10]
    B11 <- ests[11]
    B12 <- ests[12]
    S1 <- ests[13]
    S12 <- ests[14]
    S13 <- ests[15]
    S2 <- ests[16]
    S23 <- ests[17]
    S3 <- ests[18]
    ScoreMatrix[9,(d*n+1-n):(d*n)] <- (1/2)*(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/((2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2))*sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))
    ScoreMatrix[10,(d*n+1-n):(d*n)] <- (1/2)*(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))*X[,5]*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/((2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2))*sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))
    ScoreMatrix[11,(d*n+1-n):(d*n)] <- (1/2)*(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))*X[,6]*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/((2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2))*sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))
    ScoreMatrix[12,(d*n+1-n):(d*n)] <- (1/2)*(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))*X[,7]*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/((2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2))*sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))
    ScoreMatrix[13,(d*n+1-n):(d*n)] <- (1/2)*(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))*(S12*S3-S13*S23)*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/((2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2))*(S3*S2-S23^2)*sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))+(1/2)*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/((2*S2-2*S23^2/S3)*sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))
    ScoreMatrix[14,(d*n+1-n):(d*n)] <- (1/2)*(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))*(S12*S3-S13*S23)*X[,5]*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/((2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2))*(S3*S2-S23^2)*sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))+(1/2)*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)*X[,5]*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/((2*S2-2*S23^2/S3)*sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))
    ScoreMatrix[15,(d*n+1-n):(d*n)] <- (1/2)*(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))*(S12*S3-S13*S23)*X[,6]*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/((2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2))*(S3*S2-S23^2)*sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))+(1/2)*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)*X[,6]*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/((2*S2-2*S23^2/S3)*sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))
    ScoreMatrix[16,(d*n+1-n):(d*n)] <- (1/2)*(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))*(S12*S3-S13*S23)*X[,7]*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/((2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2))*(S3*S2-S23^2)*sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))+(1/2)*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)*X[,7]*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/((2*S2-2*S23^2/S3)*sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))
    ScoreMatrix[17,(d*n+1-n):(d*n)] <- (1/2)*(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))*(S13*S2-S12*S23)*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/((2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2))*(S3*S2-S23^2)*sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))+(1/2)*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)*S23*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/((2*S2-2*S23^2/S3)*S3*sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))+(1/4)*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/(S3*sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))
    ScoreMatrix[18,(d*n+1-n):(d*n)] <- (1/2)*(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))*(S13*S2-S12*S23)*X[,5]*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/((2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2))*(S3*S2-S23^2)*sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))+(1/2)*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)*S23*X[,5]*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/((2*S2-2*S23^2/S3)*S3*sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))+(1/4)*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)*X[,5]*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/(S3*sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3)) 
    ScoreMatrix[19,(d*n+1-n):(d*n)] <- (1/2)*(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))*(S13*S2-S12*S23)*X[,6]*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/((2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2))*(S3*S2-S23^2)*sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))+(1/2)*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)*S23*X[,6]*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/((2*S2-2*S23^2/S3)*S3*sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))+(1/4)*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)*X[,6]*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/(S3*sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))
    ScoreMatrix[20,(d*n+1-n):(d*n)] <- (1/2)*(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))*(S13*S2-S12*S23)*X[,7]*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/((2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2))*(S3*S2-S23^2)*sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))+(1/2)*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)*S23*X[,7]*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/((2*S2-2*S23^2/S3)*S3*sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))+(1/4)*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)*X[,7]*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/(S3*sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3)) 
    ScoreMatrix[21,(d*n+1-n):(d*n)] <- (1/2)*(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/((2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2))^2*sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))-(1/8)*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)*pi/((pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))^(3/2)*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))
    ScoreMatrix[22,(d*n+1-n):(d*n)] <- (1/4)*(-2*(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2*(2*S12*S3-2*S13*S23)/((2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2))^2*(S3*S2-S23^2))-(2*(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2)))*(S3*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)-S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/(sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))+(1/8)*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)*pi*(2*S12*S3-2*S13*S23)/((pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))^(3/2)*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3)*(S3*S2-S23^2))
    ScoreMatrix[23,(d*n+1-n):(d*n)] <- (1/4)*(-2*(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2*(-2*S12*S23+2*S13*S2)/((2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2))^2*(S3*S2-S23^2))-(2*(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2)))*(-S23*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+S2*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/(sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))+(1/8)*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)*pi*(-2*S12*S23+2*S13*S2)/((pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))^(3/2)*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3)*(S3*S2-S23^2))
    ScoreMatrix[24,(d*n+1-n):(d*n)] <- (1/4)*((X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2*(-2*S13^2/(S3*S2-S23^2)+(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))*S3/(S3*S2-S23^2)^2)/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2))^2-(2*(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2)))*(-(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)*S3/(S3*S2-S23^2)^2+S13*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2)-(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)*S3/(S3*S2-S23^2)^2)/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/(sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))+(1/2)*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/((2*S2-2*S23^2/S3)^2*sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))-(1/8)*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)*pi*(-S13^2/(S3*S2-S23^2)+(S12^2*S3-2*S13*S12*S23+S13^2*S2)*S3/(S3*S2-S23^2)^2)/((pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))^(3/2)*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))-(1/8)*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)*pi/(sqrt(pi*(S1-(S12^2*S3-2*S13*S12								*S23+S13^2*S2)/(S3*S2-S23^2)))*(pi*(S2-S23^2/S3))^(3/2)*sqrt(pi*S3))
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           ScoreMatrix[25,(d*n+1-n):(d*n)] <- (1/4)*((X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2*(-2*S13^2/(S3*S2-S23^2)+(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))*S3/(S3*S2-S23^2)^2)/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2))^2-(2*(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2)))*(-(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)*S3/(S3*S2-S23^2)^2+S13*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2)-(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)*S3/(S3*S2-S23^2)^2)/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/(sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))+(1/2)*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/((2*S2-2*S23^2/S3)^2*sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))-(1/8)*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)*pi*(-S13^2/(S3*S2-S23^2)+(S12^2*S3-2*S13*S12*S23+S13^2*S2)*S3/(S3*S2-S23^2)^2)/((pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))^(3/2)*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))-(1/8)*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)*pi/(sqrt(pi*(S1-(S12^2*S3-2*S13*S12								*S23+S13^2*S2)/(S3*S2-S23^2)))*(pi*(S2-S23^2/S3))^(3/2)*sqrt(pi*S3))
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  ScoreMatrix[26,(d*n+1-n):(d*n)] <- (1/4)*((X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2*(4*S13*S12/(S3*S2-S23^2)-(4*(S12^2*S3-2*S13*S12*S23+S13^2*S2))*S23/(S3*S2-S23^2)^2)/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2))^2-(2*(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2)))*(-S13*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(2*(S12*S3-S13*S23))*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)*S23/(S3*S2-S23^2)^2-S12*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2)+(2*(S13*S2-S12*S23))*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)*S23/(S3*S2-S23^2)^2)/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/(sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))+(1/4)*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*(-4*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2*S23/((2*S2-2*S23^2/S3)^2*S3)-(2*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3))*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/((2*S2-2*S23^2/S3)*S3))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)/(sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))-(1/8)*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)*pi*(2*S13*S12/(S3*S2-S23^2)-(2*(S12^2*S3-2*S13*S12*S23+S13^2*S2))*S23/(S3*S2-S23^2)^2)/((pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))^(3/2)*sqrt(pi*(S2-S23^2/S3))*sqrt(pi*S3))+(1/4)*exp(-(X[,2]-B1-X[,5]*B2-X[,6]*B3-X[,7]*B4+(S12*S3-S13*S23)*(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8)/(S3*S2-S23^2)+(S13*S2-S12*S23)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/(S3*S2-S23^2))^2/(2*S1-(2*(S12^2*S3-2*S13*S12								*S23+S13^2*S2))/(S3*S2-S23^2)))*exp(-(X[,3]-B5-X[,5]*B6-X[,6]*B7-X[,7]*B8+S23*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)/S3)^2/(2*S2-2*S23^2/S3))*exp(-(1/2)*(X[,4]-B9-X[,5]*B10-X[,6]*B11-X[,7]*B12)^2/S3)*sqrt(2)*pi*S23/(sqrt(pi*(S1-(S12^2*S3-2*S13*S12*S23+S13^2*S2)/(S3*S2-S23^2)))*(pi*(S2-S23^2/S3))^(3/2)*sqrt(pi*S3)*S3)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
  }
  
  IFmIN <- matrix(0,26,26)
  for(i in 1:n){
    Means <- rowMeans(ScoreMatrix[,seq(i,M*n,n)])
    for(d in 1:M){
      IFmIN <- IFmIN + (ScoreMatrix[,(n*(d-1)+i)]-Means)%*%t((ScoreMatrix[,(n*(d-1)+i)]-Means)) 
    }
  }
  
  IFmIN <- IFmIN/(M-1)
  F + (1+1/M)*F%*%IFmIN%*%F + F%*%IFmIN%*%I%*%IFmIN%*%F
}

#Calculates consistent standard error for MI with exponential kernel SNP
VarExp <- function(Timp,Cimp,Ximp,cnone,cone,par,phi,k,M,I,F) {
  
  ScoreMatrix <- matrix(0,dim(par)[[1]],M*n)
  
  for(d in 1:M){
    X <- Ximp[,((d-1)*7+1):(d*7)]
    T <- as.vector(Timp[,d])
    C <- as.vector(Cimp[,d])
    beta <- as.vector(par[1:7,d])
    sigma <- as.numeric(par[8,d])
    if(k>0) a <-acoef(phi[,d],k,EXP=TRUE) else a <- c(1,0,0)
    ests <- as.vector(par[9:dim(par)[[1]],d])
    
    S <- exp((T[1:cnone]-as.vector(X[1:cnone,]%*%beta))/exp(sigma))
    SB <- ((T[1:cnone]-as.vector(X[1:cnone,]%*%beta))/exp(2*sigma))
    PK <- (a[1]+a[2]*S+a[3]*S^2)
    V <- (X[,2]-ests[1]-ests[2]*X[,4]-ests[6]/ests[7]*(X[,3]-ests[3]-ests[4]*X[,4]))
    Sig <- (ests[5]-ests[6]^2/ests[7])
    Num <- (X[,3]-ests[3]-ests[4]*X[,4])
    
    for(j in 1:7) ScoreMatrix[j,(d*n+1-n):(d*n-cone)] <- -X[1:cnone,j]/exp(sigma)+S*X[1:cnone,j]/exp(sigma)-2/PK*(a[2]*S*X[1:cnone,j]/exp(sigma)+a[3]*S^2*2*X[1:cnone,j]/exp(sigma))
    ScoreMatrix[8,(d*n+1-n):(d*n-cone)] <- -1 - SB+S*SB - 2/PK*(a[2]*S*SB+2*a[3]*S^2*SB)
    ScoreMatrix[6,(d*n+1-n):(d*n)] <- V/Sig
    ScoreMatrix[7,(d*n+1-n):(d*n)] <- X[,4]*V/Sig
    ScoreMatrix[8,(d*n+1-n):(d*n)] <- -V/Sig*ests[6]^2/ests[7]+Num/ests[7]
    ScoreMatrix[9,(d*n+1-n):(d*n)] <- -V/Sig*X[,4]*ests[6]^2/ests[7]+X[,4]*Num/ests[7]
    ScoreMatrix[10,(d*n+1-n):(d*n)] <- -1/(2*Sig)+V^2/(2*Sig^2)
    ScoreMatrix[11,(d*n+1-n):(d*n)] <- -ests[6]^2/ests[7]^2/(2*Sig)-1/(2*ests[7])-ests[6]/ests[7]^2*V/Sig+V^2*ests[6]^2/ests[7]^2/(2*Sig^2)+Num^2/(2*ests[7]^2)
    ScoreMatrix[12,(d*n+1-n):(d*n)] <- ests[6]/ests[7]/Sig + V*1/ests[7]^2/Sig-V^2*ests[6]/ests[7]/Sig^2
    
    S2 <- exp((C[(cnone+1):(cnone+cone)]-as.vector(X[(cnone+1):(cnone+cone),]%*%beta))/exp(sigma))
    S2B <- ((T[(cnone+1):(cnone+cone)]-as.vector(X[(cnone+1):(cnone+cone),]%*%beta))/exp(2*sigma))
    for(j in 1:7) ScoreMatrix[j,(d*n+1+cnone-n):(d*n)] <- -SNPexp(S2,a)*S2*X[(cnone+1):(cnone+cone),j]/exp(sigma)
    ScoreMatrix[8,(d*n+1+cnone-n):(d*n)] <- -SNPnorm(S2,a)*S2*S2B
  }
  
  IFmIN <- matrix(0,26,26)
  for(i in 1:n){
    Means <- rowMeans(ScoreMatrix[,seq(i,M*n,n)])
    for(d in 1:M){
      IFmIN <- IFmIN + (ScoreMatrix[,(n*(d-1)+i)]-Means)%*%t((ScoreMatrix[,(n*(d-1)+i)]-Means)) 
    }
  }
  
  IFmIN <- IFmIN/(M-1)
  F + (1+1/M)*F%*%IFmIN%*%F + F%*%IFmIN%*%I%*%IFmIN%*%F
}

#functions to calculate cdf in order to be able to use *apply* function, eliminating the need for a loop which slows maximization
Int1 <- function(X,t,sig12,Lower) sadmvn(lower=Lower,upper=X[8:9],mean=(c(t[1]+t[2]*X[5]+t[3]*X[6]+t[4]*X[7],t[5]+t[6]*X[5]+t[7]*X[6]+t[8]*X[7])+as.vector(matrix(c(t[15],t[17]),2,1)*1/t[18]*(X[4]-t[9]-t[10]*X[5]-t[11]*X[6]-t[12]*X[7]))),varcov=make.positive.definite(sig12-1/t[18]*matrix(c(t[15],t[17]),2,1)%*%t(matrix(c(t[15],t[17]),2,1))),abseps=0.000000001)
Int2 <- function(X,t,sig13,Lower,Upper) sadmvn(lower=Lower,upper=Upper,mean=(c(t[1]+t[2]*X[5]+t[3]*X[6]+t[4]*X[7],t[9]+t[10]*X[5]+t[11]*X[6]+t[12]*X[7])+as.vector(matrix(c(t[14],t[17]),2,1)*1/t[16]*(X[3]-t[5]-t[6]*X[5]-t[7]*X[6]-t[8]*X[7]))),varcov=make.positive.definite(sig13-1/t[16]*matrix(c(t[14],t[17]),2,1)%*%t(matrix(c(t[14],t[17]),2,1))),abseps=0.000000001)
Int3 <- function(X,t,sig23,Lower) sadmvn(lower=Lower,upper=X[8:9],mean=(c(t[5]+t[6]*X[5]+t[7]*X[6]+t[8]*X[7],t[9]+t[10]*X[5]+t[11]*X[6]+t[12]*X[7])+as.vector(matrix(c(t[14],t[15]),2,1)*1/t[13]*(X[2]-t[1]-t[2]*X[5]-t[3]*X[6]-t[4]*X[7]))),varcov=make.positive.definite(sig23-1/t[13]*matrix(c(t[14],t[15]),2,1)%*%t(matrix(c(t[14],t[15]),2,1))),abseps=0.000000001)
Int4 <- function(X,t,sig123,Lower) sadmvn(lower=Lower,upper=X[8:10],mean=c(t[1]+t[2]*X[5]+t[3]*X[6]+t[4]*X[7],t[5]+t[6]*X[5]+t[7]*X[6]+t[8]*X[7],t[9]+t[10]*X[5]+t[11]*X[6]+t[12]*X[7]),varcov=make.positive.definite(sig123),abseps=0.000000001)

#Gets c coefficients, then a coefficients
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
  if(EXP==TRUE) for(i in 1:(k+1)) s[,i] <- (exp((T-as.real(X%*%c(t[1],t[2],t[3],t[4],t[5],t[6],t[7])))/exp(t[8])))^(i-1)
  if(EXP==FALSE) for(i in 1:(k+1)) s[,i] <- ((T-as.real(X%*%c(t[1],t[2],t[3],t[4],t[5],t[6],t[7])))/t[8])^(i-1)
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
  Beta <- InitPar[1:7]
  Rest <- InitPar[8:length(InitPar)]
  if(dist==1){if(C==0) gen <- EV(T,loc=X%*%Beta) #standard extreme value (exponential model for T)
  if(C!=0) gen <- pgumbel(-C,loc=-X%*%Beta)}
  if(dist==2){if(C==0)  gen <- EV(T,loc=X%*%Beta,scale=Rest) #Gumbel/log-weibull  (exponential model for T)
  if(C!=0) gen <- pgumbel(-C,loc=-X%*%Beta,scale=Rest)}
  if(dist==3){if(C==0) gen <- dnorm(T,mean=X%*%Beta,Rest)   #normal (log-normal model for T)
  if(C!=0) gen <- pnorm(C,mean=X%*%Beta,Rest,lower.tail=FALSE) }
  if(dist==4){if(C==0)  gen <- logistic(T,loc=X%*%Beta,scale=Rest) #logistic (log-logisitic model for T)
  if(C!=0) gen <-plogis(C,location=X%*%Beta,Rest,lower.tail=FALSE) }
  
  gen
}


#CDF fucionts used in GenTsnp below (needed for MI algorithm)
pSNPexp <- function(x,a) a[1]^2*pgamma(x,1,1,lower.tail=FALSE)+2*a[1]*a[2]*pgamma(x,2,1,lower.tail=FALSE)+2*(a[2]^2+2*a[1]*a[3])*pgamma(x,3,1,lower.tail=FALSE)+12*a[2]*a[3]*pgamma(x,4,1,lower.tail=FALSE)+24*a[3]^2*pgamma(x,5,1,lower.tail=FALSE)
pSNPnorm <-function(x,a) a[1]^2*pnorm(x,lower.tail=FALSE)+2*a[1]*a[2]*dnorm(x)+(a[2]^2+2*a[1]*a[3])*(x*dnorm(x)+pnorm(x,lower.tail=FALSE))+2*a[2]*a[3]*(x^2+2)*dnorm(x)+a[3]^2*(x^3*dnorm(x)+3*x*dnorm(x)+3*pnorm(x,lower.tail=FALSE))

#Used in MH-algorithm for imputing values when assuming an SNP model
GenTsnp <- function(InitSNP,k,a,EXP=TRUE,T,C,X){
  
  Beta <- InitSNP[1:7]
  Sigma <- InitSNP[8]
  if(k==0) a1 <- c(a,0,0)
  if(k==1) a1 <- c(a,0)
  if(k==2) a1 <- a
  if(EXP==TRUE){if(T<=C) gen <- 1/(exp(Sigma))*exp((T-as.real(X%*%Beta))/exp(Sigma))*(Pk(T,X,InitSNP,k,a,1))^2*dexp(exp((T-as.real(X%*%Beta))/exp(Sigma)))
  if(T>C) gen <- pSNPexp(min(exp((C-as.real(X%*%Beta))/exp(Sigma)),1e50),a1)
  }
  if(EXP==FALSE){if(T<=C) gen <- 1/Sigma*(Pk(T,X,InitSNP,k,a,1,EXP=FALSE))^2*dnorm((T-as.real(X%*%Beta))/Sigma)
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

#Likelihood Values at Initial Values
IV <- function(Phi,T1,C1,X1,Ins,k,EXP=TRUE,cc,mt,FINAL=FALSE) {
  BetaI <- Ins[1:7]
  TauI <- Ins[8]
  a <- c(acoef(Phi,k,EXP), rep(0, 2-k)) 
  if(FINAL==FALSE){
    if(EXP==TRUE) SigmaI <- log(TauI/sqrt(integrate(SNPlog2,lower=0,upper=Inf,stop.on.error=FALSE,a)$value-(integrate(SNPlog,lower=0,upper=Inf,stop.on.error=FALSE,a)$value)^2))
    if(EXP==FALSE) SigmaI <- TauI/sqrt(integrate(SNPx2,lower=-Inf,upper=Inf,stop.on.error=FALSE,a)$value-(integrate(SNPx,lower=-Inf,upper=Inf,stop.on.error=FALSE,a)$value)^2)
    Beta0I <- Beta0If(Phi,SigmaI,BetaI,a,EXP,FINAL)
  }
  if(FINAL==TRUE){EstPhi <- Ins[9:length(Ins)]
  a2 <- c(acoef(EstPhi,k,EXP), rep(0, 2-k)) 
  if(EXP==TRUE) SigmaI <- log(exp(TauI)*sqrt((integrate(SNPlog2,lower=0,upper=Inf,stop.on.error=FALSE,a2)$value-(integrate(SNPlog,lower=0,upper=Inf,stop.on.error=FALSE,a2)$value)^2)/(integrate(SNPlog2,lower=0,upper=Inf,stop.on.error=FALSE,a)$value-(integrate(SNPlog,lower=0,upper=Inf,stop.on.error=FALSE,a)$value)^2)))
  if(EXP==FALSE) SigmaI <- (TauI*sqrt((integrate(SNPx2,lower=-Inf,upper=Inf,stop.on.error=FALSE,a2)$value-(integrate(SNPx,lower=-Inf,upper=Inf,stop.on.error=FALSE,a2)$value)^2)/(integrate(SNPx2,lower=-Inf,upper=Inf,stop.on.error=FALSE,a)$value-(integrate(SNPx,lower=-Inf,upper=Inf,stop.on.error=FALSE,a)$value)^2)))
  Beta0I <- Beta0If(Phi,SigmaI,BetaI,a,EXP,FINAL,a2)
  }
  Init(c(Beta0I,BetaI[2:7],SigmaI,Phi),T1,C1,X1,k,cc,mt,EXP)
}

#Calculates likelihood Values at Initial Values and outputs only updated initial values
IV2 <- function(Phi,T1,C1,X1,Ins,k,EXP=TRUE,cc,mt,FINAL=FALSE) {
  BetaI <- Ins[1:7]
  TauI <- Ins[8]
  a <- c(acoef(Phi,k,EXP), rep(0, 2-k)) 
  if(FINAL==FALSE){
    if(EXP==TRUE) SigmaI <- log(TauI/sqrt(integrate(SNPlog2,lower=0,upper=Inf,stop.on.error=FALSE,a)$value-(integrate(SNPlog,lower=0,upper=Inf,stop.on.error=FALSE,a)$value)^2))
    if(EXP==FALSE) SigmaI <- TauI/sqrt(integrate(SNPx2,lower=-Inf,upper=Inf,stop.on.error=FALSE,a)$value-(integrate(SNPx,lower=-Inf,upper=Inf,stop.on.error=FALSE,a)$value)^2)
    Beta0I <- Beta0If(Phi,SigmaI,BetaI,a,EXP,FINAL)
  }
  if(FINAL==TRUE){EstPhi <- Ins[9:length(Ins)]
  a2 <- c(acoef(EstPhi,k,EXP), rep(0, 2-k)) 
  if(EXP==TRUE) SigmaI <- log(exp(TauI)*sqrt((integrate(SNPlog2,lower=0,upper=Inf,stop.on.error=FALSE,a2)$value-(integrate(SNPlog,lower=0,upper=Inf,stop.on.error=FALSE,a2)$value)^2)/(integrate(SNPlog2,lower=0,upper=Inf,stop.on.error=FALSE,a)$value-(integrate(SNPlog,lower=0,upper=Inf,stop.on.error=FALSE,a)$value)^2)))
  if(EXP==FALSE) SigmaI <- (TauI*sqrt((integrate(SNPx2,lower=-Inf,upper=Inf,stop.on.error=FALSE,a2)$value-(integrate(SNPx,lower=-Inf,upper=Inf,stop.on.error=FALSE,a2)$value)^2)/(integrate(SNPx2,lower=-Inf,upper=Inf,stop.on.error=FALSE,a)$value-(integrate(SNPx,lower=-Inf,upper=Inf,stop.on.error=FALSE,a)$value)^2)))
  Beta0I <- Beta0If(Phi,SigmaI,BetaI,a,EXP,FINAL,a2)
  }
  list(Beta0I,SigmaI)
}

Init <- function(t,T,C,X,k,cnone,cone,EXP=TRUE) {
  fc<-rep(0,(cnone+cone))
  if(EXP==TRUE){a <- acoef(t[9:length(t)],k,EXP=TRUE)
  P_k <- Pk(T[1:cnone],X[1:cnone,],t,k,a,cnone,EXP=TRUE)
  a <- c(a, rep(0, 3-length(a))) 
  if(cnone>0) fc[1:cnone] <- t[8]-(T[1:cnone]-as.vector(X[1:cnone,]%*%c(t[1],t[2],t[3],t[4],t[5],t[6],t[7])))/exp(t[8])-log((P_k[1:cnone])^2)+exp((T[1:cnone]-as.vector(X[1:cnone,]%*%c(t[1],t[2],t[3],t[4],t[5],t[6],t[7])))/exp(t[8]))
  if(cone>0) for(i in (cnone+1):(cnone+cone)) fc[i] <- -log(tryCatch(integrate(SNPexp,lower=min(exp((C[i]-as.real(X[i,]%*%c(t[1],t[2],t[3],t[4],t[5],t[6],t[7])))/exp(t[8])),1e50),upper=Inf,stop.on.error=FALSE,rel.tol = .Machine$double.eps^0.5,a)$value, error=function(...) 1))}		    
  if(EXP==FALSE){a <- acoef(t[9:length(t)],k,EXP=FALSE)
  P_k <- Pk(T[1:cnone],X[1:cnone,],t,k,a,cnone,EXP=FALSE)
  a <- c(a, rep(0, 3-length(a))) 
  if(cnone>0) fc[1:cnone] <-  -log(1/t[8]*(P_k[1:cnone])^2*dnorm(((T[1:cnone]-as.vector(X[1:cnone,]%*%c(t[1],t[2],t[3],t[4],t[5],t[6],t[7])))/t[8])))
  if(cone>0) for(i in (cnone+1):(cnone+cone)) fc[i] <- -log(tryCatch(integrate(SNPnorm,lower=min(((C[i]-as.real(X[i,]%*%c(t[1],t[2],t[3],t[4],t[5],t[6],t[7])))/t[8]),1e50),upper=Inf,stop.on.error=FALSE,rel.tol = .Machine$double.eps^0.5,a)$value, error=function(...) 1))}					    
  sum(fc)
}
