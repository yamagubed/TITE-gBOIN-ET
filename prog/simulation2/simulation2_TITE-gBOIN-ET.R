
library(Iso)

source("../function/myfunction.R")

toxeffprob <- read.csv("../../data/true_probbaility.csv")
toxprob    <- toxeffprob[(toxeffprob$Item=="Tox")&(toxeffprob$Grade!="DLT")&(toxeffprob$Grade!="ETS")&(toxeffprob$Grade!="nETS"),]
effprob    <- toxeffprob[(toxeffprob$Item=="Eff")&(toxeffprob$Grade!="ORR")&(toxeffprob$Grade!="EES")&(toxeffprob$Grade!="nEES"),]

ls1   <- length(unique(toxprob$Scenario))
lds   <- 6
dosen <- 1:lds
dose  <- paste("Dose",1:lds,sep="")

ncoh <- 3
nmax <- 36
nesc <- nmax/ncoh

lambdaeta <- boundary(phi=0.313,delta=0.583)
phi       <- lambdaeta$phi
phi1      <- phi*0.1
phi2      <- phi*1.4
delta     <- lambdaeta$delta
delta1    <- delta*0.6
lambda1   <- lambdaeta$lambda1
lambda2   <- lambdaeta$lambda2
eta1      <- lambdaeta$eta1

weight.T <- c(0,0.5,1,1.5)
weight.E <- c(0,0.25,1,3)

tau.T <- 30
tau.E <- c(45,60,90)
ls2   <- length(tau.E)

accrual <- c(5,10)
ls3     <- length(accrual)

pr.alpha <- 1
pr.beta  <- 1

alpha.T1 <- 0.5
alpha.T2 <- 0.5
alpha.E1 <- 0.5
alpha.E2 <- 0.5

lss   <- 10000
maxev <- 50

data.obs.n   <- array(0,dim=c(lss,ls1,ls2,ls3,lds))
data.obs.tox <- array(0,dim=c(lss,ls1,ls2,ls3,lds))
data.obs.eff <- array(0,dim=c(lss,ls1,ls2,ls3,lds))
data.dur     <- array(0,dim=c(lss,ls1,ls2,ls3))

set.seed(100)

for(ss in 1:lss){
for(s1 in 1:ls1){
for(s2 in 1:ls2){
for(s3 in 1:ls3){

 efftoxp <- list(effp=effprob[effprob$Scenario==s1,dose],toxp=toxprob[toxprob$Scenario==s1,dose])

 obs.n     <- numeric(lds)
 obs.tox   <- numeric(lds)
 obs.tox.n <- numeric(lds)
 obs.eff   <- numeric(lds)
 obs.eff.n <- numeric(lds)
 pe        <- numeric(lds)
 pt        <- numeric(lds)

 wk.time.tox <- array(0,dim=c(nmax,maxev))
 wk.event.T  <- array(0,dim=c(nmax,maxev))
 wk.grade    <- array(0,dim=c(nmax,maxev))
 wk.endtox   <- array(0,dim=c(nmax,maxev))

 wk.time.eff <- array(0,dim=c(nmax,maxev))
 wk.event.E  <- array(0,dim=c(nmax,maxev))
 wk.response <- array(0,dim=c(nmax,maxev))
 wk.endeff   <- array(0,dim=c(nmax,maxev))

 nevent.T <- numeric(nmax)
 nevent.E <- numeric(nmax)
 subjdose <- numeric(nmax)
 enter    <- numeric(nmax)

 grade    <- numeric(nmax)
 ETS      <- numeric(nmax)
 nETS     <- numeric(nmax)
 censor.T <- numeric(nmax)

 response <- numeric(nmax)
 EES      <- numeric(nmax)
 nEES     <- numeric(nmax)
 censor.E <- numeric(nmax)

 t.enter    <- numeric(ncoh)
 t.decision <- 0

 curdose <- 1

 for(i in 1:nesc){

  dlab <- paste("Dose",curdose,sep="")
  obs.n[curdose] <- obs.n[curdose] + ncoh

  for(j in 1:ncoh){
   if(j==1){
    t.enter[j] <- t.decision
   }else{
    t.enter[j] <- t.enter[j-1]+runif(1,0,2*accrual[s3])
  }}

  if(i==nesc){
   t.decision <- t.enter[length(t.enter)]+max(tau.T,tau.E[s2])
  }else{
   t.decision <- t.enter[length(t.enter)]+runif(1,0,2*accrual[s3])
  }

  psi.T    <- sum(efftoxp$toxp[2:4,dlab])
  zetta.T1 <- log(log(1-psi.T)/log(1-psi.T+alpha.T1*psi.T))/log(1/(1-alpha.T2))
  zetta.T2 <- tau.T/(-log(1-psi.T))^(1/zetta.T1)
  psi.E    <- sum(efftoxp$effp[2:4,dlab])
  zetta.E1 <- log(log(1-psi.E)/log(1-psi.E+alpha.E1*psi.E))/log(1/(1-alpha.E2))
  zetta.E2 <- tau.E[s2]/(-log(1-psi.E))^(1/zetta.E1)

  for(j in 1:ncoh){
   subjid <- i*ncoh-3+j
   enter[subjid]    <- t.enter[j]
   subjdose[subjid] <- curdose

   for(k in 1:maxev){
    if(k==1){
     wk.time.tox[subjid,k] <- rweibull(1,shape=zetta.T1,scale=zetta.T2)
    }else{
     wk.time.tox[subjid,k] <- wk.time.tox[subjid,k-1]+rweibull(1,shape=zetta.T1,scale=zetta.T2)
    }
    wk.event.T[subjid,k] <- as.numeric(wk.time.tox[subjid,k]<=tau.T)
    wk.grade[subjid,k]   <- wk.event.T[subjid,k]*((1:3)%*%rmultinom(1,1,efftoxp$toxp[2:4,dlab]))[1,1]+1
    if(wk.event.T[subjid,k]==0){
     nevent.T[subjid] <- k
     break
    }
   }
   wk.time.tox[subjid,nevent.T[subjid]] <- tau.T
   wk.endtox[subjid,1:nevent.T[subjid]] <- t.enter[j]+wk.time.tox[subjid,1:nevent.T[subjid]]

   for(k in 1:maxev){
    if(k==1){
     wk.time.eff[subjid,k] <- rweibull(1,shape=zetta.E1,scale=zetta.E2)
    }else{
     wk.time.eff[subjid,k] <- wk.time.eff[subjid,k-1]+rweibull(1,shape=zetta.E1,scale=zetta.E2)
    }
    wk.event.E[subjid,k]  <- as.numeric(wk.time.eff[subjid,k]<=tau.E[s2])
    wk.response[subjid,k] <- wk.event.E[subjid,k]*((1:3)%*%rmultinom(1,1,efftoxp$effp[2:4,dlab]))[1,1]+1
    if(wk.event.E[subjid,k]==0){
     nevent.E[subjid] <- k
     break
    }
   }
   wk.time.eff[subjid,nevent.E[subjid]] <- tau.E[s2]
   wk.endeff[subjid,1:nevent.E[subjid]] <- t.enter[j]+wk.time.eff[subjid,1:nevent.E[subjid]]

  }

  tite.flg <- (subjdose==curdose)
  cu.endtox   <- wk.endtox[tite.flg,]
  cu.nevent.E <- nevent.E[tite.flg]
  cu.endeff   <- wk.endeff[tite.flg,]
  cu.nevent.T <- nevent.T[tite.flg]

  cu.minendtox <- apply(as.matrix(1:sum(tite.flg)),1,function(x){min(cu.endtox[x,1:cu.nevent.T[x]])})
  cu.minendeff <- apply(as.matrix(1:sum(tite.flg)),1,function(x){min(cu.endeff[x,1:cu.nevent.E[x]])})

  gamma.T <- as.numeric(cu.minendtox<=t.decision)
  gamma.E <- as.numeric(cu.minendeff<=t.decision)

  while(mean(gamma.T*gamma.E)<0.5){
   nextdec.T  <- cu.minendtox[cu.minendtox>t.decision]
   nextdec.E  <- cu.minendeff[cu.minendeff>t.decision]
   t.decision <- min(nextdec.T,nextdec.E)
   gamma.T    <- as.numeric(cu.minendtox<=t.decision)
   gamma.E    <- as.numeric(cu.minendeff<=t.decision)
  }

  for(ii in 1:subjid){
   censor.T[ii] <- as.numeric(min(wk.endtox[ii,1:nevent.T[ii]])>t.decision)
   if(censor.T[ii]==0){
    flg.T        <- (wk.endtox[ii,]<=t.decision)&((1:maxev)<=nevent.T[ii])
    grade[ii]    <- max((wk.grade[ii,])[flg.T])
    ETS[ii]      <- weight.T[grade[ii]]
    nETS[ii]     <- ETS[ii]/max(weight.T)
   }

   censor.E[ii] <- as.numeric(min(wk.endeff[ii,1:nevent.E[ii]])>t.decision)
   if(censor.E[ii]==0){
    flg.E        <- (wk.endeff[ii,]<=t.decision)&((1:maxev)<=nevent.E[ii])
    response[ii] <- max((wk.response[ii,])[flg.E])
    EES[ii]      <- weight.E[response[ii]]
    nEES[ii]     <- EES[ii]/max(weight.E)
   }
  }

  for(ds in 1:lds){
   if(sum(subjdose==ds)>0){

    dsflg <- (subjdose==ds)

    x.nETS <- sum(nETS[dsflg&(censor.T==0)])
    n.nETS <- x.nETS+sum(1-nETS[dsflg&(censor.T==0)])+sum(t.decision-enter[dsflg&(censor.T==1)])/tau.T

    x.nEES <- sum(nEES[dsflg&(censor.E==0)])
    n.nEES <- x.nEES+sum(1-nEES[dsflg&(censor.E==0)])+sum(t.decision-enter[dsflg&(censor.E==1)])/tau.E[s2]

    obs.tox[ds]   <- x.nETS
    obs.tox.n[ds] <- n.nETS
    pt[ds]        <- x.nETS/n.nETS

    obs.eff[ds]   <- x.nEES
    obs.eff.n[ds] <- n.nEES
    pe[ds]        <- x.nEES/n.nEES

  }}

  if((pt[curdose]<=lambda1)&(pe[curdose]<=eta1)){
   nxtdose <- curdose+1
  }else if((pt[curdose]<lambda2)&(pe[curdose]>eta1)){
   nxtdose <- curdose
  }else if(pt[curdose]>=lambda2){
   nxtdose <- curdose-1

  }else if((pt[curdose]>lambda1)&(pt[curdose]<lambda2)&(pe[curdose]<=eta1)){

   if(curdose==lds){
    three   <- c(curdose-1,curdose)
    maxpe   <- max(pe[three])
    nxtdose <- sample(dosen[which((pe==maxpe)&(is.element(dosen,three)))],1)

   }else if(obs.n[curdose+1]==0){
    nxtdose <- curdose+1

   }else if(curdose==1){
    three   <- c(curdose,curdose+1)
    maxpe   <- max(pe[three])
    nxtdose <- sample(dosen[which((pe==maxpe)&(is.element(dosen,three)))],1)

   }else{
    three   <- c(curdose-1,curdose,curdose+1)
    maxpe   <- max(pe[three])
    nxtdose <- sample(dosen[which((pe==maxpe)&(is.element(dosen,three)))],1)
   }

  }

  tpo.shape1 <- pr.alpha + obs.tox
  tpo.shape2 <- pr.beta  + (obs.tox.n-obs.tox)
  tterm      <- pbeta(phi,tpo.shape1,tpo.shape2)

  epo.shape1 <- pr.alpha + obs.eff
  epo.shape2 <- pr.beta  + (obs.eff.n-obs.eff)
  eterm      <- 1-pbeta(delta1,epo.shape1,epo.shape2)

  admflg  <- !((eterm<0.05)|(tterm<0.05))
  admdose <- dosen[admflg]

  if(sum(admflg)==0){
   break
  }else{

   if(nxtdose==0){
    if(admflg[1]){
     curdose <- 1
    }else{
     break
    }

   }else if(nxtdose==(lds+1)){
    curdose <- lds

   }else if(is.element(nxtdose,admdose)){
    curdose <- nxtdose

   }else if(curdose<nxtdose){
    if(sum(admdose>=nxtdose)!=0){
     curdose <- min(admdose[admdose>=nxtdose])
    }
   }else if(curdose>=nxtdose){
    if(sum(admdose<=nxtdose)!=0){
     curdose <- max(admdose[admdose<=nxtdose])
    }else{
     break
    }
   }

  }

 }

 data.obs.n[ss,s1,s2,s3,]   <- obs.n
 data.obs.tox[ss,s1,s2,s3,] <- obs.tox
 data.obs.eff[ss,s1,s2,s3,] <- obs.eff
 data.dur[ss,s1,s2,s3]      <- t.decision

}}}}

#####

bi.pr.alpha <- 1
bi.pr.beta  <- 1

res <- array(0,dim=c(lss,ls1,ls2,ls3,4))

for(ss in 1:lss){
for(s1 in 1:ls1){
for(s2 in 1:ls2){
for(s3 in 1:ls3){

 bi.prob.pe <- numeric(lds)
 bi.prob.pt <- numeric(lds)

 obsn <- data.obs.n[ss,s1,s2,s3,]
 obst <- data.obs.tox[ss,s1,s2,s3,]
 obse <- data.obs.eff[ss,s1,s2,s3,]
 evadose <- dosen[data.obs.n[ss,s1,s2,s3,]!=0]

 obspe <- obse[evadose]/obsn[evadose]
 obspt <- obst[evadose]/obsn[evadose]

 for(i in evadose){
  po.shape1 <- bi.pr.alpha + obse[i]
  po.shape2 <- bi.pr.beta  + (obsn[i]-obse[i])
  bi.prob.pe[i] <- 1-pbeta(delta1,po.shape1,po.shape2)

  po.shape1 <- bi.pr.alpha + obst[i]
  po.shape2 <- bi.pr.beta  + (obsn[i]-obst[i])
  bi.prob.pt[i] <- pbeta(phi,po.shape1,po.shape2)
 }

 if(sum(obsn)<nmax){

   res[ss,s1,s2,s3,] <- 0

 }else if(length(evadose)==1){

  if((bi.prob.pe[evadose]>=0.05)&(bi.prob.pt[evadose]>=0.05)){
   res[ss,s1,s2,s3,] <- evadose
  }

 }else if(sum((bi.prob.pe[evadose]>=0.05)&(bi.prob.pt[evadose]>=0.05))>=1){

  iso.obspt <- pava(obspt[evadose])
  res[ss,s1,s2,s3,] <- obdselect(obspe,iso.obspt,phi,phi2,bi.prob.pe[evadose],bi.prob.pt[evadose])

 }

}}}}

sprop <- array(0,dim=c(ls1,ls2,ls3,4,7))

for(s1 in 1:ls1){
for(s2 in 1:ls2){
for(s3 in 1:ls3){

 sprop[s1,s2,s3,,1] <- apply(res[,s1,s2,s3,]==0,2,mean)*100
 sprop[s1,s2,s3,,2] <- apply(res[,s1,s2,s3,]==1,2,mean)*100
 sprop[s1,s2,s3,,3] <- apply(res[,s1,s2,s3,]==2,2,mean)*100
 sprop[s1,s2,s3,,4] <- apply(res[,s1,s2,s3,]==3,2,mean)*100
 sprop[s1,s2,s3,,5] <- apply(res[,s1,s2,s3,]==4,2,mean)*100
 sprop[s1,s2,s3,,6] <- apply(res[,s1,s2,s3,]==5,2,mean)*100
 sprop[s1,s2,s3,,7] <- apply(res[,s1,s2,s3,]==6,2,mean)*100

}}}

#####

method <- "TITE-gBOIN-ET"
nobd   <- as.numeric(apply(toxeffprob[toxeffprob$Grade=="Utility",dose],1,which.max))

obdm <- c(1,3)
ls4  <- length(obdm)

pcs.s1 <- array(0,dim=c(ls1,ls2,ls3,ls4))
pca.s1 <- array(0,dim=c(ls1,ls2,ls3,ls4))
pos.s1 <- array(0,dim=c(ls1,ls2,ls3,ls4))
poa.s1 <- array(0,dim=c(ls1,ls2,ls3,ls4))
pet.s1 <- array(0,dim=c(ls1,ls2,ls3,ls4))
dur.s1 <- array(0,dim=c(ls1,ls2,ls3,ls4))

for(s1 in 1:ls1){
for(s2 in 1:ls2){
for(s3 in 1:ls3){
for(s4 in 1:ls4){

 pcs.s1[s1,s2,s3,s4] <- mean(res[,s1,s2,s3,obdm[s4]]==nobd[s1])*100
 pca.s1[s1,s2,s3,s4] <- sum(data.obs.n[,s1,s2,s3,nobd[s1]])/sum(data.obs.n[,s1,s2,s3,])*100
 if(nobd[s1]<max(dosen)){
  pos.s1[s1,s2,s3,s4] <- mean((res[,s1,s2,s3,obdm[s4]]!=0)&(res[,s1,s2,s3,obdm[s4]]>nobd[s1]))*100
  poa.s1[s1,s2,s3,s4] <- sum(data.obs.n[,s1,s2,s3,(nobd[s1]+1):lds])/sum(data.obs.n[,s1,s2,s3,])*100
 }
 pet.s1[s1,s2,s3,s4] <- mean(res[,s1,s2,s3,obdm[s4]]==0)*100
 dur.s1[s1,s2,s3,s4] <- mean(data.dur[,s1,s2,s3])

}}}}

obdml <- c("Utility","Highest Efficacy")

out.df <- NULL
for(s2 in 1:ls2){
for(s3 in 1:ls3){
for(s4 in 1:ls4){

 out    <- cbind(pcs.s1[,s2,s3,s4],pca.s1[,s2,s3,s4],pos.s1[,s2,s3,s4],poa.s1[,s2,s3,s4],pet.s1[,s2,s3,s4],dur.s1[,s2,s3,s4])
 out.df <- rbind(out.df,data.frame(M  = method,
                                   Z1 = tau.E[s2],
                                   Z2 = accrual[s3],
                                   Z3 = obdml[s4],
                                   Z4 = 1:ls1,
                                   data.frame(out)))
}}}

write.csv(out.df,"../../data/simulation2_obd_TITE-gBOIN-ET.csv")

out.df <- NULL
for(s2 in 1:ls2){
for(s3 in 1:ls3){
for(s4 in 1:ls4){
for(s1 in 1:ls1){

 out <- rbind(sprop[s1,s2,s3,obdm[s4],],c(NA,apply(data.obs.n[,s1,s2,s3,],2,mean)))
 out.df <- rbind(out.df,data.frame(M  = method,
                                   Z1 = tau.E[s2],
                                   Z2 = accrual[s3],
                                   Z3 = obdml[s4],
                                   Z4 = s1,
                                   Z5 = c("Rec","Pats"),
                                   data.frame(out)))
}}}}

write.csv(out.df,"../../data/simulation2_dose_allocation_TITE-gBOIN-ET.csv")


