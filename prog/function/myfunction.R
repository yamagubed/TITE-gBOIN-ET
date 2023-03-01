
mypinc <- function(pi,phi,phi1,phi2,delta,delta1,lambda1,lambda2,eta1)
{
  nsub <- 100

  l1p  <- pbinom(nsub*lambda1,nsub,phi)
  l1p1 <- pbinom(nsub*lambda1,nsub,phi1)
  l1p2 <- pbinom(nsub*lambda1,nsub,phi2)

  l2p  <- pbinom(nsub*lambda2-1,nsub,phi)
  l2p1 <- pbinom(nsub*lambda2-1,nsub,phi1)
  l2p2 <- pbinom(nsub*lambda2-1,nsub,phi2)

  e1d  <- pbinom(nsub*eta1,nsub,delta)
  e1d1 <- pbinom(nsub*eta1,nsub,delta1)

  pinc <- (pi[1]*(l1p1*(1-e1d1)             + 2/3*(l2p1-l1p1)*e1d1 + (l2p1-l1p1)*(1-e1d1) + (1-l2p1))
         + pi[2]*(l1p1*e1d                  + 2/3*(l2p1-l1p1)*e1d  + (1-l2p1))
         + pi[4]*(l1p*e1d                   + 2/3*(l2p-l1p1)*e1d   + (1-l2p))
         + pi[5]*(l1p2*e1d1 + l1p2*(1-e1d1) + 2/3*(l2p2-l1p2)*e1d1 + (l2p2-l1p2)*(1-e1d1))
         + pi[6]*(l1p2*e1d  + l1p2*(1-e1d)  + 2/3*(l2p2-l1p2)*e1d  + (l2p2-l1p2)*(1-e1d)))

  return(pinc)
}

boundary <- function(phi,delta)
{
  l1  <- round(seq(phi*0.1,phi,by=0.01),digits=2)
  l1s <- length(l1)
  l2  <- round(seq(phi,phi*1.4,by=0.01),digits=2)
  l2s <- length(l2)
  e1  <- round(seq(delta*0.6,delta,by=0.01),digits=2)
  e1s <- length(e1)

  pincval <- array(0,dim=c(l1s,l2s,e1s))

  for(s1 in 1:l1s){
  for(s2 in 1:l2s){
  for(s3 in 1:e1s){

    pincval[s1,s2,s3] <- mypinc(
                          pi      = rep(1/6,6),
                          phi     = phi,
                          phi1    = phi*0.1,
                          phi2    = phi*1.4,
                          delta   = delta,
                          delta1  = delta*0.6,
                          lambda1 = l1[s1],
                          lambda2 = l2[s2],
                          eta1    = e1[s3])

  }}}

  pnum <- which(pincval==min(pincval),arr.ind=TRUE)
  ledf <- data.frame(phi=phi,delta=delta,lambda1=l1[pnum[1]],lambda2=l2[pnum[2]],eta1=e1[pnum[3]])
  return(ledf)
}


futility3 <- function(pe,pt,psi00,psi11)
{
 psi.e0t1 <- 0
 psi.e0t0 <- psi00
 psi.e1t1 <- psi11
 psi.e1t0 <- 100

 ut = (  psi.e0t1*(1-pe)*pt
       + psi.e0t0*(1-pe)*(1-pt)
       + psi.e1t1*pe    *pt
       + psi.e1t0*pe    *(1-pt))

 return(ut)
}

obdselect <- function(pevec,ptvec,ph,ph2,probe,probt)
{
 candose <- which((probe>=0.05)&(probt>=0.05))
 adddose <- which((probe>=0.05)&(probt>=0.05)&(ptvec<=ph2))

 fu31 <- futility3(pe = pevec,
                   pt = ptvec,
                   psi00=40,psi11=60)
 fu31.max <- which(fu31==max(fu31[candose]))

 re311 <- min(intersect(fu31.max,candose))

 if(length(adddose)==0){
  re312 <- 0
 }else{
  re312.max <- which(fu31==max(fu31[adddose]))
  re312 <- min(intersect(re312.max,adddose))
 }

 mdif1 <- min(abs(ptvec[candose]-ph))
 mtd1  <- max(which(abs(ptvec-ph)==mdif1))
 meff1 <- max(pevec[intersect(candose,1:mtd1)])
 deff1 <- which(pevec==meff1)
 re41  <- min(intersect(candose,deff1))

 if(length(adddose)==0){
  re42 <- 0
 }else{
  mdif2 <- min(abs(ptvec[adddose]-ph))
  mtd2  <- max(which(abs(ptvec-ph)==mdif2))
  meff2 <- max(pevec[intersect(adddose,1:mtd2)])
  deff2 <- which(pevec==meff2)
  re42  <- min(intersect(adddose,deff2))
 }

 rval <- c(re311,re312,re41,re42)

 return(rval)
}


