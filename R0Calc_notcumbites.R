#################################################
## Original Model without cumulative bites
## Calculating R0 from NGM
#################################################

library(rSymPy)

b <- Var("b")
b1  <- Var("b1")
S <- Var("S")
H <- Var("H")
tep <- Var("tep")
pep <- Var("pep")
Iep <- Var("Iep")
tpb <- Var("tpb")
ppb <- Var("ppb")
Ipb <- Var("Ipb")
tb <- Var("tb")
pb <- Var("pb")
Ib <- Var("Ib")
U <- Var("U")
alpha <- Var("alpha")
I <- Var("I")
sigma <- Var("sigma")
mu <- Var("mu")
epsilon <- Var("epsilon")
lambdaA <- Var("lambdaA")
lambdaB <- Var("lambdaB")
lambdaC <- Var("lambdaC")
muf <- Var("muf")
tau <- Var("tau")
mupb <- Var("mupb")
mub <- Var("mub")
m <- Var("m") # U/H # vector to host ratio

#save a new variable called Ltrans to the sympy workspace, which is not the main workspace (it's annoying)
sympy("Ltrans = (b*Iep*pep*tep + b*Ipb*ppb*tpb + b1*Ib*pb*tb)")

sympy("T_L_Iep = diff(Ltrans,Iep,1)")
sympy("T_L_Ipb = diff(Ltrans,Ipb,1)")
sympy("T_L_Ib = diff(Ltrans,Ib,1)")
sympy("T_Iep_I = diff(I*m*b*alpha,I,1)")

# transmission matrix, T
cat(sympy("T = Matrix([[0,0,T_L_Iep,T_L_Ipb,T_L_Ib], 
          [0,0,0,0,0],
          [0,T_Iep_I,0,0,0],
          [0,0,0,0,0],
          [0,0,0,0,0]])"),"\n") # each row of a matrix is a different List



sympy("E_L_L = -(sigma+mu)")
sympy("E_I_L = sigma")
sympy("E_I_I = -(mu+epsilon)")
sympy("E_Iep_Iep = -(lambdaA+lambdaB+lambdaC+muf)")
sympy("E_Ipb_Iep = lambdaA")
sympy("E_Ipb_Ipb = -(tau+mupb)")
sympy("E_Ib_Ipb = tau")
sympy("E_Ib_Ib = -mub")

## transition matrix E (Sigma)
cat(sympy("E = Matrix([[E_L_L,0,0,0,0],
      [E_I_L,E_I_I,0,0,0],
      [0,0,E_Iep_Iep,0,0],
      [0,0,E_Ipb_Iep,E_Ipb_Ipb,0],
      [0,0,0,E_Ib_Ipb,E_Ib_Ib]])"),"\n")


Sym("NegInvE = -E.inv()") # take the negative inverse of E


sympy("NGM = T*NegInvE")
Sym("R0 = NGM.eigenvals()")
# [1] "{-(-alpha*m*pep*sigma*tep*b**2/((-epsilon - mu)*(-mu - sigma)*(-lambdaA - lambdaB - lambdaC - muf)) + alpha*lambdaA*m*ppb*sigma*tpb*b**2/((-epsilon - mu)*(-mu - sigma)*(-mupb - tau)*(-lambdaA - lambdaB - lambdaC - muf)) + alpha*b*b1*lambdaA*m*pb*sigma*tau*tb/(mub*(-epsilon - mu)*(-mu - sigma)*(-mupb - tau)*(-lambdaA - lambdaB - lambdaC - muf)))**(1/2): 1, 0: 3, (-alpha*m*pep*sigma*tep*b**2/((-epsilon - mu)*(-mu - sigma)*(-lambdaA - lambdaB - lambdaC - muf)) + alpha*lambdaA*m*ppb*sigma*tpb*b**2/((-epsilon - mu)*(-mu - sigma)*(-mupb - tau)*(-lambdaA - lambdaB - lambdaC - muf)) + alpha*b*b1*lambdaA*m*pb*sigma*tau*tb/(mub*(-epsilon - mu)*(-mu - sigma)*(-mupb - tau)*(-lambdaA - lambdaB - lambdaC - muf)))**(1/2): 1}"

## tried both the first expression and the second (non-zero) expression, and the first one gives negative number, the second is the same positive number. So going with second.

# can't get the below expression to work unless restart session, because it wants it to be a sympy object
# .rs.restartR()




R0.origmodel.function <- function(params, m ){
  R0 <- with(as.list(params),
             (-alpha*m*pep*sigma*tep*b^2/((-epsilon - mu)*(-mu - sigma)*(-lambdaA - lambdaB - lambdaC - muf)) + alpha*lambdaA*m*ppb*sigma*tpb*b^2/((-epsilon - mu)*(-mu - sigma)*(-mupb - tau)*(-lambdaA - lambdaB - lambdaC - muf)) + alpha*b*b1*lambdaA*m*pb*sigma*tau*tb/(mub*(-epsilon - mu)*(-mu - sigma)*(-mupb - tau)*(-lambdaA - lambdaB - lambdaC - muf)))^(1/2)
             
             )
  return(R0)
}

# what m (vector to host ratio) makes R0=1?
# sympy can solve for a variable, but have to make expression equal to zero, so if I want R0=1, then expression I give is R0-1. Then solve for m... so below is solve(R0-1,m)
sympy("solve((-alpha*m*pep*sigma*tep*b**2/((-epsilon - mu)*(-mu - sigma)*(-lambdaA - lambdaB - lambdaC - muf)) + alpha*lambdaA*m*ppb*sigma*tpb*b**2/((-epsilon - mu)*(-mu - sigma)*(-mupb - tau)*(-lambdaA - lambdaB - lambdaC - muf)) + alpha*b*b1*lambdaA*m*pb*sigma*tau*tb/(mub*(-epsilon - mu)*(-mu - sigma)*(-mupb - tau)*(-lambdaA - lambdaB - lambdaC - muf)))**(1/2) - 1,m)")
# ugh it returns no solution




source("RunSIRmodel_cumulative bites.R")

R0.origmodel.function(params.rat_100CFU, m=10)
# 1.975724
# at first glance, it matches up to dynamics. Need to now play with vector to host ratio to manipulate R0 and see if the model dynamics concur. 



###############################################################################
## Vector to Host ratio and R0 for different blood/scenarios

########################## Rat blood

mseq <- seq(0,15,by=0.1)
R0seq_rat100 <- R0.origmodel.function(params.rat_100CFU, m=mseq)
plot(mseq,R0seq_rat100,xlab="vector to host ratio",ylab="Ro",main="rat blood",type="l",lwd=2)
R0seq_rat10 <- R0.origmodel.function(params.rat_10CFU, m=mseq)
lines(mseq,R0seq_rat10,lty=2)
R0seq_rat1 <- R0.origmodel.function(params.rat_1CFU, m=mseq)
lines(mseq,R0seq_rat1,lty=4)
abline(h=1,lty=3)
legend("bottomright",c("100 CFU","10 CFU","1 CFU"),lty=c(1,2,4),bty="n")

# some disagreement on whether should square what comes from the NGM - either way, where R0=1 is the same (1 squared is 1)
plot(mseq,R0seq_rat100^2,xlab="vector to host ratio",ylab="Ro",main="rat blood",type="l",lwd=2)
lines(mseq,R0seq_rat10^2,lty=2)
lines(mseq,R0seq_rat1^2,lty=4)
abline(h=1,lty=3)
legend("bottomright",c("100 CFU","10 CFU","1 CFU"),lwd=c(2,1,1),lty=c(1,2,4),bty="n")

# rat blood 10 CFU with and without early phase - assuming not squared
plot(mseq,R0seq_rat10,xlab="vector to host ratio",ylab="Ro",main="rat blood 10 CFU",type="l",lwd=2)
R0seq_rat10_BPBonly <- R0.origmodel.function(params.rat_10CFU.BPBonly, m=mseq)
lines(mseq,R0seq_rat10_BPBonly,col="red")
R0seq_rat10_EPTonly <- R0.origmodel.function(params.rat_10CFU.EPTonly, m=mseq)
lines(mseq,R0seq_rat10_EPTonly,col="green")
abline(h=1,lty=3)
legend("topleft",c("All","BPB only","EPT only"),lwd=c(2,1,1),col=c("black","red","green"),bty="n")

# rat blood 10 CFU with and without early phase - assuming squared
plot(mseq,R0seq_rat10^2,xlab="vector to host ratio",ylab="Ro",main="rat blood 10 CFU",type="l",lwd=2)
lines(mseq,R0seq_rat10_BPBonly^2,col="red")
lines(mseq,R0seq_rat10_EPTonly^2,col="green")
abline(h=1,lty=3)
legend("topleft",c("All","BPB only","EPT only"),lwd=c(2,1,1),col=c("black","red","green"),bty="n")



#using squared
R0_squared_RatBlood <- data.frame(VtoH_ratio=mseq,CFU100=R0seq_rat100^2,CFU10=R0seq_rat10^2,CFU1=R0seq_rat1^2,BPBonly_CFU10=R0seq_rat10_BPBonly^2,EPTonly_CFU10=R0seq_rat10_EPTonly^2)

# not squared
R0_RatBlood <- data.frame(VtoH_ratio=mseq,CFU100=R0seq_rat100,CFU10=R0seq_rat10,CFU1=R0seq_rat1,BPBonly_CFU10=R0seq_rat10_BPBonly,EPTonly_CFU10=R0seq_rat10_EPTonly)


########################## Mouse blood

#using squared
mseq <- seq(0,15,by=0.1)
R0seq_mouse100 <- R0.origmodel.function(params.mouse_100CFU, m=mseq)
plot(mseq,R0seq_mouse100^2,xlab="vector to host ratio",ylab="Ro",main="mouse blood",type="l",lwd=2)
R0seq_mouse10 <- R0.origmodel.function(params.mouse_10CFU, m=mseq)
lines(mseq,R0seq_mouse10^2,lty=2)
R0seq_mouse1 <- R0.origmodel.function(params.mouse_1CFU, m=mseq)
lines(mseq,R0seq_mouse1^2,lty=4)
abline(h=1,lty=3)
legend("bottomright",c("100 CFU","10 CFU","1 CFU"),lwd=c(2,1,1),lty=c(1,2,4),bty="n")

# not squared
plot(mseq,R0seq_mouse100,xlab="vector to host ratio",ylab="Ro",main="mouse blood",type="l",lwd=2)
lines(mseq,R0seq_mouse10,lty=2)
lines(mseq,R0seq_mouse1,lty=4)
abline(h=1,lty=3)
legend("bottomright",c("100 CFU","10 CFU","1 CFU"),lwd=c(2,1,1),lty=c(1,2,4),bty="n")



# mouse blood 10 CFU with and without early phase - assuming squared
plot(mseq,R0seq_mouse10^2,xlab="vector to host ratio",ylab="Ro",main="mouse blood 10 CFU",type="l",lwd=2)
R0seq_mouse10_BPBonly <- R0.origmodel.function(params.mouse_10CFU.BPBonly, m=mseq)
lines(mseq,R0seq_mouse10_BPBonly^2,col="red")
R0seq_mouse10_EPTonly <- R0.origmodel.function(params.mouse_10CFU.EPTonly, m=mseq)
lines(mseq,R0seq_mouse10_EPTonly^2,col="green")
abline(h=1,lty=3)
legend("topleft",c("All","BPB only","EPT only"),lwd=c(2,1,1),col=c("black","red","green"),bty="n")

# not squared
plot(mseq,R0seq_mouse10,xlab="vector to host ratio",ylab="Ro",main="mouse blood 10 CFU",type="l",lwd=2)
lines(mseq,R0seq_mouse10_BPBonly,col="red")
lines(mseq,R0seq_mouse10_EPTonly,col="green")
abline(h=1,lty=3)
legend("topleft",c("All","BPB only","EPT only"),lwd=c(2,1,1),col=c("black","red","green"),bty="n")



#using squared
R0_squared_mouseBlood <- data.frame(VtoH_ratio=mseq,CFU100=R0seq_mouse100^2,CFU10=R0seq_mouse10^2,CFU1=R0seq_mouse1^2,BPBonly_CFU10=R0seq_mouse10_BPBonly^2,EPTonly_CFU10=R0seq_mouse10_EPTonly^2)

#not squared
R0_mouseBlood <- data.frame(VtoH_ratio=mseq,CFU100=R0seq_mouse100,CFU10=R0seq_mouse10,CFU1=R0seq_mouse1,BPBonly_CFU10=R0seq_mouse10_BPBonly,EPTonly_CFU10=R0seq_mouse10_EPTonly)











########################################################################



SIR.rat.m10.100CFU <- print.plague.one.SIRmodel(
  model = SIR.model.fleasperhost,
  params = params.rat_100CFU, 
  xstart=c(S=9,L=0,I=1,E=0,R=0,dr=0,Id=0, U=100,Iep=0,Ipb=0,Ib=0,df=0,cumb=0),
  T=100,
  title="Rat 100CFU, V/H=10,")
# infected dead 
SIR.rat.m10.100CFU$Id[dim(SIR.rat.m10.100CFU)[1]]
# recovered
SIR.rat.m10.100CFU$R[dim(SIR.rat.m10.100CFU)[1]]

# similar R0:
SIR.rat.m4.7.10CFU <- print.plague.one.SIRmodel(
  model = SIR.model.fleasperhost,
  params = params.rat_100CFU, 
  xstart=c(S=9,L=0,I=1,E=0,R=0,dr=0,Id=0, U=47,Iep=0,Ipb=0,Ib=0,df=0,cumb=0),
  T=100,
  title="Rat 10CFU, V/H=4.7,")
# infected dead 
SIR.rat.m4.7.10CFU$Id[dim(SIR.rat.m4.7.10CFU)[1]]
# recovered 
SIR.rat.m4.7.10CFU$R[dim(SIR.rat.m4.7.10CFU)[1]]




######################### Table of different rat scenarios with vector to host ratio of 5:

Rat.Table <- data.frame(scenario=c("1 CFU, m=5","10 CFU, m=5","100 CFU, m=5","10 CFU BPBonly, m=5","10 CFU EPT only, m=5"),R0=rep(NA,5),Idead=rep(NA,5),Recovered=rep(NA,5),cumblocked=rep(NA,5))


paramslist<-list(params.rat_1CFU,params.rat_10CFU, params.rat_100CFU, params.rat_10CFU.BPBonly,params.rat_10CFU.EPTonly)

for(i in 1:length(paramslist)){
  out <- print.plague.one.SIRmodel(
    model = SIR.model.fleasperhost,
    params = paramslist[[i]], 
    xstart=c(S=9,L=0,I=1,E=0,R=0,dr=0,Id=0, U=50,Iep=0,Ipb=0,Ib=0,df=0,cumb=0),
    T=100,
    title=i)
  Rat.Table$R0[i] =  R0.origmodel.function(paramslist[[i]], m=5)
  Rat.Table$Idead[i] = out$Id[dim(out)[1]]
  Rat.Table$Recovered[i] = out$R[dim(out)[1]]
  Rat.Table$cumblocked[i] = out$cumb[dim(out)[1]]
}



# the last model should have R0<1, but still getting some transmission. About 2 in R and 2 I-dead by end

# how could I approximate R0 with simulatiins? How many fleas does 1 I infect over the average infectious period? Mult by how many S's are infected from 1 infected flea. But there are 3 categories of infected fleas, that have different potentials. Calc the average infectious period of each?

# Conceptually, R0 is the # of I's (vert host) that each I infects (over 1 generation time). IT should be equal to the number of fleas the I infects, times the number of vertebrates all those fleas infect and make it to the I class.

# R0 = X*Y*Z

# Where:
# X = number of fleas infected from 1 I
X = alpha*b*V/H

######################### Table of different mouse scenarios with vector to host ratio of 5:

Mouse.Table <- data.frame(scenario=c("1 CFU, m=5","10 CFU, m=5","100 CFU, m=5","10 CFU BPBonly, m=5","10 CFU EPT only, m=5"),R0=rep(NA,5),Idead=rep(NA,5),Recovered=rep(NA,5),cumblocked=rep(NA,5))


paramslist<-list(params.mouse_1CFU,params.mouse_10CFU, params.mouse_100CFU, params.mouse_10CFU.BPBonly,params.mouse_10CFU.EPTonly)

for(i in 1:length(paramslist)){
  out <- print.plague.one.SIRmodel(
    model = SIR.model.fleasperhost,
    params = paramslist[[i]], 
    xstart=c(S=9,L=0,I=1,E=0,R=0,dr=0,Id=0, U=50,Iep=0,Ipb=0,Ib=0,df=0,cumb=0),
    T=100,
    title=i)
  Mouse.Table$R0[i] =  R0.origmodel.function(paramslist[[i]], m=5)
  Mouse.Table$Idead[i] = out$Id[dim(out)[1]]
  Mouse.Table$Recovered[i] = out$R[dim(out)[1]]
  Mouse.Table$cumblocked[i] = out$cumb[dim(out)[1]]
}



# Y = # S's infected to L state from each infectious flea
# this is tricky because depends on the probability a flea is each infectious category

# Z = prob that each of those L survives their latent state to I
Z = sigma /(sigma+mu)


### attempt to define Y:

# bites times probability infectious for each category
prob.inf_ep = b*pep*tep
prob.inf_pb = b*ppb*tpb
prob.inf_b = b1*pb*tb

#infectious period for each category
inf.per_ep = 1/(lambdaA + lambdaB + lambdaC + muf)
inf.per_pb = 1/(tau + mupb)
inf.per_b = 1/mub

# prob survive through early phase and go to partially blocked
surv_ep = lambdaA/(lambdaA + lambdaB + lambdaC + muf)

# prob survive through both early phase and partially blocked phase
surv_pb = surv_ep * tau/(tau + mupb)

# total infectious period # not sure about this
inf.per_ep + inf.per_pb * surv_ep + inf.per_b * surv_pb

# what proportion of the time are they in each of the different periods?



