

#### Used rSymPy packages to get expression for R0 using the Next Generation Matrix
#### the resulting equation is in the function below


R0.origmodel.function <- function(params, m ){
  R0 <- with(as.list(params),
             (-alpha*m*pep*sigma*tep*b^2/((-epsilon - mu)*(-mu - sigma)*(-lambdaA - lambdaB - lambdaC - muf)) + alpha*lambdaA*m*ppb*sigma*tpb*b^2/((-epsilon - mu)*(-mu - sigma)*(-mupb - tau)*(-lambdaA - lambdaB - lambdaC - muf)) + alpha*b*b1*lambdaA*m*pb*sigma*tau*tb/(mub*(-epsilon - mu)*(-mu - sigma)*(-mupb - tau)*(-lambdaA - lambdaB - lambdaC - muf)))^(1/2)
             
  )
  return(R0)
}


# get the rat and mouse parameters from this R file
source("RunSIRmodel_fleasperhost.R")


###############################################################################
## Plots of R0 against Vector to Host ratio  for different blood/scenarios

# There is some disagreement in the literature on whether you should square 
# what comes from the NGM - either way, where R0=1 is the same (1 squared is 1)
# all those below assume not squared

library(viridis) # color scales that are supposed to be good for colorblindness and gray scale
#CFUcolors <- viridis(3) # color scale to represent different CFUS for model with all phases
modelcolors <- viridis(3)  # mako(5)[2:4] # color scale to represent different models for a certain CFU



#################################################### Rat blood
####################################################

############# plots of how vector to host ratio affects R0 for the 3 CFUs
mseq <- seq(0,20,by=0.1)
R0seq_rat100 <- R0.origmodel.function(params.rat_100CFU, m=mseq)
R0seq_rat10 <- R0.origmodel.function(params.rat_10CFU, m=mseq)
R0seq_rat1 <- R0.origmodel.function(params.rat_1CFU, m=mseq)

pdf("Ro_VtoH_Rat_byCFU.pdf", onefile = TRUE, paper = "special", height = 4, width = 4)
  plot(mseq,R0seq_rat100,xlab="vector to host ratio",ylab="Ro",main="rat blood",lwd=2,type="l",ylim=c(0,2.5), col="gray80")
  lines(mseq,R0seq_rat10,lwd=2,col="gray50")
  lines(mseq,R0seq_rat1,lwd=2,col="gray0")
  abline(h=1,lty=3)
  legend("bottomright",c("100 CFU","10 CFU","1 CFU"),col=c("gray80","gray50","gray0"),lwd=2,bty="n")
dev.off()


################ plots of how vector to host ratio affects R0 for phases
# rat blood 10 CFU with and without early phase  
R0seq_rat10_BPBonly <- R0.origmodel.function(params.rat_10CFU.BPBonly, m=mseq)
R0seq_rat10_EPTonly <- R0.origmodel.function(params.rat_10CFU.EPTonly, m=mseq)

pdf("Ro_VtoH_Rat10CFU_byPhase.pdf", onefile = TRUE, paper = "special", height = 4, width = 4)
  plot(mseq,R0seq_rat10,xlab="vector to host ratio",ylab="Ro",main="rat blood 10 CFU",col=modelcolors[1],lwd=2,type="l",ylim=c(0,2.5))
  lines(mseq,R0seq_rat10_BPBonly,col=modelcolors[2],lwd=2)
  lines(mseq,R0seq_rat10_EPTonly,col=modelcolors[3],lwd=2)
  abline(h=1,lty=3)
  legend("bottomright",c("All","BPB only","EPT only"),col=modelcolors,lwd=2,bty="n")
dev.off()

# rat blood 1 CFU with and without early phase  
R0seq_rat1_BPBonly <- R0.origmodel.function(params.rat_1CFU.BPBonly, m=mseq)
R0seq_rat1_EPTonly <- R0.origmodel.function(params.rat_1CFU.EPTonly, m=mseq)

pdf("Ro_VtoH_Rat1CFU_byPhase.pdf", onefile = TRUE, paper = "special", height = 4, width = 4)
plot(mseq,R0seq_rat1,xlab="vector to host ratio",ylab="Ro",main="rat blood 1 CFU",col=modelcolors[1],lwd=2,type="l",ylim=c(0,2.5))
lines(mseq,R0seq_rat1_BPBonly,col=modelcolors[2],lwd=2)
lines(mseq,R0seq_rat1_EPTonly,col=modelcolors[3],lwd=2)
abline(h=1,lty=3)
legend("bottomright",c("All","BPB only","EPT only"),col=modelcolors,lwd=2,bty="n")
dev.off()



# rat blood 100 CFU with and without early phase  
R0seq_rat100_BPBonly <- R0.origmodel.function(params.rat_100CFU.BPBonly, m=mseq)
R0seq_rat100_EPTonly <- R0.origmodel.function(params.rat_100CFU.EPTonly, m=mseq)

pdf("Ro_VtoH_Rat100CFU_byPhase.pdf", onefile = TRUE, paper = "special", height = 4, width = 4)
plot(mseq,R0seq_rat100,xlab="vector to host ratio",ylab="Ro",main="rat blood 100 CFU",col=modelcolors[1],lwd=2,type="l",ylim=c(0,2.5))
lines(mseq,R0seq_rat100_BPBonly,col=modelcolors[2],lwd=2)
lines(mseq,R0seq_rat100_EPTonly,col=modelcolors[3],lwd=2)
abline(h=1,lty=3)
legend("bottomright",c("All","BPB only","EPT only"),col=modelcolors,lwd=2,bty="n")
dev.off()

#################################################### Mouse blood
####################################################
 
# plots of how vector to host ratio affects R0 for the 3 CFUs
R0seq_mouse100 <- R0.origmodel.function(params.mouse_100CFU, m=mseq)
R0seq_mouse10 <- R0.origmodel.function(params.mouse_10CFU, m=mseq)
R0seq_mouse1 <- R0.origmodel.function(params.mouse_1CFU, m=mseq)

pdf("Ro_VtoH_Mouse_byCFU.pdf", onefile = TRUE, paper = "special", height = 4, width = 4)
  plot(mseq,R0seq_mouse100,xlab="vector to host ratio",ylab="Ro",main="mouse blood",type="l",lwd=2,ylim=c(0,2.5), col="gray80")
  lines(mseq,R0seq_mouse10,lwd=2, col="gray50")
  lines(mseq,R0seq_mouse1,lwd=2, col="gray0")
  abline(h=1,lty=3)
  legend("bottomright",c("100 CFU","10 CFU","1 CFU"),col=c("gray80","gray50","gray0"),lwd=2,bty="n")
dev.off()


################ plots of how vector to host ratio affects R0 for phases
 
# mouse blood 10 CFU with and without early phase  
R0seq_mouse10_BPBonly <- R0.origmodel.function(params.mouse_10CFU.BPBonly, m=mseq)
R0seq_mouse10_EPTonly <- R0.origmodel.function(params.mouse_10CFU.EPTonly, m=mseq)

pdf("Ro_VtoH_Mouse10CFU_byPhase.pdf", onefile = TRUE, paper = "special", height = 4, width = 4)
  plot(mseq,R0seq_mouse10,xlab="vector to host ratio",ylab="Ro",main="mouse blood 10 CFU",col=modelcolors[1],lwd=2,type="l",ylim=c(0,2.5))
  lines(mseq,R0seq_mouse10_BPBonly,col=modelcolors[2],lwd=2)
  lines(mseq,R0seq_mouse10_EPTonly,col=modelcolors[3],lwd=2)
  abline(h=1,lty=3)
  legend("bottomright",c("All","BPB only","EPT only"),lwd=2,col=modelcolors,bty="n")
dev.off()
  

# mouse blood 1 CFU with and without early phase  
R0seq_mouse1_BPBonly <- R0.origmodel.function(params.mouse_1CFU.BPBonly, m=mseq)
R0seq_mouse1_EPTonly <- R0.origmodel.function(params.mouse_1CFU.EPTonly, m=mseq)

pdf("Ro_VtoH_Mouse1CFU_byPhase.pdf", onefile = TRUE, paper = "special", height = 4, width = 4)
plot(mseq,R0seq_mouse1,xlab="vector to host ratio",ylab="Ro",main="mouse blood 1 CFU",col=modelcolors[1],lwd=2,type="l",ylim=c(0,2.5))
lines(mseq,R0seq_mouse1_BPBonly,col=modelcolors[2],lwd=2)
lines(mseq,R0seq_mouse1_EPTonly,col=modelcolors[3],lwd=2)
abline(h=1,lty=3)
legend("bottomright",c("All","BPB only","EPT only"),lwd=2,col=modelcolors,bty="n")
dev.off()


# mouse blood 100 CFU with and without early phase  
R0seq_mouse100_BPBonly <- R0.origmodel.function(params.mouse_100CFU.BPBonly, m=mseq)
R0seq_mouse100_EPTonly <- R0.origmodel.function(params.mouse_100CFU.EPTonly, m=mseq)

pdf("Ro_VtoH_Mouse100CFU_byPhase.pdf", onefile = TRUE, paper = "special", height = 4, width = 4)
plot(mseq,R0seq_mouse100,xlab="vector to host ratio",ylab="Ro",main="mouse blood 100 CFU",col=modelcolors[1],lwd=2,type="l",ylim=c(0,2.5))
lines(mseq,R0seq_mouse100_BPBonly,col=modelcolors[2],lwd=2)
lines(mseq,R0seq_mouse100_EPTonly,col=modelcolors[3],lwd=2)
abline(h=1,lty=3)
legend("bottomright",c("All","BPB only","EPT only"),lwd=2,col=modelcolors,bty="n")
dev.off()


############################################################
## PLots of the dynamics
############################################################



######## Rat


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


