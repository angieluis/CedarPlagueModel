
## Define parameter  values
library(deSolve)
source("SIRmodel_fleasperhost.R")

###########################
### All classes
###########################

params.mouse_1CFU <- c(alpha=1,lambdaA=0.035,lambdaB=0.20,lambdaC=0.07,b=0.4,b1=2,tau=0.39, muf=0.02, mupb=0.14,mub=0.20,
                       pep=0.03,ppb=0.11,pb=0.5, tep=1, tpb=1, tb=1, 
                       mu=0.002,sigma=0.25,gamma=0.14,epsilon=0.5)
params.rat_1CFU <- c(alpha=1,lambdaA=0.04,lambdaB=0.02,lambdaC=0.06,b=0.4,b1=2,tau=0.48, muf=0.02, mupb=0.13,mub=0.26,
                     pep=0.1,ppb=0.10,pb=0.67,tep=1, tpb=1 , tb=1 ,
                     mu=0.002,sigma=0.25,gamma=0.14,epsilon=0.5)	
params.mouse_10CFU <- c(alpha=1,lambdaA=0.035,lambdaB=0.20,lambdaC=0.07,b=0.4,b1=2,tau=0.39, muf=0.02, mupb=0.14,mub=0.20,
                        pep=0.03,ppb=0.11,pb=0.5, tep=0.001, tpb=0.5, tb=0.8, 
                        mu=0.002,sigma=0.25,gamma=0.14,epsilon=0.5)	
params.rat_10CFU <- c(alpha=1,lambdaA=0.04,lambdaB=0.02,lambdaC=0.06,b=0.4,b1=2,tau=0.48, muf=0.02, mupb=0.13,mub=0.26,
                      pep=0.1,ppb=0.10,pb=0.67,tep=0.5, tpb=1 , tb=0.8 ,
                      mu=0.002,sigma=0.25,gamma=0.14,epsilon=0.5)	
params.mouse_100CFU <- c(alpha=1,lambdaA=0.035,lambdaB=0.20,lambdaC=0.07,b=0.4,b1=2,tau=0.39, muf=0.02, mupb=0.14,mub=0.20,
                         pep=0.03,ppb=0.11,pb=0.5, tep=0.00, tpb=0, tb=0.65, 
                         mu=0.002,sigma=0.25,gamma=0.14,epsilon=0.5)	
params.rat_100CFU <- c(alpha=1,lambdaA=0.04,lambdaB=0.02,lambdaC=0.06,b=0.4,b1=2,tau=0.48, muf=0.02, mupb=0.13,mub=0.26,
                       pep=0.1,ppb=0.10,pb=0.67,tep=0, tpb=1 , tb=0.41 ,
                       mu=0.002,sigma=0.25,gamma=0.14,epsilon=0.5)		


###############################################################################
### Now with Early Phase Only
###############################################################################


##params
params.mouse_1CFU.EPTonly <-  params.mouse_1CFU
  params.mouse_1CFU.EPTonly[c("ppb","pb","tpb","tb")]<-0
params.rat_1CFU.EPTonly <- params.rat_1CFU
  params.rat_1CFU.EPTonly[c("ppb","pb","tpb","tb")] <- 0
params.mouse_10CFU.EPTonly <- params.mouse_10CFU 
  params.mouse_10CFU.EPTonly[c("ppb","pb","tpb","tb")] <- 0
params.rat_10CFU.EPTonly <- params.rat_10CFU 
  params.rat_10CFU.EPTonly[c("ppb","pb","tpb","tb")] <- 0
params.mouse_100CFU.EPTonly <- params.mouse_100CFU 
  params.mouse_100CFU.EPTonly[c("ppb","pb","tpb","tb")] <- 0
params.rat_100CFU.EPTonly <- params.rat_100CFU 
  params.rat_100CFU.EPTonly[c("ppb","pb","tpb","tb")] <- 0
  

###############################################################################
### Now with no early phase, only partially and fully blocked
###############################################################################
  
 
##params
params.mouse_1CFU.BPBonly <-  params.mouse_1CFU
  params.mouse_1CFU.BPBonly[c("pep","tep")]<-0
params.rat_1CFU.BPBonly <- params.rat_1CFU
  params.rat_1CFU.BPBonly[c("pep","tep")] <- 0
params.mouse_10CFU.BPBonly <- params.mouse_10CFU 
  params.mouse_10CFU.BPBonly[c("pep","tep")] <- 0
params.rat_10CFU.BPBonly <- params.rat_10CFU 
  params.rat_10CFU.BPBonly[c("pep","tep")] <- 0
params.mouse_100CFU.BPBonly <- params.mouse_100CFU 
  params.mouse_100CFU.BPBonly[c("pep","tep")] <- 0
params.rat_100CFU.BPBonly <- params.rat_100CFU 
  params.rat_100CFU.BPBonly[c("pep","tep")] <- 0
  

##################################################################
# using comparison function



params.mouse <- list(params.mouse_1CFU, params.mouse_1CFU.EPTonly, params.mouse_1CFU.BPBonly, 
                     params.mouse_10CFU, params.mouse_10CFU.EPTonly, params.mouse_10CFU.BPBonly,
                     params.mouse_100CFU, params.mouse_100CFU.EPTonly, params.mouse_100CFU.BPBonly)
names(params.mouse) <- c("mouse_1CFU","mouse_1CFU.EPTonly", "mouse_1CFU.BPBonly", 
                         "mouse_10CFU", "mouse_10CFU.EPTonly", "mouse_10CFU.BPBonly",
                         "mouse_100CFU", "mouse_100CFU.EPTonly", "mouse_100CFU.BPBonly")

SIRcum.mouse.comparisons <- print.plague.SIRmodel(
  model = SIR.model.fleasperhost,
  params = params.mouse, 
  xstart=c(S=9,L=0,I=1,E=0,R=0,dr=0,Id=0, 
           U=50,Iep=0,Ipb=0,Ib=0,df=0,cumb=0),
  T=100,
  plot.name="Plague SIR, m=5, mouse comparison.pdf")




params.rat <- list(params.rat_1CFU, params.rat_1CFU.EPTonly, params.rat_1CFU.BPBonly, 
                   params.rat_10CFU, params.rat_10CFU.EPTonly, params.rat_10CFU.BPBonly,
                   params.rat_100CFU, params.rat_100CFU.EPTonly, params.rat_100CFU.BPBonly)
names(params.rat) <- c("rat_1CFU","rat_1CFU.EPTonly", "rat_1CFU.BPBonly", 
                       "rat_10CFU", "rat_10CFU.EPTonly", "rat_10CFU.BPBonly",
                       "rat_100CFU", "rat_100CFU.EPTonly", "rat_100CFU.BPBonly")

SIR.rat.comparisons <- print.plague.SIRmodel(
  model = SIR.model.fleasperhost,
  params = params.rat, 
  xstart=c(S=9,L=0,I=1,E=0,R=0,dr=0,Id=0, 
           U=50,Iep=0,Ipb=0,Ib=0,df=0,cumb=0),
  T=100,
  plot.name="Plague SIR, m=5, rat comparison.pdf")
