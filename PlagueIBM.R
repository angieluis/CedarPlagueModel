####################################################################
## Code to Run individual based model & plot
####################################################################

#### think about params now 
#### In this model, they are probabilties not rates
#### Some of the params are already probabilties/proportions (t's and p's)
#### And the ones I made a Poisson variable are ok? (biting rates)
#### I converted some of the rates to probabilities via prob= 1-exp(-rate*t)
#### including death rates, rates of transition (mu, epsilon, lambdas, tau, gamma, etc)

#### Could get more detailed parameter distributions from Joe to see how 
#### the variability matters

#### Think about whether I'm modeling t's and p's right. 
#### Prob of transmission is per bite and all or nothing
#### Proportion leaving to L & I is cumulative - I added t's like they were doses,
#### but they are defined as probabilities, so I should probably do 1-product(1-t)


source("IndividualBasedModel.R")

xstart=c(S=9,L=0,I=1,E=0,R=0,dr=0,Id=0, U=100,Iep=0,Ipb=0,Ib=0,df=0)					#beginning population size
T = 100 # total number of days to simulate
n.sim = 100


###############################################################################
### Baseline scenarios with all the flea compartments
###############################################################################


##params
params.mouse_1CFU <- c(alpha=1,lambdaA=1-exp(-0.035),lambdaB=1-exp(-0.20),lambdaC=1-exp(-0.07),
                       b=0.4,b1=2,tau=1-exp(-0.39), muf=1-exp(-0.02), mupb=1-exp(-0.14),mub=1-exp(-0.20), 
                       pep=0.03,ppb=0.11,pb=0.5, tep=1, tpb=1, tb=1, mu=1-exp(-0.002),
                       sigma=1-exp(-0.25),gamma=1-exp(-0.14),epsilon=1-exp(-0.5))
params.rat_1CFU <- c(alpha=1,lambdaA=1-exp(-0.04),lambdaB=1-exp(-0.02),lambdaC=1-exp(-0.06),b=0.4,
                     b1=2,tau=1-exp(-0.48), muf=1-exp(-0.02), mupb=1-exp(-0.13),mub=1-exp(-0.26), pep=0.1,
                     ppb=0.10,pb=0.67,tep=1, tpb=1, tb=1, mu=1-exp(-0.002),sigma=1-exp(-0.25),
                     gamma=1-exp(-0.14),epsilon=1-exp(-0.5))	
params.mouse_10CFU <- c(alpha=1,lambdaA=1-exp(-0.035),lambdaB=1-exp(-0.20),lambdaC=1-exp(-0.07),b=0.4,
                        b1=2,tau=1-exp(-0.39), muf=1-exp(-0.02), mupb=1-exp(-0.14),mub=1-exp(-0.20),
                        pep=0.03,ppb=0.11,pb=0.5, tep=0.001, tpb=0.5, tb=0.8, 
                        mu=1-exp(-0.002),sigma=1-exp(-0.25),gamma=1-exp(-0.14),epsilon=1-exp(-0.5))	
params.rat_10CFU <- c(alpha=1,lambdaA=1-exp(-0.04),lambdaB=1-exp(-0.02),lambdaC=1-exp(-0.06),b=0.4,
                      b1=2,tau=1-exp(-0.48), muf=1-exp(-0.02), mupb=1-exp(-0.13),mub=1-exp(-0.26),
                      pep=0.1,ppb=0.10,pb=0.67,tep=0.5, tpb=1 , tb=0.8 ,
                      mu=1-exp(-0.002),sigma=1-exp(-0.25),gamma=1-exp(-0.14),epsilon=1-exp(-0.5))	
params.mouse_100CFU <- c(alpha=1,lambdaA=1-exp(-0.035),lambdaB=1-exp(-0.20),lambdaC=1-exp(-0.07),b=0.4,
                        b1=2,tau=1-exp(-0.39), muf=1-exp(-0.02), mupb=1-exp(-0.14),mub=1-exp(-0.20),
                        pep=0.03,ppb=0.11,pb=0.5, tep=0.00, tpb=0, tb=0.65, 
                        mu=1-exp(-0.002),sigma=1-exp(-0.25),gamma=1-exp(-0.14),epsilon=1-exp(-0.5))	
params.rat_100CFU <- c(alpha=1,lambdaA=1-exp(-0.04),lambdaB=1-exp(-0.02),lambdaC=1-exp(-0.06),b=0.4,
                      b1=2,tau=1-exp(-0.48), muf=1-exp(-0.02), mupb=1-exp(-0.13),mub=1-exp(-0.26),
                      pep=0.1,ppb=0.10,pb=0.67,tep=0, tpb=1 , tb=0.41 ,
                      mu=1-exp(-0.002),sigma=1-exp(-0.25),gamma=1-exp(-0.14),epsilon=1-exp(-0.5))	









######################################################################################################

output.mouse_1CFU <- plague.IBM(params.mouse_1CFU, # vector of params with names, see below
                  xstart, # vector of starting values with names
                  T,      # number of time steps to simulate
                  n.sim  # number of simulations
)

rodent.means.mouse_1CFU <- as.data.frame(output.mouse_1CFU$mean.rodent.ts)
plot.ts(rodent.means.mouse_1CFU$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(rodent.means.mouse_1CFU$E,col="blue",lwd=2,lty=1)
lines(rodent.means.mouse_1CFU$S,col="forestgreen",lwd=2,lty=1)
lines(rodent.means.mouse_1CFU$L,col="red",lwd=2,lty=5)
lines(rodent.means.mouse_1CFU$R,col="blue",lwd=2,lty=5)
lines(rodent.means.mouse_1CFU$Id,col="black",lwd=2,lty=1)
legend("topright",c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp =0.5, y.intersp =1)
title(main="Mouse Blood, 1 CFU")

#############
output.rat_1CFU <- plague.IBM(params.rat_1CFU, # vector of params with names, see below
                                xstart, # vector of starting values with names
                                T,      # number of time steps to simulate
                                n.sim  # number of simulations
)

rodent.means.rat_1CFU <- as.data.frame(output.rat_1CFU$mean.rodent.ts)
plot.ts(rodent.means.rat_1CFU$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(rodent.means.rat_1CFU$E,col="blue",lwd=2,lty=1)
lines(rodent.means.rat_1CFU$S,col="forestgreen",lwd=2,lty=1)
lines(rodent.means.rat_1CFU$L,col="red",lwd=2,lty=5)
lines(rodent.means.rat_1CFU$R,col="blue",lwd=2,lty=5)
lines(rodent.means.rat_1CFU$Id,col="black",lwd=2,lty=1)
legend("topright",c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp =0.5, y.intersp =1)
title(main="Rat Blood, 1 CFU")


#############
output.mouse_10CFU <- plague.IBM(params.mouse_10CFU, # vector of params with names, see below
                                xstart, # vector of starting values with names
                                T,      # number of time steps to simulate
                                n.sim  # number of simulations
)

rodent.means.mouse_10CFU <- as.data.frame(output.mouse_10CFU$mean.rodent.ts)
plot.ts(rodent.means.mouse_10CFU$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(rodent.means.mouse_10CFU$E,col="blue",lwd=2,lty=1)
lines(rodent.means.mouse_10CFU$S,col="forestgreen",lwd=2,lty=1)
lines(rodent.means.mouse_10CFU$L,col="red",lwd=2,lty=5)
lines(rodent.means.mouse_10CFU$R,col="blue",lwd=2,lty=5)
lines(rodent.means.mouse_10CFU$Id,col="black",lwd=2,lty=1)
legend("topright",c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp =0.5, y.intersp =1)
title(main="Mouse Blood, 10 CFU")

###############
output.rat_10CFU <- plague.IBM(params.rat_10CFU, # vector of params with names, see below
                              xstart, # vector of starting values with names
                              T,      # number of time steps to simulate
                              n.sim  # number of simulations
)

rodent.means.rat_10CFU <- as.data.frame(output.rat_10CFU$mean.rodent.ts)
plot.ts(rodent.means.rat_10CFU$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(rodent.means.rat_10CFU$E,col="blue",lwd=2,lty=1)
lines(rodent.means.rat_10CFU$S,col="forestgreen",lwd=2,lty=1)
lines(rodent.means.rat_10CFU$L,col="red",lwd=2,lty=5)
lines(rodent.means.rat_10CFU$R,col="blue",lwd=2,lty=5)
lines(rodent.means.rat_10CFU$Id,col="black",lwd=2,lty=1)
legend("topright",c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp =0.5, y.intersp =1)
title(main="Rat Blood, 10 CFU")

#############
output.mouse_100CFU <- plague.IBM(params.mouse_100CFU, # vector of params with names, see below
                                 xstart, # vector of starting values with names
                                 T,      # number of time steps to simulate
                                 n.sim  # number of simulations
)

rodent.means.mouse_100CFU <- as.data.frame(output.mouse_100CFU$mean.rodent.ts)
plot.ts(rodent.means.mouse_100CFU$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(rodent.means.mouse_100CFU$E,col="blue",lwd=2,lty=1)
lines(rodent.means.mouse_100CFU$S,col="forestgreen",lwd=2,lty=1)
lines(rodent.means.mouse_100CFU$L,col="red",lwd=2,lty=5)
lines(rodent.means.mouse_100CFU$R,col="blue",lwd=2,lty=5)
lines(rodent.means.mouse_100CFU$Id,col="black",lwd=2,lty=1)
legend("topright",c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp =0.5, y.intersp =1)
title(main="Mouse Blood, 100 CFU")

###############
output.rat_100CFU <- plague.IBM(params.rat_100CFU, # vector of params with names, see below
                               xstart, # vector of starting values with names
                               T,      # number of time steps to simulate
                               n.sim  # number of simulations
)

rodent.means.rat_100CFU <- as.data.frame(output.rat_100CFU$mean.rodent.ts)
plot.ts(rodent.means.rat_100CFU$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(rodent.means.rat_100CFU$E,col="blue",lwd=2,lty=1)
lines(rodent.means.rat_100CFU$S,col="forestgreen",lwd=2,lty=1)
lines(rodent.means.rat_100CFU$L,col="red",lwd=2,lty=5)
lines(rodent.means.rat_100CFU$R,col="blue",lwd=2,lty=5)
lines(rodent.means.rat_100CFU$Id,col="black",lwd=2,lty=1)
legend("topright",c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp =0.5, y.intersp =1)
title(main="Rat Blood, 100 CFU")




###############################################################################
### Now with Early Phase Only
###############################################################################


##params
params.mouse_1CFU.EPTonly <- c(alpha=1,lambdaA=1-exp(-0.035),lambdaB=1-exp(-0.20),lambdaC=1-exp(-0.07),
                       b=0.4,b1=2,tau=1-exp(-0.39), muf=1-exp(-0.02), mupb=1-exp(-0.14),mub=1-exp(-0.20), 
                       pep=0.03,ppb=0 ,pb=0 , tep=1, tpb=1, tb=1, mu=1-exp(-0.002),
                       sigma=1-exp(-0.25),gamma=1-exp(-0.14),epsilon=1-exp(-0.5))
params.rat_1CFU.EPTonly <- c(alpha=1,lambdaA=1-exp(-0.04),lambdaB=1-exp(-0.02),lambdaC=1-exp(-0.06),b=0.4,
                     b1=2,tau=1-exp(-0.48), muf=1-exp(-0.02), mupb=1-exp(-0.13),mub=1-exp(-0.26), pep=0.1,
                     ppb=0,pb=0,tep=1, tpb=1, tb=1, mu=1-exp(-0.002),sigma=1-exp(-0.25),
                     gamma=1-exp(-0.14),epsilon=1-exp(-0.5))	
params.mouse_10CFU.EPTonly <- c(alpha=1,lambdaA=1-exp(-0.035),lambdaB=1-exp(-0.20),lambdaC=1-exp(-0.07),b=0.4,
                        b1=2,tau=1-exp(-0.39), muf=1-exp(-0.02), mupb=1-exp(-0.14),mub=1-exp(-0.20),
                        pep=0.03,ppb=0 ,pb=0 , tep=0.001, tpb=0.5, tb=0.8, 
                        mu=1-exp(-0.002),sigma=1-exp(-0.25),gamma=1-exp(-0.14),epsilon=1-exp(-0.5))	
params.rat_10CFU.EPTonly <- c(alpha=1,lambdaA=1-exp(-0.04),lambdaB=1-exp(-0.02),lambdaC=1-exp(-0.06),b=0.4,
                      b1=2,tau=1-exp(-0.48), muf=1-exp(-0.02), mupb=1-exp(-0.13),mub=1-exp(-0.26),
                      pep=0.1,ppb=0 ,pb=0 ,tep=0.5, tpb=1 , tb=0.8 ,
                      mu=1-exp(-0.002),sigma=1-exp(-0.25),gamma=1-exp(-0.14),epsilon=1-exp(-0.5))	
params.mouse_100CFU.EPTonly <- c(alpha=1,lambdaA=1-exp(-0.035),lambdaB=1-exp(-0.20),lambdaC=1-exp(-0.07),b=0.4,
                         b1=2,tau=1-exp(-0.39), muf=1-exp(-0.02), mupb=1-exp(-0.14),mub=1-exp(-0.20),
                         pep=0.03,ppb=0 ,pb=0 , tep=0.00, tpb=0, tb=0.65, 
                         mu=1-exp(-0.002),sigma=1-exp(-0.25),gamma=1-exp(-0.14),epsilon=1-exp(-0.5))	
params.rat_100CFU.EPTonly <- c(alpha=1,lambdaA=1-exp(-0.04),lambdaB=1-exp(-0.02),lambdaC=1-exp(-0.06),b=0.4,
                       b1=2,tau=1-exp(-0.48), muf=1-exp(-0.02), mupb=1-exp(-0.13),mub=1-exp(-0.26),
                       pep=0.1,ppb=0,pb=0,tep=0, tpb=1 , tb=0.41 ,
                       mu=1-exp(-0.002),sigma=1-exp(-0.25),gamma=1-exp(-0.14),epsilon=1-exp(-0.5))	









######################################################################################################

output.mouse_1CFU.EPTonly <- plague.IBM(params.mouse_1CFU.EPTonly, # vector of params with names, see below
                                xstart, # vector of starting values with names
                                T,      # number of time steps to simulate
                                n.sim  # number of simulations
)

rodent.means.mouse_1CFU.EPTonly <- as.data.frame(output.mouse_1CFU.EPTonly$mean.rodent.ts)
plot.ts(rodent.means.mouse_1CFU.EPTonly$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(rodent.means.mouse_1CFU.EPTonly$E,col="blue",lwd=2,lty=1)
lines(rodent.means.mouse_1CFU.EPTonly$S,col="forestgreen",lwd=2,lty=1)
lines(rodent.means.mouse_1CFU.EPTonly$L,col="red",lwd=2,lty=5)
lines(rodent.means.mouse_1CFU.EPTonly$R,col="blue",lwd=2,lty=5)
lines(rodent.means.mouse_1CFU.EPTonly$Id,col="black",lwd=2,lty=1)
legend("topright",c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp =0.5, y.intersp =1)
title(main="Mouse Blood, 1 CFU, EPTonly")

#############
output.rat_1CFU.EPTonly <- plague.IBM(params.rat_1CFU.EPTonly, # vector of params with names, see below
                              xstart, # vector of starting values with names
                              T,      # number of time steps to simulate
                              n.sim  # number of simulations
)

rodent.means.rat_1CFU.EPTonly <- as.data.frame(output.rat_1CFU.EPTonly$mean.rodent.ts)
plot.ts(rodent.means.rat_1CFU.EPTonly$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(rodent.means.rat_1CFU.EPTonly$E,col="blue",lwd=2,lty=1)
lines(rodent.means.rat_1CFU.EPTonly$S,col="forestgreen",lwd=2,lty=1)
lines(rodent.means.rat_1CFU.EPTonly$L,col="red",lwd=2,lty=5)
lines(rodent.means.rat_1CFU.EPTonly$R,col="blue",lwd=2,lty=5)
lines(rodent.means.rat_1CFU.EPTonly$Id,col="black",lwd=2,lty=1)
legend("topright",c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp =0.5, y.intersp =1)
title(main="Rat Blood, 1 CFU, EPTonly")


#############
output.mouse_10CFU.EPTonly <- plague.IBM(params.mouse_10CFU.EPTonly, # vector of params with names, see below
                                 xstart, # vector of starting values with names
                                 T,      # number of time steps to simulate
                                 n.sim  # number of simulations
)

rodent.means.mouse_10CFU.EPTonly <- as.data.frame(output.mouse_10CFU.EPTonly$mean.rodent.ts)
plot.ts(rodent.means.mouse_10CFU.EPTonly$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(rodent.means.mouse_10CFU.EPTonly$E,col="blue",lwd=2,lty=1)
lines(rodent.means.mouse_10CFU.EPTonly$S,col="forestgreen",lwd=2,lty=1)
lines(rodent.means.mouse_10CFU.EPTonly$L,col="red",lwd=2,lty=5)
lines(rodent.means.mouse_10CFU.EPTonly$R,col="blue",lwd=2,lty=5)
lines(rodent.means.mouse_10CFU.EPTonly$Id,col="black",lwd=2,lty=1)
legend("topright",c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp =0.5, y.intersp =1)
title(main="Mouse Blood, 10 CFU,EPTonly")

###############
output.rat_10CFU.EPTonly <- plague.IBM(params.rat_10CFU.EPTonly, # vector of params with names, see below
                               xstart, # vector of starting values with names
                               T,      # number of time steps to simulate
                               n.sim  # number of simulations
)

rodent.means.rat_10CFU.EPTonly <- as.data.frame(output.rat_10CFU.EPTonly$mean.rodent.ts)
plot.ts(rodent.means.rat_10CFU.EPTonly$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(rodent.means.rat_10CFU.EPTonly$E,col="blue",lwd=2,lty=1)
lines(rodent.means.rat_10CFU.EPTonly$S,col="forestgreen",lwd=2,lty=1)
lines(rodent.means.rat_10CFU.EPTonly$L,col="red",lwd=2,lty=5)
lines(rodent.means.rat_10CFU.EPTonly$R,col="blue",lwd=2,lty=5)
lines(rodent.means.rat_10CFU.EPTonly$Id,col="black",lwd=2,lty=1)
legend("topright",c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp =0.5, y.intersp =1)
title(main="Rat Blood, 10 CFU, EPTonly")

#############
output.mouse_100CFU.EPTonly <- plague.IBM(params.mouse_100CFU.EPTonly, # vector of params with names, see below
                                  xstart, # vector of starting values with names
                                  T,      # number of time steps to simulate
                                  n.sim  # number of simulations
)

rodent.means.mouse_100CFU.EPTonly <- as.data.frame(output.mouse_100CFU.EPTonly$mean.rodent.ts)
plot.ts(rodent.means.mouse_100CFU.EPTonly$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(rodent.means.mouse_100CFU.EPTonly$E,col="blue",lwd=2,lty=1)
lines(rodent.means.mouse_100CFU.EPTonly$S,col="forestgreen",lwd=2,lty=1)
lines(rodent.means.mouse_100CFU.EPTonly$L,col="red",lwd=2,lty=5)
lines(rodent.means.mouse_100CFU.EPTonly$R,col="blue",lwd=2,lty=5)
lines(rodent.means.mouse_100CFU.EPTonly$Id,col="black",lwd=2,lty=1)
legend("topright",c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp =0.5, y.intersp =1)
title(main="Mouse Blood, 100 CFU, EPTonly")

###############
output.rat_100CFU.EPTonly <- plague.IBM(params.rat_100CFU.EPTonly, # vector of params with names, see below
                                xstart, # vector of starting values with names
                                T,      # number of time steps to simulate
                                n.sim  # number of simulations
)

rodent.means.rat_100CFU.EPTonly <- as.data.frame(output.rat_100CFU.EPTonly$mean.rodent.ts)
plot.ts(rodent.means.rat_100CFU.EPTonly$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(rodent.means.rat_100CFU.EPTonly$E,col="blue",lwd=2,lty=1)
lines(rodent.means.rat_100CFU.EPTonly$S,col="forestgreen",lwd=2,lty=1)
lines(rodent.means.rat_100CFU.EPTonly$L,col="red",lwd=2,lty=5)
lines(rodent.means.rat_100CFU.EPTonly$R,col="blue",lwd=2,lty=5)
lines(rodent.means.rat_100CFU.EPTonly$Id,col="black",lwd=2,lty=1)
legend("topright",c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp =0.5, y.intersp =1)
title(main="Rat Blood, 100 CFU, EPTonly")



###############################################################################
### Now with no early phase, only partially and fully blocked
###############################################################################


##params
params.mouse_1CFU.BPBonly <- c(alpha=1,lambdaA=1-exp(-0.035),lambdaB=1-exp(-0.20),lambdaC=1-exp(-0.07),
                       b=0.4,b1=2,tau=1-exp(-0.39), muf=1-exp(-0.02), mupb=1-exp(-0.14),mub=1-exp(-0.20), 
                       pep=0 ,ppb=0.11,pb=0.5, tep=1, tpb=1, tb=1, mu=1-exp(-0.002),
                       sigma=1-exp(-0.25),gamma=1-exp(-0.14),epsilon=1-exp(-0.5))
params.rat_1CFU.BPBonly <- c(alpha=1,lambdaA=1-exp(-0.04),lambdaB=1-exp(-0.02),lambdaC=1-exp(-0.06),b=0.4,
                     b1=2,tau=1-exp(-0.48), muf=1-exp(-0.02), mupb=1-exp(-0.13),mub=1-exp(-0.26), pep=0 ,
                     ppb=0.10,pb=0.67,tep=1, tpb=1, tb=1, mu=1-exp(-0.002),sigma=1-exp(-0.25),
                     gamma=1-exp(-0.14),epsilon=1-exp(-0.5))	
params.mouse_10CFU.BPBonly <- c(alpha=1,lambdaA=1-exp(-0.035),lambdaB=1-exp(-0.20),lambdaC=1-exp(-0.07),b=0.4,
                        b1=2,tau=1-exp(-0.39), muf=1-exp(-0.02), mupb=1-exp(-0.14),mub=1-exp(-0.20),
                        pep=0 ,ppb=0.11,pb=0.5, tep=0.001, tpb=0.5, tb=0.8, 
                        mu=1-exp(-0.002),sigma=1-exp(-0.25),gamma=1-exp(-0.14),epsilon=1-exp(-0.5))	
params.rat_10CFU.BPBonly <- c(alpha=1,lambdaA=1-exp(-0.04),lambdaB=1-exp(-0.02),lambdaC=1-exp(-0.06),b=0.4,
                      b1=2,tau=1-exp(-0.48), muf=1-exp(-0.02), mupb=1-exp(-0.13),mub=1-exp(-0.26),
                      pep=0 ,ppb=0.10,pb=0.67,tep=0.5, tpb=1 , tb=0.8 ,
                      mu=1-exp(-0.002),sigma=1-exp(-0.25),gamma=1-exp(-0.14),epsilon=1-exp(-0.5))	
params.mouse_100CFU.BPBonly <- c(alpha=1,lambdaA=1-exp(-0.035),lambdaB=1-exp(-0.20),lambdaC=1-exp(-0.07),b=0.4,
                         b1=2,tau=1-exp(-0.39), muf=1-exp(-0.02), mupb=1-exp(-0.14),mub=1-exp(-0.20),
                         pep=0 ,ppb=0.11,pb=0.5, tep=0.00, tpb=0, tb=0.65, 
                         mu=1-exp(-0.002),sigma=1-exp(-0.25),gamma=1-exp(-0.14),epsilon=1-exp(-0.5))	
params.rat_100CFU.BPBonly <- c(alpha=1,lambdaA=1-exp(-0.04),lambdaB=1-exp(-0.02),lambdaC=1-exp(-0.06),b=0.4,
                       b1=2,tau=1-exp(-0.48), muf=1-exp(-0.02), mupb=1-exp(-0.13),mub=1-exp(-0.26),
                       pep=0 ,ppb=0.10,pb=0.67,tep=0, tpb=1 , tb=0.41 ,
                       mu=1-exp(-0.002),sigma=1-exp(-0.25),gamma=1-exp(-0.14),epsilon=1-exp(-0.5))	









######################################################################################################

output.mouse_1CFU.BPBonly <- plague.IBM(params.mouse_1CFU.BPBonly, # vector of params with names, see below
                                xstart, # vector of starting values with names
                                T,      # number of time steps to simulate
                                n.sim  # number of simulations
)

rodent.means.mouse_1CFU.BPBonly <- as.data.frame(output.mouse_1CFU.BPBonly$mean.rodent.ts)
plot.ts(rodent.means.mouse_1CFU.BPBonly$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(rodent.means.mouse_1CFU.BPBonly$E,col="blue",lwd=2,lty=1)
lines(rodent.means.mouse_1CFU.BPBonly$S,col="forestgreen",lwd=2,lty=1)
lines(rodent.means.mouse_1CFU.BPBonly$L,col="red",lwd=2,lty=5)
lines(rodent.means.mouse_1CFU.BPBonly$R,col="blue",lwd=2,lty=5)
lines(rodent.means.mouse_1CFU.BPBonly$Id,col="black",lwd=2,lty=1)
legend("topright",c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp =0.5, y.intersp =1)
title(main="Mouse Blood, 1 CFU, BPBonly")

#############
output.rat_1CFU.BPBonly <- plague.IBM(params.rat_1CFU.BPBonly, # vector of params with names, see below
                              xstart, # vector of starting values with names
                              T,      # number of time steps to simulate
                              n.sim  # number of simulations
)

rodent.means.rat_1CFU.BPBonly <- as.data.frame(output.rat_1CFU.BPBonly$mean.rodent.ts)
plot.ts(rodent.means.rat_1CFU.BPBonly$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(rodent.means.rat_1CFU.BPBonly$E,col="blue",lwd=2,lty=1)
lines(rodent.means.rat_1CFU.BPBonly$S,col="forestgreen",lwd=2,lty=1)
lines(rodent.means.rat_1CFU.BPBonly$L,col="red",lwd=2,lty=5)
lines(rodent.means.rat_1CFU.BPBonly$R,col="blue",lwd=2,lty=5)
lines(rodent.means.rat_1CFU.BPBonly$Id,col="black",lwd=2,lty=1)
legend("topright",c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp =0.5, y.intersp =1)
title(main="Rat Blood, 1 CFU, BPBonly")


#############
output.mouse_10CFU.BPBonly <- plague.IBM(params.mouse_10CFU.BPBonly, # vector of params with names, see below
                                 xstart, # vector of starting values with names
                                 T,      # number of time steps to simulate
                                 n.sim  # number of simulations
)

rodent.means.mouse_10CFU.BPBonly <- as.data.frame(output.mouse_10CFU.BPBonly$mean.rodent.ts)
plot.ts(rodent.means.mouse_10CFU.BPBonly$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(rodent.means.mouse_10CFU.BPBonly$E,col="blue",lwd=2,lty=1)
lines(rodent.means.mouse_10CFU.BPBonly$S,col="forestgreen",lwd=2,lty=1)
lines(rodent.means.mouse_10CFU.BPBonly$L,col="red",lwd=2,lty=5)
lines(rodent.means.mouse_10CFU.BPBonly$R,col="blue",lwd=2,lty=5)
lines(rodent.means.mouse_10CFU.BPBonly$Id,col="black",lwd=2,lty=1)
legend("topright",c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp =0.5, y.intersp =1)
title(main="Mouse Blood, 10 CFU, BPBonly")

###############
output.rat_10CFU.BPBonly <- plague.IBM(params.rat_10CFU.BPBonly, # vector of params with names, see below
                               xstart, # vector of starting values with names
                               T,      # number of time steps to simulate
                               n.sim  # number of simulations
)

rodent.means.rat_10CFU.BPBonly <- as.data.frame(output.rat_10CFU.BPBonly$mean.rodent.ts)
plot.ts(rodent.means.rat_10CFU.BPBonly$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(rodent.means.rat_10CFU.BPBonly$E,col="blue",lwd=2,lty=1)
lines(rodent.means.rat_10CFU.BPBonly$S,col="forestgreen",lwd=2,lty=1)
lines(rodent.means.rat_10CFU.BPBonly$L,col="red",lwd=2,lty=5)
lines(rodent.means.rat_10CFU.BPBonly$R,col="blue",lwd=2,lty=5)
lines(rodent.means.rat_10CFU.BPBonly$Id,col="black",lwd=2,lty=1)
legend("topright",c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp =0.5, y.intersp =1)
title(main="Rat Blood, 10 CFU, BPBonly")

#############
output.mouse_100CFU.BPBonly <- plague.IBM(params.mouse_100CFU.BPBonly, # vector of params with names, see below
                                  xstart, # vector of starting values with names
                                  T,      # number of time steps to simulate
                                  n.sim  # number of simulations
)

rodent.means.mouse_100CFU.BPBonly <- as.data.frame(output.mouse_100CFU.BPBonly$mean.rodent.ts)
plot.ts(rodent.means.mouse_100CFU.BPBonly$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(rodent.means.mouse_100CFU.BPBonly$E,col="blue",lwd=2,lty=1)
lines(rodent.means.mouse_100CFU.BPBonly$S,col="forestgreen",lwd=2,lty=1)
lines(rodent.means.mouse_100CFU.BPBonly$L,col="red",lwd=2,lty=5)
lines(rodent.means.mouse_100CFU.BPBonly$R,col="blue",lwd=2,lty=5)
lines(rodent.means.mouse_100CFU.BPBonly$Id,col="black",lwd=2,lty=1)
legend("topright",c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp =0.5, y.intersp =1)
title(main="Mouse Blood, 100 CFU, BPBonly")

###############
output.rat_100CFU.BPBonly <- plague.IBM(params.rat_100CFU.BPBonly, # vector of params with names, see below
                                xstart, # vector of starting values with names
                                T,      # number of time steps to simulate
                                n.sim  # number of simulations
)

rodent.means.rat_100CFU.BPBonly <- as.data.frame(output.rat_100CFU.BPBonly$mean.rodent.ts)
plot.ts(rodent.means.rat_100CFU.BPBonly$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(rodent.means.rat_100CFU.BPBonly$E,col="blue",lwd=2,lty=1)
lines(rodent.means.rat_100CFU.BPBonly$S,col="forestgreen",lwd=2,lty=1)
lines(rodent.means.rat_100CFU.BPBonly$L,col="red",lwd=2,lty=5)
lines(rodent.means.rat_100CFU.BPBonly$R,col="blue",lwd=2,lty=5)
lines(rodent.means.rat_100CFU.BPBonly$Id,col="black",lwd=2,lty=1)
legend("topright",c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp =0.5, y.intersp =1)
title(main="Rat Blood, 100 CFU,BPBonly")




###############################################################################
## Compare how many rodents were dead at the end for the different scenarios
###############################################################################

RodentsDead <- data.frame(rodent=rep(c("mouse","rat"),3),scenario=rep(c("baseline","EPTonly","BPBonly"),each=2),CFU.1=rep(NA,6),CFU.10=rep(NA,6),CFU.100=rep(NA,6))

RodentsDead$CFU.1[1] <- rodent.means.mouse_1CFU$Id[T]
RodentsDead$CFU.10[1] <- rodent.means.mouse_10CFU$Id[T]
RodentsDead$CFU.100[1] <- rodent.means.mouse_100CFU$Id[T]
RodentsDead$CFU.1[2] <- rodent.means.rat_1CFU$Id[T]
RodentsDead$CFU.10[2] <- rodent.means.rat_10CFU$Id[T]
RodentsDead$CFU.100[2] <- rodent.means.rat_100CFU$Id[T]

RodentsDead$CFU.1[3] <- rodent.means.mouse_1CFU.EPTonly$Id[T]
RodentsDead$CFU.10[3] <- rodent.means.mouse_10CFU.EPTonly$Id[T]
RodentsDead$CFU.100[3] <- rodent.means.mouse_100CFU.EPTonly$Id[T]
RodentsDead$CFU.1[4] <- rodent.means.rat_1CFU.EPTonly$Id[T]
RodentsDead$CFU.10[4] <- rodent.means.rat_10CFU.EPTonly$Id[T]
RodentsDead$CFU.100[4] <- rodent.means.rat_100CFU.EPTonly$Id[T]

RodentsDead$CFU.1[5] <- rodent.means.mouse_1CFU.BPBonly$Id[T]
RodentsDead$CFU.10[5] <- rodent.means.mouse_10CFU.BPBonly$Id[T]
RodentsDead$CFU.100[5] <- rodent.means.mouse_100CFU.BPBonly$Id[T]
RodentsDead$CFU.1[6] <- rodent.means.rat_1CFU.BPBonly$Id[T]
RodentsDead$CFU.10[6] <- rodent.means.rat_10CFU.BPBonly$Id[T]
RodentsDead$CFU.100[6] <- rodent.means.rat_100CFU.BPBonly$Id[T]


RodentsDead





##################################################################
# using comparison function



params.mouse <- list(params.mouse_1CFU, params.mouse_1CFU.EPTonly, params.mouse_1CFU.BPBonly, 
                     params.mouse_10CFU, params.mouse_10CFU.EPTonly, params.mouse_10CFU.BPBonly,
                     params.mouse_100CFU, params.mouse_100CFU.EPTonly, params.mouse_100CFU.BPBonly)
names(params.mouse) <- c("params.mouse_1CFU","params.mouse_1CFU.EPTonly", "params.mouse_1CFU.BPBonly", 
                         "params.mouse_10CFU", "params.mouse_10CFU.EPTonly", "params.mouse_10CFU.BPBonly",
                         "params.mouse_100CFU", "params.mouse_100CFU.EPTonly", "params.mouse_100CFU.BPBonly")

mouse.comparisons <- print.plague.IBM(params.mouse, 
                                      xstart=c(S=9,L=0,I=1,E=0,R=0,dr=0,Id=0, U=100,Iep=0,Ipb=0,Ib=0,df=0),
                                      T=100,
                                      n.sim = 100,
                                      plot.name="Plague IBM mouse comparison.eps")

params.rat <- list(params.rat_1CFU, params.rat_1CFU.EPTonly, params.rat_1CFU.BPBonly, 
                   params.rat_10CFU, params.rat_10CFU.EPTonly, params.rat_10CFU.BPBonly,
                   params.rat_100CFU, params.rat_100CFU.EPTonly, params.rat_100CFU.BPBonly)
names(params.rat) <- c("params.rat_1CFU","params.rat_1CFU.EPTonly", "params.rat_1CFU.BPBonly", 
                       "params.rat_10CFU", "params.rat_10CFU.EPTonly", "params.rat_10CFU.BPBonly",
                       "params.rat_100CFU", "params.rat_100CFU.EPTonly", "params.rat_100CFU.BPBonly")

rat.comparisons <- print.plague.IBM(params.rat, 
                                    xstart=c(S=9,L=0,I=1,E=0,R=0,dr=0,Id=0, U=100,Iep=0,Ipb=0,Ib=0,df=0),
                                    T=100,
                                    n.sim = 100,
                                    plot.name="Plague IBM rat comparison.eps")

