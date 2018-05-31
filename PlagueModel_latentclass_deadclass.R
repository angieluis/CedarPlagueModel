###############################################################################
### Model for Plague by Cedar Mitchell Spring 2017
### Exploring the importance of early-phase in transmission
###############################################################################


library(deSolve)

SIR.model=function(t,x,params){
  b=params["b"]           # biting rate of fleas
  mu=params["mu"]		      # natural mortality rate of rodent (uninfected) (mu_r in fig)
  p_ep=params["p_ep"]       # Proportion of ep fleas that transmit
  p_pb=params["p_pb"]       # Proportion of partially blocked fleas that transmit
  p_b=params["p_b"]         # Proportion of fully blocked fleas that transmit
  gamma=params["gamma"]		# recovery rate in rodent from low-dose flea infection
  epsilon=params["epsilon"]  #disease induced mortality rate in rodent from high-dose flea infection
  sigma=params["sigma"]   # rate to become infectious from latent class
  
  alpha=params["alpha"]   # proportion of fleas infected from host
  mu_f=params["mu_f"]       # natural mortality of flea
  mu_pb=params["mu_pb"]     # mortality of partially blocked flea
  mu_b=params["mu_b"]       # mortality of blocked flea
  
  lambdaA=params["lambdaA"] # rate of developing partial blockage from EP
  lambdaB=params["lambdaB"] # rate of clearing infection in EP (back to uninfected)
  lambdaC=params["lambdaC"] # rate of leaving EP (still infected but not enough to block)
  tau=params["tau"]         # rate of developing full blockage
  
  beta=params["beta"]       # Infection rate of rodent OMITS EARLY-PHASE # this isn't specified in params but calculated below
  
  
  S=x[1]
  L=x[2]
  I=x[3]
  E=x[4]
  R=x[5]
  dr=x[6]
  
  
  U=x[7]
  Iep=x[8]
  Ipb=x[9]
  Ib=x[10]
  df=x[11]
  
  beta=(b*S*Ipb*ppb)+(b*S*Ib*pb)
  
  
  dS = -(beta + (b*S*Iep*p_ep))-mu*S   #Susceptible Rodent Pop
  dL = (beta + (b*S*0.1*Iep*p_ep)) - L*sigma  #infection of host, not yet infectious
  dI = sigma*L - I*epsilon  #infectious rodent hosts
  dE = (b*S*0.9*Iep*p_ep) - E*(mu+gamma) #exposed rodents (non-infectious)
  dR = gamma*E - R*mu  #recovery of rodent
  dr = mu*(S+R+E) + epsilon*I  #dead rodents
    
  dU = -I*(U*b*alpha) + (0.21*Iep*lambdaB) -U*muf    #Uninfected Flea pop
  dIep = (I*(U*b*alpha)) - ((0.21*Iep*lambdaB) + (0.13*Iep*lambdaC) + (0.66*Iep*lambdaA) + Iep*mu_f) #Infection of early phase fleas 
  dIpb = 0.66*Iep*lambdaA - Ipb*(tau+mu_pb) #Infection of partially blocked fleas
  dIb = tau*Ipb - mub*Ib  #Infection of fully blocked fleas
  df = m_uf*(U+Iep) + (mu_pb*Ipb) + mu_b*Ib + 0.13*Iep*lambdaC #dead fleas
  
  list(c(dS,dL,dI,dE,dR,dr, dU,dIep,dIpb,dIb,df))
  
}
          
    

times=seq(1,30,by=0.3)					#time steps to output
xstart=c(S=9, L=0, I=0, E=0, R=0, dr=0, U=90, Iep=10, Ipb=0, Ib=0, df=0)	#beginning population sizes

parms=c(         # all parameter values are specified in days
  b=.48,              # biting rate of fleas 
  mu=1/365,           # natural mortality rate of rodent (uninfected)
  p_ep=2/83,           # Proportion of ep fleas that transmit
  p_pb=5/26,           # Proportion of partially blocked fleas that transmit
  p_b=10/20,           # Proportion of fully blocked fleas that transmit
  sigma=1/2.6,        # rate to become infectious from latent class
  gamma=1/7,          # recovery rate in rodent from low-dose flea infection
  epsilon=1/2,        # disease induced mortality rate in rodent from high-dose flea infection
  mu_f=1/20.6,         # natural mortality of flea
  alpha=19/20,        # proportion of fleas infected from host
  lambdaA=1/5.05,     # rate of developing partial blockage from EP
  lambdaB=1/3,        # rate of clearing infection in EP (back to uninfected)
  lambdaC=1/3,        # rate of leaving EP (still infected but not enough to block)
  mu_pb=1/12.3,        # mortality of partially blocked flea
  tau=1/4.8,          # rate of developing full blockage
  mu_b=1/6)	          # mortality of blocked flea

out=lsoda(xstart, times, SIR.model, parms)  #run the model

time=out[,1]
S=out[,2]
L=out[,3]
I=out[,4]	
E=out[,5]
R=out[,6]
dr=out[,7]

U=out[,8]
Iep=out[,9]
Ipb=out[,10]
Ib=out[,11]
df=out[,12]

par(mfrow=c(2,1))
plot(time,S,ylab="Abundance",xlab="time",type="l",col="olivedrab3",lwd=2,ylim=c(0,max(10)),lty=2)
lines(time,L,col="ivory3",lwd=2)
lines(time,I,col="tomato3",lwd=2)
lines(time,E,col="black",lwd=2)
lines(time,R,col="steelblue",lwd=2,lty=3)
lines(time,dr,col="grey66",lwd=2, lty=4)
legend(25,10,c("S","L","I","E","R","dr"),col=c("olivedrab3","ivory3","tomato3","black","steelblue","grey66"),bty="n",lty=c(2,1,1,1,3,4),lwd=2)
title(main="Host Dynamics")


plot(time,U,col="burlywood1",lwd=2,lty=2,ylab="Abundance",xlab="time", type="l",ylim=c(0,max(100)))
lines(time,Iep,col="darkgreen",lwd=2,lty=1)
lines(time,Ipb,col="plum2",lwd=2,lty=1)
lines(time,Ib,col="orchid4",lwd=2,lty=1)
lines(time,df,col="grey66",lwd=2, lty=4)
title(main="Flea Dynamics")


legend(25,100,c("U","Iep","Ipb","Ib","df"),col=c("burlywood1","darkgreen","plum2","orchid4","grey66"),bty="n",lty=c(2,1,1,1,4),lwd=2,x.intersp = 1, y.intersp = 0.75)
out
colors() #many color options

