###############################################################################
### Model for Plague by Cedar Mitchell Spring 2017
### Exploring the importance of early-phase vs blocked fleas in transmission
###############################################################################


library(deSolve)

SIR.model=function(t,x,params){
  b = params["b"]           # biting rate of all fleas except blocked
  b_b = params["b_b"]       # biting rate of blocked fleas 
  mu = params["mu"]		      # natural mortality rate of rodent (uninfected) (mu_r in fig)
  p_ep = params["p_ep"]     # Proportion of ep fleas that transmit
  p_pb = params["p_pb"]     # Proportion of partially blocked fleas that transmit
  p_b = params["p_b"]       # Proportion of fully blocked fleas that transmit
  gamma = params["gamma"]		# recovery rate in rodent from low-dose flea infection
  epsilon = params["epsilon"]  #disease induced mortality rate in rodent from high-dose flea infection
  sigma = params["sigma"]   # rate to become infectious from latent class
  
  rho = params["rho"]       # transmisison efficiency of early phase fleas
  
  alpha = params["alpha"]     # proportion of fleas infected from host
  mu_f = params["mu_f"]       # natural mortality of flea
  mu_pb = params["mu_pb"]     # mortality of partially blocked flea
  mu_b = params["mu_b"]       # mortality of blocked flea
  
  lambdaA = params["lambdaA"] # rate of developing partial blockage from EP
  lambdaB = params["lambdaB"] # rate of clearing infection in EP (back to uninfected)
  lambdaC = params["lambdaC"] # rate of leaving EP (still infected but not enough to block)
  tau = params["tau"]         # rate of developing full blockage (from partial blockage)
  
  
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
  
  
  dS = -S*(b*p_pb*Ipb + b_b*p_b*Ib + b*p_ep*Iep) - mu*S   #Susceptible Rodent Pop
  dL = S*(b*p_pb*Ipb + b_b*p_b*Ib + b*p_ep*Iep*rho) - L*sigma  # doesn't have death here (-mu*L)   #infection of host, not yet infectious
  dI = sigma*L - I*epsilon  #infectious rodent hosts
  dE = S*b*p_ep*Iep*(1-rho) - E*(mu+gamma) #exposed rodents (non-infectious)
  dR = gamma*E - R*mu  # recovered rodents
  dr = mu*(S+R+E) + epsilon*I  #dead rodents
    
  dU = -b*alpha*U*I + (Iep*lambdaB) - U*mu_f   #Uninfected Flea pop
  dIep = b*alpha*U*I - ((Iep*lambdaB) + (Iep*lambdaC) + (Iep*lambdaA) + Iep*mu_f) #Infection of early phase fleas 
  dIpb = Iep*lambdaA - Ipb*(tau+mu_pb) #Infection of partially blocked fleas
  dIb = tau*Ipb - mu_b*Ib  #Infection of fully blocked fleas
  df = mu_f*(U+Iep) + (mu_pb*Ipb) + mu_b*Ib + Iep*lambdaC # lost fleas
  
  list(c(dS,dL,dI,dE,dR,dr, dU,dIep,dIpb,dIb,df))
  
}
          
    

times=seq(1,30,by=0.3)					#time steps to output
xstart=c(S=9, L=0, I=0, E=0, R=0, dr=0, U=90, Iep=10, Ipb=0, Ib=0, df=0)	#beginning population sizes

old.parms=c(         # all parameter values are specified in days
  b=0.48,             # biting rate of fleas (uninfected,early-phase, & partially blocked)
  b_b=0.9,            # biting rate of blocked fleas 
  mu=1/365,           # natural mortality rate of rodent (uninfected)
  p_ep=2/83,          # Proportion of ep fleas that transmit
  p_pb=5/26,          # Proportion of partially blocked fleas that transmit
  p_b=10/20,          # Proportion of fully blocked fleas that transmit
  sigma=1/2.6,        # rate to become infectious from latent class
  gamma=1/7,          # recovery rate in rodent from low-dose flea infection
  epsilon=1/2,        # disease induced mortality rate in rodent from high-dose flea infection
  rho=0.1,
  mu_f=1/20.6,        # natural mortality of flea
  alpha=19/20,        # proportion of fleas infected from host
  lambdaA=1/5.05,     # rate of developing partial blockage from EP
  lambdaB=1/3,        # rate of clearing infection in EP (back to uninfected)
  lambdaC=1/3,        # rate of leaving EP (still infected but not enough to block)
  mu_pb=1/12.3,       # mortality of partially blocked flea
  tau=1/4.8,          # rate of developing full blockage
  mu_b=1/6)	          # mortality of blocked flea

mouse.parms=c(         # all parameter values are specified in days
  b=0.39,             # biting rate of fleas (uninfected,early-phase, & partially blocked)
  b_b=0.9,            # biting rate of blocked fleas 
  mu=1/365,           # natural mortality rate of rodent (uninfected)
  p_ep=0.11,          # Proportion of ep fleas that transmit
  p_pb=0.17,          # Proportion of partially blocked fleas that transmit
  p_b=0.63,           # Proportion of fully blocked fleas that transmit
  sigma=0.5,          # rate to become infectious from latent class
  gamma=0.065,        # recovery rate in rodent from low-dose flea infection
  epsilon=0.8,        # disease induced mortality rate in rodent from high-dose flea infection
  rho=0.21,
  mu_f=0.029,        # natural mortality of flea
  alpha=0.97,        # proportion of fleas infected from host
  lambdaA=0.033,     # rate of developing partial blockage from EP
  lambdaB=0.049,     # rate of clearing infection in EP (back to uninfected)
  lambdaC=0.099,     # rate of leaving EP (still infected but not enough to block)
  mu_pb=0.19,        # mortality of partially blocked flea
  tau=0.42,          # rate of developing full blockage
  mu_b=0.22)	       # mortality of blocked flea

rat.parms=c(         # all parameter values are specified in days
  b=0.39,             # biting rate of fleas (uninfected,early-phase, & partially blocked)
  b_b=0.9,            # biting rate of blocked fleas 
  mu=1/365,           # natural mortality rate of rodent (uninfected)
  p_ep=0.09,          # Proportion of ep fleas that transmit
  p_pb=0.14,          # Proportion of partially blocked fleas that transmit
  p_b=0.73,           # Proportion of fully blocked fleas that transmit
  sigma=0.5,          # rate to become infectious from latent class
  gamma=0.065,        # recovery rate in rodent from low-dose flea infection
  epsilon=0.8,        # disease induced mortality rate in rodent from high-dose flea infection
  rho=0.21,
  mu_f=0.029,        # natural mortality of flea
  alpha=1,           # proportion of fleas infected from host
  lambdaA=0.044,     # rate of developing partial blockage from EP
  lambdaB=0.007,     # rate of clearing infection in EP (back to uninfected)
  lambdaC=0.058,     # rate of leaving EP (still infected but not enough to block)
  mu_pb=0.16,        # mortality of partially blocked flea
  tau=0.55,          # rate of developing full blockage
  mu_b=0.26)	       # mortality of blocked flea

out=lsoda(xstart, times, SIR.model, rat.parms)  #run the model

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





################################################### R0

#mouse
N_f = 1:20
N_r = 1
R0 <- with(as.list(mouse.parms), sqrt(((b*rho*p_ep + b*p_pb + b_b*p_b)*sigma/(epsilon*(sigma+mu))) * (b*alpha/(lambdaB + lambdaC + lambdaA + mu_f + tau + mu_pb + mu_b)) * N_r * N_f))

R0_ep <- with(as.list(mouse.parms), sqrt(((b*rho*p_ep + b*0 + b_b*0)*sigma/(epsilon*(sigma+mu))) * (b*alpha/(lambdaB + lambdaC + lambdaA + mu_f + tau + mu_pb + mu_b)) * N_r * N_f))

R0_pb <- with(as.list(mouse.parms), sqrt(((b*rho*0 + b*p_pb + b_b*0)*sigma/(epsilon*(sigma+mu))) * (b*alpha/(lambdaB + lambdaC + lambdaA + mu_f + tau + mu_pb + mu_b)) * N_r * N_f))

R0_b <- with(as.list(mouse.parms), sqrt(((b*rho*0 + b*0 + b_b*p_b)*sigma/(epsilon*(sigma+mu))) * (b*alpha/(lambdaB + lambdaC + lambdaA + mu_f + tau + mu_pb + mu_b)) * N_r * N_f))

par(mfrow=c(1,1))
plot(N_f,R0,type="l",ylim=c(0,2.5))
lines(N_f,R0_ep,col="blue")
lines(N_f,R0_pb,col="green")
lines(N_f,R0_b,col="red")
legend("topleft",c("full","blocked","partially blocked","earlyphase"))





#rat
N_f = 1:20
N_r = 1
R0 <- with(as.list(rat.parms), sqrt(((b*rho*p_ep + b*p_pb + b_b*p_b)*sigma/(epsilon*(sigma+mu))) * (b*alpha/(lambdaB + lambdaC + lambdaA + mu_f + tau + mu_pb + mu_b)) * N_r * N_f))

R0_ep <- with(as.list(rat.parms), sqrt(((b*rho*p_ep + b*0 + b_b*0)*sigma/(epsilon*(sigma+mu))) * (b*alpha/(lambdaB + lambdaC + lambdaA + mu_f + tau + mu_pb + mu_b)) * N_r * N_f))

R0_pb <- with(as.list(rat.parms), sqrt(((b*rho*0 + b*p_pb + b_b*0)*sigma/(epsilon*(sigma+mu))) * (b*alpha/(lambdaB + lambdaC + lambdaA + mu_f + tau + mu_pb + mu_b)) * N_r * N_f))

R0_b <- with(as.list(rat.parms), sqrt(((b*rho*0 + b*0 + b_b*p_b)*sigma/(epsilon*(sigma+mu))) * (b*alpha/(lambdaB + lambdaC + lambdaA + mu_f + tau + mu_pb + mu_b)) * N_r * N_f))

par(mfrow=c(1,1))
plot(N_f,R0,type="l",ylim=c(0,2.5))
lines(N_f,R0_ep,col="blue")
lines(N_f,R0_pb,col="green")
lines(N_f,R0_b,col="red")
legend("topleft",c("all modes","blocked","partially blocked","early-phase"),col=c("black","red","green","blue"),lty=1,bty="n")
