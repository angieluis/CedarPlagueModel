library(deSolve)
library(tidyverse)
SIR.model.EPT=function(t,x,params){
  b=params["b"] #biting rate of fleas
  b1=params["b1"]  #biting rate of blocked fleas
  tep=params["tep"]  #proprotion of INFECTIOUS bites from early-phase fleas
  tpb=params["tpb"]  #proprotion of INFECTIOUS bites from partially blocked fleas
  tb=params["tb"]  #proprotion of INFECTIOUS bites from fully blocked fleas
  mu=params["mu"]		#natural mortality rate of rodent (uninfected)
  pep=params["pep"] #Proportion of ep fleas that transmit
  ppb=params["ppb"]  #Proportion of partially blocked fleas that transmit
  pb=params["pb"]  #Proportion of fully blocked fleas that transmit
  gamma=params["gamma"]		#recovery rate in rodent from low-dose flea infection
  epsilon=params["epsilon"]  #disease induced mortality rate in rodent from high-dose flea infection
  sigma=params["sigma"]  #rate to become infectious from latent class
  
  alpha=params["alpha"] #proportion of fleas infected from host
  muf=params["muf"] #natural mortality of flea
  mupb=params["mupb"] #mortality of partially blocked flea
  mub=params["mub"] #mortality of blocked flea
  
  lambdaA=params["lambdaA"] #rate of developing partial blockage from EP
  lambdaB=params["lambdaB"] #rate of clearing infection in EP (back to uninfected)
  lambdaC=params["lambdaC"] #rate of leaving EP (still infected but not enough to block)
  tau=params["tau"] #rate of developing full blockage
  
  
  S1=x[1]
  L1=x[2]
  I1=x[3]
  E1=x[4]
  R1=x[5]
  dr1=x[6]
  Id1=x[7]
  
  U1=x[8]
  Iep1=x[9]
  Ipb1=x[10]
  Ib1=x[11]
  df1=x[12]
  
  
  dS1=-S1*((b*Iep1*pep)+(b*Ipb1*ppb)+(b1*Ib1*pb))-mu*S1   #Susceptible Rodent Pop
  dL1=(b*S1*Iep1*pep*tep + b*S1*Ipb1*ppb*tpb + b1*S1*Ib1*pb*tb) - sigma*L1 - mu*L1  #infection of host, not yet infectious
  dI1=sigma*L1 - I1*mu - I1*epsilon  #infectious rodent hosts
  dE1=S1*(Iep1*b*(1-tep)*pep + Ipb1*b*(1-tpb)*ppb + Ib1*b1*(1-tb)*pb) - E1*(mu+gamma) #exposed rodents (non-infectious)
  dR1=gamma*E1 - R1*mu  #recovery of rodent
  dr1=mu*(S1+R1+E1+L1+I1) + epsilon*I1  #dead rodents
  Id1=epsilon*I1
  
  dU1=-I1*(U1*b*alpha) + Iep1*lambdaB - U1*muf    #Uninfected Flea pop
  dIep1=I1*(U1*b*alpha) - (Iep1*lambdaB + Iep1*lambdaC + Iep1*lambdaA + Iep1*muf) #Infection of early phase fleas 
  dIpb1=Iep1*lambdaA - Ipb1*(tau+mupb) #Infection of partially blocked fleas
  dIb1=tau*Ipb1 - mub*Ib1  #Infection of fully blocked fleas
  df1=muf*(U1+Iep1) + mupb*Ipb1 + mub*Ib1 + Iep1*lambdaC #dead fleas
  
  
  
  list(c(dS1,dL1,dI1,dE1,dR1,dr1,Id1, dU1,dIep1,dIpb1,dIb1,df1))
  
}


times=seq(1,30,by=0.3)					#time steps to output
xstart=c(S1=99,L1=0,I1=1,E1=0,R1=0,dr1=0,Id1=0, U1=500,Iep1=0,Ipb1=0,Ib1=0,df1=0)					#beginning population size
parms=c(alpha=1,lambdaA=0.04,lambdaB=0.20,lambdaC=0.07,b=0.4,b1=1,tau=0.39, muf=0.02, mupb=0.14,mub=0.20,
        pep=0.05,ppb=0.11,pb=0.50, tep=1, tpb=1, tb=1, 
        mu=0.002,sigma=0.25,gamma=0.14,epsilon=0.5)	#set parameter values
out=lsoda(xstart, times, SIR.model.EPT, parms)  #run the model

time=out[,1]
S1=out[,2]
L1=out[,3]
I1=out[,4]	
E1=out[,5]
R1=out[,6]
dr1=out[,7]
Id1=out[,8]

U1=out[,9]
Iep1=out[,10]
Ipb1=out[,11]
Ib1=out[,12]
df1=out[,13]

#Verify that all fleas at end add up to starting number;
sum(df1[97]+Ib1[97]+Ipb1[97]+U1[97]+Iep1[97])
#Verify that all rodents add up to starting number at end of outbreak;
sum(S1[97], L1[97], I1[97], E1[97], R1[97], dr1[97])

#Create dataframe of infection output by time for each model state;
#mouse_SIR<-cbind.data.frame(time, S1, L1, I1, E1, R1, dr1, U1, Iep1, Ipb1, Ib1, df1)

#................................................
#Infected and Maintained on Rat Blood####
#................................................

SIR.model=function(t,x,params){
  b=params["b"] #biting rate of fleas
  b1=params["b1"]  #biting rate of blocked fleas
  tep=params["tep"]  #proprotion of INFECTIOUS bites from early-phase fleas
  tpb=params["tpb"]  #proprotion of INFECTIOUS bites from partially blocked fleas
  tb=params["tb"]  #proprotion of INFECTIOUS bites from fully blocked fleas
  mu=params["mu"]		#natural mortality rate of rodent (uninfected)
  pep=params["pep"] #Proportion of ep fleas that transmit
  ppb=params["ppb"]  #Proportion of partially blocked fleas that transmit
  pb=params["pb"]  #Proportion of fully blocked fleas that transmit
  gamma=params["gamma"]		#recovery rate in rodent from low-dose flea infection
  epsilon=params["epsilon"]  #disease induced mortality rate in rodent from high-dose flea infection
  sigma=params["sigma"]  #rate to become infectious from latent class
  
  alpha=params["alpha"] #proportion of fleas infected from host
  muf=params["muf"] #natural mortality of flea
  mupb=params["mupb"] #mortality of partially blocked flea
  mub=params["mub"] #mortality of blocked flea
  
  lambdaA=params["lambdaA"] #rate of developing partial blockage from EP
  lambdaB=params["lambdaB"] #rate of clearing infection in EP (back to uninfected)
  lambdaC=params["lambdaC"] #rate of leaving EP (still infected but not enough to block)
  tau=params["tau"] #rate of developing full blockage
  
  S=x[1]
  L=x[2]
  I=x[3]
  E=x[4]
  R=x[5]
  dr=x[6]
  Id=x[7]
  
  U=x[8]
  Iep=x[9]
  Ipb=x[10]
  Ib=x[11]
  df=x[12]
  
  
  dS=-S*(b*Iep*pep + b*Ipb*ppb + b1*Ib*pb) - mu*S   #Susceptible Rodent Pop
  dL=S*(b*Iep*pep*tep + b*Ipb*ppb*tpb + b1*Ib*pb*tb) - sigma*L - mu*L  #infection of host, not yet infectious
  dI=sigma*L - I*(mu+epsilon)  #infectious rodent hosts
  dE=S*(b*Iep*pep*(1-tep) + b*Ipb*ppb*(1-tpb) + b1*Ib*pb*(1-tb)) - E*(mu+gamma) #exposed rodents (non-infectious)
  dR=gamma*E - R*mu  #recovery of rodent
  dr=mu*(S+R+E+L+I) + epsilon*I  #dead rodents
  Id=epsilon*I
  
  dU=-I*(U*b*alpha) + Iep*lambdaB - U*muf    #Uninfected Flea pop
  dIep=I*U*b*alpha - (Iep*lambdaB + Iep*lambdaC + Iep*lambdaA + Iep*muf) #Infection of early phase fleas 
  dIpb=Iep*lambdaA - Ipb*(tau+mupb) #Infection of partially blocked fleas
  dIb=tau*Ipb - mub*Ib  #Infection of fully blocked fleas
  df=muf*(U+Iep) + mupb*Ipb + mub*Ib+Iep*lambdaC #dead fleas
  
  
  
  list(c(dS,dL,dI,dE,dR,dr,Id, dU,dIep,dIpb,dIb,df))
  
}


times=seq(1,30,by=0.3)					#time steps to output
xstart=c(S=99,L=0,I=1,E=0,R=0,dr=0,Id=0, U=500,Iep=0,Ipb=0,Ib=0,df=0)					#beginning population size
parms=c(alpha=1,lambdaA=0.04,lambdaB=0.02,lambdaC=0.06,b=0.4,b1=1,tau=0.48, muf=0.02, mupb=0.13,mub=0.26,
        pep=0.10,ppb=0.10,pb=0.67,tep=1, tpb=1 , tb=1 ,
        mu=0.002,sigma=0.25,gamma=0.14,epsilon=0.5)	#set parameter values
out=lsoda(xstart, times, SIR.model, parms)  #run the model

time=out[,1]
S=out[,2]
L=out[,3]
I=out[,4]	
E=out[,5]
R=out[,6]
dr=out[,7]
Id=out[,8]

U=out[,9]
Iep=out[,10]
Ipb=out[,11]
Ib=out[,12]
df=out[,13]

#rat_SIR<-cbind.data.frame(time, S, L, I, E, R, dr, U, Iep, Ipb, Ib, df)


#..................................................................
#Graph the model output####
#..................................................................



#HOST Compartment;
#par(mfrow=c(2,1))
#postscript("C:/Users/Cedar/Desktop/Plague/Models/Graphs/Host_mouseblood_1CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)
postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Host_mouseblood_1CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(time,I1,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(100)),lty=1)
lines(time,E1,col="blue",lwd=2,lty=1)
lines(time,S1,col="forestgreen",lwd=2,lty=1)
lines(time,L1,col="red",lwd=2,lty=5)
lines(time,R1,col="blue",lwd=2,lty=5)
lines(time,Id1,col="black",lwd=2,lty=1)
legend(23,50,c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp =0.5, y.intersp =1)
title(main="Host Dynamics with Mouse Blood Infected Fleas, 99UH1IH 500UF LD50 of 1 CFU")
dev.off()
#rat plot;
#postscript("C:/Users/Cedar/Desktop/Plague/Models/Graphs/Host_ratblood_1CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)
postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Host_ratblood_1CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(time,I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(100)),lty=1)
lines(time,E,col="blue",lwd=2,lty=1)
lines(time,S,col="forestgreen",lwd=2,lty=1)
lines(time,L,col="red",lwd=2,lty=5)
lines(time,R,col="blue",lwd=2,lty=5)
lines(time,Id,col="black",lwd=2,lty=1)
legend(23,50,c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Host Dynamics with Rat Blood Infected Fleas, 99UH1IH 500UF LD50 of 1 CFU")
dev.off()

#............................................................................
#FLEA Compartment

#par(mfrow=c(2,1))
#postscript("C:/Users/Cedar/Desktop/Plague/Models/Graphs/Flea_mouseblood_1CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)
postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Flea_mouseblood_1CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(time,U1,ylab="Abundance",xlab="Time",type="l",col="forestgreen",lwd=2,ylim=c(0,max(500)),lty=1)
lines(time,Iep1,col="blue",lwd=2,lty=1)
lines(time,Ipb1,col="red",lwd=2,lty=5)
lines(time,Ib1,col="red",lwd=2,lty=1)
lines(time,df1,col="black",lwd=2,lty=1)
legend(20,100,c("Uninfected", "Infectious-Early Phase", "Infectious-Partially Blocked","Infectious-Fully Blocked", "Cumulative Dead"),
       col=c("forestgreen","blue", "red","red", "black"),bty="n",lty=c(1,1, 5, 1, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Flea Dynamics with Mouse Blood Infected Fleas, 99UH1IH 500UF LD50 of 1 CFU")
dev.off()
#rat plot;
#postscript("C:/Users/Cedar/Desktop/Plague/Models/Graphs/Flea_ratblood_1CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)
postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Flea_ratblood_1CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(time,U,ylab="Abundance",xlab="Time",type="l",col="forestgreen",lwd=2,ylim=c(0,max(500)),lty=1)
lines(time,Iep,col="blue",lwd=2,lty=1)
lines(time,Ipb,col="red",lwd=2,lty=5)
lines(time,Ib,col="red",lwd=2,lty=1)
lines(time,df,col="black",lwd=2,lty=1)
legend(20,100,c("Uninfected", "Infectious-Early Phase", "Infectious-Partially Blocked","Infectious-Fully Blocked", "Cumulative Dead"),
       col=c("forestgreen","blue", "red","red", "black"),bty="n",lty=c(1,1, 5, 1, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Flea Dynamics with Rat Blood Infected Fleas, 99UH1IH 500UF LD50 of 1 CFU")
dev.off()

#..................................................
#Threshold of 10 CFU####
#..................................................
SIR.model.EPT=function(t,x,params){
  b=params["b"] #biting rate of fleas
  b1=params["b1"]  #biting rate of blocked fleas
  tep=params["tep"]  #proprotion of INFECTIOUS bites from early-phase fleas
  tpb=params["tpb"]  #proprotion of INFECTIOUS bites from partially blocked fleas
  tb=params["tb"]  #proprotion of INFECTIOUS bites from fully blocked fleas
  mu=params["mu"]		#natural mortality rate of rodent (uninfected)
  pep=params["pep"] #Proportion of ep fleas that transmit
  ppb=params["ppb"]  #Proportion of partially blocked fleas that transmit
  pb=params["pb"]  #Proportion of fully blocked fleas that transmit
  gamma=params["gamma"]		#recovery rate in rodent from low-dose flea infection
  epsilon=params["epsilon"]  #disease induced mortality rate in rodent from high-dose flea infection
  sigma=params["sigma"]  #rate to become infectious from latent class
  
  alpha=params["alpha"] #proportion of fleas infected from host
  muf=params["muf"] #natural mortality of flea
  mupb=params["mupb"] #mortality of partially blocked flea
  mub=params["mub"] #mortality of blocked flea
  
  lambdaA=params["lambdaA"] #rate of developing partial blockage from EP
  lambdaB=params["lambdaB"] #rate of clearing infection in EP (back to uninfected)
  lambdaC=params["lambdaC"] #rate of leaving EP (still infected but not enough to block)
  tau=params["tau"] #rate of developing full blockage
  
  S1=x[1]
  L1=x[2]
  I1=x[3]
  E1=x[4]
  R1=x[5]
  dr1=x[6]
  Id1=x[7]
  
  U1=x[8]
  Iep1=x[9]
  Ipb1=x[10]
  Ib1=x[11]
  df1=x[12]
  
  
  dS1=-S1*(b*Iep1*pep + b*Ipb1*ppb +b1*Ib1*pb) - mu*S1   #Susceptible Rodent Pop
  dL1=S1*(b*Iep1*pep*tep + b*Ipb1*ppb*tpb + b1*Ib1*pb*tb) - sigma*L1 - mu*L1  #infection of host, not yet infectious
  dI1=sigma*L1 - I1*(mu+epsilon)  #infectious rodent hosts
  dE1=S1*(b*Iep1*pep*(1-tep) + b*Ipb1*ppb*(1-tpb) + b1*Ib1*pb*(1-tb)) - E1*(mu+gamma) #exposed rodents (non-infectious)
dR1=gamma*E1 - R1*mu  #recovery of rodent
dr1=mu*(S1+R1+E1+L1+I1)+epsilon*I1  #dead rodents
Id1=epsilon*I1

dU1=-I1*(U1*b*alpha) + Iep1*lambdaB - U1*muf    #Uninfected Flea pop
dIep1=I1*(U1*b*alpha) - (Iep1*lambdaB + Iep1*lambdaC + Iep1*lambdaA + Iep1*muf) #Infection of early phase fleas 
dIpb1=Iep1*lambdaA - Ipb1*(tau+mupb) #Infection of partially blocked fleas
dIb1=tau*Ipb1 - mub*Ib1  #Infection of fully blocked fleas
df1=muf*(U1+Iep1) + mupb*Ipb1 + mub*Ib1 + Iep1*lambdaC #dead fleas



list(c(dS1,dL1,dI1,dE1,dR1,dr1,Id1, dU1,dIep1,dIpb1,dIb1,df1))

}


times=seq(1,30,by=0.3)					#time steps to output
xstart=c(S1=99,L1=0,I1=1,E1=0,R1=0,dr1=0,Id1=0, U1=500,Iep1=0,Ipb1=0,Ib1=0,df1=0)					#beginning population size
parms=c(alpha=1,lambdaA=0.04,lambdaB=0.20,lambdaC=0.07,b=0.4,b1=1,tau=0.39, muf=0.02, mupb=0.14,mub=0.20,
        pep=0.05,ppb=0.11,pb=0.50, tep=0, tpb=0.5, tb=0.80, 
        mu=0.002,sigma=0.25,gamma=0.14,epsilon=0.5)	#set parameter values
out=lsoda(xstart, times, SIR.model.EPT, parms)  #run the model

time=out[,1]
S1=out[,2]
L1=out[,3]
I1=out[,4]	
E1=out[,5]
R1=out[,6]
dr1=out[,7]
Id1=out[,8]

U1=out[,9]
Iep1=out[,10]
Ipb1=out[,11]
Ib1=out[,12]
df1=out[,13]

#Verify that all fleas at end add up to starting number;
sum(df1[97]+Ib1[97]+Ipb1[97]+U1[97]+Iep1[97])
#Verify that all rodents add up to starting number at end of outbreak;
sum(S1[97], L1[97], I1[97], E1[97], R1[97], dr1[97])

#Create dataframe of infection output by time for each model state;
#mouse_SIR<-cbind.data.frame(time, S1, L1, I1, E1, R1, dr1, U1, Iep1, Ipb1, Ib1, df1)

#................................................
#Infected and Maintained on Rat Blood####
#................................................

SIR.model=function(t,x,params){
  b=params["b"] #biting rate of fleas
  b1=params["b1"]  #biting rate of blocked fleas
  tep=params["tep"]  #proprotion of INFECTIOUS bites from early-phase fleas
  tpb=params["tpb"]  #proprotion of INFECTIOUS bites from partially blocked fleas
  tb=params["tb"]  #proprotion of INFECTIOUS bites from fully blocked fleas
  mu=params["mu"]		#natural mortality rate of rodent (uninfected)
  pep=params["pep"] #Proportion of ep fleas that transmit
  ppb=params["ppb"]  #Proportion of partially blocked fleas that transmit
  pb=params["pb"]  #Proportion of fully blocked fleas that transmit
  gamma=params["gamma"]		#recovery rate in rodent from low-dose flea infection
  epsilon=params["epsilon"]  #disease induced mortality rate in rodent from high-dose flea infection
  sigma=params["sigma"]  #rate to become infectious from latent class
  
  alpha=params["alpha"] #proportion of fleas infected from host
  muf=params["muf"] #natural mortality of flea
  mupb=params["mupb"] #mortality of partially blocked flea
  mub=params["mub"] #mortality of blocked flea
  
  lambdaA=params["lambdaA"] #rate of developing partial blockage from EP
  lambdaB=params["lambdaB"] #rate of clearing infection in EP (back to uninfected)
  lambdaC=params["lambdaC"] #rate of leaving EP (still infected but not enough to block)
  tau=params["tau"] #rate of developing full blockage
  
  
  S=x[1]
  L=x[2]
  I=x[3]
  E=x[4]
  R=x[5]
  dr=x[6]
  Id=x[7]
  
  U=x[8]
  Iep=x[9]
  Ipb=x[10]
  Ib=x[11]
  df=x[12]
  
  dS=-S*(b*Iep*pep + b*Ipb*ppb + b1*Ib*pb)-mu*S   #Susceptible Rodent Pop
  dL=S*(b*Iep*pep*tep + b*Ipb*ppb*tpb + b1*Ib*pb*tb) - sigma*L - mu*L  #infection of host, not yet infectious
  dI=sigma*L - I*(mu+epsilon)  #infectious rodent hosts
  dE=S*(b*Iep*pep*(1-tep) + b*Ipb*ppb*(1-tpb) + b1*Ib*pb*(1-tb)) - E*(mu+gamma) #exposed rodents (non-infectious)
  dR=gamma*E - R*mu  #recovery of rodent
  dr=mu*(S+R+E+L+I)+epsilon*I  #dead rodents
  Id=epsilon*I
  
  dU=-I*(U*b*alpha) + Iep*lambdaB - U*muf    #Uninfected Flea pop
  dIep=I*(U*b*alpha) - (Iep*lambdaB + Iep*lambdaC + Iep*lambdaA + Iep*muf) #Infection of early phase fleas 
  dIpb=Iep*lambdaA - Ipb*(tau+mupb) #Infection of partially blocked fleas
  dIb=tau*Ipb - mub*Ib  #Infection of fully blocked fleas
  df=muf*(U+Iep) + mupb*Ipb + mub*Ib + Iep*lambdaC #dead fleas
  
  
  
  list(c(dS,dL,dI,dE,dR,dr,Id, dU,dIep,dIpb,dIb,df))
  
}


times=seq(1,30,by=0.3)					#time steps to output
xstart=c(S=99,L=0,I=1,E=0,R=0,dr=0,Id=0, U=500,Iep=0,Ipb=0,Ib=0,df=0)					#beginning population size
parms=c(alpha=1,lambdaA=0.04,lambdaB=0.02,lambdaC=0.06,b=0.4,b1=1,tau=0.48, muf=0.02, mupb=0.13,mub=0.26,
        pep=0.10,ppb=0.10,pb=0.67,tep=0.5, tpb=1 , tb=0.80 ,
        mu=0.002,sigma=0.25,gamma=0.14,epsilon=0.5)	#set parameter values
out=lsoda(xstart, times, SIR.model, parms)  #run the model

time=out[,1]
S=out[,2]
L=out[,3]
I=out[,4]	
E=out[,5]
R=out[,6]
dr=out[,7]
Id=out[,8]

U=out[,9]
Iep=out[,10]
Ipb=out[,11]
Ib=out[,12]
df=out[,13]

#rat_SIR<-cbind.data.frame(time, S, L, I, E, R, dr, U, Iep, Ipb, Ib, df)


#..................................................................
#Graph the model output####
#..................................................................



#HOST Compartment;
#par(mfrow=c(2,1))
#postscript("C:/Users/Cedar/Desktop/Plague/Models/Graphs/Host_mouseblood_10CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)
postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Host_mouseblood_10CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(time,I1,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(100)),lty=1)
lines(time,E1,col="blue",lwd=2,lty=1)
lines(time,S1,col="forestgreen",lwd=2,lty=1)
lines(time,L1,col="red",lwd=2,lty=5)
lines(time,R1,col="blue",lwd=2,lty=5)
lines(time,Id1,col="black",lwd=2,lty=1)
legend(23,50,c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1,5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp =1)
title(main="Host Dynamics with Mouse Blood Infected Fleas, 99UH1IH 500UF LD50 of 10 CFU")
dev.off()
#rat plot;
#postscript("C:/Users/Cedar/Desktop/Plague/Models/Graphs/Host_ratblood_10CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)
postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Host_ratblood_10CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(time,I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(100)),lty=1)
lines(time,E,col="blue",lwd=2,lty=1)
lines(time,S,col="forestgreen",lwd=2,lty=1)
lines(time,L,col="red",lwd=2,lty=5)
lines(time,R,col="blue",lwd=2,lty=5)
lines(time,Id,col="black",lwd=2,lty=1)
legend(23,50,c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1,5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Host Dynamics with Rat Blood Infected Fleas, 99UH1IH 500UF LD50 of 10 CFU")
dev.off()

#............................................................................
#FLEA Compartment

#par(mfrow=c(2,1))
#postscript("C:/Users/Cedar/Desktop/Plague/Models/Graphs/Flea_mouseblood_10CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)
postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Flea_mouseblood_10CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(time,U1,ylab="Abundance",xlab="Time",type="l",col="forestgreen",lwd=2,ylim=c(0,max(500)),lty=1)
lines(time,Iep1,col="blue",lwd=2,lty=1)
lines(time,Ipb1,col="red",lwd=2,lty=5)
lines(time,Ib1,col="red",lwd=2,lty=1)
lines(time,df1,col="black",lwd=2,lty=1)
legend(20,100,c("Uninfected", "Infectious-Early Phase", "Infectious-Partially Blocked","Infectious-Fully Blocked", "Cumulative Dead"),
       col=c("forestgreen","blue", "red","red", "black"),bty="n",lty=c(1,1, 5, 1, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Flea Dynamics with Mouse Blood Infected Fleas, 99UH1IH 500UF LD50 of 10 CFU")
dev.off()
#rat plot;
#postscript("C:/Users/Cedar/Desktop/Plague/Models/Graphs/Flea_ratblood_10CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)
postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Flea_ratblood_10CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(time,U,ylab="Abundance",xlab="Time",type="l",col="forestgreen",lwd=2,ylim=c(0,max(500)),lty=1)
lines(time,Iep,col="blue",lwd=2,lty=1)
lines(time,Ipb,col="red",lwd=2,lty=5)
lines(time,Ib,col="red",lwd=2,lty=1)
lines(time,df,col="black",lwd=2,lty=1)
legend(20,100,c("Uninfected", "Infectious-Early Phase", "Infectious-Partially Blocked","Infectious-Fully Blocked", "Cumulative Dead"),
       col=c("forestgreen","blue", "red","red", "black"),bty="n",lty=c(1,1, 5, 1, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Flea Dynamics with Rat Blood Infected Fleas, 99UH1IH 500UF LD50 of 10 CFU")
dev.off()

#..................................................
#Threshold of 100 CFU####
#..................................................
SIR.model.EPT=function(t,x,params){
  b=params["b"] #biting rate of fleas
  b1=params["b1"]  #biting rate of blocked fleas
  tep=params["tep"]  #proprotion of INFECTIOUS bites from early-phase fleas
  tpb=params["tpb"]  #proprotion of INFECTIOUS bites from partially blocked fleas
  tb=params["tb"]  #proprotion of INFECTIOUS bites from fully blocked fleas
  mu=params["mu"]		#natural mortality rate of rodent (uninfected)
  pep=params["pep"] #Proportion of ep fleas that transmit
  ppb=params["ppb"]  #Proportion of partially blocked fleas that transmit
  pb=params["pb"]  #Proportion of fully blocked fleas that transmit
  gamma=params["gamma"]		#recovery rate in rodent from low-dose flea infection
  epsilon=params["epsilon"]  #disease induced mortality rate in rodent from high-dose flea infection
  sigma=params["sigma"]  #rate to become infectious from latent class
  
  alpha=params["alpha"] #proportion of fleas infected from host
  muf=params["muf"] #natural mortality of flea
  mupb=params["mupb"] #mortality of partially blocked flea
  mub=params["mub"] #mortality of blocked flea
  
  lambdaA=params["lambdaA"] #rate of developing partial blockage from EP
  lambdaB=params["lambdaB"] #rate of clearing infection in EP (back to uninfected)
  lambdaC=params["lambdaC"] #rate of leaving EP (still infected but not enough to block)
  tau=params["tau"] #rate of developing full blockage
  
  S1=x[1]
  L1=x[2]
  I1=x[3]
  E1=x[4]
  R1=x[5]
  dr1=x[6]
  Id1=x[7]
  
  U1=x[8]
  Iep1=x[9]
  Ipb1=x[10]
  Ib1=x[11]
  df1=x[12]
  
  
  dS1=-S1*(b*Iep1*pep + b*Ipb1*ppb + b1*Ib1*pb) - mu*S1   #Susceptible Rodent Pop
  dL1=S1*(b*Iep1*pep*tep + b*Ipb1*ppb*tpb + b1*Ib1*pb*tb) - sigma*L1 - mu*L1  #infection of host, not yet infectious
  dI1=sigma*L1-I1*(mu+epsilon)  #infectious rodent hosts
  dE1=S1*(b*Iep1*pep*(1-tep) + b*Ipb1*ppb*(1-tpb) + b1*Ib1*pb*(1-tb)) - E1*(mu+gamma) #exposed rodents (non-infectious)
  dR1=gamma*E1 - R1*mu  #recovery of rodent
  dr1=mu*(S1+R1+E1+L1+I1) + epsilon*I1  #dead rodents
  Id1=epsilon*I1

  dU1=-I1*(U1*b*alpha) + Iep1*lambdaB - U1*muf    #Uninfected Flea pop
  dIep1=I1*(U1*b*alpha) - (Iep1*lambdaB + Iep1*lambdaC + Iep1*lambdaA + Iep1*muf) #Infection of early phase fleas 
  dIpb1=Iep1*lambdaA - Ipb1*(tau+mupb) #Infection of partially blocked fleas
  dIb1=tau*Ipb1 - mub*Ib1  #Infection of fully blocked fleas
  df1=muf*(U1+Iep1)+ mupb*Ipb1 + mub*Ib1 + Iep1*lambdaC #dead fleas



list(c(dS1,dL1,dI1,dE1,dR1,dr1,Id1, dU1,dIep1,dIpb1,dIb1,df1))

}


times=seq(1,30,by=0.3)					#time steps to output
xstart=c(S1=99,L1=0,I1=1,E1=0,R1=0,dr1=0,Id1=0, U1=500,Iep1=0,Ipb1=0,Ib1=0,df1=0)					#beginning population size
parms=c(alpha=1,lambdaA=0.04,lambdaB=0.20,lambdaC=0.07,b=0.4,b1=1,tau=0.39, muf=0.02, mupb=0.14,mub=0.20,
        pep=0.05,ppb=0.11,pb=0.50, tep=0, tpb=0.5, tb=0.65, 
        mu=0.002,sigma=0.25,gamma=0.14,epsilon=0.5)	#set parameter values
out=lsoda(xstart, times, SIR.model.EPT, parms)  #run the model

time=out[,1]
S1=out[,2]
L1=out[,3]
I1=out[,4]	
E1=out[,5]
R1=out[,6]
dr1=out[,7]
Id1=out[,8]

U1=out[,9]
Iep1=out[,10]
Ipb1=out[,11]
Ib1=out[,12]
df1=out[,13]

#Verify that all fleas at end add up to starting number;
sum(df1[97]+Ib1[97]+Ipb1[97]+U1[97]+Iep1[97])
#Verify that all rodents add up to starting number at end of outbreak;
sum(S1[97], L1[97], I1[97], E1[97], R1[97], dr1[97])

#Create dataframe of infection output by time for each model state;
#mouse_SIR<-cbind.data.frame(time, S1, L1, I1, E1, R1, dr1, U1, Iep1, Ipb1, Ib1, df1)

#................................................
#Infected and Maintained on Rat Blood####
#................................................

SIR.model=function(t,x,params){
  b=params["b"] #biting rate of fleas
  b1=params["b1"]  #biting rate of blocked fleas
  tep=params["tep"]  #proprotion of INFECTIOUS bites from early-phase fleas
  tpb=params["tpb"]  #proprotion of INFECTIOUS bites from partially blocked fleas
  tb=params["tb"]  #proprotion of INFECTIOUS bites from fully blocked fleas
  mu=params["mu"]		#natural mortality rate of rodent (uninfected)
  pep=params["pep"] #Proportion of ep fleas that transmit
  ppb=params["ppb"]  #Proportion of partially blocked fleas that transmit
  pb=params["pb"]  #Proportion of fully blocked fleas that transmit
  gamma=params["gamma"]		#recovery rate in rodent from low-dose flea infection
  epsilon=params["epsilon"]  #disease induced mortality rate in rodent from high-dose flea infection
  sigma=params["sigma"]  #rate to become infectious from latent class
  
  alpha=params["alpha"] #proportion of fleas infected from host
  muf=params["muf"] #natural mortality of flea
  mupb=params["mupb"] #mortality of partially blocked flea
  mub=params["mub"] #mortality of blocked flea
  
  lambdaA=params["lambdaA"] #rate of developing partial blockage from EP
  lambdaB=params["lambdaB"] #rate of clearing infection in EP (back to uninfected)
  lambdaC=params["lambdaC"] #rate of leaving EP (still infected but not enough to block)
  tau=params["tau"] #rate of developing full blockage
  
  S=x[1]
  L=x[2]
  I=x[3]
  E=x[4]
  R=x[5]
  dr=x[6]
  Id=x[7]
  
  U=x[8]
  Iep=x[9]
  Ipb=x[10]
  Ib=x[11]
  df=x[12]
  
  dS=-S*(b*Iep*pep + b*Ipb*ppb + b1*Ib*pb) - mu*S   #Susceptible Rodent Pop
  dL=S*(b*Iep*pep*tep + b*Ipb*ppb*tpb + b1*Ib*pb*tb) - sigma*L - mu*L  #infection of host, not yet infectious
  dI=sigma*L - I*(mu+epsilon)  #infectious rodent hosts
  dE=S*(b*Iep*pep*(1-tep) + b*Ipb*ppb*(1-tpb) + b1*Ib*pb*(1-tb)) - E*(mu+gamma) #exposed rodents (non-infectious)
  dR=gamma*E - R*mu  #recovery of rodent
  dr=mu*(S+R+E+L+I)+epsilon*I  #dead rodents
  Id=epsilon*I
  
  dU=-I*(U*b*alpha) + Iep*lambdaB - U*muf    #Uninfected Flea pop
  dIep=I*(U*b*alpha) - (Iep*lambdaB + Iep*lambdaC + Iep*lambdaA +Iep*muf) #Infection of early phase fleas 
  dIpb=Iep*lambdaA - Ipb*(tau+mupb) #Infection of partially blocked fleas
  dIb=tau*Ipb - mub*Ib  #Infection of fully blocked fleas
  df=muf*(U+Iep) + mupb*Ipb + mub*Ib + Iep*lambdaC #dead fleas
  
  
  list(c(dS,dL,dI,dE,dR,dr,Id, dU,dIep,dIpb,dIb,df))
  
}


times=seq(1,30,by=0.3)					#time steps to output
xstart=c(S=99,L=0,I=1,E=0,R=0,dr=0,Id=0, U=500,Iep=0,Ipb=0,Ib=0,df=0)					#beginning population size
parms=c(alpha=1,lambdaA=0.04,lambdaB=0.02,lambdaC=0.06,b=0.4,b1=1,tau=0.48, muf=0.02, mupb=0.13,mub=0.26,
        pep=0.10,ppb=0.10,pb=0.67,tep=0, tpb=1 , tb=0.41 ,
        mu=0.002,sigma=0.25,gamma=0.14,epsilon=0.5)	#set parameter values
out=lsoda(xstart, times, SIR.model, parms)  #run the model

time=out[,1]
S=out[,2]
L=out[,3]
I=out[,4]	
E=out[,5]
R=out[,6]
dr=out[,7]
Id=out[,8]

U=out[,9]
Iep=out[,10]
Ipb=out[,11]
Ib=out[,12]
df=out[,13]

#rat_SIR<-cbind.data.frame(time, S, L, I, E, R, dr, U, Iep, Ipb, Ib, df)


#..................................................................
#Graph the model output####
#..................................................................

#HOST Compartment;
#par(mfrow=c(2,1))
#postscript("C:/Users/Cedar/Desktop/Plague/Models/Graphs/Host_mouseblood_100CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)
postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Host_mouseblood_100CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(time,I1,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(100)),lty=1)
lines(time,E1,col="blue",lwd=2,lty=1)
lines(time,S1,col="forestgreen",lwd=2,lty=1)
lines(time,L1,col="red",lwd=2,lty=5)
lines(time,R1,col="blue",lwd=2,lty=5)
lines(time,Id1,col="black",lwd=2,lty=1)
legend(23,50,c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1,5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp =1)
title(main="Host Dynamics with Mouse Blood Infected Fleas, 99UH1IH 500UF LD50 of 100 CFU")
dev.off()
#rat plot;
#postscript("C:/Users/Cedar/Desktop/Plague/Models/Graphs/Host_ratblood_100CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)
postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Host_ratblood_100CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(time,I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(100)),lty=1)
lines(time,E,col="blue",lwd=2,lty=1)
lines(time,S,col="forestgreen",lwd=2,lty=1)
lines(time,L,col="red",lwd=2,lty=5)
lines(time,R,col="blue",lwd=2,lty=5)
lines(time,Id,col="black",lwd=2,lty=1)
legend(23,50,c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1,5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Host Dynamics with Rat Blood Infected Fleas, 99UH1IH 500UF LD50 of 100 CFU")
dev.off()

#............................................................................
#FLEA Compartment

#par(mfrow=c(2,1))
#postscript("C:/Users/Cedar/Desktop/Plague/Models/Graphs/Flea_mouseblood_100CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)
postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Flea_mouseblood_100CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(time,U1,ylab="Abundance",xlab="Time",type="l",col="forestgreen",lwd=2,ylim=c(0,max(500)),lty=1)
lines(time,Iep1,col="blue",lwd=2,lty=1)
lines(time,Ipb1,col="red",lwd=2,lty=5)
lines(time,Ib1,col="red",lwd=2,lty=1)
lines(time,df1,col="black",lwd=2,lty=1)
legend(20,100,c("Uninfected", "Infectious-Early Phase", "Infectious-Partially Blocked","Infectious-Fully Blocked", "Cumulative Dead"),
       col=c("forestgreen","blue", "red","red", "black"),bty="n",lty=c(1,1, 5, 1, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Flea Dynamics with Mouse Blood Infected Fleas, 99UH1IH 500UF LD50 of 100 CFU")
dev.off()
#rat plot;
#postscript("C:/Users/Cedar/Desktop/Plague/Models/Graphs/Flea_ratblood_100CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)
postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Flea_ratblood_100CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(time,U,ylab="Abundance",xlab="Time",type="l",col="forestgreen",lwd=2,ylim=c(0,max(500)),lty=1)
lines(time,Iep,col="blue",lwd=2,lty=1)
lines(time,Ipb,col="red",lwd=2,lty=5)
lines(time,Ib,col="red",lwd=2,lty=1)
lines(time,df,col="black",lwd=2,lty=1)
legend(20,100,c("Uninfected", "Infectious-Early Phase", "Infectious-Partially Blocked","Infectious-Fully Blocked", "Cumulative Dead"),
       col=c("forestgreen","blue", "red","red", "black"),bty="n",lty=c(1,1, 5, 1, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Flea Dynamics with Rat Blood Infected Fleas, 99UH1IH 500UF LD50 of 100 CFU")
dev.off()
