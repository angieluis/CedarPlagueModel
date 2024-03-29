library(deSolve)
library(tidyverse)
#source("/Users/jhinnebusch/R/Plague/Models/SIRmodel.R")
source("SIRmodel.R")

#### This is using the fleas per host model

#### These simulations explore if there was no transmission from fleas in the early phase of infection. Only transmission from partially blocked and blocked fleas. (set pep=0)

#................................................
# Infected and Maintained on Mouse Blood####
# Assuming lethal dose = 1 CFU
#................................................

times=seq(1,100,by=0.3)					#time steps to output
xstart=c(S=9,L=0,I=1,E=0,R=0,dr=0,Id=0, U=100,Iep=0,Ipb=0,Ib=0,df=0)					#beginning population size
parms=c(alpha=1,lambdaA=0.035,lambdaB=0.20,lambdaC=0.07,b=0.4,b1=2,tau=0.39, muf=0.02, mupb=0.14,mub=0.20,
        pep=0,ppb=0.11,pb=0.5, tep=1, tpb=1, tb=1, 
        mu=0.002,sigma=0.25,gamma=0.14,epsilon=0.5)	#set parameter values
out=lsoda(xstart, times, SIR.model.fleasperhost, parms)  #run the model

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
BPBonly_mouse_1CFU<-cbind.data.frame(time, S=S1, L=L1, I=I1, E=E1, R=R1, Id=Id1, dr=dr1, U=U1, Iep=Iep1, Ipb=Ipb1, Ib=Ib1, df=df1)




#................................................
# Infected and Maintained on Rat Blood####
# Assuming lethal dose = 1 CFU 
#................................................


times=seq(1,100,by=0.3)					#time steps to output
xstart=c(S=9,L=0,I=1,E=0,R=0,dr=0,Id=0, U=100,Iep=0,Ipb=0,Ib=0,df=0)					#beginning population size
parms=c(alpha=1,lambdaA=0.04,lambdaB=0.02,lambdaC=0.06,b=0.4,b1=2,tau=0.48, muf=0.02, mupb=0.13,mub=0.26,
        pep=0,ppb=0.10,pb=0.67,tep=1, tpb=1 , tb=1 ,
        mu=0.002,sigma=0.25,gamma=0.14,epsilon=0.5)	#set parameter values
out=lsoda(xstart, times, SIR.model.fleasperhost, parms)  #run the model

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

BPBonly_rat_1CFU<-cbind.data.frame(time, S, L, I, E, R, dr, Id, U, Iep, Ipb, Ib, df)


#..................................................................
#Graph the model output####
#..................................................................



#HOST Compartment;
#par(mfrow=c(2,1))
#postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Host_mouseblood_1CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

#mouse plot
plot(BPBonly_mouse_1CFU$time,BPBonly_mouse_1CFU$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(BPBonly_mouse_1CFU$time,BPBonly_mouse_1CFU$E,col="blue",lwd=2,lty=1)
lines(BPBonly_mouse_1CFU$time,BPBonly_mouse_1CFU$S,col="forestgreen",lwd=2,lty=1)
lines(BPBonly_mouse_1CFU$time,BPBonly_mouse_1CFU$L,col="red",lwd=2,lty=5)
lines(BPBonly_mouse_1CFU$time,BPBonly_mouse_1CFU$R,col="blue",lwd=2,lty=5)
lines(BPBonly_mouse_1CFU$time,BPBonly_mouse_1CFU$Id,col="black",lwd=2,lty=1)
legend(23,10,c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp =0.5, y.intersp =1)
title(main="Mouse Blood, no early phase, 1 CFU")
#dev.off()

#rat plot;
#postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Host_ratblood_1CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(BPBonly_rat_1CFU$time,BPBonly_rat_1CFU$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(BPBonly_rat_1CFU$time,BPBonly_rat_1CFU$E,col="blue",lwd=2,lty=1)
lines(BPBonly_rat_1CFU$time,BPBonly_rat_1CFU$S,col="forestgreen",lwd=2,lty=1)
lines(BPBonly_rat_1CFU$time,BPBonly_rat_1CFU$L,col="red",lwd=2,lty=5)
lines(BPBonly_rat_1CFU$time,BPBonly_rat_1CFU$R,col="blue",lwd=2,lty=5)
lines(BPBonly_rat_1CFU$time,BPBonly_rat_1CFU$Id,col="black",lwd=2,lty=1)
#legend(23,10,c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Rat Blood, no early phase, 1 CFU")
#dev.off()

#............................................................................
#FLEA Compartment

par(mfrow=c(2,1))
#postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Flea_mouseblood_1CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

#mouse fleas plot
plot(BPBonly_mouse_1CFU$time,BPBonly_mouse_1CFU$U,ylab="Abundance",xlab="Time",type="l",col="forestgreen",lwd=2,ylim=c(0,max(100)),lty=1)
lines(BPBonly_mouse_1CFU$time,BPBonly_mouse_1CFU$Iep,col="blue",lwd=2,lty=1)
lines(BPBonly_mouse_1CFU$time,BPBonly_mouse_1CFU$Ipb,col="red",lwd=2,lty=5)
lines(BPBonly_mouse_1CFU$time,BPBonly_mouse_1CFU$Ib,col="red",lwd=2,lty=1)
lines(BPBonly_mouse_1CFU$time,BPBonly_mouse_1CFU$df,col="black",lwd=2,lty=1)
legend(20,50,c("Uninfected", "Infectious-Early Phase", "Infectious-Partially Blocked","Infectious-Fully Blocked", "Cumulative Dead"),
       col=c("forestgreen","blue", "red","red", "black"),bty="n",lty=c(1,1, 5, 1, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Flea Dynamics with Mouse Blood Infected Fleas, 9UH1IH 100UF LD50 of 1 CFU v6")
#dev.off()

#rat fleas plot;
#postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Flea_ratblood_1CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(BPBonly_rat_1CFU$time,BPBonly_rat_1CFU$U,ylab="Abundance",xlab="Time",type="l",col="forestgreen",lwd=2,ylim=c(0,max(100)),lty=1)
lines(BPBonly_rat_1CFU$time,BPBonly_rat_1CFU$Iep,col="blue",lwd=2,lty=1)
lines(BPBonly_rat_1CFU$time,BPBonly_rat_1CFU$Ipb,col="red",lwd=2,lty=5)
lines(BPBonly_rat_1CFU$time,BPBonly_rat_1CFU$Ib,col="red",lwd=2,lty=1)
lines(BPBonly_rat_1CFU$time,BPBonly_rat_1CFU$df,col="black",lwd=2,lty=1)
legend(20,50,c("Uninfected", "Infectious-Early Phase", "Infectious-Partially Blocked","Infectious-Fully Blocked", "Cumulative Dead"),
       col=c("forestgreen","blue", "red","red", "black"),bty="n",lty=c(1,1, 5, 1, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Flea Dynamics with Rat Blood Infected Fleas, 9UH1IH 100UF LD50 of 1 CFU v6")
#dev.off()





#..................................................
#Threshold of 10 CFU####
# Mouse blood
#..................................................

times=seq(1,100,by=0.3)					#time steps to output
xstart=c(S=9,L=0,I=1,E=0,R=0,dr=0,Id=0, U=100,Iep=0,Ipb=0,Ib=0,df=0)					#beginning population size
parms=c(alpha=1,lambdaA=0.035,lambdaB=0.20,lambdaC=0.07,b=0.4,b1=2,tau=0.39, muf=0.02, mupb=0.14,mub=0.20,
        pep=0,ppb=0.11,pb=0.5, tep=0.001, tpb=0.5, tb=0.8, 
        mu=0.002,sigma=0.25,gamma=0.14,epsilon=0.5)	#set parameter values
out=lsoda(xstart, times, SIR.model.fleasperhost, parms)  #run the model

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
BPBonly_mouse_10CFU<-cbind.data.frame(time, S=S1, L=L1, I=I1, E=E1, R=R1, Id=Id1, dr=dr1, U=U1, Iep=Iep1, Ipb=Ipb1, Ib=Ib1, df=df1)



#................................................
# 10 CFU
# Infected and Maintained on Rat Blood####
#................................................

times=seq(1,100,by=0.3)					#time steps to output
xstart=c(S=9,L=0,I=1,E=0,R=0,dr=0,Id=0, U=100,Iep=0,Ipb=0,Ib=0,df=0)					#beginning population size
parms=c(alpha=1,lambdaA=0.04,lambdaB=0.02,lambdaC=0.06,b=0.4,b1=2,tau=0.48, muf=0.02, mupb=0.13,mub=0.26,
        pep=0,ppb=0.10,pb=0.67,tep=0.5, tpb=1 , tb=0.8 ,
        mu=0.002,sigma=0.25,gamma=0.14,epsilon=0.5)	#set parameter values
out=lsoda(xstart, times, SIR.model.fleasperhost, parms)  #run the model

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

BPBonly_rat_10CFU<-cbind.data.frame(time, S, L, I, E, R, dr, Id, U, Iep, Ipb, Ib, df)


#..................................................................
#Graph the model output####
#..................................................................


par(mfrow=c(2,1))

#mouse plot
plot(BPBonly_mouse_10CFU$time,BPBonly_mouse_10CFU$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(BPBonly_mouse_10CFU$time,BPBonly_mouse_10CFU$E,col="blue",lwd=2,lty=1)
lines(BPBonly_mouse_10CFU$time,BPBonly_mouse_10CFU$S,col="forestgreen",lwd=2,lty=1)
lines(BPBonly_mouse_10CFU$time,BPBonly_mouse_10CFU$L,col="red",lwd=2,lty=5)
lines(BPBonly_mouse_10CFU$time,BPBonly_mouse_10CFU$R,col="blue",lwd=2,lty=5)
lines(BPBonly_mouse_10CFU$time,BPBonly_mouse_10CFU$Id,col="black",lwd=2,lty=1)
#legend(23,10,c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp =0.5, y.intersp =1)
title(main="Mouse Blood, no early phase, 10 CFU")
#dev.off()

#rat plot;
#postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Host_ratblood_10CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(BPBonly_rat_10CFU$time,BPBonly_rat_10CFU$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(BPBonly_rat_10CFU$time,BPBonly_rat_10CFU$E,col="blue",lwd=2,lty=1)
lines(BPBonly_rat_10CFU$time,BPBonly_rat_10CFU$S,col="forestgreen",lwd=2,lty=1)
lines(BPBonly_rat_10CFU$time,BPBonly_rat_10CFU$L,col="red",lwd=2,lty=5)
lines(BPBonly_rat_10CFU$time,BPBonly_rat_10CFU$R,col="blue",lwd=2,lty=5)
lines(BPBonly_rat_10CFU$time,BPBonly_rat_10CFU$Id,col="black",lwd=2,lty=1)
#legend(23,10,c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Rat Blood, no early phase, 10 CFU")
#dev.off()

#............................................................................
#FLEA Compartment

##par(mfrow=c(2,1))
#postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Flea_mouseblood_10CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

#mouse fleas plot
plot(BPBonly_mouse_10CFU$time,BPBonly_mouse_10CFU$U,ylab="Abundance",xlab="Time",type="l",col="forestgreen",lwd=2,ylim=c(0,max(100)),lty=1)
lines(BPBonly_mouse_10CFU$time,BPBonly_mouse_10CFU$Iep,col="blue",lwd=2,lty=1)
lines(BPBonly_mouse_10CFU$time,BPBonly_mouse_10CFU$Ipb,col="red",lwd=2,lty=5)
lines(BPBonly_mouse_10CFU$time,BPBonly_mouse_10CFU$Ib,col="red",lwd=2,lty=1)
lines(BPBonly_mouse_10CFU$time,BPBonly_mouse_10CFU$df,col="black",lwd=2,lty=1)
legend(20,50,c("Uninfected", "Infectious-Early Phase", "Infectious-Partially Blocked","Infectious-Fully Blocked", "Cumulative Dead"),
       col=c("forestgreen","blue", "red","red", "black"),bty="n",lty=c(1,1, 5, 1, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Flea Dynamics with Mouse Blood Infected Fleas, 9UH1IH 100UF LD50 of 10 CFU v6")
#dev.off()

#rat fleas plot;
#postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Flea_ratblood_10CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(BPBonly_rat_10CFU$time,BPBonly_rat_10CFU$U,ylab="Abundance",xlab="Time",type="l",col="forestgreen",lwd=2,ylim=c(0,max(100)),lty=1)
lines(BPBonly_rat_10CFU$time,BPBonly_rat_10CFU$Iep,col="blue",lwd=2,lty=1)
lines(BPBonly_rat_10CFU$time,BPBonly_rat_10CFU$Ipb,col="red",lwd=2,lty=5)
lines(BPBonly_rat_10CFU$time,BPBonly_rat_10CFU$Ib,col="red",lwd=2,lty=1)
lines(BPBonly_rat_10CFU$time,BPBonly_rat_10CFU$df,col="black",lwd=2,lty=1)
legend(20,50,c("Uninfected", "Infectious-Early Phase", "Infectious-Partially Blocked","Infectious-Fully Blocked", "Cumulative Dead"),
       col=c("forestgreen","blue", "red","red", "black"),bty="n",lty=c(1,1, 5, 1, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Flea Dynamics with Rat Blood Infected Fleas, 9UH1IH 100UF LD50 of 10 CFU v6")
#dev.off()





#..................................................
# Threshold of 100 CFU####
# Mouse blood
#..................................................

times=seq(1,100,by=0.3)					#time steps to output
xstart=c(S=9,L=0,I=1,E=0,R=0,dr=0,Id=0, U=100,Iep=0,Ipb=0,Ib=0,df=0)					#beginning population size
parms=c(alpha=1,lambdaA=0.035,lambdaB=0.20,lambdaC=0.07,b=0.4,b1=2,tau=0.39, muf=0.02, mupb=0.14,mub=0.20,
        pep=0,ppb=0.11,pb=0.5, tep=0.00, tpb=0, tb=0.65, 
        mu=0.002,sigma=0.25,gamma=0.14,epsilon=0.5)	#set parameter values
out=lsoda(xstart, times, SIR.model.fleasperhost, parms)  #run the model

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
BPBonly_mouse_100CFU<-cbind.data.frame(time, S=S1, L=L1, I=I1, E=E1, R=R1, Id=Id1, dr=dr1, U=U1, Iep=Iep1, Ipb=Ipb1, Ib=Ib1, df=df1)

#................................................
#Infected and Maintained on Rat Blood####
#................................................


times=seq(1,100,by=0.3)					#time steps to output
xstart=c(S=9,L=0,I=1,E=0,R=0,dr=0,Id=0, U=100,Iep=0,Ipb=0,Ib=0,df=0)					#beginning population size
parms=c(alpha=1,lambdaA=0.04,lambdaB=0.02,lambdaC=0.06,b=0.4,b1=2,tau=0.48, muf=0.02, mupb=0.13,mub=0.26,
        pep=0,ppb=0.10,pb=0.67,tep=0, tpb=1 , tb=0.41 ,
        mu=0.002,sigma=0.25,gamma=0.14,epsilon=0.5)	#set parameter values
out=lsoda(xstart, times, SIR.model.fleasperhost, parms)  #run the model

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

BPBonly_rat_100CFU<-cbind.data.frame(time, S, L, I, E, R, dr, Id, U, Iep, Ipb, Ib, df)


#..................................................................
#Graph the model output####
#..................................................................

#HOST Compartment;
#par(mfrow=c(2,1))
#postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Host_mouseblood_100CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

#mouse
plot(BPBonly_mouse_100CFU$time,BPBonly_mouse_100CFU$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(BPBonly_mouse_100CFU$time,BPBonly_mouse_100CFU$E,col="blue",lwd=2,lty=1)
lines(BPBonly_mouse_100CFU$time,BPBonly_mouse_100CFU$S,col="forestgreen",lwd=2,lty=1)
lines(BPBonly_mouse_100CFU$time,BPBonly_mouse_100CFU$L,col="red",lwd=2,lty=5)
lines(BPBonly_mouse_100CFU$time,BPBonly_mouse_100CFU$R,col="blue",lwd=2,lty=5)
lines(BPBonly_mouse_100CFU$time,BPBonly_mouse_100CFU$Id,col="black",lwd=2,lty=1)
#legend(23,10,c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1,5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp =1)
title(main="Mouse Blood, no early phase, 100 CFU")
#dev.off()

#rat plot;
#postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Host_ratblood_100CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(BPBonly_rat_100CFU$time,BPBonly_rat_100CFU$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(BPBonly_rat_100CFU$time,BPBonly_rat_100CFU$E,col="blue",lwd=2,lty=1)
lines(BPBonly_rat_100CFU$time,BPBonly_rat_100CFU$S,col="forestgreen",lwd=2,lty=1)
lines(BPBonly_rat_100CFU$time,BPBonly_rat_100CFU$L,col="red",lwd=2,lty=5)
lines(BPBonly_rat_100CFU$time,BPBonly_rat_100CFU$R,col="blue",lwd=2,lty=5)
lines(BPBonly_rat_100CFU$time,BPBonly_rat_100CFU$Id,col="black",lwd=2,lty=1)
#legend(23,10,c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1,5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Rat Blood, no early phase, 100 CFU")
#dev.off()

#............................................................................
#FLEA Compartment

#par(mfrow=c(2,1))
#postscript("C:/Users/Cedar/Desktop/Plague/Models/Graphs/Flea_mouseblood_100CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)
postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Flea_mouseblood_100CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(BPBonly_mouse_100CFU$time,BPBonly_mouse_100CFU$U,ylab="Abundance",xlab="Time",type="l",col="forestgreen",lwd=2,ylim=c(0,max(100)),lty=1)
lines(BPBonly_mouse_100CFU$time,BPBonly_mouse_100CFU$Iep,col="blue",lwd=2,lty=1)
lines(BPBonly_mouse_100CFU$time,BPBonly_mouse_100CFU$Ipb,col="red",lwd=2,lty=5)
lines(BPBonly_mouse_100CFU$time,BPBonly_mouse_100CFU$Ib,col="red",lwd=2,lty=1)
lines(BPBonly_mouse_100CFU$time,BPBonly_mouse_100CFU$df,col="black",lwd=2,lty=1)
legend(20,50,c("Uninfected", "Infectious-Early Phase", "Infectious-Partially Blocked","Infectious-Fully Blocked", "Cumulative Dead"),
       col=c("forestgreen","blue", "red","red", "black"),bty="n",lty=c(1,1, 5, 1, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Flea Dynamics with Mouse Blood Infected Fleas, 9UH1IH 100UF LD50 of 100 CFU v6")
dev.off()
#rat plot;
#postscript("C:/Users/Cedar/Desktop/Plague/Models/Graphs/Flea_ratblood_100CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)
postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Flea_ratblood_100CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(BPBonly_rat_100CFU$time,BPBonly_rat_100CFU$U,ylab="Abundance",xlab="Time",type="l",col="forestgreen",lwd=2,ylim=c(0,max(100)),lty=1)
lines(BPBonly_rat_100CFU$time,BPBonly_rat_100CFU$Iep,col="blue",lwd=2,lty=1)
lines(BPBonly_rat_100CFU$time,BPBonly_rat_100CFU$Ipb,col="red",lwd=2,lty=5)
lines(BPBonly_rat_100CFU$time,BPBonly_rat_100CFU$Ib,col="red",lwd=2,lty=1)
lines(BPBonly_rat_100CFU$time,BPBonly_rat_100CFU$df,col="black",lwd=2,lty=1)
legend(20,50,c("Uninfected", "Infectious-Early Phase", "Infectious-Partially Blocked","Infectious-Fully Blocked", "Cumulative Dead"),
       col=c("forestgreen","blue", "red","red", "black"),bty="n",lty=c(1,1, 5, 1, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Flea Dynamics with Rat Blood Infected Fleas, 9UH1IH 100UF LD50 of 100 CFU v6")
dev.off()



###############################################################################
## Compare how many rodents were dead at the end for the different scenarios
###############################################################################

RodentsDead.BPBonly <- data.frame(rodent=c("mouse","rat"),scenario=rep("BPBonly",2),CFU.1=rep(NA,2),CFU.10=rep(NA,2),CFU.100=rep(NA,2))

RodentsDead.BPBonly$CFU.1[which(RodentsDead.BPBonly$rodent=="mouse")] <- BPBonly_mouse_1CFU$dr[length(times)]
RodentsDead.BPBonly$CFU.10[which(RodentsDead.BPBonly$rodent=="mouse")] <- BPBonly_mouse_10CFU$dr[length(times)]
RodentsDead.BPBonly$CFU.100[which(RodentsDead.BPBonly$rodent=="mouse")] <- BPBonly_mouse_100CFU$dr[length(times)]
RodentsDead.BPBonly$CFU.1[which(RodentsDead.BPBonly$rodent=="rat")] <- BPBonly_rat_1CFU$dr[length(times)]
RodentsDead.BPBonly$CFU.10[which(RodentsDead.BPBonly$rodent=="rat")] <- BPBonly_rat_10CFU$dr[length(times)]
RodentsDead.BPBonly$CFU.100[which(RodentsDead.BPBonly$rodent=="rat")] <- BPBonly_rat_100CFU$dr[length(times)]

