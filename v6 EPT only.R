library(deSolve)
library(tidyverse)
#source("/Users/jhinnebusch/R/Plague/Models/SIRmodel.R")
source("SIRmodel.R")

#................................................
# Infected and Maintained on Mouse Blood####
# Assuming lethal dose = 1 CFU
#................................................

times=seq(1,30,by=0.3)					#time steps to output
xstart=c(S=9,L=0,I=1,E=0,R=0,dr=0,Id=0, U=100,Iep=0,Ipb=0,Ib=0,df=0)					#beginning population size
parms=c(alpha=1,lambdaA=0.035,lambdaB=0.20,lambdaC=0.07,b=0.4,b1=2,tau=0.39, muf=0.02, mupb=0.14,mub=0.20,
        pep=0.03,ppb=0,pb=0, tep=1, tpb=1, tb=1, 
        mu=0.002,sigma=0.25,gamma=0.14,epsilon=0.5)	#set parameter values
out=lsoda(xstart, times, SIR.model, parms)  #run the model

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
EPTonly_mouse_1CFU<-cbind.data.frame(time, S=S1, L=L1, I=I1, E=E1, R=R1, Id=Id1, dr=dr1, U=U1, Iep=Iep1, Ipb=Ipb1, Ib=Ib1, df=df1)





#................................................
# Infected and Maintained on Rat Blood####
# Assuming lethal dose = 1 CFU 
#................................................


times=seq(1,30,by=0.3)					#time steps to output
xstart=c(S=9,L=0,I=1,E=0,R=0,dr=0,Id=0, U=100,Iep=0,Ipb=0,Ib=0,df=0)					#beginning population size
parms=c(alpha=1,lambdaA=0.04,lambdaB=0.02,lambdaC=0.06,b=0.4,b1=2,tau=0.48, muf=0.02, mupb=0.13,mub=0.26,
        pep=0.1,ppb=0,pb=0,tep=1, tpb=1 , tb=1 ,
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

EPTonly_rat_1CFU<-cbind.data.frame(time, S, L, I, E, R, dr, Id, U, Iep, Ipb, Ib, df)


#..................................................................
#Graph the model output####
#..................................................................



#HOST Compartment;
#par(mfrow=c(2,1))
#postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Host_mouseblood_1CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

#mouse plot
plot(EPTonly_mouse_1CFU$time,EPTonly_mouse_1CFU$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(EPTonly_mouse_1CFU$time,EPTonly_mouse_1CFU$E,col="blue",lwd=2,lty=1)
lines(EPTonly_mouse_1CFU$time,EPTonly_mouse_1CFU$S,col="forestgreen",lwd=2,lty=1)
lines(EPTonly_mouse_1CFU$time,EPTonly_mouse_1CFU$L,col="red",lwd=2,lty=5)
lines(EPTonly_mouse_1CFU$time,EPTonly_mouse_1CFU$R,col="blue",lwd=2,lty=5)
lines(EPTonly_mouse_1CFU$time,EPTonly_mouse_1CFU$Id,col="black",lwd=2,lty=1)
legend("topright",c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp =0.5, y.intersp =1)
title(main="Mouse Blood, Early Phase only, 1 CFU")
#title(main="Host Dynamics with Mouse Blood Infected Fleas, 9UH1IH 100UF LD50 of 1 CFU  v6")
#dev.off()

#rat plot;
#postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Host_ratblood_1CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(EPTonly_rat_1CFU$time,EPTonly_rat_1CFU$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(EPTonly_rat_1CFU$time,EPTonly_rat_1CFU$E,col="blue",lwd=2,lty=1)
lines(EPTonly_rat_1CFU$time,EPTonly_rat_1CFU$S,col="forestgreen",lwd=2,lty=1)
lines(EPTonly_rat_1CFU$time,EPTonly_rat_1CFU$L,col="red",lwd=2,lty=5)
lines(EPTonly_rat_1CFU$time,EPTonly_rat_1CFU$R,col="blue",lwd=2,lty=5)
lines(EPTonly_rat_1CFU$time,EPTonly_rat_1CFU$Id,col="black",lwd=2,lty=1)
legend("topright",c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Rat Blood, Early Phase only, 1 CFU")
#title(main="Host Dynamics with Rat Blood Infected Fleas, 9UH1IH 100UF LD50 of 1 CFU v6")
#dev.off()

#............................................................................
#FLEA Compartment

par(mfrow=c(2,1))
#postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Flea_mouseblood_1CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

#mouse fleas plot
plot(EPTonly_mouse_1CFU$time,EPTonly_mouse_1CFU$U,ylab="Abundance",xlab="Time",type="l",col="forestgreen",lwd=2,ylim=c(0,max(100)),lty=1)
lines(EPTonly_mouse_1CFU$time,EPTonly_mouse_1CFU$Iep,col="blue",lwd=2,lty=1)
lines(EPTonly_mouse_1CFU$time,EPTonly_mouse_1CFU$Ipb,col="red",lwd=2,lty=5)
lines(EPTonly_mouse_1CFU$time,EPTonly_mouse_1CFU$Ib,col="red",lwd=2,lty=1)
lines(EPTonly_mouse_1CFU$time,EPTonly_mouse_1CFU$df,col="black",lwd=2,lty=1)
legend(20,50,c("Uninfected", "Infectious-Early Phase", "Infectious-Partially Blocked","Infectious-Fully Blocked", "Cumulative Dead"),
       col=c("forestgreen","blue", "red","red", "black"),bty="n",lty=c(1,1, 5, 1, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Flea Dynamics with Mouse Blood Infected Fleas, 9UH1IH 100UF LD50 of 1 CFU v6")
#dev.off()

#rat fleas plot;
#postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Flea_ratblood_1CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(EPTonly_rat_1CFU$time,EPTonly_rat_1CFU$U,ylab="Abundance",xlab="Time",type="l",col="forestgreen",lwd=2,ylim=c(0,max(100)),lty=1)
lines(EPTonly_rat_1CFU$time,EPTonly_rat_1CFU$Iep,col="blue",lwd=2,lty=1)
lines(EPTonly_rat_1CFU$time,EPTonly_rat_1CFU$Ipb,col="red",lwd=2,lty=5)
lines(EPTonly_rat_1CFU$time,EPTonly_rat_1CFU$Ib,col="red",lwd=2,lty=1)
lines(EPTonly_rat_1CFU$time,EPTonly_rat_1CFU$df,col="black",lwd=2,lty=1)
legend(20,50,c("Uninfected", "Infectious-Early Phase", "Infectious-Partially Blocked","Infectious-Fully Blocked", "Cumulative Dead"),
       col=c("forestgreen","blue", "red","red", "black"),bty="n",lty=c(1,1, 5, 1, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Flea Dynamics with Rat Blood Infected Fleas, 9UH1IH 100UF LD50 of 1 CFU v6")
#dev.off()





#..................................................
#Threshold of 10 CFU####
# Mouse blood
#..................................................

times=seq(1,30,by=0.3)					#time steps to output
xstart=c(S=9,L=0,I=1,E=0,R=0,dr=0,Id=0, U=100,Iep=0,Ipb=0,Ib=0,df=0)					#beginning population size
parms=c(alpha=1,lambdaA=0.035,lambdaB=0.20,lambdaC=0.07,b=0.4,b1=2,tau=0.39, muf=0.02, mupb=0.14,mub=0.20,
        pep=0.03,ppb=0,pb=0, tep=0.001, tpb=0.5, tb=0.8, 
        mu=0.002,sigma=0.25,gamma=0.14,epsilon=0.5)	#set parameter values
out=lsoda(xstart, times, SIR.model, parms)  #run the model

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
EPTonly_mouse_10CFU<-cbind.data.frame(time, S=S1, L=L1, I=I1, E=E1, R=R1, Id=Id1, dr=dr1, U=U1, Iep=Iep1, Ipb=Ipb1, Ib=Ib1, df=df1)



#................................................
# 10 CFU
# Infected and Maintained on Rat Blood####
#................................................

times=seq(1,30,by=0.3)					#time steps to output
xstart=c(S=9,L=0,I=1,E=0,R=0,dr=0,Id=0, U=100,Iep=0,Ipb=0,Ib=0,df=0)					#beginning population size
parms=c(alpha=1,lambdaA=0.04,lambdaB=0.02,lambdaC=0.06,b=0.4,b1=2,tau=0.48, muf=0.02, mupb=0.13,mub=0.26,
        pep=0.1,ppb=0,pb=0,tep=0.5, tpb=1 , tb=0.8 ,
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

EPTonly_rat_10CFU<-cbind.data.frame(time, S, L, I, E, R, dr, Id, U, Iep, Ipb, Ib, df)


#..................................................................
#Graph the model output####
#..................................................................



#mouse plot
plot(EPTonly_mouse_10CFU$time,EPTonly_mouse_10CFU$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(EPTonly_mouse_10CFU$time,EPTonly_mouse_10CFU$E,col="blue",lwd=2,lty=1)
lines(EPTonly_mouse_10CFU$time,EPTonly_mouse_10CFU$S,col="forestgreen",lwd=2,lty=1)
lines(EPTonly_mouse_10CFU$time,EPTonly_mouse_10CFU$L,col="red",lwd=2,lty=5)
lines(EPTonly_mouse_10CFU$time,EPTonly_mouse_10CFU$R,col="blue",lwd=2,lty=5)
lines(EPTonly_mouse_10CFU$time,EPTonly_mouse_10CFU$Id,col="black",lwd=2,lty=1)
legend("topright",c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp =0.5, y.intersp =1)
title(main="Mouse Blood, Early Phase only, 10 CFU")
#dev.off()

#rat plot;
#postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Host_ratblood_10CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(EPTonly_rat_10CFU$time,EPTonly_rat_10CFU$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(EPTonly_rat_10CFU$time,EPTonly_rat_10CFU$E,col="blue",lwd=2,lty=1)
lines(EPTonly_rat_10CFU$time,EPTonly_rat_10CFU$S,col="forestgreen",lwd=2,lty=1)
lines(EPTonly_rat_10CFU$time,EPTonly_rat_10CFU$L,col="red",lwd=2,lty=5)
lines(EPTonly_rat_10CFU$time,EPTonly_rat_10CFU$R,col="blue",lwd=2,lty=5)
lines(EPTonly_rat_10CFU$time,EPTonly_rat_10CFU$Id,col="black",lwd=2,lty=1)
legend("topright",c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Rat Blood, Early Phase only, 10 CFU")
#dev.off()

#............................................................................
#FLEA Compartment

##par(mfrow=c(2,1))
#postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Flea_mouseblood_10CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

#mouse fleas plot
plot(EPTonly_mouse_10CFU$time,EPTonly_mouse_10CFU$U,ylab="Abundance",xlab="Time",type="l",col="forestgreen",lwd=2,ylim=c(0,max(100)),lty=1)
lines(EPTonly_mouse_10CFU$time,EPTonly_mouse_10CFU$Iep,col="blue",lwd=2,lty=1)
lines(EPTonly_mouse_10CFU$time,EPTonly_mouse_10CFU$Ipb,col="red",lwd=2,lty=5)
lines(EPTonly_mouse_10CFU$time,EPTonly_mouse_10CFU$Ib,col="red",lwd=2,lty=1)
lines(EPTonly_mouse_10CFU$time,EPTonly_mouse_10CFU$df,col="black",lwd=2,lty=1)
legend(20,50,c("Uninfected", "Infectious-Early Phase", "Infectious-Partially Blocked","Infectious-Fully Blocked", "Cumulative Dead"),
       col=c("forestgreen","blue", "red","red", "black"),bty="n",lty=c(1,1, 5, 1, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Flea Dynamics with Mouse Blood Infected Fleas, 9UH1IH 100UF LD50 of 10 CFU v6")
#dev.off()

#rat fleas plot;
#postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Flea_ratblood_10CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(EPTonly_rat_10CFU$time,EPTonly_rat_10CFU$U,ylab="Abundance",xlab="Time",type="l",col="forestgreen",lwd=2,ylim=c(0,max(100)),lty=1)
lines(EPTonly_rat_10CFU$time,EPTonly_rat_10CFU$Iep,col="blue",lwd=2,lty=1)
lines(EPTonly_rat_10CFU$time,EPTonly_rat_10CFU$Ipb,col="red",lwd=2,lty=5)
lines(EPTonly_rat_10CFU$time,EPTonly_rat_10CFU$Ib,col="red",lwd=2,lty=1)
lines(EPTonly_rat_10CFU$time,EPTonly_rat_10CFU$df,col="black",lwd=2,lty=1)
legend(20,50,c("Uninfected", "Infectious-Early Phase", "Infectious-Partially Blocked","Infectious-Fully Blocked", "Cumulative Dead"),
       col=c("forestgreen","blue", "red","red", "black"),bty="n",lty=c(1,1, 5, 1, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Flea Dynamics with Rat Blood Infected Fleas, 9UH1IH 100UF LD50 of 10 CFU v6")
#dev.off()





#..................................................
# Threshold of 100 CFU####
# Mouse blood
#..................................................

times=seq(1,30,by=0.3)					#time steps to output
xstart=c(S=9,L=0,I=1,E=0,R=0,dr=0,Id=0, U=100,Iep=0,Ipb=0,Ib=0,df=0)					#beginning population size
parms=c(alpha=1,lambdaA=0.035,lambdaB=0.20,lambdaC=0.07,b=0.4,b1=2,tau=0.39, muf=0.02, mupb=0.14,mub=0.20,
        pep=0.03,ppb=0,pb=0, tep=0.00, tpb=0, tb=0.65, 
        mu=0.002,sigma=0.25,gamma=0.14,epsilon=0.5)	#set parameter values
out=lsoda(xstart, times, SIR.model, parms)  #run the model

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
EPTonly_mouse_100CFU<-cbind.data.frame(time, S=S1, L=L1, I=I1, E=E1, R=R1, Id=Id1, dr=dr1, U=U1, Iep=Iep1, Ipb=Ipb1, Ib=Ib1, df=df1)

#................................................
#Infected and Maintained on Rat Blood####
#................................................


times=seq(1,30,by=0.3)					#time steps to output
xstart=c(S=9,L=0,I=1,E=0,R=0,dr=0,Id=0, U=100,Iep=0,Ipb=0,Ib=0,df=0)					#beginning population size
parms=c(alpha=1,lambdaA=0.04,lambdaB=0.02,lambdaC=0.06,b=0.4,b1=2,tau=0.48, muf=0.02, mupb=0.13,mub=0.26,
        pep=0.1,ppb=0,pb=0,tep=0, tpb=1 , tb=0.41 ,
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

EPTonly_rat_100CFU<-cbind.data.frame(time, S, L, I, E, R, dr, Id, U, Iep, Ipb, Ib, df)


#..................................................................
#Graph the model output####
#..................................................................

#HOST Compartment;
#par(mfrow=c(2,1))
#postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Host_mouseblood_100CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

#mouse
plot(EPTonly_mouse_100CFU$time,EPTonly_mouse_100CFU$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(EPTonly_mouse_100CFU$time,EPTonly_mouse_100CFU$E,col="blue",lwd=2,lty=1)
lines(EPTonly_mouse_100CFU$time,EPTonly_mouse_100CFU$S,col="forestgreen",lwd=2,lty=1)
lines(EPTonly_mouse_100CFU$time,EPTonly_mouse_100CFU$L,col="red",lwd=2,lty=5)
lines(EPTonly_mouse_100CFU$time,EPTonly_mouse_100CFU$R,col="blue",lwd=2,lty=5)
lines(EPTonly_mouse_100CFU$time,EPTonly_mouse_100CFU$Id,col="black",lwd=2,lty=1)
legend("topright",c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1,5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp =1)
title(main="Mouse Blood, Early Phase only, 100 CFU")
#dev.off()

#rat plot;
#postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Host_ratblood_100CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(EPTonly_rat_100CFU$time,EPTonly_rat_100CFU$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(EPTonly_rat_100CFU$time,EPTonly_rat_100CFU$E,col="blue",lwd=2,lty=1)
lines(EPTonly_rat_100CFU$time,EPTonly_rat_100CFU$S,col="forestgreen",lwd=2,lty=1)
lines(EPTonly_rat_100CFU$time,EPTonly_rat_100CFU$L,col="red",lwd=2,lty=5)
lines(EPTonly_rat_100CFU$time,EPTonly_rat_100CFU$R,col="blue",lwd=2,lty=5)
lines(EPTonly_rat_100CFU$time,EPTonly_rat_100CFU$Id,col="black",lwd=2,lty=1)
legend(23,10,c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1,5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Rat Blood,Early Phase only, 100 CFU")
#dev.off()

#............................................................................
#FLEA Compartment

#par(mfrow=c(2,1))
#postscript("C:/Users/Cedar/Desktop/Plague/Models/Graphs/Flea_mouseblood_100CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)
postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Flea_mouseblood_100CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(EPTonly_mouse_100CFU$time,EPTonly_mouse_100CFU$U,ylab="Abundance",xlab="Time",type="l",col="forestgreen",lwd=2,ylim=c(0,max(100)),lty=1)
lines(EPTonly_mouse_100CFU$time,EPTonly_mouse_100CFU$Iep,col="blue",lwd=2,lty=1)
lines(EPTonly_mouse_100CFU$time,EPTonly_mouse_100CFU$Ipb,col="red",lwd=2,lty=5)
lines(EPTonly_mouse_100CFU$time,EPTonly_mouse_100CFU$Ib,col="red",lwd=2,lty=1)
lines(EPTonly_mouse_100CFU$time,EPTonly_mouse_100CFU$df,col="black",lwd=2,lty=1)
legend(20,50,c("Uninfected", "Infectious-Early Phase", "Infectious-Partially Blocked","Infectious-Fully Blocked", "Cumulative Dead"),
       col=c("forestgreen","blue", "red","red", "black"),bty="n",lty=c(1,1, 5, 1, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Flea Dynamics with Mouse Blood, Early phase only, LD50 of 100 CFU v6")
dev.off()
#rat plot;
#postscript("C:/Users/Cedar/Desktop/Plague/Models/Graphs/Flea_ratblood_100CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)
postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Flea_ratblood_100CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(EPTonly_rat_100CFU$time,EPTonly_rat_100CFU$U,ylab="Abundance",xlab="Time",type="l",col="forestgreen",lwd=2,ylim=c(0,max(100)),lty=1)
lines(EPTonly_rat_100CFU$time,EPTonly_rat_100CFU$Iep,col="blue",lwd=2,lty=1)
lines(EPTonly_rat_100CFU$time,EPTonly_rat_100CFU$Ipb,col="red",lwd=2,lty=5)
lines(EPTonly_rat_100CFU$time,EPTonly_rat_100CFU$Ib,col="red",lwd=2,lty=1)
lines(EPTonly_rat_100CFU$time,EPTonly_rat_100CFU$df,col="black",lwd=2,lty=1)
legend(20,50,c("Uninfected", "Infectious-Early Phase", "Infectious-Partially Blocked","Infectious-Fully Blocked", "Cumulative Dead"),
       col=c("forestgreen","blue", "red","red", "black"),bty="n",lty=c(1,1, 5, 1, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Flea Dynamics with Rat Blood Infected Fleas, Early Phase only, of 100 CFU v6")
dev.off()
