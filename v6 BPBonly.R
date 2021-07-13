library(deSolve)
library(tidyverse)
source("/Users/jhinnebusch/R/Plague/Models/SIRmodel.R")



times=seq(1,30,by=0.3)					#time steps to output
xstart=c(S=9,L=0,I=1,E=0,R=0,dr=0,Id=0, U=100,Iep=0,Ipb=0,Ib=0,df=0)					#beginning population size
parms=c(alpha=1,lambdaA=0.04,lambdaB=0.20,lambdaC=0.07,b=0.4,b1=2,tau=0.39, muf=0.02, mupb=0.14,mub=0.20,
        pep=0,ppb=0.11,pb=0.75, tep=0, tpb=1, tb=1, 
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
#mouse_SIR<-cbind.data.frame(time, S1, L1, I1, E1, R1, dr1, U1, Iep1, Ipb1, Ib1, df1)

#................................................
#Infected and Maintained on Rat Blood####
#................................................


times=seq(1,30,by=0.3)					#time steps to output
xstart=c(S=9,L=0,I=1,E=0,R=0,dr=0,Id=0, U=100,Iep=0,Ipb=0,Ib=0,df=0)					#beginning population size
parms=c(alpha=1,lambdaA=0.04,lambdaB=0.02,lambdaC=0.06,b=0.4,b1=2,tau=0.48, muf=0.02, mupb=0.13,mub=0.26,
        pep=0,ppb=0.10,pb=0.89,tep=0, tpb=1 , tb=1 ,
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

plot(time,I1,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(time,E1,col="blue",lwd=2,lty=1)
lines(time,S1,col="forestgreen",lwd=2,lty=1)
lines(time,L1,col="red",lwd=2,lty=5)
lines(time,R1,col="blue",lwd=2,lty=5)
lines(time,Id1,col="black",lwd=2,lty=1)
legend(23,10,c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp =0.5, y.intersp =1)
title(main="Host Dynamics with Mouse Blood Infected Fleas, 9UH1IH 100UF LD50 of 1 CFU  v6 BPB only")
dev.off()
#rat plot;
#postscript("C:/Users/Cedar/Desktop/Plague/Models/Graphs/Host_ratblood_1CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)
postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Host_ratblood_1CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(time,I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(time,E,col="blue",lwd=2,lty=1)
lines(time,S,col="forestgreen",lwd=2,lty=1)
lines(time,L,col="red",lwd=2,lty=5)
lines(time,R,col="blue",lwd=2,lty=5)
lines(time,Id,col="black",lwd=2,lty=1)
legend(23,10,c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Host Dynamics with Rat Blood Infected Fleas, 9UH1IH 100UF LD50 of 1 CFU v6 BPB only")
dev.off()

#............................................................................
#FLEA Compartment

#par(mfrow=c(2,1))
#postscript("C:/Users/Cedar/Desktop/Plague/Models/Graphs/Flea_mouseblood_1CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)
postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Flea_mouseblood_1CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(time,U1,ylab="Abundance",xlab="Time",type="l",col="forestgreen",lwd=2,ylim=c(0,max(100)),lty=1)
lines(time,Iep1,col="blue",lwd=2,lty=1)
lines(time,Ipb1,col="red",lwd=2,lty=5)
lines(time,Ib1,col="red",lwd=2,lty=1)
lines(time,df1,col="black",lwd=2,lty=1)
legend(20,50,c("Uninfected", "Infectious-Early Phase", "Infectious-Partially Blocked","Infectious-Fully Blocked", "Cumulative Dead"),
       col=c("forestgreen","blue", "red","red", "black"),bty="n",lty=c(1,1, 5, 1, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Flea Dynamics with Mouse Blood Infected Fleas, 9UH1IH 100UF LD50 of 1 CFU v6 BPB only")
dev.off()
#rat plot;
#postscript("C:/Users/Cedar/Desktop/Plague/Models/Graphs/Flea_ratblood_1CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)
postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Flea_ratblood_1CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(time,U,ylab="Abundance",xlab="Time",type="l",col="forestgreen",lwd=2,ylim=c(0,max(100)),lty=1)
lines(time,Iep,col="blue",lwd=2,lty=1)
lines(time,Ipb,col="red",lwd=2,lty=5)
lines(time,Ib,col="red",lwd=2,lty=1)
lines(time,df,col="black",lwd=2,lty=1)
legend(20,50,c("Uninfected", "Infectious-Early Phase", "Infectious-Partially Blocked","Infectious-Fully Blocked", "Cumulative Dead"),
       col=c("forestgreen","blue", "red","red", "black"),bty="n",lty=c(1,1, 5, 1, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Flea Dynamics with Rat Blood Infected Fleas, 9UH1IH 100UF LD50 of 1 CFU v6 BPB only")
dev.off()



#..................................................
#Threshold of 10 CFU####
#..................................................

times=seq(1,30,by=0.3)					#time steps to output
xstart=c(S=9,L=0,I=1,E=0,R=0,dr=0,Id=0, U=100,Iep=0,Ipb=0,Ib=0,df=0)					#beginning population size
parms=c(alpha=1,lambdaA=0.04,lambdaB=0.20,lambdaC=0.07,b=0.4,b1=2,tau=0.39, muf=0.02, mupb=0.14,mub=0.20,
        pep=0,ppb=0.11,pb=0.75, tep=0, tpb=0.5, tb=0.96, 
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
#mouse_SIR<-cbind.data.frame(time, S1, L1, I1, E1, R1, dr1, U1, Iep1, Ipb1, Ib1, df1)



#................................................
#Infected and Maintained on Rat Blood####
#................................................

times=seq(1,30,by=0.3)					#time steps to output
xstart=c(S=9,L=0,I=1,E=0,R=0,dr=0,Id=0, U=100,Iep=0,Ipb=0,Ib=0,df=0)					#beginning population size
parms=c(alpha=1,lambdaA=0.04,lambdaB=0.02,lambdaC=0.06,b=0.4,b1=2,tau=0.48, muf=0.02, mupb=0.13,mub=0.26,
        pep=0,ppb=0.10,pb=0.89,tep=0, tpb=1 , tb=0.96 ,
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

plot(time,I1,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(time,E1,col="blue",lwd=2,lty=1)
lines(time,S1,col="forestgreen",lwd=2,lty=1)
lines(time,L1,col="red",lwd=2,lty=5)
lines(time,R1,col="blue",lwd=2,lty=5)
lines(time,Id1,col="black",lwd=2,lty=1)
legend(23,10,c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1,5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp =1)
title(main="Host Dynamics with Mouse Blood Infected Fleas, 9UH1IH 100UF LD50 of 10 CFU v6 BPB only")
dev.off()
#rat plot;
#postscript("C:/Users/Cedar/Desktop/Plague/Models/Graphs/Host_ratblood_10CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)
postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Host_ratblood_10CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(time,I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(time,E,col="blue",lwd=2,lty=1)
lines(time,S,col="forestgreen",lwd=2,lty=1)
lines(time,L,col="red",lwd=2,lty=5)
lines(time,R,col="blue",lwd=2,lty=5)
lines(time,Id,col="black",lwd=2,lty=1)
legend(23,10,c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1,5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Host Dynamics with Rat Blood Infected Fleas, 9UH1IH 100UF LD50 of 10 CFU v6 BPB only")
dev.off()

#............................................................................
#FLEA Compartment

#par(mfrow=c(2,1))
#postscript("C:/Users/Cedar/Desktop/Plague/Models/Graphs/Flea_mouseblood_10CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)
postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Flea_mouseblood_10CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(time,U1,ylab="Abundance",xlab="Time",type="l",col="forestgreen",lwd=2,ylim=c(0,max(100)),lty=1)
lines(time,Iep1,col="blue",lwd=2,lty=1)
lines(time,Ipb1,col="red",lwd=2,lty=5)
lines(time,Ib1,col="red",lwd=2,lty=1)
lines(time,df1,col="black",lwd=2,lty=1)
legend(20,50,c("Uninfected", "Infectious-Early Phase", "Infectious-Partially Blocked","Infectious-Fully Blocked", "Cumulative Dead"),
       col=c("forestgreen","blue", "red","red", "black"),bty="n",lty=c(1,1, 5, 1, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Flea Dynamics with Mouse Blood Infected Fleas, 9UH1IH 100UF LD50 of 10 CFU v6 BPB only")
dev.off()
#rat plot;
#postscript("C:/Users/Cedar/Desktop/Plague/Models/Graphs/Flea_ratblood_10CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)
postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Flea_ratblood_10CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(time,U,ylab="Abundance",xlab="Time",type="l",col="forestgreen",lwd=2,ylim=c(0,max(100)),lty=1)
lines(time,Iep,col="blue",lwd=2,lty=1)
lines(time,Ipb,col="red",lwd=2,lty=5)
lines(time,Ib,col="red",lwd=2,lty=1)
lines(time,df,col="black",lwd=2,lty=1)
legend(20,50,c("Uninfected", "Infectious-Early Phase", "Infectious-Partially Blocked","Infectious-Fully Blocked", "Cumulative Dead"),
       col=c("forestgreen","blue", "red","red", "black"),bty="n",lty=c(1,1, 5, 1, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Flea Dynamics with Rat Blood Infected Fleas, 9UH1IH 100UF LD50 of 10 CFU v6 BPB only")
dev.off()



#..................................................
#Threshold of 100 CFU####
#..................................................

times=seq(1,30,by=0.3)					#time steps to output
xstart=c(S=9,L=0,I=1,E=0,R=0,dr=0,Id=0, U=100,Iep=0,Ipb=0,Ib=0,df=0)					#beginning population size
parms=c(alpha=1,lambdaA=0.04,lambdaB=0.20,lambdaC=0.07,b=0.4,b1=2,tau=0.39, muf=0.02, mupb=0.14,mub=0.20,
        pep=0,ppb=0.11,pb=0.75, tep=0.00, tpb=0.5, tb=0.88, 
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
#mouse_SIR<-cbind.data.frame(time, S1, L1, I1, E1, R1, dr1, U1, Iep1, Ipb1, Ib1, df1)

#................................................
#Infected and Maintained on Rat Blood####
#................................................


times=seq(1,30,by=0.3)					#time steps to output
xstart=c(S=9,L=0,I=1,E=0,R=0,dr=0,Id=0, U=100,Iep=0,Ipb=0,Ib=0,df=0)					#beginning population size
parms=c(alpha=1,lambdaA=0.04,lambdaB=0.02,lambdaC=0.06,b=0.4,b1=2,tau=0.48, muf=0.02, mupb=0.13,mub=0.26,
        pep=0,ppb=0.10,pb=0.89,tep=0, tpb=1 , tb=0.65 ,
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

plot(time,I1,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(time,E1,col="blue",lwd=2,lty=1)
lines(time,S1,col="forestgreen",lwd=2,lty=1)
lines(time,L1,col="red",lwd=2,lty=5)
lines(time,R1,col="blue",lwd=2,lty=5)
lines(time,Id1,col="black",lwd=2,lty=1)
legend(23,10,c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1,5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp =1)
title(main="Host Dynamics with Mouse Blood Infected Fleas, 9UH1IH 100UF LD50 of 100 CFU v6 BPB only")
dev.off()
#rat plot;
#postscript("C:/Users/Cedar/Desktop/Plague/Models/Graphs/Host_ratblood_100CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)
postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Host_ratblood_100CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(time,I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(time,E,col="blue",lwd=2,lty=1)
lines(time,S,col="forestgreen",lwd=2,lty=1)
lines(time,L,col="red",lwd=2,lty=5)
lines(time,R,col="blue",lwd=2,lty=5)
lines(time,Id,col="black",lwd=2,lty=1)
legend(23,10,c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1,5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Host Dynamics with Rat Blood Infected Fleas, 9UH1IH 100UF LD50 of 100 CFU v6 BPB only")
dev.off()

#............................................................................
#FLEA Compartment

#par(mfrow=c(2,1))
#postscript("C:/Users/Cedar/Desktop/Plague/Models/Graphs/Flea_mouseblood_100CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)
postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Flea_mouseblood_100CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(time,U1,ylab="Abundance",xlab="Time",type="l",col="forestgreen",lwd=2,ylim=c(0,max(100)),lty=1)
lines(time,Iep1,col="blue",lwd=2,lty=1)
lines(time,Ipb1,col="red",lwd=2,lty=5)
lines(time,Ib1,col="red",lwd=2,lty=1)
lines(time,df1,col="black",lwd=2,lty=1)
legend(20,50,c("Uninfected", "Infectious-Early Phase", "Infectious-Partially Blocked","Infectious-Fully Blocked", "Cumulative Dead"),
       col=c("forestgreen","blue", "red","red", "black"),bty="n",lty=c(1,1, 5, 1, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Flea Dynamics with Mouse Blood Infected Fleas, 9UH1IH 100UF LD50 of 100 CFU v6 BPB only")
dev.off()
#rat plot;
#postscript("C:/Users/Cedar/Desktop/Plague/Models/Graphs/Flea_ratblood_100CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)
postscript("/Users/jhinnebusch/R/Plague/Models/Graphs/Flea_ratblood_100CFU.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)

plot(time,U,ylab="Abundance",xlab="Time",type="l",col="forestgreen",lwd=2,ylim=c(0,max(100)),lty=1)
lines(time,Iep,col="blue",lwd=2,lty=1)
lines(time,Ipb,col="red",lwd=2,lty=5)
lines(time,Ib,col="red",lwd=2,lty=1)
lines(time,df,col="black",lwd=2,lty=1)
legend(20,50,c("Uninfected", "Infectious-Early Phase", "Infectious-Partially Blocked","Infectious-Fully Blocked", "Cumulative Dead"),
       col=c("forestgreen","blue", "red","red", "black"),bty="n",lty=c(1,1, 5, 1, 1),lwd=2,seg.len=2.0,x.intersp = 0.5, y.intersp = 1)
title(main="Flea Dynamics with Rat Blood Infected Fleas, 9UH1IH 100UF LD50 of 100 CFU v6 BPB only")
dev.off()
