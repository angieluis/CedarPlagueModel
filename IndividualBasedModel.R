###################################################################################
## Individual-based discrete stochastic model for Plague 
##  based on Cedar & Joe's model and data
###################################################################################

##params

b=0.4  # mean of 0.4 bites per day from a Poisson distribution (for non-blocked fleas)
b1=2
muf=0.02
mupb=0.13
mub=0.26
tau=0.48
lambdaA=0.04
lambdaB=0.02
lambdaC=0.06
pep=0.03
ppb=0.11
pb=0.5
tep=0.001
tpb=0.5
tb=0.8
alpha=1
mu=0.002 #rodent mortality
gamma=0.14	 #recovery rate in rodent from low-dose flea infection
epsilon=0.5  #disease induced mortality rate in rodent from high-dose flea infection
sigma=0.25  #rate to become infectious from latent class


# Total number of Rodent Hosts to start with
S=9  # susceptible rodents
L=0  # rodents latently infected (pre-infecitous)
I=1  # rodents infected (bacteremic)
E=0  # rodents exposed (non-infectious)
R=0  # rodents recovered and immune
dr=0 # rodents that died naturally
Id=0 # rodents that died of infection

H = sum(c(S,L,E,I,R)) # total hosts alive 

U=100    # uninfected fleas
Iep=0  # early phase fleas
Ipb=0 # partially blocked fleas
Ib=0  # fully blocked fleas
df=0  # dead fleas

F= sum(U,Iep,Ipb,Ib) # total fleas alive

T = 100 # total number of days to simulate

# create matrices to hold each individual's status over time
rodents <- matrix(NA,T,H)
fleas <- matrix(NA,T,F)

rodents[1,] <- sample(c(rep("S",S),rep("L",L),rep("I",I),rep("E",E),rep("R",R))) 
fleas[1,] <- sample( c(rep("U",U),rep("Iep",Iep),rep("Ipb",Ipb),rep("Ib",Ib)))

dose=matrix(0,T,H)

for(t in 2:T){
  # set up a vector to hold cumulative dose from infected fleas to link flea model to rodent model
  #dose <- rep(0, H)
###################################################################################
## Flea Model
###################################################################################

# Simulate each flea separately 
  
  for(f in 1:F){
    bites <-rpois(1,ifelse(fleas[t-1,f]=="Ib",b1,b)) # how many bites does this flea make today? Draw from a Poisson distribution. If blocked use mean of b1, otherwise b
  
    if(bites>0){
      adH <- c(which(rodents[t-1,]=="dr"),which(rodents[t-1,]=="Id")) # dead rodents so fleas wont feed on these
      # remove dead hosts from possibilty
      if(length(adH)>0){
        possible_rodents <- (1:H)[-adH]
      }else{  
        possible_rodents <-  1:H
      }
      who_bites <- sample(possible_rodents,bites,replace=TRUE) # which rodents are these bites on? Draw from those that are alive - assumes the same flea could bite the same host more than once
      status_of_bitten <- rodents[t-1,who_bites] # what class(es) are the bitten rodents in
  
      # bites are useless unless by an infected flea on uninfected host, or from an uninfected flea on Infected host
      for(i in 1:length(status_of_bitten)){ # for each of the bites
        if(status_of_bitten[i]=="S"){ # if rodent was uninfected, it could get infected but depends on cumulative dose, so keep track of that for the rodent model below
          if(fleas[t-1,f]=="Iep"){
            dose[t,who_bites[i]] <- min(dose[t,who_bites[i]] + rbinom(1,1,pep) * tep ,1)#if it transmits, it adds dose of tep
          }
          if(fleas[t-1,f]=="Ipb"){
            dose[t,who_bites[i]] <- min(dose[t,who_bites[i]] + rbinom(1,1,ppb) * tpb ,1)#if it transmits, it adds dose of tpb
          }
          if(fleas[t-1,f]=="Ib"){
            dose[t,who_bites[i]] <- min(dose[t,who_bites[i]] + rbinom(1,1,pb) * tb ,1) #if it transmits, it adds dose of tb
          }
        }
        if(status_of_bitten[i]=="I"){ # if rodent was infected, it could transmit to uninfected flea
          if(fleas[t-1,f]=="U"){
            # if already infected this time step from another bite, stay infected, otherwise see if this bite leads to infection
            # already infected from another bite? if so stay infected, if not can become infected
            if(!is.na(fleas[t,f])){
              if(fleas[t,f]=="Iep"){
                fleas[t,f] <- "Iep"
              } 
            }
            if(is.na(fleas[t,f])){
              fleas[t,f] <- ifelse(rbinom(1,1,alpha), "Iep", "U")
            }
          }
        }
        if(status_of_bitten[i]!="I"){ # if didn't bite an infected, then stay where you are
          fleas[t,f] <- fleas[t-1,f] 
        }
         # if not uninfected, stay where you are, for now (just modeling becoming infected here, other transitions can happen below)
        if(fleas[t-1,f]!="U"){ 
          fleas[t,f] <- fleas[t-1,f] 
        }
        
      } #i, end of each bites loop

    } # end of if bites>0

    if(bites==0){
      # remain in the same class unless die or move in next lines
      fleas[t,f] <- fleas[t-1,f]
    }
    
    ## need to also add stay in same class if bite but not on U or I.... 
    
    # simulate deaths and other movements 
    # assume these deaths and movements happen after any bites
  
    # simulates deaths
    prob.die <- ifelse(fleas[t-1,f]=="Ipb",mupb,ifelse(fleas[t-1,f]=="Ib",mub,muf))
    fleas[t,f] <- ifelse(rbinom(1,1,prob=prob.die), "df", fleas[t,f])
  
    # if alive, then allow the other few movements
    if(fleas[t,f] == "Iep"){ # could move back to U (lambdaB) or become partially blocked (lambdaA), or leave by lambda C (consider dead here, i don't like this- need to ask Joe)
      fleas[t,f] <-  sample(c("U","Ipb","dr","Iep"),1,prob=c(lambdaB,lambdaA,lambdaC,1-(lambdaA+lambdaB+lambdaC))) # definitions say rate for lambda so need to see if ok to use probs 
    }
    if(fleas[t,f] == "Ipb"){ # could move into fully blocked
      fleas[t,f] <- ifelse(rbinom(1,1,prob=tau),"Ib","Ipb") # right now using probability but may want to look at how long have been here in this class instead
    }
  
  } # f,  end of flea loop



###################################################################################
## Rodent Model
###################################################################################

  # ok to use rates (gamma, epsilon, sigma) as probabilities? how variable are these?
  
#### Now Simulate each rodent separately using info from flea bites above
  for (h in 1:H){
    if(rodents[t-1,h]=="dr"){ # if dead, stay dead
      rodents[t,h] <- "dr"
    }
    if(rodents[t-1,h]=="Id"){ # if dead, stay dead
      rodents[t,h] <- "Id"
    }
    if(rodents[t-1,h]=="S"){ # if was susceptible, check that it survived, then check for cumulative infectious dose from bites
      die <- rbinom(1,1,mu)
      if(die==1){
        rodents[t,h] <- "dr"
      }else{
        if(dose[t,h]==0){
          rodents[t,h] <- "S"
        }
        if(dose[t,h]>0){
          rodents[t,h] <- ifelse(rbinom(1,1,prob=dose[t,h]),"L", "E")
        }
      }
    } # end S

    if(rodents[t-1,h]=="L"){ # if was latent, check that it survived, then could move to I class
      die <- rbinom(1,1,mu)
      if(die==1){
        rodents[t,h] <- "dr"
      }else{
        rodents[t,h] <- ifelse(rbinom(1,1,prob=sigma), "I", "L")
      }
    } # end L
    
    if(rodents[t-1,h]=="I"){ # if was infectious, check if it survived, 
      rodents[t,h] <- ifelse(rbinom(1,1,prob=epsilon), "Id", "I")
    } # end I
    
    if(rodents[t-1,h]=="E"){ # if was exposed, check that it survived, then could move to R class
      die <- rbinom(1,1,mu)
      if(die==1){
        rodents[t,h] <- "dr"
      }else{
        rodents[t,h] <- ifelse(rbinom(1,1,prob=gamma), "R", "E")
      }
    } # end E

    if(rodents[t-1,h]=="R"){ # if was recovered, check if it survived 
      rodents[t,h] <- ifelse(rbinom(1,1,prob=mu), "dr", "R")
    } # end R
  
  }# end rodent model  
    
} # t 



rodent.time.series <- data.frame(
  S = apply(rodents,1,function(x){length(which(x=="S"))}),
  L = apply(rodents,1,function(x){length(which(x=="L"))}),
  I = apply(rodents,1,function(x){length(which(x=="I"))}),
  E = apply(rodents,1,function(x){length(which(x=="E"))}),
  R = apply(rodents,1,function(x){length(which(x=="R"))}),
  dr = apply(rodents,1,function(x){length(which(x=="dr"))}),
  Id = apply(rodents,1,function(x){length(which(x=="Id"))})
)  

flea.time.series <- data.frame(
  U = apply(fleas,1,function(x){length(which(x=="U"))}),
  Iep = apply(fleas,1,function(x){length(which(x=="Iep"))}),
  Ipb = apply(fleas,1,function(x){length(which(x=="Ipb"))}),
  Ib = apply(fleas,1,function(x){length(which(x=="Ib"))}),
  df = apply(fleas,1,function(x){length(which(x=="df"))})
)  

  
plot.ts(rodent.time.series$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,max(10)),lty=1)
lines(rodent.time.series$E,col="blue",lwd=2,lty=1)
lines(rodent.time.series$S,col="forestgreen",lwd=2,lty=1)
lines(rodent.time.series$L,col="red",lwd=2,lty=5)
lines(rodent.time.series$R,col="blue",lwd=2,lty=5)
lines(rodent.time.series$Id,col="black",lwd=2,lty=1)
legend("topright",c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp =0.5, y.intersp =1)
title(main="Mouse Blood, 1 CFU")
