###################################################################################
## Individual-based discrete stochastic model for Plague 
##  based on Cedar & Joe's model and data
###################################################################################


######## To work on :
#######     # when all rodents are dead, the fleas have no one to bite on
            # error at sample(possible_rodents, bites, replace = TRUE)
            # return other summaries - time series


##params
params=c(alpha=1,lambdaA=0.035,lambdaB=0.20,lambdaC=0.07,b=0.4,b1=2,tau=0.39, muf=0.02, mupb=0.14,mub=0.20,
        pep=0.03,ppb=0.11,pb=0.5, tep=1, tpb=1, tb=1, 
        mu=0.002,sigma=0.25,gamma=0.14,epsilon=0.5)	#set parameter values
xstart=c(S=9,L=0,I=1,E=0,R=0,dr=0,Id=0, U=100,Iep=0,Ipb=0,Ib=0,df=0)					#beginning population size
T = 100 # total number of days to simulate
n.sim = 100


plague.IBM <- function(params, # vector of params with names, see below
                       xstart, # vector of starting values with names
                       T,      # number of time steps to simulate
                       n.sim  # number of simulations
                       ){

  b=params["b"]          #biting rate of fleas
  b1=params["b1"]        #biting rate of blocked fleas
  tep=params["tep"]      #proportion of INFECTIOUS bites from early-phase fleas
  tpb=params["tpb"]      #proportion of INFECTIOUS bites from partially blocked fleas
  tb=params["tb"]        #proportion of INFECTIOUS bites from fully blocked fleas
  mu=params["mu"]		     #natural mortality rate of rodent (uninfected)
  pep=params["pep"]      #Proportion of ep fleas that transmit
  ppb=params["ppb"]      #Proportion of partially blocked fleas that transmit
  pb=params["pb"]        #Proportion of fully blocked fleas that transmit
  gamma=params["gamma"]	 #recovery rate in rodent from low-dose flea infection
  epsilon=params["epsilon"]  #disease induced mortality rate in rodent from high-dose flea infection
  sigma=params["sigma"]  #rate to become infectious from latent class
  alpha=params["alpha"]  #proportion of fleas infected from host
  muf=params["muf"]      #natural mortality of flea
  mupb=params["mupb"]    #mortality of partially blocked flea
  mub=params["mub"]      #mortality of blocked flea
  lambdaA=params["lambdaA"] #rate of developing partial blockage from EP
  lambdaB=params["lambdaB"] #rate of clearing infection in EP (back to uninfected)
  lambdaC=params["lambdaC"] #rate of leaving EP (still infected but not enough to block)
  tau=params["tau"]         #rate of developing full blockage
  
  S = xstart["S"]  # susceptible rodents
  L = xstart["L"]  # rodents latently infected (pre-infecitous)
  I = xstart["I"]  # rodents infected (bacteremic)
  E = xstart["E"]  # rodents exposed (non-infectious)
  R = xstart["R"]  # rodents recovered and immune
  dr= xstart["dr"] # all rodents that died
  Id= xstart["Id"] # rodents that died of infection
  
  U = xstart["U"]    # uninfected fleas
  Iep = xstart["Iep"]  # early phase fleas
  Ipb = xstart["Ipb"] # partially blocked fleas
  Ib = xstart["Ib"]  # fully blocked fleas
  df = xstart["df"]  # dead fleas
  
  H = sum(c(S,L,E,I,R)) # total hosts alive 

  F= sum(U,Iep,Ipb,Ib) # total fleas alive


  # create matrices to hold each individual's status over time
  rodents <- array(data=NA, dim=c(T,H,n.sim))
  fleas   <- array(data=NA, dim=c(T,F,n.sim))

  rodents[1,,] <- sample(c(rep("S",S),rep("L",L),rep("I",I),rep("E",E),rep("R",R))) 
  fleas[1,,] <- sample( c(rep("U",U),rep("Iep",Iep),rep("Ipb",Ipb),rep("Ib",Ib)))


  for(s in 1:n.sim){  
    for(t in 2:T){
      dose <- rep(0, H)

      ###################################################################################
      ## Flea Model
      ###################################################################################

      # Simulate each flea separately 
  
      for(f in 1:F){
        bites <-rpois(1,ifelse(fleas[t-1,f,s]=="Ib",b1,b)) # how many bites does this flea make today? Draw from a Poisson distribution. If blocked use mean of b1, otherwise b
  
        if(bites>0){
          adH <- c(which(rodents[t-1,,s]=="dr"),which(rodents[t-1,,s]=="Id")) # dead rodents so fleas wont feed on these
          # remove dead hosts from possibilty
          if(length(adH)>0){
            possible_rodents <- (1:H)[-adH]
          }else{  
            possible_rodents <-  1:H
          }
          who_bites <- sample(possible_rodents,bites,replace=TRUE) # which rodents are these bites on? Draw from those that are alive - assumes the same flea could bite the same host more than once
          status_of_bitten <- rodents[t-1,who_bites,s] # what class(es) are the bitten rodents in
  
          # bites are useless unless by an infected flea on uninfected host, or from an uninfected flea on Infected host
          for(i in 1:length(status_of_bitten)){ # for each of the bites
            if(status_of_bitten[i]=="S"){ # if rodent was uninfected, it could get infected but depends on cumulative dose, so keep track of that for the rodent model below
              if(fleas[t-1,f,s]=="Iep"){
                dose[who_bites[i]] <- min(dose[who_bites[i]] + rbinom(1,1,pep) * tep ,1)#if it transmits, it adds dose of tep
              }
              if(fleas[t-1,f,s]=="Ipb"){
                dose[who_bites[i]] <- min(dose[who_bites[i]] + rbinom(1,1,ppb) * tpb ,1)#if it transmits, it adds dose of tpb
              }
              if(fleas[t-1,f,s]=="Ib"){
                dose[who_bites[i]] <- min(dose[who_bites[i]] + rbinom(1,1,pb) * tb ,1) #if it transmits, it adds dose of tb
              }
            }
            if(status_of_bitten[i]=="I"){ # if rodent was infected, it could transmit to uninfected flea
              if(fleas[t-1,f,s]=="U"){
                # if already infected this time step from another bite, stay infected, otherwise see if this bite leads to infection
                # already infected from another bite? if so stay infected, if not can become infected
                if(!is.na(fleas[t,f,s])){
                  if(fleas[t,f,s]=="Iep"){
                    fleas[t,f,s] <- "Iep"
                  } 
                }
                if(is.na(fleas[t,f,s])){
                  fleas[t,f,s] <- ifelse(rbinom(1,1,alpha), "Iep", "U")
                }
              }
            }
            if(status_of_bitten[i]!="I"){ # if didn't bite an infected, then stay where you are
              fleas[t,f,s] <- fleas[t-1,f,s] 
            }
            # if not uninfected, stay where you are, for now (just modeling becoming infected here, other transitions can happen below)
            if(fleas[t-1,f,s]!="U"){ 
              fleas[t,f,s] <- fleas[t-1,f,s] 
            }
        
          } #i, end of each bites loop

        } # end of if bites>0

        if(bites==0){
          # remain in the same class unless die or move in next lines
          fleas[t,f,s] <- fleas[t-1,f,s]
        }
    
        ## need to also add stay in same class if bite but not on U or I.... 
    
        # simulate deaths and other movements 
        # assume these deaths and movements happen after any bites
  
        # simulates deaths
        prob.die <- ifelse(fleas[t-1,f,s]=="Ipb",mupb,ifelse(fleas[t-1,f,s]=="Ib",mub,muf))
        fleas[t,f,s] <- ifelse(rbinom(1,1,prob=prob.die), "df", fleas[t,f,s])
  
        # if alive, then allow the other few movements
        if(fleas[t,f,s] == "Iep"){ # could move back to U (lambdaB) or become partially blocked (lambdaA), or leave by lambda C (consider dead here, i don't like this- need to ask Joe)
          fleas[t,f,s] <-  sample(c("U","Ipb","dr","Iep"),1,prob=c(lambdaB,lambdaA,lambdaC,1-(lambdaA+lambdaB+lambdaC))) # definitions say rate for lambda so need to see if ok to use probs 
        }
        if(fleas[t,f,s] == "Ipb"){ # could move into fully blocked
          fleas[t,f,s] <- ifelse(rbinom(1,1,prob=tau),"Ib","Ipb") # right now using probability but may want to look at how long have been here in this class instead
        }
  
     } # f,  end of flea loop



    ###################################################################################
    ## Rodent Model
    ###################################################################################

      # ok to use rates (gamma, epsilon, sigma) as probabilities? how variable are these?
  
    #### Now Simulate each rodent separately using info from flea bites above
      for (h in 1:H){
        if(rodents[t-1,h,s]=="dr"){ # if dead, stay dead
          rodents[t,,s] <- "dr"
        }
        if(rodents[t-1,h,s]=="Id"){ # if dead, stay dead
          rodents[t,h,s] <- "Id"
        }
        if(rodents[t-1,h,s]=="S"){ # if was susceptible, check that it survived, then check for cumulative infectious dose from bites
          die <- rbinom(1,1,mu)
          if(die==1){
            rodents[t,h,s] <- "dr"
          }else{
            if(dose[h]==0){
              rodents[t,h,s] <- "S"
            }
            if(dose[h]>0){
              rodents[t,h,s] <- ifelse(rbinom(1,1,prob=dose[h]),"L", "E")
            }
          }
        } # end S

        if(rodents[t-1,h,s]=="L"){ # if was latent, check that it survived, then could move to I class
          die <- rbinom(1,1,mu)
          if(die==1){
            rodents[t,h,s] <- "dr"
          }else{
            rodents[t,h,s] <- ifelse(rbinom(1,1,prob=sigma), "I", "L")
          }
        } # end L
    
        if(rodents[t-1,h,s]=="I"){ # if was infectious, check if it survived, 
          rodents[t,h,s] <- ifelse(rbinom(1,1,prob=epsilon), "Id", "I")
        } # end I
    
        if(rodents[t-1,h,s]=="E"){ # if was exposed, check that it survived, then could move to R class
          die <- rbinom(1,1,mu)
          if(die==1){
            rodents[t,h,s] <- "dr"
          }else{
            rodents[t,h,s] <- ifelse(rbinom(1,1,prob=gamma), "R", "E")
          }
        } # end E

        if(rodents[t-1,h,s]=="R"){ # if was recovered, check if it survived 
          rodents[t,h,s] <- ifelse(rbinom(1,1,prob=mu), "dr", "R")
        } # end R
  
      }# end rodent model  
    
    } # t 
  } # s


  return(list(rodents, fleas))
} #end of function



out <- plague.IBM(params, # vector of params with names, see below
                              xstart, # vector of starting values with names
                              T,      # number of time steps to simulate
                              n.sim  # number of simulations
)




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
