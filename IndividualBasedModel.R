###################################################################################
## Individual-based discrete stochastic model for Plague 
##  based on Cedar & Joe's model and data
###################################################################################

### add calculation of R0. How many fleas (& who) does the first I transmit to. Then track
### those fleas and see how many hosts they infect (that make it to I).

# If who_bites==1 & status_of_bitten == "I" & fleas[t,f,s] <- "Iep", then keep track
# of that flea [f] & who it infects unless it dies or moves back to "U".
# keep vector of fleas.for.R0 to watch, but need to remove those that go back to U from Iep
# keep vector of L.for.R0 to track who those infect and if they survive to I


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

  rodents[1,,] <- c(rep("I",I),rep("S",S),rep("L",L),rep("E",E),rep("R",R)) # infected first
  fleas[1,,] <- sample( c(rep("U",U),rep("Iep",Iep),rep("Ipb",Ipb),rep("Ib",Ib)))

  R0 <- rep(0,n.sim)

  for(s in 1:n.sim){  
    
    fleas.for.R0 <-numeric() # keep track of the fleas that the first I infects
    S.to.watch <- numeric() # keep track of S's bitten by fleas to watch
    L.for.R0 <- numeric() # keep track of hosts that become L from those fleas
    
    for(t in 2:T){
      prob.notI <- rep(1, H) # the probability it doesn't move into L then I (will be a fxn of all bites)

      ###################################################################################
      ## Flea Model
      ###################################################################################

      # Simulate each flea separately 
  
      for(f in 1:F){
        bites <-rpois(1,ifelse(fleas[t-1,f,s]=="Ib",b1,b)) # how many bites does this flea make today? Draw from a Poisson distribution. If blocked use mean of b1, otherwise b
        adH <- c(which(rodents[t-1,,s]=="dr"),which(rodents[t-1,,s]=="Id")) # dead rodents so fleas wont feed on these
        # only want to evaluate bites if not all the rodents are dead 
        # need expression to determine if there are any dead rodents, and if so are they all of them
        eval.bites <- TRUE
        if(length(adH)>0){
          eval.bites <- ifelse(length(adH)<H,TRUE,FALSE)
        }
        if(bites>0 & eval.bites){ # if there are bites and rodents left to bite
          # remove dead hosts from possibilty
          if(length(adH)>0){
            possible_rodents <- (1:H)[-adH]
          }else{  
            possible_rodents <-  1:H
          }
          who_bites <- sample(possible_rodents,bites,replace=TRUE) # which rodents are these bites on? Draw from those that are alive - assumes the same flea could bite the same host more than once
          status_of_bitten <- rodents[t-1,who_bites,s] # what class(es) are the bitten rodents in
  
          # if who_bites is #1 and status_of_bitten is I, then need to track it for R0.
          
          # bites are useless unless by an infected flea on uninfected host, or from an uninfected flea on Infected host
          for(i in 1:length(status_of_bitten)){ # for each of the bites
            if(status_of_bitten[i]=="S"){ # if rodent was uninfected, it could get infected but depends on cumulative number of bites from different categories, so keep track of that for the rodent model below
              prob.t <- 0
              if(fleas[t-1,f,s]=="Iep"){
                prob.t <- (rbinom(1,1,pep) * tep)
                prob.notI[who_bites[i]] <- prob.notI[who_bites[i]] * (1-prob.t) # prob it doesn't move to infectious trajectory is product of result of previous bites and this bite (with prob of moving to I is tep, so prob of not moving there to 1-tep)
              }
              if(fleas[t-1,f,s]=="Ipb"){
                prob.t <- (rbinom(1,1,ppb) * tpb)
                prob.notI[who_bites[i]] <- prob.notI[who_bites[i]] * (1-prob.t)
              }
              if(fleas[t-1,f,s]=="Ib"){
                prob.t <- (rbinom(1,1,pb) * tb)
                prob.notI[who_bites[i]] <- prob.notI[who_bites[i]] * (1-prob.t)
              }
              # if some prob of transmission & if this flea is in the vector to watch, keep track for R0
              if(prob.t>0 & length(which(fleas.for.R0==f))>0){
                S.to.watch <- c(S.to.watch, who_bites[i])
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
                  # if the flea just moved to Iep phase and was bitten by first infected rodent, start tracking it for R0 (if not already)
                  if(fleas[t,f,s]=="Iep" & who_bites[i]==1){
                    fleas.for.R0 <- unique(c(fleas.for.R0,f))
                  }
                  
                }
              } #if "U"
            } # if bit "I"
            if(status_of_bitten[i]!="I"){ # if didn't bite an infected, then stay where you are
              fleas[t,f,s] <- fleas[t-1,f,s] 
            }
            # if not uninfected, stay where you are, for now (just modeling becoming infected here, other transitions can happen below)
            if(fleas[t-1,f,s]!="U"){ 
              fleas[t,f,s] <- fleas[t-1,f,s] 
            }
        
          } #i, end of each bites loop

        } else { # if not evaluating bites
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
        if(fleas[t,f,s] == "Iep"){ # could move back to U (lambdaB+lambdaC) or become partially blocked (lambdaA)
          fleas[t,f,s] <-  sample(c("U","Ipb","U","Iep"),1,prob=c(lambdaB,lambdaA,lambdaC,1-(lambdaA+lambdaB+lambdaC))) # definitions say rate for lambda so need to see if ok to use probs 
        
          if(fleas[t,f,s]=="U"){ # if the flea moved back to U
            if(length(which(fleas.for.R0==f))>0){ # & if the flea was already in the fleas Ro vector, remove it
              fleas.for.R0 <- fleas.for.R0[which(fleas.for.R0==f)]
            }
          }
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
          rodents[t,h,s] <- "dr"
        }
        if(rodents[t-1,h,s]=="Id"){ # if dead, stay dead
          rodents[t,h,s] <- "Id"
        }
        if(rodents[t-1,h,s]=="S"){ # if was susceptible, check that it survived, then check for cumulative infectious prob from bites
          die <- rbinom(1,1,mu)
          if(die==1){
            rodents[t,h,s] <- "dr"
          }else{
            if(prob.notI[h]==1){
              rodents[t,h,s] <- "S"
            }
            if(prob.notI[h]<1){
              rodents[t,h,s] <- ifelse(rbinom(1,1,prob=(1-prob.notI[h])),"L", "E")
              if(length(which(S.to.watch==h))>0 & rodents[t,h,s]=="L"){ # keep track for R0 if in watch vector
                L.for.R0 <- unique(c(L.for.R0, h))
              }
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
          if(rodents[t,h,s]=="I" & length(which(L.for.R0==h))>0){
            R0[s] <- R0[s]+1 #add to R0
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

  ####### Summaries to output
  
  rodent.ts.list <- list()
  flea.ts.list <- list()
  for(n in 1:n.sim){
    rodent.time.series <- data.frame(
      S = apply(rodents[,,n],1,function(x){length(which(x=="S"))}),
      L = apply(rodents[,,n],1,function(x){length(which(x=="L"))}),
      I = apply(rodents[,,n],1,function(x){length(which(x=="I"))}),
      E = apply(rodents[,,n],1,function(x){length(which(x=="E"))}),
      R = apply(rodents[,,n],1,function(x){length(which(x=="R"))}),
      dr = apply(rodents[,,n],1,function(x){length(which(x=="dr"))}),
      Id = apply(rodents[,,n],1,function(x){length(which(x=="Id"))})
    )  
    
    flea.time.series <- data.frame(
      U = apply(fleas[,,n],1,function(x){length(which(x=="U"))}),
      Iep = apply(fleas[,,n],1,function(x){length(which(x=="Iep"))}),
      Ipb = apply(fleas[,,n],1,function(x){length(which(x=="Ipb"))}),
      Ib = apply(fleas[,,n],1,function(x){length(which(x=="Ib"))}),
      df = apply(fleas[,,n],1,function(x){length(which(x=="df"))})
    )  
    rodent.ts.list[[n]] <- rodent.time.series
    flea.ts.list[[n]] <- flea.time.series
    
  }
  
  rodent.ts.array<-array(as.numeric(unlist(rodent.ts.list)),dim=c(100,7,100),dimnames=list(paste("time",1:100,sep=""),c("S","L","I","E","R","dr","Id"),paste("sim",1:100,sep="")))
  
  flea.ts.array<-array(as.numeric(unlist(flea.ts.list)),dim=c(100,5,100),dimnames=list(paste("time",1:100,sep=""),c("U","Iep","Ipb","Ib","df"),paste("sim",1:100,sep="")))
  
  
  mean.rodent.ts <- apply(rodent.ts.array,c(1,2),mean)
  mean.flea.ts <- apply(flea.ts.array,c(1,2),mean)
  
  
  return(list(rodent.status.array=rodents, flea.status.array=fleas, rodent.ts.array=rodent.ts.array,flea.ts.array=flea.ts.array,mean.rodent.ts=mean.rodent.ts,mean.flea.ts=mean.flea.ts,R0=R0))

} #end of function




######################### function to run and compare several scenarios as once


print.plague.IBM <- function(params, # list of vectors of params with labels for each model
                             xstart, 
                             T, 
                             n.sim, 
                             plot=TRUE, # do you want comparison plots?, if TRUE, spit out pdf of plots in dir
                             plot.name=NULL){ #optional name for plot pdf file
  out <- list()
  for(i in 1:length(params)){
    out[[i]] <- plague.IBM(params[[i]], xstart, T, n.sim)
  }
  names(out) <- names(params)
  
  if(plot==TRUE){
    if(length(plot.name)==0){
      plot.name="plague.IBM.plot"
    }
    postscript(plot.name, horizontal = FALSE, onefile = FALSE, paper = "special", height = 7, width = 11)
    par(mfrow=c(length(params)/3,3))
    
    for(i in 1:length(params)){
      dat <- as.data.frame(out[[i]]$mean.rodent.ts)
      plot.ts(dat$I,ylab="Abundance",xlab="Time",type="l",col="red",lwd=2,ylim=c(0,sum(dat[1,])),lty=1)
      lines(dat$E,col="blue",lwd=2,lty=1)
      lines(dat$S,col="forestgreen",lwd=2,lty=1)
      lines(dat$L,col="red",lwd=2,lty=5)
      lines(dat$R,col="blue",lwd=2,lty=5)
      lines(dat$Id,col="black",lwd=2,lty=1)
      legend("topright",c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp =0.5, y.intersp =1)
      title(main=names(params)[i])
    }
    
    dev.off()
  }
  
  out[[length(params)+1]] <- unlist(lapply(out,function(y){as.data.frame(y$mean.rodent.ts)$Id[dim(y$mean.rodent.ts)[1]]}))
  names(out)[length(params)+1] <- "Infected.dead.rodents"
    
  return(out)
}


