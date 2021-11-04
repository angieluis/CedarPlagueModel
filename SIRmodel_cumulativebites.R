##################################################################
## Model Formulation
##################################################################

### divide transmission terms by number of vertebrate hosts to match Ross-MacDonald model
### calculate probability move from S to L,I vs S to E,R based on cumulative probabilty across all bites

SIR.model.fleasperhost.cumbites=function(t,x,params){
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
  
  S=x[1]  # susceptible rodents
  L=x[2]  # rodents latently infected (pre-infecitous)
  I=x[3]  # rodents infected (bacteremic)
  E=x[4]  # rodents exposed (non-infectious)
  R=x[5]  # rodents recovered and immune
  dr=x[6] # all rodents that died
  Id=x[7] # rodents that died of infection
  
  U=x[8]    # uninfected fleas
  Iep=x[9]  # early phase fleas
  Ipb=x[10] # partially blocked fleas
  Ib=x[11]  # fully blocked fleas
  df=x[12]  # dead fleas

  H  = S+R+E+L+I # total number of rodent hosts
  
  # probability move to L and I when S is bitten : cumulative probabilities across all bites
  prob.infectious <- 1 - ((1-tep)^(Iep*b*S/H))*((1-tpb)^(Ipb*b*S/H))*((1-tb)^(Ib*b1*S/H))
  
  
  dS = -S*(b*Iep*pep + b*Ipb*ppb + b1*Ib*pb)/H - mu*S   #Susceptible Rodent Pop
  dL = prob.infectious*S*(b*Iep*pep + b*Ipb*ppb + b1*Ib*pb)/H - sigma*L - mu*L  #infection of host, not yet infectious
  dI = sigma*L - I*(mu+epsilon)  #infectious rodent hosts
  dE = (1-prob.infectious)*S*(b*Iep*pep + b*Ipb*ppb + b1*Ib*pb)/H - E*(mu+gamma) #exposed rodents (non-infectious)
  dR = gamma*E - R*mu  #recovery of rodent
  dr = mu*H + epsilon*I  #dead rodents
  Id = epsilon*I
  
  dU=-I*(U*b*alpha)/H + Iep*(lambdaB + lambdaC) - U*muf    #Uninfected Flea pop
  dIep=I*U*b*alpha/H - (Iep*lambdaB + Iep*lambdaC + Iep*lambdaA + Iep*muf) #Infection of early phase fleas 
  dIpb=Iep*lambdaA - Ipb*(tau+mupb) #Infection of partially blocked fleas
  dIb=tau*Ipb - mub*Ib  #Infection of fully blocked fleas
  df=muf*(U+Iep) + mupb*Ipb + mub*Ib  #dead fleas
  
  
  
  list(c(dS,dL,dI,dE,dR,dr,Id, dU,dIep,dIpb,dIb,df))
  
}





######################### function to run and compare several scenarios at once


print.plague.SIRmodel <- function(model = SIR.model.fleasperhost.cumbites,
                             params, # list of vectors of params with labels for each model
                             xstart, 
                             T, # number of days to simulate
                             plot=TRUE, # do you want comparison plots?, if TRUE, spit out pdf of plots in dir
                             plot.name=NULL){ #optional name for plot pdf file
  times=seq(0,T,by=0.3)					#time steps to output
  
  out <- list()
  for(i in 1:length(params)){
    out [[i]] <- as.data.frame(lsoda(xstart, times, model, params[[i]]))  #run the model
    
  }
  names(out) <- names(params)
  
  if(plot==TRUE){
    if(length(plot.name)==0){
      plot.name="plague.IBM.plot"
    }
    pdf(plot.name, onefile = TRUE, paper = "special", height = 7, width = 11)
    par(mfrow=c(length(params)/3,3))
    
    for(i in 1:length(params)){
      dat <- as.data.frame(out[[i]])
      plot(dat$time,dat$I,ylab="Abundance",xlab="Time (days)",type="l",col="red",lwd=2,ylim=c(0,sum(dat[1,2:8])),lty=1)
      lines(dat$time,dat$E,col="blue",lwd=2,lty=1)
      lines(dat$time,dat$S,col="forestgreen",lwd=2,lty=1)
      lines(dat$time,dat$L,col="red",lwd=2,lty=5)
      lines(dat$time,dat$R,col="blue",lwd=2,lty=5)
      lines(dat$time,dat$Id,col="black",lwd=2,lty=1)
      legend("topright",c("Susceptible", "Latent", "Infectious","Exposed", "Recovered", "Infected-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1),lwd=2,seg.len=2.0,x.intersp =0.5, y.intersp =1)
      title(main=names(params)[i])
    }
    
    par(mfrow=c(length(params)/3,3))
    
    for(i in 1:length(params)){
      dat <- as.data.frame(out[[i]])
      plot(dat$time,dat$U,ylab="Abundance",xlab="Time (days)",type="l",col="blue",lwd=2,ylim=c(0,sum(dat[1,9:13])),lty=1)
      lines(dat$time,dat$Iep,col="forestgreen",lwd=2,lty=1)
      lines(dat$time,dat$Ipb,col="magenta",lwd=2,lty=1)
      lines(dat$time,dat$Ib,col="red",lwd=2,lty=1)
      lines(dat$time,dat$df,col="black",lwd=2,lty=5)
      legend("topright",c("U", "Iep", "Ipb","Ib", "dead"),col=c("blue","forestgreen","magenta", "red", "black"),bty="n",lty=c(1, 1, 1, 1, 5),lwd=2,seg.len=2.0,x.intersp =0.5, y.intersp =1)
      title(main=names(params)[i])
    }
    
    dev.off()
  }
  
  #out[[length(params)+1]] <- unlist(lapply(out,function(y){as.data.frame(y$mean.rodent.ts)$Id[dim(y$mean.rodent.ts)[1]]}))
  #names(out)[length(params)+1] <- "Infected.dead.rodents"
  
  return(out)
}


