##################################################################
## Model Formulation
##################################################################

### divide transmission terms by number of vertebrate hosts to match Ross-MacDonald model

SIR.model.fleasperhost=function(t,x,params){
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
  
  dS = -S*(b*Iep*pep + b*Ipb*ppb + b1*Ib*pb)/H - mu*S   #Susceptible Rodent Pop
  dL = S*(b*Iep*pep*tep + b*Ipb*ppb*tpb + b1*Ib*pb*tb)/H - sigma*L - mu*L  #infection of host, not yet infectious
  dI = sigma*L - I*(mu+epsilon)  #infectious rodent hosts
  dE = S*(b*Iep*pep*(1-tep) + b*Ipb*ppb*(1-tpb) + b1*Ib*pb*(1-tb))/H - E*(mu+gamma) #exposed rodents (non-infectious)
  dR = gamma*E - R*mu  #recovery of rodent
  dr = mu*H + epsilon*I  #dead rodents
  Id = epsilon*I
  
  dU=-I*(U*b*alpha)/H + Iep*lambdaB + Iep*lambdaC - U*muf    #Uninfected Flea pop
  dIep=I*U*b*alpha/H - (Iep*lambdaB + Iep*lambdaC + Iep*lambdaA + Iep*muf) #Infection of early phase fleas 
  dIpb=Iep*lambdaA - Ipb*(tau+mupb) #Infection of partially blocked fleas
  dIb=tau*Ipb - mub*Ib  #Infection of fully blocked fleas
  df=muf*(U+Iep) + mupb*Ipb + mub*Ib #dead fleas
  cumb= tau*Ipb #dead blocked fleas
  
  
  list(c(dS,dL,dI,dE,dR,dr,Id, dU,dIep,dIpb,dIb,df,cumb))
  
}



#######################################################################
## Function to calculate R0 based on this model
#######################################################################

R0.function <- function(params, m ){
  R0 <- with(as.list(params),
             (-alpha*m*pep*sigma*tep*b^2/((-epsilon - mu)*(-mu - sigma)*(-lambdaA - lambdaB - lambdaC - muf)) + alpha*lambdaA*m*ppb*sigma*tpb*b^2/((-epsilon - mu)*(-mu - sigma)*(-mupb - tau)*(-lambdaA - lambdaB - lambdaC - muf)) + alpha*b*b1*lambdaA*m*pb*sigma*tau*tb/(mub*(-epsilon - mu)*(-mu - sigma)*(-mupb - tau)*(-lambdaA - lambdaB - lambdaC - muf)))^(1/2)
             
  )
  return(R0)
}




#######################################################################
## Function to run and compare several scenarios at once
#######################################################################

print.plague.SIRmodel <- function(model = SIR.model.fleasperhost,
                                  params, # list of vectors of params with labels for each model
                                  xstart, 
                                  T, # number of days to simulate
                                  plot=TRUE, # do you want comparison plots?, if TRUE, spit out pdf of plots in dir
                                  plot.name=NULL){ #optional name for plot pdf file
  times=seq(0,T,by=0.3)					#time steps to output
  
  out <- list()
  for(i in 1:length(params)){
    out [[i]] <- as.data.frame(lsoda(xstart, times, model, params[[i]], hmax = 0.001))  #run the model
    
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
      lines(dat$time,dat$cumb,col="orange",lwd=2,lty=5)
      legend("topright",c("U", "Iep", "Ipb","Ib", "dead","cumBlocked"),col=c("blue","forestgreen","magenta", "red", "black","orange"),bty="n",lty=c(1, 1, 1, 1, 5, 5),lwd=2,seg.len=2.0,x.intersp =0.5, y.intersp =1)
      title(main=names(params)[i])
    }
    
    dev.off()
  }
  
  #out[[length(params)+1]] <- unlist(lapply(out,function(y){as.data.frame(y$mean.rodent.ts)$Id[dim(y$mean.rodent.ts)[1]]}))
  #names(out)[length(params)+1] <- "Infected.dead.rodents"
  
  ## table of summaries 
  
  summary.table <- data.frame(scenario=names(params),R0=rep(NA,length(params)),Idead=rep(NA,length(params)),Recovered=rep(NA,length(params)),cumblocked=rep(NA,length(params)))
  
  for(i in 1:length(out)){
    dat <- as.data.frame(out[[i]])
    summary.table$R0[i] = R0.function(params[[i]], m=xstart[8]/sum(xstart[1:7]))
    summary.table$Idead[i] = dat$Id[dim(dat)[1]]
    summary.table$Recovered[i] = dat$R[dim(dat)[1]]
    summary.table$cumblocked[i] = dat$cumb[dim(dat)[1]]
    
  }
  
  print(summary.table)
  return(out)
}




print.plague.one.SIRmodel <- function(model = SIR.model.fleasperhost,
                                      params, # one vector of params
                                      xstart, 
                                      T, # number of days to simulate
                                      plot=TRUE, # do you want host and vector plots? (not pdfs)
                                      title=NULL){ 
  times=seq(0,T,by=0.3)					#time steps to output
  
  out <- as.data.frame(lsoda(xstart, times, model, params, hmax = 0.001))  #run the model
  
  par(mfrow=c(2,1))
  
  dat <- as.data.frame(out)
  plot(dat$time,dat$I,ylab="Abundance",xlab="Time (days)",type="l",col="red",lwd=2,ylim=c(0,sum(dat[1,2:8])),lty=1)
  lines(dat$time,dat$E,col="blue",lwd=2,lty=1)
  lines(dat$time,dat$S,col="forestgreen",lwd=2,lty=1)
  lines(dat$time,dat$L,col="red",lwd=2,lty=5)
  lines(dat$time,dat$R,col="blue",lwd=2,lty=5)
  lines(dat$time,dat$Id,col="black",lwd=2,lty=1)
  legend("topright",c("S", "L", "I","E", "R", "I-dead"),col=c("forestgreen","red", "red","blue", "blue", "black"),bty="n",lty=c(1, 5, 1, 1, 5, 1))
  title(main=paste(title,"host"))
  
  
  plot(dat$time,dat$U,ylab="Abundance",xlab="Time (days)",type="l",col="blue",lwd=2,ylim=c(0,sum(dat[1,9:13])),lty=1)
  lines(dat$time,dat$Iep,col="forestgreen",lwd=2,lty=1)
  lines(dat$time,dat$Ipb,col="magenta",lwd=2,lty=1)
  lines(dat$time,dat$Ib,col="red",lwd=2,lty=1)
  lines(dat$time,dat$df,col="black",lwd=2,lty=5)
  lines(dat$time,dat$cumb,col="orange",lwd=2,lty=5)
  legend("topright",c("U", "Iep", "Ipb","Ib", "dead","cumBlocked"),col=c("blue","forestgreen","magenta", "red", "black","orange"),bty="n",lty=c(1, 1, 1, 1, 5, 5))
  title(main=paste(title,"fleas"))
  
  
  
  
  
  return(out)
}


#######################################################################
## Run the function to compare several scenarios at once
#######################################################################


# rat blood at vector to host ratio of 5
#SIR.rat.comparisons <- print.plague.SIRmodel(
#    model = SIR.model.fleasperhost,
#    params = params.rat, 
#    xstart=c(S=9,L=0,I=1,E=0,R=0,dr=0,Id=0, 
#             U=50,Iep=0,Ipb=0,Ib=0,df=0,cumb=0),
#    T=100,
#    plot.name="Plague SIR m=5, rat comparison.pdf")
  