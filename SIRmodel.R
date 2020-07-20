##################################################################
## Model Formulation
##################################################################

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
