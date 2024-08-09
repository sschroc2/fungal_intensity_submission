### Integral Projection Model & code for Figures S5, S6, S7, S8 (IPM model outputs)
### Step 1: Load Sierra dataset
### Step 2: Set IPM functions
### Step 3: Set parameters
### Step 4: Set functions for aggregation metrics
### Step 5: Set functions to simulate data & extract values
### Step 6: S7 plot
### Step 7: S5 plot
### Step 8: S6 plot
### Step 9: S8 plot


#Clear the environment to be safe
rm(list=ls())

#Load necessary packages
library(data.table)
library(ggplot2)
library(tictoc)
library(pracma)
library(gridExtra)
library(cowplot)
library(grid)

###########################################
#### Step 1: Load data & set time step ####
###########################################

#Read in just Sierra data
sierra_dat = fread("data/formatted/analysis_aggregation_dataset.csv")[dataset == "sierra_nevada"]

#Time step is one day everywhere
deltat=1


#########################################
#### Step 2: Establish IPM functions ####
#########################################

#Create bins of load
set_discretized_values = function(min_size, max_size, bins){

  # Calculates the necessary parameters to use the midpoint rule to evaluate
  # the IPM model

  # Parameters
  # ----------
  # min_size : The lower bound of the integral
  # max_size : The upper bound of the integral
  # bins : The number of bins in the discretized matrix

  # Returns
  # -------
  # list
  # min_size, max_size, bins, bnd (edges of discretized kernel), y (midpoints),
  # h (width of cells)


  # Set the edges of the discretized kernel
  bnd = min_size+c(0:bins)*(max_size-min_size) / bins

  # Set the midpoints of the discretizing kernel. Using midpoint rule for evaluation
  y = 0.5 * (bnd[1:bins] + bnd[2:(bins + 1)])

  # Width of cells
  h = y[2] - y[1]

  return(list(min_size=min_size, max_size=max_size, bins=bins, bnd=bnd, y=y,
                h=h))

}

#Write a function for growth from all bins to all bins
growth_function = function(x_lower, x_upper, x_now, params){
	# Bd growth function on a host
	#
	# Parameters
	# ----------
	# x_upper: Upper bound on Bd load bin at time t + 1
	# x_lower: Lower bound on Bd load bin at time t + 1
	# x_now: Bd load at time t
	# params : list
	#
	# Returns
	# -------
	# : Probability of increasing from x at t to between x_upper and x_lower at time t + 1

	deltat = params$deltat
	a = params$growth_int
	b = params$growth_size
	mu = (a / (1 - b))*(1 - b^deltat) + b^deltat * x_now
	sigma = sqrt((params$growth_sd^2 / (1 - b^2)) * (1 - b^(2*deltat)))
	prob = pnorm(x_upper, mean=mu, sd=sigma) - pnorm(x_lower, mean=mu, sd=sigma)
	return(prob)

}

#Calculate survival probability as a function of load
survival_function = function(x_now, params){
	# Probability of surviving to the next time step with Bd load of x_now at time t
	#
	# Parameters
	# ----------
	# x_now : Bd load at time t
	#
	# Returns
	# -------
	# : Survival probability over a time step

	deltat = params$deltat
	prob=1-pnorm((x_now-params$surv_LD50)/params$surv_sigma2,mean=0,sd=1)
	return(prob^deltat)
}

#Fixed probability of infection loss
loss_function = function(x_now, params){
	# Probability of losing infection given infection of x_now by the next time step
	#
	# Parameters
	# ----------
	# x_now : Bd load at time t
	#
	# Returns
	# -------
	# : Survival probability over a time step

	deltat = params$deltat
	#Parametrized for loss probability in a 3 day window
	lossprob3day=exp(params$loss_LD50+params$loss_sigma2*x_now)/(1+exp(params$loss_LD50+params$loss_sigma2*x_now))
	#Can do the conversion
	#(1-lossprob3day)=(1-lossprob1day)^3
	#(1-lossprob3day)^(1/3)=1-lossprob1day
	#lossprob1day=1-(1-lossprob3day)^(1/3)
	prob=1-(1-lossprob3day)^(1/3)
	return(prob^deltat)
}

initial_inf_function = function(x_lower, x_upper, params){
	# Initial infection function gives the probability of a new infection starting
	# beetween x_upper and x_lower
	#
	# Parameters
	# -----------
	# x_upper: Upper bound on Bd load bin at time t + 1
	# x_lower: Lower bound on Bd load bin at time t + 1
	#
	# Returns
	# -------
	# : Probability of new Bd load being between x_upper and x_lower

	mu = params$init_inf_int
	sigma = params$init_inf_sd
	prob = pnorm(x_upper, mean=mu, sd=sigma) - pnorm(x_lower, mean=mu, sd=sigma)
	return(prob)
}

transmission_function = function(Z, params){
	# Transmission function dependent on the zoospore pool
	#
	# Parameters
	# ----------
	# Z: Zoospore density at time t
	# params : list of parameters, should include trans_beta which is beta'*delta_t (delta t is baked in!)
	#
	# Returns
	# -------
	# transmission probability over a time step

	deltat = params$deltat
	prob = 1 - exp(-Z*params$trans_beta*deltat)
	return(prob)
}


build_P_matrix = function(lower, upper, bins, params){
	# Transition matrix from I to I
	#
	# Parameters
	# ----------
	# lower : lower bound on kernel
	# upper : upper bound on kernel
	# bins : number of bins in kernel
	# params : list of IPM parameters
	#
	# Return
	# ------
	# : I to I transition matrix

	ipm_bnds = set_discretized_values(lower, upper, bins)
	x_lower = ipm_bnds$bnd[1:(bins)]
	x_upper = ipm_bnds$bnd[2:(bins + 1)]
	x = ipm_bnds$y # midpoints

	#For each value of x_now based on midpoints, 
  #this calculates a vector of probabilities corresponding
	#to the pairs of x_lower and x_upper that define each bin. Column i corresponds
  #to a value of x_now and row j corresponds to growth into bin j.
	G = sapply(x, function(i) growth_function(x_lower, x_upper, i, params))

	S = survival_function(x, params)
	L = loss_function(x, params)
  
	#%*% is matrix multiplication
	#diag turns a vector into a diagonal matrix, in this case
	
	#P is matrix that returns the joint probability of survival based on x_now,
	#loss based on x_now, and growth from x_now into a bin with column i corresponding
	#to a value of x_now and row j corresponding to the bin in the next time step
	P = G %*% diag(S * (1 - L))
	return(P)
}

build_full_P_matrix = function(P, Z, lower, upper, bins, params){
	# Build the full transition matrix assuming S - I - Z classes
	# Shedding into a zoospore pool and transmission from the Z pool

	ipm_bnds = set_discretized_values(lower, upper, bins)
	x_lower = ipm_bnds$bnd[1:(bins)]
	x_upper = ipm_bnds$bnd[2:(bins + 1)]
	x = ipm_bnds$y # midpoints

	fullP = array(0, dim=c(bins + 2, bins + 2))
	deltat = params$deltat

	inf_prob = transmission_function(Z, params)
	StoS = params$s0^deltat*(1 - inf_prob)

	#A vector of probabilities of initial infection going into each bin
	G0_init = initial_inf_function(x_lower, x_upper, params)
	
	#A vector of probabilities of survival and infection going into each bin
	StoI = params$s0^deltat*inf_prob*G0_init

	S = survival_function(x, params)
	L = loss_function(x, params)
	
	#A vector of probabilites of survival and then infection loss from each midpoint
	ItoS = S*L

	#A vector of the amount of Z shed at each midpoint per time unit (day in my case)
	ItoZ = params$lambda*deltat * exp(x)
	
	#Probability of survival of spores for a time unit (day in my case)
	ZtoZ = 0 # Make this zero to make it easier to count zoospores outside of the
	#genotype specific function

	# Fill in the matrix
	fullP[2:(bins + 1), 2:(bins + 1)] = P
	fullP[1, 1] = StoS
	fullP[2:(bins + 1), 1] = StoI
	fullP[1, 2:(bins + 1)] = ItoS
	fullP[bins + 2, 2:(bins + 1)] = ItoZ
	fullP[bins + 2, bins + 2] = ZtoZ
  
	#Quantity of flow from column class/bin to row class/bin. For host to host, 
	#this is the proportion that will transition from column class/bin to row 
	#class/bin. For last row, it is the amount of Z produced by different columns
	#per day.
	return(fullP)
}


build_full_P_matrix_singlegen = function(P, Z, lower, upper, bins, params){
  # Build the full transition matrix assuming S - I - Z classes
  # Shedding into a zoospore pool and transmission from the Z pool
  
  ipm_bnds = set_discretized_values(lower, upper, bins)
  x_lower = ipm_bnds$bnd[1:(bins)]
  x_upper = ipm_bnds$bnd[2:(bins + 1)]
  x = ipm_bnds$y # midpoints
  
  fullP = array(0, dim=c(bins + 2, bins + 2))
  deltat = params$deltat
  
  inf_prob = transmission_function(Z, params)
  StoS = params$s0^deltat*(1 - inf_prob)
  
  #A vector of probabilities of initial infection going into each bin
  G0_init = initial_inf_function(x_lower, x_upper, params)
  
  #A vector of probabilities of survival and infection going into each bin
  StoI = params$s0^deltat*inf_prob*G0_init
  
  S = survival_function(x, params)
  L = loss_function(x, params)
  
  #A vector of probabilites of survival and then infection loss from each midpoint
  ItoS = S*L
  
  #A vector of the amount of Z shed at each midpoint per time unit (day in my case)
  ItoZ = params$lambda*deltat * exp(x)
  
  #Probability of survival of spores for a time unit (day in my case)
  ZtoZ = params$s_z^deltat
  
  # Fill in the matrix
  # Technically, I only need to do this once and then update the first rows
  fullP[2:(bins + 1), 2:(bins + 1)] = P
  fullP[1, 1] = StoS
  fullP[2:(bins + 1), 1] = StoI
  fullP[1, 2:(bins + 1)] = ItoS
  fullP[bins + 2, 2:(bins + 1)] = ItoZ
  fullP[bins + 2, bins + 2] = ZtoZ
  
  #Quantity of flow from column class/bin to row class/bin. For host to host, 
  #this is the proportion that will transition from column class/bin to row 
  #class/bin. For last row, it is the amount of Z produced by different columns
  #per day.
  return(fullP)
}


####################################
#### Step 3: Set IPM Parameters ####
####################################

# Set IPM approximation parameters
lower = -10
upper = 30
bins = 150

# Set IPM parameters
params = list()

params$growth_int = 0.6640405
params$growth_size = 0.9279308
params$growth_sd = 1.124049
params$surv_LD50 = 15.94864
params$surv_sigma2 = 2.952637
params$loss_LD50=1.213-0.151*20
params$loss_sigma2=-0.472
params$lambda = 0.06
params$trans_beta = 4.035e-5
params$s_z = 0.009481621
params$s0 = 0.9997718
params$init_inf_int = 3.382
params$init_inf_sd = 2.707928
params$deltat = 1
params$r=0.1711157
params$q=0.9984626
params$mu=5e-6
params_focal=params

#Then overwriting some
params$init_inf_int=3.382
params$init_inf_sd=1.5
params$growth_sd=0.5
params$growth_size=0.96

params$growth_int/(1-params$growth_size)
params$surv_LD50=25
params$surv_sigma2=2.952637
params$lambda=0.015

params_temp=params

num_gens=10

rs=rep(0.1711157,num_gens)*(1-0.999*linspace(0,1,num_gens)^10)
muss=linspace(params$surv_LD50,params$surv_LD50+10,num_gens)

#######################################################
#### Step 4: Set functions for aggregation metrics ####
#######################################################

#Write a function for getting Poulin's D from density of infected hosts in each bin
PD=function(N,Is,bins,x,min_thresh){
  #Vector of number of infected hosts in each class
  Inf_vec=round(N*Is/(sum(Is)))
  
  #Fix rounding issues which should be inconsequential but can cause errors
  if((N-sum(Inf_vec))!=0){
    Inf_vec[which.max(Inf_vec)]=Inf_vec[which.max(Inf_vec)]+(N-sum(Inf_vec))
  }
  
  num_pars=rep(NA,N)
  incrementing_index=1
  for(bindex in 1:bins){
    if(incrementing_index<N){
      num_pars[incrementing_index:(incrementing_index+max(c(Inf_vec[bindex]-1),0))]=exp(x[bindex]) 
    }
    #print(incrementing_index:(incrementing_index+max(c(Inf_vec[bindex]-1),0)))
    incrementing_index=incrementing_index+Inf_vec[bindex]
    #Sys.sleep(0.5)
  }
  
  num_pars_final=num_pars[which(num_pars>min_thresh)]
  
  Poulin_sum=0
  for(i in 1:length(num_pars_final)){
    for(j in 1:i){
      Poulin_sum=Poulin_sum+num_pars_final[j]
    }
  }
  
  return(1-2*Poulin_sum/(mean(num_pars_final)*length(num_pars_final)*(length(num_pars_final)+1)))
}

#McVinish formulation of Poulin's D
PD_mcvinish=function(N,Is,bins,x,min_thresh){
  #Vector of number of infected hosts in each class
  Inf_vec=round(N*Is/(sum(Is)))
  
  #Fix rounding issues which should be inconsequential but can cause errors
  if((N-sum(Inf_vec))!=0){
    Inf_vec[which.max(Inf_vec)]=Inf_vec[which.max(Inf_vec)]+(N-sum(Inf_vec))
  }
  
  num_pars=rep(NA,N)
  incrementing_index=1
  for(bindex in 1:bins){
    if(incrementing_index<N){
      num_pars[incrementing_index:(incrementing_index+max(c(Inf_vec[bindex]-1),0))]=exp(x[bindex]) 
    }
    #print(incrementing_index:(incrementing_index+max(c(Inf_vec[bindex]-1),0)))
    incrementing_index=incrementing_index+Inf_vec[bindex]
    #Sys.sleep(0.5)
  }
  
  num_pars_final=num_pars[which(num_pars>exp(min_thresh))]
  
  Poulin_sum=0
  for(i in 1:length(num_pars_final)){
    for(j in 1:length(num_pars_final)){
      Poulin_sum=Poulin_sum+abs(num_pars_final[j]-num_pars_final[i])
    }
  }
  
  return(Poulin_sum/(2*mean(num_pars_final)*length(num_pars_final)^2))
}


##############################################################
#### Step 5: Set functions simulate data & extract values ####
##############################################################

#Choose a number of time steps
steps = 5*365

#Run a multiple genotype simulation
multi_genotype_runner=function(params_10,init,steps){
  
  #Column i lists density of S, I in each bin, then density of Z
  #at time i
results = array(NA, dim=c(bins + 2, steps + 1,num_gens))
#Midpoints for estimations
x=set_discretized_values(lower,upper,bins)$y
x_lower=x-(upper-lower)/(2*bins)
x_upper=x+(upper-lower)/(2*bins)
nonwt_factor=1e-6
init_conds = array(0,dim=c(bins+2,num_gens))

#This works as long as starting Z is zero. Otherwise, there can be a discrepancy
#between WT and others
init_conds[,1]=init
init_conds[,2:num_gens]=init_conds[,1]*nonwt_factor

#Put initial conditions in the first column, all layers
results[, 1,] = init_conds
  
# Just build this once
P = array(NA,c(bins,bins,num_gens))
for(g in 1:num_gens){
  P[,,g]=build_P_matrix(lower, upper, bins, params_10[[g]])
}

fullP=array(NA,dim=c(bins+2,bins+2,num_gens))

n_next=array(NA,dim=c(bins+2,num_gens))

#Iterate time steps starting from 1 all the way to the number of steps
for(i in 1:steps){
  for(g in 1:num_gens){  
  n_now = results[,i,g]
    
  Z_now = n_now[bins+2]
  fullP[,,g] = build_full_P_matrix(P[,,g], Z_now, lower, upper, bins, params_10[[g]])
    
  n_next[,g] = fullP[,,g] %*% n_now
  }
  
  #Add in fecundity as well
  if((1-params_10[[1]]$q*sum(results[1:(bins+1),i,]))>0){
    crowd_mod=(1-params_10[[1]]$q*sum(results[1:(bins+1),i,]))
  }else{
    crowd_mod=0
  }
  
  for(g in 1:num_gens){
    n_next[1,g]=n_next[1,g]+(params_10[[g]]$r*(sum(results[1:(bins+1),i,g]))*(1-params_10[[g]]$mu-params_10[[g]]$mu/(num_gens-1))+params_10[[g]]$mu/(num_gens-1)*sum(rs*colSums(results[1:(bins+1),i,])))*crowd_mod    
  }
  
  
  #This sums shed and persisting zoospores
  Z_next=sum(n_next[bins+2,])+params_10[[1]]$s_z*results[bins+2,i,1]
  n_next[bins+2,]=Z_next
  results[, i + 1,] = n_next
}
  return(results)
}

summary_finder=function(results,thresh,PDsize){
  steps=dim(results)[2]-1
  # Calculate means and variances
  ipm_bnds = set_discretized_values(lower, upper, bins)
  x = ipm_bnds$y # midpoints
  
  thresh_index=min(which(x>thresh))
  
  x_threshed=x[thresh_index:bins]
  
  bins_threshed=length(x_threshed)
  
  trun = 1
  
  #And calculate overall mean and variation
  Ns = colSums(results[(thresh_index+1):(bins+1),trun:(steps+1)])
  norm_dist = t(t(results[(thresh_index+1):(bins+1),trun:(steps+1)])/Ns)
  
  norm_means_everyone=rep(NA,steps+1)
  norm_vars_everyone=rep(NA,steps+1)
  nat_means_everyone=rep(NA,steps+1)
  nat_vars_everyone=rep(NA,steps+1)
  log_cv_everyone=rep(NA,steps+1)
  PD_everyone=rep(NA,steps+1)
  PD_everyone_Mc=rep(NA,steps+1)
  
  for(i in 1:(steps+1)){
    #What is the mean on the log scale from sum of freq of bins times midpoint of bins?
    norm_means_everyone[i] = sum(norm_dist[, i] * x_threshed)
    
    #What is the variance on the log scale using freq of bins and midpoints of bins?
    norm_vars_everyone[i] = sum(norm_dist[, i] * x_threshed^2) - norm_means_everyone[i]^2
    
    nat_means_everyone[i] = sum(norm_dist[,i]*exp(x_threshed))
    
    nat_vars_everyone[i] = sum(norm_dist[,i]*exp(x_threshed)^2)-nat_means_everyone[i]^2
    
    if(Ns[i]>0){
      #PD_everyone[i]=PD(PDsize,PDsize*norm_dist[,i],bins_threshed,x_threshed,-Inf)      
      PD_everyone_Mc[i]=PD_mcvinish(1e3,1e3*norm_dist[,i],bins_threshed,x_threshed,-Inf)      
    }
  }
  
  log_cv_everyone=log10(sqrt(nat_vars_everyone)/nat_means_everyone)
  
  return(cbind(norm_means_everyone,sqrt(norm_vars_everyone),nat_means_everyone,sqrt(nat_vars_everyone),log_cv_everyone,PD_everyone,PD_everyone_Mc))
}

summary_finder_multi_sparse=function(results,thresh,reduce_factor){
  time_indices=(1:floor((dim(results)[2]-1)/reduce_factor))*reduce_factor
  # Calculate means and variances
  ipm_bnds = set_discretized_values(lower, upper, bins)
  x = ipm_bnds$y # midpoints
  
  thresh_index=min(which(x>thresh))
  
  x_threshed=x[thresh_index:bins]
  
  bins_threshed=length(x_threshed)
  
  trun = 1
  
  #And calculate overall mean and variation
  Ns = apply(results[(thresh_index+1):(bins+1),time_indices,],MARGIN=2,FUN=sum)
  norm_dist = t(t(apply(results[(thresh_index+1):(bins+1),time_indices,],MARGIN=c(1,2),FUN=sum))/Ns)
  
  norm_means_everyone=rep(NA,length(time_indices))
  norm_vars_everyone=rep(NA,length(time_indices))
  nat_means_everyone=rep(NA,length(time_indices))
  nat_vars_everyone=rep(NA,length(time_indices))
  log_cv_everyone=rep(NA,length(time_indices))
  PD_everyone=rep(NA,length(time_indices))
  PD_everyone_Mc=rep(NA,length(time_indices))
  
  for(i in 1:length(time_indices)){
    #What is the mean on the log scale from sum of freq of bins times midpoint of bins?
    norm_means_everyone[i] = sum(norm_dist[, i] * x_threshed)
    
    #What is the variance on the log scale using freq of bins and midpoints of bins?
    norm_vars_everyone[i] = sum(norm_dist[, i] * x_threshed^2) - norm_means_everyone[i]^2
    
    nat_means_everyone[i] = sum(norm_dist[,i]*exp(x_threshed))
    
    nat_vars_everyone[i] = sum(norm_dist[,i]*exp(x_threshed)^2)-nat_means_everyone[i]^2
    
    if(Ns[i]>0){
      #PD_everyone[i]=PD(1e3,1e3*norm_dist[,i],bins_threshed,x_threshed,-Inf)      
      PD_everyone_Mc[i]=PD_mcvinish(1e3,1e3*norm_dist[,i],bins_threshed,x_threshed,-Inf)      
    }
  }
  
  log_cv_everyone=log10(sqrt(nat_vars_everyone)/nat_means_everyone)
  
  return(cbind(norm_means_everyone,sqrt(norm_vars_everyone),nat_means_everyone,sqrt(nat_vars_everyone),log_cv_everyone,PD_everyone,PD_everyone_Mc))
}

multi_genotype_plotter=function(output,reduce_factor){
  timevec=(0:(dim(output)[2]-1))/365
  Hs=apply(output[1:(bins+1),,],MARGIN = 2, FUN = sum)
  prev=apply(output[2:(bins+1),,],MARGIN = 2, FUN = sum)/Hs
  time_indices=(1:floor((dim(output)[2]-1)/reduce_factor))*reduce_factor
  summary=summary_finder_multi_sparse(output,-Inf,reduce_factor)
  par(mfrow=c(2,2))
  plot(timevec,Hs,type="l",xlab="Year",ylab="Host density")
  plot(timevec,prev,type="l",xlab="Year",ylab="Prevalence")
  plot(timevec[time_indices],summary[,1],type="l",xlab="Year",ylab="Mean log load")
 plot(timevec[time_indices],summary[,2],type="l",xlab="Year",ylab="SD log load")
}


params_temp=params
growth_ints=linspace(params_temp$growth_int,0.1*params_temp$growth_int,num_gens)

params_tol=list()
params_a=list()

for(g in 1:num_gens){
  params_tol[[g]]=params_temp
  params_tol[[g]]$r=rs[g]
  params_tol[[g]]$surv_LD50=muss[g]
  
  params_a[[g]]=params_temp
  params_a[[g]]$r=rs[g]
  params_a[[g]]$growth_int=growth_ints[g]
}

runner1=function(params,init,steps){
  
  #Column i lists density of S, I in each bin, then density of Z
  #at time i
  results = array(NA, dim=c(bins + 2, steps + 1))
  #Midpoints for estimations
  x=set_discretized_values(lower,upper,bins)$y
  x_lower=x-(upper-lower)/(2*bins)
  x_upper=x+(upper-lower)/(2*bins)
  
  
  #Put initial conditions in the first column, all layers
  results[, 1] = init
  
  # Just build this once
  P=build_P_matrix(lower, upper, bins, params)
  
  n_next=array(NA,dim=c(bins+2))
  
  #Iterate time steps starting from 1 all the way to the number of steps
  for(i in 1:steps){
    n_now = results[,i]
    
    Z_now = n_now[bins+2]
    fullP = build_full_P_matrix_singlegen(P[,], Z_now, lower, upper, bins, params)
    
    n_next = fullP %*% n_now
    
    #Add in fecundity as well
    if((1-params$q*sum(results[1:(bins+1),i]))>0){
      crowd_mod=(1-params$q*sum(results[1:(bins+1),i]))
    }else{
      crowd_mod=0
    }
    
    n_next[1]=n_next[1]+(params$r*(sum(results[1:(bins+1),i]))*crowd_mod)    
    
    results[, i + 1] = n_next
  }
  return(results)
}


###########################
#### Step 6: Figure S7 ####
###########################

temp=runner1(params_temp,c(1,rep(0,bins),1),365)
summary_temp=summary_finder(temp,-Inf,1e4) 

tol_results=multi_genotype_runner(params_tol,temp[,366],30*365)
a_results=multi_genotype_runner(params_a,temp[,366],60*365)

phase_planes_tol=summary_finder_multi_sparse(tol_results,-Inf,60)
phase_planes_a=summary_finder_multi_sparse(a_results,-Inf,60)

eco_df=data.frame(summary_temp) 
fwrite(eco_df, "data/model/eco_df.csv") #Used for Fig. 6 in main manuscript
evo_df=data.frame(phase_planes_a) 
fwrite(evo_df, "data/model/evo_df.csv") #Used for Fig. 6 in main manuscript
evo_tol_df=data.frame(phase_planes_tol)

pd_model=ggplot(data=eco_df)+geom_path(data=eco_df, aes(x=norm_means_everyone*log10(5)/log(5), y=PD_everyone_Mc),color="black")  + ylab("Poulin's D") + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) +
  theme_classic() + theme_bw() + theme(legend.title = element_text(size=14),legend.text = element_text(size=12),axis.title = element_text(size=14),legend.position = "none")+
  geom_path(data=evo_df, aes(x=norm_means_everyone*log10(5)/log(5), y=PD_everyone_Mc),color="blue")+
  geom_path(data=evo_tol_df, aes(x=norm_means_everyone*log10(5)/log(5), y=PD_everyone_Mc),color="red")+
  labs(tag="A")

lcv_model=ggplot(data=eco_df)+geom_path(data=eco_df, aes(x=norm_means_everyone*log10(5)/log(5), y=log_cv_everyone),color="black")  + ylab(expression(paste(Log[10],"(CV [natural scale])",sep=""))) + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) +
  theme_classic() + theme_bw() + theme(legend.title = element_text(size=14),legend.text = element_text(size=12),axis.title = element_text(size=14),legend.position = "none")+
  geom_path(data=evo_df, aes(x=norm_means_everyone*log10(5)/log(5), y=log_cv_everyone),color="blue")+
  geom_path(data=evo_tol_df, aes(x=norm_means_everyone*log10(5)/log(5), y=log_cv_everyone),color="red")+
  labs(tag="B")

cv_model=ggplot(data=eco_df)+geom_path(data=eco_df, aes(x=norm_means_everyone*log10(5)/log(5), y=V4/nat_means_everyone),color="black")  + ylab("CV [natural scale]") + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) +
  theme_classic() + theme_bw() + theme(legend.title = element_text(size=14),legend.text = element_text(size=12),axis.title = element_text(size=14),legend.position = "none")+
  geom_path(data=evo_df, aes(x=norm_means_everyone*log10(5)/log(5), y=V4/nat_means_everyone),color="blue")+
  geom_path(data=evo_tol_df, aes(x=norm_means_everyone*log10(5)/log(5), y=V4/nat_means_everyone),color="red")+
  labs(tag="C")

cvl_model=ggplot(data=eco_df)+geom_path(data=eco_df, aes(x=norm_means_everyone*log10(5)/log(5), y=V2/norm_means_everyone),color="black")  + ylab(expression(paste("CV [",log[10]," scale]",sep=""))) + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) +
  theme_classic() + theme_bw() + theme(legend.title = element_text(size=14),legend.text = element_text(size=12),axis.title = element_text(size=14),legend.position = "none")+
  geom_path(data=evo_df, aes(x=norm_means_everyone*log10(5)/log(5), y=V2/norm_means_everyone),color="blue")+
  geom_path(data=evo_tol_df, aes(x=norm_means_everyone*log10(5)/log(5), y=V2/norm_means_everyone),color="red")+
  labs(tag="D")

png("results/plots/Tolerance_supplment.png",width=6,height=6,units="in",res=300)
agg_plots_all = grid.arrange(
  grobs = list(pd_model, lcv_model,cv_model, cvl_model),
  layout_matrix = rbind(c(1, 2),
                        c(3, 4))
)
dev.off()

###########################
#### Step 7: Figure S5 ####
###########################

temp_focal=runner1(params_focal,c(1,rep(0,bins),1),365)
summary_temp_focal=summary_finder(temp_focal,-Inf,1e4)

eco_focal_df=data.frame(summary_temp_focal)

pd_focal=ggplot(data=eco_focal_df)+geom_path(data=eco_focal_df, aes(x=norm_means_everyone*log10(5)/log(5), y=PD_everyone_Mc),color="black")  + ylab("Poulin's D") + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) +
  theme_classic() + theme_bw() + theme(legend.title = element_text(size=14),legend.text = element_text(size=12),axis.title = element_text(size=14),legend.position = "none")+
  scale_y_continuous(limit=range(summary_temp[,7],na.rm=T))+
  labs(tag="A")

lcv_focal=ggplot(data=eco_focal_df)+geom_path(data=eco_focal_df, aes(x=norm_means_everyone*log10(5)/log(5), y=log_cv_everyone),color="black")  + ylab(expression(paste(Log[10],"(CV [natural scale])",sep=""))) + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) +
  theme_classic() + theme_bw() + theme(legend.title = element_text(size=14),legend.text = element_text(size=12),axis.title = element_text(size=14),legend.position = "none")+
  labs(tag="B")

cv_focal=ggplot(data=eco_focal_df)+geom_path(data=eco_focal_df, aes(x=norm_means_everyone*log10(5)/log(5), y=V4/nat_means_everyone),color="black")  + ylab("CV [natural scale]") + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) +
  theme_classic() + theme_bw() + theme(legend.title = element_text(size=14),legend.text = element_text(size=12),axis.title = element_text(size=14),legend.position = "none")+
  labs(tag="C")

cvl_focal=ggplot(data=eco_focal_df)+geom_path(data=eco_focal_df, aes(x=norm_means_everyone*log10(5)/log(5), y=V2/norm_means_everyone),color="black")  + ylab(expression(paste("CV [",log[10]," scale]",sep=""))) + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) +
  theme_classic() + theme_bw() + theme(legend.title = element_text(size=14),legend.text = element_text(size=12),axis.title = element_text(size=14),legend.position = "none")+
  labs(tag="D")

png("results/plots/Focal_supplement.png",width=6,height=6,units="in",res=300)
focal_plots_all = grid.arrange(
  grobs = list(pd_focal, lcv_focal,cv_focal, cvl_focal),
  layout_matrix = rbind(c(1, 2),
                        c(3, 4))
)
dev.off()

###########################
#### Step 8: Figure S6 ####
###########################

params_nearfocal=params_focal
params_nearfocal$init_inf_sd=params$init_inf_sd
params_nearfocal$growth_size=params$growth_size
temp_nearfocal=runner1(params_nearfocal,c(1,rep(0,bins),1),3*365)
summary_temp_nearfocal=summary_finder(temp_nearfocal,-Inf,1e3)

eco_nearfocal_df=data.frame(summary_temp_nearfocal)

pd_nearfocal=ggplot(data=eco_nearfocal_df)+geom_path(data=eco_nearfocal_df, aes(x=norm_means_everyone*log10(5)/log(5), y=PD_everyone_Mc),color="black")  + ylab("Poulin's D") + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) +
  theme_classic() + theme_bw() + theme(legend.title = element_text(size=14),legend.text = element_text(size=12),axis.title = element_text(size=14),legend.position = "none")+
  scale_y_continuous(limit=range(summary_temp[,7],na.rm=T))+
  labs(tag="A")

lcv_nearfocal=ggplot(data=eco_nearfocal_df)+geom_path(data=eco_nearfocal_df, aes(x=norm_means_everyone*log10(5)/log(5), y=log_cv_everyone),color="black")  + ylab(expression(paste(Log[10],"(CV [natural scale])",sep=""))) + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) +
  theme_classic() + theme_bw() + theme(legend.title = element_text(size=14),legend.text = element_text(size=12),axis.title = element_text(size=14),legend.position = "none")+
  labs(tag="B")

cv_nearfocal=ggplot(data=eco_nearfocal_df)+geom_path(data=eco_nearfocal_df, aes(x=norm_means_everyone*log10(5)/log(5), y=V4/nat_means_everyone),color="black")  + ylab("CV [natural scale]") + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) +
  theme_classic() + theme_bw() + theme(legend.title = element_text(size=14),legend.text = element_text(size=12),axis.title = element_text(size=14),legend.position = "none")+
  labs(tag="C")

cvl_nearfocal=ggplot(data=eco_nearfocal_df)+geom_path(data=eco_nearfocal_df, aes(x=norm_means_everyone*log10(5)/log(5), y=V2/norm_means_everyone),color="black")  + ylab(expression(paste("CV [",log[10]," scale]",sep=""))) + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) +
  theme_classic() + theme_bw() + theme(legend.title = element_text(size=14),legend.text = element_text(size=12),axis.title = element_text(size=14),legend.position = "none")+
  labs(tag="D")

png("results/plots/nearfocal_supplement.png",width=6,height=6,units="in",res=300)
nearfocal_plots_all = grid.arrange(
  grobs = list(pd_nearfocal, lcv_nearfocal,cv_nearfocal, cvl_nearfocal),
  layout_matrix = rbind(c(1, 2),
                        c(3, 4))
)
dev.off()


###########################
#### Step 9: Figure S8 ####
###########################

params_hiLD50=params
params_hiLD50$surv_LD50=1e4
temp_hiLD50=runner1(params_hiLD50,c(1,rep(0,bins),1),3*365)
summary_temp_hiLD50=summary_finder(temp_hiLD50,-Inf,1e3)

eco_hiLD50_df=data.frame(summary_temp_hiLD50)

pd_hiLD50=ggplot(data=eco_hiLD50_df)+geom_path(data=eco_hiLD50_df, aes(x=norm_means_everyone*log10(5)/log(5), y=PD_everyone_Mc),color="black")  + ylab("Poulin's D") + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) +
  theme_classic() + theme_bw() + theme(legend.title = element_text(size=14),legend.text = element_text(size=12),axis.title = element_text(size=14),legend.position = "none")+
  scale_y_continuous(limit=range(summary_temp[,7],na.rm=T))+
  labs(tag="A")

lcv_hiLD50=ggplot(data=eco_hiLD50_df)+geom_path(data=eco_hiLD50_df, aes(x=norm_means_everyone*log10(5)/log(5), y=log_cv_everyone),color="black")  + ylab(expression(paste(Log[10],"(CV [natural scale])",sep=""))) + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) +
  theme_classic() + theme_bw() + theme(legend.title = element_text(size=14),legend.text = element_text(size=12),axis.title = element_text(size=14),legend.position = "none")+
  labs(tag="B")

cv_hiLD50=ggplot(data=eco_hiLD50_df)+geom_path(data=eco_hiLD50_df, aes(x=norm_means_everyone*log10(5)/log(5), y=V4/nat_means_everyone),color="black")  + ylab("CV [natural scale]") + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) +
  theme_classic() + theme_bw() + theme(legend.title = element_text(size=14),legend.text = element_text(size=12),axis.title = element_text(size=14),legend.position = "none")+
  labs(tag="C")

cvl_hiLD50=ggplot(data=eco_hiLD50_df)+geom_path(data=eco_hiLD50_df, aes(x=norm_means_everyone*log10(5)/log(5), y=V2/norm_means_everyone),color="black")  + ylab(expression(paste("CV [",log[10]," scale]",sep=""))) + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) +
  theme_classic() + theme_bw() + theme(legend.title = element_text(size=14),legend.text = element_text(size=12),axis.title = element_text(size=14),legend.position = "none")+
  labs(tag="D")

png("results/plots/hiLD50_supplement.png",width=6,height=6,units="in",res=300)
hiLD50_plots_all = grid.arrange(
  grobs = list(pd_hiLD50, lcv_hiLD50,cv_hiLD50, cvl_hiLD50),
  layout_matrix = rbind(c(1, 2),
                        c(3, 4))
)
dev.off()
