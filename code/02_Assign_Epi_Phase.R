### Calculate the epizoological states for the Sierra Nevada data
###
### Step 1: Set function to calculate Poulin's D
### Step 2: Load the summarized data
### Step 3: For each site by year, look at the prevalence data and abundance data
###  to decide whether a site is in a invasion, epizootic/decline, or post-decline state.
###  Invasion: Prevalence <= cp n_inf > 5
###  Epizootic/decline: Prevalence >= cutoff_decline, but abundance is not less than 25% of previous max
###  Post-decline: Prevalence > cutoff_decline and abundance < max 
### ***NOTE: Cannot provide raw data. File will not run.

library(data.table)
library(ggplot2)
library(patchwork)
library(brms)
library(rstan)
library(ggnewscale)

#############################################
#### Step 1: Set function for Poulin's D ####
#############################################

poulin_D = function(x, drop_zeros=FALSE, log_data=FALSE){
  # Poulin's D as a measure of aggregation
  #
  # Parameters
  # ---------
  # x: vector, the distribution of fungal loads
  # 
  # Notes
  # -----
  # From McVinish and Lester 2021. Slightly different equation then in Poulin
  
  
  if(drop_zeros){
    x = x[x > 0]
  }
  
  n = length(x)
  
  if(n == 0){
    return(NA)
  }
  
  if(log_data == TRUE){
    x = log10(x)
  }
  
  num = 0
  for(i in 1:n){
    for(j in 1:n){
      num = num + abs(x[i] - x[j])
    }
  }
  D = num / (mean(x)*2*n^2)
  return(D)
}

######################################
#### Step 2: Load and format data ####
######################################

sierra_dat = fread("data/formatted/summary_aggregation_data.csv")[region == "california"]

# Magnitude of decline data
decline_dat = fread("data/Archival/SierraNevadaFieldData/analysisII_magnitudeofdecline_withoutsubadults.csv")

# We know these sites have experienced invasion and decline. Helpful to just look at this sites 
decline_sites = decline_dat$site_id

# Load the full sierra data WITH VES. Drop tadpole data.  Using this to compare answers
full_dat = fread("data/Archival/SierraNevadaFieldData/combined_Bd_data_for_analysis_allVES.csv")[capture_life_stage != "tadpole"]
abund_and_prev = full_dat[, 
                          .(bd_present=ifelse(all(is.na(bd_load)), as.integer(-1), as.integer(any(bd_load > 0))), 
                            num_infected=ifelse(all(is.na(bd_load)), as.integer(0), sum(bd_load > 0)), 
                            num_samps=ifelse(all(is.na(bd_load)), as.integer(0), length(bd_load)),
                            adult_abundance=ifelse(all(is.na(adult)), as.integer(-1), max(adult, na.rm=T)),
                            prev=ifelse(all(is.na(bd_load)), as.numeric(-1), sum(bd_load > 0) / length(bd_load)),
                            mean_bd=ifelse(all(is.na(bd_load)), as.numeric(-1), ifelse(any(bd_load > 0), mean(log10(bd_load[bd_load > 0])), -1)),
                            poulin_D=ifelse(all(is.na(bd_load)), as.numeric(-1), ifelse(any(bd_load > 0), poulin_D(bd_load, drop_zeros=TRUE, log_data=FALSE), -1))
                          ), by=.(site_id, year, capture_life_stage)][year > 2004]


fwrite(abund_and_prev, "data/formatted/abund_and_prev.csv")

########################################
#### Step 3: Assign invasion stages ####
########################################

# In Wilber et al. 2022, Journal of Animal Ecology, the "invasion phase" is defined with prevalence is
# less than rho, with rho being between 0.25 and 0.5
#
# From Knapp et al. 2016, we know all Yosemite sites are generally enzootic. These are sites that
# start with a "7".

p_cut = 0.5 # Use the upper prevalence criteria established in Wilber et al. 2022. Lower criteria is 0.25.  Try both
#p_cut = 0.25

# Assign stages to summarizes data
sierra_dat[, "prev":=num_infected / N]
sierra_dat$phase = character(nrow(sierra_dat))
sierra_dat[prev <= p_cut, "phase":="invasion\n(Southern Sierra)"]
sierra_dat[prev > p_cut, "phase":="post-invasion\n(Southern Sierra)"]
sierra_dat[site %like% "7....", "phase":="enzootic\n(Northern Sierra)"]

sierra_dat$site_id = sierra_dat$site
sierra_dat$var_natural_with_zeros = as.numeric(sierra_dat$var_natural_with_zeros)
sierra_dat$var_natural_without_zeros = as.numeric(sierra_dat$var_natural_without_zeros)

# Also split by location: Yosemite and everything else
sierra_dat[, "location":=ifelse(site %like% "7....", "Northern Sierra", "Southern Sierra")]


fwrite(sierra_dat, "data/formatted/summary_sierra_with_epi_phase.csv")


#Load and combine summary datasets
three_datasets = fread("data/formatted/summary_aggregation_data.csv")[region != "california" & dataset != "panama"]
analysis_aggregation_dataset = rbind(three_datasets,sierra_dat,fill=TRUE)

fwrite(analysis_aggregation_dataset, "data/formatted/analysis_aggregation_dataset.csv") #This dataset will be used for the analyses

