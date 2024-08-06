### Test the significance of life stage on the relationship between log mean and log variance of fungal intensity
### Step 1: Load summary file & keep only groups with at least 3 Bd-positive records (for meaningful mean and variance calculations)
### Step 2: Remove unknowns (and try without larval stage)
### Step 3: Run linear regression models with mean fungal load as the predictor variable and variance as the response

library(data.table)
library(lme4)
library(dplyr) 

######################################
#### Step 1: Load and format data ####
######################################

#Load and combine summary datasets
three_dat = fread("data/formatted/summary_aggregation_full.csv")[region != "california" & dataset != "panama"]
sierra_dat = fread("data/formatted/summary_sierra_with_epi_phase.csv")
full_dat = rbind(three_dat,sierra_dat,fill=TRUE)

# Keep only datasets with at least 3 bd-positive records to calculate mean and variance
dat = full_dat[full_dat$num_infected>=3,]

#################################
#### Step 2: Remove unknowns ####
#################################

#Check the lifestage and see if there are any overlapping categories
unique(dat$life_stage)

#Remove unknowns
data = dat[dat$life_stage != "unknown" & dat$life_stage != "",]
#data = dat[dat$life_stage != "unknown" & dat$life_stage != "" & dat$life_stage != "larva",] #Tried with removing larva stage without significant change

#convert variance and mean to numeric
data$var_natural_with_zeros = as.numeric(data$var_natural_with_zeros)
data$mean_natural_with_zeros = as.numeric(data$mean_natural_with_zeros)

############################################
#### Step 3: Run Models and compare AIC ####
############################################

mod = lm(log10(var_natural_with_zeros) ~ log10(mean_natural_with_zeros) , data=data) #Fixed effect: mean; Random effect: NA

mod_ls_slp_rdm = lmer(log10(var_natural_with_zeros) ~ log10(mean_natural_with_zeros) +  (1 + log10(mean_natural_with_zeros) | life_stage), data=data) #Fixed effect: mean, region; Random effect: life_stage (slp)

mod_reg_slp_rdm = lmer(log10(var_natural_with_zeros) ~ log10(mean_natural_with_zeros)  + (1 + log10(mean_natural_with_zeros) | region) , data=data) #Fixed effect: mean, region; Random effect: subregion (slp)

mod_ls_slp_reg_slp_rdm = lmer(log10(var_natural_with_zeros) ~ log10(mean_natural_with_zeros)  + (1 + log10(mean_natural_with_zeros) | region) +  (1 + log10(mean_natural_with_zeros) | life_stage), data=data)  #Fixed effect: mean, region; Random effect: subregion (slp), life_stage (slp)


AIC(mod) #1063.3
AIC(mod_ls_slp_rdm) #828.0
AIC(mod_reg_slp_rdm) #985.8
AIC(mod_ls_slp_reg_slp_rdm) #794.6

AIC(mod_reg_slp_rdm) - AIC(mod_ls_slp_reg_slp_rdm) #191.2


