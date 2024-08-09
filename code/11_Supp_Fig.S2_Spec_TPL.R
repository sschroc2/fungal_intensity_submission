### Reproduce Figure S2 (Histogram of species-specific TPL slopes)
### Step 1: Load summary file & keep only groups with at least 3 Bd-positive records (for meaningful mean and variance calculations)
### Step 2: Examine species list for duplicate species labels and combine and remove unknowns
### Step 3: Run linear regression models with mean fungal load as the predictor variable and variance as the response
### Step 4: Plot

library(data.table)
library(lme4)
library(dplyr)
library(ggplot2)


######################################
#### Step 1: Load and format data ####
######################################

#Load summary dataset
full_dat = fread("data/formatted/analysis_aggregation_dataset.csv")

# Keep only datasets with at least 3 bd-positive records to calculate mean and variance
dat = full_dat[full_dat$num_infected>=3,]

#convert variance and mean to numeric
dat$var_natural_with_zeros = as.numeric(dat$var_natural_with_zeros)
dat$mean_natural_with_zeros = as.numeric(dat$mean_natural_with_zeros)

#####################################################################
#### Step 2: Combine duplicate species names and remove unknowns ####
#####################################################################

#Check the species and see if there are any overlapping categories
unique(dat$species) 

dat$species = ifelse(dat$species=="RACA","rana_catesbeiana",dat$species) #combine RACA and rana_catesbeiana as a single species
data = dat[dat$species!="unknown_species",]
unique(data$species) #36


############################################
#### Step 3: Run Models and compare AIC ####
############################################

mod_null = lm(log10(var_natural_with_zeros) ~ log10(mean_natural_with_zeros) , data=data) #Fixed effect: mean; Random effect: NA

mod_spec_slp_rdm = lmer(log10(var_natural_with_zeros) ~ log10(mean_natural_with_zeros) +  (1 + log10(mean_natural_with_zeros) | species), data=data) #Fixed effect: mean; Random effect: species (slp)

mod_reg_slp_rdm = lmer(log10(var_natural_with_zeros) ~ log10(mean_natural_with_zeros)  + (1 + log10(mean_natural_with_zeros) | region) , data=data) #Fixed effect: mean; Random effect: subregion (slp)

mod_spec_slp_reg_slp_rdm = lmer(log10(var_natural_with_zeros) ~ log10(mean_natural_with_zeros)  + (1 + log10(mean_natural_with_zeros) | region) +  (1 + log10(mean_natural_with_zeros) | species), data=data)  #Fixed effect: mean; Random effect: subregion (slp), species (slp)

AIC(mod_null) #1196.2  #no effect from species or region
AIC(mod_spec_slp_rdm) #1034.2  #species --> random effect on slope
AIC(mod_reg_slp_rdm) #1041.8  #region --> random effect on slope
AIC(mod_spec_slp_reg_slp_rdm) #1028.2  #species --> random effect on slope; region --> random effect on slope

AIC(mod_reg_slp_rdm) - AIC(mod_spec_slp_reg_slp_rdm) #13.56

#Species is a better predictor than null relationship, but it is not a better predictor than region. And it does not benefit the model to include it with region.


spec_coef_slp = coef(mod_spec_slp_reg_slp_rdm)$species


################################
####  Step 4: Plot Fig. S2  ####
################################

spec_slope_histo = ggplot() + geom_histogram(aes(x=spec_coef_slp$`log10(mean_natural_with_zeros)`),bins=20,color="cyan4",fill="azure") + geom_hline(yintercept = 0,color="black") + theme_bw() + xlab("TPL Slope") + ylab("Number of Species") + xlim(1.45,2.3) + geom_rect(aes(xmin=1.48, xmax=1.62, ymin=-Inf, ymax=Inf),fill="pink3",alpha=0.1) + geom_vline(aes(xintercept=1.55),color="red") + geom_vline(aes(xintercept=1.48),color="pink3",linetype="longdash") + geom_vline(aes(xintercept=1.62),color="pink3",linetype="longdash") + annotate("text",x=1.53,y=12,label="Average Macroparasite Slope",angle=90,color="red",size=4.5) + theme(axis.title = element_text(size=14),axis.text = element_text(size = 12))


ggsave("results/plots/spec_slope_histo.png",spec_slope_histo, width=180,height=116,units="mm",dpi=600)

