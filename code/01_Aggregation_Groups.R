### Aggregate records into analysis groups
###
### Step 1: Create functions to calculate Poulin's D/Gini Index & coefficient of variation metrics
### Step 2: Load the cleaned combined datasets
### Step 3: Group records based on dataset, site, year, species, life stage, and temporal group (season) and remove groups with fewer than 10 records
### Step 4: Run loop and create summary dataset
###    - calculates natural mean & variance and log10 mean & variance with and without zero values
###    - fits distributions: 
###       + discretize data by rounding up to a whole number and fit Poisson and Negative Binomial
###       + fit lognormal, Weibull, exponential, & gamma to Bd-positive records
###       + output AIC values to compare fit
###       + fit normal distribution to log-transformed data and run goodness of fit tests
###    - applies function to calculate Poulin's D for each group
### Step 5: Assign dataset names and save summary file
### ***NOTE: Cannot provide raw data. File will not run.

library(dplyr)
library(fitdistrplus)
library(goft)
library(stats)
library(purrr)
library(data.table)
library(lubridate)
library(actuar)
library(ggplot2)
library(KScorrect)
library(nortest)


##########################################################
#### Step 1: Create functions for aggregation metrics ####
##########################################################

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

# poulin_D_other = function(x){
#   # Poulin's D as a measure of aggregation
#   #
#   # Parameters
#   # ---------
#   # x: vector, the distribution of fungal loads
#   # 

#   n = length(x)
#   x = ceiling(x)
#   x = sort(x)

#   num = 0
#   for(i in 1:n){
#     for(j in 1:i){
#       num = num + x[j]
#     }
#   }
#   D = 1 - ((2*num) / (mean(x)*n*(n + 1)))
#   return(D)
# }


coef.var = function(x, drop_zeros=FALSE, log_data=FALSE, log_cv=FALSE){
  # Coefficient of variation as a measure of aggregation (std dev / mean)
  # CV of natural data, CV of log-transformed data, and log(CV) of natural data
  #
  # Parameters
  # ---------
  # x: vector, the distribution of fungal loads

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
  
  if(log_cv == TRUE){
    CV = log10(sd(x)/mean(x))
  } else{
    CV = sd(x)/mean(x)
  }
  
  return(CV)
}


###########################
#### Step 2: Load data ####
###########################

#Split data set into a list of data frames based on criteria

bd = fread("data/cleaned/full_aggregation_data.csv")
bd$bd_load = bd$target_quant_per_swab

# Drop all blank swab ids
bd = bd[bd_swab_id != ""]

########################################
#### Step 3: Split data into groups ####
########################################

split.dat = split(as.data.frame(bd), ~location + region + site + species_capture + life_stage + temporal_group, drop=TRUE)


# Remove sets that have fewer than 10 observations and combine with the non-bd records
dat.list = purrr::discard(split.dat, ~nrow(.) < 10)


#####################################################
#### Step 4: Run loop and create summary dataset ####
#####################################################

# Define continuous distributions to fit
dists = c("gamma", "lnorm", "weibull", "exp")
num_params = c(2, 2, 2, 1)
name_params = list(c("shape", "rate"), c("meanlog", "sdlog"), c("shape", 'scale'), c("rate"))
memp <- function(x, order) mean(x^order)

summary_dat = list()

# Loop through data
for(d in 1:length(dat.list)){
  
  # Extract bd load
  tdat = dat.list[[d]]
  
  my_res = data.frame(site=unique(tdat$site), 
                      location=unique(tdat$location), 
                      region=unique(tdat$region))
  
  ### Calculate summary information for each dataset ###
  
  my_res$species = unique(tdat$species_capture)
  my_res$life_stage = unique(tdat$life_stage)
  my_res$temporal_group = unique(tdat$temporal_group)
  my_res$season = unique(tdat$season)
  my_res$year = unique(tdat$temporal_group)
  my_res$N = nrow(tdat)
  my_res$num_infected = sum(tdat$bd_load > 0)
  my_res$num_infected_discrete = sum(floor(tdat$bd_load) > 0) #Number of infected after discretizing data by rounding down
  my_res$mean_natural_with_zeros = mean(tdat$bd_load)
  my_res$var_natural_with_zeros = var(tdat$bd_load)
  my_res$mean_natural_without_zeros = with(tdat, ifelse(any(bd_load > 0), mean(bd_load[bd_load != 0]), NA))
  my_res$var_natural_without_zeros = with(tdat, ifelse(any(bd_load > 0), var(bd_load[bd_load != 0]), NA))
  my_res$mean_log10_without_zeros = with(tdat, ifelse(any(bd_load > 0), mean(log10(bd_load[bd_load != 0])), NA))
  my_res$var_log10_without_zeros = with(tdat, ifelse(any(bd_load > 0), var(log10(bd_load[bd_load != 0])), NA))
  my_res$mean_discrete_ceiling_with_zeros = mean(ceiling(tdat$bd_load))
  my_res$var_discrete_ceiling_with_zeros = var(ceiling(tdat$bd_load))
  my_res$mean_discrete_floor_with_zeros = mean(floor(tdat$bd_load))
  my_res$var_discrete_floor_with_zeros = var(floor(tdat$bd_load))
  
  ### Fit discrete distributions ###
  
  y_discrete = ceiling(tdat$bd_load)
  
  if(my_res$num_infected >= 3){
    
    fitnbd = tryCatch(
      expr = {
        fitnbd = fitdist(y_discrete, "nbinom", method="mle")
      },
      error = function(e){ 
        # Use method of moments if MLE doesn't work
        fitnbd = fitdist(y_discrete, "nbinom", method="mme")
      }
    )
    
    my_res$nbd_k = fitnbd$estimate[1]
    my_res$nbd_mu = fitnbd$estimate[2]
    my_res$nbd_aic = fitnbd$aic
    
    fitpois = tryCatch(
      expr = {
        fitpois = fitdist(y_discrete, "pois", method="mle")
      },
      error = function(e){ 
        fitpois = fitdist(y_discrete, "pois", method="mme")
      }
    )
    
    my_res$pois_lam = fitpois$estimate[1]
    my_res$pois_aic = fitpois$aic
    
  } else{
    my_res$nbd_k = NA
    my_res$nbd_mu = NA
    my_res$nbd_aic = NA
    my_res$pois_lam = NA
    my_res$pois_aic = NA
  }
  
  ### Fit continuous distributions ###
  
  y = tdat$bd_load
  
  if(my_res$num_infected >= 10){
    
    y_no_0 = y[y != 0]
    
    # Fit continuous distributions
    for(k in 1:length(dists)){
      
      dist = dists[k]
      fitcont = tryCatch(
        expr = {
          fitcont = fitdist(y_no_0, dist, method="mle")
        },
        error = function(e){ 
          # Use method of moments if MLE doesn't work
          fitcont = fitdist(y_no_0, dist, method="mme", order=1:num_params[k], memp=memp)
        }
      )
      
      my_res[[paste0(dist, "_aic")]] = fitcont$aic
      
      # Save the estimated parameters
      for(p in 1:num_params[k]){
        my_res[[paste0(dist, "_", name_params[[k]][p])]] = fitcont$estimate[p]
      }
      
    }
    
    # Calculate goodness of fit of normal distribution to log-transformed data using 4 different tests
    y_no_0_log = log(y_no_0)
    sw = shapiro.test(y_no_0_log)
    my_res[["norm_gof_pvalue_sw"]] = sw$p.value
    my_res[['norm_gof_statistic_sw']] = sw$statistic
    ks = ks.test(y_no_0_log,"pnorm")
    my_res[["norm_gof_pvalue_ks"]] = ks$p.value
    my_res[['norm_gof_statistic_ks']] = ks$statistic
    kslc = LcKS(y_no_0_log,"pnorm")
    my_res[["norm_gof_pvalue_kslc"]] = kslc$p.value
    my_res[['norm_gof_statistic_kslc']] = kslc$D.obs
    ad = ad.test(y_no_0_log)
    my_res[["norm_gof_pvalue_ad"]] = ad$p.value
    my_res[['norm_gof_statistic_ad']] = ad$statistic
    
    
  } else{
    
    # Assign NAs
    for(k in 1:length(dists)){
      
      dist = dists[k]
      my_res[[paste0(dist, "_aic")]] = NA 
      
      for(p in 1:num_params[k]){
        my_res[[paste0(dist, "_", name_params[[k]][p])]] = NA
      }
      
    }
    
    my_res[["norm_gof_pvalue_sw"]] = NA
    my_res[['norm_gof_statistic_sw']] = NA
    my_res[["norm_gof_pvalue_ks"]] = NA
    my_res[['norm_gof_statistic_ks']] = NA
    my_res[["norm_gof_pvalue_kslc"]] = NA
    my_res[['norm_gof_statistic_kslc']] = NA
    my_res[["norm_gof_pvalue_ad"]] = NA
    my_res[['norm_gof_statistic_ad']] = NA
    
    
  }
  
  # Add in Poulin's D/Gini index as a measure of aggregation 
  my_res[["poulin_D"]] = poulin_D(y, drop_zeros=TRUE, log_data=FALSE)
  my_res[["poulin_D_log"]] = poulin_D(y, drop_zeros=TRUE, log_data=TRUE)
  
  # Add coefficient of variation
  # my_res[["cv_log"]]=sqrt(var(y_no_0_log)/mean(y_no_0_log))
  # my_res[["log_cv"]]=log10(sqrt(var(y_no_0)/mean(y_no_0)))
  # my_res[["cv"]]=sqrt(var(y_no_0)/mean(y_no_0))
  
  # Add coefficient of variation
  my_res[["cv_log"]]=coef.var(y, drop_zeros=TRUE, log_data=TRUE, log_cv = FALSE)
  my_res[["log_cv"]]=coef.var(y, drop_zeros=TRUE, log_data=FALSE, log_cv = TRUE)
  my_res[["cv"]]=coef.var(y, drop_zeros=TRUE, log_data=FALSE, log_cv = FALSE)

  
  summary_dat[[d]] = my_res
  
}

# Combine dataset
summary_aggregation_data = data.table(do.call(rbind, summary_dat))


######################################
#### Step 5: Assign dataset names ####
######################################

# Relabel for plotting
summary_aggregation_data$dataset = NA
summary_aggregation_data$dataset[summary_aggregation_data$location == "panama"] = "panama"
summary_aggregation_data$dataset[summary_aggregation_data$location == "brazil"] = "brazil"
summary_aggregation_data$dataset[summary_aggregation_data$region == "california"] = "sierra_nevada"
summary_aggregation_data$dataset[summary_aggregation_data$region == "east_bay"] = "east_bay"
summary_aggregation_data$dataset[summary_aggregation_data$region %in% c("vermont", "new_mexico", "pennsylvania", "tennessee", "louisiana")] = "serdp"


#Save summary dataset
fwrite(summary_aggregation_data, "data/formatted/summary_aggregation_data.csv")




