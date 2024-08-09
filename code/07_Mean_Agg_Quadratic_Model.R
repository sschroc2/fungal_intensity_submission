### Test the relationship between mean fungal load and Poulin's D
### Step 1: Load Sierra summary file and remove records with fewer than 2 Bd-positive and larval groups (for consistent comparisons)
### Step 2: Observe if there is a significant quadratic effect using frequentist (a) and Bayesian (b) approaches
### Step 3: Observe is a there is a significant effect of life stage (Adult and subadult)
### Step 4: Observe is a there is a significant effect of location (Northern vs Southern Sierra)


library(brms)
library(data.table)
library(ggplot2)
library(betareg)

######################################
#### Step 1: Load and format data ####
######################################

sierra_dat = fread("data/formatted/analysis_aggregation_dataset.csv")[dataset == "sierra_nevada"]

# Beta regression 1: Does aggregation vary with epizoological stage?
sierra_dat[, "log_mean":=mean_log10_without_zeros]

# Fit varying models to see if life stage or location affect the relationship
fit_dat = sierra_dat[life_stage != "larva" & num_infected >= 2 & poulin_D > 0]


###################################################################################
#### Step 2a: Compare model with quadratic effect to only linear (Frequentist) ####
###################################################################################

lin_mod = betareg(poulin_D ~ log_mean, data=fit_dat)
quad_mod = betareg(poulin_D ~ I(log_mean^2) + log_mean, data=fit_dat)

AIC(lin_mod) #-125.16
AIC(quad_mod) #-275.29
delta_AIC = AIC(lin_mod) - AIC(quad_mod) #150.13


ggplot() + geom_point(aes(x=log_mean,y=poulin_D,color=phase),data=fit_dat)  +
  geom_line(aes(y = predict(lin_mod, fit_dat),x=fit_dat$log_mean)) +
  geom_line(aes(y = predict(quad_mod, fit_dat),x=fit_dat$log_mean))

#################################################################################
#### Step 2b: Compare model with quadratic effect to only linear (Bayesian) #####
#################################################################################

mod_fxn = function(eq){
  brm(bf(eq,
         phi ~ 1),
      data=fit_dat,
      family=Beta(),
      prior=c(set_prior("cauchy(2, 2.5)", class="Intercept", dpar="phi")),
      chains=4, iter=2000, warmup=1000, cores=4)
}


mod_null = mod_fxn(poulin_D ~ log_mean)
mod_quad = mod_fxn(poulin_D ~ log_mean + I(log_mean^2))


# Perform model comparison
diff1 = loo(mod_null, mod_quad)

# Model comparisons:
#   elpd_diff se_diff
# mod_quad   0.0       0.0  
# mod_null -74.6      11.7


# Get posterior predictions overall
newdata = fit_dat[,c("poulin_D","log_mean")]
pred = posterior_epred(mod_quad, newdata=newdata, re.form=NA)
newdata$med = apply(pred, 2, function(x) quantile(x, .50))
newdata$lower = apply(pred, 2, function(x) quantile(x, .025))
newdata$upper = apply(pred, 2, function(x) quantile(x, .975))

ggplot() + geom_point(aes(x=log_mean, y=poulin_D),data=fit_dat) + 
  geom_line(data=newdata, aes(x=log_mean, y=med)) +
  geom_ribbon(data=newdata, aes(x=log_mean, ymin=lower, ymax=upper), alpha=0.5) + ylim(c(0, 1)) + ylab("Poulin's D") + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) +  theme_bw()


###################################################
#### Step 3: Observe the effect of life stage #####
###################################################

mod_quad_ls = mod_fxn(poulin_D ~ life_stage*(log_mean + I(log_mean^2)))

# Perform model comparison
diff2 = loo(mod_quad, mod_quad_ls)

# Model comparisons:
#   elpd_diff se_diff
# mod_quad_ls  0.0       0.0   
# mod_quad    -1.1       3.3 

### No significant effect

# Get posterior predictions for lifestage
newdata = fit_dat[,c("poulin_D","log_mean","life_stage")]
pred = posterior_epred(mod_quad_ls, newdata=newdata, re.form=NA)
newdata$med = apply(pred, 2, function(x) quantile(x, .50))
newdata$lower = apply(pred, 2, function(x) quantile(x, .025))
newdata$upper = apply(pred, 2, function(x) quantile(x, .975))

ggplot() + geom_point(aes(x=log_mean, y=poulin_D,colour = life_stage),data=fit_dat) + 
  geom_line(data=newdata, aes(x=log_mean, y=med,colour = life_stage)) +
  geom_ribbon(data=newdata, aes(x=log_mean, ymin=lower, ymax=upper,colour = life_stage,fill=life_stage), alpha=0.5) + ylim(c(0, 1)) + ylab("Poulin's D") + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) +  theme_bw()


#################################################
#### Step 4: Observe the effect of location #####
#################################################

mod_quad_loc = mod_fxn(poulin_D ~ location*(log_mean + I(log_mean^2))) 

# Perform model comparison
diff3 = loo(mod_quad, mod_quad_loc)

# Model comparisons:
#   elpd_diff se_diff
# mod_quad_loc  0.0       0.0   
# mod_quad     -0.4       2.8 

### No significant effect

# Get posterior predictions for location
newdata = fit_dat[,c("poulin_D","log_mean","location")]
pred = posterior_epred(mod_quad_loc, newdata=newdata, re.form=NA)
newdata$med = apply(pred, 2, function(x) quantile(x, .50))
newdata$lower = apply(pred, 2, function(x) quantile(x, .025))
newdata$upper = apply(pred, 2, function(x) quantile(x, .975))

ggplot() + geom_point(aes(x=log_mean, y=poulin_D,colour = location),data=fit_dat) + 
  geom_line(data=newdata, aes(x=log_mean, y=med,colour = location)) +
  geom_ribbon(data=newdata, aes(x=log_mean, ymin=lower, ymax=upper,colour = location,fill=location), alpha=0.5) + ylim(c(0, 1)) + ylab("Poulin's D") + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) +  theme_bw()


