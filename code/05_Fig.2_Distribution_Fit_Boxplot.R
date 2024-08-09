### Recreate Fig. 2 and calculate statistics of the best fits
### Step 1: Load aggregation dataset
### Step 2: Subset data by those populations with at least 10 infected individuals and values for all fitted distributions 
### Step 3: Calculate the difference between AICs of the fit distributions and rank which distributions had the lowest value
### Step 4: Create datasets of delta AICs and distribution rank
### Step 5: Plot Fig. 2 (Boxplot of delta AICs for each distribution)
### Step 6: Calculate proportions of the datasets that fit particular distributions best
### Step 7: Similar to 6, but includes multiple best fits (within 2 AIC of lowest value)

library(ggplot2)

###########################
#### Step 1: Load data ####
###########################

full_dat = fread("data/formatted/analysis_aggregation_dataset.csv")
dat = full_dat[full_dat$num_infected>=3,]


#############################
#### Step 2: Subset data ####
#############################

#Populations with at least 10 infected for distribution analyses
dat_distr = dat[dat$num_infected>9,]

#List vars and labels
aic_cont <- c( "gamma_aic", "lnorm_aic", "weibull_aic", "exp_aic") #Untransformed data 
distr_labs = c("gamma","log normal","weibull","exponential")

df.subset = subset(dat_distr,select=c("num_infected","dataset","gamma_aic", "lnorm_aic", "weibull_aic", "exp_aic")) #subset only vars needed
cmplt = df.subset[complete.cases(df.subset)] #include only those records with values (1000 out of 1001)

#################################################################################
#### Step 3: Find delta AIC for each distribution and rank based on best fit ####
#################################################################################

delta_aic = cmplt[, 3:(2+length(aic_cont))] - apply(cmplt[, 3:(2+length(aic_cont))], 1, min) #find difference in aic of lowest aic value
ranked_aic = data.table(t(apply(delta_aic, 1, rank))) #rank the lowest aic 


########################################################
#### Step 4: Create datasets of delta AICs and rank ####
########################################################

#Create dataframes of datasets and distribution rankings
delta_aic[, c("num_infected", "dataset"):=.(cmplt$num_infected, cmplt$dataset)]
ranked_aic[, c("num_infected", "dataset"):=.(cmplt$num_infected, cmplt$dataset)]
melted_aic = melt(delta_aic, id.vars=c("num_infected", "dataset"))
melted_aic = data.table(melted_aic)
melted_rank = melt(ranked_aic, id.vars=c("num_infected","dataset"))
melted_rank = data.table(melted_rank)

melted_aic$variable <- ordered(melted_aic$variable,
                               levels = c("gamma_aic","lnorm_aic","weibull_aic","exp_aic"),
                               labels = c("gamma", "log normal", "weibull","exponential"))

############################
#### Step 5: Plot Fig 2 ####
############################

distr_plot = ggplot(melted_aic) + geom_boxplot(aes(x=variable, y=log10(value + 1),color=melted_aic$dataset))  +   
  ylab(expression(paste(Log[10],"(",Delta," AIC + 1)", sep=""))) + xlab("Distribution") + labs(color="Dataset")  + theme_bw() +  
  theme(axis.text.x = element_text(angle=45, hjust=1, size=12),axis.title = element_text(size=14),legend.text = element_text(size=12),legend.title = element_text(size=12)) + scale_color_manual(values=c("dodgerblue", "seagreen", "darkorange", "darkorchid"),labels=c( "brazil" = "Brazil","east_bay" = "Eastbay", "serdp" = "SERDP","sierra_nevada" = "Sierra")) + scale_x_discrete(labels=c("log normal" = "Lognormal","weibull" = "Weibull","exponential"="Exponential","gamma" = "Gamma"))

ggsave("results/plots/distr_plot.png",distr_plot, width=180,height=115,units="mm",dpi=600)

#Add data points to plot 
# +  geom_point(aes(x=variable, y=log10(value + 1),color=melted_aic$dataset),alpha=0.5,position=position_dodge(width=0.75))

##############################################################
#### Step 6: Find the percentage of best fit distribution ####
##############################################################

#Identify which distribution is the best fit
delta_aic_update = delta_aic
delta_aic_update$best_fit = ifelse(delta_aic_update$lnorm_aic==0,"lnorm",ifelse(delta_aic_update$weibull_aic==0,"weibull",ifelse(delta_aic_update$gamma_aic==0,"gamma","exp")))

#Create dataframe with proportions of overall best fit 
best.fit = as.data.frame(delta_aic_update %>%
                           count(best_fit))
total = as.data.frame(delta_aic_update %>%
                      summarize(n_tot=n()))

combined = merge(best.fit,total)
combined$perc = (combined$n/combined$n_tot)*100
combined


#Create dataframe with proportions of best fit distributions
best_fit = as.data.frame(delta_aic_update %>%
                           group_by(dataset) %>%
                           count(best_fit))

tot = as.data.frame(delta_aic_update %>%
                      group_by(dataset) %>%
                      summarize(n_tot=n()))

cmb = merge(best_fit,tot)
cmb$perc = (cmb$n/cmb$n_tot)*100
cmb


##########################################################################
#### Step 7: Find the percentage of best fit distribution (w/i 2 AIC) ####
##########################################################################

delta_aic_update$exp_under2 = ifelse(delta_aic_update$exp_aic <= 2, 1, 0)
delta_aic_update$gamma_under2 = ifelse(delta_aic_update$gamma_aic <= 2, 1, 0)
delta_aic_update$lnorm_under2 = ifelse(delta_aic_update$lnorm_aic <= 2, 1, 0)
delta_aic_update$weibull_under2 = ifelse(delta_aic_update$weibull_aic <= 2, 1, 0)

delta_aic_update$exp_cat = ifelse(delta_aic_update$exp_aic <= 2, "E", "0")
delta_aic_update$gamma_cat = ifelse(delta_aic_update$gamma_aic <= 2, "G", "0")
delta_aic_update$lnorm_cat = ifelse(delta_aic_update$lnorm_aic <= 2, "L", "0")
delta_aic_update$weibull_cat = ifelse(delta_aic_update$weibull_aic <= 2, "W", "0")
delta_aic_update$all_fit = paste0(delta_aic_update$exp_cat,delta_aic_update$gamma_cat,delta_aic_update$lnorm_cat,delta_aic_update$weibull_cat)

delta_aic_update$lnorm_best_fit = ifelse(delta_aic_update$best_fit!="lnorm",NA,ifelse(delta_aic_update$exp_under2==0 & delta_aic_update$gamma_under2==0 & delta_aic_update$weibull_under2==0,1,0))

delta_aic_update$exp_best_fit = ifelse(delta_aic_update$best_fit!="exp",NA,ifelse(delta_aic_update$lnorm_under2==0 & delta_aic_update$gamma_under2==0 & delta_aic_update$weibull_under2==0,1,0))

delta_aic_update$weibull_best_fit = ifelse(delta_aic_update$best_fit!="weibull",NA,ifelse(delta_aic_update$exp_under2==0 & delta_aic_update$gamma_under2==0 & delta_aic_update$lnorm_under2==0,1,0))

delta_aic_update$gamma_best_fit = ifelse(delta_aic_update$best_fit!="gamma",NA,ifelse(delta_aic_update$exp_under2==0 & delta_aic_update$lnorm_under2==0 & delta_aic_update$weibull_under2==0,1,0))

#Proportion of groups that were best described by a single distribution
ln_best = nrow(delta_aic_update[delta_aic_update$lnorm_best_fit==1,])/nrow(delta_aic_update) #38.0%
exp_best = nrow(delta_aic_update[delta_aic_update$exp_best_fit==1,])/nrow(delta_aic_update) #0%
wb_best = nrow(delta_aic_update[delta_aic_update$weibull_best_fit==1,])/nrow(delta_aic_update) #4.4%
gam_best = nrow(delta_aic_update[delta_aic_update$gamma_best_fit==1,])/nrow(delta_aic_update) #).4%

#Proportion of each distribution that is within 2 AIC points of the lowest value
tot = nrow(delta_aic_update)
exp_perc = sum(delta_aic_update$exp_under2)/tot
gamma_perc = sum(delta_aic_update$gamma_under2)/tot
lnorm_perc = sum(delta_aic_update$lnorm_under2)/tot
weibull_perc = sum(delta_aic_update$weibull_under2)/tot

exp_perc #25.7%
gamma_perc #35.2%
lnorm_perc #76.7%
weibull_perc #58.8%

1-(ln_best + exp_best + wb_best + gam_best) #57.3

#Percent based on best fit category
diff_under2 = as.data.frame(delta_aic_update %>%
                              group_by(best_fit) %>%
                              summarise(tot = n(),exp = sum(exp_under2),gamma = sum(gamma_under2),lnorm = sum(lnorm_under2),weibull = sum(weibull_under2)))
diff_under2$exp_perc = diff_under2$exp/diff_under2$tot
diff_under2$gamma_perc = diff_under2$gamma/diff_under2$tot
diff_under2$lnorm_perc = diff_under2$lnorm/diff_under2$tot
diff_under2$weibull_perc = diff_under2$weibull/diff_under2$tot

diff_under2


#Interpretation: 66% of those with best fit of exponential also fit gamma within 2 AIC points; 37% fit log normal; 100% fit weibull
#Interpretation: 89.5% of those with best fit of gamma also fit weibull within 2 AIC points; 26% fit exponential; 31.6% fit log normal
#Interpretation: lnorm distribution with the fewest other fits (5.5% exp, 13.3% gamma, 32% weibull)
#Interpretation: 61% of those with weibull best fit also fit gamma and 49% lnorm and 21% exponential

#Proportion that fit multiple distributions
### Letters represent all distributions with =< 2 AIC of best fit
### Ex. EGLW: all distributons; 0G00: only Gamma; 0GL0: Gamma and lognormal
all_fit = as.data.frame(delta_aic_update %>%
                          group_by(all_fit) %>%
                          summarize(n=n()))
all_fit$total = sum(all_fit$n)
all_fit$perc = all_fit$n/all_fit$total

all_fit

#Interpretation: 38% of populations only fit log normal; 13.5% fit lnorm and weibull; 11.1% fit all except exp; 9.7% fit all; 7.6% fit all except log norm 

