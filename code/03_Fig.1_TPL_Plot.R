### Reproduce Figure 1 (Slope of Taylor's Power Law for each dataset)
### Step 1: Load summary file & keep only groups with at least 3 Bd-positive records (for meaningful mean and variance calculations)
### Step 2: Set function to specify number of decimal places to include (for plot aesthetics)
### Step 3: Run linear regression models with mean fungal load as the predictor variable and variance as the response
### Step 4: Create plot labels by pulling from model output
### Step 5: Plot 
### Optional --> Run based on discrete data (use var_discrete_ceiling_with_zeros & mean_discrete_ceiling_with_zeros)


library(data.table)
library(lme4)
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


##################################################################
#### Step 2: Create function to specify decimal place display ####
##################################################################

specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

##########################################
#### Step 3: Linear regression models ####
##########################################

#Run linear models for each dataset
mod_eastbay = lm(log10(var_natural_with_zeros) ~ log10(mean_natural_with_zeros), data = dat[dat$dataset=="east_bay",])
mod_serdp = lm(log10(var_natural_with_zeros) ~ log10(mean_natural_with_zeros), data = dat[dat$dataset=="serdp",])
mod_sierra = lm(log10(var_natural_with_zeros) ~ log10(mean_natural_with_zeros), data = dat[dat$dataset=="sierra_nevada",])
mod_brazil = lm(log10(var_natural_with_zeros) ~ log10(mean_natural_with_zeros), data = dat[dat$dataset=="brazil",])

####################################
#### Step 4: Create plot labels ####
####################################

#Assign legend text, including confidence interval
l_eastbay = paste("Eastbay=",specify_decimal(coef(mod_eastbay)[2],2),"  [",specify_decimal(confint(mod_eastbay)[2],2),"-",specify_decimal(confint(mod_eastbay)[4],2),"]")
l_serdp = paste("Eastern US=",specify_decimal(coef(mod_serdp)[2],2),"  [",specify_decimal(confint(mod_serdp)[2],2),"-",specify_decimal(confint(mod_serdp)[4],2),"]")
l_sierra = paste("Sierra=",specify_decimal(coef(mod_sierra)[2],2),"  [", specify_decimal(confint(mod_sierra)[2],2),"-",specify_decimal(confint(mod_sierra)[4],2),"]")
l_brazil = paste("Brazil=",specify_decimal(coef(mod_brazil)[2],2),"  [", specify_decimal(confint(mod_brazil)[2],2),"-",specify_decimal(confint(mod_brazil)[4],2),"]")
#l_macro = "Macroparasites= 1.55 [ 1.48 - 1.62 ]"

#join the legend text
stat_lab = rbind(l_brazil,l_eastbay,l_serdp,l_sierra)

legend_lines <- c("Lognormal = 2" = "dotted", "Poisson = 1" = "solid", "Macroparasites = 1.55 [ 1.48 - 1.62 ]" = "dashed")

###############################
#### Step 5: Plot Figure 1 ####
###############################

tpl_plot = ggplot()  + 
  geom_abline(aes(intercept = 1.098, slope =1.551,linetype = "A"), color="black", linewidth=1) + 
  geom_abline(aes(intercept = 0, slope =1,linetype = "B"), color="black", linewidth=1) + 
  geom_abline(aes(intercept = 0, slope =2,linetype = "C"), color="black", linewidth=1) + 
  geom_point(alpha=0.35,size=2.5,shape=21, aes(x=log10(mean_natural_with_zeros), y=log10(var_natural_with_zeros), color=dataset, fill=dataset),data=dat) + 
  geom_smooth(aes(x=log10(mean_natural_with_zeros),y=log10(var_natural_with_zeros), color=dataset),method="lm",se=TRUE,linewidth=1.5,data=dat) + 
  scale_linetype_manual(name="",values = c("solid","dotted","dashed"),labels = c("Macroparasites = 1.55 [ 1.48 - 1.62 ]", "Poisson = 1", "Lognormal = 2")) + 
  scale_fill_manual(name="",values = c("lightskyblue1", "palegreen", "goldenrod1", "plum1"),labels = c(stat_lab[1:nrow(stat_lab)])) +  
  scale_color_manual(name="",values = c("dodgerblue", "seagreen", "darkorange", "darkorchid"),labels = c(stat_lab[1:nrow(stat_lab)])) + guides(color=guide_legend(override.aes=list(fill="lightgray"))) + 
  xlim(-2,8) +  ylim(-5,17.5) + xlab(expression(paste(Log[10],"(Mean Bd Load)", sep=""))) + ylab(expression(paste(Log[10],"(Variance Bd Load)", sep=""))) + theme_bw()  +   theme( legend.text=element_text(size=10), legend.title=element_text(size=14), title = element_text(size=18), axis.title = element_text(size=14),legend.spacing.y = unit(0.1, "cm")) + guides(linetype=guide_legend(title="TPL Slopes"))

### Options for making "null" intercepts fit the data more closely
# geom_abline(aes(intercept = 1.5, slope =1,linetype = "B"), color="black", linewidth=0.5) + 
# geom_abline(aes(intercept = 0.8, slope =2,linetype = "C"), color="black", linewidth=0.5)

#Save plot
ggsave("results/plots/tpl_lin_reg.png",tpl_plot, width=180,height=98,units="mm",dpi=600)


