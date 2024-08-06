### Check the goodness of fit (GOF) of the normal distribution on the log-transformed data & recreate Fig. S1
### Step 1: Load aggregation dataset
### Step 2: Adjust p-values for multiple tests
### Step 3: Check GOF by life stage
### Step 4: Check GOF by epizoological phase
### Step 5: Check GOF by dataset
### Step 6: Plot Fig. S1

library(ggplot2)

###########################
#### Step 1: Load data ####
###########################

all_data = fread("data/formatted/analysis_aggregation_dataset.csv")

####################################################
#### Step 2: Adjust p-values for multiple tests ####
####################################################

all_data$p.adj = p.adjust(all_data$norm_gof_pvalue_sw,method="fdr")
#Check the percentage that are not statistically different from a normal distribution
1-(nrow(all_data[all_data$p.adj<=0.05,])/nrow(all_data)) #96.7%

#########################################
#### Step 3: Check GOF by life stage ####
#########################################

all_dat_adult = all_data[all_data$life_stage=="adult",]
all_dat_larva = all_data[all_data$life_stage=="larva",]
all_dat_subadult = all_data[all_data$life_stage=="subadult",]

all_dat_adult$p_adj = p.adjust(all_dat_adult$norm_gof_pvalue_sw,method="fdr")
all_dat_larva$p_adj = p.adjust(all_dat_larva$norm_gof_pvalue_sw,method="fdr")
all_dat_subadult$p_adj = p.adjust(all_dat_subadult$norm_gof_pvalue_sw,method="fdr")


ggplot() + geom_boxplot(aes(x=all_dat_adult$life_stage,y=all_dat_adult$p_adj)) + geom_boxplot(aes(x=all_dat_larva$life_stage,y=all_dat_larva$p_adj)) + geom_boxplot(aes(x=all_dat_subadult$life_stage,y=all_dat_subadult$p_adj)) + geom_abline(intercept = 0.05,slope=0,color="red")+ xlab("Life Stage") + ylab("Goodness of Fit (p-value)") + ggtitle("Eastbay: Log-Normal Goodness of Fit: Shapiro-Wilk")


nrow(all_dat_adult[all_dat_adult$p_adj<=0.05,])/nrow(all_dat_adult) #3.7%
nrow(all_dat_larva[all_dat_larva$p_adj<=0.05,])/nrow(all_dat_larva) #6.2%
nrow(all_dat_subadult[all_dat_subadult$p_adj<=0.05,])/nrow(all_dat_subadult) #1.3%

########################################
#### Step 4: Check GOF by epi phase ####
########################################

dat_enzo = all_data[all_data$phase=="enzootic\n(Northern Sierra)",]
dat_inv = all_data[all_data$phase=="invasion\n(Southern Sierra)",]
dat_decl = all_data[all_data$phase=="post-invasion\n(Southern Sierra)",]

dat_enzo$p_adj = p.adjust(dat_enzo$norm_gof_pvalue_sw,method="fdr")
dat_inv$p_adj = p.adjust(dat_inv$norm_gof_pvalue_sw,method="fdr")
dat_decl$p_adj = p.adjust(dat_decl$norm_gof_pvalue_sw,method="fdr")

ggplot() + geom_boxplot(aes(x=dat_enzo$phase,y=dat_enzo$p_adj)) + geom_boxplot(aes(x=dat_inv$phase,y=dat_inv$p_adj)) + geom_boxplot(aes(x=dat_decl$phase,y=dat_decl$p_adj)) + geom_abline(intercept = 0.05,slope=0,color="red")+ xlab("Life Stage") + ylab("Goodness of Fit (p-value)") 

nrow(dat_enzo[dat_enzo$p_adj<=0.05,])/nrow(dat_enzo) #4.2%
nrow(dat_inv[dat_inv$p_adj<=0.05,])/nrow(dat_inv) #0.4% | 1.1% (change in prev cutoff? [50%])
nrow(dat_decl[dat_decl$p_adj<=0.05,])/nrow(dat_decl) #19.4% | 20.8% (change in prev cutoff? [50%])

#########################################
#### Step 5: Check GOF by dataset ####
#########################################

dat_distr = all_data[all_data$num_infected>9,]

dat_eastbay = dat_distr[dat_distr$dataset=="east_bay",]
dat_serdp = dat_distr[dat_distr$dataset=="serdp",]
dat_sierra = dat_distr[dat_distr$dataset=="sierra_nevada",]
dat_brazil = dat_distr[dat_distr$dataset=="brazil",]

dat_eastbay$p_adj = p.adjust(dat_eastbay$norm_gof_pvalue_sw,method="fdr")
dat_serdp$p_adj = p.adjust(dat_serdp$norm_gof_pvalue_sw,method="fdr")
dat_sierra$p_adj = p.adjust(dat_sierra$norm_gof_pvalue_sw,method="fdr")
dat_brazil$p_adj = p.adjust(dat_brazil$norm_gof_pvalue_sw,method="fdr")

nrow(dat_eastbay[dat_eastbay$p_adj<=0.05,])/nrow(dat_eastbay) #4.9% 
nrow(dat_serdp[dat_serdp$p_adj<=0.05,])/nrow(dat_serdp) #0.0% 
nrow(dat_sierra[dat_sierra$p_adj<=0.05,])/nrow(dat_sierra) #17.0% 
nrow(dat_brazil[dat_brazil$p_adj<=0.05,])/nrow(dat_brazil) #20.0% 

#########################
#### Step 6: Plot S1 ####
#########################

gof_plot = ggplot() + geom_boxplot(aes(x=dat_eastbay$dataset,y=dat_eastbay$p_adj),color="seagreen") + geom_boxplot(aes(x=dat_serdp$dataset,y=dat_serdp$p_adj),color="darkorange") + geom_boxplot(aes(x=dat_sierra$dataset,y=dat_sierra$p_adj),color="darkorchid") + geom_boxplot(aes(x=dat_brazil$dataset,y=dat_brazil$p_adj),color="dodgerblue") + geom_abline(intercept = 0.05,slope=0,color="red")+ xlab("Dataset") + ylab("Goodness of Fit (adjusted p-value)") + theme_bw() + theme( axis.title = element_text(size=14),axis.text.x=element_text(size=12))+scale_x_discrete(labels=c("east_bay" = "Eastbay", "serdp" = "Eastern US","sierra_nevada" = "Sierra", "brazil" = "Brazil")) + geom_text(aes(y=0.06, label="\np=0.05", x=0.6), colour="red", angle=0)

ggsave("results/plots/gof_plot.png",gof_plot, width=180,height=116,units="mm",dpi=600)


