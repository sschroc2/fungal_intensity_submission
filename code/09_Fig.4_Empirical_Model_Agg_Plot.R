### Reproduce Figure 4 (Mean-aggregation plot for empirical data and model predictions)
### Step 1: Load Sierra summary file and remove records with fewer than 2 Bd-positive and larval groups (for consistent comparisons)
### Step 2: Plot empirical by log mean fungal load against aggregation metrics, grouped by epizoological stage
### Step 3: Plot the results of the integral projection model
### Step 4: Combine plots

library(ggplot2)
library(gridExtra)
library(cowplot)
library(grid)
library(data.table)


######################################
#### Step 1: Load and format data ####
######################################

sierra_dat = fread("data/formatted/analysis_aggregation_dataset.csv")[dataset == "sierra_nevada"]
sierra_dat[, "log_mean":=mean_log10_without_zeros]
sierra_no_larva = sierra_dat[life_stage != "larva" & num_infected >= 2 & poulin_D > 0]

######################################
#### Step 2: Plot empirical data  ####
######################################

#Reorder levels of this factor to have chronological order in the legend
sierra_no_larva$phase<- factor(sierra_no_larva$phase, levels = c("invasion\n(Southern Sierra)" , "post-invasion\n(Southern Sierra)", "enzootic\n(Northern Sierra)"))

pd=ggplot(data=sierra_no_larva) + geom_point(size=2.5,aes(x=log_mean, y=poulin_D, color=phase)) + scale_color_manual(name="Phase",values = c("goldenrod3",  "darkturquoise","maroon4"),labels = c("Invasion \n(Southern Sierra)","Post-Invasion \n(Southern Sierra)","Enzootic \n(Northern Sierra)")) + 
  geom_smooth(data=sierra_no_larva, aes(x=log_mean, y=poulin_D),color="black")  + ylab("Poulin's D") + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) +
  labs(tag="A")+
  theme_classic() + theme_bw() + theme(legend.title = element_text(size=14),legend.text = element_text(size=12),axis.title = element_text(size=14),legend.position = "none")

cvl=ggplot(data=sierra_no_larva) + geom_point(size=2.5,aes(x=log_mean, y=cv_log, color=phase)) + scale_color_manual(name="Phase",values = c("goldenrod3",  "darkturquoise","maroon4"),labels = c("Invasion \n(Southern Sierra)","Post-Invasion \n(Southern Sierra)","Enzootic \n(Northern Sierra)")) + 
  geom_smooth(data=sierra_no_larva, aes(x=log_mean, y=cv_log),color="black") + ylab(expression(paste("CV [",log[10]," scale]",sep=""))) + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) +
  labs(tag="D")+
  theme_classic() + theme_bw() + theme(legend.title = element_text(size=14),legend.text = element_text(size=12),axis.title = element_text(size=14),legend.position = "none")

lcv=ggplot(data=sierra_no_larva) + geom_point(size=2.5,aes(x=log_mean, y=log_cv, color=phase)) + scale_color_manual(name="Phase",values = c("goldenrod3",  "darkturquoise","maroon4"),labels = c("Invasion \n(Southern Sierra)","Post-Invasion \n(Southern Sierra)","Enzootic \n(Northern Sierra)")) + 
  geom_smooth(data=sierra_no_larva, aes(x=log_mean, y=log_cv),color="black") + ylab(expression(paste(Log[10],"(CV [natural scale])",sep="")))  + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) +
  labs(tag="B")+
  theme_classic() + theme_bw() + theme(legend.title = element_text(size=14),legend.text = element_text(size=12),axis.title = element_text(size=14),legend.position = "none")

cv=ggplot(data=sierra_no_larva) + geom_point(size=2.5,aes(x=log_mean, y=cv, color=phase)) + scale_color_manual(name="Phase",values = c("goldenrod3",  "darkturquoise","maroon4"),labels = c("Invasion \n(Southern Sierra)","Post-Invasion \n(Southern Sierra)","Enzootic \n(Northern Sierra)")) + 
  geom_smooth(data=sierra_no_larva, aes(x=log_mean, y=cv),color="black")  + ylab("CV [natural scale]") + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) +
  labs(tag="C")+
  theme_classic() + theme_bw() + theme(legend.title = element_text(size=14),legend.text = element_text(size=12),axis.title = element_text(size=14),legend.position = "none")


# plot with legend 
full_plt = ggplot(data=sierra_no_larva) + geom_point(size=2.5,aes(x=log_mean, y=cv_log, color=phase)) + scale_color_manual(name="Phase",values = c("goldenrod3",  "darkturquoise","maroon4"),labels = c("Invasion \n(Southern Sierra)","Post-Invasion \n(Southern Sierra)","Enzootic \n(Northern Sierra)")) + 
  geom_smooth(data=sierra_no_larva, aes(x=log_mean, y=cv_log),color="black")  + ylim(c(0, 1)) + ylab(expression(paste("CV [",log[10]," scale]",sep=""))) + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) +
  theme_classic() + theme_bw() + theme(legend.title = element_text(size=14),legend.text = element_text(size=12),axis.title = element_text(size=14))

# function to extract legend from plot 
get_only_legend <- function(plot) { 
  plot_table <- ggplot_gtable(ggplot_build(plot)) 
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box") 
  legend <- plot_table$grobs[[legend_plot]] 
  return(legend) 
} 

# extract legend from using above function 
legend <- get_only_legend(full_plt)    

#Arranged figure
agg_plots = grid.arrange(
  grobs = list(pd, lcv, cv, cvl,legend),
  widths = c(2, 2, 1),
  layout_matrix = rbind(c(1, 2, NA),
                        c(3, 4, 5))
)

####################################
#### Step 3: Plot Model Results ####
#################################### 

eco_df=fread("data/model/eco_df.csv")
evo_df=fread("data/model/evo_df.csv")

both_df=rbind(eco_df[c(7,365),],evo_df[365,])
both_df$phase=unique(sierra_dat$phase)[c(1,3,2)]

pd_model=ggplot(data=eco_df)+geom_path(data=eco_df, aes(x=norm_means_everyone*log10(5)/log(5), y=PD_everyone_Mc),color="black")  + ylab("Poulin's D") + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) +
  theme_classic() + theme_bw() + theme(legend.title = element_text(size=14),legend.text = element_text(size=12),axis.title = element_text(size=14),legend.position = "none")+
  geom_path(data=evo_df, aes(x=norm_means_everyone*log10(5)/log(5), y=PD_everyone_Mc),color="black",linetype="dashed")+
  geom_point(size=2.5,data=both_df,aes(x=norm_means_everyone*log10(5)/log(5), y=PD_everyone_Mc, color=phase))+
  geom_path(data=both_df,aes(x=norm_means_everyone*log10(5)/log(5), y=PD_everyone_Mc),col="gray",linetype="dashed")+
  labs(tag="E")+
  scale_color_manual(name="Phase",values = c("maroon4","goldenrod3","darkturquoise"),labels = c("Invasion \n(Southern Sierra)","Post-Invasion \n(Southern Sierra)","Enzootic \n(Northern Sierra)"))

cvl_model=ggplot(data=eco_df)+geom_path(data=eco_df, aes(x=norm_means_everyone*log10(5)/log(5), y=V2/norm_means_everyone),color="black")  + ylab(expression(paste("CV [", log[10], " scale]",sep="")))  + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) +
  theme_classic() + theme_bw() + theme(legend.title = element_text(size=14),legend.text = element_text(size=12),axis.title = element_text(size=14),legend.position = "none")+
  geom_path(data=evo_df, aes(x=norm_means_everyone*log10(5)/log(5), y=V2/norm_means_everyone),color="black",linetype="dashed")+
  geom_point(size=2.5,data=both_df,aes(x=norm_means_everyone*log10(5)/log(5), y=V2/norm_means_everyone, color=phase))+
  geom_path(data=both_df,aes(x=norm_means_everyone*log10(5)/log(5), y=V2/norm_means_everyone),col="gray",linetype="dashed")+
  labs(tag="H")+
  scale_color_manual(name="Phase",values = c("maroon4","goldenrod3","darkturquoise"),labels = c("Invasion \n(Southern Sierra)","Post-Invasion \n(Southern Sierra)","Enzootic \n(Northern Sierra)"))

lcv_model=ggplot(data=eco_df)+geom_path(data=eco_df, aes(x=norm_means_everyone*log10(5)/log(5), y=log_cv_everyone),color="black")  + ylab(expression(paste(Log[10],"(CV [natural scale])",sep=""))) + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) +
  theme_classic() + theme_bw() + theme(legend.title = element_text(size=14),legend.text = element_text(size=12),axis.title = element_text(size=14),legend.position = "none")+
  geom_path(data=evo_df, aes(x=norm_means_everyone*log10(5)/log(5), y=log_cv_everyone),color="black",linetype="dashed")+
  geom_point(size=2.5,data=both_df,aes(x=norm_means_everyone*log10(5)/log(5), y=log_cv_everyone, color=phase))+
  geom_path(data=both_df,aes(x=norm_means_everyone*log10(5)/log(5), y=log_cv_everyone),col="gray",linetype="dashed")+
  labs(tag="F")+
  scale_color_manual(name="Phase",values = c("maroon4","goldenrod3","darkturquoise"),labels = c("Invasion \n(Southern Sierra)","Post-Invasion \n(Southern Sierra)","Enzootic \n(Northern Sierra)"))

cv_model=ggplot(data=eco_df)+geom_path(data=eco_df, aes(x=norm_means_everyone*log10(5)/log(5), y=V4/nat_means_everyone),color="black")  + ylab("CV [natural scale]") + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) +
  theme_classic() + theme_bw() + theme(legend.title = element_text(size=14),legend.text = element_text(size=12),axis.title = element_text(size=14),legend.position = "none")+
  geom_path(data=evo_df, aes(x=norm_means_everyone*log10(5)/log(5), y=V4/nat_means_everyone),color="black",linetype="dashed")+
  geom_point(size=2.5,data=both_df,aes(x=norm_means_everyone*log10(5)/log(5), y=V4/nat_means_everyone, color=phase))+
  geom_path(data=both_df,aes(x=norm_means_everyone*log10(5)/log(5), y=V4/nat_means_everyone),col="gray",linetype="dashed")+
  labs(tag="G")+
  scale_color_manual(name="Phase",values = c("maroon4","goldenrod3","darkturquoise"),labels = c("Enzootic \n(Northern Sierra)","Invasion \n(Southern Sierra)","Post-Invasion \n(Southern Sierra)"))


#################################################
#### Step 4: Combine empirical & model plots ####
#################################################   

png("results/plots/agg_emp_mod.png",width=12,height=6,units="in",res=300)
agg_plots_all = grid.arrange(
  grobs = list(pd, lcv, pd_model, lcv_model, cv, cvl,legend, cv_model, cvl_model),
  top=textGrob(c("Empirical results","Model results"),x=c(0.22,0.82),gp=gpar(fontface="bold"),hjust=0.5),
  layout_matrix = rbind(c(1, 2, NA, 3, 4),
                        c(5, 6, 7, 8, 9))
)
dev.off()


