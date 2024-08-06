### Reproduce Figure S2 (Pennsylvania green frog mean-aggregation plot)
### Step 1: Load, format, and select relevant data
### Step 2: Plot empirical data by log mean fungal load against Poulin's D
### Step 3: Remove population seemingly in the "invasion space" and re-plot
### Step 4: Combine plots


library(ggplot2)
library(gridExtra)
library(cowplot)
library(grid)
library(dplyr)
library(data.table)
library(gridtext)


######################################
#### Step 1: Load and format data ####
######################################

#Load and combine summary datasets
three_dat = fread("data/formatted/summary_aggregation_full.csv")[region != "california" & dataset != "panama"]
full_dat = rbind(three_dat,sierra_dat,fill=TRUE)
dat = full_dat[full_dat$num_infected>=3,]

dat$var_natural_with_zeros = as.numeric(dat$var_natural_with_zeros)
dat$mean_natural_with_zeros = as.numeric(dat$mean_natural_with_zeros)

dat_penn = dat[dat$region=="pennsylvania" & dat$species=="rana_clamitans",]

##############################################
#### Step 2: Create mean-aggregation plot ####
##############################################

#Plot PA green frogs in mean load-agg space
penn_green_plot = ggplot(data=dat_penn) + geom_point(size=2.5,aes(x=mean_log10_without_zeros, y=poulin_D,color=region),show.legend = TRUE) + scale_color_manual(values="firebrick",labels="Pennsylvania\n green frogs",name="") +
  geom_smooth(data=dat_penn, aes(x=mean_log10_without_zeros, y=poulin_D),color="black") + ylab("Poulin's D") + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) + ylim(-0.1,1)+ xlim(1.4,4.5) + theme_bw() + theme(legend.title = element_text(size=14),legend.text = element_text(size=12),axis.title = element_text(size=14),legend.position = "none")+
  labs(tag="A")


#########################################################
#### Step 3: Plot data without "invasion" population ####
#########################################################

#Which population appears to be in the invasion phase (i.e. Poulin's D under 0.25)?
x=dat_penn[dat_penn$poulin_D < 0.25,]

#Plot without invasion population (see if pattern holds)
penn_green_trun_plot = ggplot(data=dat_penn[dat_penn$poulin_D>0.25,]) + geom_point(size=2.5,aes(x=mean_log10_without_zeros, y=poulin_D,color=region),show.legend = TRUE) + scale_color_manual(values="firebrick",labels="Pennsylvania\n green frogs",name="") +
  geom_smooth(aes(x=mean_log10_without_zeros, y=poulin_D),color="black") + ylab("") + ylim(-0.1,1) + xlim(1.4,4.5) + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) + theme_bw() + theme(legend.title = element_text(size=14),legend.text = element_text(size=12),axis.title = element_text(size=14),legend.position = "none")+
  labs(tag="B")

###############################
#### Step 4: Combine plots ####
###############################

green_w_leg = ggplot(data=dat_penn[dat_penn$poulin_D>0.25,]) + geom_point(size=2.5,aes(x=mean_log10_without_zeros, y=poulin_D,color=region),show.legend = TRUE) + scale_color_manual(values="firebrick",labels="Pennsylvania\n green frogs",name="") +
  geom_smooth(aes(x=mean_log10_without_zeros, y=poulin_D),color="black") + ylab("") + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) + theme_bw() + theme(legend.title = element_text(size=14),legend.text = element_text(size=12),axis.title = element_text(size=14),legend.position = "right")

# extract legend from using above function 
legend <- get_only_legend(green_w_leg)    

full_green_plot = grid.arrange(
  grobs = list(penn_green_plot,penn_green_trun_plot,legend),
  widths = c(2,2,1),
  layout_matrix = rbind(c(1, 2, 3))
)

ggsave("C:/Users/sarss/OneDrive/Documents/Univ of TN/BD_Database/fungal_intensity/results/plots/penn_green_full_plot.png",full_green_plot, width=180,height=105,units="mm",dpi=600)
