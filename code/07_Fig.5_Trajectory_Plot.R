### Reproduce Figure 5 (Trajectory plot & mean-aggregation plot for 7 well-studied sites)
### Step 1: Load Sierra summary file and remove records with fewer than 2 Bd-positive and larval groups (for consistent comparisons)
### Step 2: Select specific sites, relabel, and plot time vs. host abundance
### Step 3: Plot these groups by log mean fungal load against Poulin's D aggregation metric
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

sierra_dat = fread("data/formatted/summary_sierra_with_epi_phase.csv")
sierra_dat[, "log_mean":=mean_log10_without_zeros]
sierra_no_larva = sierra_dat[life_stage != "larva" & num_infected >= 2 & poulin_D > 0]

abund_and_prev = fread("data/formatted/abund_and_prev.csv")

######################################
####    Step 2: Trajectory Plot   ####
######################################

trajectory_sites = c("10100", "10101", "10102", "11858", "12590", "10090", "12618")
y=sierra_dat[site_id %in% trajectory_sites,]

sierra_dat$Year <- plyr::mapvalues(sierra_dat$year, 
                                   from=c(1:16), 
                                   to=c(2004:2019))
sierra_dat$Site <- plyr::mapvalues(sierra_dat$site, 
                                   from=c(10100, 10101, 10102, 11858, 12590, 10090, 12618), 
                                   to=c("Site A","Site B","Site C","Site D","Site E","Site F","Site G"))
abund_and_prev$Site <- plyr::mapvalues(abund_and_prev$site_id, 
                                       from=c(10100, 10101, 10102, 11858, 12590, 10090, 12618), 
                                       to=c("Site A","Site B","Site C","Site D","Site E","Site F","Site G"))


abunds = merge(sierra_dat[num_infected >= 2 & site_id %in% trajectory_sites & life_stage == "adult"], 
               abund_and_prev[capture_life_stage == "adult", .(site_id, year, adult_abundance)], 
               by.x=c("site_id", "Year"), by.y=c("site_id", 'year'))[, .(site_id, Year, adult_abundance, phase, Site)]

abund_trun = abund_and_prev[site_id %in% trajectory_sites & capture_life_stage == "adult" & adult_abundance != -1]

#Abundance graph
p1 = ggplot(abunds) +
  geom_line(data=abund_trun, aes(x=year, y=log10(adult_abundance)),linewidth = 0.5) + geom_point(aes(x=Year, y=log10(adult_abundance), color=phase), size=3.5) +
  geom_point(data=abund_trun, aes(x=year, y=log10(adult_abundance)),size=2) +
  scale_color_manual(values=c('goldenrod3','darkturquoise'))  + theme_bw() + ylab(expression(paste(Log[10],"(Adult Abundance)", sep="")))  +
  facet_wrap(~Site)+
  xlab("Year") + 
  guides(color="none") + labs(tag="A") + theme(title=element_text(size=14),axis.title = element_text(size=14),axis.text.x = element_text(angle = 45,margin = margin(t = 5)),legend.title = element_text(size=14),legend.text = element_text(size=12),plot.tag = element_text(size=14),axis.text = element_text(size=10))



######################################
####Step 3: Mean-Aggregation Plot ####
######################################

#Load-Aggregation plot
p2 = ggplot(sierra_dat[num_infected >= 2 & site_id %in% trajectory_sites & life_stage == "adult"][order(Site, year)]) +
  geom_path(aes(x=mean_log10_without_zeros, y=poulin_D, group=Site),linewidth = 0.5) + 
  geom_point(aes(x=mean_log10_without_zeros, y=poulin_D, color=phase), size=3.5) +
  facet_wrap(~Site)  +
  scale_color_manual(name="Phase",values=c('goldenrod3','darkturquoise'),labels=c("Invasion \n(Southern Sierra)","Post-Invasion \n(Southern Sierra)","Enzootic \n(Northern Sierra)"))  + theme_bw() + ylab("Poulin's D") + xlab(expression(paste("Mean ",Log[10],"(Bd Load)", sep=""))) + labs(tag="B")  + theme(axis.title = element_text(size=14),legend.title = element_text(size=14),legend.text = element_text(size=12),plot.tag = element_text(size=14),axis.text = element_text(size=10),legend.position = c(0.7, 0.12))


######################################
####     Step 4: Combine Plots    ####
######################################

#Arranged figure
full_plot = grid.arrange(
  grobs = list(p1, p2),
  widths = c(1, 1)
)

#save figure
ggsave("results/plots/trajectory_plot.png",full_plot, width=180,height=105,units="mm",dpi=600)


