### Extract data for Table S1


library(dplyr)
library(data.table)

bd = fread("data/formatted/full_aggregation_data.csv")[location != "panama" & region !="new_mexico"]
bd$bd_load = bd$target_quant_per_swab 
bd_pos=bd[bd$target_quant_per_swab>0,]

full_dat = fread("data/formatted/summary_aggregation_full.csv")[region != "california" & dataset != "panama"]
sierra_dat = fread("data/formatted/summary_sierra_with_epi_phase.csv")
all_data = rbind(full_dat,sierra_dat,fill=TRUE)

#Statistics on database
nrow(bd) #Number of records (56,912) 
nrow(bd_pos) #Number of positive records (26,451) 

unique(bd$region)
unique(all_data$region)

bd%>%
  group_by(location) %>%
  count(region)
#Brazil: 4365 records 
#Sierra: 29,600 records
#East Bay: 10,490 records
#SERDP: 12,457 (2663+5937+2817+1040)


bd%>%
  count(region)

nrow(bd[bd$region=="louisiana" | bd$region=="vermont" | bd$region=="pennsylvania" | bd$region=="tennessee",])

bd[bd$region=="louisiana" | bd$region=="vermont" | bd$region=="pennsylvania" | bd$region=="tennessee",] %>%
  count(species_capture)

bd[bd$location=="brazil",] %>%
  count(species_capture)

unique(bd$species_capture)

#Find number of unique species of all samples
genus_level = unique(bd$species_capture[like(bd$species_capture, "spp") | endsWith(bd$species_capture, "_sp") | endsWith(bd$species_capture, "_sp.") | like(bd$species_capture, "_sp._") | like(bd$species_capture, "_sp_")])
unknowns = unique(bd$species_capture[like(bd$species_capture,"unknown") | like(bd$species_capture,"see_notes") | like(bd$species_capture,"possibly")])
repeats = unique(bd$species_capture[like(bd$species_capture,"RACA") | like(bd$species_capture,"lithobates_sylvaticus") | like(bd$species_capture,"bufo_americanus")]) #[lithobates_sylvaticus/rana_sylvatica, bufo_americanus/anaxyrus_americanus, RACA/rana_catesbeiana] 

unique_spec = bd$species_capture[!(bd$species_capture %in% genus_level) & !(bd$species_capture %in% unknowns) & !(bd$species_capture %in% repeats)]
length(unique(unique_spec)) #number of unique species (93)

#Find number of unique species of all bd-positive
unique_spec_pos = bd_pos$species_capture[!(bd_pos$species_capture %in% genus_level) & !(bd_pos$species_capture %in% unknowns) & !(bd_pos$species_capture %in% repeats)]
length(unique(unique_spec_pos)) #number of unique species (78)

### Statistics on datset populations
all = all_data %>%
  count(all_data$dataset)
#brazil: 109 pops 
#eastbay: 714 pops
#serdp: 391 pops
#sierra: 647 pops

#Datasets limited by number of infected
three_infected = all_data[all_data$num_infected>2,]
more_than_two = three_infected %>%
  count(three_infected$dataset)
#brazil: 49 
#eastbay: 267
#serdp: 231
#sierra: 414

sum(more_than_two$n) #961

unique(three_infected$species[!(three_infected$species %in% genus_level) & !(three_infected$species %in% unknowns) & !(three_infected$species %in% repeats)])
#number of unique species (35) with 3+ infected

#Datasets limited by number of infected
two_infected = all_data[all_data$num_infected>1,]
more_than_one = two_infected %>%
  count(two_infected$dataset)
#brazil: 59 
#eastbay: 326
#serdp: 257
#sierra: 444

sum(more_than_one$n)

mean(two_infected$num_infected[two_infected$dataset=="brazil"]) #5.36
sd(two_infected$num_infected[two_infected$dataset=="brazil"]) #3.44
mean(two_infected$num_infected[two_infected$dataset=="east_bay"]) #6.63
sd(two_infected$num_infected[two_infected$dataset=="east_bay"]) #4.56
mean(two_infected$num_infected[two_infected$dataset=="serdp"]) #10.28
sd(two_infected$num_infected[two_infected$dataset=="serdp"]) #8.40
mean(two_infected$num_infected[two_infected$dataset=="sierra_nevada"]) #35.23
sd(two_infected$num_infected[two_infected$dataset=="sierra_nevada"]) #64.59


unique(two_infected$species[!(two_infected$species %in% genus_level) & !(two_infected$species %in% unknowns) & !(two_infected$species %in% repeats)])
#number of unique species (41) with 2+ infected


ten_infected = all_data[all_data$num_infected>9,]
more_than_nine = ten_infected %>%
  count(ten_infected$dataset)
#brazil: 5
#eastbay: 81
#serdp: 110
#sierra: 329

sum(more_than_nine$n)

unique(ten_infected$species[!(ten_infected$species %in% genus_level) & !(ten_infected$species %in% unknowns) & !(ten_infected$species %in% repeats)])
#number of unique species (22)

#Compare number of records across these subsets
all_tab = cbind(all,more_than_one,more_than_nine)
all_df = as.data.frame(all_tab[,c(1,2,4,6)])
colnames(all_df) = c("dataset","all","at_least_2","at_least_10")


ls_list = two_infected %>%
  group_by(two_infected$dataset) %>%
  count(life_stage)

ls_list$n_tot = ifelse(ls_list$`two_infected$dataset`=="brazil",59,ifelse(ls_list$`two_infected$dataset`=="east_bay",326,ifelse(ls_list$`two_infected$dataset`=="serdp",257,444)))

ls_list$perc = (ls_list$n / ls_list$n_tot)*100

two_infected %>%
  group_by(two_infected$dataset) %>%
  count(phase)

szn_list = two_infected %>%
  group_by(two_infected$dataset) %>%
  count(season)

szn_list$n_tot = ifelse(szn_list$`two_infected$dataset`=="brazil",59,ifelse(szn_list$`two_infected$dataset`=="east_bay",326,ifelse(szn_list$`two_infected$dataset`=="serdp",257,444)))

szn_list$perc = (szn_list$n / szn_list$n_tot)*100

spec_list = two_infected %>%
  group_by(two_infected$dataset) %>%
  count(species)


spec_list$n_tot = ifelse(spec_list$`two_infected$dataset`=="brazil",59,ifelse(spec_list$`two_infected$dataset`=="east_bay",326,ifelse(spec_list$`two_infected$dataset`=="serdp",257,444)))

spec_list$perc = (spec_list$n / spec_list$n_tot)*100


#average and sd of fungal load for each dataset
ds_sum = as.data.frame(dat %>%
                         group_by(dataset) %>%
                         summarise(n_pop=n(),ave=mean(N),sd=sd(N)))

#number of populations for each life stage and dataset
ls_sum = as.data.frame(dat %>%
                         group_by(dataset) %>%
                         count(life_stage))

merge_sum = merge(ds_sum,ls_sum,by="dataset") 
merge_sum$perc = (merge_sum$n/merge_sum$n_pop)*100
merge_sum

#number of populations for each phase
phz_sum = as.data.frame(dat %>%
                          group_by(dataset) %>%
                          count(phase))
merge_phz = merge(ds_sum,phz_sum,by="dataset")
merge_phz$perc = (merge_phz$n/merge_phz$n_pop)*100
merge_phz


#number of populations for each season
szn_sum = as.data.frame(dat %>%
                          group_by(dataset) %>%
                          count(season))
merge_szn = merge(ds_sum,szn_sum,by="dataset")
merge_szn$perc = (merge_szn$n/merge_szn$n_pop)*100
merge_szn


#number of populations for each species
spp_sum = as.data.frame(dat %>%
                          group_by(dataset) %>%
                          count(species))
merge_spp = merge(ds_sum,spp_sum,by="dataset")
merge_spp$perc = (merge_spp$n/merge_spp$n_pop)*100
merge_spp

morethan2 = all_data[all_data$num_infected>2,]
morethan2 %>%
  group_by(dataset) %>%
  summarize(sum = sum(N))

all_data %>%
  group_by(dataset) %>%
  summarize(sum = sum(N))
