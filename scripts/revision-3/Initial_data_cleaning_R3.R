rm(list=ls())

#install.packages("TMB", type = "source")
#install.packages("sjPlot")
#install.packages("plyr")
#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("lmerTest")
#install.packages("ggforce")

library(plyr)
library(dplyr)
library(ggplot2)
library(lmerTest)
library(sjPlot)
library(ggforce)
library(glmmTMB)

#set home directory
#homewd= "/Users/carabrook/Developer/git-transit-time"
homewd= "/Users/katherinemcferrin/Developer/git-transit-time"
#homewd <- "/Users/gavindehnert/Desktop/GitHub_repos/git-transit-time"
setwd(homewd)

#load the GIT transit data:
#dat <- read.csv(file = paste0(homewd, "/data/final-GIT-transit-database-april-2021.csv"), header = T, stringsAsFactors = F )
dat <- read.csv(file = paste0(homewd, "/data/revision-3/RuhsMcFerrin_database_R3_v3.csv"), header = T, stringsAsFactors = F )
#View(dat)

# R3: this is still the same even after Katherines additions of non-flying MRT
length(unique(dat$Retention.Citation)) #127 unique papers -> 791 papers
length(unique(dat$Genus.species)) #144 unique species -> 434 species
sort(unique(dat$Genus.species))

#We also added the collection of MRT, but never plotted it
#I'm pulling that out here to include in the actual paper - starts on line 561

#dat = only the columns we will use for modeling and plotting and renaming columns
#original code: dat <- dplyr::select(dat, can_fly, Class, Order, Family, Genus.species, Common.Name, Typical.Diet, Mass..mean.median.from.study., Feeding.Trial.Food.Cat, N_measured, N_trials, Minimum..min., Median..min., Mean..min.,Maximum..min.)
dat <- dplyr::select(dat, can_fly, Class, Order, Family, Genus.species, Common.Name, Mass.g..mean.median.from.study., Feeding.Trial.Food.Cat, N_individuals, N_trials, Minimum..min., Median..min., Mean..min.,Maximum..min., mrt..min.,mrt_sd, mrt_se)
head(dat)
#rename columns
names(dat) <- c("fly", "class", "order", "family", "genus.species", "common.name", "mass", "diet","N_individuals", "N_trials", "min", "median", "mean", "max", "MRT_min", "MRT_sd", "MRT_se")
#View(dat)

#R3: KEEPING THOSE IN WITH UNKOWN DIETS
#dat <- subset(dat, diet!="unknown") # these are the citations that only had label-solution as a diet and no other info
#length(unique(dat$genus.species)) #142 unique species

#choose one mean/median to report
dat.plot <- dat
#View(dat.plot)
dat.plot$transit <- dat.plot$mean #copying the mean git transit values to new column transit
dat.plot$transit[is.na(dat.plot$transit)]<- dat.plot$median[is.na(dat.plot$transit)] #if transit doesn't have a mean, fill it with median

dat.plot <- dplyr::select(dat.plot, -(min), -(median), -(mean), -(max)) #simplifying by taking out min, median, mean and max
head(dat.plot)

#summarize by species
dat.plot$N_individuals[dat.plot$N_individuals=="not reported"] <- 1 #assigning 1 individual if total # not reported
dat.plot$N_individuals[is.na(dat.plot$N_individuals)] <- 1 #assigning 1 individual if #individuals is NA
dat.plot$N_individuals[dat.plot$N_individuals=="unknown"] <- 1
dat.plot$N_individuals[dat.plot$N_individuals=="unk"] <- 1
dat.plot$N_individuals[dat.plot$N_individuals==""] <- 1 #Cara added this
unique(dat.plot$N_individuals)
dat.plot$N_individuals = as.numeric(dat.plot$N_individuals)
dat.plot$total_transit = dat.plot$transit*dat.plot$N_individuals #multiplying transit time x #individuals

#and also add for MRT
sort(unique(dat.plot$MRT_min))
dat.plot$MRT_min = as.numeric(dat.plot$MRT_min)
dat.plot$total_MRT = dat.plot$MRT_min*dat.plot$N_individuals

# create new dataset that only has MRT and remove NAs
dat.plot2 <- dat.plot
dat.plot2 <- dplyr::select(dat.plot2, -(total_transit), -(transit)) 
dat.plot2 <- dat.plot2 %>% filter(!is.na(total_MRT))

# remove entries that are NA for git transit
dat.plot <- dplyr::select(dat.plot, -(total_MRT), -(MRT_min), -(MRT_sd), -(MRT_se)) 
dat.plot <- dat.plot %>% filter(!is.na(total_transit))

#what are the numbers on non-flying birds?
nonflyingbirds_git <- subset(dat.plot, dat.plot$fly=="N" & dat.plot$class=="Aves")
length(unique(nonflyingbirds_git$genus.species)) #n=5
nonflyingbirds_git <- nonflyingbirds_git %>% filter(!is.na(mass))
length(unique(nonflyingbirds_git$genus.species)) #n=5

nonflyingbirds_mrt <- subset(dat.plot2, dat.plot2$fly=="N" & dat.plot2$class=="Aves")
length(unique(nonflyingbirds_mrt$genus.species)) #n=7
nonflyingbirds_mrt <- nonflyingbirds_mrt %>% filter(!is.na(mass))
length(unique(nonflyingbirds_mrt$genus.species)) #n=7

################## CLEANING FOR GIT TRANSIT ONLY ##########################
#now get one entry for each species for GIT_transit
dat.split <- dlply(dat.plot, .(fly,class,order,family,genus.species, common.name, diet))
#making a function

str(dat)

summarise.dat <- function(dat){
  
  dat$mass = as.numeric(dat$mass)
  dat$N_individuals = as.numeric(dat$N_individuals)
  dat$total_transit = as.numeric(dat$total_transit)
  dat2 <- ddply(dat, .(fly, class, order, family, genus.species, common.name), summarise, 
                sum_mass = sum(mass, na.rm = TRUE), 
                trial.diet = unique(diet),  
                N_tot = sum(N_individuals, na.rm = TRUE), 
                total_transit = if (all(is.na(total_transit))) NA else sum(total_transit, na.rm = TRUE))
  
  return(dat2)
}

dat.out <- lapply(dat.split, summarise.dat)
dat.sum.tot <- data.table::rbindlist(dat.out)
head(dat.sum.tot)
dat.sum.tot$transit <- dat.sum.tot$total_transit/dat.sum.tot$N_tot
#dat.sum.tot$MRT <- dat.sum.tot$total_MRT/dat.sum.tot$N_tot
dat.sum.tot$avg_mass <- dat.sum.tot$sum_mass/dat.sum.tot$N_tot
nrow(subset(dat.sum.tot, !is.na(transit))) #313. Cara: I get 326



#now categorize -- vertebrates are signified by class but mammals by order
dat.sum.tot$re_class <- dat.sum.tot$class
dat.sum.tot$re_class[dat.sum.tot$re_class=="Mammalia"] <- dat.sum.tot$order[dat.sum.tot$re_class=="Mammalia"]
dat.sum.tot$re_class[dat.sum.tot$re_class=="Aves" & dat.sum.tot$fly=="Y"] <- "Flying Birds"
dat.sum.tot$re_class[dat.sum.tot$re_class=="Aves" & dat.sum.tot$fly=="N"] <- "Non-Flying Birds"

unique(dat.sum.tot$re_class)
#find out how many entries per each "cat"
# dat.simp.big <- ddply(dat.sum.tot, .(re_class, typical.diet), summarize, N=length(re_class))
# dat.simp.big

paper.dat <- ddply(dat.sum.tot, .(re_class), summarise, N_species = length(unique(genus.species)))
paper.dat

# re_class N_species
# 1          Amphibia         1
# 2      Artiodactyla        17
# 3         Carnivora        26
# 4           Cetacea         1
# 5        Chiroptera        37
# 6    Chondrichthyes         2
# 7    Dasyuromorphia         1
# 8        Dermoptera         1
# 9     Diprotodontia         3
# 10     Flying Birds        62
# 11       Lagomorpha         1
# 12    Macroscelidea         1
# 13 Non-Flying Birds         5
# 14  Peramelemorphia         2
# 15   Perissodactyla         7
# 16         Primates        60
# 17      Proboscidea         2
# 18         Reptilia        19
# 19         Rodentia        15
# 20          Sirenia         2
# 21     Soricomorpha         3

write.csv(dat.sum.tot, "data/revision-3/dat_sum_tot_GIT_R3.csv")

#and group
#remove any non-mammalian classes and any mammalian orders with < =4 entries
# also removing marsupials because they are weird (and all < 4 spp)
#Amphibia, Cetacea, Chondrichtythes, Dermoptera, Pilosa, Sirenia

dat.sum.tot = subset(dat.sum.tot, re_class !="Amphibia" & 
                       re_class != "Cetacea" & 
                       re_class != "Chondrichthyes" & 
                       re_class != "Cetartiodactyla" &
                       re_class != "Dermoptera" & 
                       re_class != "Dasyuromorphia" &
                       re_class != "Diprotodontia" & #marsupials
                       re_class != "Lagomorpha" &
                       re_class != "Peramelemorphia" &
                       re_class != "Macroscelidea" &
                       re_class != "Pilosa" &
                       re_class != "Proboscidea" &
                       #re_class != "Reptilia" &
                       re_class != "Soricomorpha" & 
                       re_class != "Sirenia") 


unique(dat.sum.tot$re_class)
# > unique(dat.sum.tot$re_class)
# [1] "Flying Birds"     "Carnivora"        "Artiodactyla"     "Primates"         "Chiroptera"      
# [6] "Reptilia"         "Non-Flying Birds" "Rodentia"         "Perissodactyla"  

#and write over some of the others
dat.sum.tot$re_class[dat.sum.tot$re_class=="Reptilia"] <- "Reptiles"
dat.sum.tot$re_class[dat.sum.tot$re_class=="Rodentia"] <- "Rodents"
dat.sum.tot$re_class[dat.sum.tot$re_class=="Chiroptera"] <- "Bats"
dat.sum.tot$re_class[dat.sum.tot$re_class=="Artiodactyla"] <- "Even-toed Ungulates"
dat.sum.tot$re_class[dat.sum.tot$re_class=="Perissodactyla"] <- "Odd-toed Ungulates"
dat.sum.tot$re_class[dat.sum.tot$re_class=="Carnivora"] <- "Carnivores"
#dat.sum.tot$re_class[dat.sum.tot$re_class=="Non-Flying Birds"] <- "Non-Flying Birds"
#dat.sum.tot$re_class[dat.sum.tot$re_class=="Perissodactyla" | dat.sum.tot$re_class=="Artiodactyla"] <- "Ungulates"

#dat.sum.tot$re_class[dat.sum.tot$re_class=="Diprotodontia" | dat.sum.tot$re_class=="Peramelemorphia"] <- "Marsupials"

#dat.sum.tot$re_class[dat.sum.tot$re_class=="Lagomorpha" ] <- "Lagomorphs"

unique(dat.sum.tot$re_class)
# [1] "Flying Birds"        "Carnivores"          "Even-toed Ungulates" "Primates"            "Bats"               
# [6] "Reptiles"            "Non-Flying Birds"    "Rodents"             "Odd-toed Ungulates" 

dat.simp <- ddply(dat.sum.tot, .(re_class), summarize, N=length(unique(genus.species)))
dat.simp
# re_class  N
# 1                Bats 37
# 2          Carnivores 26
# 3 Even-toed Ungulates 17
# 4        Flying Birds 62
# 5    Non-Flying Birds  5
# 6  Odd-toed Ungulates  7
# 7            Primates 60
# 8            Reptiles 19
# 9             Rodents 15



#convert transit time to hours
dat.sum.tot$transit_hrs = dat.sum.tot$transit/60
#dat.sum.tot$MRT_hrs = dat.sum.tot$MRT/60
#convert mass to kg
dat.sum.tot$mass_kg = as.numeric(dat.sum.tot$avg_mass)/1000

#write.csv(dat.sum.tot, "data/R_cleaned_data_R2.csv")

#check how transit time varies by mass
#generally, we see that transit is longer at higher mass
#all taxa agree with this trend, though reptiles and bats are a little longer in time
#than we might otherwise predict -- this is likely due to an effect of temperature that is not
#reported here.

write.csv(dat.sum.tot, "data/revision-3/dat_sum_tot_clean_GIT_R3.csv")
# dat.sum.tot_clean <- read.csv(file = paste0(homewd, "data/dat.sum.tot.csv"), header = T, stringsAsFactors = F )
#dat.sum.tot_clean <- dat.sum.tot










################## CLEANING FOR MRT ONLY ##########################
#now get one entry for each species for GIT_transit
dat.split <- dlply(dat.plot2, .(fly,class,order,family,genus.species, common.name, diet))
#making a function


summarise.dat2 <- function(dat){
  dat$mass = as.numeric(dat$mass)
  dat$N_individuals = as.numeric(dat$N_individuals)
  dat$total_MRT = as.numeric(dat$total_MRT)
  dat2 <- ddply(dat, .(fly, class, order, family, genus.species, common.name), summarise,
                sum_mass = sum(mass),
                trial.diet = unique(diet),
                N_tot = sum(N_individuals),
                total_MRT = if (all(is.na(total_MRT))) NA else sum(total_MRT, na.rm = TRUE))
  return(dat2)
}


dat.out <- lapply(dat.split, summarise.dat2)
dat.sum.tot2 <- data.table::rbindlist(dat.out)
head(dat.sum.tot2)
#dat.sum.tot$transit <- dat.sum.tot$total_transit/dat.sum.tot$N_tot
dat.sum.tot2$MRT <- dat.sum.tot2$total_MRT/dat.sum.tot2$N_tot
dat.sum.tot2$avg_mass <- dat.sum.tot2$sum_mass/dat.sum.tot2$N_tot
#subset(dat.sum.tot, !is.na(transit)) 
subset(dat.sum.tot2, !is.na(total_MRT)) 


#now categorize -- vertebrates are signified by class but mammals by order
dat.sum.tot2$re_class <- dat.sum.tot2$class
dat.sum.tot2$re_class[dat.sum.tot2$re_class=="Mammalia"] <- dat.sum.tot2$order[dat.sum.tot2$re_class=="Mammalia"]
dat.sum.tot2$re_class[dat.sum.tot2$re_class=="Aves" & dat.sum.tot2$fly=="Y"] <- "Flying Birds"
dat.sum.tot2$re_class[dat.sum.tot2$re_class=="Aves" & dat.sum.tot2$fly=="N"] <- "Non-Flying Birds"

unique(dat.sum.tot2$re_class)
# [1] "Flying Birds"     "Artiodactyla"     "Primates"         "Non-Flying Birds" "Diprotodontia"   
# [6] "Carnivora"        "Rodentia"         "Eulipotyphla"     "Cetartiodactyla"  "Pilosa"          
# [11] "Perissodactyla"   "Dermoptera"       "Reptilia"         "Sirenia"          "Proboscidea"     
# [16] "Chiroptera"       "Cetacea"          "Peramelemorphia"  "Lagomorpha"       "Amphibia"        
# [21] "Dasyuromorphia"   "Soricomorpha"   

paper.dat <- ddply(dat.sum.tot2, .(re_class), summarise, N_species = length(unique(genus.species)))
paper.dat

# re_class N_species
# 1          Amphibia         1
# 2      Artiodactyla        55
# 3         Carnivora        18
# 4           Cetacea         1
# 5   Cetartiodactyla         2
# 6        Chiroptera        13
# 7    Dasyuromorphia         1
# 8        Dermoptera         1
# 9     Diprotodontia        14
# 10     Eulipotyphla         1
# 11     Flying Birds        69
# 12       Lagomorpha         1
# 13 Non-Flying Birds         7
# 14  Peramelemorphia         2
# 15   Perissodactyla        10
# 16           Pilosa         3
# 17         Primates        40
# 18      Proboscidea         2
# 19         Reptilia         3
# 20         Rodentia        27
# 21          Sirenia         2
# 22     Soricomorpha         1

write.csv(dat.sum.tot2, "data/revision-3/dat_sum_tot_MRT_R2.csv")

#and group
#remove any non-mammalian classes and any mammalian orders with <= 4 entries
#Amphibia, Cetacea, Chondrichtythes, Dermoptera, Pilosa, Sirenia

dat.sum.tot2 = subset(dat.sum.tot2, re_class !="Amphibia" & 
                       re_class != "Cetacea" &
                        re_class != "Cetartiodactyla" &
                       re_class != "Dermoptera" & 
                       re_class != "Diprotodontia" & #marsupial
                        re_class != "Dasyuromorphia" &
                       re_class != "Eulipotyphla" &
                       re_class != "Lagomorpha" &
                       re_class != "Peramelemorphia" & #marsupial
                       #re_class != "Perissodactyla" &
                       re_class != "Pilosa" & 
                        re_class != "Proboscidea" & 
                       re_class != "Reptilia" &
                        re_class != "Soricomorpha" &
                       re_class != "Sirenia") 


unique(dat.sum.tot2$re_class)
# [1] "Flying Birds"     "Artiodactyla"     "Primates"         "Non-Flying Birds" "Carnivora"       
# [6] "Rodentia"         "Perissodactyla"   "Chiroptera" 

#and write over some of the others
#dat.sum.tot2$re_class[dat.sum.tot2$re_class=="Reptilia"] <- "Reptiles"
dat.sum.tot2$re_class[dat.sum.tot2$re_class=="Rodentia"] <- "Rodents"
dat.sum.tot2$re_class[dat.sum.tot2$re_class=="Chiroptera"] <- "Bats"
dat.sum.tot2$re_class[dat.sum.tot2$re_class=="Artiodactyla"] <- "Even-toed Ungulates"
dat.sum.tot2$re_class[dat.sum.tot2$re_class=="Perissodactyla"] <- "Odd-toed Ungulates"
#dat.sum.tot2$re_class[dat.sum.tot$re_class=="Carnivora"] <- "Carnivores"
#dat.sum.tot$re_class[dat.sum.tot$re_class=="Perissodactyla" | dat.sum.tot$re_class=="Artiodactyla"] <- "Ungulates"



unique(dat.sum.tot2$re_class)
# [1] "Flying Birds"        "Even-toed Ungulates" "Primates"            "Non-Flying Birds"    "Carnivora"          
# [6] "Rodents"             "Odd-toed Ungulates"  "Bats"   

dat.simp <- ddply(dat.sum.tot2, .(re_class), summarize, N=length(re_class))
dat.simp
# re_class  N
# 1                Bats 13
# 2           Carnivora 26
# 3 Even-toed Ungulates 69
# 4        Flying Birds 82
# 5    Non-Flying Birds 11
# 6  Odd-toed Ungulates 15
# 7            Primates 50
# 8             Rodents 28

#convert transit time to hours
#dat.sum.tot$transit_hrs = dat.sum.tot$transit/60
dat.sum.tot2$MRT_hrs = dat.sum.tot2$MRT/60
#convert mass to kg
dat.sum.tot2$mass_kg = as.numeric(dat.sum.tot2$avg_mass)/1000

#write.csv(dat.sum.tot, "data/R_cleaned_data.csv")

#check how transit time varies by mass
#generally, we see that transit is longer at higher mass
#all taxa agree with this trend, though reptiles and bats are a little longer in time
#than we might otherwise predict -- this is likely due to an effect of temperature that is not
#reported here.

write.csv(dat.sum.tot2, "data/revision-3/dat_sum_tot_clean_MRT_R2.csv") #I don't know why there are some NA so manually fixing and re-uploading
# dat.sum.tot_clean <- read.csv(file = paste0(homewd, "data/dat.sum.tot.csv"), header = T, stringsAsFactors = F )
#dat.sum.tot_clean <- dat.sum.tot
