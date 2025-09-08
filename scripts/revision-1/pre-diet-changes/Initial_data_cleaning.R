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
#homewd= "/Users/katherinemcferrin/Developer/git-transit-time"
homewd <- "/Users/emilyruhs/Desktop/UChi_Brook_Lab/GitHub_repos/git-transit-time/"
setwd(homewd)

#load the GIT transit data:
#dat <- read.csv(file = paste0(homewd, "/data/final-GIT-transit-database-april-2021.csv"), header = T, stringsAsFactors = F )
dat <- read.csv(file = paste0(homewd, "/data/McFerrin_database_R1.csv"), header = T, stringsAsFactors = F )
#View(dat)

length(unique(dat$Retention.Citation)) #127 unique papers
length(unique(dat$Genus.species)) #144 unique species
sort(unique(dat$Genus.species))

#We also added the collection of MRT, but never plotted it
#I'm pulling that out here to include in the actual paper - starts on line 561


#dat = only the columns we will use for modeling and plotting and renaming columns
#original code: dat <- dplyr::select(dat, can_fly, Class, Order, Family, Genus.species, Common.Name, Typical.Diet, Mass..mean.median.from.study., Feeding.Trial.Food.Cat, N_measured, N_trials, Minimum..min., Median..min., Mean..min.,Maximum..min.)
dat <- dplyr::select(dat, can_fly, Class, Order, phylo_dist,Family, Genus.species, Common.Name, Typical_Diet, Mass.g..mean.median.from.study.,  N_individuals, N_trials, Minimum..min., Median..min., Mean..min.,Maximum..min., mrt..min.,mrt_sd, mrt_se)
head(dat)
#rename columns
names(dat) <- c("fly", "class", "order", "phylo_dist","family", "genus.species", "common.name", "typical.diet", "mass", "N_individuals", "N_trials", "min", "median", "mean", "max", "MRT_min", "MRT_sd", "MRT_se")
#View(dat)


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
dat.plot$N_individuals = as.numeric(dat.plot$N_individuals)
dat.plot$total_transit = dat.plot$transit*dat.plot$N_individuals #multiplying transit time x #individuals

#and also add for MRT
dat.plot$total_MRT = dat.plot$MRT_min*dat.plot$N_individuals


# create new dataset that only has MRT and remove NAs
dat.plot2 <- dat.plot
dat.plot2 <- dplyr::select(dat.plot2, -(total_transit), -(transit)) 
dat.plot2 <- dat.plot2 %>% filter(!is.na(total_MRT))


# remove entries that are NA for git transit
dat.plot <- dplyr::select(dat.plot, -(total_MRT), -(MRT_min), -(MRT_sd), -(MRT_se)) 
dat.plot <- dat.plot %>% filter(!is.na(total_transit))


################## CLEANING FOR GIT TRANSIT ONLY ##########################
#now get one entry for each species for GIT_transit
dat.split <- dlply(dat.plot, .(fly,class,order,family,genus.species, common.name, typical.diet))
#making a function


summarise.dat <- function(dat){
  dat2 <- ddply(dat, .(fly, class, order, family, genus.species, common.name), summarise, 
                sum_mass = sum(mass), 
                typical.diet = unique(typical.diet),  
                N_tot = sum(N_individuals), 
                total_transit = if (all(is.na(total_transit))) NA else sum(total_transit, na.rm = TRUE))
                #total_MRT = if (all(is.na(total_MRT))) NA else sum(total_MRT, na.rm = TRUE))
  return(dat2)
}


dat.out <- lapply(dat.split, summarise.dat)
dat.sum.tot <- data.table::rbindlist(dat.out)
head(dat.sum.tot)
dat.sum.tot$transit <- dat.sum.tot$total_transit/dat.sum.tot$N_tot
#dat.sum.tot$MRT <- dat.sum.tot$total_MRT/dat.sum.tot$N_tot
dat.sum.tot$avg_mass <- dat.sum.tot$sum_mass/dat.sum.tot$N_tot
subset(dat.sum.tot, !is.na(transit)) #131
#subset(dat.sum.tot, !is.na(total_MRT)) #75


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
# 2      Artiodactyla         7
# 3         Carnivora         6
# 4           Cetacea         1
# 5        Chiroptera        36
# 6    Chondrichthyes         2
# 7        Dermoptera         1
# 8     Diprotodontia         3
# 9      Flying Birds        14
# 10       Lagomorpha         1
# 11 Non-Flying Birds         2
# 12  Peramelemorphia         2
# 13   Perissodactyla         1
# 14         Primates        21
# 15         Reptilia        19
# 16         Rodentia        13
# 17          Sirenia         1

write.csv(dat.sum.tot, "data/dat_sum_tot_GIT.csv")

#and group
#remove any non-mammalian classes and any mammalian orders with < 2 entries
#Amphibia, Cetacea, Chondrichtythes, Dermoptera, Pilosa, Sirenia

dat.sum.tot = subset(dat.sum.tot, re_class !="Amphibia" & 
                       re_class != "Cetacea" & 
                       re_class != "Chondrichthyes" & 
                       re_class != "Dermoptera" & 
                       re_class != "Diprotodontia" & 
                       re_class != "Lagomorpha" &
                       re_class != "Non-Flying Birds" &
                       re_class != "Peramelemorphia" &
                       re_class != "Perissodactyla" &
                       re_class != "Pilosa" & 
                       re_class != "Sirenia") 


unique(dat.sum.tot$re_class)
#"Artiodactyla" "Carnivora"    "Primates"     "Rodentia"     "Reptilia"     "Flying Birds" "Chiroptera" 

#and write over some of the others
dat.sum.tot$re_class[dat.sum.tot$re_class=="Reptilia"] <- "Reptiles"
dat.sum.tot$re_class[dat.sum.tot$re_class=="Rodentia"] <- "Rodents"
dat.sum.tot$re_class[dat.sum.tot$re_class=="Chiroptera"] <- "Bats"
dat.sum.tot$re_class[dat.sum.tot$re_class=="Artiodactyla"] <- "Ungulates"
dat.sum.tot$re_class[dat.sum.tot$re_class=="Carnivora"] <- "Carnivores"
#dat.sum.tot$re_class[dat.sum.tot$re_class=="Perissodactyla" | dat.sum.tot$re_class=="Artiodactyla"] <- "Ungulates"

#dat.sum.tot$re_class[dat.sum.tot$re_class=="Diprotodontia" | dat.sum.tot$re_class=="Peramelemorphia"] <- "Marsupials"

#dat.sum.tot$re_class[dat.sum.tot$re_class=="Lagomorpha" ] <- "Lagomorphs"

unique(dat.sum.tot$re_class)
# [1] "Ungulates"    "Carnivores"   "Primates"     "Rodents"      "Reptiles"     "Flying Birds" "Bats" 

dat.simp <- ddply(dat.sum.tot, .(re_class), summarize, N=length(re_class))
dat.simp
# re_class  N
# 6      Rodents 13
# 3 Flying Birds 14
# 1         Bats 36
# 2   Carnivores  6
# 4     Primates 21
# 7    Ungulates  7
# 5     Reptiles 19



#convert transit time to hours
dat.sum.tot$transit_hrs = dat.sum.tot$transit/60
#dat.sum.tot$MRT_hrs = dat.sum.tot$MRT/60
#convert mass to kg
dat.sum.tot$mass_kg = as.numeric(dat.sum.tot$avg_mass)/1000

#write.csv(dat.sum.tot, "data/R_cleaned_data.csv")

#check how transit time varies by mass
#generally, we see that transit is longer at higher mass
#all taxa agree with this trend, though reptiles and bats are a little longer in time
#than we might otherwise predict -- this is likely due to an effect of temperature that is not
#reported here.

write.csv(dat.sum.tot, "data/dat_sum_tot_clean_GIT.csv") #I don't know why there are some NA so manually fixing and re-uploading
# dat.sum.tot_clean <- read.csv(file = paste0(homewd, "data/dat.sum.tot.csv"), header = T, stringsAsFactors = F )
#dat.sum.tot_clean <- dat.sum.tot










################## CLEANING FOR MRT ONLY ##########################
#now get one entry for each species for GIT_transit
dat.split <- dlply(dat.plot2, .(fly,class,order,family,genus.species, common.name, typical.diet))
#making a function


summarise.dat2 <- function(dat){
  dat2 <- ddply(dat, .(fly, class, order, family, genus.species, common.name), summarise, 
                sum_mass = sum(mass), 
                typical.diet = unique(typical.diet),  
                N_tot = sum(N_individuals), 
                #total_transit = if (all(is.na(total_transit))) NA else sum(total_transit, na.rm = TRUE))
                total_MRT = if (all(is.na(total_MRT))) NA else sum(total_MRT, na.rm = TRUE))
  return(dat2)
}


dat.out <- lapply(dat.split, summarise.dat2)
dat.sum.tot2 <- data.table::rbindlist(dat.out)
head(dat.sum.tot2)
#dat.sum.tot$transit <- dat.sum.tot$total_transit/dat.sum.tot$N_tot
dat.sum.tot2$MRT <- dat.sum.tot2$total_MRT/dat.sum.tot2$N_tot
dat.sum.tot2$avg_mass <- dat.sum.tot2$sum_mass/dat.sum.tot2$N_tot
#subset(dat.sum.tot, !is.na(transit)) #131
subset(dat.sum.tot2, !is.na(total_MRT)) #74


#now categorize -- vertebrates are signified by class but mammals by order
dat.sum.tot2$re_class <- dat.sum.tot2$class
dat.sum.tot2$re_class[dat.sum.tot2$re_class=="Mammalia"] <- dat.sum.tot2$order[dat.sum.tot2$re_class=="Mammalia"]
dat.sum.tot2$re_class[dat.sum.tot2$re_class=="Aves" & dat.sum.tot2$fly=="Y"] <- "Flying Birds"
dat.sum.tot2$re_class[dat.sum.tot2$re_class=="Aves" & dat.sum.tot2$fly=="N"] <- "Non-Flying Birds"

unique(dat.sum.tot2$re_class)
# [1] "Amphibia"         "Non-Flying Birds" "Artiodactyla"     "Carnivora"        "Dermoptera"       "Diprotodontia"    "Lagomorpha"       "Peramelemorphia" 
#[9] "Perissodactyla"   "Pilosa"           "Primates"         "Rodentia"         "Sirenia"          "Reptilia"         "Flying Birds"     "Chiroptera"   

paper.dat <- ddply(dat.sum.tot2, .(re_class), summarise, N_species = length(unique(genus.species)))
paper.dat

# re_class N_species
# 1          Amphibia         1
# 2      Artiodactyla        12
# 3         Carnivora         1
# 4        Chiroptera        12
# 5        Dermoptera         1
# 6     Diprotodontia         3
# 7      Flying Birds         1
# 8        Lagomorpha         1
# 9  Non-Flying Birds         3
# 10  Peramelemorphia         2
# 11   Perissodactyla         1
# 12           Pilosa         1
# 13         Primates        19
# 14         Reptilia         1
# 15         Rodentia        14
# 16          Sirenia         1

write.csv(dat.sum.tot2, "data/dat_sum_tot_MRT.csv")

#and group
#remove any non-mammalian classes and any mammalian orders with < 2 entries
#Amphibia, Cetacea, Chondrichtythes, Dermoptera, Pilosa, Sirenia

dat.sum.tot2 = subset(dat.sum.tot2, re_class !="Amphibia" & 
                       re_class != "Carnivora" &
                       re_class != "Dermoptera" & 
                       re_class != "Diprotodontia" & #marsupial
                       re_class != "Flying Birds" &
                       re_class != "Lagomorpha" &
                       re_class != "Peramelemorphia" & #marsupial
                       re_class != "Perissodactyla" &
                       re_class != "Pilosa" & 
                       re_class != "Reptilia" &
                       re_class != "Sirenia") 


unique(dat.sum.tot2$re_class)
#"Non-Flying Birds" "Artiodactyla"     "Primates"         "Rodentia"         "Flying Birds"     "Chiroptera" 

#and write over some of the others
#dat.sum.tot2$re_class[dat.sum.tot2$re_class=="Reptilia"] <- "Reptiles"
dat.sum.tot2$re_class[dat.sum.tot2$re_class=="Rodentia"] <- "Rodents"
dat.sum.tot2$re_class[dat.sum.tot2$re_class=="Chiroptera"] <- "Bats"
dat.sum.tot2$re_class[dat.sum.tot2$re_class=="Artiodactyla"] <- "Ungulates"
#dat.sum.tot2$re_class[dat.sum.tot$re_class=="Carnivora"] <- "Carnivores"
#dat.sum.tot$re_class[dat.sum.tot$re_class=="Perissodactyla" | dat.sum.tot$re_class=="Artiodactyla"] <- "Ungulates"



unique(dat.sum.tot2$re_class)
# "Non-Flying Birds"  "Ungulates"        "Primates"         "Rodents"    "Bats"  

dat.simp <- ddply(dat.sum.tot2, .(re_class), summarize, N=length(re_class))
dat.simp
# re_class  N
# 1             Bats 12
# 2     Flying Birds  1
# 3 Non-Flying Birds  3
# 4         Primates 19
# 5          Rodents 14
# 6        Ungulates 12

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

write.csv(dat.sum.tot2, "data/dat_sum_tot_clean_MRT.csv") #I don't know why there are some NA so manually fixing and re-uploading
# dat.sum.tot_clean <- read.csv(file = paste0(homewd, "data/dat.sum.tot.csv"), header = T, stringsAsFactors = F )
#dat.sum.tot_clean <- dat.sum.tot
