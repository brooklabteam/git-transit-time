rm(list=ls())

#install.packages("TMB", type = "source")
#install.packages("glmmTMB")
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
#homewd <- "/Users/emilyruhs/Desktop/UChi_Brook_Lab/GitHub_repos/git-transit-time/"
setwd(homewd)

#load the GIT transit data:
#dat <- read.csv(file = paste0(homewd, "/data/final-GIT-transit-database-april-2021.csv"), header = T, stringsAsFactors = F )

#this csv is the database used in the first submission + updated with the papers recommended by R2
dat <- read.csv(file = paste0(homewd, "/data/revision-1/McFerrin_database_R1.csv"), header = T, stringsAsFactors = F )
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



#now get one entry for each species/food category combination
dat.split <- dlply(dat.plot, .(fly,class,order,family,genus.species, common.name, typical.diet))
#making a function


summarise.dat <- function(dat){
  dat2 <- ddply(dat, .(fly, class, order, family, genus.species, common.name), summarise, 
                sum_mass = sum(mass), 
                typical.diet = unique(typical.diet),  
                N_tot = sum(N_individuals), 
                total_transit = if (all(is.na(total_transit))) NA else sum(total_transit, na.rm = TRUE),
                total_MRT = if (all(is.na(total_MRT))) NA else sum(total_MRT, na.rm = TRUE))
  return(dat2)
}


dat.out <- lapply(dat.split, summarise.dat)
dat.sum.tot <- data.table::rbindlist(dat.out)
head(dat.sum.tot)
dat.sum.tot$transit <- dat.sum.tot$total_transit/dat.sum.tot$N_tot
dat.sum.tot$MRT <- dat.sum.tot$total_MRT/dat.sum.tot$N_tot
dat.sum.tot$avg_mass <- dat.sum.tot$sum_mass/dat.sum.tot$N_tot
subset(dat.sum.tot, !is.na(transit)) #129
subset(dat.sum.tot, !is.na(total_MRT)) #75


#now categorize -- vertebrates are signified by class but mammals by order
dat.sum.tot$re_class <- dat.sum.tot$class
dat.sum.tot$re_class[dat.sum.tot$re_class=="Mammalia"] <- dat.sum.tot$order[dat.sum.tot$re_class=="Mammalia"]
dat.sum.tot$re_class[dat.sum.tot$re_class=="Aves" & dat.sum.tot$fly=="Y"] <- "Flying Birds"
dat.sum.tot$re_class[dat.sum.tot$re_class=="Aves" & dat.sum.tot$fly=="N"] <- "Non-Flying Birds"

unique(dat.sum.tot$re_class)
#find out how many entries per each "cat"
# dat.simp.big <- ddply(dat.sum.tot, .(re_class, typical.diet), summarize, N=length(re_class))
# dat.simp.big
dat.simp <- ddply(dat.sum.tot, .(re_class), summarize, N=length(re_class))
dat.simp

paper.dat <- ddply(dat.sum.tot, .(re_class), summarise, N_species = length(unique(genus.species)))
paper.dat

# re_class N_species
# 1          Amphibia         1
# 2      Artiodactyla        16
# 3         Carnivora         6
# 4           Cetacea         1
# 5        Chiroptera        36
# 6    Chondrichthyes         2
# 7        Dermoptera         1
# 8     Diprotodontia         3
# 9      Flying Birds        14
# 10       Lagomorpha         1
# 11 Non-Flying Birds         4
# 12  Peramelemorphia         2
# 13   Perissodactyla         1
# 14           Pilosa         1
# 15         Primates        21
# 16         Reptilia        19
# 17         Rodentia        14
# 18          Sirenia         1

write.csv(dat.sum.tot, "data/dat_sum_tot.csv")

#and group
#remove any non-mammalian classes and any mammalian orders with < 4 entries
#Amphibia, Cetacea, Chondrichtythes, Dermoptera, Pilosa, Sirenia

dat.sum.tot = subset(dat.sum.tot, re_class !="Amphibia" & 
                       re_class != "Cetacea" & 
                       re_class != "Chondrichthyes" & 
                       re_class != "Dermoptera" & 
                       re_class != "Diprotodontia" & 
                       re_class != "Lagomorpha" &
                       re_class != "Peramelemorphia" &
                       re_class != "Perissodactyla" &
                       re_class != "Pilosa" & 
                       re_class != "Sirenia") 

#&re_class != "Non-Flying Birds")

unique(dat.sum.tot$re_class)
#[1] "Primates"         "Rodentia"         "Non-Flying Birds" "Reptilia"         "Artiodactyla"     "Chiroptera"       "Carnivora"        "Flying Birds"   

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
dat.simp <- ddply(dat.sum.tot, .(re_class), summarize, N=length(re_class))
dat.simp
# re_class  N
# 1             Bats 36
# 2       Carnivores  6
# 3     Flying Birds 14
# 4 Non-Flying Birds  4
# 5         Primates 21
# 6         Reptiles 19
# 7          Rodents 14
# 8        Ungulates 16

#convert transit time to hours
dat.sum.tot$transit_hrs = dat.sum.tot$transit/60
dat.sum.tot$MRT_hrs = dat.sum.tot$MRT/60
#convert mass to kg
dat.sum.tot$mass_kg = as.numeric(dat.sum.tot$avg_mass)/1000

#write.csv(dat.sum.tot, "data/R_cleaned_data.csv")

#check how transit time varies by mass
#generally, we see that transit is longer at higher mass
#all taxa agree with this trend, though reptiles and bats are a little longer in time
#than we might otherwise predict -- this is likely due to an effect of temperature that is not
#reported here.

write.csv(dat.sum.tot, "data/dat_sum_tot_clean.csv") #I don't know why there are some NA so manually fixing and re-uploading
# dat.sum.tot_clean <- read.csv(file = paste0(homewd, "data/dat.sum.tot.csv"), header = T, stringsAsFactors = F )
dat.sum.tot_clean <- dat.sum.tot
