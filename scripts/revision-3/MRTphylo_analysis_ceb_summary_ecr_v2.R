rm(list=ls())
###Phylo analysis
# load the packages
library(lme4)
library(nlme)
library(ape)
library(dplyr)
library(plyr)
library(phytools)
library(ggplot2)
library(ggforce)

# load data
#homewd <- "/Users/carabrook/Developer/git-transit-time/"
#homewd <- "/Users/katherinemcferrin/Developer/git-transit-time/"
homewd <- "/Users/gavindehnert/Desktop/GitHub_repos/git-transit-time"

setwd(homewd)

#load the MRT transit data with phylo name
dat_MRT <- read.csv(file = paste0(homewd, "/data/dat_sum_tot_clean_MRT_R3.csv"), header = T, stringsAsFactors = F )

mass_col    <- "avg_mass"        # body-mass column
species_col <- "genus.species"  # tip / phylogeny-name column

# species missing a mass in either dataset
miss <- function(d) d[[species_col]][is.na(d[[mass_col]]) | trimws(d[[mass_col]]) == "0"]
no_mass <- sort(unique(c(miss(dat_MRT)))) 

# report which ones
cat("Species without mass (", length(no_mass), "):\n", sep = "")
print(no_mass) #0

# remove them from both datasets
dat_MRT <- dat_MRT[!dat_MRT[[species_col]] %in% no_mass, ]

#### species counts ####
#table 1 in the manuscript

length(unique(dat_MRT$genus.species)) #242 unique species for transit time
tt.pivot <-ddply(dat_MRT, .(re_class), summarise, N_species_transit = length(unique(genus.species)))
print(tt.pivot)
# re_class N_species_transit
# 1                Bats                13
# 2           Carnivora                18
# 3 Even-toed Ungulates                55
# 4        Flying Birds                72
# 5    Non-Flying Birds                 7
# 6  Odd-toed Ungulates                10
# 7            Primates                40
# 8             Rodents                27


#dat_MRT <- subset(dat_MRT, dat_MRT$re_class != "Non-Flying Birds")
#dat_MRT <- subset(dat_MRT, dat_MRT$re_class != "Reptiles")
#288

#load the tree - pulled from timetree so should be in principle an ultrametric tree
tree <- read.tree(file = paste0(homewd, "/data/Book3.nwk"))
# is.ultrametric(tree) #check for ultrametric - is FALSE
# 

tree_ultra2 <- force.ultrametric(tree, method = "nnls") #no warnings

is.ultrametric(tree_ultra2) # TRUE

### LET'S GO WITH THIS force.ultrametric APPROACH AS IT DOES NOT HAVE WARNINGS





################### MRT phylogenetic signal ##########################
dat_MRT$phylo_name <- as.character(dat_MRT$phylo_name)

# tree <- read.tree(file = paste0(homewd, "data/Book3.nwk"))
# tree_ultra2 <- force.ultrametric(tree, method = "nnls") #no warnings
# is.ultrametric(tree_ultra2)

#for MRT
MRT <- dat_MRT[c("phylo_name", "MRT_hrs")]
#MRT <- na.omit(MRT)

#make sure they match
# If you have a data frame GIT with the trait
MRT_hrs <- setNames(MRT$MRT_hrs, MRT$phylo_name)

# Then only keep the overlapping species
common_species <- intersect(tree_ultra2$tip.label, names(MRT_hrs)) #222
tree <- drop.tip(tree_ultra2, setdiff(tree_ultra2$tip.label, common_species)) 
MRT_hrs <- MRT_hrs[common_species] #223

#for MRT
lambda_gs2<-phylosig(tree, MRT_hrs,
                     method="lambda",test=TRUE)
lambda_gs2

# Phylogenetic signal lambda : 0.775571 
# logL(lambda) : -950.042 
# LR(lambda=0) : 139.153 
# P-value (based on LR test) : 4.07696e-32 



############# plotting phylo signal by tree/phenogram ############


contMap(tree,log10(MRT_hrs), fsize =0.6)




#########################################  MGT. #################################

########### MGT data cleaning #######

#_________________________________________________________________________________________________________________
#I think this code runs PGLS and also provides lambda estimates as part of the model


# Assuming that we already have things loaded from above
# Load your data 
# dat_MRT <- read.csv(file = paste0(homewd, "/data/dat_sum_tot_clean_MRT_R3.csv"), header = T, stringsAsFactors = F )
# names(dat_MRT)


# clean dataset with only the columns you need
MRT <- dat_MRT[c("phylo_name", "MRT_hrs", "trial.diet", "re_class", "mass_kg")]
#MRT <- na.omit(MRT$transit_hrs)

tree <- read.tree(paste0(homewd, "/data/Book3.nwk"))
tree_ultra2 <- force.ultrametric(tree, method = "nnls") #no warnings
is.ultrametric(tree_ultra2)

#make sure they match
# Check if the species names match the tree tips
common_species <- intersect(tree_ultra2$tip.label, MRT$phylo_name)
tree_pruned <- drop.tip(tree_ultra2, setdiff(tree_ultra2$tip.label, common_species))
MRT_pruned <- MRT[MRT$phylo_name %in% common_species, ]

# there are some duplicated taxa with different diets:
MRT_pruned[duplicated(MRT_pruned$phylo_name),]



species1 <- unique(dat_MRT$phylo_name)
species2 <- unique(MRT_pruned$phylo_name)

# Species in df1 but not in df2
missing_in_df2 <- setdiff(species1, species2)
missing_in_df2 

#[1] "Gazella"           "Nectarinia_famosa" ""

#sort(tree$tip.label)

###################### MGT - model 1 ################


# Fit PGLS model
#ML=maximum likelihood

# factor the groups so that rodent is the comparasion
unique(MRT_pruned$re_class)
# [1] "Flying Birds"        "Carnivora"           "Even-toed Ungulates" "Primates"            "Rodents"            
# [6] "Odd-toed Ungulates"  "Bats"   

length(unique(MRT_pruned$phylo_name)) #230


# MRT_pruned$re_class <- factor(MRT_pruned$re_class, levels = c("Rodents","Bats", "Non-Flying Birds", 
#                                                               "Primates", "Ungulates"))

#transform data
MRT_pruned$log_MRT_hrs <- log10(MRT_pruned$MRT_hrs)


# add the "flyer" parameter
MRT_pruned$flyer <- 0
MRT_pruned$flyer[MRT_pruned$re_class=="Bats"] <- 1

MRT_pruned$flyer <- as.factor(MRT_pruned$flyer)

# Create correlation structure
# this is full Brownian motion assumption of phylogentic effects:
cor_phylo_fixed1 <- corPagel(1, phy = tree_pruned, fixed = TRUE, form = ~phylo_name) #assumes phylogeny equivalent to brownian motion

# this is no phylogentic effects:
cor_phylo_fixed0 <- corPagel(0, phy = tree_pruned, fixed = TRUE, form = ~phylo_name) #assumes no phylogenetic signal

# this is what was estimated from your data above - removing because lambda is 1 in the above model
cor_phylo_fixed2 <- corPagel(lambda_gs2$lambda, phy = tree_pruned, fixed = TRUE, form = ~phylo_name) 



# Cara: now test the effect of flyer on MRT transit without accounting for mass or diet
# All text below is from Cara

#with full phylogenetic effects:
gls_model_MRT_1 <- gls(log_MRT_hrs ~ flyer, 
                         data = MRT_pruned, 
                        correlation = cor_phylo_fixed1,
                         method = "ML")
summary(gls_model_MRT_1)$tTable
# Value Std.Error   t-value    p-value
# (Intercept)  0.6604288 0.4035738  1.636451 0.10281518
# flyer1      -0.5484253 0.2756196 -1.989790 0.04754016

#with no phylogenetic effects:
gls_model_MRT_0 <- gls(log_MRT_hrs ~ flyer, 
                        data = MRT_pruned, 
                        #correlation = cor_phylo_fixed0,
                        method = "ML")
summary(gls_model_MRT_0)$tTable
# Value  Std.Error   t-value      p-value
# (Intercept)  1.0169502 0.04298513 23.658186 5.043775e-70
# flyer1      -0.5519944 0.20511263 -2.691177 7.527674e-03

#with estimated phylogenetic effects:
gls_model_MRT_2 <- gls(log_MRT_hrs ~ flyer,
                        data = MRT_pruned,
                        correlation = cor_phylo_fixed2,
                        method = "ML")
summary(gls_model_MRT_2)$tTable
# Value Std.Error   t-value    p-value
# (Intercept)  0.6796533 0.3652269  1.860907 0.06375544
# flyer1      -0.5955751 0.2903126 -2.051496 0.04110444

AIC(gls_model_MRT_0, gls_model_MRT_1, gls_model_MRT_2)
# Emily: The model with full brownian phylogenetic effects is the best fit.
# df        AIC
# gls_model_MRT_0  3   652.0917
# gls_model_MRT_1  3 -2521.3884
# gls_model_MRT_2  3 -2104.7056



summary(gls_model_MRT_1)$tTable
# Value Std.Error   t-value    p-value
# (Intercept)  0.6604288 0.4035738  1.636451 0.10281518
# flyer1      -0.5484253 0.2756196 -1.989790 0.04754016

# Emily: flyer is a marginally significant term which says that flight is significantly
# negatively related to MRT transit time.
# This is on top of a model that already includes significant 
# phylogenetic clustering in transit time, suggesting that flight affects transit
# time independent of phylogeny.



# now, let's also account for diet
# first try as a fixed effect:
MRT_pruned$trial.diet <- as.factor(MRT_pruned$trial.diet)

# Cara: with full phylogenetic effects:
gls_model_MRT_1_diet <- gls(log_MRT_hrs ~ flyer + trial.diet, 
                        data = MRT_pruned, 
                        correlation = cor_phylo_fixed1,
                        method = "ML")
summary(gls_model_MRT_1_diet)$tTable
# Value  Std.Error    t-value    p-value
# (Intercept)                    0.66605423 0.40665139  1.6378998 0.10253466
# flyer1                        -0.60778820 0.27510048 -2.2093316 0.02793848
# trial.dietfruit/nectar/pollen -0.05727265 0.08362692 -0.6848590 0.49398369
# trial.dietmeat                 0.35573397 0.13043309  2.7273292 0.00677686
# trial.dietmixed                0.07023159 0.33187695  0.2116194 0.83255368
# trial.dietprotein              0.12845194 0.11269389  1.1398305 0.25530398
# trial.dietunk                  0.19216605 0.47080970  0.4081608 0.68345893
# trial.dietunknown              0.49175591 0.24176613  2.0340149 0.04286635



# With no phylogenetic effects:
gls_model_MRT_0_diet <- gls(log_MRT_hrs ~ flyer + trial.diet, 
                        data = MRT_pruned, 
                        #correlation = cor_phylo_fixed0,
                        method = "ML")
summary(gls_model_MRT_0_diet)$tTable
# Value  Std.Error    t-value      p-value
# (Intercept)                    1.12735456 0.05506415 20.4734768 4.037345e-58
# flyer1                        -0.36571946 0.23001371 -1.5899898 1.129342e-01
# trial.dietfruit/nectar/pollen -0.41290758 0.10757484 -3.8383287 1.523051e-04
# trial.dietmeat                -0.12499968 0.13026829 -0.9595557 3.380837e-01
# trial.dietmixed               -0.07043189 0.50316702 -0.1399772 8.887759e-01
# trial.dietprotein             -0.29060165 0.16301632 -1.7826537 7.569596e-02
# trial.dietunk                  0.71774348 0.70945192  1.0116873 3.125368e-01
# trial.dietunknown              0.40270659 0.35791696  1.1251397 2.614666e-01


# With estimated phylogenetic effects:
gls_model_MRT_2_diet <- gls(log_MRT_hrs ~ flyer + trial.diet,
                        data = MRT_pruned,
                        correlation = cor_phylo_fixed2,
                        method = "ML")
summary(gls_model_MRT_2_diet)$tTable
# Value  Std.Error     t-value     p-value
# (Intercept)                    0.68191203 0.36910653  1.84746669 0.065704960
# flyer1                        -0.61915229 0.28946755 -2.13893506 0.033281848
# trial.dietfruit/nectar/pollen -0.07335471 0.08517432 -0.86123035 0.389827636
# trial.dietmeat                 0.35290954 0.12902113  2.73528476 0.006619321
# trial.dietmixed                0.02624841 0.32888095  0.07981129 0.936442777
# trial.dietprotein              0.09662541 0.11613324  0.83202197 0.406085890
# trial.dietunk                  0.30857524 0.48945305  0.63044912 0.528900706
# trial.dietunknown              0.47511909 0.24686532  1.92460849 0.055263338

AIC(gls_model_MRT_0, gls_model_MRT_1, gls_model_MRT_2, 
    gls_model_MRT_0_diet, gls_model_MRT_1_diet, gls_model_MRT_2_diet)
# df        AIC
# gls_model_MRT_0       3   652.0917
# gls_model_MRT_1       3 -2521.3884
# gls_model_MRT_2       3 -2104.7056
# gls_model_MRT_0_diet  9   644.9016
# gls_model_MRT_1_diet  9 -2523.4310
# gls_model_MRT_2_diet  9 -2106.6616


summary(gls_model_MRT_1_diet)$tTable
# Value  Std.Error    t-value    p-value
# (Intercept)                    0.66605423 0.40665139  1.6378998 0.10253466
# flyer1                        -0.60778820 0.27510048 -2.2093316 0.02793848
# trial.dietfruit/nectar/pollen -0.05727265 0.08362692 -0.6848590 0.49398369
# trial.dietmeat                 0.35573397 0.13043309  2.7273292 0.00677686
# trial.dietmixed                0.07023159 0.33187695  0.2116194 0.83255368
# trial.dietprotein              0.12845194 0.11269389  1.1398305 0.25530398
# trial.dietunk                  0.19216605 0.47080970  0.4081608 0.68345893
# trial.dietunknown              0.49175591 0.24176613  2.0340149 0.04286635

# Emily: model gls_model_MRT_1_diet has a close AIC to the model without diet
# Plus only meat comes out as significant when you do the summary output. 
# So prefer the simpler model?



# Cara: let's try it also as a random effect - here we use the lme function which
# is also supported in the nlme package. we can compare these since all were fit via ML:

gls_model_MRT_1_diet_re <- lme(log_MRT_hrs ~ flyer,
                            random = ~1 | trial.diet,
                            correlation = cor_phylo_fixed1,
                            data = MRT_pruned,
                            method = "ML")
summary(gls_model_MRT_1_diet_re)$tTable
# Value Std.Error  DF   t-value      p-value
# (Intercept)  0.9876054 0.1958407 288  5.042902 8.124496e-07
# flyer1      -0.7264127 0.2370835 288 -3.063952 2.391159e-03

gls_model_MRT_0_diet_re <- lme(log_MRT_hrs ~ flyer,
                                random = ~1 | trial.diet,
                                #correlation = cor_phylo_fixed0,
                                data = MRT_pruned,
                                method = "ML")
summary(gls_model_MRT_0_diet_re)$tTable
# Value Std.Error  DF   t-value      p-value
# (Intercept)  0.9840423 0.0862259 288 11.412374 3.880157e-25
# flyer1      -0.4382938 0.2176507 288 -2.013749 4.496575e-02

gls_model_MRT_2_diet_re <- lme(log_MRT_hrs ~ flyer,
                                random = ~1 | trial.diet,
                                correlation = cor_phylo_fixed2,
                                data = MRT_pruned,
                                method = "ML")
summary(gls_model_MRT_2_diet_re)$tTable
# Value Std.Error  DF   t-value      p-value
# (Intercept)  0.9180366 0.1737722 288  5.282987 2.511246e-07
# flyer1      -0.6805407 0.2167721 288 -3.139429 1.868579e-03

AIC(gls_model_MRT_0, gls_model_MRT_0_diet, gls_model_MRT_0_diet_re,
    gls_model_MRT_1, gls_model_MRT_1_diet, gls_model_MRT_1_diet_re,
    gls_model_MRT_2, gls_model_MRT_2_diet, gls_model_MRT_2_diet_re)

# df        AIC
# gls_model_MRT_0          3   652.0917
# gls_model_MRT_0_diet     9   644.9016
# gls_model_MRT_0_diet_re  4   647.5806
# gls_model_MRT_1          3 -2521.3884
# gls_model_MRT_1_diet     9 -2523.4310
# gls_model_MRT_1_diet_re  4 -1471.6609
# gls_model_MRT_2          3 -2104.7056
# gls_model_MRT_2_diet     9 -2106.6616
# gls_model_MRT_2_diet_re  4 -1140.6559

#Emily: in all cases, the fit with diet as a fixed effect was better supported than the fit as
# a random effect. 
# Cara: The models with the lambda=1, which incorporate phylogenetic clustering
# of MRT transit time, are best supported. The top fit model, gls_model_MRT_1,
# offers the best fit overall.




###################
########## model gls_model_MRT_1 is best for reporting.




# now, let's also consider the effects of mass and, later, mass and diet
 
#with full phylogenetic effects:
gls_model_MRT_1_mass <- gls(log_MRT_hrs ~ log10(mass_kg)*flyer, 
                        data = MRT_pruned, 
                        correlation = cor_phylo_fixed1,
                        method = "ML")
summary(gls_model_MRT_1_mass)$tTable
# Value Std.Error   t-value     p-value
# (Intercept)            2.2910150 1.1879302  1.928577 0.054752610
# log10(mass_kg)         0.5816735 0.3577441  1.625949 0.105039299
# flyer1                -2.1598011 0.6749336 -3.200020 0.001525552
# log10(mass_kg):flyer1 -0.6939963 0.2561255 -2.709595 0.007134686


#with no phylogenetic effects:
gls_model_MRT_0_mass <- gls(log_MRT_hrs ~ log10(mass_kg)*flyer, 
                        data = MRT_pruned, 
                        #correlation = cor_phylo_fixed0,
                        method = "ML")
summary(gls_model_MRT_0_mass)$tTable
# Value  Std.Error   t-value      p-value
# (Intercept)            0.9982601 0.03504252 28.487107 2.724862e-86
# log10(mass_kg)         0.2919655 0.02386620 12.233430 4.618833e-28
# flyer1                -1.4117490 0.72219206 -1.954811 5.155934e-02
# log10(mass_kg):flyer1 -0.6395551 0.27903411 -2.292032 2.261525e-02

#with estimated phylogenetic effects:
gls_model_MRT_2_mass <- gls(log_MRT_hrs ~ log10(mass_kg)*flyer,
                        data = MRT_pruned,
                        correlation = cor_phylo_fixed2,
                        method = "ML")
summary(gls_model_MRT_2_mass)$tTable
# Value Std.Error   t-value     p-value
# (Intercept)            2.1286830 1.0254755  2.075801 0.038787357
# log10(mass_kg)         0.5332326 0.3102336  1.718810 0.086708834
# flyer1                -2.1195855 0.6892104 -3.075382 0.002301379
# log10(mass_kg):flyer1 -0.6913428 0.2703146 -2.557549 0.011046487



AIC(gls_model_MRT_0_mass, gls_model_MRT_1_mass, gls_model_MRT_2_mass)
#Emily: model 1 is the best, by a lot
#Cara: This shows that transit time is still influenced by phylogeny, not
# exclusively by mass (e.g. the model with phylogeny AND mass performed better
# than the model with mass alone)

# df        AIC
# gls_model_MRT_0_mass  5   532.5672
# gls_model_MRT_1_mass  5 -2526.9941 **
# gls_model_MRT_2_mass  5 -2109.4875

summary(gls_model_MRT_1_mass)$tTable

# Value Std.Error   t-value     p-value
# (Intercept)            2.2910150 1.1879302  1.928577 0.054752610
# log10(mass_kg)         0.5816735 0.3577441  1.625949 0.105039299
# flyer1                -2.1598011 0.6749336 -3.200020 0.001525552
# log10(mass_kg):flyer1 -0.6939963 0.2561255 -2.709595 0.007134686

# Emily: flyer is significant (negative) term here (p=0.005)
# We can test whether it would be best to drop the interaction which 
# might clarify the variables that matter most.


gls_model_MRT_1_mass_simple <- gls(log_MRT_hrs ~ flyer + log10(mass_kg), 
                             data = MRT_pruned, 
                             correlation = cor_phylo_fixed1,
                             method = "ML")
summary(gls_model_MRT_1_mass_simple)$tTable
# Value Std.Error   t-value    p-value
# (Intercept)     2.3533870 1.2004921  1.960352 0.05090206
# flyer1         -0.4895454 0.2778380 -1.761982 0.07911518
# log10(mass_kg)  0.5408144 0.3612736  1.496966 0.13547851


AIC(gls_model_MRT_1_mass, gls_model_MRT_1_mass_simple)

# df       AIC
# gls_model_MRT_1_mass         5 -2526.994
# gls_model_MRT_1_mass_simple  4 -2521.644


# Emily: there is no difference in AIC,really, though the complex model is slightly better

summary(gls_model_MRT_1_mass_simple)$tTable
# Value Std.Error   t-value    p-value
# (Intercept)     2.3533870 1.2004921  1.960352 0.05090206
# flyer1         -0.4895454 0.2778380 -1.761982 0.07911518
# log10(mass_kg)  0.5408144 0.3612736  1.496966 0.13547851


# Emily: note that flight is significant in the complex model.

# now, let's see if adding diet improves this model:


# I'll just illustrate with the best fit model
# first, as a fixed effect
gls_model_MRT_1_mass_simple_diet <- gls(log_MRT_hrs ~ flyer + log10(mass_kg) +
                                          trial.diet, 
                                    data = MRT_pruned, 
                                    correlation = cor_phylo_fixed1,
                                    method = "ML")
summary(gls_model_MRT_1_mass_simple_diet)$tTable
# Value  Std.Error    t-value     p-value
# (Intercept)                    2.49585577 1.18697299  2.1027065 0.036362231
# flyer1                        -0.54650667 0.27682947 -1.9741637 0.049321801
# log10(mass_kg)                 0.58506942 0.35669879  1.6402338 0.102052071
# trial.dietfruit/nectar/pollen -0.05654948 0.08338375 -0.6781835 0.498201752
# trial.dietmeat                 0.36475619 0.13016826  2.8021899 0.005421111
# trial.dietmixed                0.06804594 0.33090995  0.2056328 0.837223547
# trial.dietprotein              0.13368176 0.11240985  1.1892353 0.235330040
# trial.dietunk                  0.18448182 0.46945746  0.3929681 0.694634742
# trial.dietunknown              0.48696436 0.24107744  2.0199499 0.044317542

# then, as a random effect
gls_model_MRT_1_mass_simple_diet_re <- lme(log_MRT_hrs ~ flyer + log10(mass_kg),
                                            random = ~1 | trial.diet,
                                            data = MRT_pruned, 
                                            correlation = cor_phylo_fixed1,
                                            method = "ML")
summary(gls_model_MRT_1_mass_simple_diet_re)$tTable
# Value  Std.Error  DF   t-value      p-value
# (Intercept)     0.99955555 0.17345347 287  5.762672 2.133193e-08
# flyer1         -0.62399048 0.24564121 287 -2.540252 1.160448e-02
# log10(mass_kg)  0.05330623 0.03594279 287  1.483086 1.391490e-01

# then, compare
AIC(gls_model_MRT_1_mass_simple, gls_model_MRT_1_mass_simple_diet, gls_model_MRT_1_mass_simple_diet_re)

# df       AIC
# gls_model_MRT_1_mass_simple          4 -2521.644
# gls_model_MRT_1_mass_simple_diet    10 -2524.193
# gls_model_MRT_1_mass_simple_diet_re  5 -1471.657

# here, Adding diet does nothing really for the fit. So only have mass in the model is still the best

summary(gls_model_MRT_1_mass_simple)$tTable

# Value Std.Error   t-value    p-value
# (Intercept)     2.3533870 1.2004921  1.960352 0.05090206
# flyer1         -0.4895454 0.2778380 -1.761982 0.07911518
# log10(mass_kg)  0.5408144 0.3612736  1.496966 0.13547851

# Emily: flyer is not significantly negatively related to transit, in the simple model
# after accounting for mass (and phylogeny). 


# Cara: what if we compare our two best fit models so far? 
AIC(gls_model_MRT_1, gls_model_MRT_1_diet, gls_model_MRT_1_mass_simple)

# df       AIC
# gls_model_MRT_1              3 -2521.388
# gls_model_MRT_1_diet         9 -2523.431
# gls_model_MRT_1_mass_simple  4 -2521.644


gls_model_MRT_1_mass_only <- gls(log_MRT_hrs ~ log10(mass_kg), 
                                   data = MRT_pruned, 
                                   correlation = cor_phylo_fixed1,
                                   method = "ML")
summary(gls_model_MRT_1_mass_only)$tTable
# Value Std.Error  t-value    p-value
# (Intercept)    2.4694245 1.2029671 2.052778 0.04097899
# log10(mass_kg) 0.6309302 0.3589129 1.757892 0.07980647

gls_model_MRT_0_mass_only <- gls(log_MRT_hrs ~ log10(mass_kg), 
                                 data = MRT_pruned, 
                                 #correlation = cor_phylo_fixed0,
                                 method = "ML")
summary(gls_model_MRT_0_mass_only)$tTable
# Value Std.Error  t-value    p-value
# (Intercept)    2.2597950 1.0358953 2.181490 0.02993917
# log10(mass_kg) 0.5571227 0.3086132 1.805246 0.07205918

# mass not significant

gls_model_MRT_2_mass_only <- gls(log_MRT_hrs ~ log10(mass_kg), 
                                 data = MRT_pruned, 
                                 correlation = cor_phylo_fixed2,
                                 method = "ML")
summary(gls_model_MRT_2_mass_only)$tTable
# Value Std.Error  t-value    p-value
# (Intercept)    2.2597950 1.0358953 2.181490 0.02993917
# log10(mass_kg) 0.5571227 0.3086132 1.805246 0.07205918

AIC(gls_model_MRT_0_mass_only, gls_model_MRT_1_mass_only, gls_model_MRT_2_mass_only)
# 
# df        AIC
# gls_model_MRT_0_mass_only  3   535.0061
# gls_model_MRT_1_mass_only  3 -2520.5238
# gls_model_MRT_2_mass_only  3 -2103.7614








###### ------ PLOTTING FOR FIGURE S2
library(nlme)
library(ggplot2)
library(dplyr)

unique(MRT_pruned$re_class)
# [1] "Flying Birds"        "Carnivora"           "Even-toed Ungulates" "Primates"            "Non-Flying Birds"   
# [6] "Rodents"             "Odd-toed Ungulates"  "Bats"             

colz <- c(
  "Flying\nBirds"= "#F8766D",
  "Bats" = "#C49A00",
  "Non-Flying\nBirds" = "#edf8b1",
  "Rodents" = "navy",
  "Primates" = "#00B6EB",
  "Even-toed\nUngulates" = "#A58AFF",
  "Reptiles" = "#FB61D7",
  "Carnivores" = "#00C094","Marsupials"="#E08B45", "Shrews/Moles"="#c51b8c", "Odd-toed\nUngulates"="thistle1"
)

unique(MRT_pruned$trial.diet)
# [1] meat                fruit/nectar/pollen fiber/foliage       protein             unknown             mixed              
# Levels: fiber/foliage fruit/nectar/pollen meat mixed protein unknown

shapez <- c(
  "meat" = 25,
  "fruit/nectar/pollen" = 22,
  "fiber/foliage" = 23,
  "protein" = 21,
  "unknown" = 24,
  "mixed" = 8
)

# --- Fix factor levels ---
# Convert to character to safely modify
MRT_pruned$re_class <- as.character(MRT_pruned$re_class)

# Replace the value
MRT_pruned$re_class[MRT_pruned$re_class == "Flying Birds"] <- "Flying\nBirds"
MRT_pruned$re_class[MRT_pruned$re_class == "Non-Flying Birds"] <- "Non-Flying\nBirds"
MRT_pruned$re_class[MRT_pruned$re_class == "Even-toed Ungulates"] <- "Even-toed\nUngulates"
MRT_pruned$re_class[MRT_pruned$re_class == "Odd-toed Ungulates"] <- "Odd-toed\nUngulates"
MRT_pruned$re_class[MRT_pruned$re_class == "Carnivora"] <- "Carnivores"

unique(MRT_pruned$re_class)
# [1] "Flying Birds"         "Even-toed\nUngulates" "Primates"             "Marsupials"           "Carnivores"          
# [6] "Rodents"              "Odd-toed\nUngulates"  "Bats"   

MRT_pruned$re_class <- factor(MRT_pruned$re_class, levels = c("Flying\nBirds", "Bats", "Rodents", "Non-Flying\nBirds", 
                                                              "Carnivores","Primates",
                                                              "Even-toed\nUngulates", "Odd-toed\nUngulates"))

# --- Plot ---
p1 <- ggplot(MRT_pruned) +
  geom_boxplot(aes(x = re_class, y = log10(MRT_hrs), fill = re_class), show.legend = FALSE, outlier.shape = NA) +
  geom_point(aes(x = re_class, y = log10(MRT_hrs), shape = trial.diet),
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
             size = 3) +
  scale_shape_manual(values = shapez, name = "Trial diet") +
  scale_fill_manual(values = colz) +
  ylab(expression(Log[10] ~ "Mean retention time (hrs)")) +
  coord_cartesian(ylim = c(-2, 2.95)) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid = element_blank(),
    legend.position = c(0.1, 0.85),
    legend.background = element_rect(color = "black"),
    plot.margin = unit(c(0.2, 0.2, 0.8, 0.2), "cm")
  )

print(p1)


ggsave(file = paste0(homewd,"/figures_r4/Fig_S1A.png"),
       units="mm",  
       width=120, 
       height=80, 
       scale=3, 
       dpi=300)




#### panel B - mass - raw data
# 
# # --- Get predictions from the model ---
# library(nlme)
# library(ggplot2)
# library(dplyr)
# 
# # --- Prepare dataset ---
# MRT_pruned <- MRT_pruned %>%
#   mutate(log_mass = log10(mass_kg))
# 
# # --- Fit GLS model ---
# gls_model <- gls(
#   log_MRT_hrs ~ flyer + log_mass,
#   data = MRT_pruned,
#   correlation = cor_phylo_fixed1,
#   method = "ML"
# )
# 
# # --- Build prediction grid for lines ---
# # Use full range of log_mass for each flyer level
# newdat <- expand.grid(
#   log_mass = seq(min(MRT_pruned$log_mass, na.rm = TRUE),
#                  max(MRT_pruned$log_mass, na.rm = TRUE),
#                  length.out = 100),
#   flyer = unique(MRT_pruned$flyer)
# )
# 
# # --- Get predicted values ---
# newdat$fit <- predict(gls_model, newdata = newdat)
# 
# # --- Plot ---
# ggplot(MRT_pruned, aes(x = log_mass, y = log_MRT_hrs)) +
#   geom_point(size = 3, alpha = 0.8, aes(color=flyer)) +
#   # lines colored by flyer to show vertical shift
#   geom_line(data = newdat,
#             aes(x = log_mass, y = fit, color = flyer),
#             linewidth = 1.2) +
#   scale_color_manual(values = c("darkblue", "#E08929")) + #flyer is redish orange and non-flyer is blue
#   theme_classic(base_size = 14) +
#   labs(
#     x = expression("Log"[10] * " body mass (kg)"),
#     y = "Log Mean Retention Time (hrs)\n",
#     shape = "Diet",
#     color = "Flyer"
#   ) +
#   theme(legend.position = "none")







##### panel b with mass - mean centered mass 
library(nlme)
library(ggplot2)
library(dplyr)

# --- Prepare dataset with mean-centered log_mass ---
MRT_pruned <- MRT_pruned %>%
  mutate(
    log_mass = log10(mass_kg),
    log_mass_c = log_mass - mean(log_mass, na.rm = TRUE)
  )

# --- Fit GLS model ---
gls_model <- gls(
  log_MRT_hrs ~ flyer + log_mass_c,
  data = MRT_pruned,
  correlation = cor_phylo_fixed1,
  method = "ML"
)

# --- Build prediction grid ---
newdat <- expand.grid(
  log_mass_c = seq(min(MRT_pruned$log_mass_c, na.rm = TRUE),
                   max(MRT_pruned$log_mass_c, na.rm = TRUE),
                   length.out = 100),
  flyer = unique(MRT_pruned$flyer)
)

# --- Get predicted values ---
newdat$fit <- predict(gls_model, newdata = newdat)

# --- Plot ---
p2 <-
ggplot(MRT_pruned, aes(x = log_mass_c, y = log_MRT_hrs)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  # Optional: annotate vertical difference between flyer lines
  # annotate("segment",
  #          x = 0, xend = 0,
  #          y = min(newdat$fit[newdat$flyer=="yes"]),
  #          yend = max(newdat$fit[newdat$flyer=="no"]),
  #          color = "black", size = 0.5,
  #          arrow = arrow(length = unit(0.15, "cm"))) +
  geom_point(aes(color=flyer),size = 3, alpha = 0.8) +
  # GLS lines
  # geom_line(data = newdat,
  #           aes(x = log_mass_c, y = fit, color = flyer),
  #           linewidth = 1.2) +
  # Vertical line at mean mass
  scale_color_manual(values = c("darkblue", "#E08929")) +
  theme_classic(base_size = 14) +
  labs(
    x = "\nMean-centered Log10 body mass (kg)",
    y = "Log Mean Retention Time (hrs)\n",
    shape = "Diet",
    color = "Flyer"
  ) + theme(legend.position= "none")


# --- Save the plot ---
ggsave("figures_r4/MRT_plot2b.png", plot = p2, width = 7, height = 6, dpi = 300)




##### cowplot the figures
out.plot2 <- cowplot::plot_grid(p1, p2, nrow=1, ncol=2, labels=c("(A)", "(B)"), rel_widths = c(1.2,1.2))

print(out.plot2)

ggsave(file = paste0(homewd,"/figures_r4/FigS2_TwoPanel_R4.jpeg"),
       units="mm",  
       width=170, 
       height=70, 
       scale=3, 
       dpi=300)
