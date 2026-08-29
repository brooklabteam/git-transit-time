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
dat_MRT <- read.csv(file = paste0(homewd, "/data/revision-3/dat_sum_tot_clean_MRT_R3_v2.csv"), header = T, stringsAsFactors = F )

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
tree <- read.tree(file = paste0(homewd, "/data/revision-3/Book3.nwk"))
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
length(common_species) # CARA: I get 230 here, not 222, EMILY: SAME
tree <- drop.tip(tree_ultra2, setdiff(tree_ultra2$tip.label, common_species)) 
MRT_hrs <- MRT_hrs[common_species] #223
length(MRT_hrs) #CARA: I get 230 here as well
#for MRT
lambda_gs2<-phylosig(tree, MRT_hrs,
                     method="lambda",test=TRUE)
lambda_gs2

#EMILY GETS:
# SAME NOW


#CARA gets:
# Phylogenetic signal lambda : 0.77132 
# logL(lambda) : -956.715 
# LR(lambda=0) : 140.083 
# P-value (based on LR test) : 2.55322e-32 

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

tree <- read.tree(paste0(homewd, "/data/revision-3/Book3.nwk"))
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
# [1] "Flying Birds"        "Carnivora"           "Even-toed Ungulates" "Primates"            "Non-Flying Birds"   
# [6] "Rodents"             "Odd-toed Ungulates"  "Bats" 

#CARA I get this (includes non-flying birds):
# [1] "Rodents"             "Bats"                "Flying Birds"        "Primates"            "Even-toed Ungulates"
# [6] "Carnivora"           "Non-Flying Birds"    "Odd-toed Ungulates" 

length(unique(MRT_pruned$phylo_name)) #230. CARA: I agree


# MRT_pruned$re_class <- factor(MRT_pruned$re_class, levels = c("Rodents","Bats", "Non-Flying Birds", 
#                                                               "Primates", "Ungulates"))

#transform data
MRT_pruned$log_MRT_hrs <- log10(MRT_pruned$MRT_hrs)


# add the "flyer" parameter
MRT_pruned$flyer <- 0
MRT_pruned$flyer[MRT_pruned$re_class=="Bats"] <- 1
#CARA: I added the next line below:
MRT_pruned$flyer[MRT_pruned$re_class=="Flying Birds"] <- 1

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

#CARA: I get this
# Value Std.Error   t-value      p-value
# (Intercept)  1.3267590 0.3702254  3.583652 0.0003964841
# flyer1      -0.6448687 0.1839075 -3.506483 0.0005249530


#with no phylogenetic effects:
gls_model_MRT_0 <- gls(log_MRT_hrs ~ flyer, 
                        data = MRT_pruned, 
                        #correlation = cor_phylo_fixed0,
                        method = "ML")
summary(gls_model_MRT_0)$tTable


# Value  Std.Error   t-value       p-value
# (Intercept)  1.3192035 0.04013526  32.86894 1.841621e-100
# flyer1      -0.9861522 0.06975235 -14.13791  5.576318e-35

#CARA: I get the same

#with estimated phylogenetic effects:
gls_model_MRT_2 <- gls(log_MRT_hrs ~ flyer,
                        data = MRT_pruned,
                        correlation = cor_phylo_fixed2,
                        method = "ML")
summary(gls_model_MRT_2)$tTable

#CARA: I get this:
# Value Std.Error   t-value      p-value
# (Intercept)  1.4813660 0.2984518  4.963501 1.174555e-06
# flyer1      -0.5947408 0.1618582 -3.674455 2.831908e-04


AIC(gls_model_MRT_0, gls_model_MRT_1, gls_model_MRT_2)
# Emily: The model with full brownian phylogenetic effects is the best fit.
#CARA: I agree with conclusion, though values are different

#CARA: I get this
# df        AIC
# gls_model_MRT_0  3   505.7558
# gls_model_MRT_1  3 -3346.0869
# gls_model_MRT_2  3 -2888.3107

#EMILY:
# df        AIC
# gls_model_MRT_0  3   505.7558
# gls_model_MRT_1  3 -2607.2022
# gls_model_MRT_2  3 -2211.8379


summary(gls_model_MRT_1)$tTable

#CARA: I get this, EMILY: same
# Value Std.Error   t-value      p-value
# (Intercept)  1.3267590 0.3702254  3.583652 0.0003964841
# flyer1      -0.6448687 0.1839075 -3.506483 0.0005249530

# Emily: flyer is a significant term which says that flight is significantly
# negatively related to MRT transit time.
# This is on top of a model that already includes significant 
# phylogenetic clustering in transit time, suggesting that flight affects transit
# time independent of phylogeny.

#CARA: I agree with conclusion


# now, let's also account for diet
# first try as a fixed effect:
MRT_pruned$trial.diet <- as.factor(MRT_pruned$trial.diet)

# Cara: with full phylogenetic effects:
gls_model_MRT_1_diet <- gls(log_MRT_hrs ~ flyer + trial.diet, 
                        data = MRT_pruned, 
                        correlation = cor_phylo_fixed1,
                        method = "ML")
summary(gls_model_MRT_1_diet)$tTable

#CARA: I get this. EMILY: same
# Value  Std.Error     t-value      p-value
# (Intercept)                    1.339701142 0.36806143  3.63988460 0.0003230573
# flyer1                        -0.635685290 0.18503128 -3.43555586 0.0006780549
# trial.dietfruit/nectar/pollen -0.016624469 0.08194591 -0.20287126 0.8393784143
# trial.dietmeat                 0.376447653 0.13121206  2.86900193 0.0044213120
# trial.dietmixed                0.078159067 0.36760384  0.21261766 0.8317750876
# trial.dietprotein              0.009424897 0.11618347  0.08112081 0.9354020276
# trial.dietunknown              0.074813791 0.19385207  0.38593238 0.6998306927

# With no phylogenetic effects:
gls_model_MRT_0_diet <- gls(log_MRT_hrs ~ flyer + trial.diet, 
                        data = MRT_pruned, 
                        #correlation = cor_phylo_fixed0,
                        method = "ML")
summary(gls_model_MRT_0_diet)$tTable
# Value  Std.Error     t-value      p-value
# (Intercept)                    1.32783581 0.04527590  29.3276532 1.291034e-88
# flyer1                        -1.00240624 0.07265751 -13.7963196 1.307279e-33
# trial.dietfruit/nectar/pollen -0.18404776 0.08536267  -2.1560685 3.190342e-02
# trial.dietmeat                 0.27694179 0.10563402   2.6217102 9.212383e-03
# trial.dietmixed               -0.27091314 0.39210477  -0.6909203 4.901702e-01
# trial.dietprotein             -0.06414133 0.11381402  -0.5635627 5.734888e-01
# trial.dietunknown              0.46571397 0.25003447   1.8625990 6.353327e-02

#CARA: I get the same



# With estimated phylogenetic effects:
gls_model_MRT_2_diet <- gls(log_MRT_hrs ~ flyer + trial.diet,
                        data = MRT_pruned,
                        correlation = cor_phylo_fixed2,
                        method = "ML")
summary(gls_model_MRT_2_diet)$tTable

#CARA: I get this: EMILY, same:
# Value  Std.Error     t-value      p-value
# (Intercept)                    1.460816049 0.29522260  4.94818498 1.273993e-06
# flyer1                        -0.591898721 0.16268143 -3.63839141 3.248527e-04
# trial.dietfruit/nectar/pollen -0.003657727 0.07659246 -0.04775571 9.619439e-01
# trial.dietmeat                 0.373756371 0.11743144  3.18276235 1.617853e-03
# trial.dietmixed                0.257646320 0.34784395  0.74069513 4.594797e-01
# trial.dietprotein              0.027303431 0.10784541  0.25317194 8.003153e-01
# trial.dietunknown              0.241778753 0.19817459  1.22002905 2.234485e-01

AIC(gls_model_MRT_0, gls_model_MRT_1, gls_model_MRT_2, 
    gls_model_MRT_0_diet, gls_model_MRT_1_diet, gls_model_MRT_2_diet)

#CARA: I get this
# df        AIC
# gls_model_MRT_0       3   505.7558
# gls_model_MRT_1       3 -3346.0869
# gls_model_MRT_2       3 -2888.3107
# gls_model_MRT_0_diet  8   495.8800
# gls_model_MRT_1_diet  8 -3345.9899
# gls_model_MRT_2_diet  8 -2890.6239

#EMILY
# df        AIC
# gls_model_MRT_0       3   505.7558
# gls_model_MRT_1       3 -2607.2022
# gls_model_MRT_2       3 -2211.8379
# gls_model_MRT_0_diet  8   495.8800
# gls_model_MRT_1_diet  8 -2607.1052
# gls_model_MRT_2_diet  8 -2214.1510

summary(gls_model_MRT_1_diet)$tTable


#CARA: I get this, EMILY SAME
# Value  Std.Error     t-value      p-value
# (Intercept)                    1.339701142 0.36806143  3.63988460 0.0003230573
# flyer1                        -0.635685290 0.18503128 -3.43555586 0.0006780549
# trial.dietfruit/nectar/pollen -0.016624469 0.08194591 -0.20287126 0.8393784143
# trial.dietmeat                 0.376447653 0.13121206  2.86900193 0.0044213120
# trial.dietmixed                0.078159067 0.36760384  0.21261766 0.8317750876
# trial.dietprotein              0.009424897 0.11618347  0.08112081 0.9354020276
# trial.dietunknown              0.074813791 0.19385207  0.38593238 0.6998306927

# Emily: model gls_model_MRT_1_diet has a close AIC to the model without diet
# Plus only meat comes out as significant when you do the summary output. 
# So prefer the simpler model?
# CARA: I think we can say this. typically only AIC changes of 4 or more are 
# considered worth the additional complexity. here, only 1. in most recent version,
# simpler model is 1 AIC point lower anyhow in my version



# Cara: let's try it also as a random effect - here we use the lme function which
# is also supported in the nlme package. we can compare these since all were fit via ML:

gls_model_MRT_1_diet_re <- lme(log_MRT_hrs ~ flyer,
                            random = ~1 | trial.diet,
                            correlation = cor_phylo_fixed1,
                            data = MRT_pruned,
                            method = "ML")
summary(gls_model_MRT_1_diet_re)$tTable
# Value Std.Error  DF   t-value      p-value
# (Intercept)  1.1716953 0.1924782 289  6.087419 3.641363e-09
# flyer1      -0.6495386 0.1541441 289 -4.213841 3.359595e-05

#CARA: I get the same


gls_model_MRT_0_diet_re <- lme(log_MRT_hrs ~ flyer,
                                random = ~1 | trial.diet,
                                #correlation = cor_phylo_fixed0,
                                data = MRT_pruned,
                                method = "ML")
summary(gls_model_MRT_0_diet_re)$tTable
# Value  Std.Error  DF   t-value      p-value
# (Intercept)  1.3548978 0.08652659 289  15.65875 1.929431e-40
# flyer1      -0.9983666 0.07189340 289 -13.88676 6.134079e-34

#CARA: I get the same


gls_model_MRT_2_diet_re <- lme(log_MRT_hrs ~ flyer,
                                random = ~1 | trial.diet,
                                correlation = cor_phylo_fixed2,
                                data = MRT_pruned,
                                method = "ML")
summary(gls_model_MRT_2_diet_re)$tTable
# Value Std.Error  DF   t-value      p-value
# (Intercept)  1.2124287 0.1608278 289  7.538674 6.162535e-13
# flyer1      -0.6536117 0.1356990 289 -4.816628 2.359685e-06

#Cara: I get the same

AIC(gls_model_MRT_0, gls_model_MRT_0_diet, gls_model_MRT_0_diet_re,
    gls_model_MRT_1, gls_model_MRT_1_diet, gls_model_MRT_1_diet_re,
    gls_model_MRT_2, gls_model_MRT_2_diet, gls_model_MRT_2_diet_re)
#EMILY
# df        AIC
# gls_model_MRT_0          3   505.7558
# gls_model_MRT_0_diet     8   495.8800
# gls_model_MRT_0_diet_re  4   501.8596
# gls_model_MRT_1          3 -2607.2022
# gls_model_MRT_1_diet     8 -2607.1052
# gls_model_MRT_1_diet_re  4 -1536.0665
# gls_model_MRT_2          3 -2211.8379
# gls_model_MRT_2_diet     8 -2214.1510
# gls_model_MRT_2_diet_re  4 -1196.6698

#CARA: I get this
# df        AIC
# gls_model_MRT_0          3   505.7558
# gls_model_MRT_0_diet     8   495.8800
# gls_model_MRT_0_diet_re  4   501.8596
# gls_model_MRT_1          3 -3346.0869 *
# gls_model_MRT_1_diet     8 -3345.9899 *
# gls_model_MRT_1_diet_re  4 -1557.1889
# gls_model_MRT_2          3 -2888.3107
# gls_model_MRT_2_diet     8 -2890.6239
# gls_model_MRT_2_diet_re  4 -1226.5116

#Emily: in all cases, the fit with diet as a fixed effect was better supported than the fit as
# a random effect. 
# Cara: The models with the lambda=1, which incorporate phylogenetic clustering
# of MRT transit time, are best supported. The top fit model, gls_model_MRT_1,
# offers the best fit overall.

#CARA: I agree with conclusion, also, diet prob not necessary to include here



###################
########## model gls_model_MRT_1 is best for reporting.




# now, let's also consider the effects of mass and, later, mass and diet
 
#with full phylogenetic effects:
gls_model_MRT_1_mass <- gls(log_MRT_hrs ~ log10(mass_kg)*flyer, 
                        data = MRT_pruned, 
                        correlation = cor_phylo_fixed1,
                        method = "ML")
summary(gls_model_MRT_1_mass)$tTable

#CARA: I get this
# Value  Std.Error    t-value     p-value
# (Intercept)            1.1797290 0.36129339  3.2652937 0.001223687
# log10(mass_kg)         0.1291901 0.04508749  2.8653212 0.004468380
# flyer1                -0.1506248 0.21140315 -0.7125003 0.476723961
# log10(mass_kg):flyer1  0.1484832 0.08765610  1.6939289 0.091345117

#with no phylogenetic effects:
gls_model_MRT_0_mass <- gls(log_MRT_hrs ~ log10(mass_kg)*flyer, 
                        data = MRT_pruned, 
                        #correlation = cor_phylo_fixed0,
                        method = "ML")
summary(gls_model_MRT_0_mass)$tTable
# Value  Std.Error   t-value      p-value
# (Intercept)            1.25802238 0.04420326 28.459945 3.344981e-86
# log10(mass_kg)         0.08852655 0.03201528  2.765135 6.052202e-03
# flyer1                -0.48653978 0.11402962 -4.266784 2.682666e-05
# log10(mass_kg):flyer1  0.19493542 0.06638343  2.936507 3.583426e-03

#CARA: I get the same


#with estimated phylogenetic effects:
gls_model_MRT_2_mass <- gls(log_MRT_hrs ~ log10(mass_kg)*flyer,
                        data = MRT_pruned,
                        correlation = cor_phylo_fixed2,
                        method = "ML")
summary(gls_model_MRT_2_mass)$tTable

#CARA: I get this
# Value  Std.Error    t-value      p-value
# (Intercept)            1.37138106 0.29543633  4.6418836 5.222790e-06
# log10(mass_kg)         0.11114006 0.04403690  2.5237942 1.214005e-02
# flyer1                -0.26061231 0.19768522 -1.3183196 1.884298e-01
# log10(mass_kg):flyer1  0.06207325 0.08605438  0.7213259 4.712862e-01



AIC(gls_model_MRT_0_mass, gls_model_MRT_1_mass, gls_model_MRT_2_mass)
#Emily: model 1 is the best, by a lot
#Cara: This shows that transit time is still influenced by phylogeny, not
# exclusively by mass (e.g. the model with phylogeny AND mass performed better
# than the model with mass alone)

#CARA: I agree with conclusion

#EMILY
# df        AIC
# gls_model_MRT_0_mass  5   479.5181
# gls_model_MRT_1_mass  5 -2625.3342
# gls_model_MRT_2_mass  5 -2219.8868

#CARA: I get this
# df        AIC
# gls_model_MRT_0_mass  5   479.5181
# gls_model_MRT_1_mass  5 -3364.2189
# gls_model_MRT_2_mass  5 -2896.3597

summary(gls_model_MRT_1_mass)$tTable

#CARA: I get this
# Value  Std.Error    t-value     p-value
# (Intercept)            1.1797290 0.36129339  3.2652937 0.001223687
# log10(mass_kg)         0.1291901 0.04508749  2.8653212 0.004468380
# flyer1                -0.1506248 0.21140315 -0.7125003 0.476723961
# log10(mass_kg):flyer1  0.1484832 0.08765610  1.6939289 0.091345117

# Emily: flyer is not significant term here (p=0.41)
# We can test whether it would be best to drop the interaction which 
# might clarify the variables that matter most.
#CARA: flyer is borderline NOT sig in mine either.
#however, the interaction of mass and flyer is significant and positive
#suggesting that bigger flyers have longer GIT transit. 

gls_model_MRT_1_mass_simple <- gls(log_MRT_hrs ~ flyer + log10(mass_kg), 
                             data = MRT_pruned, 
                             correlation = cor_phylo_fixed1,
                             method = "ML")
summary(gls_model_MRT_1_mass_simple)$tTable

#CARA: I get this:
# Value Std.Error   t-value      p-value
# (Intercept)     1.1366997 0.3615472  3.143987 0.0018376244
# flyer1         -0.2893619 0.1955137 -1.480009 0.1399453774
# log10(mass_kg)  0.1698441 0.0382908  4.435637 0.0000130024

AIC(gls_model_MRT_1_mass, gls_model_MRT_1_mass_simple)
#EMILY
# df       AIC
# gls_model_MRT_1_mass         5 -2625.334
# gls_model_MRT_1_mass_simple  4 -2624.440

#CARA: I get this
# df       AIC
# gls_model_MRT_1_mass         4 -3363.324
# gls_model_MRT_1_mass_simple  4 -3363.324

# Emily: there is no difference in AIC
#CARA: I agree with conclusion, no difference at all

summary(gls_model_MRT_1_mass_simple)$tTable

#CARA: I get this
# Value Std.Error   t-value      p-value
# (Intercept)     1.1366997 0.3615472  3.143987 0.0018376244
# flyer1         -0.2893619 0.1955137 -1.480009 0.1399453774
# log10(mass_kg)  0.1698441 0.0382908  4.435637 0.0000130024

# Emily: note that flight is significant in the complex model.
#CARA: No, I don't think it is...

# now, let's see if adding diet improves this model:


# I'll just illustrate with the best fit model
# first, as a fixed effect
gls_model_MRT_1_mass_simple_diet <- gls(log_MRT_hrs ~ flyer + log10(mass_kg) +
                                          trial.diet, 
                                    data = MRT_pruned, 
                                    correlation = cor_phylo_fixed1,
                                    method = "ML")
summary(gls_model_MRT_1_mass_simple_diet)$tTable

#CARA: I get this
# Value  Std.Error    t-value      p-value
# (Intercept)                    1.172184268 0.35997225  3.2563184 1.263528e-03
# flyer1                        -0.314655089 0.19525832 -1.6114811 1.081702e-01
# log10(mass_kg)                 0.166595345 0.03943597  4.2244520 3.217123e-05
# trial.dietfruit/nectar/pollen  0.009095385 0.07988949  0.1138496 9.094364e-01
# trial.dietmeat                 0.340314443 0.12783378  2.6621638 8.200297e-03
# trial.dietmixed                0.133017954 0.35757254  0.3720027 7.101644e-01
# trial.dietprotein              0.114579878 0.11564907  0.9907549 3.226369e-01
# trial.dietunknown              0.142084836 0.18910942  0.7513366 4.530636e-01


# then, as a random effect
gls_model_MRT_1_mass_simple_diet_re <- lme(log_MRT_hrs ~ flyer + log10(mass_kg),
                                            random = ~1 | trial.diet,
                                            data = MRT_pruned, 
                                            correlation = cor_phylo_fixed1,
                                            method = "ML")
summary(gls_model_MRT_1_mass_simple_diet_re)$tTable
# Value  Std.Error  DF   t-value      p-value
# (Intercept)     1.11373914 0.19125249 288  5.823397 1.538261e-08
# flyer1         -0.46896803 0.16488611 288 -2.844194 4.771407e-03
# log10(mass_kg)  0.09119103 0.03190991 288  2.857765 4.577483e-03

#CARA: I get the same

# then, compare
AIC(gls_model_MRT_1_mass_simple, gls_model_MRT_1_mass_simple_diet, gls_model_MRT_1_mass_simple_diet_re)
#EMILY
# df       AIC
# gls_model_MRT_1_mass_simple          4 -2624.440
# gls_model_MRT_1_mass_simple_diet     9 -2622.901
# gls_model_MRT_1_mass_simple_diet_re  5 -1542.204

#CARA: I get this
# df       AIC
# gls_model_MRT_1_mass_simple          4 -3363.324
# gls_model_MRT_1_mass_simple_diet     9 -3361.786
# gls_model_MRT_1_mass_simple_diet_re  5 -1563.326

# here, Adding diet does nothing really for the fit. So only have mass in the model is still the best
#CARA: I agree with conclusion

summary(gls_model_MRT_1_mass_simple)$tTable

#CARA: I get this
# Value Std.Error   t-value      p-value
# (Intercept)     1.1366997 0.3615472  3.143987 0.0018376244
# flyer1         -0.2893619 0.1955137 -1.480009 0.1399453774
# log10(mass_kg)  0.1698441 0.0382908  4.435637 0.0000130024


# Emily: flyer is not significantly negatively related to transit, in the simple model
# after accounting for mass (and phylogeny). 

#CARA: I agree with above statement


# Cara: what if we compare our two best fit models so far? 
AIC(gls_model_MRT_1, gls_model_MRT_1_diet, gls_model_MRT_1_mass_simple)
#EMILY
# df       AIC
# gls_model_MRT_1              3 -2607.202
# gls_model_MRT_1_diet         8 -2607.105
# gls_model_MRT_1_mass_simple  4 -2624.440

#CARA: I get this
# df       AIC
# gls_model_MRT_1              3 -3346.087
# gls_model_MRT_1_diet         8 -3345.990
# gls_model_MRT_1_mass_simple  4 -3363.324


gls_model_MRT_1_mass_only <- gls(log_MRT_hrs ~ log10(mass_kg), 
                                   data = MRT_pruned, 
                                   correlation = cor_phylo_fixed1,
                                   method = "ML")
summary(gls_model_MRT_1_mass_only)$tTable

#CARA: I get this, EMILY SAME
# Value  Std.Error  t-value      p-value
# (Intercept)    1.0178135 0.35322361 2.881499 4.249339e-03
# log10(mass_kg) 0.1930754 0.03499625 5.517032 7.565699e-08

gls_model_MRT_0_mass_only <- gls(log_MRT_hrs ~ log10(mass_kg), 
                                 data = MRT_pruned, 
                                 #correlation = cor_phylo_fixed0,
                                 method = "ML")
summary(gls_model_MRT_0_mass_only)$tTable
# Value  Std.Error  t-value      p-value
# (Intercept)    1.0065686 0.03450639 29.17050 8.680808e-89
# log10(mass_kg) 0.2783932 0.02247669 12.38586 1.231486e-28

#CARA: I get the same

# mass hugely significant

gls_model_MRT_2_mass_only <- gls(log_MRT_hrs ~ log10(mass_kg), 
                                 data = MRT_pruned, 
                                 correlation = cor_phylo_fixed2,
                                 method = "ML")
summary(gls_model_MRT_2_mass_only)$tTable

#CARA: I get this
# Value  Std.Error  t-value      p-value
# (Intercept)    1.2321946 0.28675963 4.296960 2.356533e-05
# log10(mass_kg) 0.1582657 0.03361975 4.707523 3.867915e-06

AIC(gls_model_MRT_0_mass_only, gls_model_MRT_1_mass_only, gls_model_MRT_2_mass_only)

#EMILY
# df        AIC
# gls_model_MRT_0_mass_only  3   535.0061
# gls_model_MRT_1_mass_only  3 -2624.2350
# gls_model_MRT_2_mass_only  3 -2220.0579

#CARA: I get this
# df        AIC
# gls_model_MRT_0_mass_only  3   535.0061
# gls_model_MRT_1_mass_only  3 -3363.1197
# gls_model_MRT_2_mass_only  3 -2896.5308




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


ggsave(file = paste0(homewd,"/figures/revision-3/Fig_S1A.png"),
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
ggsave("figures/revision-3/MRT_plot2b.png", plot = p2, width = 7, height = 6, dpi = 300)




##### cowplot the figures
out.plot2 <- cowplot::plot_grid(p1, p2, nrow=1, ncol=2, labels=c("(A)", "(B)"), rel_widths = c(1.2,1.2))

print(out.plot2)

ggsave(file = paste0(homewd,"/figures/revision-3/FigS2_TwoPanel_R4.jpeg"),
       units="mm",  
       width=170, 
       height=70, 
       scale=3, 
       dpi=300)

