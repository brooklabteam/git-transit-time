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

length(unique(dat_MRT$genus.species)) #240 unique species for transit time
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

# Phylogenetic signal lambda : 0.820087 
# logL(lambda) : -967.56 
# LR(lambda=0) : 137.925 
# P-value (based on LR test) : 7.56859e-32 



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
# Value Std.Error   t-value      p-value
# (Intercept)  1.106883 0.3658369  3.025619 0.0027005942
# flyer1      -0.979569 0.2897306 -3.380966 0.0008200439

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
# Value Std.Error   t-value      p-value
# (Intercept)  1.2704204 0.2990841  4.247703 2.901354e-05
# flyer1      -0.9643594 0.2523240 -3.821910 1.616340e-04

AIC(gls_model_MRT_0, gls_model_MRT_1, gls_model_MRT_2)
# Emily: The model with full brownian phylogenetic effects is the best fit.
# df        AIC
# gls_model_MRT_0  3   652.0917
# gls_model_MRT_1  3 -2603.3888
# gls_model_MRT_2  3 -2356.8415



summary(gls_model_MRT_1)$tTable
# Value Std.Error   t-value      p-value
# (Intercept)  1.106883 0.3658369  3.025619 0.0027005942
# flyer1      -0.979569 0.2897306 -3.380966 0.0008200439

# Emily: flyer is a significant term which says that flight is significantly
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
# Value  Std.Error    t-value     p-value
# (Intercept)                    1.11354193 0.36249869  3.0718509 0.002330760
# flyer1                        -0.96150947 0.29563591 -3.2523433 0.001280691
# trial.dietfruit/nectar/pollen -0.02454781 0.08267623 -0.2969149 0.766745472
# trial.dietmeat                 0.39162042 0.13203398  2.9660578 0.003269156
# trial.dietmixed                0.04147622 0.34789285  0.1192212 0.905183191
# trial.dietprotein              0.03075949 0.11749517  0.2617937 0.793667639
# trial.dietunk                  0.67686643 0.52385479  1.2920879 0.197362526
# trial.dietunknown             -0.04481134 0.20539463 -0.2181719 0.827449670



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
# Value  Std.Error    t-value      p-value
# (Intercept)                    1.255792696 0.29520449  4.2539756 2.842666e-05
# flyer1                        -0.965816895 0.25754984 -3.7500195 2.137276e-04
# trial.dietfruit/nectar/pollen -0.008484854 0.07601623 -0.1116190 9.112033e-01
# trial.dietmeat                 0.397883091 0.11895486  3.3448243 9.324894e-04
# trial.dietmixed                0.092743820 0.31435772  0.2950264 7.681863e-01
# trial.dietprotein              0.040895331 0.10842718  0.3771686 7.063261e-01
# trial.dietunk                  0.604989346 0.48199768  1.2551707 2.104343e-01
# trial.dietunknown              0.119854681 0.20885876  0.5738552 5.665137e-01

AIC(gls_model_MRT_0, gls_model_MRT_1, gls_model_MRT_2, 
    gls_model_MRT_0_diet, gls_model_MRT_1_diet, gls_model_MRT_2_diet)
# df        AIC
# gls_model_MRT_0       3   652.0917
# gls_model_MRT_1       3 -2603.3888
# gls_model_MRT_2       3 -2356.8415
# gls_model_MRT_0_diet  9   644.9016
# gls_model_MRT_1_diet  9 -2604.0126
# gls_model_MRT_2_diet  9 -2359.1106


summary(gls_model_MRT_1_diet)$tTable
# Value  Std.Error    t-value     p-value
# (Intercept)                    1.11354193 0.36249869  3.0718509 0.002330760
# flyer1                        -0.96150947 0.29563591 -3.2523433 0.001280691
# trial.dietfruit/nectar/pollen -0.02454781 0.08267623 -0.2969149 0.766745472
# trial.dietmeat                 0.39162042 0.13203398  2.9660578 0.003269156
# trial.dietmixed                0.04147622 0.34789285  0.1192212 0.905183191
# trial.dietprotein              0.03075949 0.11749517  0.2617937 0.793667639
# trial.dietunk                  0.67686643 0.52385479  1.2920879 0.197362526
# trial.dietunknown             -0.04481134 0.20539463 -0.2181719 0.827449670

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
# (Intercept)  1.0020987 0.1776598 288  5.640547 4.048481e-08
# flyer1      -0.8692119 0.2400730 288 -3.620615 3.471541e-04

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
# (Intercept)  1.0186786 0.1489817 288  6.837611 4.816690e-11
# flyer1      -0.8699707 0.2076528 288 -4.189544 3.720722e-05

AIC(gls_model_MRT_0, gls_model_MRT_0_diet, gls_model_MRT_0_diet_re,
    gls_model_MRT_1, gls_model_MRT_1_diet, gls_model_MRT_1_diet_re,
    gls_model_MRT_2, gls_model_MRT_2_diet, gls_model_MRT_2_diet_re)

# df        AIC
# gls_model_MRT_0          3   652.0917
# gls_model_MRT_0_diet     9   644.9016
# gls_model_MRT_0_diet_re  4   647.5806
# gls_model_MRT_1          3 -2603.3888
# gls_model_MRT_1_diet     9 -2604.0126 **
# gls_model_MRT_1_diet_re  4 -1529.7376
# gls_model_MRT_2          3 -2356.8415
# gls_model_MRT_2_diet     9 -2359.1106
# gls_model_MRT_2_diet_re  4 -1219.6939

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
# Value  Std.Error   t-value      p-value
# (Intercept)            1.0238383 0.35026088  2.923074 3.737047e-03
# log10(mass_kg)         0.1969964 0.03792113  5.194899 3.850322e-07
# flyer1                -1.7562799 0.61182358 -2.870566 4.396946e-03
# log10(mass_kg):flyer1 -0.6041602 0.24590960 -2.456839 1.459830e-02


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
# Value  Std.Error   t-value      p-value
# (Intercept)            1.2483943 0.29035285  4.299576 2.335257e-05
# log10(mass_kg)         0.1528385 0.03682937  4.149907 4.367707e-05
# flyer1                -1.7715347 0.58042807 -3.052118 2.481383e-03
# log10(mass_kg):flyer1 -0.5486613 0.23180124 -2.366947 1.858772e-02



AIC(gls_model_MRT_0_mass, gls_model_MRT_1_mass, gls_model_MRT_2_mass)
#Emily: model 1 is the best
#Cara: This shows that transit time is still influenced by phylogeny, not
# exclusively by mass (e.g. the model with phylogeny AND mass performed better
# than the model with mass alone)

# df        AIC
# gls_model_MRT_0_mass  5   532.5672
# gls_model_MRT_1_mass  5 -2627.7650 *
# gls_model_MRT_2_mass  5 -2372.4915

summary(gls_model_MRT_1_mass)$tTable

# Value  Std.Error   t-value      p-value
# (Intercept)            1.0238383 0.35026088  2.923074 3.737047e-03
# log10(mass_kg)         0.1969964 0.03792113  5.194899 3.850322e-07
# flyer1                -1.7562799 0.61182358 -2.870566 4.396946e-03
# log10(mass_kg):flyer1 -0.6041602 0.24590960 -2.456839 1.459830e-02

# Emily: flyer is significant (negative) term here (p=0.03)
# We can test whether it would be best to drop the interaction which 
# might clarify the variables that matter most.


gls_model_MRT_1_mass_simple <- gls(log_MRT_hrs ~ flyer + log10(mass_kg), 
                             data = MRT_pruned, 
                             correlation = cor_phylo_fixed1,
                             method = "ML")
summary(gls_model_MRT_1_mass_simple)$tTable
# Value  Std.Error   t-value      p-value
# (Intercept)     1.0330241 0.35323805  2.924442 3.720190e-03
# flyer1         -0.4441066 0.30100972 -1.475390 1.411816e-01
# log10(mass_kg)  0.1802153 0.03762011  4.790399 2.647871e-06


AIC(gls_model_MRT_1_mass, gls_model_MRT_1_mass_simple)

# df       AIC
# gls_model_MRT_1_mass         5 -2627.765
# gls_model_MRT_1_mass_simple  4 -2623.709


# Emily: there is no difference in AIC,really, though the complex model is slightly better

summary(gls_model_MRT_1_mass_simple)$tTable
# Value  Std.Error   t-value      p-value
# (Intercept)     1.0330241 0.35323805  2.924442 3.720190e-03
# flyer1         -0.4441066 0.30100972 -1.475390 1.411816e-01
# log10(mass_kg)  0.1802153 0.03762011  4.790399 2.647871e-06


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
# Value  Std.Error    t-value      p-value
# (Intercept)                    1.05610780 0.35163515  3.0034193 0.0029049121
# flyer1                        -0.50953955 0.30431780 -1.6743666 0.0951479084
# log10(mass_kg)                 0.17171835 0.03889695  4.4146994 0.0000143375
# trial.dietfruit/nectar/pollen  0.00972326 0.08051873  0.1207577 0.9039674216
# trial.dietmeat                 0.34545688 0.12841588  2.6901414 0.0075605840
# trial.dietmixed                0.12080555 0.33771438  0.3577151 0.7208194225
# trial.dietprotein              0.12976388 0.11608284  1.1178559 0.2645633825
# trial.dietunk                  0.58132795 0.50826866  1.1437415 0.2536839293
# trial.dietunknown              0.06074774 0.20053347  0.3029307 0.7621621627

# then, as a random effect
gls_model_MRT_1_mass_simple_diet_re <- lme(log_MRT_hrs ~ flyer + log10(mass_kg),
                                            random = ~1 | trial.diet,
                                            data = MRT_pruned, 
                                            correlation = cor_phylo_fixed1,
                                            method = "ML")
summary(gls_model_MRT_1_mass_simple_diet_re)$tTable
# Value  Std.Error  DF   t-value      p-value
# (Intercept)     0.9985459 0.17498721 287  5.706394 2.873297e-08
# flyer1         -0.6130110 0.24987085 287 -2.453311 1.474951e-02
# log10(mass_kg)  0.1005826 0.03170947 287  3.172006 1.678036e-03

# then, compare
AIC(gls_model_MRT_1_mass_simple, gls_model_MRT_1_mass_simple_diet, gls_model_MRT_1_mass_simple_diet_re)

# df       AIC
# gls_model_MRT_1_mass_simple          4 -2623.709
# gls_model_MRT_1_mass_simple_diet    10 -2621.460
# gls_model_MRT_1_mass_simple_diet_re  5 -1537.732

# here, Adding diet does nothing really for the fit. So only have mass in the model is still the best

summary(gls_model_MRT_1_mass_simple)$tTable

# Value  Std.Error   t-value      p-value
# (Intercept)     1.0330241 0.35323805  2.924442 3.720190e-03
# flyer1         -0.4441066 0.30100972 -1.475390 1.411816e-01
# log10(mass_kg)  0.1802153 0.03762011  4.790399 2.647871e-06

# Emily: flyer is not significantly negatively related to transit, 
# after accounting for mass (and phylogeny). 


# Cara: what if we compare our two best fit models so far? 
AIC(gls_model_MRT_1, gls_model_MRT_1_diet, gls_model_MRT_1_mass_simple)

# df       AIC
# gls_model_MRT_1              3 -2603.389
# gls_model_MRT_1_diet         9 -2604.013
# gls_model_MRT_1_mass_simple  4 -2623.709


gls_model_MRT_1_mass_only <- gls(log_MRT_hrs ~ log10(mass_kg), 
                                   data = MRT_pruned, 
                                   correlation = cor_phylo_fixed1,
                                   method = "ML")
summary(gls_model_MRT_1_mass_only)$tTable
# Value  Std.Error  t-value      p-value
# (Intercept)    1.0118063 0.35365084 2.861032 4.525474e-03
# log10(mass_kg) 0.2008265 0.03499993 5.737913 2.382712e-08

gls_model_MRT_0_mass_only <- gls(log_MRT_hrs ~ log10(mass_kg), 
                                 data = MRT_pruned, 
                                 #correlation = cor_phylo_fixed0,
                                 method = "ML")
summary(gls_model_MRT_0_mass_only)$tTable
# Value  Std.Error  t-value      p-value
# (Intercept)    1.0065686 0.03450639 29.17050 8.680808e-89
# log10(mass_kg) 0.2783932 0.02247669 12.38586 1.231486e-28

gls_model_MRT_2_mass_only <- gls(log_MRT_hrs ~ log10(mass_kg), 
                                 data = MRT_pruned, 
                                 correlation = cor_phylo_fixed2,
                                 method = "ML")
summary(gls_model_MRT_2_mass_only)$tTable
# Value  Std.Error  t-value      p-value
# (Intercept)    1.232982 0.29405351 4.193054 3.646249e-05
# log10(mass_kg) 0.168163 0.03362204 5.001570 9.796621e-07
# mass is always signficant

AIC(gls_model_MRT_0_mass_only, gls_model_MRT_1_mass_only, gls_model_MRT_2_mass_only)
# 
# df        AIC
# gls_model_MRT_0_mass_only  3   535.0061
# gls_model_MRT_1_mass_only  3 -2623.5177
# gls_model_MRT_2_mass_only  3 -2366.6604








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
