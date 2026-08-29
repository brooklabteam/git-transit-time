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

#load the GIT transit data with phylo name
dat_GIT <- read.csv(file = paste0(homewd, "/data/revision-3/dat_sum_tot_clean_GIT_R3_v2.csv"), header = T, stringsAsFactors = F )
length(unique(dat_GIT$phylo_name)) #249
dat_MRT <- read.csv(file = paste0(homewd, "/data/revision-3/dat_sum_tot_clean_MRT_R3_v2.csv"), header = T, stringsAsFactors = F )
length(unique(dat_MRT$phylo_name)) #233

#remove the two species that don't have masses earlier on - While they have data on diet, it will be
# easier to explain if we just remove them from the offset.
# --- names to match your columns ---
mass_col    <- "avg_mass"        # body-mass column
species_col <- "genus.species"  # tip / phylogeny-name column

# species missing a mass in either dataset
miss <- function(d) d[[species_col]][is.na(d[[mass_col]]) | trimws(d[[mass_col]]) == "0"]
no_mass <- sort(unique(c(miss(dat_GIT), miss(dat_MRT))))

# report which ones
cat("Species without mass (", length(no_mass), "):\n", sep = "")
print(no_mass)
#[1] "Atheris nitscheri"    "Uromacer oxyrhynchus"

# flightless = c("Struthio camelus", "Rhea americana", "Dromaius novaehollandiae",
#                 "Casuaris casuaris", "Pygoscelis papua")

# remove them from both datasets
dat_GIT <- dat_GIT[!dat_GIT[[species_col]] %in% no_mass, ]
dat_MRT <- dat_MRT[!dat_MRT[[species_col]] %in% no_mass, ]

#### species counts ####
#table 1 in the manuscript

length(unique(dat_GIT$genus.species)) #257 unique species for transit time
tt.pivot <-ddply(dat_GIT, .(re_class), summarise, N_species_transit = length(unique(genus.species)))
print(tt.pivot)
# re_class N_species_transit
# 1                Bats                37
# 2          Carnivores                26
# 3 Even-toed Ungulates                17
# 4        Flying Birds                73
# 5    Non-Flying Birds                 5
# 6  Odd-toed Ungulates                 7
# 7            Primates                60
# 8            Reptiles                17
# 9             Rodents                15

# dat_GIT <- subset(
#   dat_GIT,
#   re_class != "Non-Flying Birds"
# )


length(unique(dat_MRT$genus.species)) 
(mrt.pivot <-ddply(dat_MRT, .(re_class), summarise, N_species_mrt = length(unique(genus.species))))
# re_class N_species_mrt
# 1                Bats            13
# 2           Carnivora            18
# 3 Even-toed Ungulates            55
# 4        Flying Birds            72
# 5    Non-Flying Birds             7
# 6  Odd-toed Ungulates            10
# 7            Primates            40
# 8             Rodents            27


#load the tree - pulled from timetree so should be in principle an ultrametric tree
tree <- read.tree(file = paste0(homewd, "/data/revision-3/Book3.nwk"))
# is.ultrametric(tree) #check for ultrametric - is FALSE

# "false convergence" happens when the optimization routine thinks it has converged, but really just got stuck in a flat likelihood surface.
# Usually this means the fit didn’t fully optimize branch lengths, but you still get an ultrametric tree out.
# 
# For comparative methods (like gls() or pgls): Usually yes, it’s fine. You just need an ultrametric tree, and the exact timing doesn’t matter much if your question isn’t about divergence dating.

#try to make the warnings less by fitting lambda

# tree_ultra <- chronos(tree, lambda=0.38)
# is.ultrametric(tree_ultra) #still give false convergence

# does this way give less warnings?
tree_ultra2 <- force.ultrametric(tree, method = "nnls") #no warnings

is.ultrametric(tree_ultra2) # TRUE

### LET'S GO WITH THIS force.ultrametric APPROACH AS IT DOES NOT HAVE WARNINGS






##### plot the phylogeny for GIT and MRT #######

# Plot the tree with enhancements
# plot(tree_ultra2,
#      type = "phylogram",    # or "fan", "cladogram", "unrooted"
#      cex = 0.6,             # Tip label size
#      no.margin = TRUE,      # Remove extra white space
#      edge.width = 1.5,      # Thicker branches
#      edge.color = "darkgray", # Branch color
#      label.offset = 0.001,  # Space between tip and label
#      font = 3,              # Italic tip labels
#      main = "")
# 
# # Add a scale bar with adjusted size and location
# add.scale.bar(length = 0.05, lwd = 2, cex = 0.7, col = "black")

# 
# 
# 




##################  GIT phylogenetic signal ##########
dat_GIT$phylo_name <- as.character(dat_GIT$phylo_name)

#for GIT
GIT <- dat_GIT[c("phylo_name", "transit_hrs")] #275. CARA: I get a different number (see below)
nrow(GIT) #296
#GIT <- na.omit(GIT)

length(unique(GIT$phylo_name)) #247. CARA: I get the same here
#three don't have matches in the tree

#make sure they match
# If you have a data frame GIT with the trait
transit_hrs <- setNames(GIT$transit_hrs, GIT$phylo_name) 

# Then only keep the overlapping species
common_species <- intersect(tree_ultra2$tip.label, names(transit_hrs)) #243
length(common_species)#243

########
#### What is missing and not matching?
# in your data but NOT on the tree
in_data_not_tree <- setdiff(names(transit_hrs), tree_ultra2$tip.label)

# on the tree but NOT in your data
in_tree_not_data <- setdiff(tree_ultra2$tip.label, names(transit_hrs))

cat("In data but not on tree (", length(in_data_not_tree), "):\n", sep = ""); print(sort(in_data_not_tree))

# In data but not on tree (4):
#   [1] ""                "Bycanistes"      "Cercopithecinae" "Pycnonotidae" 
#cat("\nOn tree but not in data (", length(in_tree_not_data), "):\n", sep = ""); print(sort(in_tree_not_data))






tree <- drop.tip(tree_ultra2, setdiff(tree_ultra2$tip.label, common_species))
#tree <- drop.tip(tree_ultra, setdiff(tree$tip.label, names(transit_hrs)))
transit_hrs <- transit_hrs[common_species]
setdiff(names(transit_hrs), tree$tip.label)#none
setdiff(tree$tip.label, names(transit_hrs))#none

lambda_gs1<-phylosig(tree, transit_hrs,
                    method="lambda",test=TRUE)
lambda_gs1 

# Emily: lambda is now so close to one that we will use lambda = 1 in the model, and not need model 2
# CARA 9=(on 8/28): agree!

# Cara: a little research suggests that values >1 are pretty rare and since
# this is basically 1, that seems fine to me

#CARA... I get this: EMILY, same
# Phylogenetic signal lambda : 1.00145 
# logL(lambda) : -1552.62 
# LR(lambda=0) : 152.375 
# P-value (based on LR test) : 5.24537e-35 

#tree$tip.label

contMap(tree,log10(transit_hrs), fsize =0.5)
# save as MGT_contmap
# 


#########################################  MGT. #################################

########### MGT data cleaning #######

#_________________________________________________________________________________________________________________
#I think this code runs PGLS and also provides lambda estimates as part of the model


# if you haven't continued from above
# Load your data
# dat_GIT <- read.csv(file = paste0(homewd, "/data/dat_GIT_foranalysis.csv"), header = T, stringsAsFactors = F )
# # add dat_MRT

# clean dataset with only the columns you need
GIT <- dat_GIT[c("phylo_name", "transit_hrs", "trial.diet", "re_class", "mass_kg")]
#GIT <- na.omit(GIT$transit_hrs)

tree <- read.tree(paste0(homewd, "/data/revision-3/Book3.nwk"))
tree_ultra2 <- force.ultrametric(tree, method = "nnls") #no warnings
is.ultrametric(tree_ultra2)

#make sure they match
# Check if the species names match the tree tips
common_species <- intersect(tree_ultra2$tip.label, GIT$phylo_name)
tree_pruned <- drop.tip(tree_ultra2, setdiff(tree_ultra2$tip.label, common_species))
GIT_pruned <- GIT[GIT$phylo_name %in% common_species, ]

# there are some duplicated taxa with different diets:
GIT_pruned[duplicated(GIT_pruned$phylo_name),]


species1 <- unique(dat_GIT$phylo_name)
species2 <- unique(GIT_pruned$phylo_name)

# Species in df1 but not in df2
missing_in_df2 <- setdiff(species1, species2)
missing_in_df2 

#[1] ""                "Pycnonotidae"    "Bycanistes"      "Cercopithecinae"    # we knew that

#sort(tree$tip.label)

###################### MGT - model 1 ################


# Fit PGLS model
#ML=maximum likelihood

# factor the groups so that rodent is the comparasion
unique(GIT_pruned$re_class)
# [1] "Flying Birds"        "Carnivores"          "Even-toed Ungulates" "Primates"            "Bats"               
# [6] "Reptiles"            "Non-Flying Birds"    "Rodents"             "Odd-toed Ungulates" 

length(unique(GIT_pruned$phylo_name)) #243


#GIT_pruned$re_class <- factor(GIT_pruned$re_class, levels = c("Rodents","Flying Birds","Bats","Carnivores",  "Non-Flying Birds",  "Primates", "Ungulates", "Reptiles"))

#transform data
GIT_pruned$log_transit_hrs <- log10(GIT_pruned$transit_hrs)


# add the "flyer" parameter
GIT_pruned$flyer <- 0
GIT_pruned$flyer[GIT_pruned$re_class=="Flying Birds" |GIT_pruned$re_class=="Bats"] <- 1

GIT_pruned$flyer <- as.factor(GIT_pruned$flyer)

# Create correlation structure
# this is full Brownian motion assumption of phylogentic effects:
cor_phylo_fixed1 <- corPagel(1, phy = tree_pruned, fixed = TRUE, form = ~phylo_name) #assumes phylogeny equivalent to brownian motion

# this is no phylogentic effects:
cor_phylo_fixed0 <- corPagel(0, phy = tree_pruned, fixed = TRUE, form = ~phylo_name) #assumes no phylogenetic signal



# Cara: now test the effect of flyer on GIT transit without accounting for mass or diet
# All text below is from Cara

#with full phylogenetic effects:
gls_model_GIT_1 <- gls(log_transit_hrs ~ flyer, 
                         data = GIT_pruned, 
                        correlation = cor_phylo_fixed1,
                         method = "ML")

summary(gls_model_GIT_1)$tTable

#CARA: I get this. EMILY, same
# Value Std.Error   t-value      p-value
# (Intercept)  1.654959 0.2416697  6.848022 4.550490e-11
# flyer1      -1.148172 0.1531520 -7.496942 8.173648e-13

# removing correlation structure, as its the same and the lambda=0 breaks the matrix
# with no phylogenetic effects:
gls_model_GIT_0 <- gls(log_transit_hrs ~ flyer, 
                        data = GIT_pruned, 
                        #correlation = cor_phylo_fixed0,
                        method = "ML")
summary(gls_model_GIT_0)$tTable
# Value  Std.Error   t-value      p-value
# (Intercept)  1.200419 0.04133458  29.04151 2.006484e-87
# flyer1      -1.324838 0.06468766 -20.48053 4.440990e-58

#CARA: I get the same



AIC(gls_model_GIT_0, gls_model_GIT_1)

#CARA: I get this
# df       AIC
# gls_model_GIT_0  3   468.540
# gls_model_GIT_1  3 -2285.136

# Emily: flyer is a significant term which says that flight is significantly
# negatively related to GIT transit time.
# R3: This is on top of a model that already includes significant 
# phylogenetic clustering in transit time, suggesting that flight affects transit
# time independent of phylogeny. 
# CARA: I agree



# now, let's also account for diet
# first try as a fixed effect:
GIT_pruned$trial.diet <- as.factor(GIT_pruned$trial.diet)

# Cara: with full phylogenetic effects:
gls_model_GIT_1_diet <- gls(log_transit_hrs ~ flyer + trial.diet, 
                        data = GIT_pruned, 
                        correlation = cor_phylo_fixed1,
                        method = "ML")
summary(gls_model_GIT_1_diet)$tTable


#CARA: I get this
# Value  Std.Error    t-value      p-value
# (Intercept)                    1.69576552 0.24065823  7.0463640 1.411006e-11
# flyer1                        -1.11432166 0.15347891 -7.2604219 3.768194e-12

# trial.dietfruit/nectar/pollen -0.11983615 0.07092630 -1.6895870 9.221217e-02
# trial.dietmeat                -0.03425177 0.09341499 -0.3666624 7.141458e-01
# trial.dietmixed               -0.18716642 0.23077848 -0.8110220 4.180370e-01
# trial.dietprotein             -0.10804018 0.08057619 -1.3408450 1.810498e-01
# trial.dietunknown              0.26001617 0.11503338  2.2603541 2.456190e-02

# Cara: with no phylogenetic effects:
gls_model_GIT_0_diet <- gls(log_transit_hrs ~ flyer + trial.diet, 
                        data = GIT_pruned, 
                        #correlation = cor_phylo_fixed0,
                        method = "ML")
summary(gls_model_GIT_0_diet)$tTable

# Value  Std.Error     t-value      p-value
# (Intercept)                    1.23140059 0.06328795  19.4571087 4.697948e-54
# flyer1                        -1.35614045 0.07220339 -18.7822264 1.314605e-51
# trial.dietfruit/nectar/pollen -0.11466745 0.09216882  -1.2441024 2.144950e-01
# trial.dietmeat                 0.09355365 0.10496697   0.8912675 3.735455e-01
# trial.dietmixed               -0.16921235 0.31552316  -0.5362914 5.921800e-01
# trial.dietprotein             -0.08527729 0.09517346  -0.8960197 3.710065e-01
# trial.dietunknown              0.25025965 0.14189221   1.7637307 7.886014e-02

#CARA: I get the same


AIC(gls_model_GIT_0, gls_model_GIT_1, gls_model_GIT_0_diet, gls_model_GIT_1_diet)

#CARA: I get this
# df        AIC
# gls_model_GIT_0       3   468.5400
# gls_model_GIT_1       3 -2285.1360
# gls_model_GIT_0_diet  8   467.0779
# gls_model_GIT_1_diet  8 -2291.5858

# CARA: conclusion is still the same: gls_model_GIT_1_diet is best model but 
# gls_model_GIT_1 is also pretty good. either way, seems like including flight
# is important even after phylogeny is accounted for



summary(gls_model_GIT_1_diet)$tTable

#CARA: I get this
# Value  Std.Error    t-value      p-value
# (Intercept)                    1.69576552 0.24065823  7.0463640 1.411006e-11
# flyer1                        -1.11432166 0.15347891 -7.2604219 3.768194e-12
# trial.dietfruit/nectar/pollen -0.11983615 0.07092630 -1.6895870 9.221217e-02
# trial.dietmeat                -0.03425177 0.09341499 -0.3666624 7.141458e-01
# trial.dietmixed               -0.18716642 0.23077848 -0.8110220 4.180370e-01
# trial.dietprotein             -0.10804018 0.08057619 -1.3408450 1.810498e-01
# trial.dietunknown              0.26001617 0.11503338  2.2603541 2.456190e-02

# Emily: model gls_model_GIT_1_diet has the best fit.
# CARA: I agree with your conclusion


# Cara: let's try it also as a random effect - here we use the lme function which
# is also supported in the nlme package. we can compare these since all were fit via ML:

gls_model_GIT_1_diet_re <- lme(log_transit_hrs ~ flyer,
                            random = ~1 | trial.diet,
                            correlation = cor_phylo_fixed1,
                            data = GIT_pruned,
                            method = "ML")

summary(gls_model_GIT_1_diet_re)$tTable
# Value Std.Error  DF   t-value      p-value
# (Intercept)  1.457103 0.1867034 282  7.804375 1.174627e-13
# flyer1      -1.044916 0.1601138 282 -6.526085 3.124348e-10

#CARA: I get the same

gls_model_GIT_0_diet_re <- lme(log_transit_hrs ~ flyer,
                                random = ~1 | trial.diet,
                                #correlation = cor_phylo_fixed0,
                                data = GIT_pruned,
                                method = "ML")
summary(gls_model_GIT_0_diet_re)$tTable
# Value  Std.Error  DF   t-value      p-value
# (Intercept)  1.216819 0.05492749 282  22.15318 1.108390e-63
# flyer1      -1.327690 0.06787143 282 -19.56184 1.964575e-54

#CARA : I get the same


AIC(gls_model_GIT_0, gls_model_GIT_0_diet, gls_model_GIT_0_diet_re,
    gls_model_GIT_1, gls_model_GIT_1_diet, gls_model_GIT_1_diet_re)


#CARA: I get this
# df        AIC
# gls_model_GIT_0          3   468.5400
# gls_model_GIT_0_diet     8   467.0779
# gls_model_GIT_0_diet_re  4   469.2664
# gls_model_GIT_1          3 -2285.1360
# gls_model_GIT_1_diet     8 -2291.5858
# gls_model_GIT_1_diet_re  4  -490.7875


#Emily: in all cases, the fit with diet as a fixed effect was better supported than the fit as
# a random effect. 

#CARA: conclusion is the same as yours (fixed effect better than random)


# Cara: The models with the lambda=1, which incorporate phylogenetic clustering
# of GIT transit time, are best supported. The top fit model, gls_model_GIT_1_diet,
# offers the best fit overall.

# you see that the significance of "flyer" is slightly improved when diet is 
# also modeled--making it significant** by most standards. 

# Because the flyer p-value improved with the addition of diet,
# I am guessing there are some fruit/nectar/pollen diet non-flyers with fairly fast 
# transit, so this additional variable helps clarify the what traits are correlated 
# with rapid transit and when. Thus, to conclude, flight has independently significant
# negative effect on GIT transit even after accounting for both phylogeny and diet.


###################
########## model gls_model_GIT_1_diet is best for reporting!


# now, let's also consider the effects of mass and, later, mass and diet


 
#with full phylogenetic effects:
gls_model_GIT_1_mass <- gls(log_transit_hrs ~ log10(mass_kg)*flyer, 
                        data = GIT_pruned, 
                        correlation = cor_phylo_fixed1,
                        method = "ML")
summary(gls_model_GIT_1_mass)$tTable

#CARA: I get this
# Value  Std.Error   t-value      p-value
# (Intercept)            1.66308317 0.23933097  6.948884 2.505479e-11
# log10(mass_kg)         0.12024435 0.03848077  3.124790 1.962812e-03
# flyer1                -1.03144545 0.17572566 -5.869635 1.211922e-08
# log10(mass_kg):flyer1 -0.08382952 0.06564203 -1.277071 2.026167e-01


#with no phylogenetic effects:
gls_model_GIT_0_mass <- gls(log_transit_hrs ~ log10(mass_kg)*flyer, 
                        data = GIT_pruned, 
                        #correlation = cor_phylo_fixed0,
                        method = "ML")
summary(gls_model_GIT_0_mass)$tTable
# Value  Std.Error    t-value      p-value
# (Intercept)            1.18023808 0.04172500  28.286116 1.001072e-84
# log10(mass_kg)         0.09058320 0.03507337   2.582677 1.030250e-02
# flyer1                -1.26415593 0.10776374 -11.730810 3.415045e-26
# log10(mass_kg):flyer1 -0.06842639 0.05878842  -1.163943 2.454204e-01

#CARA: I get the same



AIC(gls_model_GIT_0_mass, gls_model_GIT_1_mass)
#Emily: model 1 is the best
#Cara: This shows that transit time is still influenced by phylogeny, not
# exclusively by mass (e.g. the model with phylogeny AND mass performed better
# than the model with mass alone)

#CARA: I get this
# df        AIC
# gls_model_GIT_0_mass  5   465.6357
# gls_model_GIT_1_mass  5 -2291.3548

#CARA: conclusion is the same as yours

summary(gls_model_GIT_1_mass)$tTable

#CARA: I get this
# Value  Std.Error   t-value      p-value
# (Intercept)            1.66308317 0.23933097  6.948884 2.505479e-11
# log10(mass_kg)         0.12024435 0.03848077  3.124790 1.962812e-03
# flyer1                -1.03144545 0.17572566 -5.869635 1.211922e-08
# log10(mass_kg):flyer1 -0.08382952 0.06564203 -1.277071 2.026167e-01
#CARA: I agree with conclusion- values different

gls_model_GIT_1_mass_simple <- gls(log_transit_hrs ~ flyer + log10(mass_kg), 
                             data = GIT_pruned, 
                             correlation = cor_phylo_fixed1,
                             method = "ML")
summary(gls_model_GIT_1_mass_simple)$tTable

# Value  Std.Error   t-value      p-value
# (Intercept)     1.63554928 0.23862062  6.854182 4.410320e-11
# flyer1         -0.95377199 0.16504606 -5.778823 1.963752e-08
# log10(mass_kg)  0.09124249 0.03109935  2.933904 3.618423e-03

#CARA: I get the same

# mass is now significant

AIC(gls_model_GIT_1_mass, gls_model_GIT_1_mass_simple)

#CARA: I get this
# df       AIC
# gls_model_GIT_1_mass         5 -2291.355
# gls_model_GIT_1_mass_simple  4 -2291.706

# identical




#so I'm going to write up the simple model because they are so close and significance is across the board
# (no interaction term), shown here:
summary(gls_model_GIT_1_mass_simple)$tTable


#CARA: I get this
# Value  Std.Error   t-value      p-value
# (Intercept)     1.63554928 0.23862062  6.854182 4.410320e-11
# flyer1         -0.95377199 0.16504606 -5.778823 1.963752e-08
# log10(mass_kg)  0.09124249 0.03109935  2.933904 3.618423e-03

# flyer is significantly negatively related to transit suggesting 
# that flyers have rapid transit independent of small body size; 
# however, there is also significance across the board, and its significant in the complex model too
#CARA : I agree

#

# now, let's see if adding diet improves this model:


# first, as a fixed effect
gls_model_GIT_1_mass_simple_diet <- gls(log_transit_hrs ~ flyer + log10(mass_kg) + trial.diet, 
                                    data = GIT_pruned, 
                                    correlation = cor_phylo_fixed1,
                                    method = "ML")
summary(gls_model_GIT_1_mass_simple_diet)$tTable

#CARA: I get this
# Value  Std.Error    t-value      p-value
# (Intercept)                    1.68901518 0.23737493  7.1153898 9.306649e-12
# flyer1                        -0.92282513 0.16445144 -5.6115357 4.804738e-08
# log10(mass_kg)                 0.09548708 0.03204079  2.9801726 3.132814e-03
# trial.dietfruit/nectar/pollen -0.09102744 0.07062021 -1.2889715 1.984681e-01
# trial.dietmeat                -0.07715611 0.09325430 -0.8273732 4.087276e-01
# trial.dietmixed               -0.20764524 0.22772330 -0.9118313 3.626391e-01
# trial.dietprotein             -0.07877702 0.08007758 -0.9837587 3.260804e-01
# trial.dietunknown              0.29222516 0.11397241  2.5639991 1.086761e-02



# then, as a random effect

gls_model_GIT_1_mass_simple_diet_re <- lme(log_transit_hrs ~ flyer + log10(mass_kg),
                                            random = ~1 | trial.diet,
                                            data = GIT_pruned, 
                                            correlation = cor_phylo_fixed1,
                                            method = "ML")
summary(gls_model_GIT_1_mass_simple_diet_re)$tTable
# Value  Std.Error  DF   t-value      p-value
# (Intercept)     1.4492430 0.18716737 281  7.743032 1.766444e-13
# flyer1         -0.8066961 0.16996020 281 -4.746383 3.302343e-06
# log10(mass_kg)  0.1115418 0.03060434 281  3.644640 3.188962e-04

#CARA: I get the same


# then, compare
AIC(gls_model_GIT_1_mass_simple, gls_model_GIT_1_mass_simple_diet, gls_model_GIT_1_mass_simple_diet_re)
# EMILY
# df        AIC
# gls_model_GIT_1_mass_simple          4 -1880.5901
# gls_model_GIT_1_mass_simple_diet     9 -1887.4631
# gls_model_GIT_1_mass_simple_diet_re  5  -511.7503

#CARA: I get this:
# df        AIC
# gls_model_GIT_1_mass_simple          4 -2291.7057
# gls_model_GIT_1_mass_simple_diet     9 -2298.5787
# gls_model_GIT_1_mass_simple_diet_re  5  -501.9003

# here, mass and diet does improve model fit when diet is modeled as a fixed effect! 
#CARA: conclusion of mine is the same as yours
summary(gls_model_GIT_1_mass_simple_diet)$tTable

#CARA: I get this
# Value  Std.Error    t-value      p-value
# (Intercept)                    1.68901518 0.23737493  7.1153898 9.306649e-12
# flyer1                        -0.92282513 0.16445144 -5.6115357 4.804738e-08
# log10(mass_kg)                 0.09548708 0.03204079  2.9801726 3.132814e-03
# trial.dietfruit/nectar/pollen -0.09102744 0.07062021 -1.2889715 1.984681e-01
# trial.dietmeat                -0.07715611 0.09325430 -0.8273732 4.087276e-01
# trial.dietmixed               -0.20764524 0.22772330 -0.9118313 3.626391e-01
# trial.dietprotein             -0.07877702 0.08007758 -0.9837587 3.260804e-01
# trial.dietunknown              0.29222516 0.11397241  2.5639991 1.086761e-02

# Emily: flyer is still significantly negatively related to transit, 
# after accounting for mass and diet (and phylogeny). 
# The unknown diets are significant but nothing else -
# but their inclusion does improve the overall model, so it is best to keep them.

#CARA: I agree. my output is differe. note tht mass is also sig here

# R3: I'll also note that including diet here again
# improved the significance of flyer but had no effect on mass, supporting my 
# prior hypothesis that there are probably some non-flyers with fruit/nectar/pollen diets
# that had fairly fast transit time and meat eaters with slow transit time

# Cara: what if we compare our two best fit models so far? 
AIC(gls_model_GIT_1_diet, gls_model_GIT_1_mass_simple_diet)

# df       AIC
# gls_model_GIT_1_diet              8 -1880.470
# gls_model_GIT_1_mass_simple_diet  9 -1887.463

#CARA: I get this
# df       AIC
# gls_model_GIT_1_diet              8 -2291.586
# gls_model_GIT_1_mass_simple_diet  9 -2298.579


gls_model_GIT_1_mass_only <- gls(log_transit_hrs ~ log10(mass_kg), 
                                   data = GIT_pruned, 
                                   correlation = cor_phylo_fixed1,
                                   method = "ML")

gls_model_GIT_0_mass_only <- gls(log_transit_hrs ~ log10(mass_kg), 
                                 data = GIT_pruned, 
                                 #correlation = cor_phylo_fixed0,
                                 method = "ML")



summary(gls_model_GIT_0_mass_only)$tTable
#mass is very significant and positively related to gut transit

# Value  Std.Error  t-value      p-value
# (Intercept)    0.8596155 0.04385958 19.59927 7.109775e-55
# log10(mass_kg) 0.3256692 0.02695168 12.08345 1.881124e-27

#CARA: I get the same as above. I agree with conclusion

summary(gls_model_GIT_1_mass_only)$tTable
# mass is hugely significant and positively related to gut transit

#EMILY: mass is hugely significant again

#CARA: I get this
# Value  Std.Error  t-value      p-value
# (Intercept)    1.4368124 0.24909960 5.768024 2.073395e-08
# log10(mass_kg) 0.1633924 0.03004766 5.437773 1.155222e-07

#CARA: conclusion still the same
AIC(gls_model_GIT_0_mass_only, gls_model_GIT_1_mass_only)
# incorporating phylogeny with mass offers a much better fit to the data however!
# mass is positively correlated with git transit

# df        AIC
# gls_model_GIT_0_mass_only  3   610.0041
# gls_model_GIT_1_mass_only  3 -1850.6741

#CARA: I get this
# df        AIC
# gls_model_GIT_0_mass_only  3   610.0041
# gls_model_GIT_1_mass_only  3 -2261.7897


#CARA: I agree with conclusions. including phylogeny is a better fit...
#mass positively correlated with GIT transit


###### ------ PLOTTING FOR FIGURE 2
library(nlme)
library(ggplot2)
library(dplyr)

# --- Colors and shapes ---
# [1] "Flying Birds"        "Carnivores"          "Even-toed Ungulates" "Primates"            "Bats"               
# [6] "Reptiles"            "Rodents"             "Odd-toed Ungulates"  

colz <- c(
  "Flying\nBirds"= "#F8766D",
  "Bats" = "#C49A00",
  "Non-Flying\nBirds" = "#edf8b1",
  "Rodents" = "navy",
  "Primates" = "#00B6EB",
  "Even-toed\nUngulates" = "#A58AFF",
  "Reptiles" = "#FB61D7",
  "Carnivores" = "#00C094","Odd-toed\nUngulates"="thistle1"
)

unique(GIT_pruned$trial.diet)

shapez <- c(
  "protein" = 21,
  "meat" = 25,
  "fruit/nectar/pollen" = 22,
  "unknown" = 24,
  "mixed" = 8,
  "fiber/foliage" = 23
)

# --- Fix factor levels ---
# Convert to character to safely modify
GIT_pruned$re_class2 <- as.character(GIT_pruned$re_class)

unique(GIT_pruned$re_class2)
# [1] "Flying Birds"        "Carnivores"          "Even-toed Ungulates" "Primates"            "Bats"               
# [6] "Reptiles"            "Non-Flying Birds"    "Rodents"             "Odd-toed Ungulates"   

# Replace the value
GIT_pruned$re_class2[GIT_pruned$re_class2 == "Flying Birds"] <- "Flying\nBirds"
GIT_pruned$re_class2[GIT_pruned$re_class2 == "Non-Flying Birds"] <- "Non-Flying\nBirds"
#GIT_pruned$re_class2[GIT_pruned$re_class2 == "Flying Birds"] <- "Flying\nBirds"
GIT_pruned$re_class2[GIT_pruned$re_class2 == "Even-toed Ungulates"] <- "Even-toed\nUngulates"
GIT_pruned$re_class2[GIT_pruned$re_class2 == "Odd-toed Ungulates"] <- "Odd-toed\nUngulates"

GIT_pruned$re_class2 <- factor(GIT_pruned$re_class2, levels = c("Flying\nBirds","Bats", "Rodents", "Carnivores",
                                                                "Non-Flying\nBirds",
                                                              "Primates", "Even-toed\nUngulates", "Odd-toed\nUngulates", "Reptiles"))

# --- Plot ---
p1 <- ggplot(GIT_pruned) +
  geom_boxplot(aes(x = re_class2, y = log10(transit_hrs), fill = re_class2), 
               outlier.shape = NA, show.legend = FALSE) +
  geom_point(aes(x = re_class2, y = log10(transit_hrs), shape = trial.diet),
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
             size = 3) +
  scale_fill_manual(values = colz) +
  scale_shape_manual(values = shapez, name = "Trial diet") +
  ylab(expression(Log[10] ~ "GIT transit time (hrs)")) +
  coord_cartesian(ylim = c(-1, 3.5)) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid = element_blank(),
    legend.position = c(0.15, 0.8),
    legend.background = element_rect(color = "black"),
    plot.margin = unit(c(0.2, 0.2, 0.8, 0.2), "cm")
  )

print(p1)


ggsave(file = paste0(homewd,"/figures/revision-3/Fig_2A.png"),
       units="mm",  
       width=120, 
       height=80, 
       scale=3, 
       dpi=300)




##### panel b with mass - mean centered mass 
library(nlme)
library(ggplot2)
library(dplyr)

# --- Prepare dataset with mean-centered log_mass ---
GIT_pruned <- GIT_pruned %>%
  mutate(
    log_mass = log10(mass_kg),
    log_mass_c = log_mass - mean(log_mass, na.rm = TRUE)
  )

# --- Fit GLS model ---
gls_model <- gls(
  log_transit_hrs ~ flyer +log_mass_c,
  data = GIT_pruned,
  correlation = cor_phylo_fixed1,
  method = "ML"
)

summary(gls_model)
# Coefficients:
#   Value  Std.Error   t-value p-value
# (Intercept)  1.6666752 0.24539620  6.791773  0.0000
# flyer1      -0.9014954 0.16841113 -5.352944  0.0000
# log_mass_c   0.0897846 0.03178465  2.824780  0.0051



# --- Build prediction grid ---
newdat <- expand.grid(
  log_mass_c = seq(min(GIT_pruned$log_mass_c, na.rm = TRUE),
                   max(GIT_pruned$log_mass_c, na.rm = TRUE),
                   length.out = 100),
  flyer = unique(GIT_pruned$flyer)
)

# --- Get predicted values ---
newdat$fit <- predict(gls_model, newdata = newdat)

# --- Plot ---
p2 <-
  ggplot(GIT_pruned, aes(x = log_mass_c, y = log_transit_hrs)) +
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
  #geom_line(data = newdat,
            # aes(x = log_mass_c, y = fit, color = flyer),
            # linewidth = 1.2) +
  # Vertical line at mean mass
  scale_color_manual(values = c("darkblue", "#E08929")) +
  theme_classic(base_size = 14) +
  labs(
    x = "\nMean-centered Log10 body mass (kg)",
    y = "Log GIT transit time (hrs)\n",
    shape = "Diet",
    color = "Flyer"
  ) + theme(legend.position= "none")

p2

# --- Save the plot ---
ggsave("figures/revision-3/GIT_plot2b.png", plot = p2, width = 7, height = 6, dpi = 300)




##### cowplot the figures
out.plot2 <- cowplot::plot_grid(p1, p2, nrow=1, ncol=2, labels=c("(A)", "(B)"), rel_widths = c(1.2,1.2))

print(out.plot2)

ggsave(file = paste0(homewd,"/figures/revision-3/Fig2_TwoPanel_R2.jpeg"),
       units="mm",  
       width=170, 
       height=70, 
       scale=3, 
       dpi=300)

