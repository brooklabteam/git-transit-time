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
homewd <- "/Users/emilyruhs/Desktop/GitHub_repos/git-transit-time/"
homewd <- "/Users/katherinemcferrin/Developer/git-transit-time/"

setwd(homewd)

#load the GIT transit data with phylo name
dat_GIT <- read.csv(file = paste0(homewd, "data/revision-1/dat_sum_tot_clean_GIT_R2.csv"), header = T, stringsAsFactors = F )
dat_MRT <- read.csv(file = paste0(homewd, "data/revision-1/dat_sum_tot_clean_MRT_R2.csv"), header = T, stringsAsFactors = F )

#remove the two species that don't have masses earlier on - While they have data on diet, it will be
# easier to explain if we just remove them from the offset.
dat_GIT <- subset(dat_GIT, phylo_name!="Uromacer_oxyrhynchus")
dat_GIT <- subset(dat_GIT, phylo_name!="Atheris_squamigera")


#### species counts ####
#table 1 in the manuscript

length(unique(dat_GIT$genus.species)) #112 unique species for transit time, this is down 2 species due to unknown diet
(tt.pivot <-ddply(dat_GIT, .(re_class), summarise, N_species_transit = length(unique(genus.species))))
# re_class N_species_transit
# 1         Bats                36
# 2   Carnivores                 6
# 3 Flying Birds                13
# 4     Primates                21
# 5     Reptiles                16 #used to be 18
# 6      Rodents                13
# 7    Ungulates                 7

length(unique(dat_MRT$genus.species)) #60 unique species for mean retention time
(mrt.pivot <-ddply(dat_MRT, .(re_class), summarise, N_species_mrt = length(unique(genus.species))))
# re_class N_species_mrt
# 1             Bats            12
# 2 Non-Flying Birds             3
# 3         Primates            19
# 4          Rodents            14
# 5        Ungulates            12



#load the tree - pulled from timetree so should be in principle an ultrametric tree
tree <- read.tree(file = paste0(homewd, "data/revision-1/timetree_names.nwk"))
# is.ultrametric(tree) #check for ultrametric - is FALSE
# 
# tree_ultra <- chronos(tree)
# is.ultrametric(tree_ultra)

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
## Emily: I will remove the species above and those without diets for the final paper tree
plot(tree_ultra2,
     type = "phylogram",    # or "fan", "cladogram", "unrooted"
     cex = 0.6,             # Tip label size
     no.margin = TRUE,      # Remove extra white space
     edge.width = 1.5,      # Thicker branches
     edge.color = "darkgray", # Branch color
     label.offset = 0.001,  # Space between tip and label
     font = 3,              # Italic tip labels
     main = "")

# Add a scale bar with adjusted size and location
add.scale.bar(length = 0.05, lwd = 2, cex = 0.7, col = "black")

# 
# 
# 



 
##################  MGT phylogenetic signal ##########
dat_GIT$phylo_name <- as.character(dat_GIT$phylo_name)

#for GIT
GIT <- dat_GIT[c("phylo_name", "transit_hrs")]
#GIT <- na.omit(GIT)

#make sure they match
# If you have a data frame GIT with the trait
transit_hrs <- setNames(GIT$transit_hrs, GIT$phylo_name)

# Then only keep the overlapping species
common_species <- intersect(tree_ultra2$tip.label, names(transit_hrs)) #107 here. there are 120 entries for transit hours and 124 tips in tree_ultra2
tree <- drop.tip(tree_ultra2, setdiff(tree$tip.label, common_species)) #from Cara: it's the same here but to be fully accurate, shouldn't you do setdiff(tree_ultra2$tip.label, ...) rather than tree$tip.label?
#tree <- drop.tip(tree_ultra, setdiff(tree$tip.label, names(transit_hrs)))
transit_hrs <- transit_hrs[common_species]
setdiff(names(transit_hrs), tree$tip.label)#none
setdiff(tree$tip.label, names(transit_hrs))#none

lambda_gs1<-phylosig(tree, transit_hrs,
                    method="lambda",test=TRUE)
lambda_gs1 
# Cara: my understanding is that this is the number (1.00745) you'd want 
# to put in your phylogenetic model (rather than 1) but I could be wrong here, 
# Emily: lambda is now so close to one that we will use lambda = 1 in the model, and not need model 2
# Cara: a little research suggests that values >1 are pretty rare and since
# this is basically 1, that seems fine to me

# Phylogenetic signal lambda : 1.00745
# logL(lambda) : -731.72
# LR(lambda=0) : 54.7237 
# P-value (based on LR test) : 1.38721e-13 

#tree$tip.label

contMap(tree,log10(transit_hrs), fsize =0.6)
# save as MGT_contmap
# 





################### MRT phylogenetic signal ##########################
dat_MRT$phylo_name <- as.character(dat_MRT$phylo_name)

tree <- read.tree(file = paste0(homewd, "data/revision-1/timetree_names.nwk"))
tree_ultra2 <- force.ultrametric(tree, method = "nnls") #no warnings
is.ultrametric(tree_ultra2)

#for MRT
MRT <- dat_MRT[c("phylo_name", "MRT_hrs")]
#MRT <- na.omit(MRT)

#make sure they match
# If you have a data frame GIT with the trait
MRT_hrs <- setNames(MRT$MRT_hrs, MRT$phylo_name)

# Then only keep the overlapping species
common_species <- intersect(tree_ultra2$tip.label, names(MRT_hrs))
tree <- drop.tip(tree_ultra2, setdiff(tree$tip.label, common_species)) #Cara: same as above, to keep it clean, I would write as setdiff(tree_ultra2$tip.label,...)
MRT_hrs <- MRT_hrs[common_species]

#for MRT
lambda_gs2<-phylosig(tree, MRT_hrs,
                     method="lambda",test=TRUE)
lambda_gs2
### put the 0.89 into the model instead of 1? 
#Cara: I think so -- this suggests the phylogenetic signal is a bit weaker for MRT vs MGT
#which makes sense because the data are fewer
# 
# Phylogenetic signal lambda : 0.890595 
# logL(lambda) : -234.723 
# LR(lambda=0) : 22.4312 
# P-value (based on LR test) : 2.17806e-06 



############# plotting phylo signal by tree/phenogram ############


contMap(tree,log10(MRT_hrs), fsize =0.6)


#phenogram(tree, log10(MRT_hrs), spread.labels = T, fsize=0.6, colors = "black")





#########################################  MGT. #################################

########### MGT data cleaning #######

#_________________________________________________________________________________________________________________
#I think this code runs PGLS and also provides lambda estimates as part of the model



# Load your data
dat_GIT <- read.csv(file = paste0(homewd, "/data/revision-1/dat_sum_tot_clean_GIT_R2.csv"), header = T, stringsAsFactors = F )
names(dat_GIT)

#remove these
dat_GIT <- subset(dat_GIT, phylo_name!="Uromacer_oxyrhynchus")
dat_GIT <- subset(dat_GIT, phylo_name!="Atheris_squamigera")


# clean dataset with only the columns you need
GIT <- dat_GIT[c("phylo_name", "transit_hrs", "trial.diet", "re_class", "mass_kg")]
#GIT <- na.omit(GIT$transit_hrs)

tree <- read.tree(paste0(homewd, "/data/revision-1/timetree_names.nwk"))
tree_ultra2 <- force.ultrametric(tree, method = "nnls") #no warnings
is.ultrametric(tree_ultra2)

#make sure they match
# Check if the species names match the tree tips
common_species <- intersect(tree_ultra2$tip.label, GIT$phylo_name)
tree_pruned <- drop.tip(tree_ultra2, setdiff(tree$tip.label, common_species))
GIT_pruned <- GIT[GIT$phylo_name %in% common_species, ]

# there are some duplicated taxa with different diets:
GIT_pruned[duplicated(GIT_pruned$phylo_name),]

# check to see that the names match
# Species in common_species but not in tree_pruned$tip.label
# missing_from_tree <- setdiff(common_species, tree_pruned$tip.label)
# # Print the missing species
# missing_from_tree 


species1 <- unique(dat_GIT$phylo_name)
species2 <- unique(GIT_pruned$phylo_name)

# Species in df1 but not in df2
missing_in_df2 <- setdiff(species1, species2)
missing_in_df2 

# [1] "Cercopithecus_albogularis" "Cercopithecus_mitis"       "Cercopithecus_pogonias"   
# [4] "Acanthagenys_rufogularis"  "Eptesicus_innoxius"     

#sort(tree$tip.label)

###################### MGT - model 1 ################


# Fit PGLS model
#ML=maximum likelihood

# factor the groups so that rodent is the comparasion
unique(GIT_pruned$re_class)
# "Primates"     "Bats"         "Reptiles"     "Ungulates"    "Flying Birds" "Rodents"      "Carnivores" 

length(unique(GIT_pruned$phylo_name)) #107


GIT_pruned$re_class <- factor(GIT_pruned$re_class, levels = c("Rodents","Flying Birds","Bats", "Carnivores", "Primates", "Ungulates", "Reptiles"))

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

# this is what was estimated from your data above - removing because lambda is 1 in the above model
#cor_phylo_fixed2 <- corPagel(lambda_gs1$lambda, phy = tree_pruned, fixed = TRUE, form = ~phylo_name) 


# Cara: now test the effect of flyer on GIT transit without accounting for mass or diet
# All text below is from Cara

#with full phylogenetic effects:
gls_model_GIT_1 <- gls(log_transit_hrs ~ flyer, 
                         data = GIT_pruned, 
                        correlation = cor_phylo_fixed1,
                         method = "ML")

#with no phylogenetic effects:
gls_model_GIT_0 <- gls(log_transit_hrs ~ flyer, 
                        data = GIT_pruned, 
                        correlation = cor_phylo_fixed0,
                        method = "ML")

#with estimated phylogenetic effects:
# gls_model_GIT_2 <- gls(log_transit_hrs ~ flyer, 
#                         data = GIT_pruned, 
#                         correlation = cor_phylo_fixed2,
#                         method = "ML")

AIC(gls_model_GIT_0, gls_model_GIT_1)
# Emily: the estimated phylogenetic effects (which happens to be lambda=1)
# provide a better fit than the no-phylo effects model. 

# df       AIC
# gls_model_GIT_0  3 -160.4626
# gls_model_GIT_1  3 -309.6786 ******top

summary(gls_model_GIT_1) 
# Coefficients:
#               Value    Std.Error   t-value  p-value
# (Intercept)  1.6179687 0.2872801  5.632025  0.0000
# flyer1      -0.9935241 0.4029689 -2.465511  0.0152

# Emily: flyer is a highly significant term which says that flight is significantly
# negatively related to GIT transit time.
# Cara: This is on top of a model that already includes significant 
# phylogenetic clustering in transit time, suggesting that flight affects transit
# time independent of phylogeny. I will note that I would not say it is "highly
# significant" but moderately (or even weakly*) so (I like the scale, p<0.001***, p<0.01**, p<0.1*)



# now, let's also account for diet
# first try as a fixed effect:
GIT_pruned$trial.diet <- as.factor(GIT_pruned$trial.diet)

# Cara: with full phylogenetic effects:
gls_model_GIT_1_diet <- gls(log_transit_hrs ~ flyer + trial.diet, 
                        data = GIT_pruned, 
                        correlation = cor_phylo_fixed1,
                        method = "ML")

# Cara: with no phylogenetic effects:
gls_model_GIT_0_diet <- gls(log_transit_hrs ~ flyer + trial.diet, 
                        data = GIT_pruned, 
                        correlation = cor_phylo_fixed0,
                        method = "ML")

# Cara: with estimated phylogenetic effects:
# gls_model_GIT_2_diet <- gls(log_transit_hrs ~ flyer + typical.diet, 
#                         data = GIT_pruned, 
#                         correlation = cor_phylo_fixed2,
#                         method = "ML")

AIC(gls_model_GIT_0, gls_model_GIT_1, gls_model_GIT_0_diet, gls_model_GIT_1_diet)
# df       AIC
# gls_model_GIT_0       3 -160.4626
# gls_model_GIT_1       3 -309.6786
# gls_model_GIT_0_diet  7 -184.2452
# gls_model_GIT_1_diet  7 -314.6670  ******top

# Emily: model gls_model_GIT_1_diet has the best fit, but only marginally. 
# Cara: an AIC difference of 4 or more is considered significant, so this
# is actually a sizerable improvement over the phylogenetic model without diet



# Cara: let's try it also as a random effect - here we use the lme function which
# is also supported in the nlme package. we can compare these since all were fit via ML:

gls_model_GIT_1_diet_re <- lme(log_transit_hrs ~ flyer,
                            random = ~1 | trial.diet,
                            correlation = cor_phylo_fixed1,
                            data = GIT_pruned,
                            method = "ML")

gls_model_GIT_0_diet_re <- lme(log_transit_hrs ~ flyer,
                                random = ~1 | trial.diet,
                                correlation = cor_phylo_fixed0,
                                data = GIT_pruned,
                                method = "ML")

# gls_model_GIT_2_diet_re <- lme(log_transit_hrs ~ flyer,
#                                 random = ~1 | typical.diet,
#                                 correlation = cor_phylo_fixed2,
#                                 data = GIT_pruned,
#                                 method = "ML")

AIC(gls_model_GIT_0, gls_model_GIT_0_diet, gls_model_GIT_0_diet_re,
    gls_model_GIT_1, gls_model_GIT_1_diet, gls_model_GIT_1_diet_re)
    #gls_model_GIT_2, gls_model_GIT_2_diet, gls_model_GIT_2_diet_re)

# df       AIC
# gls_model_GIT_0          3 -160.4626
# gls_model_GIT_0_diet     7 -184.2452
# gls_model_GIT_0_diet_re  4  194.8100
# gls_model_GIT_1          3 -309.6786
# gls_model_GIT_1_diet     7 -314.6670 **** top
# gls_model_GIT_1_diet_re  4  197.3335

#Emily: in all cases, the fit as a fixed effect was better supported than the fit as
# a random effect. 
# Cara: The models with the lambda=1, which incorporate phylogenetic clustering
# of GIT transit time, are best supported. The top fit model, gls_model_GIT_1_diet,
# offers the best fit overall. When you look at the results and compare:
summary(gls_model_GIT_1_diet) #flyer p=0.003
summary(gls_model_GIT_1)# flyer p=0.015
# you see that the significance of "flyer" is actually improved when diet is 
# also modeled--making it significant** by most standards (p=0.003). You also see that
# diet type "fruit/nectar/pollen" is weakly significant* (p=0.04) by some standards.
# This suggests that both flight and fruit/nectar/pollen diet are associated with
# faster GIT transit. Because the flyer p-value improved with the addition of diet,
# I am guessing there are some fruit/nectar/pollen diet non-flyers with fairly fast 
# transit, so this additional variable helps clarify the what traits are correlated 
# with rapid transit and when. Thus, to conclude, flight has independently significant
# negative effect on GIT transit even after accounting for both phylogeny and diet.


###################
########## model gls_model_GIT_1_diet is best for reporting.


#???????

# now, let's also consider the effects of mass and, later, mass and diet

# below won't run because the correlation structure is set by a 
# phylogeny of 109 species but the dataset has 117 entries due to some taxa being 
# reeated multiple times due tp differences in diet or mass reported each time
 
#with full phylogenetic effects:
gls_model_GIT_1_mass <- gls(log_transit_hrs ~ log10(mass_kg)*flyer, 
                        data = GIT_pruned, 
                        correlation = cor_phylo_fixed1,
                        method = "ML")

#with no phylogenetic effects:
gls_model_GIT_0_mass <- gls(log_transit_hrs ~ log10(mass_kg)*flyer, 
                        data = GIT_pruned, 
                        correlation = cor_phylo_fixed0,
                        method = "ML")

#with estimated phylogenetic effects:
# gls_model_GIT_2_mass <- gls(log_transit_hrs ~ log10(mass_kg)*flyer, 
#                         data = GIT_pruned, 
#                         correlation = cor_phylo_fixed2,
#                         method = "ML")

AIC(gls_model_GIT_0_mass, gls_model_GIT_1_mass)
#Emily: model 1 is the best
#Cara: This shows that transit time is still influenced by phylogeny, not
# exclusively by mass (e.g. the model with phylogeny AND mass performed better
# than the model with mass alone)

# df       AIC
# gls_model_GIT_0_mass  5 -160.8576
# gls_model_GIT_1_mass  5 -308.9503 ***** top

summary(gls_model_GIT_1_mass)

# Coefficients:
#                        Value      Std.Error   t-value p-value
# (Intercept)            1.6183264 0.2864556  5.649484  0.0000
# log10(mass_kg)         0.0648412 0.0507509  1.277636  0.2040
# flyer1                -0.7625409 0.4633963 -1.645548  0.1027
# log10(mass_kg):flyer1  0.0942537 0.1328935  0.709242  0.4797

# Cara: flyer is the only significant term here at all but weakly so.
# We can test whether it would be best to drop the interaction which 
# might clarify the variables that matter most.

gls_model_GIT_1_mass_simple <- gls(log_transit_hrs ~ flyer + log10(mass_kg), 
                             data = GIT_pruned, 
                             correlation = cor_phylo_fixed1,
                             method = "ML")

AIC(gls_model_GIT_1_mass, gls_model_GIT_1_mass_simple)

# df       AIC
# gls_model_GIT_1_mass         5 -308.9503
# gls_model_GIT_1_mass_simple  4 -310.4303

# there is no difference in AIC,really, though the simple model is slightly better
#so we prefer the simplest model 
# (no interaction term), shown here:
summary(gls_model_GIT_1_mass_simple)

# Coefficients:
#   Value Std.Error   t-value p-value
# (Intercept)     1.6312368 0.2852417  5.718788  0.0000
# flyer1         -0.9247574 0.4021237 -2.299684  0.0233
# log10(mass_kg)  0.0778024 0.0472413  1.646915  0.1024

# flyer is significantly negatively related to transit suggesting 
# that flyers have rapid transit independent of small body size; 
# however, there is no effect of mass alone on transit

# Cara: good! additionally, note that the significance of flyer improved
# in the simple model which I think clarifies the relationship.

# now, let's see if adding diet improves this model:


# I'll just illustrate with the best fit model
# first, as a fixed effect
gls_model_GIT_1_mass_simple_diet <- gls(log_transit_hrs ~ flyer + log10(mass_kg) + trial.diet, 
                                    data = GIT_pruned, 
                                    correlation = cor_phylo_fixed1,
                                    method = "ML")

# then, as a random effect

gls_model_GIT_1_mass_simple_diet_re <- lme(log_transit_hrs ~ flyer + log10(mass_kg),
                                            random = ~1 | trial.diet,
                                            data = GIT_pruned, 
                                            correlation = cor_phylo_fixed1,
                                            method = "ML")

# then, compare
AIC(gls_model_GIT_1_mass_simple, gls_model_GIT_1_mass_simple_diet, gls_model_GIT_1_mass_simple_diet_re)

# df       AIC
# gls_model_GIT_1_mass_simple          4 -310.4303
# gls_model_GIT_1_mass_simple_diet     8 -313.0891 *** top 
# gls_model_GIT_1_mass_simple_diet_re  5  192.0029

# here, mass and diet does improve model fit when diet is modeled as a fixed effect! however,
# it is extremely close if we consider the penality of AIC +2
summary(gls_model_GIT_1_mass_simple_diet)

#                                Value     Std.Error   t-value p-value
# (Intercept)                    1.9070381 0.4436281  4.298731  0.0000
# flyer1                        -1.2358346 0.4247817 -2.909341  0.0044
# log10(mass_kg)                 0.0318619 0.0505564  0.630225  0.5299
# trial.dietfiber/foliage       -0.4175332 0.3235691 -1.290399  0.1997
# trial.dietfruit/nectar/pollen -0.5795950 0.2996161 -1.934459  0.0557
# trial.dietmeat                -0.0826795 0.3425893 -0.241337  0.8098
# trial.dietprotein             -0.3997218 0.3327219 -1.201369  0.2322

# Emily: flyer is still significantly negatively related to transit, 
# after accounting for mass and diet (and phylogeny). 
#None of the individual effects of diet are significant, 
# but their inclusion does improve the overall model, so it is best to keep them.
# the fruit/nectar/pollen diet has a near-significant negative relationship
# with transit (independent of flying status and mass) which makes intuitive sense.

# Cara: yes, good interpretation. I'll also note that including diet here again
# improved the significance of flyer but had no effect on mass, supporting my 
# prior hypothesis that there are probably some non-flyers with fruit/nectar/pollen diets
# that had fairly fast transit time.

# Cara: what if we compare our two best fit models so far? 
AIC(gls_model_GIT_1_diet, gls_model_GIT_1_mass_simple_diet)
# These are really very close; however, the simpler model which includes
# diet but not mass still has a better AIC, and mass was not significantly
# related to transit in our models, so the model without it
# (gls_model_GIT_1_diet) is probably the best overall.


## Emily and Katherine questions: 
# (1) previously we had structured the methods/results as though there were two separate models being run - one # for the diet and one for the mass. Do you still view that as true, or are we presenting these above models 
# as all part of the AIC model comparison? 

# (2) just a note, we removed the 'slim' models because we no longer had to fit separate datasets, as we 
# removed those two species at the offset. (these species had a NA for mass)
# (3) second note, lambda is now basically 1 (set by our data), so we no longer need the model twos where the 
# data were deciding the lambda.