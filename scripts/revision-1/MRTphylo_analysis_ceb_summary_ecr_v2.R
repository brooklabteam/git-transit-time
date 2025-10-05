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
# dat_GIT <- read.csv(file = paste0(homewd, "data/dat_sum_tot_clean_GIT_R2.csv"), header = T, stringsAsFactors = F )
dat_MRT <- read.csv(file = paste0(homewd, "data/revision-1/dat_sum_tot_clean_MRT_R2.csv"), header = T, stringsAsFactors = F )

#remove the two species that don't have masses earlier on - While they have data on diet, it will be
# easier to explain if we just remove them from the offset.
# dat_GIT <- subset(dat_GIT, phylo_name!="Uromacer_oxyrhynchus")
# dat_GIT <- subset(dat_GIT, phylo_name!="Atheris_squamigera")


#### species counts ####
#table 1 in the manuscript

# length(unique(dat_GIT$genus.species)) #112 unique species for transit time, this is down 2 species due to unknown diet
# (tt.pivot <-ddply(dat_GIT, .(re_class), summarise, N_species_transit = length(unique(genus.species))))
# # re_class N_species_transit
# # 1         Bats                36
# # 2   Carnivores                 6
# # 3 Flying Birds                13
# # 4     Primates                21
# # 5     Reptiles                16 #used to be 18
# # 6      Rodents                13
# # 7    Ungulates                 7

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
tree_ultra2 <- force.ultrametric(tree, method = "nnls") #no warnings

is.ultrametric(tree_ultra2) # TRUE



##### plot the phylogeny for GIT and MRT #######

# Plot the tree with enhancements

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
# dat_GIT$phylo_name <- as.character(dat_GIT$phylo_name)
# 
# #for GIT
# GIT <- dat_GIT[c("phylo_name", "transit_hrs")]
# #GIT <- na.omit(GIT)
# 
# #make sure they match
# # If you have a data frame GIT with the trait
# transit_hrs <- setNames(GIT$transit_hrs, GIT$phylo_name)
# 
# # Then only keep the overlapping species
# common_species <- intersect(tree_ultra2$tip.label, names(transit_hrs)) #107 here. there are 120 entries for transit hours and 124 tips in tree_ultra2
# 
# tree <- drop.tip(tree_ultra2, setdiff(tree_ultra2$tip.label, common_species)) 
# #tree <- drop.tip(tree_ultra, setdiff(tree$tip.label, names(transit_hrs)))
# transit_hrs <- transit_hrs[common_species]
# setdiff(names(transit_hrs), tree$tip.label)#none
# setdiff(tree$tip.label, names(transit_hrs))#none
# 
# lambda_gs1<-phylosig(tree, transit_hrs,
#                     method="lambda",test=TRUE)
# lambda_gs1 
# # Cara: my understanding is that this is the number (1.00745) you'd want 
# # to put in your phylogenetic model (rather than 1) but I could be wrong here, 
# # Emily: lambda is now so close to one that we will use lambda = 1 in the model, and not need model 2
# 
# # Cara: a little research suggests that values >1 are pretty rare and since
# # this is basically 1, that seems fine to me
# 
# # Phylogenetic signal lambda : 1.00745
# # logL(lambda) : -731.72
# # LR(lambda=0) : 54.7237 
# # P-value (based on LR test) : 1.38721e-13 
# 
# #tree$tip.label
# 
# contMap(tree,log10(transit_hrs), fsize =0.6)
# # save as MGT_contmap
# # 





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
tree <- drop.tip(tree_ultra2, setdiff(tree_ultra2$tip.label, common_species)) 
MRT_hrs <- MRT_hrs[common_species]

#for MRT
lambda_gs2<-phylosig(tree, MRT_hrs,
                     method="lambda",test=TRUE)
lambda_gs2
### put the 0.89 into the model instead of 1. this suggests the phylogenetic
# signal is a bit weaker for MRT vs MGT which makes sense because the data are fewer
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
dat_MRT <- read.csv(file = paste0(homewd, "/data/revision-1/dat_sum_tot_clean_MRT_R2.csv"), header = T, stringsAsFactors = F )
names(dat_MRT)


# clean dataset with only the columns you need
MRT <- dat_MRT[c("phylo_name", "MRT_hrs", "trial.diet", "re_class", "mass_kg")]
#MRT <- na.omit(MRT$transit_hrs)

tree <- read.tree(paste0(homewd, "/data/revision-1/timetree_names.nwk"))
tree_ultra2 <- force.ultrametric(tree, method = "nnls") #no warnings
is.ultrametric(tree_ultra2)

#make sure they match
# Check if the species names match the tree tips
common_species <- intersect(tree_ultra2$tip.label, MRT$phylo_name)
tree_pruned <- drop.tip(tree_ultra2, setdiff(tree_ultra2$tip.label, common_species))
MRT_pruned <- MRT[MRT$phylo_name %in% common_species, ]

# there are some duplicated taxa with different diets:
MRT_pruned[duplicated(MRT_pruned$phylo_name),]

# check to see that the names match
# Species in common_species but not in tree_pruned$tip.label
# missing_from_tree <- setdiff(common_species, tree_pruned$tip.label)
# # Print the missing species
# missing_from_tree 


species1 <- unique(dat_MRT$phylo_name)
species2 <- unique(MRT_pruned$phylo_name)

# Species in df1 but not in df2
missing_in_df2 <- setdiff(species1, species2)
missing_in_df2 

# [1] "Cercopithecus_mitis"       "Cercopithecus_pogonias"   

#sort(tree$tip.label)

###################### MGT - model 1 ################


# Fit PGLS model
#ML=maximum likelihood

# factor the groups so that rodent is the comparasion
unique(MRT_pruned$re_class)
#[1] "Non-Flying Birds" "Ungulates"        "Primates"         "Rodents"          "Bats" 

length(unique(MRT_pruned$phylo_name)) #58


MRT_pruned$re_class <- factor(MRT_pruned$re_class, levels = c("Rodents","Bats", "Non-Flying Birds", "Primates", "Ungulates"))

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


# now test the effect of flyer on MRT transit without accounting for mass or diet


#with full phylogenetic effects:
gls_model_MRT_1 <- gls(log_MRT_hrs ~ flyer, 
                         data = MRT_pruned, 
                        correlation = cor_phylo_fixed1,
                         method = "ML")

#with no phylogenetic effects:
gls_model_MRT_0 <- gls(log_MRT_hrs ~ flyer, 
                        data = MRT_pruned, 
                        correlation = cor_phylo_fixed0,
                        method = "ML")

#with estimated phylogenetic effects:
gls_model_MRT_2 <- gls(log_MRT_hrs ~ flyer,
                        data = MRT_pruned,
                        correlation = cor_phylo_fixed2,
                        method = "ML")

AIC(gls_model_MRT_0, gls_model_MRT_1, gls_model_MRT_2)
# The model with full brownian phylogenetic effects is the best fit.
# df        AIC
# gls_model_MRT_0  3  -42.67912
# gls_model_MRT_1  3 -111.97139 ****
# gls_model_MRT_2  3  -92.43836


summary(gls_model_MRT_1) 
# Coefficients:
#   Value Std.Error   t-value p-value
# (Intercept)  1.637944 0.3066152  5.342018  0.0000
# flyer1      -1.435929 0.5262160 -2.728782  0.0084

# Flyer is a significant term which says that flight is significantly
# negatively related to MRT transit time.
# This is on top of a model that already includes significant 
# phylogenetic clustering in transit time, suggesting that flight affects transit
# time independent of phylogeny. 



# now, let's also account for diet
# first try as a fixed effect:
MRT_pruned$trial.diet <- as.factor(MRT_pruned$trial.diet)

#with full phylogenetic effects:
gls_model_MRT_1_diet <- gls(log_MRT_hrs ~ flyer + trial.diet, 
                        data = MRT_pruned, 
                        correlation = cor_phylo_fixed1,
                        method = "ML")

#with no phylogenetic effects:
gls_model_MRT_0_diet <- gls(log_MRT_hrs ~ flyer + trial.diet, 
                        data = MRT_pruned, 
                        correlation = cor_phylo_fixed0,
                        method = "ML")

#  with estimated phylogenetic effects:
gls_model_MRT_2_diet <- gls(log_MRT_hrs ~ flyer + trial.diet,
                        data = MRT_pruned,
                        correlation = cor_phylo_fixed2,
                        method = "ML")

AIC(gls_model_MRT_0, gls_model_MRT_1, gls_model_MRT_2, 
    gls_model_MRT_0_diet, gls_model_MRT_1_diet, gls_model_MRT_2_diet)
# df       AIC
# gls_model_MRT_0       3  -42.67912
# gls_model_MRT_1       3 -111.97139 *****
# gls_model_MRT_2       3  -92.43836
# gls_model_MRT_0_diet  5  -40.03282
# gls_model_MRT_1_diet  5 -109.59129 *****
# gls_model_MRT_2_diet  5  -89.55585


summary(gls_model_MRT_1_diet)
# Coefficients:
#   Value Std.Error   t-value p-value
# (Intercept)                    1.5559643 0.3199279  4.863484  0.0000
# flyer1                        -1.3399004 0.5580493 -2.401043  0.0197
# trial.dietfruit/nectar/pollen  0.0522815 0.0554737  0.942455  0.3500
# trial.dietprotein              0.2320196 0.2713630  0.855016  0.3962

#  model gls_model_MRT_1_diet has a close AIC to the model without diet
# Plus diet doesn't come out as significant when you do the summary output. 
# So prefer the simpler model (see above for the summary output for gls_model_MRT_1)


# let's try it also as a random effect - here we use the lme function which
# is also supported in the nlme package. we can compare these since all were fit via ML:

gls_model_MRT_1_diet_re <- lme(log_MRT_hrs ~ flyer,
                            random = ~1 | trial.diet,
                            correlation = cor_phylo_fixed1,
                            data = MRT_pruned,
                            method = "ML")

gls_model_MRT_0_diet_re <- lme(log_MRT_hrs ~ flyer,
                                random = ~1 | trial.diet,
                                correlation = cor_phylo_fixed0,
                                data = MRT_pruned,
                                method = "ML")

gls_model_MRT_2_diet_re <- lme(log_MRT_hrs ~ flyer,
                                random = ~1 | trial.diet,
                                correlation = cor_phylo_fixed2,
                                data = MRT_pruned,
                                method = "ML")

AIC(gls_model_MRT_0, gls_model_MRT_0_diet, gls_model_MRT_0_diet_re,
    gls_model_MRT_1, gls_model_MRT_1_diet, gls_model_MRT_1_diet_re,
    gls_model_MRT_2, gls_model_MRT_2_diet, gls_model_MRT_2_diet_re)

# df        AIC
# gls_model_MRT_0          3  -42.67912
# gls_model_MRT_0_diet     5  -40.03282
# gls_model_MRT_0_diet_re  4   33.93345
# gls_model_MRT_1          3 -111.97139 ***
# gls_model_MRT_1_diet     5 -109.59129 **
# gls_model_MRT_1_diet_re  4   20.79840
# gls_model_MRT_2          3  -92.43836
# gls_model_MRT_2_diet     5  -89.55585
# gls_model_MRT_2_diet_re  4   22.04532

# in all cases, the fit with diet as a fixed effect was better supported than the fit as
# a random effect. The models with the lambda=1, which incorporate phylogenetic clustering
# of MRT transit time, are best supported. The top fit model, gls_model_MRT_1,
# offers the best fit overall. When you look at the results and compare:
summary(gls_model_MRT_1_diet) #flyer p=0.0197
summary(gls_model_MRT_1)# flyer p=0.0084

# you see that the significance of "flyer" is declines when diet is also modeled
# none of the trial.diets end up being significant



###################
########## model gls_model_MRT_1 is best for reporting.




# now, let's also consider the effects of mass and, later, mass and diet
 
#with full phylogenetic effects:
gls_model_MRT_1_mass <- gls(log_MRT_hrs ~ log10(mass_kg)*flyer, 
                        data = MRT_pruned, 
                        correlation = cor_phylo_fixed1,
                        method = "ML")

#with no phylogenetic effects:
gls_model_MRT_0_mass <- gls(log_MRT_hrs ~ log10(mass_kg)*flyer, 
                        data = MRT_pruned, 
                        correlation = cor_phylo_fixed0,
                        method = "ML")

#with estimated phylogenetic effects:
gls_model_MRT_2_mass <- gls(log_MRT_hrs ~ log10(mass_kg)*flyer,
                        data = MRT_pruned,
                        correlation = cor_phylo_fixed2,
                        method = "ML")

AIC(gls_model_MRT_0_mass, gls_model_MRT_1_mass, gls_model_MRT_2_mass)
# model 1 is the best. This shows that MRT is still influenced by phylogeny, not
# exclusively by mass (e.g. the model with phylogeny AND mass performed better
# than the model with mass alone)

#                      df     AIC
# gls_model_MRT_0_mass  5  -58.67539
# gls_model_MRT_1_mass  5 -110.38487 ***
# gls_model_MRT_2_mass  5  -93.33385

summary(gls_model_MRT_1_mass)

# Coefficients:
#                         Value    Std.Error   t-value  p-value
# (Intercept)            1.6110127 0.3066086  5.254298  0.0000
# log10(mass_kg)         0.0538528 0.0436853  1.232744  0.2228
# flyer1                -1.7505317 0.5900825 -2.966588  0.0044
# log10(mass_kg):flyer1 -0.2889297 0.2699299 -1.070388  0.2890

# flyer is the only significant term here 
# We can test whether it would be best to drop the interaction which 
# might clarify the variables that matter most.


gls_model_MRT_1_mass_simple <- gls(log_MRT_hrs ~ flyer + log10(mass_kg), 
                             data = MRT_pruned, 
                             correlation = cor_phylo_fixed1,
                             method = "ML")

AIC(gls_model_MRT_1_mass, gls_model_MRT_1_mass_simple)

# df       AIC
# gls_model_MRT_1_mass         5 -110.3849
# gls_model_MRT_1_mass_simple  4 -111.1697

#  there is no difference in AIC,really, though the simple model is slightly better
#so we prefer the simplest model 
# (no interaction term), shown here:
summary(gls_model_MRT_1_mass_simple)

# Coefficients:
#                 Value     Std.Error   t-value p-value
# (Intercept)     1.6147972 0.3069798  5.260271  0.0000
# flyer1         -1.4631818 0.5261519 -2.780912  0.0073
# log10(mass_kg)  0.0462851 0.0431644  1.072298  0.2881

# flyer is significantly negatively related to transit suggesting 
# that flyers have rapid MRT independent of small body size; 
# however, there is no effect of mass alone on transit

# additionally, note that the significance of flyer improved
# in the simple model which I think clarifies the relationship.

# now, let's see if adding diet improves this model:


# I'll just illustrate with the best fit model
# first, as a fixed effect
gls_model_MRT_1_mass_simple_diet <- gls(log_MRT_hrs ~ flyer + log10(mass_kg) +
                                          trial.diet, 
                                    data = MRT_pruned, 
                                    correlation = cor_phylo_fixed1,
                                    method = "ML")

# then, as a random effect

gls_model_MRT_1_mass_simple_diet_re <- lme(log_MRT_hrs ~ flyer + log10(mass_kg),
                                            random = ~1 | trial.diet,
                                            data = MRT_pruned, 
                                            correlation = cor_phylo_fixed1,
                                            method = "ML")

# then, compare
AIC(gls_model_MRT_1_mass_simple, gls_model_MRT_1_mass_simple_diet, gls_model_MRT_1_mass_simple_diet_re)

# df        AIC
# gls_model_MRT_1_mass_simple          4 -111.16969 ***
# gls_model_MRT_1_mass_simple_diet     6 -109.17679
# gls_model_MRT_1_mass_simple_diet_re  5   15.35875

# here, Adding diet does nothing really for the fit. So only have mass in the model is still the best

summary(gls_model_MRT_1_mass_simple)

# Value Std.Error   t-value p-value
# (Intercept)     1.6147972 0.3069798  5.260271  0.0000
# flyer1         -1.4631818 0.5261519 -2.780912  0.0073
# log10(mass_kg)  0.0462851 0.0431644  1.072298  0.2881

#  flyer is still significantly negatively related to MRT, 
# after accounting for mass (and phylogeny). 


#  what if we compare our two best fit models so far? 
AIC(gls_model_MRT_1, gls_model_MRT_1_mass_simple)

# df       AIC
# gls_model_MRT_1              3 -111.9714
# gls_model_MRT_1_mass_simple  4 -111.1697

# These are pretty much the same

# This is pretty surprising since we know from the literature that mass is
# strongly correlated with digesta passage rates, so we might want to also test if this
# is all due to phylogeny

gls_model_MRT_1_mass_only <- gls(log_MRT_hrs ~ log10(mass_kg), 
                                   data = MRT_pruned, 
                                   correlation = cor_phylo_fixed1,
                                   method = "ML")

gls_model_MRT_0_mass_only <- gls(log_MRT_hrs ~ log10(mass_kg), 
                                 data = MRT_pruned, 
                                 correlation = cor_phylo_fixed0,
                                 method = "ML")


gls_model_MRT_2_mass_only <- gls(log_MRT_hrs ~ log10(mass_kg), 
                                 data = MRT_pruned, 
                                 correlation = cor_phylo_fixed2,
                                 method = "ML")


summary(gls_model_MRT_0_mass_only) 
#mass is significant and positively related to MRT
summary(gls_model_MRT_1_mass_only)
#mass is not significantly related to MRT when accounting for phylogeny too
summary(gls_model_MRT_2_mass_only)
#mass is not significantly related to MRT when accounting for phylogeny too

AIC(gls_model_MRT_0_mass_only, gls_model_MRT_1_mass_only, gls_model_MRT_2_mass_only)
# 
# df        AIC
# gls_model_MRT_0_mass_only  3  -25.15709
# gls_model_MRT_1_mass_only  3 -105.53606***
# gls_model_MRT_2_mass_only  3  -82.03153






#table 1 final counts used in the analysis
length(unique(MRT_pruned$phylo_name)) #58 unique species for MRT
(tt.pivot <-ddply(MRT_pruned, .(re_class), summarise, N_species_transit = length(unique(phylo_name))))
# re_class N_species_transit
# 1          Rodents                14
# 2             Bats                12
# 3 Non-Flying Birds                 3
# 4         Primates                17
# 5        Ungulates                12





###### ------ PLOTTING FOR FIGURE 2
library(nlme)
library(ggplot2)
library(dplyr)

# --- Colors and shapes ---
colz <- c(
  "Bats" = "#C49A00",
  "Non-Flying\nBirds" = "#edf8b1",
  "Rodents" = "navy",
  "Primates" = "#00B6EB",
  "Ungulates" = "#A58AFF"
)

shapez <- c(
  "protein" = 21,
  "fruit/nectar/pollen" = 22,
  "fiber/foliage" = 23,
  "insectivore" = 24,
  "omnivore" = 25
)

# --- Fix factor levels ---
# Convert to character to safely modify
MRT_pruned$re_class <- as.character(MRT_pruned$re_class)

# Replace the value
MRT_pruned$re_class[MRT_pruned$re_class == "Non-Flying Birds"] <- "Non-Flying\nBirds"

MRT_pruned$re_class <- factor(MRT_pruned$re_class, levels = c("Bats", "Non-Flying\nBirds", "Rodents", "Primates", "Ungulates"))

# --- Plot ---
p1 <- ggplot(MRT_pruned) +
  geom_boxplot(aes(x = re_class, y = log10(MRT_hrs), fill = re_class), show.legend = FALSE) +
  geom_point(aes(x = re_class, y = log10(MRT_hrs), shape = trial.diet),
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
             size = 3) +
  scale_fill_manual(values = colz) +
  scale_shape_manual(values = shapez, name = "Trial diet") +
  ylab(expression(Log[10] ~ "Mean retention time (hrs)")) +
  coord_cartesian(ylim = c(-0.8, 2.95)) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid = element_blank(),
    legend.position = c(0.85, 0.15),
    legend.background = element_rect(color = "black"),
    plot.margin = unit(c(0.2, 0.2, 0.8, 0.2), "cm")
  )

print(p1)


ggsave(file = paste0(homewd,"/figures_R2/Fig_S1A.png"),
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
  geom_line(data = newdat,
            aes(x = log_mass_c, y = fit, color = flyer),
            linewidth = 1.2) +
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
ggsave("figures_r2/MRT_plot2b.png", plot = p2, width = 7, height = 6, dpi = 300)




##### cowplot the figures
out.plot2 <- cowplot::plot_grid(p1, p2, nrow=1, ncol=2, labels=c("(A)", "(B)"), rel_widths = c(1.2,1.2))

print(out.plot2)

ggsave(file = paste0(homewd,"/figures_R2/FigS1_TwoPanel_R1.jpeg"),
       units="mm",  
       width=170, 
       height=70, 
       scale=3, 
       dpi=300)
