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
homewd <- "/Users/carabrook/Developer/git-transit-time/"
#homewd <- "/Users/emilyruhs/Desktop/UChi_Brook_Lab/GitHub_repos/git-transit-time/"
#homewd <- "/Users/katherinemcferrin/Developer/git-transit-time/"

#load the GIT transit data with phylo name
dat_GIT <- read.csv(file = paste0(homewd, "data/revision-1/dat_sum_tot_clean_GIT.csv"), header = T, stringsAsFactors = F )
dat_MRT <- read.csv(file = paste0(homewd, "data/revision-1/dat_sum_tot_clean_MRT.csv"), header = T, stringsAsFactors = F )


#### species counts ####
#table 1 in the manuscript

length(unique(dat_GIT$genus.species)) #116 unique species for transit time
(tt.pivot <-ddply(dat_GIT, .(re_class), summarise, N_species_transit = length(unique(genus.species))))
# re_class N_species_transit
# 1         Bats                36
# 2   Carnivores                 6
# 3 Flying Birds                14
# 4     Primates                21
# 5     Reptiles                19
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



#load the tree
tree <- read.tree(file = paste0(homewd, "data/revision-1/timetree_names.nwk"))



##### plot the phylogeny for GIT and MRT #######

# Plot the tree with enhancements
plot(tree,
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
common_species <- intersect(tree$tip.label, names(transit_hrs))
tree <- drop.tip(tree, setdiff(tree$tip.label, common_species))
transit_hrs <- transit_hrs[common_species]


lambda_gs1<-phylosig(tree, transit_hrs,
                    method="lambda",test=TRUE)
lambda_gs1 
# Cara: my understanding is that this is the number (0.38) you'd want 
# to put in your phylogenetic model (rather than 1) but I could be wrong here
# lambda_gs1
# Phylogenetic signal lambda : 0.380229 
# logL(lambda) : -775.403 
# LR(lambda=0) : 19.8794 
# P-value (based on LR test) : 8.24865e-06 

#tree$tip.label

contMap(tree,log10(transit_hrs), fsize =0.6)
# save as MGT_contmap
# 





################### MRT phylogenetic signal ##########################
dat_MRT$phylo_name <- as.character(dat_MRT$phylo_name)

tree <- read.tree(file = paste0(homewd, "data/revision-1/timetree_names.nwk"))

#for MRT
MRT <- dat_MRT[c("phylo_name", "MRT_hrs")]
#MRT <- na.omit(MRT)

#make sure they match
# If you have a data frame GIT with the trait
MRT_hrs <- setNames(MRT$MRT_hrs, MRT$phylo_name)

# Then only keep the overlapping species
common_species <- intersect(tree$tip.label, names(MRT_hrs))
tree <- drop.tip(tree, setdiff(tree$tip.label, common_species))
MRT_hrs <- MRT_hrs[common_species]

#for MRT
lambda_gs2<-phylosig(tree, MRT_hrs,
                     method="lambda",test=TRUE)
lambda_gs2
# 
# Phylogenetic signal lambda : 0.90359 
# logL(lambda) : -227.842 
# LR(lambda=0) : 24.439 
# P-value (based on LR test) : 7.66987e-07 



############# plotting phylo signal by tree/phenogram ############


contMap(tree,log10(MRT_hrs), fsize =0.6)


#phenogram(tree, log10(MRT_hrs), spread.labels = T, fsize=0.6, colors = "black")


#########################################  MGT. #################################

########### MGT data cleaning #######

#_________________________________________________________________________________________________________________
#I think this code runs PGLS and also provides lambda estimates as part of the model



# Load your data
dat_GIT <- read.csv(file = paste0(homewd, "/data/revision-1/dat_sum_tot_clean_GIT.csv"), header = T, stringsAsFactors = F )
names(dat_GIT)

# clean dataset with only the columns you need
GIT <- dat_GIT[c("phylo_name", "transit_hrs", "typical.diet", "re_class", "mass_kg")]
#GIT <- na.omit(GIT$transit_hrs)

tree <- read.tree(paste0(homewd, "/data/revision-1/timetree_names.nwk"))


#make sure they match
# Check if the species names match the tree tips
common_species <- intersect(tree$tip.label, GIT$phylo_name)
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, common_species))
GIT_pruned <- GIT[GIT$phylo_name %in% common_species, ]

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

# [1] "Eptesicus_innoxius"        "Acanthagenys_rufogularis"  "Cercopithecus_albogularis"
# [4] "Cercopithecus_mitis"       "Cercopithecus_pogonias"   

#sort(tree$tip.label)

###################### MGT - model 1 ################


# Fit PGLS model
#ML=maximum likelihood

# factor the groups so that rodent is the comparasion
unique(GIT_pruned$re_class)
# "Primates"     "Bats"         "Reptiles"     "Ungulates"    "Flying Birds" "Rodents"      "Carnivores" 

#Cara: keeping this here for plotting but not going to test since we already include
#phylogeny in the analysis
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
# this is what was estimated from your data above:
cor_phylo_fixed2 <- corPagel(lambda_gs1$lambda, phy = tree_pruned, fixed = TRUE, form = ~phylo_name) #assumes no phylogenetic signal


# Cara: now test the effect of flyer on GIT transit without accounting for mass or diet
# All text below is from Cara

#with full phylogenetic effects:
pgls_model_GIT_1 <- gls(log_transit_hrs ~ flyer, 
                         data = GIT_pruned, 
                        correlation = cor_phylo_fixed1,
                         method = "ML")

#with no phylogenetic effects:
pgls_model_GIT_0 <- gls(log_transit_hrs ~ flyer, 
                        data = GIT_pruned, 
                        correlation = cor_phylo_fixed0,
                        method = "ML")

#with estimated phylogenetic effects:
pgls_model_GIT_2 <- gls(log_transit_hrs ~ flyer, 
                        data = GIT_pruned, 
                        correlation = cor_phylo_fixed2,
                        method = "ML")

AIC(pgls_model_GIT_0, pgls_model_GIT_1, pgls_model_GIT_2)
# Cara: both the full phylogenetic effects and the estimated phylogenetic effects
# provide a better fit than the no-phylo effects model. unsurprisingly, the one
# with the data-estimated lambda (model 2) is the best

summary(pgls_model_GIT_2) 
# Cara: flyer is a highly significant term which says that flight is significantly
# negatively related to GIT transit time irregardless of phylogeny.

summary(pgls_model_GIT_0)
# Cara: unsurprisngly, if you don't account for phylogeny, then the effects of flight 
# are even stronger.

summary(pgls_model_GIT_1)
# Cara: and if you over-account for phylogeny, flight is still significant but slightly
# less so.

# now, let's also account for diet
# first try as a fixed effect:
GIT_pruned$typical.diet <- as.factor(GIT_pruned$typical.diet)

# Cara: with full phylogenetic effects:
pgls_model_GIT_1_diet <- gls(log_transit_hrs ~ flyer + typical.diet, 
                        data = GIT_pruned, 
                        correlation = cor_phylo_fixed1,
                        method = "ML")

# Cara: with no phylogenetic effects:
pgls_model_GIT_0_diet <- gls(log_transit_hrs ~ flyer + typical.diet, 
                        data = GIT_pruned, 
                        correlation = cor_phylo_fixed0,
                        method = "ML")

# Cara: with estimated phylogenetic effects:
pgls_model_GIT_2_diet <- gls(log_transit_hrs ~ flyer + typical.diet, 
                        data = GIT_pruned, 
                        correlation = cor_phylo_fixed2,
                        method = "ML")

AIC(pgls_model_GIT_0, pgls_model_GIT_0_diet, pgls_model_GIT_1, pgls_model_GIT_1_diet, pgls_model_GIT_2, pgls_model_GIT_2_diet)
# Cara: interestingly, diet improves the fit if phylogeny is ignored (probably because
# diet aligns with phylogeny) but it detracts in all other cases when modeled as 
# a fixed effect

# Cara: let's try it also as a random effect - here we use the lme function which
# is also supported in the nlme package. we can compare these since all were fit via ML:

pgls_model_GIT_1_diet_re <- lme(log_transit_hrs ~ flyer,
                            random = ~1 | typical.diet,
                            correlation = cor_phylo_fixed1,
                            data = GIT_pruned,
                            method = "ML")

pgls_model_GIT_0_diet_re <- lme(log_transit_hrs ~ flyer,
                                random = ~1 | typical.diet,
                                correlation = cor_phylo_fixed0,
                                data = GIT_pruned,
                                method = "ML")

pgls_model_GIT_2_diet_re <- lme(log_transit_hrs ~ flyer,
                                random = ~1 | typical.diet,
                                correlation = cor_phylo_fixed2,
                                data = GIT_pruned,
                                method = "ML")

AIC(pgls_model_GIT_0, pgls_model_GIT_0_diet, pgls_model_GIT_0_diet_re,
    pgls_model_GIT_1, pgls_model_GIT_1_diet, pgls_model_GIT_1_diet_re,
    pgls_model_GIT_2, pgls_model_GIT_2_diet, pgls_model_GIT_2_diet_re)

#Cara: in all cases, the fit as a fixed effect was more better supported than the fit as
# a random effect. however, still, the addition of diet only improved the models without
# any phylogenetic structure. We conclude this does not matter for this analysis.

# model pgls_model_GIT_2 is best for reporting.

# now, let's also consider the effects of mass and, later, mass and diet
GIT_pruned <- na.omit(GIT_pruned)
#two species don't have masses and are removed
 
#with full phylogenetic effects:
pgls_model_GIT_1_mass <- gls(log_transit_hrs ~ log10(mass_kg)*flyer, 
                        data = GIT_pruned, 
                        correlation = cor_phylo_fixed1,
                        method = "ML")

#with no phylogenetic effects:
pgls_model_GIT_0_mass <- gls(log_transit_hrs ~ log10(mass_kg)*flyer, 
                        data = GIT_pruned, 
                        correlation = cor_phylo_fixed0,
                        method = "ML")

#with estimated phylogenetic effects:
pgls_model_GIT_2_mass <- gls(log_transit_hrs ~ log10(mass_kg)*flyer, 
                        data = GIT_pruned, 
                        correlation = cor_phylo_fixed2,
                        method = "ML")

AIC(pgls_model_GIT_0_mass, pgls_model_GIT_1_mass, pgls_model_GIT_2_mass)
#Cara: again, model 2 is best.

#Cara: let's compare with a version with no mass - 
# we have to refit because the dataset is smaller

pgls_model_GIT_1_slim <- gls(log_transit_hrs ~ flyer, 
                             data = GIT_pruned, 
                             correlation = cor_phylo_fixed1,
                             method = "ML")

#with no phylogenetic effects:
pgls_model_GIT_0_slim <- gls(log_transit_hrs ~ flyer, 
                             data = GIT_pruned, 
                             correlation = cor_phylo_fixed0,
                             method = "ML")

#with estimated phylogenetic effects:
pgls_model_GIT_2_slim <- gls(log_transit_hrs ~ flyer, 
                             data = GIT_pruned, 
                             correlation = cor_phylo_fixed2,
                             method = "ML")

AIC(pgls_model_GIT_0_slim, pgls_model_GIT_0_mass, pgls_model_GIT_1_slim, 
    pgls_model_GIT_1_mass, pgls_model_GIT_2_slim, pgls_model_GIT_2_mass)

# Cara: including mass improves all versions of the models. best version is 
# with fitted lambda value (model 2)

summary(pgls_model_GIT_2_mass)

# Cara: mass is significantly positively related to transit
# but flyer is significantly negatively related to transit suggesting 
# that flyers have rapid transit independent of small body size
# the interaction of flyer and mass is not significant... we can test
# whether it would be best to drop that term

pgls_model_GIT_2_mass_simple <- gls(log_transit_hrs ~ flyer + log10(mass_kg), 
                             data = GIT_pruned, 
                             correlation = cor_phylo_fixed2,
                             method = "ML")

AIC(pgls_model_GIT_2_mass, pgls_model_GIT_2_mass_simple)
# there is no difference in AIC, so we prefer the simplest model 
# (no interaction term), shown here:
summary(pgls_model_GIT_2_mass_simple)

# as above, mass is significantly positively related to transit
# but flyer is significantly negatively related to transit suggesting 
# that flyers have rapid transit independent of small body size


# now, let's see if adding diet improves this model
# you could test the addition of diet to all three here but I would
# suggest dropping the interaction term of mass and flyer in the others
# as well to do that

# I'll just illustrate with the best fit model
# first, as a fixed effect
pgls_model_GIT_2_mass_simple_diet <- gls(log_transit_hrs ~ flyer + log10(mass_kg) + typical.diet, 
                                    data = GIT_pruned, 
                                    correlation = cor_phylo_fixed2,
                                    method = "ML")

# then, as a random effect

pgls_model_GIT_2_mass_simple_diet_re <- lme(log_transit_hrs ~ flyer + log10(mass_kg),
                                            random = ~1 | typical.diet,
                                            data = GIT_pruned, 
                                            correlation = cor_phylo_fixed2,
                                            method = "ML")

# then, compare
AIC(pgls_model_GIT_2_mass_simple, pgls_model_GIT_2_mass_simple_diet, pgls_model_GIT_2_mass_simple_diet_re)

# here, diet does improve model fit when modeled as a fixed effect!
summary(pgls_model_GIT_2_mass_simple_diet)

# Cara: flyer is still significantly negatively related to transit, 
# after accounting for mass and diet (and phylogeny). Mass is still 
# positively related to transit, after accounting for flyer and diet
# (and phylogeny). None of the individual effects of diet are significant, 
# but their inclusion does improve the overall model, so it is best to keep them.
# the frugivore/nectarivore diet has a near-significant negative relationship
# with transit (independent of flying status and mass) which makes intuitive sense.

# I would probably re-do this with experimental diet rather than typical diet
# and of course redo for MRT as well.

# and make sure to convert the tree to ultrametric first.
