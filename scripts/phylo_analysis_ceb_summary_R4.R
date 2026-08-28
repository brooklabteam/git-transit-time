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
homewd <- "/Users/katherinemcferrin/Developer/git-transit-time/"
#homewd <- "/Users/gavindehnert/Desktop/GitHub_repos/git-transit-time"


setwd(homewd)

#load the GIT transit data with phylo name
dat_GIT <- read.csv(file = paste0(homewd, "/data/dat_sum_tot_clean_GIT_R3.csv"), header = T, stringsAsFactors = F )
length(unique(dat_GIT$phylo_name)) #245
dat_MRT <- read.csv(file = paste0(homewd, "/data/dat_sum_tot_clean_MRT_R3.csv"), header = T, stringsAsFactors = F )
length(unique(dat_MRT$phylo_name)) #242

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
print(no_mass) #16
# [1] "Alopochen aegyptiaca"     "Anser anser"              "Atheris nitscheri"        "Casuaris casuaris"       
# [5] "Dromaius novaehollandiae" "Gulosus aristotelis"      "Macronectes giganteus"    "Martes melampus"         
# [9] "Meleagris gallopavo"      "Mitu salvini"             "Plectopterus gambensis"   "Pygoscelis papua"        
# [13] "Rhea americana"           "Stercorarius skua"        "Struthio camelus"         "Uromacer oxyrhynchus"    

# flightless = c("Struthio camelus", "Rhea americana", "Dromaius novaehollandiae",
#                 "Casuaris casuaris", "Pygoscelis papua")

# remove them from both datasets
dat_GIT <- dat_GIT[!dat_GIT[[species_col]] %in% no_mass, ]
dat_MRT <- dat_MRT[!dat_MRT[[species_col]] %in% no_mass, ]

#### species counts ####
#table 1 in the manuscript

length(unique(dat_GIT$genus.species)) #247 unique species for transit time
tt.pivot <-ddply(dat_GIT, .(re_class), summarise, N_species_transit = length(unique(genus.species)))
print(tt.pivot)
# re_class N_species_transit
# 1                 Bats                37
# 2           Carnivores                25
# 3  Even-toed Ungulates                17
# 4         Flying Birds                62
# 5           Marsupials                 3
# 6     Non-Flying Birds                 1 #### Need to remove due to lack of mass
# 7   Odd-toed Ungulates                 7
# 8             Primates                60
# 9             Reptiles                17
# 10             Rodents                15
# 11        Shrews/Moles                 3

dat_GIT <- subset(
  dat_GIT,
  re_class != "Non-Flying Birds"
)

dat_GIT <- subset(
  dat_GIT,
  re_class != "Shrews/Moles"
)

dat_GIT <- subset(
  dat_GIT,
  re_class != "Marsupials"
)

#post filtering counts
length(unique(dat_GIT$genus.species)) #240 unique species for transit time
tt.pivot <-ddply(dat_GIT, .(re_class), summarise, N_species_transit = length(unique(genus.species)))
print(tt.pivot)

#table 1 for GIT
# re_class N_species_transit
# 1                Bats                37
# 2          Carnivores                25
# 3 Even-toed Ungulates                17
# 4        Flying Birds                62
# 5  Odd-toed Ungulates                 7
# 6            Primates                60
# 7            Reptiles                17
# 8             Rodents                15


length(unique(dat_MRT$genus.species)) #238 unique species for mean retention time (minus 2 non-flying bird species)
(mrt.pivot <-ddply(dat_MRT, .(re_class), summarise, N_species_mrt = length(unique(genus.species))))
# re_class N_species_mrt
# 1                Bats            13
# 2          Carnivores            16
# 3        Flying Birds            61
# 4          Marsupials            14
# 6  Odd-toed Ungulates            10
# 7            Primates            40
# 8            Reptiles             3
# 9             Rodents            27
# 10          Ungulates            54 #these should be called even toed ungulates!

dat_MRT <- subset(
  dat_MRT,
  re_class != "Non-Flying Birds"
)


#Counts for the methods section 

#counts of unique species across both TT and MRT
unique_species <- c(unique(dat_MRT$genus.species),unique(dat_GIT$genus.species))
unique_species <- unique(unique_species) #380

#counts of species with unknown diets
unknown_diet_TT <- dat_GIT %>% filter(trial.diet == "unknown")
unknown_diet_MRT <- dat_GIT %>% filter(trial.diet == "unknown")
unknown_diet <- bind_rows(unknown_diet_TT, unknown_diet_MRT)
length(unique(unknown_diet$genus.species)) #21






#load the tree - pulled from timetree so should be in principle an ultrametric tree
tree <- read.tree(file = paste0(homewd, "/data/Book3.nwk"))
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



 
##################  GIT phylogenetic signal ##########
dat_GIT$phylo_name <- as.character(dat_GIT$phylo_name)

#for GIT
GIT <- dat_GIT[c("phylo_name", "transit_hrs")] #275
#GIT <- na.omit(GIT)

length(unique(GIT$phylo_name)) #231
#three don't have matches in the tree

#make sure they match
# If you have a data frame GIT with the trait
transit_hrs <- setNames(GIT$transit_hrs, GIT$phylo_name) #275

# Then only keep the overlapping species
common_species <- intersect(tree_ultra2$tip.label, names(transit_hrs)) #228

########
#### What is missing and not matching?
# in your data but NOT on the tree
in_data_not_tree <- setdiff(names(transit_hrs), tree_ultra2$tip.label)

# on the tree but NOT in your data
in_tree_not_data <- setdiff(tree_ultra2$tip.label, names(transit_hrs))

cat("In data but not on tree (", length(in_data_not_tree), "):\n", sep = ""); print(sort(in_data_not_tree))
# In data but not on tree (3):
#   [1] ""                "Cercopithecinae" "Pycnonotidae"  
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

# Cara: a little research suggests that values >1 are pretty rare and since
# this is basically 1, that seems fine to me

#### NEW VALUES R3
# Phylogenetic signal lambda : 1.00164 
# logL(lambda) : -1460.12 
# LR(lambda=0) : 151.493 
# P-value (based on LR test) : 8.17584e-35 

#tree$tip.label

contMap(tree,log10(transit_hrs), fsize =0.5)
# save as MGT_contmap
# 





################### MRT phylogenetic signal ##########################
dat_MRT$phylo_name <- as.character(dat_MRT$phylo_name)

tree <- read.tree(file = paste0(homewd, "/data/Book3.nwk"))
tree_ultra2 <- force.ultrametric(tree, method = "nnls") #no warnings
is.ultrametric(tree_ultra2)

#for MRT
MRT <- dat_MRT[c("phylo_name", "MRT_hrs")] #292
#MRT <- na.omit(MRT)

length(unique(MRT$phylo_name)) #227

#make sure theMRT#make sure they match
# If you have a data frame GIT with the trait
MRT_hrs <- setNames(MRT$MRT_hrs, MRT$phylo_name) #292

# Then only keep the overlapping species
common_species <- intersect(tree_ultra2$tip.label, names(MRT_hrs)) #225
tree <- drop.tip(tree_ultra2, setdiff(tree_ultra2$tip.label, common_species)) 
MRT_hrs <- MRT_hrs[common_species]

#for MRT
lambda_gs2<-phylosig(tree, MRT_hrs,
                     method="lambda",test=TRUE)
lambda_gs2
# Phylogenetic signal lambda : 0.943984 
# logL(lambda) : -982.221 
# LR(lambda=0) : 148.399 
# P-value (based on LR test) : 3.88009e-34 



############# plotting phylo signal by tree/phenogram ############


contMap(tree,log10(MRT_hrs), fsize =0.5)


#phenogram(tree, log10(MRT_hrs), spread.labels = T, fsize=0.6, colors = "black")





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

tree <- read.tree(paste0(homewd, "/data/Book3.nwk"))
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

#[1] ""                "Pycnonotidae"    "Cercopithecinae"    # we knew that

#sort(tree$tip.label)

###################### MGT - model 1 ################


# Fit PGLS model
#ML=maximum likelihood

# factor the groups so that rodent is the comparasion
unique(GIT_pruned$re_class)
# [1] "Flying Birds"        "Carnivores"          "Even-toed Ungulates" "Primates"            "Bats"               
# [6] "Reptiles"            "Rodents"             "Odd-toed Ungulates" 

length(unique(GIT_pruned$phylo_name)) #228


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

# the model above is already 1, so no need to run this
#cor_phylo_fixed2 <- corPagel(lambda_gs1$lambda, phy = tree_pruned, fixed = TRUE, form = ~phylo_name) 


# Cara: now test the effect of flyer on GIT transit without accounting for mass or diet
# All text below is from Cara

#with full phylogenetic effects:
gls_model_GIT_1 <- gls(log_transit_hrs ~ flyer, 
                         data = GIT_pruned, 
                        correlation = cor_phylo_fixed1,
                         method = "ML")

summary(gls_model_GIT_1)$tTable
# Value Std.Error   t-value      p-value
# (Intercept)  1.752164 0.2548899  6.874198 4.340035e-11
# flyer1      -1.230929 0.1935776 -6.358841 8.681753e-10


# removing correlation structure, as its the same and the lambda=0 breaks the matrix
# with no phylogenetic effects:
gls_model_GIT_0 <- gls(log_transit_hrs ~ flyer, 
                        data = GIT_pruned, 
                        #correlation = cor_phylo_fixed0,
                        method = "ML")
summary(gls_model_GIT_0)$tTable
# Value  Std.Error   t-value      p-value
# (Intercept)  1.215590 0.04311631  28.19329 2.723093e-82
# flyer1      -1.349293 0.06798497 -19.84694 1.304523e-54

#with estimated phylogenetic effects:
# gls_model_GIT_2 <- gls(log_transit_hrs ~ flyer,
#                         data = GIT_pruned,
#                         correlation = cor_phylo_fixed2,
#                         method = "ML")
# summary(gls_model_GIT_2)$tTable

AIC(gls_model_GIT_0, gls_model_GIT_1)
# df       AIC
# gls_model_GIT_0  3   447.827
# gls_model_GIT_1  3 -1856.137


# Emily: flyer is a significant term which says that flight is significantly
# negatively related to GIT transit time.
# R3: This is on top of a model that already includes significant 
# phylogenetic clustering in transit time, suggesting that flight affects transit
# time independent of phylogeny. 



# now, let's also account for diet
# first try as a fixed effect:
GIT_pruned$trial.diet <- as.factor(GIT_pruned$trial.diet)

# Cara: with full phylogenetic effects:
gls_model_GIT_1_diet <- gls(log_transit_hrs ~ flyer + trial.diet, 
                        data = GIT_pruned, 
                        correlation = cor_phylo_fixed1,
                        method = "ML")
summary(gls_model_GIT_1_diet)$tTable
# Value  Std.Error    t-value      p-value
# (Intercept)                    1.64063845 0.25226899  6.5035280 3.909114e-10
# flyer1                        -1.22094616 0.19114456 -6.3875537 7.568209e-10
# trial.dietfruit/nectar/pollen -0.05688240 0.06566166 -0.8662955 3.871148e-01
# trial.dietmeat                 0.01666727 0.09496347  0.1755124 8.608116e-01
# trial.dietmixed               -0.22455578 0.23448017 -0.9576749 3.391029e-01
# trial.dietprotein              0.07380284 0.04977917  1.4826050 1.393721e-01
# trial.dietunknown              0.39808656 0.11449717  3.4768245 5.931223e-04

# Cara: with no phylogenetic effects:
gls_model_GIT_0_diet <- gls(log_transit_hrs ~ flyer + trial.diet, 
                        data = GIT_pruned, 
                        #correlation = cor_phylo_fixed0,
                        method = "ML")
summary(gls_model_GIT_0_diet)$tTable
logLik(gls_model_GIT_0_diet) #'log Lik.' -213.9687 (df=8)

# Value  Std.Error     t-value      p-value
# (Intercept)                    1.26798783 0.06581813  19.2650244 3.013578e-52
# flyer1                        -1.37435400 0.07583771 -18.1223031 3.036179e-48
# trial.dietfruit/nectar/pollen -0.15445751 0.09569615  -1.6140410 1.077131e-01
# trial.dietmeat                 0.09359074 0.11372912   0.8229268 4.112923e-01
# trial.dietmixed               -0.36004994 0.31907741  -1.1284094 2.601715e-01
# trial.dietprotein             -0.10897356 0.09693681  -1.1241712 2.619615e-01
# trial.dietunknown              0.24093170 0.14833843   1.6242028 1.055258e-01

gls_model_GIT_0_diet$sigma            # residual SD — 0.533


# Cara: with estimated phylogenetic effects:
# gls_model_GIT_2_diet <- gls(log_transit_hrs ~ flyer + trial.diet,
#                         data = GIT_pruned,
#                         correlation = cor_phylo_fixed2,
#                         method = "ML")
# summary(gls_model_GIT_2_diet)$tTable

AIC(gls_model_GIT_0, gls_model_GIT_1, gls_model_GIT_0_diet, gls_model_GIT_1_diet)
# df        AIC
# gls_model_GIT_0       3   447.8270
# gls_model_GIT_1       3 -1856.1368
# gls_model_GIT_0_diet  8   443.9375
# gls_model_GIT_1_diet  8 -1868.2564

summary(gls_model_GIT_1_diet)$tTable
# Value  Std.Error    t-value      p-value
# (Intercept)                    1.64063845 0.25226899  6.5035280 3.909114e-10
# flyer1                        -1.22094616 0.19114456 -6.3875537 7.568209e-10
# trial.dietfruit/nectar/pollen -0.05688240 0.06566166 -0.8662955 3.871148e-01
# trial.dietmeat                 0.01666727 0.09496347  0.1755124 8.608116e-01
# trial.dietmixed               -0.22455578 0.23448017 -0.9576749 3.391029e-01
# trial.dietprotein              0.07380284 0.04977917  1.4826050 1.393721e-01
# trial.dietunknown              0.39808656 0.11449717  3.4768245 5.931223e-04

# Emily: model gls_model_GIT_1_diet has the best fit.


# Cara: let's try it also as a random effect - here we use the lme function which
# is also supported in the nlme package. we can compare these since all were fit via ML:

gls_model_GIT_1_diet_re <- lme(log_transit_hrs ~ flyer,
                            random = ~1 | trial.diet,
                            correlation = cor_phylo_fixed1,
                            data = GIT_pruned,
                            method = "ML")

summary(gls_model_GIT_1_diet_re)$tTable
# Value Std.Error  DF   t-value     p-value
# (Intercept)  1.446842 0.1980906 264  7.303942 3.30579e-12
# flyer1      -1.101615 0.1866870 264 -5.900867 1.10537e-08

gls_model_GIT_0_diet_re <- lme(log_transit_hrs ~ flyer,
                                random = ~1 | trial.diet,
                                #correlation = cor_phylo_fixed0,
                                data = GIT_pruned,
                                method = "ML")
summary(gls_model_GIT_0_diet_re)$tTable
# Value  Std.Error  DF   t-value      p-value
# (Intercept)  1.234750 0.06211718 264  19.87775 2.232168e-54
# flyer1      -1.348107 0.07182450 264 -18.76946 1.624744e-50

# gls_model_GIT_2_diet_re <- lme(log_transit_hrs ~ flyer,
#                                 random = ~1 | trial.diet,
#                                 correlation = cor_phylo_fixed2,
#                                 data = GIT_pruned,
#                                 method = "ML")
# summary(gls_model_GIT_2_diet_re)$tTable

AIC(gls_model_GIT_0, gls_model_GIT_0_diet, gls_model_GIT_0_diet_re,
    gls_model_GIT_1, gls_model_GIT_1_diet, gls_model_GIT_1_diet_re)

# df        AIC
# gls_model_GIT_0          3   447.8270
# gls_model_GIT_0_diet     8   443.9375
# gls_model_GIT_0_diet_re  4   447.5670
# gls_model_GIT_1          3 -1856.1368
# gls_model_GIT_1_diet     8 -1868.2564 **
# gls_model_GIT_1_diet_re  4  -467.5784

#Emily: in all cases, the fit with diet as a fixed effect was better supported than the fit as
# a random effect. 
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

# below won't run because the correlation structure is set by a 
# phylogeny of 109 species but the dataset has 117 entries due to some taxa being 
# reeated multiple times due tp differences in diet or mass reported each time

# Emily: It ran for me...
 
#with full phylogenetic effects:
gls_model_GIT_1_mass <- gls(log_transit_hrs ~ log10(mass_kg)*flyer, 
                        data = GIT_pruned, 
                        correlation = cor_phylo_fixed1,
                        method = "ML")
summary(gls_model_GIT_1_mass)$tTable

# Value  Std.Error   t-value      p-value
# (Intercept)            1.8013166 0.25270159  7.128236 9.502517e-12
# log10(mass_kg)         0.1303174 0.04027735  3.235500 1.367364e-03
# flyer1                -1.1769934 0.22870669 -5.146301 5.155149e-07
# log10(mass_kg):flyer1 -0.1307296 0.07274327 -1.797137 7.344383e-02

#with no phylogenetic effects:
gls_model_GIT_0_mass <- gls(log_transit_hrs ~ log10(mass_kg)*flyer, 
                        data = GIT_pruned, 
                        #correlation = cor_phylo_fixed0,
                        method = "ML")
summary(gls_model_GIT_0_mass)$tTable
# Value  Std.Error    t-value      p-value
# (Intercept)            1.19490442 0.04346391  27.491876 7.680851e-80
# log10(mass_kg)         0.09343323 0.03579736   2.610059 9.563785e-03
# flyer1                -1.32014336 0.12455814 -10.598612 3.882855e-22
# log10(mass_kg):flyer1 -0.08910565 0.06429447  -1.385899 1.669342e-01

#with estimated phylogenetic effects: - again, lambda is already near/at 1
# gls_model_GIT_2_mass <- gls(log_transit_hrs ~ log10(mass_kg)*flyer,
#                         data = GIT_pruned,
#                         correlation = cor_phylo_fixed2,
#                         method = "ML")
# summary(gls_model_GIT_2_mass)$tTable

AIC(gls_model_GIT_0_mass, gls_model_GIT_1_mass)
#Emily: model 1 is the best
#Cara: This shows that transit time is still influenced by phylogeny, not
# exclusively by mass (e.g. the model with phylogeny AND mass performed better
# than the model with mass alone)

# df        AIC
# gls_model_GIT_0_mass  5   444.9928
# gls_model_GIT_1_mass  5 -1862.5617 ***

summary(gls_model_GIT_1_mass)$tTable

# Value  Std.Error   t-value      p-value
# (Intercept)            1.8013166 0.25270159  7.128236 9.502517e-12
# log10(mass_kg)         0.1303174 0.04027735  3.235500 1.367364e-03
# flyer1                -1.1769934 0.22870669 -5.146301 5.155149e-07
# log10(mass_kg):flyer1 -0.1307296 0.07274327 -1.797137 7.344383e-02

# all significant, except interaction

gls_model_GIT_1_mass_simple <- gls(log_transit_hrs ~ flyer + log10(mass_kg), 
                             data = GIT_pruned, 
                             correlation = cor_phylo_fixed1,
                             method = "ML")
summary(gls_model_GIT_1_mass_simple)$tTable

# Value  Std.Error   t-value      p-value
# (Intercept)     1.7482888 0.25201485  6.937245 2.996258e-11
# flyer1         -1.0063641 0.20893150 -4.816718 2.444259e-06
# log10(mass_kg)  0.0891445 0.03326432  2.679883 7.820818e-03

AIC(gls_model_GIT_1_mass, gls_model_GIT_1_mass_simple)

# df       AIC
# gls_model_GIT_1_mass         5 -1862.562
# gls_model_GIT_1_mass_simple  4 -1861.303

# there is no difference in AIC,really, though the complex model is slightly better
#so I'm going to write up the simple model because they are so close and significance is across the board
# (no interaction term), shown here:
summary(gls_model_GIT_1_mass_simple)$tTable

# Value  Std.Error   t-value      p-value
# (Intercept)     1.7482888 0.25201485  6.937245 2.996258e-11
# flyer1         -1.0063641 0.20893150 -4.816718 2.444259e-06
# log10(mass_kg)  0.0891445 0.03326432  2.679883 7.820818e-03

# flyer is significantly negatively related to transit suggesting 
# that flyers have rapid transit independent of small body size; 
# however, there is also significance across the board, and its significant in the complex model too

# now, let's see if adding diet improves this model:


# first, as a fixed effect
gls_model_GIT_1_mass_simple_diet <- gls(log_transit_hrs ~ flyer + log10(mass_kg) + trial.diet, 
                                    data = GIT_pruned, 
                                    correlation = cor_phylo_fixed1,
                                    method = "ML")
summary(gls_model_GIT_1_mass_simple_diet)$tTable
# Value  Std.Error    t-value      p-value
# (Intercept)                    1.63605499 0.24792279  6.5990503 2.268057e-10
# flyer1                        -0.96076368 0.20452372 -4.6975661 4.246778e-06
# log10(mass_kg)                 0.10806051 0.03359509  3.2165568 1.459781e-03
# trial.dietfruit/nectar/pollen -0.03676489 0.06483173 -0.5670817 5.711424e-01
# trial.dietmeat                -0.04141964 0.09505700 -0.4357348 6.633866e-01
# trial.dietmixed               -0.21254253 0.23046690 -0.9222258 3.572560e-01
# trial.dietprotein              0.09623812 0.04941547  1.9475302 5.253528e-02
# trial.dietunknown              0.41918219 0.11271367  3.7190003 2.444767e-04

# then, as a random effect

gls_model_GIT_1_mass_simple_diet_re <- lme(log_transit_hrs ~ flyer + log10(mass_kg),
                                            random = ~1 | trial.diet,
                                            data = GIT_pruned, 
                                            correlation = cor_phylo_fixed1,
                                            method = "ML")
summary(gls_model_GIT_1_mass_simple_diet_re)$tTable
# Value  Std.Error  DF   t-value      p-value
# (Intercept)     1.4373073 0.20422277 263  7.037939 1.693145e-11
# flyer1         -0.8134608 0.19339332 263 -4.206251 3.564002e-05
# log10(mass_kg)  0.1316382 0.03091108 263  4.258608 2.866041e-05

# then, compare
AIC(gls_model_GIT_1_mass_simple, gls_model_GIT_1_mass_simple_diet, gls_model_GIT_1_mass_simple_diet_re)

# df        AIC
# gls_model_GIT_1_mass_simple          4 -1861.3033
# gls_model_GIT_1_mass_simple_diet     9 -1876.7130 **
# gls_model_GIT_1_mass_simple_diet_re  5  -483.2718

# here, mass and diet does improve model fit when diet is modeled as a fixed effect! 
summary(gls_model_GIT_1_mass_simple_diet)$tTable

# Value  Std.Error    t-value      p-value
# (Intercept)                    1.63605499 0.24792279  6.5990503 2.268057e-10
# flyer1                        -0.96076368 0.20452372 -4.6975661 4.246778e-06 *
# log10(mass_kg)                 0.10806051 0.03359509  3.2165568 1.459781e-03 *
# trial.dietfruit/nectar/pollen -0.03676489 0.06483173 -0.5670817 5.711424e-01
# trial.dietmeat                -0.04141964 0.09505700 -0.4357348 6.633866e-01
# trial.dietmixed               -0.21254253 0.23046690 -0.9222258 3.572560e-01
# trial.dietprotein              0.09623812 0.04941547  1.9475302 5.253528e-02 *marginal
# trial.dietunknown              0.41918219 0.11271367  3.7190003 2.444767e-04 *

# Emily: flyer is still significantly negatively related to transit, 
# after accounting for mass and diet (and phylogeny). 
# The protein and unknown diets are significant but nothing else -
# but their inclusion does improve the overall model, so it is best to keep them.
# the meat diet has a significant positive relationship
# with transit (independent of flying status and mass) which makes intuitive sense.

# R3: I'll also note that including diet here again
# improved the significance of flyer but had no effect on mass, supporting my 
# prior hypothesis that there are probably some non-flyers with fruit/nectar/pollen diets
# that had fairly fast transit time and meat eaters with slow transit time

# Cara: what if we compare our two best fit models so far? 
AIC(gls_model_GIT_1_diet, gls_model_GIT_1_mass_simple_diet)

# df       AIC
# gls_model_GIT_1_diet              8 -1868.256
# gls_model_GIT_1_mass_simple_diet  9 -1876.713



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
# (Intercept)    0.8899226 0.04567371 19.48435 2.449235e-53
# log10(mass_kg) 0.3316948 0.02728241 12.15782 2.130761e-27

summary(gls_model_GIT_1_mass_only)$tTable
# mass is hugely significant and positively related to gut transit

# Value  Std.Error  t-value      p-value
# (Intercept)    1.5264122 0.25779088 5.921126 9.725480e-09
# log10(mass_kg) 0.1534062 0.03170416 4.838678 2.204503e-06

AIC(gls_model_GIT_0_mass_only, gls_model_GIT_1_mass_only)
# incorporating phylogeny with mass offers a much better fit to the data however!
# 
# df        AIC
# gls_model_GIT_0_mass_only  3   573.5679
# gls_model_GIT_1_mass_only  3 -1840.8033




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
  #"Non-Flying\nBirds" = "#edf8b1",
  "Rodents" = "navy",
  "Primates" = "#00B6EB",
  "Even-toed\nUngulates" = "#A58AFF",
  "Reptiles" = "#FB61D7",
  "Carnivores" = "#00C094","Marsupials"="#E08B45", "Shrews/Moles"="#c51b8c", "Odd-toed\nUngulates"="thistle1"
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
# [6] "Reptiles"            "Rodents"             "Odd-toed Ungulates"      

# Replace the value
#GIT_pruned$re_class[GIT_pruned$re_class == "Flying Birds"] <- "Flying\nBirds"
#GIT_pruned$re_class[GIT_pruned$re_class == "Non-Flying Birds"] <- "Non-Flying\nBirds"
GIT_pruned$re_class2[GIT_pruned$re_class2 == "Flying Birds"] <- "Flying\nBirds"
GIT_pruned$re_class2[GIT_pruned$re_class2 == "Even-toed Ungulates"] <- "Even-toed\nUngulates"
GIT_pruned$re_class2[GIT_pruned$re_class2 == "Odd-toed Ungulates"] <- "Odd-toed\nUngulates"

GIT_pruned$re_class2 <- factor(GIT_pruned$re_class2, levels = c("Flying\nBirds","Bats", "Rodents",
                                                              "Carnivores",
                                                              "Marsupials", "Lagomorphs",
                                                              "Primates", "Even-toed\nUngulates", "Odd-toed\nUngulates",
                                                              "Shrews/Moles",
                                                              "Reptiles"))

# --- Plot ---
p1 <- ggplot(GIT_pruned) +
  geom_boxplot(aes(x = re_class2, y = log10(transit_hrs), fill = re_class2), show.legend = FALSE) +
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


ggsave(file = paste0(homewd,"/figures_r4/Fig_2A.png"),
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
# Value  Std.Error   t-value p-value
# (Intercept)  1.6899592 0.25307741  6.677637  0.0000
# flyer1      -1.0063641 0.20893150 -4.816718  0.0000
# log_mass_c   0.0891445 0.03326432  2.679883  0.0078



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
ggsave("figures_r4/GIT_plot2b.png", plot = p2, width = 7, height = 6, dpi = 300)




##### cowplot the figures
out.plot2 <- cowplot::plot_grid(p1, p2, nrow=1, ncol=2, labels=c("(A)", "(B)"), rel_widths = c(1.2,1.2))

print(out.plot2)

ggsave(file = paste0(homewd,"/figures_r4/Fig2_TwoPanel_R2.jpeg"),
       units="mm",  
       width=170, 
       height=70, 
       scale=3, 
       dpi=300)
