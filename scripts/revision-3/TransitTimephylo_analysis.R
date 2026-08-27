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
#homewd <- "/Users/katherinemcferrin/Developer/git-transit-time/"
#homewd <- "/Users/gavindehnert/Desktop/GitHub_repos/git-transit-time"


setwd(homewd)

#load the GIT transit data with phylo name
dat_GIT <- read.csv(file = paste0(homewd, "/data/revision-3/dat_sum_tot_clean_GIT_R3.csv"), header = T, stringsAsFactors = F )
length(unique(dat_GIT$phylo_name)) #245
dat_MRT <- read.csv(file = paste0(homewd, "/data/revision-3/dat_sum_tot_clean_MRT_R3.csv"), header = T, stringsAsFactors = F )
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

length(unique(dat_GIT$phylo_name)) #238 (7 removed)
length(unique(dat_MRT$phylo_name)) #229 42 (13 removed)

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

length(unique(dat_MRT$genus.species)) #240 unique species for mean retention time
(mrt.pivot <-ddply(dat_MRT, .(re_class), summarise, N_species_mrt = length(unique(genus.species))))

# re_class N_species_mrt
# 1                Bats            13
# 2          Carnivores            16
# 3        Flying Birds            61
# 4          Marsupials            14
# 5    Non-Flying Birds             2
# 6  Odd-toed Ungulates            10
# 7            Primates            40
# 8            Reptiles             3
# 9             Rodents            27
# 10          Ungulates            54

dat_MRT <- subset(
  dat_MRT,
  re_class != "Non-Flying Birds"
)

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
GIT <- dat_GIT[c("phylo_name", "transit_hrs")] #285
#GIT <- na.omit(GIT)

length(unique(GIT$phylo_name)) #237
#three don't have matches in the tree

#make sure they match
# If you have a data frame GIT with the trait
transit_hrs <- setNames(GIT$transit_hrs, GIT$phylo_name) #285

# Then only keep the overlapping species
common_species <- intersect(tree_ultra2$tip.label, names(transit_hrs)) #234

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






tree <- drop.tip(tree_ultra2, setdiff(tree_ultra2$tip.label, common_species)) #234
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

#### NEW VALUES R3
# Phylogenetic signal lambda : 1.00144 
# logL(lambda) : -1498.2 
# LR(lambda=0) : 150.254 
# P-value (based on LR test) : 1.52597e-34 

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
# [6] "Reptiles"            "Rodents"             "Shrews/Moles"        "Odd-toed Ungulates"  "Marsupials"  

length(unique(GIT_pruned$phylo_name)) #234


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
# (Intercept)  1.707836 0.2606986  6.550999 2.745122e-10
# flyer1      -1.111084 0.1999771 -5.556053 6.439758e-08


# removing correlation structure, as its the same and the lambda=0 breaks the matrix
# with no phylogenetic effects:
gls_model_GIT_0 <- gls(log_transit_hrs ~ flyer, 
                        data = GIT_pruned, 
                        #correlation = cor_phylo_fixed0,
                        method = "ML")
summary(gls_model_GIT_0)$tTable
# Value  Std.Error   t-value      p-value
# (Intercept)  1.137517 0.04681320  24.29907 7.935726e-71
# flyer1      -1.271220 0.07516372 -16.91268 1.180153e-44

#with estimated phylogenetic effects:
# gls_model_GIT_2 <- gls(log_transit_hrs ~ flyer,
#                         data = GIT_pruned,
#                         correlation = cor_phylo_fixed2,
#                         method = "ML")
# summary(gls_model_GIT_2)$tTable

AIC(gls_model_GIT_0, gls_model_GIT_1)
# df        AIC
# gls_model_GIT_0  3   527.2687
# gls_model_GIT_1  3 -1872.1434


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
# Value  Std.Error     t-value      p-value
# (Intercept)                    1.603231829 0.25902368  6.18951836 2.193966e-09
# flyer1                        -1.113159289 0.19911900 -5.59042228 5.472127e-08
# trial.dietfruit/nectar/pollen -0.032738474 0.06910322 -0.47376194 6.360469e-01
# trial.dietmeat                 0.006531261 0.09811557  0.06656702 9.469749e-01
# trial.dietmixed               -0.208655518 0.24678181 -0.84550608 3.985661e-01
# trial.dietprotein              0.073684204 0.05209857  1.41432289 1.584024e-01
# trial.dietunknown              0.410504755 0.12034570  3.41104617 7.446848e-04

# Cara: with no phylogenetic effects:
gls_model_GIT_0_diet <- gls(log_transit_hrs ~ flyer + trial.diet, 
                        data = GIT_pruned, 
                        #correlation = cor_phylo_fixed0,
                        method = "ML")
summary(gls_model_GIT_0_diet)$tTable
logLik(gls_model_GIT_0_diet) #'log Lik.' Inf (df=9)

# Value  Std.Error     t-value      p-value
# (Intercept)                    1.16120970 0.07115383  16.3197065 2.703367e-42
# flyer1                        -1.31830021 0.08506992 -15.4966674 2.506670e-39
# trial.dietfruit/nectar/pollen -0.07850896 0.10604584  -0.7403305 4.597334e-01
# trial.dietmeat                 0.12342959 0.12539922   0.9842931 3.258398e-01
# trial.dietmixed               -0.27195641 0.35942498  -0.7566430 4.499138e-01
# trial.dietprotein             -0.08562395 0.10583852  -0.8090055 4.192139e-01
# trial.dietunknown              0.29699449 0.16612674   1.7877585 7.492001e-02

gls_model_GIT_0_diet$sigma            # residual SD — 0.595


# Cara: with estimated phylogenetic effects:
# gls_model_GIT_2_diet <- gls(log_transit_hrs ~ flyer + trial.diet,
#                         data = GIT_pruned,
#                         correlation = cor_phylo_fixed2,
#                         method = "ML")
# summary(gls_model_GIT_2_diet)$tTable

AIC(gls_model_GIT_0, gls_model_GIT_1, gls_model_GIT_0_diet, gls_model_GIT_1_diet)
# df        AIC
# gls_model_GIT_0       3   527.2687
# gls_model_GIT_1       3 -1872.1434
# gls_model_GIT_0_diet  9   523.1859
# gls_model_GIT_1_diet  9 -1882.9853

summary(gls_model_GIT_1_diet)$tTable
# Value  Std.Error     t-value      p-value
# (Intercept)                    1.603231829 0.25902368  6.18951836 2.193966e-09
# flyer1                        -1.113159289 0.19911900 -5.59042228 5.472127e-08
# trial.dietfruit/nectar/pollen -0.032738474 0.06910322 -0.47376194 6.360469e-01
# trial.dietmeat                 0.006531261 0.09811557  0.06656702 9.469749e-01
# trial.dietmixed               -0.208655518 0.24678181 -0.84550608 3.985661e-01
# trial.dietprotein              0.073684204 0.05209857  1.41432289 1.584024e-01
# trial.dietunknown              0.410504755 0.12034570  3.41104617 7.446848e-04

# Emily: model gls_model_GIT_1_diet has the best fit.


# Cara: let's try it also as a random effect - here we use the lme function which
# is also supported in the nlme package. we can compare these since all were fit via ML:

gls_model_GIT_1_diet_re <- lme(log_transit_hrs ~ flyer,
                            random = ~1 | trial.diet,
                            correlation = cor_phylo_fixed1,
                            data = GIT_pruned,
                            method = "ML")

summary(gls_model_GIT_1_diet_re)$tTable
# Value Std.Error  DF   t-value      p-value
# (Intercept)  1.409541 0.1911813 274  7.372795 1.985164e-12
# flyer1      -1.033279 0.1931381 274 -5.349949 1.859049e-07

gls_model_GIT_0_diet_re <- lme(log_transit_hrs ~ flyer,
                                random = ~1 | trial.diet,
                                #correlation = cor_phylo_fixed0,
                                data = GIT_pruned,
                                method = "ML")
summary(gls_model_GIT_0_diet_re)$tTable
# Value  Std.Error  DF   t-value      p-value
# (Intercept)  1.146042 0.05349968 274  21.42147 1.761078e-60
# flyer1      -1.272709 0.07731297 274 -16.46178 8.306119e-43

# gls_model_GIT_2_diet_re <- lme(log_transit_hrs ~ flyer,
#                                 random = ~1 | trial.diet,
#                                 correlation = cor_phylo_fixed2,
#                                 data = GIT_pruned,
#                                 method = "ML")
# summary(gls_model_GIT_2_diet_re)$tTable

AIC(gls_model_GIT_0, gls_model_GIT_0_diet, gls_model_GIT_0_diet_re,
    gls_model_GIT_1, gls_model_GIT_1_diet, gls_model_GIT_1_diet_re)

# df        AIC
# gls_model_GIT_0          3   527.2687
# gls_model_GIT_0_diet     8   527.8280
# gls_model_GIT_0_diet_re  4   529.1113
# gls_model_GIT_1          3 -1872.1434
# gls_model_GIT_1_diet     8 -1881.2743 ** top
# gls_model_GIT_1_diet_re  4  -617.6787

#Emily: in all cases, the fit with diet as a fixed effect was better supported than the fit as
# a random effect. 
# Cara: The models with the lambda=1, which incorporate phylogenetic clustering
# of GIT transit time, are best supported. The top fit model, gls_model_GIT_1_diet,
# offers the best fit overall. When you look at the results and compare:
summary(gls_model_GIT_1_diet) #flyer p=6.439758e-08
summary(gls_model_GIT_1)# flyer p=1.859049e-07

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
# (Intercept)            1.7879830 0.25699992  6.957134 2.505491e-11
# log10(mass_kg)         0.1520478 0.03964739  3.835001 1.555048e-04
# flyer1                -1.1057437 0.23552817 -4.694741 4.203755e-06
# log10(mass_kg):flyer1 -0.1629897 0.07374201 -2.210269 2.790399e-02

#with no phylogenetic effects:
gls_model_GIT_0_mass <- gls(log_transit_hrs ~ log10(mass_kg)*flyer, 
                        data = GIT_pruned, 
                        #correlation = cor_phylo_fixed0,
                        method = "ML")
summary(gls_model_GIT_0_mass)$tTable
# Value  Std.Error   t-value      p-value
# (Intercept)            1.1160321 0.04572375 24.408147 5.397274e-71
# log10(mass_kg)         0.1595073 0.03668588  4.347920 1.932736e-05
# flyer1                -1.2412710 0.13586039 -9.136372 1.410833e-17
# log10(mass_kg):flyer1 -0.1551797 0.06908062 -2.246356 2.546918e-02

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
# gls_model_GIT_0_mass  5   512.7122
# gls_model_GIT_1_mass  5 -1882.7087 ** mass is way stronger

summary(gls_model_GIT_1_mass)$tTable

# Value  Std.Error   t-value      p-value
# (Intercept)            1.7879830 0.25699992  6.957134 2.505491e-11
# log10(mass_kg)         0.1520478 0.03964739  3.835001 1.555048e-04
# flyer1                -1.1057437 0.23552817 -4.694741 4.203755e-06
# log10(mass_kg):flyer1 -0.1629897 0.07374201 -2.210269 2.790399e-02

# all significant

gls_model_GIT_1_mass_simple <- gls(log_transit_hrs ~ flyer + log10(mass_kg), 
                             data = GIT_pruned, 
                             correlation = cor_phylo_fixed1,
                             method = "ML")
summary(gls_model_GIT_1_mass_simple)$tTable

# Value  Std.Error   t-value      p-value
# (Intercept)     1.7165445 0.25673487  6.686059 1.253521e-10
# flyer1         -0.8703255 0.21153166 -4.114398 5.119077e-05
# log10(mass_kg)  0.1048389 0.03363492  3.116967 2.018691e-03

AIC(gls_model_GIT_1_mass, gls_model_GIT_1_mass_simple)

# df       AIC
# gls_model_GIT_1_mass         5 -1882.709
# gls_model_GIT_1_mass_simple  4 -1879.796

# there is no difference in AIC,really, though the complex model is slightly better
#so I'm going to write up the simple model because they are so close and significance is across the board
# (no interaction term), shown here:
summary(gls_model_GIT_1_mass_simple)$tTable

# Value  Std.Error   t-value      p-value
# (Intercept)     1.7165445 0.25673487  6.686059 1.253521e-10
# flyer1         -0.8703255 0.21153166 -4.114398 5.119077e-05
# log10(mass_kg)  0.1048389 0.03363492  3.116967 2.018691e-03

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
# (Intercept)                    1.60872886 0.25306651  6.3569410 8.606947e-10
# flyer1                        -0.83242737 0.20844452 -3.9935200 8.375764e-05
# log10(mass_kg)                 0.12728405 0.03394380  3.7498473 2.160343e-04
# trial.dietfruit/nectar/pollen -0.01028957 0.06777772 -0.1518134 8.794462e-01
# trial.dietmeat                -0.05486068 0.09724550 -0.5641462 5.731178e-01
# trial.dietmixed               -0.19264083 0.24113996 -0.7988756 4.250570e-01
# trial.dietprotein              0.10196209 0.05145512  1.9815731 4.853003e-02
# trial.dietunknown              0.43686686 0.11778594  3.7089898 2.520584e-04

# then, as a random effect

gls_model_GIT_1_mass_simple_diet_re <- lme(log_transit_hrs ~ flyer + log10(mass_kg),
                                            random = ~1 | trial.diet,
                                            data = GIT_pruned, 
                                            correlation = cor_phylo_fixed1,
                                            method = "ML")
summary(gls_model_GIT_1_mass_simple_diet_re)$tTable
# Value  Std.Error  DF   t-value      p-value
# (Intercept)     1.4122715 0.20047312 273  7.044693 1.513850e-11
# flyer1         -0.7193753 0.19630612 273 -3.664558 2.976293e-04
# log10(mass_kg)  0.1527351 0.03090193 273  4.942575 1.347166e-06

# then, compare
AIC(gls_model_GIT_1_mass_simple, gls_model_GIT_1_mass_simple_diet, gls_model_GIT_1_mass_simple_diet_re)

# df        AIC
# gls_model_GIT_1_mass_simple          4 -1879.7960
# gls_model_GIT_1_mass_simple_diet     9 -1893.3873 *****
# gls_model_GIT_1_mass_simple_diet_re  5  -639.2439

# here, mass and diet does improve model fit when diet is modeled as a fixed effect! 
summary(gls_model_GIT_1_mass_simple_diet)$tTable

# Value  Std.Error    t-value      p-value
# (Intercept)                    1.60872886 0.25306651  6.3569410 8.606947e-10
# flyer1                        -0.83242737 0.20844452 -3.9935200 8.375764e-05 
# log10(mass_kg)                 0.12728405 0.03394380  3.7498473 2.160343e-04 *
# trial.dietfruit/nectar/pollen -0.01028957 0.06777772 -0.1518134 8.794462e-01
# trial.dietmeat                -0.05486068 0.09724550 -0.5641462 5.731178e-01
# trial.dietmixed               -0.19264083 0.24113996 -0.7988756 4.250570e-01
# trial.dietprotein              0.10196209 0.05145512  1.9815731 4.853003e-02 *
# trial.dietunknown              0.43686686 0.11778594  3.7089898 2.520584e-04 *

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
# gls_model_GIT_1_diet              8 -1881.274
# gls_model_GIT_1_mass_simple_diet  9 -1893.387



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
# (Intercept)    0.8699079 0.04613644 18.85511 1.057902e-51
# log10(mass_kg) 0.3334617 0.02755288 12.10261 2.195480e-27

summary(gls_model_GIT_1_mass_only)$tTable
# mass is hugely significant and positively related to gut transit

# Value  Std.Error  t-value      p-value
# (Intercept)    1.5453020 0.26047005 5.932743 8.802754e-09
# log10(mass_kg) 0.1553713 0.03219378 4.826130 2.295368e-06

AIC(gls_model_GIT_0_mass_only, gls_model_GIT_1_mass_only)
# incorporating phylogeny with mass offers a much better fit to the data however!
# 
# df        AIC
# gls_model_GIT_0_mass_only  3   606.9863
# gls_model_GIT_1_mass_only  3 -1865.1858




###### ------ PLOTTING FOR FIGURE 2
library(nlme)
library(ggplot2)
library(dplyr)

# --- Colors and shapes ---
# [1] "Flying Birds"        "Carnivores"          "Even-toed Ungulates" "Primates"            "Bats"               
# [6] "Reptiles"            "Rodents"             "Shrews/Moles"        "Odd-toed Ungulates"  "Marsupials"   

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
# [6] "Reptiles"            "Rodents"             "Shrews/Moles"        "Odd-toed Ungulates"  "Marsupials"       

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
  coord_cartesian(ylim = c(-0.8, 3.5)) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid = element_blank(),
    legend.position = c(0.15, 0.85),
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
# (Intercept)  1.6456490 0.25749376  6.391025   0e+00
# flyer1      -0.8703255 0.21153166 -4.114398   1e-04
# log_mass_c   0.1048389 0.03363492  3.116967   2e-03



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
