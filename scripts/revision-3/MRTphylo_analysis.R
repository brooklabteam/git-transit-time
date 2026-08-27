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

#load the MRT transit data with phylo name
dat_MRT <- read.csv(file = paste0(homewd, "/data/revision-3/dat_sum_tot_clean_MRT_R3.csv"), header = T, stringsAsFactors = F )

mass_col    <- "avg_mass"        # body-mass column
species_col <- "genus.species"  # tip / phylogeny-name column

# species missing a mass in either dataset
miss <- function(d) d[[species_col]][is.na(d[[mass_col]]) | trimws(d[[mass_col]]) == "0"]
no_mass <- sort(unique(c(miss(dat_MRT))))

# report which ones
cat("Species without mass (", length(no_mass), "):\n", sep = "")
print(no_mass) #13
# [1] "Alopochen aegyptiaca"     "Anser anser"              "Dromaius novaehollandiae" "Gulosus aristotelis"     
# [5] "Macronectes giganteus"    "Martes melampus"          "Meleagris gallopavo"      "Mitu salvini"            
# [9] "Plectopterus gambensis"   "Pygoscelis papua"         "Rhea americana"           "Stercorarius skua"       
# [13] "Struthio camelus"     

# remove them from both datasets
dat_MRT <- dat_MRT[!dat_MRT[[species_col]] %in% no_mass, ]

#### species counts ####
#table 1 in the manuscript

length(unique(dat_MRT$genus.species)) #240 unique species for transit time
tt.pivot <-ddply(dat_MRT, .(re_class), summarise, N_species_transit = length(unique(genus.species)))
print(tt.pivot)

# re_class N_species_transit
# 1                Bats                13
# 2          Carnivores                16
# 3        Flying Birds                61
# 4          Marsupials                14
# 5    Non-Flying Birds                 2
# 6  Odd-toed Ungulates                10
# 7            Primates                40
# 8            Reptiles                 3
# 9             Rodents                27
# 10          Ungulates                54


dat_MRT <- subset(dat_MRT, dat_MRT$re_class != "Non-Flying Birds")
dat_MRT <- subset(dat_MRT, dat_MRT$re_class != "Reptiles")
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
MRT_hrs <- MRT_hrs[common_species] #222

#for MRT
lambda_gs2<-phylosig(tree, MRT_hrs,
                     method="lambda",test=TRUE)
lambda_gs2

# Phylogenetic signal lambda : 0.941198 
# logL(lambda) : -966.376 
# LR(lambda=0) : 134.609 
# P-value (based on LR test) : 4.01964e-31 



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

#[1] "Gazella"           "Nectarinia_famosa"

#sort(tree$tip.label)

###################### MGT - model 1 ################


# Fit PGLS model
#ML=maximum likelihood

# factor the groups so that rodent is the comparasion
unique(MRT_pruned$re_class)
# [1] "Flying Birds"       "Ungulates"          "Primates"           "Marsupials"         "Carnivores"         "Rodents"           
# [7] "Odd-toed Ungulates" "Bats" 

length(unique(MRT_pruned$phylo_name)) #222


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
# Phylogenetic signal lambda : 0.942855


# Cara: now test the effect of flyer on MRT transit without accounting for mass or diet
# All text below is from Cara

#with full phylogenetic effects:
gls_model_MRT_1 <- gls(log_MRT_hrs ~ flyer, 
                         data = MRT_pruned, 
                        correlation = cor_phylo_fixed1,
                         method = "ML")
summary(gls_model_MRT_1)$tTable
# Value Std.Error   t-value     p-value
# (Intercept)  0.9189334 0.3804043  2.415676 0.016337307
# flyer1      -0.9032493 0.3026176 -2.984788 0.003084689

#with no phylogenetic effects:
gls_model_MRT_0 <- gls(log_MRT_hrs ~ flyer, 
                        data = MRT_pruned, 
                        #correlation = cor_phylo_fixed0,
                        method = "ML")
summary(gls_model_MRT_0)$tTable
# Value  Std.Error   t-value      p-value
# (Intercept)  1.0584276 0.04461767 23.722163 2.474435e-69
# flyer1      -0.5934719 0.20927542 -2.835841 4.899079e-03

#with estimated phylogenetic effects:
gls_model_MRT_2 <- gls(log_MRT_hrs ~ flyer,
                        data = MRT_pruned,
                        correlation = cor_phylo_fixed2,
                        method = "ML")
summary(gls_model_MRT_2)$tTable
# Value Std.Error   t-value     p-value
# (Intercept)  1.0454412 0.3351832  3.119015 0.002001056
# flyer1      -0.9060088 0.2736555 -3.310765 0.001050674

AIC(gls_model_MRT_0, gls_model_MRT_1, gls_model_MRT_2)
# Emily: The model with full brownian phylogenetic effects is the best fit.
# df        AIC
# gls_model_MRT_0  3   641.2291
# gls_model_MRT_1  3 -2536.6349 **
# gls_model_MRT_2  3 -2337.8360



summary(gls_model_MRT_1)$tTable
# Value Std.Error   t-value     p-value
# (Intercept)  0.9189334 0.3804043  2.415676 0.016337307
# flyer1      -0.9032493 0.3026176 -2.984788 0.003084689

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
# (Intercept)                    0.87122645 0.38061330  2.2890068 0.022826232
# flyer1                        -0.93943401 0.30933013 -3.0369949 0.002615539
# trial.dietfruit/nectar/pollen -0.01033575 0.09025405 -0.1145184 0.908909254
# trial.dietmeat                 0.31635315 0.14526655  2.1777426 0.030262873
# trial.dietmixed                0.16245137 0.49938818  0.3253008 0.745197243
# trial.dietprotein              0.07246705 0.12758207  0.5680034 0.570489499
# trial.dietunknown              0.16803373 0.25660246  0.6548407 0.513109834



# With no phylogenetic effects:
gls_model_MRT_0_diet <- gls(log_MRT_hrs ~ flyer + trial.diet, 
                        data = MRT_pruned, 
                        #correlation = cor_phylo_fixed0,
                        method = "ML")
summary(gls_model_MRT_0_diet)$tTable
# Value  Std.Error    t-value      p-value
# (Intercept)                    1.1723418 0.05555751 21.1014105 9.752837e-60
# flyer1                        -0.3914769 0.23656279 -1.6548542 9.907851e-02
# trial.dietfruit/nectar/pollen -0.4165126 0.11245012 -3.7039766 2.558304e-04
# trial.dietmeat                -0.1773151 0.14539883 -1.2195088 2.236811e-01
# trial.dietmixed                0.1018161 0.72438128  0.1405559 8.883222e-01
# trial.dietprotein             -0.3096478 0.16986671 -1.8228867 6.939069e-02
# trial.dietunknown              0.3533181 0.36537246  0.9670081 3.343775e-01


# With estimated phylogenetic effects:
gls_model_MRT_2_diet <- gls(log_MRT_hrs ~ flyer + trial.diet,
                        data = MRT_pruned,
                        correlation = cor_phylo_fixed2,
                        method = "ML")
summary(gls_model_MRT_2_diet)$tTable
# Value Std.Error    t-value      p-value
# (Intercept)                    0.99017419 0.3341985  2.9628323 0.0033109489
# flyer1                        -0.95173034 0.2792512 -3.4081508 0.0007504681
# trial.dietfruit/nectar/pollen  0.01131968 0.0827610  0.1367755 0.8913068911
# trial.dietmeat                 0.34994613 0.1330169  2.6308400 0.0089908380
# trial.dietmixed                0.19007558 0.4504509  0.4219674 0.6733737894
# trial.dietprotein              0.08325051 0.1172196  0.7102101 0.4781670742
# trial.dietunknown              0.17332156 0.2346984  0.7384864 0.4608397034

AIC(gls_model_MRT_0, gls_model_MRT_1, gls_model_MRT_2, 
    gls_model_MRT_0_diet, gls_model_MRT_1_diet, gls_model_MRT_2_diet)
# df        AIC
# gls_model_MRT_0       3   641.2291
# gls_model_MRT_1       3 -2536.6349 *
# gls_model_MRT_2       3 -2337.8360
# gls_model_MRT_0_diet  8   634.4242
# gls_model_MRT_1_diet  8 -2532.5178 *
# gls_model_MRT_2_diet  8 -2335.7647 *

# but remember that lambda is 0.94 for model2 

summary(gls_model_MRT_1_diet)$tTable
# Value  Std.Error    t-value     p-value
# (Intercept)                    0.87122645 0.38061330  2.2890068 0.022826232
# flyer1                        -0.93943401 0.30933013 -3.0369949 0.002615539 *
# trial.dietfruit/nectar/pollen -0.01033575 0.09025405 -0.1145184 0.908909254
# trial.dietmeat                 0.31635315 0.14526655  2.1777426 0.030262873 *
# trial.dietmixed                0.16245137 0.49938818  0.3253008 0.745197243
# trial.dietprotein              0.07246705 0.12758207  0.5680034 0.570489499
# trial.dietunknown              0.16803373 0.25660246  0.6548407 0.513109834

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
# (Intercept)  0.8321836 0.2133704 279  3.900183 1.205009e-04
# flyer1      -1.3387866 0.2577807 279 -5.193509 3.985442e-07

gls_model_MRT_0_diet_re <- lme(log_MRT_hrs ~ flyer,
                                random = ~1 | trial.diet,
                                #correlation = cor_phylo_fixed0,
                                data = MRT_pruned,
                                method = "ML")
summary(gls_model_MRT_0_diet_re)$tTable
# Value  Std.Error  DF   t-value      p-value
# (Intercept)  1.0026173 0.08888188 279 11.280334 1.456172e-24
# flyer1      -0.4648389 0.22255133 279 -2.088682 3.764272e-02

gls_model_MRT_2_diet_re <- lme(log_MRT_hrs ~ flyer,
                                random = ~1 | trial.diet,
                                correlation = cor_phylo_fixed2,
                                data = MRT_pruned,
                                method = "ML")
summary(gls_model_MRT_2_diet_re)$tTable
# Value Std.Error  DF   t-value      p-value
# (Intercept)  0.9135057 0.1702606 279  5.365339 1.699433e-07
# flyer1      -1.0783676 0.2178704 279 -4.949582 1.288623e-06

AIC(gls_model_MRT_0, gls_model_MRT_0_diet, gls_model_MRT_0_diet_re,
    gls_model_MRT_1, gls_model_MRT_1_diet, gls_model_MRT_1_diet_re,
    gls_model_MRT_2, gls_model_MRT_2_diet, gls_model_MRT_2_diet_re)

# df        AIC
# gls_model_MRT_0          3   641.2291
# gls_model_MRT_0_diet     8   634.4242
# gls_model_MRT_0_diet_re  4   636.9633
# gls_model_MRT_1          3 -2536.6349 * --- best supported thus far
# gls_model_MRT_1_diet     8 -2532.5178
# gls_model_MRT_1_diet_re  4 -1331.1706
# gls_model_MRT_2          3 -2337.8360 * 
# gls_model_MRT_2_diet     8 -2335.7647 *
# gls_model_MRT_2_diet_re  4 -1242.2332

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
# (Intercept)            0.9636360 0.36506005  2.639664 8.760712e-03
# log10(mass_kg)         0.2110837 0.04154097  5.081338 6.831644e-07
# flyer1                -1.3562368 0.63800296 -2.125753 3.439274e-02
# log10(mass_kg):flyer1 -0.4678003 0.25691434 -1.820841 6.969066e-02


#with no phylogenetic effects:
gls_model_MRT_0_mass <- gls(log_MRT_hrs ~ log10(mass_kg)*flyer, 
                        data = MRT_pruned, 
                        #correlation = cor_phylo_fixed0,
                        method = "ML")
summary(gls_model_MRT_0_mass)$tTable
# Value  Std.Error   t-value      p-value
# (Intercept)            1.0411798 0.03700275 28.137903 7.418260e-84
# log10(mass_kg)         0.2845213 0.02493923 11.408582 4.841155e-25
# flyer1                -1.4546686 0.74911485 -1.941850 5.315100e-02
# log10(mass_kg):flyer1 -0.6321109 0.28943972 -2.183912 2.979258e-02

#with estimated phylogenetic effects:
gls_model_MRT_2_mass <- gls(log_MRT_hrs ~ log10(mass_kg)*flyer,
                        data = MRT_pruned,
                        correlation = cor_phylo_fixed2,
                        method = "ML")
summary(gls_model_MRT_2_mass)$tTable
# Value  Std.Error   t-value      p-value
# (Intercept)            1.0831822 0.32689521  3.313546 1.041547e-03
# log10(mass_kg)         0.1614776 0.04051237  3.985884 8.568092e-05
# flyer1                -1.3948146 0.59963246 -2.326116 2.072131e-02
# log10(mass_kg):flyer1 -0.4189896 0.24124694 -1.736767 8.352036e-02



AIC(gls_model_MRT_0_mass, gls_model_MRT_1_mass, gls_model_MRT_2_mass)
#Emily: model 1 is the best
#Cara: This shows that transit time is still influenced by phylogeny, not
# exclusively by mass (e.g. the model with phylogeny AND mass performed better
# than the model with mass alone)

# df        AIC
# gls_model_MRT_0_mass  5   535.6873
# gls_model_MRT_1_mass  5 -2558.3849 *
# gls_model_MRT_2_mass  5 -2350.4346

summary(gls_model_MRT_1_mass)$tTable

# Value  Std.Error   t-value      p-value
# (Intercept)            0.9636360 0.36506005  2.639664 8.760712e-03
# log10(mass_kg)         0.2110837 0.04154097  5.081338 6.831644e-07
# flyer1                -1.3562368 0.63800296 -2.125753 3.439274e-02
# log10(mass_kg):flyer1 -0.4678003 0.25691434 -1.820841 6.969066e-02

# Emily: flyer is significant (negative) term here (p=0.03)
# We can test whether it would be best to drop the interaction which 
# might clarify the variables that matter most.


gls_model_MRT_1_mass_simple <- gls(log_MRT_hrs ~ flyer + log10(mass_kg), 
                             data = MRT_pruned, 
                             correlation = cor_phylo_fixed1,
                             method = "ML")
summary(gls_model_MRT_1_mass_simple)$tTable
# Value  Std.Error   t-value      p-value
# (Intercept)     0.9627713 0.36655013  2.626575 9.094682e-03
# flyer1         -0.3435944 0.31393642 -1.094471 2.746789e-01 - not significant
# log10(mass_kg)  0.1966148 0.04094033  4.802471 2.543771e-06


AIC(gls_model_MRT_1_mass, gls_model_MRT_1_mass_simple)

# df       AIC
# gls_model_MRT_1_mass         5 -2558.385
# gls_model_MRT_1_mass_simple  4 -2557.042


# Emily: there is no difference in AIC,really, though the complex model is slightly better
#so we prefer the simplest model 
# (no interaction term), shown here:
summary(gls_model_MRT_1_mass_simple)$tTable
# Value  Std.Error   t-value      p-value
# (Intercept)     0.9627713 0.36655013  2.626575 9.094682e-03
# flyer1         -0.3435944 0.31393642 -1.094471 2.746789e-01
# log10(mass_kg)  0.1966148 0.04094033  4.802471 2.543771e-06


# flyer is not significantly negatively related to transit suggesting 
# that flyers have rapid MRT independent of small body size; 

# Emily: note that flight is significant in the complex model. so this is really a toss up.

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
# (Intercept)                                   0.92136859 0.36696315  2.5107932 1.261585e-02
# flyer1                                       -0.41040091 0.31785767 -1.2911468 1.977292e-01
# log10(mass_kg)                                0.20376741 0.04241481  4.8041575 2.549326e-06
# trial.dietfruit/nectar/pollen                 0.02335313 0.08750141  0.2668886 7.897535e-01
# trial.dietfruit/nectar/pollen; fiber/foliage  0.12500050 0.47665148  0.2622472 7.933258e-01
# trial.dietfruit/nectar/pollen; meat           0.15418008 0.48130226  0.3203394 7.489526e-01
# trial.dietmeat                                0.26612168 0.14038877  1.8956052 5.905311e-02
# trial.dietprotein                             0.18371095 0.12509107  1.4686177 1.430712e-01
# trial.dietunknown                             0.32024746 0.24938868  1.2841299 2.001695e-01

# then, as a random effect
gls_model_MRT_1_mass_simple_diet_re <- lme(log_MRT_hrs ~ flyer + log10(mass_kg),
                                            random = ~1 | trial.diet,
                                            data = MRT_pruned, 
                                            correlation = cor_phylo_fixed1,
                                            method = "ML")
summary(gls_model_MRT_1_mass_simple_diet_re)$tTable
# Value  Std.Error  DF   t-value      p-value
# (Intercept)     0.93318698 0.21187154 277  4.404494 1.516434e-05
# flyer1         -1.07190391 0.27496359 277 -3.898349 1.215604e-04
# log10(mass_kg)  0.09943869 0.03718593 277  2.674094 7.938326e-03

# then, compare
AIC(gls_model_MRT_1_mass_simple, gls_model_MRT_1_mass_simple_diet, gls_model_MRT_1_mass_simple_diet_re)

# df       AIC
# gls_model_MRT_1_mass_simple          4 -2557.042
# gls_model_MRT_1_mass_simple_diet    10 -2551.696
# gls_model_MRT_1_mass_simple_diet_re  5 -1294.370

# here, Adding diet does nothing really for the fit. So only have mass in the model is still the best

summary(gls_model_MRT_1_mass_simple)$tTable

# Value  Std.Error   t-value      p-value
# (Intercept)     0.9627713 0.36655013  2.626575 9.094682e-03
# flyer1         -0.3435944 0.31393642 -1.094471 2.746789e-01
# log10(mass_kg)  0.1966148 0.04094033  4.802471 2.543771e-06

# Emily: flyer is not significantly negatively related to transit, 
# after accounting for mass (and phylogeny). 


# Cara: what if we compare our two best fit models so far? 
AIC(gls_model_MRT_1, gls_model_MRT_1_diet, gls_model_MRT_1_mass_simple)

# df       AIC
# gls_model_MRT_1              3 -2536.635
# gls_model_MRT_1_diet         8 -2532.518
# gls_model_MRT_1_mass_simple  4 -2557.042 **


gls_model_MRT_1_mass_only <- gls(log_MRT_hrs ~ log10(mass_kg), 
                                   data = MRT_pruned, 
                                   correlation = cor_phylo_fixed1,
                                   method = "ML")
summary(gls_model_MRT_1_mass_only)$tTable
# Value  Std.Error  t-value      p-value
# (Intercept)    0.9564153 0.36663178 2.608654 9.571165e-03
# log10(mass_kg) 0.2132477 0.03802842 5.607588 4.862342e-08

gls_model_MRT_0_mass_only <- gls(log_MRT_hrs ~ log10(mass_kg), 
                                 data = MRT_pruned, 
                                 #correlation = cor_phylo_fixed0,
                                 method = "ML")
summary(gls_model_MRT_0_mass_only)$tTable
# Value  Std.Error  t-value      p-value
# (Intercept)    1.0470566 0.03635364 28.80197 2.995939e-86
# log10(mass_kg) 0.2737247 0.02343028 11.68252 5.191835e-26

gls_model_MRT_2_mass_only <- gls(log_MRT_hrs ~ log10(mass_kg), 
                                 data = MRT_pruned, 
                                 correlation = cor_phylo_fixed2,
                                 method = "ML")
summary(gls_model_MRT_2_mass_only)$tTable
# Value  Std.Error  t-value      p-value
# (Intercept)    1.0747564 0.32904382 3.266302 1.223277e-03
# log10(mass_kg) 0.1737762 0.03680466 4.721583 3.681160e-06

# mass is always signficant

AIC(gls_model_MRT_0_mass_only, gls_model_MRT_1_mass_only, gls_model_MRT_2_mass_only)
# 
# df        AIC
# gls_model_MRT_0_mass_only  3   536.9815
# gls_model_MRT_1_mass_only  3 -2557.8340 **
# gls_model_MRT_2_mass_only  3 -2348.6181








###### ------ PLOTTING FOR FIGURE S2
library(nlme)
library(ggplot2)
library(dplyr)

unique(MRT_pruned$re_class)
# [1] "Flying Birds"       "Ungulates"          "Primates"           "Marsupials"         "Carnivores"         "Rodents"           
# [7] "Odd-toed Ungulates" "Bats"        

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
MRT_pruned$re_class[MRT_pruned$re_class == "Ungulates"] <- "Even-toed\nUngulates"
MRT_pruned$re_class[MRT_pruned$re_class == "Odd-toed Ungulates"] <- "Odd-toed\nUngulates"

unique(MRT_pruned$re_class)
# [1] "Flying Birds"         "Even-toed\nUngulates" "Primates"             "Marsupials"           "Carnivores"          
# [6] "Rodents"              "Odd-toed\nUngulates"  "Bats"   

MRT_pruned$re_class <- factor(MRT_pruned$re_class, levels = c("Bats", "Flying\nBirds", "Rodents", "Primates", 
                                                              "Marsupials", "Carnivores",
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
