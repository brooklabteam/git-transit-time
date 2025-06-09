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
homewd <- "/Users/emilyruhs/Desktop/UChi_Brook_Lab/GitHub_repos/git-transit-time/"
homewd <- "/Users/katherinemcferrin/Developer/git-transit-time/"

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
library(caper)

#caper automatically builds the correlation structure using the comparative data object
# Load your tree
#tree <- read.tree("Genus_species_list_2.nwk")

# Load your data
dat_GIT <- read.csv(file = paste0(homewd, "/data/dat_sum_tot_clean_GIT.csv"), header = T, stringsAsFactors = F )
names(dat_GIT)

# clean dataset with only the columns you need
GIT <- dat_GIT[c("phylo_name", "transit_hrs", "typical.diet", "re_class", "mass_kg")]
#GIT <- na.omit(GIT$transit_hrs)

tree <- read.tree("timetree_names.nwk")


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

GIT_pruned$re_class <- factor(GIT_pruned$re_class, levels = c("Rodents","Flying Birds","Bats", "Carnivores", "Primates", "Ungulates", "Reptiles"))

#transform data
GIT_pruned$log_transit_hrs <- log10(GIT_pruned$transit_hrs)

# Create correlation structure
cor_phylo_fixed1 <- corPagel(1, phy = tree_pruned, fixed = TRUE, form = ~phylo_name)


#converges fine with 'transit_hrs' but not log_transit_hrs
pgls_model_GIT_1 <- gls(log_transit_hrs ~ re_class, 
                         data = GIT_pruned, 
                        correlation = cor_phylo_fixed1,
                         method = "ML",
                         weights = varIdent(form = ~1|typical.diet))


#check that residuals are normal
hist(resid(pgls_model_GIT_1), main = "Residual Histogram", xlab = "Residuals")

# View the results
summary(pgls_model_GIT_1)

# Value Std.Error   t-value p-value
# (Intercept)           1.0820500 0.3358128  3.222182  0.0017
# re_classFlying Birds -0.8439406 0.4802604 -1.757256  0.0818
# re_classBats         -1.0087676 0.3814457 -2.644590  0.0094
# re_classCarnivores    0.0976089 0.4781655  0.204132  0.8386
# re_classPrimates      0.2547125 0.3530821  0.721398  0.4723
# re_classUngulates     0.4491604 0.4895050  0.917581  0.3610
# re_classReptiles      1.0378403 0.5994547  1.731307  0.0864


pgls_model_GIT_1_reduced <- gls(log_transit_hrs ~ re_class, 
                                data = GIT_pruned, 
                                correlation = cor_phylo_fixed1,
                                method = "ML")


AIC(pgls_model_GIT_1, pgls_model_GIT_1_reduced)
# df      AIC
# pgls_model_GIT_1         12 160.3381
# pgls_model_GIT_1_reduced  8 179.1278

# delta AIC is 18.7897



############### Plot MGT - model 1 ###############

# Add residuals and fitted values to your data
GIT_pruned$resid <- resid(pgls_model_GIT_1, type = "normalized")
GIT_pruned$fitted <- fitted(pgls_model_GIT_1)


D <- ggplot(GIT_pruned, aes(x = fitted, y = resid, shape = typical.diet)) +
  geom_point(size=3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_shape_manual(values=c(21,23,22,24,25)) +
  theme_bw() +
  labs(title = "Residuals vs Fitted by Diet Type",
       y = "Normalized Residuals\n\n", x = "\nFitted Values",
       shape = "Natural Diet") +
  theme(axis.text.x=element_text(size=11),
        axis.text.y=element_text(size=11))+
  theme(
    legend.text = element_text(size=8),
    legend.title = element_text(size=10),
    legend.position = c(0.8,0.8),
    legend.background = element_rect(fill="white", color="black"))

ggsave(file = paste0(homewd,"/figures/MGT_residualVSfitted.jpeg"),
       plot=D,
       units="mm",  
       width=60, 
       height=70, 
       scale=3, 
       dpi=300)





summary(pgls_model_GIT_1)$modelStruct$varStruct

# Extract variance weights
var_weights <- coef(pgls_model_GIT_1$modelStruct$varStruct, unconstrained = FALSE)
var_weights_df <- data.frame(
  typical_diet = names(var_weights),
  relative_variance = var_weights
)


C <- ggplot(GIT_pruned, aes(x = typical.diet, y = resid)) +
  geom_boxplot() +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x=element_text(size=11),
        axis.text.y=element_text(size=11))+
  labs(title = "Residual Spread by Natural Diet",
       y = "\nNormalized Residuals",
       x="")

ggsave(file = paste0(homewd,"/figures/MGT_resid.jpeg"),
       plot=C,
       units="mm",  
       width=70, 
       height=70, 
       scale=3, 
       dpi=300)









#and reshape order
GIT_pruned$re_class <- as.character(GIT_pruned$re_class)
GIT_pruned$re_class[GIT_pruned$re_class=="Flying Birds"] <- "Flying\nBirds"
#GIT_pruned$re_class[GIT_pruned$re_class=="Non-Flying Birds"] <- "Non-Flying\nBirds"
#dat.sum.tot_clean$re_class[dat.sum.tot_clean$re_class=="Non-Flying Birds"] <- "Non-Flying\nBirds"
GIT_pruned$re_class <- factor(GIT_pruned$re_class, levels = c( "Flying\nBirds","Bats","Rodents", "Carnivores","Primates", "Ungulates", "Reptiles"))
#GIT_pruned$typical_diet = factor(GIT_pruned$typical_diet, levels= c("ommnivore", "insectivore", "frugivore/nectarivore", "folivore", "carnivore")) #this is messing up and adding NAs



#scales::show_col(scales::hue_pal()(9))
colz =  c( "Flying\nBirds"= "#F8766D", "Non-Flying\nBirds"="#edf8b1","Bats" = "#C49A00", "Rodents"= "navy","Carnivores" = "#00C094", "Primates"= "#00B6EB", "Ungulates"= "#A58AFF",   "Reptiles" = "#FB61D7")

#shapez = c("fiber/foliage" = 21, "fruit/nectar/pollen" = 22, "label-solution" = 23, "protein" = 24)
shapez = c("carnivore" = 21, "frugivore/nectarivore" = 22, "folivore" = 23, "insectivore" = 24, "omnivore" = 25)


#and the version that is 2panel
label.dat <- cbind.data.frame(re_class=c("Flying\nBirds","Bats", "Rodents","Carnivores", "Primates", "Ungulates", "Reptiles"), label = c("", "**", "", "","","",""))

# Value Std.Error   t-value p-value
# (Intercept)           1.0820500 0.3358128  3.222182  0.0017
# re_classFlying Birds -0.8439406 0.4802604 -1.757256  0.0818
# re_classBats         -1.0087676 0.3814457 -2.644590  0.0094
# re_classCarnivores    0.0976089 0.4781655  0.204132  0.8386
# re_classPrimates      0.2547125 0.3530821  0.721398  0.4723
# re_classUngulates     0.4491604 0.4895050  0.917581  0.3610
# re_classReptiles      1.0378403 0.5994547  1.731307  0.0864

label.dat$re_class <- factor(label.dat$re_class, levels=c("Flying\nBirds","Bats","Rodents","Carnivores", "Primates", "Ungulates", "Reptiles"))


#and all together
pA <- ggplot(data=GIT_pruned)+ geom_boxplot(aes(x=re_class, y=log10(transit_hrs), fill=re_class), show.legend = F) + theme_bw() + 
  theme(axis.title.x = element_blank(), panel.grid = element_blank(), legend.position = c(.85,.15), legend.background = element_rect(color="black"),
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=12), axis.title.y = element_text(size=15),
        plot.margin = unit(c(.2,.2,.8,.2), "cm")) + scale_fill_manual(values=colz) + 
  geom_point(aes(x=re_class, y=log10(transit_hrs),shape=typical.diet), position = position_jitterdodge(jitter.width = 0), size=3) +
  scale_shape_manual(values=shapez, name="Natural diet") +  ylab("GIT transit time, hrs") +
  geom_hline(yintercept=0.6684442, linetype=2) +
  scale_y_continuous(breaks=c(0,1,2), labels = c("1", "10", "100")) +
  coord_cartesian(ylim=c(-.8,2.95)) + geom_label(data=label.dat, aes(x=re_class, y=3.0, label=label), label.size = NA,size=6) 
  #scale_x_discrete(labels = c("Flying\nBirds","Bats","Non-Flying\nBirds", "Rodents", "Carnivores", "Marsupials","Lagomorphs","Primates", "Ungulates", "Reptiles"))
print(pA)


ggsave(file = paste0(homewd,"/figures/Fig_1A.png"),
       units="mm",  
       width=120, 
       height=80, 
       scale=3, 
       dpi=300)




###################### MGT - model 2 - mass ################
GIT_pruned <- na.omit(GIT_pruned)
pgls_model_GIT_2 <- gls(log_transit_hrs ~ log10(mass_kg):re_class, 
                        data = GIT_pruned, 
                        correlation = cor_phylo_fixed1,
                        method = "ML",
                        weights = varIdent(form = ~1|typical.diet))

# check the residuals
hist(resid(pgls_model_GIT_2), main = "Residual Histogram", xlab = "Residuals")

# View the results
summary(pgls_model_GIT_2)

# Value Std.Error   t-value p-value
# (Intercept)                           1.3778763 0.12781705 10.780067  0.0000
# log10(mass_kg):re_classFlying\nBirds  0.4635773 0.15470984  2.996431  0.0034
# log10(mass_kg):re_classBats           0.3298203 0.07959178  4.143898  0.0001
# log10(mass_kg):re_classRodents        0.3608113 0.06710976  5.376435  0.0000
# log10(mass_kg):re_classCarnivores     0.1017976 0.28554113  0.356508  0.7222
# log10(mass_kg):re_classPrimates       0.0184652 0.09294120  0.198676  0.8429
# log10(mass_kg):re_classUngulates     -0.2159061 0.10978986 -1.966539  0.0520
# log10(mass_kg):re_classReptiles       0.1240345 0.14717984  0.842741  0.4014



pgls_model_GIT_2_reduced <- gls(log_transit_hrs ~ log10(mass_kg):re_class, 
                        data = GIT_pruned, 
                        correlation = cor_phylo_fixed1,
                        method = "ML")

AIC(pgls_model_GIT_2,pgls_model_GIT_2_reduced)
# df      AIC
# pgls_model_GIT_2         13 146.3067
# pgls_model_GIT_2_reduced  9 165.8505

# delta AIC = 19.3067




colz =  c( "Flying\nBirds"= "#F8766D", "Bats" = "#C49A00", "Non-Flying Birds"="#edf8b1", "Rodents"= "navy","Carnivores" = "#00C094","Primates"= "#00B6EB", "Ungulates"= "#A58AFF",   "Reptiles" = "#FB61D7")

#shapez = c("fiber/foliage" = 21, "fruit/nectar/pollen" = 22, "label-solution" = 23, "protein" = 24)
shapez = c("carnivore" = 21, "frugivore/nectarivore" = 22, "folivore" = 23, "insectivore" = 24, "omnivore" = 25)


pB <- ggplot(data = GIT_pruned, aes(x = log10(mass_kg), y = log10(transit_hrs))) + 
  geom_mark_ellipse(expand = 0, radius = 0, aes(fill = re_class, color = re_class), size = 0.1) +
  #geom_line(aes(x = log10(mass_kg), y = log10(predicted_transit), color = re_class, linetype = re_class), 
            #alpha = 0.6, size = 0.7, show.legend = TRUE) +
  scale_linetype_manual(values = c("solid", "dotted","solid", "solid", "solid", "solid", "dotted", "dotted", "solid", "dotted"), guide="none") +
  geom_point(aes(x = log10(mass_kg), y = log10(transit_hrs), fill = re_class, shape = typical.diet), 
             size = 3, show.legend = FALSE) + 
  ylab("GIT transit time, hrs") + xlab("Mass (kg)") + 
  scale_shape_manual(name = "Food Type", values = shapez) + 
  scale_fill_manual(name = "Taxonomic Group", values = colz) +  # Adjusted legend name for re_class
  scale_color_manual(name = "Taxonomic Group", values = colz) +  # Adjusted legend name for re_class
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        legend.position = c(0.85, 0.2),  # Adjusted legend position if needed
        legend.background = element_rect(color = "black"),
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 15),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")) +
  scale_y_continuous(breaks = c(0, 1, 2), labels = c("1", "10", "100")) +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2, 3), labels = c(".01", ".1", "1", "10", "100", "1000")) +
  coord_cartesian(ylim = c(-0.8, 2.9)) +
  guides(color = guide_legend(title = "Taxonomic Group"),  # Adjusted legend title for re_class
         fill = guide_legend(title = "Taxonomic Group"),  # Adjusted legend title for re_class
         shape = "none")  # Remove legend for shape (food.cat)

print(pB)

ggsave(file = paste0(homewd,"/figures/Fig_1B.png"),
       units="mm",  
       width=120, 
       height=80, 
       scale=3, 
       dpi=300)



out.plot2 <- cowplot::plot_grid(pA, pB, nrow=1, ncol=2, labels=c("(a)", "(b)"), rel_widths = c(1.2,1.2))

print(out.plot2)

ggsave(file = paste0(homewd,"/figures/Fig1_TwoPanel_R1.jpeg"),
       units="mm",  
       width=170, 
       height=70, 
       scale=3, 
       dpi=300)







##########################################################  MRT. #################################


############################ MRT Cleaning ###########################
tree <- read.tree("timetree_names.nwk")

# clean dataset with only the columns you need
MRT <- dat_MRT[c("phylo_name", "MRT_hrs", "typical.diet", "re_class", "mass_kg")]
#MRT <- na.omit(MRT)

#make sure they match
# Check if the species names match the tree tips
common_species <- intersect(tree$tip.label, MRT$phylo_name)
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, common_species))
MRT_pruned <- MRT[MRT$phylo_name %in% common_species, ]

# check to see that the names match
# Species in common_species but not in tree_pruned$tip.label
# missing_from_tree <- setdiff(common_species, tree_pruned$tip.label)
# # # Print the missing species
# missing_from_tree #there are zero


species1 <- unique(dat_MRT$phylo_name)
species2 <- unique(MRT_pruned$phylo_name)

# Species in df1 but not in df2
missing_in_df2 <- setdiff(species1, species2)
missing_in_df2
#[1] "Cercopithecus_mitis"    "Cercopithecus_pogonias"



# what is in the database
unique(MRT_pruned$re_class)
# [1] "Non-Flying Birds"    "Ungulates"        "Primates"         "Rodents"          "Bats" 

# factor the groups so that rodent is the comparasion
MRT_pruned$re_class <- factor(MRT_pruned$re_class, levels = c("Rodents","Non-Flying Birds","Bats", "Primates", "Ungulates"))

#transform data
MRT_pruned$log_MRT_hrs <- log10(MRT_pruned$MRT_hrs)





############################# MRT - model 1 ###############

# Create correlation structure
cor_phylo_fixed1 <- corPagel(1, phy = tree_pruned, fixed = TRUE, form = ~phylo_name)


pgls_model_MRT_1 <- gls(log_MRT_hrs ~ re_class, 
                        data = MRT_pruned, 
                        correlation = cor_phylo_fixed1,
                        method = "ML",
                        weights = varIdent(form = ~1|typical.diet))


#check that residuals are normal
hist(resid(pgls_model_MRT_1), main = "Residual Histogram", xlab = "Residuals")

# View the results
summary(pgls_model_MRT_1)

# Value Std.Error   t-value p-value
# (Intercept)               1.2752458 0.1578501  8.078841  0.0000
# re_classNon-Flying Birds -0.3729681 0.8072160 -0.462043  0.6459
# re_classBats             -0.7031048 0.3109098 -2.261443  0.0279
# re_classPrimates          0.2770018 0.1471707  1.882180  0.0653
# re_classUngulates         0.3745205 0.2716690  1.378591  0.1738




pgls_model_MRT_1_reduced <- gls(log_MRT_hrs ~ re_class, 
                        data = MRT_pruned, 
                        correlation = cor_phylo_fixed1,
                        method = "ML")

summary(pgls_model_MRT_1_reduced)

AIC(pgls_model_MRT_1,pgls_model_MRT_1_reduced)
# df       AIC
# pgls_model_MRT_1          9 -46.87001
# pgls_model_MRT_1_reduced  6 -39.26213

# delta AIC = 7.60788



######### MRT model 2 - mass ######


pgls_model_MRT_2 <- gls(log_MRT_hrs ~ log10(mass_kg):re_class, 
                        data = MRT_pruned, 
                        correlation = cor_phylo_fixed1,
                        method = "ML",
                        weights = varIdent(form = ~1|typical.diet))


#check that residuals are normal
hist(resid(pgls_model_MRT_2), main = "Residual Histogram", xlab = "Residuals")

# View the results
summary(pgls_model_MRT_2)

# Value Std.Error  t-value p-value
# (Intercept)                              1.4935506 0.0348951 42.80118  0.0000
# log10(mass_kg):re_classRodents           0.1631381 0.0263710  6.18627  0.0000
# log10(mass_kg):re_classNon-Flying Birds -0.4856223 0.5635864 -0.86166  0.3928
# log10(mass_kg):re_classBats              0.1682470 0.0801067  2.10029  0.0406
# log10(mass_kg):re_classPrimates         -0.1724393 0.0772188 -2.23313  0.0299
# log10(mass_kg):re_classUngulates        -0.0376446 0.0332968 -1.13058  0.2634



pgls_model_MRT_2_reduced <- gls(log_MRT_hrs ~ log10(mass_kg):re_class, 
                        data = MRT_pruned, 
                        correlation = cor_phylo_fixed1,
                        method = "ML")

AIC(pgls_model_MRT_2, pgls_model_MRT_2_reduced)
# df       AIC
# pgls_model_MRT_2         10 -62.07957
# pgls_model_MRT_2_reduced  7 -40.68041

# delta AIC = 21.39916
