###Phylo analysis
# load the packages
library(lme4)
library(nlme)
library(ape)
library(dplyr)
library(phytools)
library(ggplot2)
library(ggforce)

# load data
#homewd <- "/Users/emilyruhs/Desktop/UChi_Brook_Lab/GitHub_repos/git-transit-time/"
homewd <- "/Users/katherinemcferrin/Developer/git-transit-time/"

#load the GIT transit data with phylo name
dat <- read.csv(file = paste0(homewd, "/data/revision-1/dat_sum_tot_clean.csv"), header = T, stringsAsFactors = F )

#load the tree
tree <- read.tree(file = paste0(homewd, "/data/revision-1/Genus_species_list_2.nwk"))



##### plot the phylogeny for GIT and MRT #######
# plot(tree,
#      cex = 0.6,         # Tip label size
#      no.margin = FALSE,  # Remove extra white space
#      main = "Phylogenetic Tree")
# 
# # Add a scale bar (optional)
# add.scale.bar()
# 
# 
# 
# 
# # Ladderize for cleaner visualization
# tree <- ladderize(tree)
# 
# # Basic styled plot
# plot(tree,
#      show.tip.label = TRUE,   # Show tip labels (default is TRUE, but made explicit)
#      cex = 0.6,               # Font size for tip labels
#      main = "Phylogenetic Tree with Tip Labels",
#      no.margin = TRUE, 
#      edge.width = 1.5)
# 
# # Add a scale bar
# add.scale.bar()







 
####transit time phylogenetic signal ####

#for GIT
GIT <- dat[c("phylo_name", "transit_hrs")]
GIT <- na.omit(GIT)

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
# Phylogenetic signal lambda : 0.458034 
# logL(lambda) : -688.179 
# LR(lambda=0) : 21.2645 
# P-value (based on LR test) : 4.00066e-06 



contMap(tree,log10(transit_hrs), fsize =0.7)
# save as MGT_contmap, 550 width 1000 height

#phenogram(tree, log10(transit_hrs), spread.labels = T, fsize=0.4, colors = "black",sublabel.angle=45)
# save as MGT_phenogram





################### MRT phylogenetic signal ##########################
#for MRT
MRT <- dat[c("phylo_name", "MRT_hrs")]
MRT <- na.omit(MRT)

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
# Phylogenetic signal lambda : 0.999927 
# logL(lambda) : -173.223 
# LR(lambda=0) : 77.0561 
# P-value (based on LR test) : 1.66168e-18 



############# plotting phylo signal by tree/phenogram ############


contMap(tree,log10(MRT_hrs), fsize =0.7)


#phenogram(tree, log10(MRT_hrs), spread.labels = T, fsize=0.6, colors = "black")


#########################################  MGT. #################################

########### Transit time data cleaning #######

#_________________________________________________________________________________________________________________
#I think this code runs PGLS and also provides lambda estimates as part of the model
library(caper)

#caper automatically builds the correlation structure using the comparative data object
# Load your tree
tree <- read.tree(file = paste0(homewd, "/data/revision-1/Genus_species_list_2.nwk"))

# Load your data
dat <- read.csv(file = paste0(homewd, "/data/revision-1/dat_sum_tot_clean.csv"), header = T, stringsAsFactors = F )
names(dat)

# clean dataset with only the columns you need
GIT <- dat[c("phylo_name", "transit_hrs", "typical.diet", "re_class", "mass_kg")]
GIT <- na.omit(GIT)

#make sure they match
# Check if the species names match the tree tips
common_species <- intersect(tree$tip.label, GIT$phylo_name)
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, common_species))
GIT_pruned <- GIT[GIT$phylo_name %in% common_species, ]

# check to see that the names match
# Species in common_species but not in tree_pruned$tip.label
missing_from_tree <- setdiff(common_species, tree_pruned$tip.label)
# Print the missing species
missing_from_tree #there are zero



###################### Transit time - model 1 ################


# Fit PGLS model
#ML=maximum likelihood

# factor the groups so that rodent is the comparasion
unique(GIT_pruned$re_class)

GIT_pruned$re_class <- factor(GIT_pruned$re_class, levels = c("Rodents","Flying Birds","Non-Flying Birds","Bats", "Carnivores", "Primates", "Ungulates", "Reptiles"))

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
# (Intercept)               1.1051483 0.3253658  3.396633  0.0010
# re_classFlying Birds     -0.9278879 0.4952769 -1.873473  0.0643
# re_classNon-Flying Birds -0.3183586 0.6617106 -0.481115  0.6316
# re_classBats             -1.0078472 0.3764212 -2.677446  0.0088
# re_classCarnivores        0.0504731 0.5014241  0.100659  0.9200
# re_classPrimates          0.1524649 0.3526723  0.432313  0.6666
# re_classUngulates         0.4928203 0.5106141  0.965152  0.3371
# re_classReptiles          0.8950368 0.6314527  1.417425  0.1599


pgls_model_GIT_1_reduced <- gls(log_transit_hrs ~ re_class, 
                                data = GIT_pruned, 
                                correlation = cor_phylo_fixed1,
                                method = "ML")


AIC(pgls_model_GIT_1, pgls_model_GIT_1_reduced)
# df      AIC
# pgls_model_GIT_1         13 147.6720
# pgls_model_GIT_1_reduced  9 164.0248





############### Plot MGT - model 1 ###############

# Add residuals and fitted values to your data
GIT_pruned$resid <- resid(pgls_model_GIT_1, type = "normalized")
GIT_pruned$fitted <- fitted(pgls_model_GIT_1)

# Basic plot
# library(ggplot2)
#shapez = c("carnivore" = 21, "frugivore/nectarivore" = 22, "folivore" = 23, "insectivore" = 24, "omnivore" = 25)

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

D

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

# ggplot(var_weights_df, aes(x = reorder(typical_diet, -relative_variance), y = relative_variance)) +
#   geom_col(fill = "steelblue") +
#   coord_flip() +
#   labs(title = "Estimated Relative Variance by Diet Type",
#        x = "Typical Diet",
#        y = "Relative Residual Variance") +
#   theme_minimal()

C <- ggplot(GIT_pruned, aes(x = typical.diet, y = resid)) +
  geom_boxplot() +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x=element_text(size=11),
        axis.text.y=element_text(size=11))+
  labs(title = "Residual Spread by Natural Diet",
       y = "\nNormalized Residuals",
       x="")
C

ggsave(file = paste0(homewd,"/figures/MGT_resid.jpeg"),
       plot=C,
       units="mm",  
       width=70, 
       height=70, 
       scale=3, 
       dpi=300)


# plot(fitted(pgls_model_GIT_1), residuals(pgls_model_GIT_1),
#      xlab = "Fitted values",
#      ylab = "Residuals",
#      main="Residuals vs. Fitted")
# abline(h=0,col='red')
# 
# plot(GIT_pruned$log_transit_hrs, fitted(pgls_model_GIT_1),
#      xlab = "Observed",
#      ylab = "Fitted",
#      main = "Observed vs. Fitted")
# abline(0,1,col="blue", lwd=2)






#and reshape order
GIT_pruned$re_class <- as.character(GIT_pruned$re_class)
GIT_pruned$re_class[GIT_pruned$re_class=="Flying Birds"] <- "Flying\nBirds"
GIT_pruned$re_class[GIT_pruned$re_class=="Non-Flying Birds"] <- "Non-Flying\nBirds"
#dat.sum.tot_clean$re_class[dat.sum.tot_clean$re_class=="Non-Flying Birds"] <- "Non-Flying\nBirds"
GIT_pruned$re_class <- factor(GIT_pruned$re_class, levels = c( "Flying\nBirds","Bats","Rodents", "Non-Flying\nBirds", "Carnivores","Primates", "Ungulates", "Reptiles"))
#GIT_pruned$typical_diet = factor(GIT_pruned$typical_diet, levels= c("ommnivore", "insectivore", "frugivore/nectarivore", "folivore", "carnivore")) #this is messing up and adding NAs



#scales::show_col(scales::hue_pal()(9))
colz =  c( "Flying\nBirds"= "#F8766D", "Non-Flying\nBirds"="#edf8b1","Bats" = "#C49A00", "Rodents"= "navy","Carnivores" = "#00C094", "Primates"= "#00B6EB", "Ungulates"= "#A58AFF",   "Reptiles" = "#FB61D7")

#shapez = c("fiber/foliage" = 21, "fruit/nectar/pollen" = 22, "label-solution" = 23, "protein" = 24)
shapez = c("carnivore" = 21, "frugivore/nectarivore" = 22, "folivore" = 23, "insectivore" = 24, "omnivore" = 25)


#and the version that is 2panel
label.dat <- cbind.data.frame(re_class=c("Flying\nBirds","Bats", "Rodents","Non-Flying\nBirds", "Carnivores", "Primates", "Ungulates", "Reptiles"), label = c("", "**", "", "","","","","**"))

label.dat$re_class <- factor(label.dat$re_class, levels=c("Flying\nBirds","Bats","Rodents","Non-Flying\nBirds", "Carnivores", "Primates", "Ungulates", "Reptiles"))


#and all together
pA <- ggplot(data=GIT_pruned)+ geom_boxplot(aes(x=re_class, y=log10(transit_hrs), fill=re_class), show.legend = F) + theme_bw() + 
  theme(axis.title.x = element_blank(), panel.grid = element_blank(), legend.position = c(.85,.15), legend.background = element_rect(color="black"),
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=12), axis.title.y = element_text(size=15),
        plot.margin = unit(c(.2,.2,.8,.2), "cm")) + scale_fill_manual(values=colz) + 
  geom_point(aes(x=re_class, y=log10(transit_hrs),shape=typical.diet), position = position_jitterdodge(jitter.width = 0), size=3) +
  scale_shape_manual(values=shapez, name="Natural diet") +  ylab("GIT transit time, hrs") +
  geom_hline(yintercept=0.6684442, linetype=2) +
  scale_y_continuous(breaks=c(0,1,2), labels = c("1", "10", "100")) +
  coord_cartesian(ylim=c(-.8,2.95)) + geom_label(data=label.dat, aes(x=re_class, y=3.0, label=label), label.size = NA,size=5) 
  #scale_x_discrete(labels = c("Flying\nBirds","Bats","Non-Flying\nBirds", "Rodents", "Carnivores", "Marsupials","Lagomorphs","Primates", "Ungulates", "Reptiles"))
print(pA)


ggsave(file = paste0(homewd,"/figures/Fig_1A.png"),
       units="mm",  
       width=120, 
       height=80, 
       scale=3, 
       dpi=300)




###################### Transit time - model 2 - mass ################
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
# (Intercept)                               1.1098810 0.1775645  6.250580  0.0000
# log10(mass_kg):re_classRodents            0.3419527 0.0754266  4.533580  0.0000
# log10(mass_kg):re_classFlying\nBirds      0.5165710 0.1558750  3.314008  0.0013
# log10(mass_kg):re_classBats               0.3139096 0.0789419  3.976464  0.0001
# log10(mass_kg):re_classNon-Flying\nBirds -0.3109468 0.4828060 -0.644041  0.5212
# log10(mass_kg):re_classCarnivores         0.7227498 0.6442551  1.121838  0.2650
# log10(mass_kg):re_classPrimates           0.1503524 0.1047483  1.435368  0.1548

# log10(mass_kg):re_classUngulates         -0.2114904 0.1149872 -1.839251  0.0693
# log10(mass_kg):re_classReptiles           0.1809061 0.1635589  1.106060  0.2717

pgls_model_GIT_2_reduced <- gls(log_transit_hrs ~ log10(mass_kg):re_class, 
                        data = GIT_pruned, 
                        correlation = cor_phylo_fixed1,
                        method = "ML")

AIC(pgls_model_GIT_2,pgls_model_GIT_2_reduced)
# df      AIC
# pgls_model_GIT_2         14 135.6655
# pgls_model_GIT_2_reduced 10 153.1407


#the trend lines in panel b are based off of this linear mixed effect model (not the PGLS)
m2c <- lmer(log10(transit_hrs)~ log10(mass_kg):re_class + (1|re_class), data=GIT_pruned)

GIT_pruned$predicted_transit[!is.na(GIT_pruned$transit_hrs) & !is.na(GIT_pruned$mass_kg) & !is.na(GIT_pruned$re_class)]   <- 10^(predict(m2c))
summary(m2c)




colz =  c( "Flying Birds"= "#F8766D", "Bats" = "#C49A00", "Non-Flying Birds"="#edf8b1", "Rodents"= "navy","Carnivores" = "#00C094","Primates"= "#00B6EB", "Ungulates"= "#A58AFF",   "Reptiles" = "#FB61D7")

#shapez = c("fiber/foliage" = 21, "fruit/nectar/pollen" = 22, "label-solution" = 23, "protein" = 24)
shapez = c("carnivore" = 21, "frugivore/nectarivore" = 22, "folivore" = 23, "insectivore" = 24, "omnivore" = 25)


pB <- ggplot(data = GIT_pruned, aes(x = log10(mass_kg), y = log10(transit_hrs))) + 
  geom_mark_ellipse(expand = 0, radius = 0, aes(fill = re_class, color = re_class), size = 0.1) +
  geom_line(aes(x = log10(mass_kg), y = log10(predicted_transit), color = re_class, linetype = re_class), 
            alpha = 0.6, size = 0.7, show.legend = TRUE) +
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
# clean dataset with only the columns you need
MRT <- dat[c("phylo_name", "MRT_hrs", "typical.diet", "re_class", "mass_kg")]
MRT <- na.omit(MRT)

#make sure they match
# Check if the species names match the tree tips
common_species <- intersect(tree$tip.label, MRT$phylo_name)
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, common_species))
MRT_pruned <- MRT[MRT$phylo_name %in% common_species, ]

# check to see that the names match
# Species in common_species but not in tree_pruned$tip.label
missing_from_tree <- setdiff(common_species, tree_pruned$tip.label)
# Print the missing species
missing_from_tree #there are zero


# factor the groups so that rodent is the comparasion
MRT_pruned$re_class <- factor(MRT_pruned$re_class, levels = c("Rodents","Flying Birds","Non-Flying Birds","Bats", "Carnivores", "Primates", "Ungulates", "Reptiles"))

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
# (Intercept)               1.0706518 0.0964935 11.095581  0.0000
# re_classNon-Flying Birds  0.2283571 0.3541948  0.644722  0.5230
# re_classBats             -1.1232313 0.3516257 -3.194395  0.0028
# re_classPrimates          0.1279635 0.1454538  0.879754  0.3845
# re_classUngulates         0.3384721 0.4637528  0.729855  0.4700
# re_classReptiles          1.2157485 0.3541948  3.432429  0.0015


pgls_model_MRT_1_reduced <- gls(log_MRT_hrs ~ re_class, 
                        data = MRT_pruned, 
                        correlation = cor_phylo_fixed1,
                        method = "ML")

summary(pgls_model_MRT_1_reduced)

AIC(pgls_model_MRT_1,pgls_model_MRT_1_reduced)
# df       AIC
# pgls_model_MRT_1         10  9.998707
# pgls_model_MRT_1_reduced  7 14.425907





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
# (Intercept)                              1.249133 0.0357769 34.91452  0.0000
# log10(mass_kg):re_classRodents           0.118489 0.0259190  4.57151  0.0001
# log10(mass_kg):re_classNon-Flying Birds  0.043879 0.1308871  0.33524  0.7393
# log10(mass_kg):re_classBats              0.108696 0.1281041  0.84849  0.4016
# log10(mass_kg):re_classPrimates          0.155993 0.0266529  5.85275  0.0000
# log10(mass_kg):re_classUngulates        -0.050290 0.1127025 -0.44622  0.6580
# log10(mass_kg):re_classReptiles         -3.445726 0.4942224 -6.97201  0.0000

pgls_model_MRT_2_reduced <- gls(log_MRT_hrs ~ log10(mass_kg):re_class, 
                        data = MRT_pruned, 
                        correlation = cor_phylo_fixed1,
                        method = "ML")

AIC(pgls_model_MRT_2, pgls_model_MRT_2_reduced)
# df       AIC
# pgls_model_MRT_2         11  5.246758
# pgls_model_MRT_2_reduced  8 20.434000


