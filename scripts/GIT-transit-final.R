rm(list=ls())

#install.packages("TMB", type = "source")
#install.packages("sjPlot")
#install.packages("plyr")
#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("lmerTest")
#install.packages("ggforce")

library(plyr)
library(dplyr)
library(ggplot2)
library(lmerTest)
library(sjPlot)
library(ggforce)
library(glmmTMB)

#set home directory
#homewd= "/Users/carabrook/Developer/git-transit-time"
#homewd= "/Users/katherinemcferrin/Developer/git-transit-time"
homewd <- "/Users/emilyruhs/Desktop/UChi_Brook_Lab/GitHub_repos/git-transit-time/"


#load the GIT transit data:
#dat <- read.csv(file = paste0(homewd, "/data/final-GIT-transit-database-april-2021.csv"), header = T, stringsAsFactors = F )
dat <- read.csv(file = paste0(homewd, "/data/McFerrin_database.csv"), header = T, stringsAsFactors = F )
#View(dat)

length(unique(dat$Retention.Citation)) #107 unique papers
length(unique(dat$Genus.species)) #112 unique species
unique(dat$Genus.species)


#dat = only the columns we will use for modeling and plotting and renaming columns
#original code: dat <- dplyr::select(dat, can_fly, Class, Order, Family, Genus.species, Common.Name, Typical.Diet, Mass..mean.median.from.study., Feeding.Trial.Food.Cat, N_measured, N_trials, Minimum..min., Median..min., Mean..min.,Maximum..min.)
dat <- dplyr::select(dat, can_fly, Class, Order, Family, Genus.species, Common.Name, Typical.Diet, Mass.g..mean.median.from.study., Feeding.Trial.Food.Cat, N_individuals, N_trials, Minimum..min., Median..min., Mean..min.,Maximum..min.)
head(dat)
#rename columns
names(dat) <- c("fly", "class", "order", "family", "genus.species", "common.name", "typical.diet", "mass", "food.cat", "N_individuals", "N_trials", "min", "median", "mean", "max")
#View(dat)


#choose one mean/median to report
dat.plot <- dat
#View(dat.plot)
dat.plot$transit <- dat.plot$mean #copying the mean git transit values to new column transit
dat.plot$transit[is.na(dat.plot$transit)]<- dat.plot$median[is.na(dat.plot$transit)] #if transit doesn't have a mean, fill it with median

dat.plot <- dplyr::select(dat.plot, -(min), -(median), -(mean), -(max)) #simplifying by taking out min, median, mean and max
head(dat.plot)

#summarize by species
dat.plot$N_individuals[dat.plot$N_individuals=="not reported"] <- 1 #assigning 1 individual if total # not reported
dat.plot$N_individuals[is.na(dat.plot$N_individuals)] <- 1 #assigning 1 individual if #individuals is NA
dat.plot$N_individuals = as.numeric(dat.plot$N_individuals)
dat.plot$total_transit = dat.plot$transit*dat.plot$N_individuals #multiplying transit time x #individuals
unique(dat.plot$food.cat)
#View(dat.plot)


#now get one entry for each species/food category combination
dat.split <- dlply(dat.plot, .(fly,class,order,family,genus.species, common.name, typical.diet,food.cat))
#making a function
summarise.dat <- function(dat1){
  
  dat2 <- ddply(dat1, .(fly,class,order,family,genus.species, common.name, typical.diet,food.cat), summarise, mass=unique(mass), food=unique(food.cat),  N_tot =sum(N_individuals), total_transit=sum(total_transit))
  return(dat2)
}
dat.out <- lapply(dat.split, summarise.dat)
dat.sum.tot <- data.table::rbindlist(dat.out)
head(dat.sum.tot)
dat.sum.tot$transit <- dat.sum.tot$total_transit/dat.sum.tot$N_tot
subset(dat.sum.tot, is.na(transit)) #0

#now categorize -- vertebrates are signified by class but mammals by order
dat.sum.tot$re_class <- dat.sum.tot$class
dat.sum.tot$re_class[dat.sum.tot$re_class=="Mammalia"] <- dat.sum.tot$order[dat.sum.tot$re_class=="Mammalia"]
dat.sum.tot$re_class[dat.sum.tot$re_class=="Aves" & dat.sum.tot$fly=="Y"] <- "Flying Birds"
dat.sum.tot$re_class[dat.sum.tot$re_class=="Aves" & dat.sum.tot$fly=="N"] <- "Non-Flying Birds"

unique(dat.sum.tot$re_class)
#find out how many entries per each "cat"
dat.simp.big <- ddply(dat.sum.tot, .(re_class, food.cat), summarize, N=length(re_class))
dat.simp <- ddply(dat.sum.tot, .(re_class), summarize, N=length(re_class))
dat.simp

#and group
#remove any non-mammalian classes and any mammalian orders with < 4 entries
#Amphibia, Chondrichtythes, Dermoptera, Lagomorpha, Non-Flying Birds 
#Peramelemorphia, Sirenia
dat.sum.tot = subset(dat.sum.tot, class!="Chondrichthyes" & class !="Amphibia")
#and write over some of the others
dat.sum.tot$re_class[dat.sum.tot$re_class=="Reptilia"] <- "Reptiles"
#dat.sum.tot$re_class[dat.sum.tot$re_class=="Lagomorpha"] <- "Rodents and Lagomorphs"
dat.sum.tot$re_class[dat.sum.tot$re_class=="Rodentia"] <- "Rodents"
dat.sum.tot$re_class[dat.sum.tot$re_class=="Chiroptera"] <- "Bats"
#dat.sum.tot$re_class[dat.sum.tot$re_class=="Artiodactyla"] <- "Ungulates"
dat.sum.tot$re_class[dat.sum.tot$re_class=="Carnivora"] <- "Carnivores"
dat.sum.tot$re_class[dat.sum.tot$re_class=="Perissodactyla" | dat.sum.tot$re_class=="Artiodactyla"] <- "Ungulates"
dat.sum.tot$re_class[dat.sum.tot$re_class=="Lagomorpha" | dat.sum.tot$re_class=="Eulipotyphla" | dat.sum.tot$re_class=="Peramelemorphia" | dat.sum.tot$re_class=="Sirenia" |  dat.sum.tot$re_class=="Dermoptera" ] <- "Other Mammals"

dat.sum.tot = subset(dat.sum.tot, re_class!="Non-Flying Birds" & re_class!= "Other Mammals")
#dat.sum.tot = subset(dat.sum.tot,  re_class!= "Other Mammals")
unique(dat.sum.tot$re_class)

#convert transit time to hours
dat.sum.tot$transit_hrs = dat.sum.tot$transit/60
#convert mass to kg
dat.sum.tot$mass_kg = as.numeric(dat.sum.tot$mass)/1000

#check how transit time varies by mass
#generally, we see that transit is longer at higher mass
#all taxa agree with this trend, though reptiles and bats are a little longer in time
#than we might otherwise predict -- this is likely due to an effect of temperature that is not
#reported here.

#and check how transit varies by food category
#taxa appear to have a stronger effect than food type
#reclass.colz =  c("Rodents and Lagomorphs" = "black", "Bats"= "black", "Flying Birds"= "black", "Reptiles"= "black",  "Ungulates"= "black", "Carnivores"= "black", "Primates"= "black")
#reclass.colz =  c("Rodents" = "black", "Bats"= "black", "Flying Birds"= "black", "Reptiles"= "black",  "Ungulates"= "black", "Carnivores"= "black", "Primates"= "black")
#and replot the mass by transit with categories
dat.sum.tot$re_class <- factor(dat.sum.tot$re_class, levels = c( "Flying Birds","Bats", "Rodents", "Carnivores", "Primates", "Ungulates",   "Reptiles"))


#now some stats
dat.sum.tot$re_class = factor(dat.sum.tot$re_class, levels= c( "Rodents",  "Flying Birds","Bats","Carnivores", "Primates", "Ungulates",  "Reptiles"))
dat.sum.tot$food.cat = factor(dat.sum.tot$food.cat, levels= c("label-solution", "protein", "fiber/foliage", "fruit/nectar/pollen"))

m1 <- lmer(log10(transit_hrs)~ re_class + (1|food.cat), data=dat.sum.tot, REML = F)
summary(m1)
# Fixed effects:
#                         Estimate  Std. Error       df     t value     Pr(>|t|)    
# (Intercept)            0.5473     0.1353      40.5218     4.045     0.000228 ***
# re_classFlying Birds  -0.5873     0.1492      139.7038    -3.936    0.000130 ***
# re_classBats          -0.4273     0.1387      135.8771    -3.080    0.002505 ** 
# re_classCarnivores     0.4952     0.1635      133.0883     3.028    0.002955 ** 
# re_classPrimates       0.8392     0.1435      139.8450    5.849     3.34e-08 ***
# re_classUngulates      1.0173     0.1487      139.5672    6.840     2.29e-10 ***
# re_classReptiles       1.4621     0.1488      139.9711    9.827     < 2e-16 ***
rand(m1) 
# Model:
#   log10(transit_hrs) ~ re_class + (1 | food.cat)
# npar  logLik    AIC    LRT Df Pr(>Chisq)  
# <none>            9 -48.649 115.30                       
# (1 | food.cat)    8 -51.205 118.41 5.1113  1    0.02377 *

#for plotting later
m1b <- lme4::lmer(log10(transit_hrs)~ re_class + (1|food.cat), data=dat.sum.tot, REML = F)
summary(m1b)
# Fixed effects:
# Estimate Std. Error t value
# (Intercept)            0.5473     0.1353   4.045
# re_classFlying Birds  -0.5873     0.1492  -3.936
# re_classBats          -0.4273     0.1387  -3.080
# re_classCarnivores     0.4952     0.1635   3.028
# re_classPrimates       0.8392     0.1435   5.849
# re_classUngulates      1.0173     0.1487   6.840
# re_classReptiles       1.4621     0.1488   9.827


#and here plot the raw partial effects of taxonomic grouping for Figure S2:


#plot fixed effects
pS2 <- plot_model(m1b, type="est", vline.color = "black",
                  axis.labels = rev(c("Flying Birds", "Bats", "Carnivores", "Primates", "Ungulates", "Reptiles")),
                  title = "Effects of taxon on GIT transit") + theme_bw() + 
  theme(panel.grid = element_blank(), axis.text = element_text(size=14), axis.title.x = element_blank()) 

print(pS2)

ggsave(file = paste0(homewd,"/figures/FigS2_Taxon_effects.png"),
       units="mm",  
       width=50, 
       height=40, 
       scale=3, 
       dpi=300)

m2 <- lmer(log10(transit_hrs)~ log10(mass_kg):re_class + (1|re_class) +(1|food.cat), data=dat.sum.tot, REML = F)
# Fixed effects:
# Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)                           0.89560    0.25182   7.31854   3.557  0.00861 ** 
# log10(mass_kg):re_classRodents        0.25758    0.11039 133.56186   2.333  0.02112 *  
# log10(mass_kg):re_classFlying Birds   0.48209    0.08776  81.20270   5.493 4.42e-07 ***
# log10(mass_kg):re_classBats          -0.12193    0.09147  88.65304  -1.333  0.18594    
# log10(mass_kg):re_classCarnivores     0.49250    0.19771 136.26238   2.491  0.01394 *  
# log10(mass_kg):re_classPrimates       0.28920    0.13624 133.33183   2.123  0.03563 *  
# log10(mass_kg):re_classUngulates      0.04259    0.09621 135.44017   0.443  0.65867    
# log10(mass_kg):re_classReptiles       0.13418    0.06606 132.48521   2.031  0.04423 *  
m2 <- lmer(log10(transit_hrs)~ log10(mass_kg):re_class + (1|re_class), data=dat.sum.tot, REML = F)
summary(m2)
# Fixed effects:
# Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)                           0.89115    0.25511   7.13648   3.493  0.00978 ** 
# log10(mass_kg):re_classFlying Birds   0.50171    0.08364 136.43675   5.998 1.69e-08 ***
# log10(mass_kg):re_classBats          -0.16267    0.08661 139.23271  -1.878  0.06244 . 
# log10(mass_kg):re_classRodents        0.27406    0.11004 133.81241   2.491  0.01397 *  
# log10(mass_kg):re_classCarnivores     0.48411    0.19952 139.91864   2.426  0.01652 *  
# log10(mass_kg):re_classPrimates       0.29365    0.13470 137.56043   2.180  0.03095 *  
# log10(mass_kg):re_classUngulates      0.03099    0.09638 135.84425   0.322  0.74826    
# log10(mass_kg):re_classReptiles       0.12151    0.06620 132.74398   1.836  0.06865 .
rand(m2) 
# Model:
# log10(transit_hrs) ~ (1 | re_class) + log10(mass_kg):re_class
# npar   logLik    AIC    LRT Df Pr(>Chisq)    
# <none>           10  -44.803 109.61                         
# (1 | re_class)    9 -108.287 234.57 126.97  1  < 2.2e-16 ***

m2b <- lme4::lmer(log10(transit_hrs)~ log10(mass_kg):re_class + (1|re_class)+(1|food.cat), data=dat.sum.tot, REML = F)
summary(m2b)

pSm2 <- plot_model(m2b, type="est", vline.color = "black",
                   axis.labels = rev(c("Rodents", "Flying Birds", "Bats", "Carnivores", "Primates", "Ungulates", "Reptiles")),
                   title = "Effects of mass on GIT transit") + theme_bw() + 
  theme(panel.grid = element_blank(), axis.text = element_text(size=14), axis.title.x = element_blank()) 
print(pSm2)

library(gridExtra)
pSm2b <- plot_model(m2b, type="re", vline.color = "black", facet.grid=FALSE) + theme_bw() + 
  theme(panel.grid = element_blank(), axis.text = element_text(size=14)) 
grid.arrange(pSm2b[[1]], pSm2b[[2]])

mass.mod.plot <- cowplot::plot_grid(pSm2,pSm2b, nrow=1, ncol=2, labels=c("(a)", "(b)"), rel_widths = c(1,1))
print(mass.mod.plot)




m3 <- lmer(log10(transit_hrs)~ log10(mass_kg):re_class + (1|food.cat), data=dat.sum.tot, REML = F)
summary(m3)
# Fixed effects:
# Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)                           1.13781    0.12686   6.41860   8.969 7.29e-05 ***
# log10(mass_kg):re_classRodents        0.36725    0.16867 137.36160   2.177  0.03117 *  
# log10(mass_kg):re_classFlying Birds   0.78490    0.10169 139.80400   7.718 2.03e-12 ***
# log10(mass_kg):re_classBats           0.61208    0.06877 138.23840   8.900 2.80e-15 ***
# log10(mass_kg):re_classCarnivores    -0.19167    0.18480 139.78473  -1.037  0.30145    
# log10(mass_kg):re_classPrimates       0.40018    0.14368 139.98917   2.785  0.00609 ** 
# log10(mass_kg):re_classUngulates      0.19175    0.06816 135.12359   2.813  0.00564 ** 
# log10(mass_kg):re_classReptiles       0.30837    0.10756 136.36347   2.867  0.00480 ** 
rand(m3)
# Model:
# log10(transit_hrs) ~ (1 | food.cat) + log10(mass_kg):re_class
# npar  logLik    AIC    LRT Df Pr(>Chisq)    
# <none>           10 -101.93 223.87                         
# (1 | food.cat)    9 -108.29 234.57 12.709  1  0.0003639 ***

m3b <- lme4::lmer(log10(transit_hrs)~ log10(mass_kg):re_class + (1|food.cat), data=dat.sum.tot, REML = F)
summary(m3b)

pSm3 <- plot_model(m3b, type="est", vline.color = "black",
                  axis.labels = rev(c("Rodents", "Flying Birds", "Bats", "Carnivores", "Primates", "Ungulates", "Reptiles")),
                  title = "Effects of mass on GIT transit") + theme_bw() + 
  theme(panel.grid = element_blank(), axis.text = element_text(size=14), axis.title.x = element_blank()) 
print(pSm3)

pSm3b <- plot_model(m3b, type="re", vline.color = "black") + theme_bw() + 
  theme(panel.grid = element_blank(), axis.text = element_text(size=14)) 
print(pSm3b)

mass.mod.plot <- cowplot::plot_grid(pSm3,pSm3b, nrow=1, ncol=2, labels=c("(a)", "(b)"), rel_widths = c(1,1))
print(mass.mod.plot)

ggsave(file = paste0(homewd,"/figures/Fig3_Randomeffects.png"),
       units="mm",  
       width=120, 
       height=50, 
       scale=3, 
       dpi=300)




dat.sum.tot$predicted_transit[!is.na(dat.sum.tot$transit_hrs) & !is.na(dat.sum.tot$mass_kg) & !is.na(dat.sum.tot$re_class)]   <- 10^(predict(m2))



#and reshape order
dat.sum.tot$re_class <- as.character(dat.sum.tot$re_class)
dat.sum.tot$re_class[dat.sum.tot$re_class=="Flying Birds"] <- "Flying\nBirds"
dat.sum.tot$re_class <- factor(dat.sum.tot$re_class, levels = c( "Flying\nBirds","Bats", "Rodents", "Carnivores", "Ungulates", "Primates",  "Reptiles"))
dat.sum.tot$food.cat = factor(dat.sum.tot$food.cat, levels= c("fiber/foliage", "fruit/nectar/pollen","protein","label-solution"))


scales::show_col(scales::hue_pal()(7))
colz =  c( "Flying\nBirds"= "#F8766D", "Bats" = "#C49A00", "Rodents"= "navy","Carnivores" = "#00C094",  "Primates"= "#00B6EB", "Ungulates"= "#A58AFF",   "Reptiles" = "#FB61D7")
shapez = c("fiber/foliage" = 21, "fruit/nectar/pollen" = 22, "label-solution" = 23, "protein" = 24)


#and the version that is 2panel
label.dat <- cbind.data.frame(re_class=c("Flying\nBirds","Bats", "Rodents", "Carnivores", "Primates", "Ungulates",   "Reptiles"), label = c("***", "**", "", "**", "***", "***", "***"))
label.dat$re_class <- factor(label.dat$re_class, levels=c( "Flying\nBirds","Bats", "Rodents", "Carnivores", "Primates", "Ungulates", "Reptiles"))

#and all together
pA <- ggplot(data=dat.sum.tot)+ geom_boxplot(aes(x=re_class, y=log10(transit_hrs), fill=re_class), show.legend = F) + theme_bw() + 
  theme(axis.title.x = element_blank(), panel.grid = element_blank(), legend.position = c(.85,.15), legend.background = element_rect(color="black"),
        axis.text.y = element_text(size=13), axis.text.x = element_text(size=13), axis.title.y = element_text(size=16),
        plot.margin = unit(c(.2,.2,.8,.2), "cm")) + scale_fill_manual(values=colz) + 
  geom_point(aes(x=re_class, y=log10(transit_hrs),shape=food.cat), position = position_jitterdodge(jitter.width = 0), size=3) +  
  scale_shape_manual(values=shapez, name="Food type") +  ylab("GIT transit time, hrs") +
  geom_hline(aes(yintercept=log10(quantile(subset(dat.sum.tot, re_class=="Rodents")$transit_hrs)["50%"])), linetype=2) +
  scale_y_continuous(breaks=c(0,1,2), labels = c("1", "10", "100")) +
  coord_cartesian(ylim=c(-.8,2.95)) + geom_label(data=label.dat, aes(x=re_class, y=2.9, label=label), label.size = NA,size=5)
print(pA)

# pA_poster <- ggplot(data=dat.sum.tot)+ geom_boxplot(aes(x=re_class, y=log10(transit_hrs), fill=re_class), show.legend = F) + theme_bw() + 
#   theme(axis.title.x = element_blank(), panel.grid = element_blank(), legend.position = c(.85,.15), 
#         legend.title = element_text(size=23), legend.text = element_text(size=20),legend.background = element_rect(color="black"),
#         axis.text.y = element_text(size=25), axis.text.x = element_text(size=25), axis.title.y = element_text(size=28),
#         plot.margin = unit(c(.2,.2,.8,.2), "cm")) + scale_fill_manual(values=colz) + 
#   geom_point(aes(x=re_class, y=log10(transit_hrs),shape=food.cat), position = position_jitterdodge(jitter.width = 0), size=4) +  
#   scale_shape_manual(values=shapez, name="Food Type") +  ylab("GIT Transit Time (hrs)") +
#   geom_hline(aes(yintercept=log10(quantile(subset(dat.sum.tot, re_class=="Rodents")$transit_hrs)["50%"])), linetype=2) +
#   scale_y_continuous(breaks=c(0,1,2), labels = c("1", "10", "100")) +
#   coord_cartesian(ylim=c(-.8,2.95)) + geom_label(data=label.dat, aes(x=re_class, y=2.9, label=label), label.size = NA,size=8)
# print(pA_poster)

ggsave(file = paste0(homewd,"/figures/Fig_transit_tax_poster.png"),
       units="mm",  
       width=120, 
       height=80, 
       scale=3, 
       dpi=300)



#now plot with mass
#bats are the only taxon for which the mass slope from model 2 is negative (and significantly so)

head(dat.sum.tot)

pB <- ggplot(data=dat.sum.tot, aes(x=log10(mass_kg), y=log10(transit_hrs))) + 
  geom_mark_ellipse(expand=0,radius=0,aes(fill=re_class, color=re_class), size=.1)+
  geom_line(aes(x=log10(mass_kg), y=log10(predicted_transit), color=re_class), alpha=.6, size=.7) +
  #geom_ribbon(aes(x=log10(mass_kg), ymin=predicted_transit_lci, ymax=predicted_transit_uci, fill=re_class), alpha=.1) +
  geom_point(aes(x=log10(mass_kg), y=log10(transit_hrs), fill= re_class, shape=food.cat), size=4, show.legend = F) + 
  geom_point(aes(x=log10(mass_kg), y=log10(transit_hrs), fill= re_class, color= re_class, shape=food.cat), size=3, show.legend = F) + 
  ylab("GIT transit time (hrs)") + xlab("Mass (kg)") + scale_shape_manual(name = "Food Type", values=shapez) + 
  scale_fill_manual(name= "taxonomic\nclass", values=colz) +
  scale_color_manual(name= "taxonomic\nclass", values=colz) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = c(.85,.2),
                     axis.text = element_text(size=13), axis.title = element_text(size=16),
                     plot.margin = unit(c(.2,.2,.2,.2), "cm")) +
  scale_y_continuous(breaks = c(0,1,2), labels= c("1", "10", "100")) +
  scale_x_continuous(breaks = c(-2,-1,0,1,2,3), labels = c(".01", ".1", "1", "10", "100", "1000")) +
  coord_cartesian(ylim=c(-.8,2.9)) +
  guides(color = "none", fill="none")
print(pB)

# pB_poster <- ggplot(data=dat.sum.tot, aes(x=log10(mass_kg), y=log10(transit_hrs))) + 
#   geom_mark_ellipse(expand=0,radius=0,aes(fill=re_class, color=re_class), size=.1)+
#   geom_line(aes(x=log10(mass_kg), y=log10(predicted_transit), color=re_class), alpha=.6, size=.7) +
#   #geom_ribbon(aes(x=log10(mass_kg), ymin=predicted_transit_lci, ymax=predicted_transit_uci, fill=re_class), alpha=.1) +
#   geom_point(aes(x=log10(mass_kg), y=log10(transit_hrs), fill= re_class, shape=food.cat), size=5, show.legend = F) + 
#   geom_point(aes(x=log10(mass_kg), y=log10(transit_hrs), fill= re_class, color= re_class, shape=food.cat), size=4, show.legend = F) + 
#   ylab("GIT Transit Time (hrs)") + xlab("Mass (kg)") + scale_shape_manual(name = "Food Type", values=shapez) + 
#   scale_fill_manual(name= "taxonomic\nclass", values=colz) +
#   scale_color_manual(name= "taxonomic\nclass", values=colz) +
#   theme_bw() + theme(panel.grid = element_blank(), legend.position = c(.85,.2),
#                      axis.text = element_text(size=25), axis.title = element_text(size=28),
#                      legend.title = element_text(size=23), legend.text = element_text(size=20),
#                      plot.margin = unit(c(.4,.4,.4,.4), "cm")) +
#   scale_y_continuous(breaks = c(0,1,2), labels= c("1", "10", "100")) +
#   scale_x_continuous(breaks = c(-2,-1,0,1,2,3), labels = c(".01", ".1", "1", "10", "100", "1000")) +
#   coord_cartesian(ylim=c(-.8,2.9)) +
#   guides(color = "none", fill="none")
#print(pB_poster)

ggsave(file = paste0(homewd,"/figures/Fig_transit_mass_tax_poster.png"),
       units="mm",  
       width=120, 
       height=80, 
       scale=3, 
       dpi=300)




#and the supplementary figures

#plot random effects
pS3 <- plot_model(m1b, type="re", vline.color = "black") + theme_bw() + 
        theme(panel.grid = element_blank(), axis.text = element_text(size=14)) 
        
print(pS3)

ggsave(file = paste0(homewd,"/figures/FigS3_Random_Effects.png"),
       units="mm",  
       width=50, 
       height=40, 
       scale=3, 
       dpi=300)


#combine partial effects into a two panel plot
out.plot <- cowplot::plot_grid(pS2,pS3, nrow=1, ncol=2, labels=c("(a)", "(b)"), rel_widths = c(1,1))

print(out.plot)


ggsave(file = paste0(homewd,"/figures/Fig2_TwoPanel.jpeg"),
       units="mm",  
       width=120, 
       height=50, 
       scale=3, 
       dpi=300)


fam.col = c("Molossidae" = "firebrick", "Mormoopidae" = "navy", "Phyllostomidae" = "darkmagenta", "Pteropodidae" = "forestgreen", "Vespertilionidae" = "lightcoral")

#and look at just bats - reverses the direction
pS4 <- ggplot(data=subset(dat.sum.tot, re_class=="Bats")) + 
  geom_point(aes(x=log10(mass), y=log10(transit), fill= family, shape=food.cat), size=4, show.legend = F) + 
  geom_point(aes(x=log10(mass), y=log10(transit), fill= family, color= family, shape=food.cat), size=3) + 
  ylab("GIT transit time (min)") + xlab("Mass (g)") + coord_cartesian(xlim=c(.5,3), ylim=c(.5,3))+
  theme_bw() + theme(panel.grid = element_blank(), #legend.position = c(.87,.83), 
                     #legend.background = element_rect(color="black"),
                     axis.text = element_text(size=14), axis.title = element_text(size=18)) +
  scale_shape_manual(values=shapez, name= "food type") + 
  scale_y_continuous(breaks = c(1,2,3), labels= c("10", "100", "1000")) + 
  scale_fill_manual(name="bat family", values=fam.col) + scale_color_manual(name="bat family", values=fam.col) +
  scale_x_continuous(breaks = c(1,2,3), labels = c("10", "100", "1000")) #+
  #guides(shape=FALSE)#, nrow=1, colour = guide_legend(nrow = 1))

print(pS4)
ggsave(file = paste0(homewd,"/figures/FigS4_bat_only.png"),
       units="mm",  
       width=80, 
       height=60, 
       scale=3, 
       dpi=300)

# pS4_poster <- ggplot(data=subset(dat.sum.tot, re_class=="Bats")) + 
#   geom_point(aes(x=log10(mass), y=log10(transit), fill= family, shape=food.cat), size=5, show.legend = F) + 
#   geom_point(aes(x=log10(mass), y=log10(transit), fill= family, color= family, shape=food.cat), size=4) + 
#   ylab("GIT Transit Time (min)") + xlab("Mass (g)") + coord_cartesian(xlim=c(.5,3), ylim=c(.5,3))+
#   theme_bw() + theme(panel.grid = element_blank(),
#                      axis.text = element_text(size=25), axis.title = element_text(size=28), 
#                      legend.title = element_text(size=23), legend.text = element_text(size=20)) +
#   scale_shape_manual(values=shapez, name= "Food Type") + 
#   scale_y_continuous(breaks = c(1,2,3), labels= c("10", "100", "1000")) + 
#   scale_fill_manual(name="Bat Family", values=fam.col) + scale_color_manual(name="Bat Family", values=fam.col) +
#   scale_x_continuous(breaks = c(1,2,3), labels = c("10", "100", "1000"))
# 
# print(pS4_poster)

ggsave(file = paste0(homewd,"/figures/FigS4_bat_only_poster.png"),
       units="mm",  
       width=100, 
       height=80, 
       scale=3, 
       dpi=300)

out.plot2 <- cowplot::plot_grid(pA, pB,pS4, nrow=1, ncol=3, labels=c("(a)", "(b)", "(c)"), rel_widths = c(1,1, 1.2))

print(out.plot2)

ggsave(file = paste0(homewd,"/figures/Fig1_ThreePanel.jpeg"),
       units="mm",  
       width=180, 
       height=60, 
       scale=3, 
       dpi=300)

