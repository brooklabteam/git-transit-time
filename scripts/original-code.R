rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(lmerTest)
library(sjPlot)
library(ggforce)


#set your home directory
#homewd= "/Users/carabrook/Developer/git-transit-time"
homewd= "/Users/katherinemcferrin/Developer/git-transit-time"

#load the bat data:
dat <- read.csv(file = paste0(homewd, "/data/final-GIT-transit-database-april-2021.csv"), header = T, stringsAsFactors = F )
head(dat)

length(unique(dat$Retention.Citation)) #97 unique papers
length(unique(dat$Genus.species)) #109 unique species

#select only the columns we will use for modleing and plotting

dat <- dplyr::select(dat, can_fly, Class, Order, Family, Genus.species, Common.Name, Typical.Diet, Mass..mean.median.from.study., Feeding.Trial.Food.Cat, N_measured, N_trials, Minimum..min., Median..min., Mean..min.,Maximum..min.)
head(dat)
names(dat) <- c("fly", "class", "order", "family", "genus.species", "common.name", "typical.diet", "mass", "food.cat", "N_individuals", "N_trials", "min", "median", "mean", "max")
head(dat)



#choose one mean/median to report
dat.plot <- dat
dat.plot$transit <- dat.plot$mean
dat.plot$transit[is.na(dat.plot$transit)]<- dat.plot$median[is.na(dat.plot$transit)]

dat.plot <- dplyr::select(dat.plot, -(min), -(median), -(mean), -(max))
head(dat.plot)

#summarize by species
dat.plot$N_individuals[dat.plot$N_individuals=="not reported"] <- 1
dat.plot$N_individuals[is.na(dat.plot$N_individuals)] <- 1
dat.plot$N_individuals = as.numeric(dat.plot$N_individuals)
dat.plot$total_transit = dat.plot$transit*dat.plot$N_individuals
head(dat.plot)
unique(dat.plot$food.cat)


#now get one entry for each species/food category combination
dat.split <- dlply(dat.plot, .(fly,class,order,family,genus.species, common.name, typical.diet,food.cat))
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
#Peramelemorphia, Perissodactyla, Sirenia
dat.sum.tot = subset(dat.sum.tot, class!="Chondrichthyes" & class !="Amphibia")
#and write over some of the otheres
dat.sum.tot$re_class[dat.sum.tot$re_class=="Reptilia"] <- "Reptiles"
#dat.sum.tot$re_class[dat.sum.tot$re_class=="Lagomorpha"] <- "Rodents and Lagomorphs"
dat.sum.tot$re_class[dat.sum.tot$re_class=="Rodentia"] <- "Rodents"
dat.sum.tot$re_class[dat.sum.tot$re_class=="Chiroptera"] <- "Bats"
dat.sum.tot$re_class[dat.sum.tot$re_class=="Artiodactyla"] <- "Even-Toed Ungulates"
dat.sum.tot$re_class[dat.sum.tot$re_class=="Carnivora"] <- "Carnivores"
#dat.sum.tot$re_class[dat.sum.tot$re_class=="Perissodactyla" | dat.sum.tot$re_class=="Artiodactyla"] <- "Ungulates"
dat.sum.tot$re_class[dat.sum.tot$re_class=="Perissodactyla" | dat.sum.tot$re_class=="Lagomorpha" | dat.sum.tot$re_class=="Eulipotyphla" | dat.sum.tot$re_class=="Peramelemorphia" | dat.sum.tot$re_class=="Sirenia" |  dat.sum.tot$re_class=="Dermoptera" ] <- "Other Mammals"

dat.sum.tot = subset(dat.sum.tot, re_class!="Non-Flying Birds" & re_class!= "Other Mammals")
#dat.sum.tot = subset(dat.sum.tot,  re_class!= "Other Mammals")
unique(dat.sum.tot$re_class)

#convert transit time to hours
dat.sum.tot$transit_hrs = dat.sum.tot$transit/60
dat.sum.tot$mass_kg = as.numeric(dat.sum.tot$mass)/1000

#check how transit time varies by mass
#generally, we see that transit is longer at higher mass
#all taxa agree with this trend, though reptiles and bats are a little longer in time
#than we might otherwise predict -- this is likely due to an effect of temperature that is not
#reported here.

#and check how transit varies by food category
#taxa appear to have a stronger effect than food type
#reclass.colz =  c("Rodents and Lagomorphs" = "black", "Bats"= "black", "Flying Birds"= "black", "Reptiles"= "black",  "Ungulates"= "black", "Carnivores"= "black", "Primates"= "black")
#reclass.colz =  c("Rodents" = "black", "Bats"= "black", "Flying Birds"= "black", "Reptiles"= "black",  "Even-Toed Ungulates"= "black", "Carnivores"= "black", "Primates"= "black")
#and replot the mass by transit with categories
dat.sum.tot$re_class <- factor(dat.sum.tot$re_class, levels = c( "Flying Birds","Bats", "Rodents", "Carnivores", "Primates", "Even-Toed Ungulates",   "Reptiles"))


#now some stats
dat.sum.tot$re_class = factor(dat.sum.tot$re_class, levels= c( "Rodents",  "Flying Birds","Bats","Carnivores", "Primates", "Even-Toed Ungulates",  "Reptiles"))
dat.sum.tot$food.cat = factor(dat.sum.tot$food.cat, levels= c("label-solution", "protein", "fiber/foliage", "fruit/nectar/pollen"))

m1 <- lmer(log10(transit_hrs)~ re_class + (1|food.cat), data=dat.sum.tot, REML = F)
summary(m1)
rand(m1) 


#for plotting later
m1b <- lme4::lmer(log10(transit_hrs)~ re_class + (1|food.cat), data=dat.sum.tot, REML = F)
summary(m1b)


#and here plot the raw partial effects of taxonomic grouping for Figure S2:


#plot fixed effects
pS2 <- plot_model(m1b, type="est", vline.color = "black",
                  axis.labels = rev(c("Flying Birds", "Bats", "Carnivores", "Primates", "Even-Toed Ungulates", "Reptiles")),
                  title = "Effects of taxon on GIT transit") + theme_bw() + 
  theme(panel.grid = element_blank(), axis.text = element_text(size=14), axis.title.x = element_blank()) 

print(pS2)

ggsave(file = paste0(homewd,"/figures/April2021_database/FigS2_Taxon_effects.png"),
       units="mm",  
       width=50, 
       height=40, 
       scale=3, 
       dpi=300)



m2 <- lmer(log10(transit_hrs)~ log10(mass_kg):re_class + (1|re_class), data=dat.sum.tot, REML = F)
summary(m2)
rand(m2) 


dat.sum.tot$predicted_transit[!is.na(dat.sum.tot$transit_hrs) & !is.na(dat.sum.tot$mass_kg) & !is.na(dat.sum.tot$re_class)]   <- 10^(predict(m2))



#and reshape order
dat.sum.tot$re_class <- as.character(dat.sum.tot$re_class)
dat.sum.tot$re_class[dat.sum.tot$re_class=="Flying Birds"] <- "Flying\nBirds"
dat.sum.tot$re_class[dat.sum.tot$re_class=="Even-Toed Ungulates"] <- "Even-Toed\nUngulates"
dat.sum.tot$re_class <- factor(dat.sum.tot$re_class, levels = c( "Flying\nBirds","Bats", "Rodents", "Carnivores", "Even-Toed\nUngulates", "Primates",  "Reptiles"))
dat.sum.tot$food.cat = factor(dat.sum.tot$food.cat, levels= c("fiber/foliage", "fruit/nectar/pollen","protein","label-solution"))


scales::show_col(scales::hue_pal()(7))
colz =  c( "Flying\nBirds"= "#F8766D", "Bats" = "#C49A00", "Rodents"= "navy","Carnivores" = "#00C094",  "Primates"= "#00B6EB", "Even-Toed\nUngulates"= "#A58AFF",   "Reptiles" = "#FB61D7")
shapez = c("fiber/foliage" = 21, "fruit/nectar/pollen" = 22, "label-solution" = 23, "protein" = 24)


#and the version that is 2panel
label.dat <- cbind.data.frame(re_class=c("Flying\nBirds","Bats", "Rodents", "Carnivores", "Primates", "Even-Toed\nUngulates",   "Reptiles"), label = c("***", "**", "", "*", "***", "***", "***"))
label.dat$re_class <- factor(label.dat$re_class, levels=c( "Flying\nBirds","Bats", "Rodents", "Carnivores", "Primates", "Even-Toed\nUngulates", "Reptiles"))

#and all together
pA <- ggplot(data=dat.sum.tot)+ geom_boxplot(aes(x=re_class, y=log10(transit_hrs), fill=re_class), show.legend = F) + theme_bw() + 
  theme(axis.title.x = element_blank(), panel.grid = element_blank(), legend.position = c(.85,.15), legend.background = element_rect(color="black"),
        axis.text.y = element_text(size=13), axis.text.x = element_text(size=13), axis.title.y = element_text(size=16),
        plot.margin = unit(c(.2,.2,.8,.2), "cm")) + scale_fill_manual(values=colz) + 
  geom_point(aes(x=re_class, y=log10(transit_hrs),shape=food.cat), position = position_jitterdodge(jitter.width = 0), size=3) +  
  scale_shape_manual(values=shapez, name="food type") +  ylab("gut transit time, hrs") +
  geom_hline(aes(yintercept=log10(quantile(subset(dat.sum.tot, re_class=="Rodents")$transit_hrs)["50%"])), linetype=2) +
  scale_y_continuous(breaks=c(0,1,2), labels = c("1", "10", "100")) +
  coord_cartesian(ylim=c(-.8,2.95)) + geom_label(data=label.dat, aes(x=re_class, y=2.9, label=label), label.size = NA,size=5)
print(pA)



#now plot with mass
#bats are the only taxon for which the mass slope from model 2 is negative (and significantly so)

head(dat.sum.tot)

pB <- ggplot(data=dat.sum.tot, aes(x=log10(mass_kg), y=log10(transit_hrs))) + 
  geom_mark_ellipse(expand=0,radius=0,aes(fill=re_class, color=re_class), size=.1)+
  geom_line(aes(x=log10(mass_kg), y=log10(predicted_transit), color=re_class), alpha=.6, size=.7) +
  #geom_ribbon(aes(x=log10(mass_kg), ymin=predicted_transit_lci, ymax=predicted_transit_uci, fill=re_class), alpha=.1) +
  geom_point(aes(x=log10(mass_kg), y=log10(transit_hrs), fill= re_class, shape=food.cat), size=4, show.legend = F) + 
  geom_point(aes(x=log10(mass_kg), y=log10(transit_hrs), fill= re_class, color= re_class, shape=food.cat), size=3, show.legend = F) + 
  ylab("gut transit time (hrs)") + xlab("mass (kg)") + scale_shape_manual(name = "food type", values=shapez) + 
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




out.plot <- cowplot::plot_grid(pA, pB, nrow=1, ncol=2, labels=c("(a)", "(b)"))

print(out.plot)

ggsave(file = paste0(homewd,"/figures/April2021_database/Fig4_TwoPanel.jpeg"),
       units="mm",  
       width=130, 
       height=60, 
       scale=3, 
       dpi=300)


#and the supplementary figures

#plot random effects
pS3 <- plot_model(m1b, type="re", vline.color = "black") + theme_bw() + 
  theme(panel.grid = element_blank(), axis.text = element_text(size=14)) 

print(pS3)

ggsave(file = paste0(homewd,"/figures/April2021_database/FigS3_Random_Effects.png"),
       units="mm",  
       width=50, 
       height=40, 
       scale=3, 
       dpi=300)

fam.col = c("Molossidae" = "firebrick", "Mormoopidae" = "navy", "Phyllostomidae" = "darkmagenta", "Pteropodidae" = "forestgreen", "Vespertilionidae" = "lightcoral")

#and look at just bats - reverses the direction
pS4 <- ggplot(data=subset(dat.sum.tot, re_class=="Bats")) + 
  geom_point(aes(x=log10(mass), y=log10(transit), fill= family, shape=food.cat), size=4, show.legend = F) + 
  geom_point(aes(x=log10(mass), y=log10(transit), fill= family, color= family, shape=food.cat), size=3) + 
  ylab("gut transit time (min)") + xlab("mass (g)") + coord_cartesian(xlim=c(.5,3), ylim=c(.5,3))+
  theme_bw() + theme(panel.grid = element_blank(), #legend.position = c(.87,.83), 
                     #legend.background = element_rect(color="black"),
                     axis.text = element_text(size=14), axis.title = element_text(size=18)) +
  scale_shape_manual(values=shapez, name= "food type") + 
  scale_y_continuous(breaks = c(1,2,3), labels= c("10", "100", "1000")) + 
  scale_fill_manual(name="bat family", values=fam.col) + scale_color_manual(name="bat family", values=fam.col) +
  scale_x_continuous(breaks = c(1,2,3), labels = c("10", "100", "1000")) #+
#guides(shape=FALSE)#, nrow=1, colour = guide_legend(nrow = 1))

print(pS4)
ggsave(file = paste0(homewd,"/figures/April2021_database/FigS3_bat_only.png"),
       units="mm",  
       width=80, 
       height=60, 
       scale=3, 
       dpi=300)