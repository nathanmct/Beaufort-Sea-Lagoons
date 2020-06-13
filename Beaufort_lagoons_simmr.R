###################
##HARRIS, c.M, N.D. MCTIGUE, J.W. MCLELLAND, & K.H. DUNTON (2018). DO HIGH ARCTIC COASTAL FOOD WEBS RELY ON A TERRESTRIAL CARBON SUBSIDY?
###################

#Some of the below script has been borrowed from other sources: 
#1) Passoti et al. (2015) Benthic trophic interactions in an Antarctic shallow water ecosystem affected by recent glacial retreat. PLoS ONE 10(11):e0141742. doi:10.1371/journal.pone.0141742
#2) The vignette for the package simmr available at https://cran.r-project.org/web/packages/simmr/vignettes/simmr.html
#3) From Andrew Parnell's GitHub for simmr (https://github.com/andrewcparnell/simmr/tree/master/R). The functionality is the same, but I have altered some of the ggplot code to change the output visualization.
#4) McTigue & Dunton (2017) Trophodynamics of the Hanna Shoal Ecosystem (Chukchi Sea, Alaska): Connecting multiple end-members to a rich food web. Deep-Sea Research II: 144:175-189. http://dx.doi.org/10.1016/j.dsr2.2017.08.010 (see https://github.com/nathanmct/Chukchi-Sea-trophodynamics)

#Much of it has been modified for our own purposes, but we acknowledge and thank these sources for their open-source code.
#The script does not produce every figure and table in the order they appear in the manuscript, but rather follows our logic for data analysis. However, each product's label corresponds with the manuscript.

#JAGS must be installed on your computer to run some of these scripts. Install it here: https://sourceforge.net/projects/mcmc-jags/files/
#It does not need to be open when running R

#Be sure to put all helper functions in working directory

###################
#This script generates Fig. 2-6 from Harris et al. (2018) along with much of the analysis within.
#Figure generation is not necessarily in numerical order, but is provided in logical flow with the below explanations
#The code is divided into three main sections: 
#1) Data wrangling, tidying, and exploration (Fig. 2 & 3)
#2) Running simmr and producing boxplots of posterior distributions (Fig. 5 & 6)
#3) Visualizing stable isotope biplots (Fig 4)
#
#

##########################################
#Data wrangling, tidying, and exploration#
##########################################

rm(list = ls())
setwd("P://Your//Directory//Here")

library(rjags)
library(simmr)
library(siar)
library(grid)
library(plyr)
library(tidyverse)
library(SIBER)
library(xlsx)
library(egg)
library(R2jags)
library(pander)
library(reshape2)
library(mcmcplots)

#we used this color palette, modified from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
cbpalette <- (c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#F781BF", "#999999", "#32CD32", "#00008B", "#A65628", "#FFD900", "#8B6969"))

###############
#Facet palette#
###############
facet_palette <- c(AN = "#E41A1C", #red
                   BE = "#A65628", #brown
                   DE = "#00008B", #navy
                   DP = "#FFD900", #yellow
                   HU = "#4DAF4A", #green
                   JA = "#FF7F00", #orange
                   KA = "#984EA3", #purple
                   NU = "#377EB8") #light blue


palette_5color <- (c("Ss/De" = "#E41A1C",
                    "Ep/Om" = "#377EB8",
                    "Fish" = "#4DAF4A",
                    "Mam/Carn" = "#984EA3",
                    "Su/FF" = "#FF7F00"))

#import data
#must be sorted alphabetically by genus!
#iso <- read.csv("BeaufortSeaLagoons_isotope_data.csv", header=T)
iso <- read.csv("BeaufortSeaLagoons_isotope_data.csv", header = TRUE)
iso <- iso[order(iso$Genus),]

#Count the number of each Genus and add a new column called 'idcount' with the count. Used in filter() in next step. 
Genus.count <- iso %>% dplyr::count(Genus)
iso <- iso %>% dplyr::mutate(idcount = (rep(Genus.count$n, Genus.count$n))) 
iso$Year <- as.factor(iso$Year)

#Filter dataset so n>4 for genera to be run in simmr
simmriso <- dplyr::filter(iso, idcount>4)
simmriso <- droplevels(simmriso)

#Add corrections (or trophic enrichment factors) for simmr
corrections.HTL <- read.csv("corrections.higherTL.csv", header=T)
corrections <- read.csv("corrections.csv")
corrections.mamm <- read.csv("corrections.mammTL.csv")

#Add food web end-members for simmr
sources <- read.csv("source.averages.csv",header=T)
sources.mamm <- read.csv("source.averages.mammals.csv")

#filter 'consumers' to separate by trophic guild for some graphing applications in last section
susp <- dplyr::filter(simmriso, Feeding.Type=="Su/FF")
depo <- dplyr::filter(simmriso, Feeding.Type=="Ss/De")
omni <- dplyr::filter(simmriso, Feeding.Type=="Ep/Om")
fish <- dplyr::filter(simmriso, Feeding.Type=="Fish")
mamm <- dplyr::filter(iso, Feeding.Type == "Mam/Carn")
beluga <- dplyr::filter(iso, Genus == "Delphinapterus")

####################################################################
##Add 'Codes' to each Genus. Codes are required by simmr. 
##Each Genus requires a different code, as a sequential integer. 
##The same genera have the same code.
####################################################################

Genus.count.susp <- susp %>% dplyr::count(Genus) #counts the number of each genus. Will be used in rep() as the number of times to repeat the code.
susp <- susp %>% droplevels() %>% #use piping to overwrite the original object. Drop levels here from previous filter()
  dplyr::mutate(Code = (rep(1:length(unique(Genus)), Genus.count.susp$n))) #use mutate() to add a column called 'Code'

Genus.count.depo <- depo %>% dplyr::count(Genus)
depo <- depo %>% droplevels() %>%
  dplyr::mutate(Code = (rep(1:length(unique(Genus)), Genus.count.depo$n)))

Genus.count.omni <- omni %>% dplyr::count(Genus)
omni <- omni %>% droplevels() %>%
  dplyr::mutate(Code = (rep(1:length(unique(Genus)), Genus.count.omni$n)))

Genus.count.fish <- fish %>% dplyr::count(Genus)
fish <- fish %>% droplevels() %>%
  dplyr::mutate(Code = (rep(1:length(unique(Genus)), Genus.count.fish$n)))

#Mammalian/carnivore group different. Need to "cheat" and count all genera as the same for simmr input
mamm <- mamm %>% droplevels() %>%
  dplyr::mutate(Code = (rep(1, length(mamm$d15N))))

#Make summary table
iso.table <- dplyr::group_by(iso, Genus) %>%
  summarise(avg.d15N = mean(d15N),
            sd.d15N = sd(d15N),
            avg.d13C = mean(d13C),
            sd.13C = sd(d13C),
            avg.cn = mean(molar.cn),
            sd.cn = sd(molar.cn),
            n = max(idcount)
  )

write.xlsx(iso.table, file = "Stable Isotope Summary Table.xlsx")

#Jitterplots
d13C.jitter <- ggplot(data = iso, aes(x = Feeding.Type, y = d13C, color = Feeding.Type))+
  geom_jitter(size = 2, alpha = 0.7, width = 0.25)+
  scale_color_manual(values = palette_5color)+
  ylab(expression(paste(delta^{13},'C (\u2030)')))+
  xlab("")+
  theme_bw()+
  theme(text = element_text(size = 16),
        axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18),
        legend.position = "none",
        panel.grid.major = element_blank(),#removes panel grids
        panel.grid.minor = element_blank())+
  coord_flip()+
  scale_y_continuous(breaks = seq(-27, -15, 3))

d13C.jitter

#ggsave(filename="Figure_2a.png", plot = d13C.jitter, device = "png", width = 8, height = 5, units = "in", dpi=330)

d15N.jitter <- ggplot(data = iso, aes(x = Feeding.Type, y = d15N, color = Feeding.Type))+
  geom_jitter(size = 2, alpha = 0.7, width = 0.25)+
  scale_color_manual(values = palette_5color)+
  ylab(expression(paste(delta^{15},'N (\u2030)')))+
  xlab("")+
  theme_bw()+
  theme(text = element_text(size = 16),
        axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 18),
        legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank())+
  coord_flip()+
  scale_y_continuous(breaks = seq(6, 22, 3))

d15N.jitter

#g1 <- ggplotGrob(d13C.jitter)
#g2 <- ggplotGrob(d15N.jitter)

#source("gtable_frame.r")
#fg1 <- gtable_frame(g1)
#fg2 <- gtable_frame(g2)

#grid.newpage()
#combined <- cbind(fg1, fg2)
#grid.draw(combined)

#library(egg)
#library(grid)
grid.newpage()
png("P:/My Documents/Pubs/Harris et al/simmr/jitterplot.png", height = 6, width = 8, units="in", res = 600)
grid.draw(ggarrange(d13C.jitter, d15N.jitter, ncol = 2))
dev.off()

#Supplemental figure for data exploration
ggplot(data = iso, aes(x = d13C, y = d15N, color = Feeding.Type))+
  geom_point(size = 2.5, alpha = 0.75)+
  scale_color_manual(values = cbpalette,
                     name = "Feeding Type")+
  ylab(expression(paste(delta^{15},'N (\u2030)')))+
  xlab(expression(paste(delta^{13},'C (\u2030)')))+
  theme_bw()+
  theme(axis.text = element_text(size=14, color = "black"),
        legend.text = element_text(size=16),
        legend.title = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(color="black"),
        panel.grid.major = element_blank(),#removes panel grids
        panel.grid.minor = element_blank(),
        axis.title = element_text(size=20),
        strip.text = element_text(size=18)
  )+
  coord_fixed(ratio = 1)+
  scale_x_continuous(breaks = seq(-27, -17, 5))+
  facet_wrap(~Location, ncol = 4)

ggsave(filename = "Station_facet_by_FeedingType.png", plot = last_plot(), device = "png", width = 7, height = 5, units = "in", dpi=300)

#Supplemental Figure 2 for data exploration
ggplot(data = iso, aes(x = d13C, y = d15N, color = Location))+
  geom_point(size = 3, alpha = 0.8)+
  scale_color_manual(values = cbpalette,
                     name = "Feeding Type")+
  ylab(expression(paste(delta^{15},'N (\u2030)')))+
  xlab(expression(paste(delta^{13},'C (\u2030)')))+
  theme_bw()+
  theme(text = element_text(size = 12),
        legend.text = element_text(size=20),
        legend.title = element_blank(),
        legend.position = "right",
        legend.direction = "vertical",
        legend.background = element_rect(color="black"),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 6)
  )+
  coord_fixed(ratio = 1)+
  scale_x_continuous(breaks = seq(-27, -17, 5))+
  facet_wrap(~Genus, ncol = 7)

ggsave(filename = "Genus_facet_by_Location.png", plot = last_plot(), device = "png", width = 10, height = 8, units = "in", dpi=600)





####################################################
##  SIBER used in earlier versions of manuscript  ##
##  These are not present in final version        ##
####################################################

#These are from Pasotti et al. (2015)
source("SEA_fitting.r")
source("BUGS_models.r")
source("comp_SEAb.r")
source("extract_SEAc.r")

#Parameters and priors for SIBER model
parms <- list()
parms$n.iter <- 2 * 10^5
parms$n.burning <- 1 * 10^4
parms$n.thin <- 10
parms$n.chains <- 2

priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

#Run SIBER for all species in dataset, subset by trophic guild
SEA.Jags.Trophic <- with(iso, siber.ellipses.test(d13C, d15N, Feeding.Type, method = "IWJAGS", parms, priors))

#Find standard ellipse paths for Trophic guilds
ell.Trophic <-dlply(iso,.(Feeding.Type),function(x) standard.ellipse(x$d13C,x$d15N,steps = 1))

e2.Trophic <- ldply(ell.Trophic, extract_SEAc) #this is object read by ggplot for SEA paths
names(e2.Trophic)[2:4] <-c("SEAc", "xSEAc", "ySEAc")

##########
#Figure 3#
##########

#plot ALL species, but only Trophic Guild SEAs
ggplot(data = iso, aes(x = d13C, y = d15N))+
  geom_point(size = 1.8, alpha = 0.6, aes(color = Feeding.Type))+
  geom_path(data = e2.Trophic, aes(x = xSEAc, y = ySEAc, color = Feeding.Type), linetype = 1, size = 2, alpha = 1)+
  scale_color_manual(values = palette_5color,
                     name = "Trophic Guild",
                     guide = guide_legend(override.aes = list(shape = NA)))+
  ylab(expression(paste(delta^{15},'N (\u2030)')))+
  xlab(expression(paste(delta^{13},'C (\u2030)')))+
  theme_bw()+
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 12, color="black"),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  coord_fixed(ratio = 1)+
  #scale_y_continuous(limits = c(5,20))+
  scale_x_continuous(breaks = seq(-28, -15, 4))

ggsave(filename = "all_SIBER_biplot.png", plot = last_plot(), device = "png", width = 7, height = 5, units = "in", dpi=600)


#############
#SIBER Table#
#############

#compare probability that g1 > g2 in isospace
comp_SEAb(SEA.Jags.Trophic$SEA.B, iso$Feeding.Type)

#make table with statistics for trophic guilds
met2 <- data.frame(SEA.Jags.Trophic$SEA)
names(met2) <- levels(iso$Feeding.Type)
met2 <- melt(met2) #throws warning about no id variables; proceed anyway
names(met2) <- c("Feeding.Type", "value")
hh <- function(x){
  h <- hdr(x$value,h = bw.nrd0(x$value))
  data.frame(mean = mean(x$value), sd = sd(x$value), mode = h$mode, lo95 = h$hdr[2,1],hi95 = h$hdr[2,2])
}
tSEAb <- ddply(met2,.(Feeding.Type), hh)

tSEAb <- cbind(tSEAb, ddply(e2.Trophic,.(Feeding.Type), summarize, unique(SEAc))[,2])
names(tSEAb)[7] <- "SEAc"

knitr::kable(tSEAb)

write.xlsx(tSEAb, file = "SIBER_stats_trophc_guild.xlsx")

#Make table of overlapping SEAc isospace
s <- as.character(levels(iso$Feeding.Type))

ovrlpSEAc <- function(x) {
  ov1 <- iso[iso$Feeding.Type==s[x[1]],]
  ov2 <- iso[iso$Feeding.Type==s[x[2]],]
  data.frame(g1=s[x[1]], g2=s[x[2]],overlap(ov1$d13C,ov1$d15N,ov2$d13C,ov2$d15N,steps=1))
}

comp1 <- do.call(rbind, apply(combn(seq_along(s), 2), 2, ovrlpSEAc))

knitr::kable(comp1)

write.xlsx(comp1, file="TrophicGuild_SEA_overlap.xlsx")

######################################################################################
#Need similar table with statistics for all species that will be analyzed (n > 9)
#Re-do above procedure, but now: 
#(1) subset n > 9, and  (2) do not divide by trophic descriptor
######################################################################################

iso.SIBER.genus<-dplyr::filter(iso, idcount > 9) %>% droplevels()

SEA.Jags.iso <- with(iso.SIBER.genus, siber.ellipses.test(d13C, d15N, Genus, method = "IWJAGS", parms, priors))

#Find standard ellipse paths for all trophic guilds and save as Excel sheet
ell.iso <-dlply(iso.SIBER.genus,.(Genus), function(x) standard.ellipse(x$d13C, x$d15N, steps=1))

extract_SEAc <- function(e) {
  data.frame(e$SEAc,e$xSEAc, e$ySEAc)
}

e2.iso <- ldply(ell.iso, extract_SEAc)
names(e2.iso)[2:4] <-c("SEAc", "xSEAc", "ySEAc")

met2 <- data.frame(SEA.Jags.iso$SEA)
names(met2) <- levels(iso.SIBER.genus$Genus)
met2 <- melt(met2)
names(met2) <- c("Genus", "value")
hh <- function(x){
  h <- hdr(x$value,h = bw.nrd0(x$value))
  data.frame(mean = mean(x$value), sd=sd(x$value), mode = h$mode,lo95 = h$hdr[2,1], hi95 = h$hdr[2,2])
}
tSEAb <- ddply(met2,.(Genus), hh)

tSEAb <- cbind(tSEAb, ddply(e2.iso,.(Genus), summarize, unique(SEAc))[,2])
names(tSEAb)[7]<-"SEAc"

knitr::kable(tSEAb)

write.xlsx(tSEAb, file = "Genus_SEA_Table.xlsx")

#Make table of overlapping SEAc isospace
s <- as.character(levels(iso.SIBER.genus$Genus))

ovrlpSEAc <- function(x) {
  ov1 <- iso.SIBER.genus[iso.SIBER.genus$Genus==s[x[1]],]
  ov2 <- iso.SIBER.genus[iso.SIBER.genus$Genus==s[x[2]],]
  data.frame(g1=s[x[1]], g2=s[x[2]],overlap(ov1$d13C,ov1$d15N,ov2$d13C,ov2$d15N,steps=1))
}

comp1 <- do.call(rbind, apply(combn(seq_along(s), 2), 2, ovrlpSEAc))

knitr::kable(comp1)

write.xlsx(comp1, file="Genera_SEA_overlap.xlsx")

#Biplot of genera that underwent SIBER
ggplot(data = iso.SIBER.genus, aes(x = d13C, y = d15N))+
  geom_path(data = e2.iso, aes(x = xSEAc, y = ySEAc, color = Genus), linetype = 1, size = 2)+
  #geom_point(size = 2, aes(color = Genus))+
  scale_color_manual(values = cbpalette,
                     name = "Genera")+
  ylab(expression(paste(delta^{15},'N (\u2030)')))+
  xlab(expression(paste(delta^{13},'C (\u2030)')))+ #turned off for gridding
  #xlab(expression(""))+
  theme_bw()+
  theme(text = element_text(size=16),
        axis.text = element_text(size=16, color="black"),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16),
        panel.grid.major = element_blank(),#removes panel grids
        panel.grid.minor = element_blank())+
  coord_fixed(ratio = 1)+
  scale_y_continuous(limits = c(5, 23), breaks = seq(5, 20, 5))+
  scale_x_continuous(limits = c(-28, -15), breaks = seq(-28, -15, 4))
  #annotate("text", x = -25.5, y = 17.5, label = "a", size = 7)

ggsave(filename = "Genus_SEA_biplot.png", plot = last_plot(), device = "png", width = 7, height = 7, units = "in", dpi=600)





####################################################################
####################################################################
###simmr (stable isotope mixing model in R)                        #
###See Parnell et al. (2013) in Envirometrics                      #
#https://cran.r-project.org/web/packages/simmr/vignettes/simmr.html#
####################################################################
####################################################################

#preen data for use with simmr
#mix=mixture, i.e., the consumers
#s=sources, i.e., the endmembers
#c=corrections, i.e., the TEFs for each endmember
#see the simmr vignette for more detail
mix.susp <- as.matrix.data.frame(susp[,c(6:5)], ncol=2, nrow=length(susp$d13C)) #in these df's d13C is the 6th column, and d15N is the 5th
mix.depo <- as.matrix.data.frame(depo[,c(6:5)], ncol=2, nrow=length(depo$d13C))
mix.omni <- as.matrix.data.frame(omni[,c(6:5)], ncol=2, nrow=length(omni$d13C))
mix.fish <- as.matrix.data.frame(fish[,c(6:5)], ncol=2, nrow=length(fish$d13C))
mix.mamm <- as.matrix.data.frame(mamm[,c(6:5)], ncol=2, nrow=length(mamm$d13C))
mix.beluga <- as.matrix.data.frame(beluga[,c(6:5)], ncol=2, nrow=length(beluga$d13C))

s_names <- levels(sources$Source)
s_means <- as.matrix.data.frame(sources[,c(2,4)], ncol=2, nrow=length(sources$Source))
s_sds <- as.matrix.data.frame(sources[,c(3,5)], ncol=2, nrow=length(sources$Source))

s_mamm_names <- levels(sources.mamm$Source)
s_mamm_means <- as.matrix.data.frame(sources.mamm[,c(2,4)], ncol=2, nrow=length(sources.mamm$Source))
s_mamm_sds <- as.matrix.data.frame(sources.mamm[,c(3,5)], ncol=2, nrow=length(sources.mamm$Source))

c_means <- as.matrix.data.frame(corrections[,c(2,4)], ncol=2, nrow=length(sources$Source))
c_sds <- as.matrix.data.frame(corrections[,c(3,5)], ncol=2, nrow=length(sources$Source))
c_means_HTL <- as.matrix.data.frame(corrections.HTL[,c(2,4)], ncol=2, nrow=length(sources$Source))
c_sds_HTL <- as.matrix.data.frame(corrections.HTL[,c(3,5)], ncol=2, nrow=length(sources$Source))
c_means_mamm <- as.matrix.data.frame(corrections.mamm[,c(2,4)], ncol=2, nrow=length(sources$Source))
c_sds_mamm <- as.matrix.data.frame(corrections.mamm[,c(3,5)], ncol=2, nrow=length(sources$Source))

grp.susp <- as.integer(as.matrix(susp$Code))
grp.depo <- as.integer(as.matrix(depo$Code))
grp.omni <- as.integer(as.matrix(omni$Code))
grp.fish <- as.integer(as.matrix(fish$Code))
grp.mamm <- as.integer(as.matrix(mamm$Code))
grp.beluga <- as.integer(as.matrix(rep(1, length(beluga$d13C))))

#load data so that simmr can use it
#create a different object for each trophic guild

simmr_groups_susp <- simmr_load(mixtures=mix.susp,
                               source_names=s_names,
                               source_means=s_means,
                               source_sds=s_sds,
                               correction_means=c_means,
                               correction_sds=c_sds,
                               group=grp.susp)

simmr_groups_depo <- simmr_load(mixtures=mix.depo,
                               source_names=s_names,
                               source_means=s_means,
                               source_sds=s_sds,
                               correction_means=c_means,
                               correction_sds=c_sds,
                               group=grp.depo)

simmr_groups_omni <- simmr_load(mixtures=mix.omni,
                               source_names=s_names,
                               source_means=s_means,
                               source_sds=s_sds,
                               correction_means=c_means,
                               correction_sds=c_sds,
                               group=grp.omni)

simmr_groups_fish <- simmr_load(mixtures=mix.fish,
                               source_names=s_names,
                               source_means=s_means,
                               source_sds=s_sds,
                               correction_means=c_means_HTL,
                               correction_sds=c_sds_HTL,
                               group=grp.fish)

simmr_groups_mamm <- simmr_load(mixtures=mix.mamm,
                               source_names=s_mamm_names,
                               source_means=s_mamm_means,
                               source_sds=s_mamm_sds,
                               correction_means=c_means_mamm,
                               correction_sds=c_sds_mamm,
                               group=grp.mamm)

simmr_groups_beluga <- simmr_load(mixtures=mix.beluga,
                                source_names=s_mamm_names,
                                source_means=s_mamm_means,
                                source_sds=s_mamm_sds,
                                correction_means=c_means_mamm,
                                correction_sds=c_sds_mamm,
                                group=grp.beluga)

#run MCMC for each group
#need JAGS loaded on computer
#should take <1 minute for each line
susp_out = simmr_mcmc(simmr_groups_susp)
depo_out = simmr_mcmc(simmr_groups_depo)
omni_out = simmr_mcmc(simmr_groups_omni)
fish_out = simmr_mcmc(simmr_groups_fish)
mamm_out = simmr_mcmc(simmr_groups_mamm)
beluga_out = simmr_mcmc(simmr_groups_beluga)



#################################################################################
##  Visually assess convergence of chains and make note of Gelman diagnostics  ##
#################################################################################

summary(susp_out, type = "diagnostics")
summary(depo_out, type = "diagnostics")
summary(omni_out, type = "diagnostics")
summary(fish_out, type = "diagnostics")
summary(mamm_out, type = "diagnostics")

#Plot OM source correlations
plot(susp_out, type='matrix', group=1:max(grp.susp))
plot(depo_out, type='matrix', group=1:max(grp.depo))
plot(omni_out, type='matrix', group=1:max(grp.omni))
plot(fish_out, type='matrix', group=1:max(grp.fish))
plot(mamm_out, type='matrix', group=1:max(grp.mamm))
                
#can also use type = “isospace”, “histogram”, “density”, “boxplot” for different visualizations

#obtain 95% credible intervals...
#save each as its own object so we can pull out values later to make a new df
#the new df will be used for graphing
#can save each obejct as xlsx and read in later so you don't have to re-run simmr
susp.summ <- summary(susp_out, type=c("quantiles"), group=c(1:max(susp$Code)))
#write.xlsx(susp.summ, "simmr_susp_summary.xlsx")

depo.summ <- summary(depo_out, type=c("quantiles"), group=c(1:max(depo$Code)))
#write.xlsx(depo.summ, "simmr_depo_summary.xlsx")

omni.summ <- summary(omni_out, type=c("quantiles"), group=c(1:max(omni$Code)))
#write.xlsx(omni.summ, "simmr_omni_summary.xlsx")

fish.summ <- summary(fish_out, type = c("quantiles"), group = c(1:max(fish$Code)))
#write.xlsx(fish.summ, "simmr_fish_summary.xlsx")

mamm.summ <- summary(mamm_out, type = c("quantiles"), group = 1)
#write.xlsx(mamm.summ, "simmr_mamm_summary.xlsx")

beluga.summ <- summary(beluga_out, type = c("quantiles"), group = 1)




############################################################################
##  extract data from summary.simmr_output lists for customized plotting  ##
############################################################################

###########
##  MPB  ##
###########

#First, make empty df to put values from simmr posteriors list (objects just made in previous step) using for loop

#Deposit Feeders
depo.MPB <- data.frame(matrix(NA, nrow = max(depo$Code), ncol = 5))
bpnames <-  c("min", "low", "mid", "top", "max")
colnames(depo.MPB) <- bpnames

#location in list where quantiles are stored
#selects quantiles for each genus and puts in new df
for (i in 1:max(depo$Code)) {
  depo.MPB[i,] <- (depo.summ[[2]][[i]][2,])
}

depo.names <- data.frame(matrix(levels(depo$Genus), ncol = 1))
colnames(depo.names) <- "Genus"

depo.type <- data.frame(matrix(rep("Ss/De", max(depo$Code))))
colnames(depo.type) <- "Feeding.Type"

depo.MPB <- cbind(depo.names, depo.type, depo.MPB)

#Suspension feeders
susp.MPB <- data.frame(matrix(NA, nrow = max(susp$Code), ncol = 5))
colnames(susp.MPB) <- bpnames

for (i in 1:max(susp$Code)) {
  susp.MPB[i,] <- (susp.summ[[2]][[i]][2,])
}

susp.names <- data.frame(matrix(levels(susp$Genus), ncol = 1))
colnames(susp.names) <- "Genus"

susp.type <- data.frame(matrix(rep("Su/FF", max(susp$Code))))
colnames(susp.type) <- "Feeding.Type"

susp.MPB <- cbind(susp.names, susp.type, susp.MPB)


#Ep/Om
omni.MPB <- data.frame(matrix(NA, nrow = max(omni$Code), ncol = 5))
colnames(omni.MPB) <- bpnames

for (i in 1:max(omni$Code)) {
  omni.MPB[i,] <- (omni.summ[[2]][[i]][2,])
}

omni.names <- data.frame(matrix(levels(omni$Genus), ncol = 1))
colnames(omni.names) <- "Genus"

omni.type <- data.frame(matrix(rep("Ep/Om", max(omni$Code))))
colnames(omni.type) <- "Feeding.Type"

omni.MPB <- cbind(omni.names, omni.type, omni.MPB)


#Fish.MPB
fish.MPB <- data.frame(matrix(NA, nrow = max(fish$Code), ncol = 5))
colnames(fish.MPB) <- bpnames

for (i in 1:max(fish$Code)) {
  fish.MPB[i,] <- (fish.summ[[2]][[i]][2,])
}

fish.names <- data.frame(matrix(levels(fish$Genus), ncol = 1))
colnames(fish.names) <- "Genus"

fish.type <- data.frame(matrix(rep("Fish", max(fish$Code))))
colnames(fish.type) <- "Feeding.Type"

fish.MPB <- cbind(fish.names, fish.type, fish.MPB)


simmr.boxplots.MPB.all <- rbind(susp.MPB, depo.MPB, omni.MPB, fish.MPB)
write.xlsx(simmr.boxplots.MPB.all, "simmr.boxplots.MPB.all.xlsx")
#after xlsx saved, can load it next time without having to ru-run simmr
#simmr.boxplots.MPB.all <- read.csv("simmr.boxplots.MPB.all.csv", header = TRUE)

#Make palette to match previous biplots: Ss/De=red, Ep/Om=blue, Fish=Green, Su/FF=orange

boxpalette <- (c("Ss/De" = "#E41A1C",
                     "Ep/Om" = "#377EB8",
                     "Fish" = "#4DAF4A",
                     "Su/FF" = "#FF7F00"))


boxplot.MPB <- ggplot(data = simmr.boxplots.MPB.all, aes(x = Genus, fill = Feeding.Type))+
  geom_boxplot(aes(ymin = min, 
                   lower = low, 
                   middle = mid, 
                   upper = top, 
                   ymax = max),
               stat = "identity",
               color = "black")+
  scale_fill_manual(values = boxpalette,
                     name = "Trophic Guild",
                    breaks = c("Su/FF", "Ss/De", "Ep/Om", "Fish"))+
  coord_flip()+        
  scale_x_discrete(limits = rev(simmr.boxplots.MPB.all[,1]),
                   labels = rev(simmr.boxplots.MPB.all[,1]))+
  theme_bw()+
  theme(text = element_text(size = 14, color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(size = 16),
        legend.position = "none",
           panel.grid.major=element_blank(),
           panel.grid.minor=element_blank()
           )+
  ylab("MPB")

boxplot.MPB

############
##Shelf POM#
############

#Deposit Feeders

depo.Shelf <- data.frame(matrix(NA, nrow = max(depo$Code), ncol = 5))
bpnames <-  c("min", "low", "mid", "top", "max")
colnames(depo.Shelf) <- bpnames

for (i in 1:max(depo$Code)) {
  depo.Shelf[i,] <- (depo.summ[[2]][[i]][3,])
}

depo.names <- data.frame(matrix(levels(depo$Genus), ncol = 1))
colnames(depo.names) <- "Genus"

depo.type <- data.frame(matrix(rep("Ss/De", max(depo$Code))))
colnames(depo.type) <- "Feeding.Type"

depo.Shelf <- cbind(depo.names, depo.type, depo.Shelf)

#Suspension feeders
susp.Shelf <- data.frame(matrix(NA, nrow = max(susp$Code), ncol = 5))
colnames(susp.Shelf) <- bpnames

for (i in 1:max(susp$Code)) {
  susp.Shelf[i,] <- (susp.summ[[2]][[i]][3,])
}

susp.names <- data.frame(matrix(levels(susp$Genus), ncol = 1))
colnames(susp.names) <- "Genus"

susp.type <- data.frame(matrix(rep("Su/FF", max(susp$Code))))
colnames(susp.type) <- "Feeding.Type"

susp.Shelf <- cbind(susp.names, susp.type, susp.Shelf)


#Ep/Om
omni.Shelf <- data.frame(matrix(NA, nrow = max(omni$Code), ncol = 5))
colnames(omni.Shelf) <- bpnames

for (i in 1:max(omni$Code)) {
  omni.Shelf[i,] <- (omni.summ[[2]][[i]][3,])
}

omni.names <- data.frame(matrix(levels(omni$Genus), ncol = 1))
colnames(omni.names) <- "Genus"

omni.type <- data.frame(matrix(rep("Ep/Om", max(omni$Code))))
colnames(omni.type) <- "Feeding.Type"

omni.Shelf <- cbind(omni.names, omni.type, omni.Shelf)


#Fish.Shelf
fish.Shelf <- data.frame(matrix(NA, nrow = max(fish$Code), ncol = 5))
colnames(fish.Shelf) <- bpnames

for (i in 1:max(fish$Code)) {
  fish.Shelf[i,] <- (fish.summ[[2]][[i]][3,])
}

fish.names <- data.frame(matrix(levels(fish$Genus), ncol = 1))
colnames(fish.names) <- "Genus"

fish.type <- data.frame(matrix(rep("Fish", max(fish$Code))))
colnames(fish.type) <- "Feeding.Type"

fish.Shelf <- cbind(fish.names, fish.type, fish.Shelf)


simmr.boxplots.Shelf.all <- rbind(susp.Shelf, depo.Shelf, omni.Shelf, fish.Shelf)
write.xlsx(simmr.boxplots.Shelf.all, "simmr.boxplots.Shelf.all.xlsx")
#simmr.boxplots.Shelf.all <- read.csv("simmr.boxplots.Shelf.all.csv", header = TRUE)

boxplot.shelf <- ggplot(data = simmr.boxplots.Shelf.all, aes(x = Genus, fill = Feeding.Type))+
  geom_boxplot(aes(ymin = min, 
                   lower = low, 
                   middle = mid, 
                   upper = top, 
                   ymax = max),
               stat = "identity",
               color = "black")+
  scale_fill_manual(values = boxpalette,
                    name = "Trophic Guild")+
  coord_flip()+        
  scale_x_discrete(limits = rev(simmr.boxplots.Shelf.all[,1]),
                   labels = rev(simmr.boxplots.Shelf.all[,1]))+
  theme_bw()+
  theme(text = element_text(size = 14, color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(size = 16),
        axis.text.y = element_blank(),
        legend.position = "none",
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()
        )+
  ylab("Shelf POM")+
  xlab("")

boxplot.shelf

##################
##Terrestrial POM#
##################

#Deposit Feeders
depo.Terrestrial <- data.frame(matrix(NA, nrow = max(depo$Code), ncol = 5))
bpnames <-  c("min", "low", "mid", "top", "max")
colnames(depo.Terrestrial) <- bpnames

for (i in 1:max(depo$Code)) {
  depo.Terrestrial[i,] <- (depo.summ[[2]][[i]][4,])
}

depo.names <- data.frame(matrix(levels(depo$Genus), ncol = 1))
colnames(depo.names) <- "Genus"

depo.type <- data.frame(matrix(rep("Ss/De", max(depo$Code))))
colnames(depo.type) <- "Feeding.Type"

depo.Terrestrial <- cbind(depo.names, depo.type, depo.Terrestrial)

#Suspension feeders
susp.Terrestrial <- data.frame(matrix(NA, nrow = max(susp$Code), ncol = 5))
colnames(susp.Terrestrial) <- bpnames

for (i in 1:max(susp$Code)) {
  susp.Terrestrial[i,] <- (susp.summ[[2]][[i]][4,])
}

susp.names <- data.frame(matrix(levels(susp$Genus), ncol = 1))
colnames(susp.names) <- "Genus"

susp.type <- data.frame(matrix(rep("Su/FF", max(susp$Code))))
colnames(susp.type) <- "Feeding.Type"

susp.Terrestrial <- cbind(susp.names, susp.type, susp.Terrestrial)


#Ep/Om
omni.Terrestrial <- data.frame(matrix(NA, nrow = max(omni$Code), ncol = 5))
colnames(omni.Terrestrial) <- bpnames

for (i in 1:max(omni$Code)) {
  omni.Terrestrial[i,] <- (omni.summ[[2]][[i]][4,])
}

omni.names <- data.frame(matrix(levels(omni$Genus), ncol = 1))
colnames(omni.names) <- "Genus"

omni.type <- data.frame(matrix(rep("Ep/Om", max(omni$Code))))
colnames(omni.type) <- "Feeding.Type"

omni.Terrestrial <- cbind(omni.names, omni.type, omni.Terrestrial)


#Fish.Terrestrial
fish.Terrestrial <- data.frame(matrix(NA, nrow = max(fish$Code), ncol = 5))
colnames(fish.Terrestrial) <- bpnames

for (i in 1:max(fish$Code)) {
  fish.Terrestrial[i,] <- (fish.summ[[2]][[i]][4,])
}

fish.names <- data.frame(matrix(levels(fish$Genus), ncol = 1))
colnames(fish.names) <- "Genus"

fish.type <- data.frame(matrix(rep("Fish", max(fish$Code))))
colnames(fish.type) <- "Feeding.Type"

fish.Terrestrial <- cbind(fish.names, fish.type, fish.Terrestrial)


simmr.boxplots.Terrestrial.all <- rbind(susp.Terrestrial, depo.Terrestrial, omni.Terrestrial, fish.Terrestrial)
write.xlsx(simmr.boxplots.Terrestrial.all, "simmr.boxplots.Terrestrial.all.xlsx")
#simmr.boxplots.Terrestrial.all <- read.csv("simmr.boxplots.Terrestrial.all.csv", header = TRUE)

boxplot.terrestrial <- ggplot(data = simmr.boxplots.Terrestrial.all, aes(x = Genus, fill = Feeding.Type))+
  geom_boxplot(aes(ymin = min, 
                   lower = low, 
                   middle = mid, 
                   upper = top, 
                   ymax = max),
               stat = "identity",
               color = "black")+
  scale_fill_manual(values = boxpalette,
                    name = "Trophic Guild",
                    breaks = c("Su/FF", "Ss/De", "Ep/Om", "Fish"))+
  coord_flip()+        
  scale_x_discrete(limits = rev(simmr.boxplots.Terrestrial.all[,1]),
                   labels = rev(simmr.boxplots.Terrestrial.all[,1]))+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 0.75, 0.25))+
  theme_bw()+
  theme(text = element_text(size=14, color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(size = 16),
        axis.text.y = element_blank(),
        legend.text = element_text(size=18),
        legend.title = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )+
  ylab("Terrestrial POM")+
  xlab("")


boxplot.terrestrial

#library(egg)
grid.newpage()
png("P:/My Documents/Pubs/Harris et al/simmr/boxplot_panels.png", height = 6, width = 10, units="in", res = 600)
grid.draw(ggarrange(boxplot.MPB, boxplot.shelf, boxplot.terrestrial, ncol = 3))
dev.off()


#Mammals.MPB
#mamm.MPB <- data.frame(matrix(NA, nrow = max(mamm$Code), ncol = 5))
#colnames(mamm.MPB) <- bpnames
#for (i in 1:max(mamm$Code)) {
#  mamm.MPB[i,] <- (mamm.summ[[2]][[i]][2,])
#}
#mamm.names <- data.frame(matrix(levels(mamm$Genus), ncol = 1))
#colnames(mamm.names) <- "Genus"
#mamm.type <- data.frame(matrix(rep("Mam/Carn", max(mamm$Code))))
#colnames(mamm.type) <- "Feeding Type"
#mamm.MPB <- cbind(mamm.names, mamm.type, mamm.MPB)


beluga.all <- as.data.frame(beluga.summ$quantiles[1:3, 1:5])
colnames(beluga.all) <- bpnames
mamm.endmembers <- data.frame(matrix(c("Ice Algae", "Marine SPOM", "Terrestrial OM"), ncol = 1))
colnames(mamm.endmembers) <- "Endmember"

beluga.all <- cbind(mamm.endmembers, beluga.all)

ggplot(data = beluga.all, aes(x = Endmember))+
  geom_boxplot(aes(ymin = min, 
                   lower = low, 
                   middle = mid, 
                   upper = top, 
                   ymax = max),
               stat = "identity")+
  coord_flip()+        
  theme_bw()+
  theme(text = element_text(size=20, color = "black"),
        axis.text = element_text(color = "black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()
  )+
  ylab("Proportion assimilated")+
  xlab("End-member")


ggsave("beluga_boxplot.png", plot = last_plot(), device = "png", width = 8, height = 8, units = "in", dpi = 600)





###################################################################################################
###################################################################################################
###  Make stable isotope biplots for all genera analyzed in simmr                               ###
###  Code below is for all organisms, each trophic guild, and then faceted panel in manuscript  ###
###################################################################################################
###################################################################################################


###################
# Exclude mammals #
###################
simmriso_facet<-dplyr::filter(simmriso, Genus != "Delphinapterus")
simmriso_facet<-droplevels(simmriso_facet)

#shapes for facet plot Fig 4
shape_years2 <- c("2011" = 21,
                 "2012" = 24,
                 "2013" = 22,
                 "2014" = 3,
                 "2015" = 7,
                 "2016" = 23)

#Visualize each end-member's isospace
#Min and max extent of end-member polygons for Su/FF, Ss/De, and Ep/Om
sources.xmin <- ((sources$Meand13C - sources$SDd13C) + (corrections$Mean13c - corrections$sd13c))
sources.xmax <- ((sources$Meand13C + sources$SDd13C) + (corrections$Mean13c + corrections$sd13c))
sources.ymin <- ((sources$Meand15N - sources$SDd15N) + (corrections$Mean15N - corrections$sd15N))
sources.ymax <- ((sources$Meand15N + sources$SDd15N) + (corrections$Mean15N + corrections$sd15N))

#Min and max extent for end-member polygons for Fish
sourcesHTL.xmin <- ((sources$Meand13C - sources$SDd13C) + (corrections.HTL$Mean13c - corrections.HTL$sd13c))
sourcesHTL.xmax <- ((sources$Meand13C + sources$SDd13C) + (corrections.HTL$Mean13c + corrections.HTL$sd13c))
sourcesHTL.ymin <- ((sources$Meand15N - sources$SDd15N) + (corrections.HTL$Mean15N - corrections.HTL$sd15N))
sourcesHTL.ymax <- ((sources$Meand15N + sources$SDd15N) + (corrections.HTL$Mean15N + corrections.HTL$sd15N))

#Min and max extent for end-member polygons for Mam/Carn
sources.mamm.xmin <- ((sources.mamm$Meand13C - sources.mamm$SDd13C) + (corrections.mamm$Mean13c - corrections.mamm$sd13c))
sources.mamm.xmax <- ((sources.mamm$Meand13C + sources.mamm$SDd13C) + (corrections.mamm$Mean13c + corrections.mamm$sd13c))
sources.mamm.ymin <- ((sources.mamm$Meand15N - sources.mamm$SDd15N) + (corrections.mamm$Mean15N - corrections.mamm$sd15N))
sources.mamm.ymax <- ((sources.mamm$Meand15N + sources.mamm$SDd15N) + (corrections.mamm$Mean15N + corrections.mamm$sd15N))

#These will map out end-member polygons and labels for each facet
#Cheated and arranged these in Excel, but could extract them from above objects and simmriso_facet[,1]
MPB_polygons <- read.csv("MPB_polygons.csv", header = TRUE)
Shelf_polygons <- read.csv("Shelf_polygons.csv", header = TRUE)
Terrestrial_polygons <- read.csv("Terrestrial_polygons.csv", header = TRUE)
MPB_Label <- read.csv("MPB_Label.csv", header = TRUE)
POM_Label <- read.csv("POM_Label.csv", header = TRUE)
Terrestrial_Label <- read.csv("Terrestrial_Label.csv", header = TRUE)
Guild_Label <- read.csv("Guild_Label.csv", header = TRUE)

#rm(p)

p <- ggplot(data = simmriso_facet, aes(x = d13C, y = d15N, fill = Location, shape = Year))

p <- p + geom_point(size = 3, alpha = 0.8)+
  ylab(expression(paste(delta^{15},'N (\u2030)'))) +
  xlab(expression(paste(delta^{13},'C (\u2030)')))+
  theme_bw()+
  theme(axis.text = element_text(size=10, colour = "black"),
        axis.title = element_text(size = 16, color = "black"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 6.5))+
  coord_fixed(ratio=1)+
  scale_fill_manual(values = facet_palette)+
  scale_shape_manual(values = shape_years2)+
  guides(fill = guide_legend(override.aes = list(shape = 21)))+
  scale_x_continuous(limits = c(-32, -10), breaks = seq(-30, -10, 5))+
  scale_y_continuous(limits = c(2.5, 20), breaks = seq(5, 15, 5))+
  facet_wrap(~Genus, ncol = 4)

#add end-member polygons specific to each trophic guild
#polygon object has a column with same name and levels as faceted factor in above... e.g., facet_wrap(~Genus)
p <- p + geom_rect(data = MPB_polygons, inherit.aes = FALSE, 
                   aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
                   color = "black", fill = "gray", alpha = 0.2, size = 0.3)

p <- p + geom_rect(data = Shelf_polygons, inherit.aes = FALSE, 
                   aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
                   color = "black", fill = "gray", alpha = 0.2, size = 0.3)

p <- p + geom_rect(data = Terrestrial_polygons, inherit.aes = FALSE, 
                   aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
                   color = "black", fill = "gray", alpha = 0.2, size = 0.3)

#same idea as polygons, but now with text labels
p <- p + geom_text(data = MPB_Label, inherit.aes = FALSE, mapping = aes(x = x, y = y, label = Label.MPB), size = 2.5)

p <- p + geom_text(data = POM_Label, inherit.aes = FALSE, mapping = aes(x = x, y = y, label = Label.POM), size = 2.5)

p <- p + geom_text(data = Terrestrial_Label, inherit.aes = FALSE, mapping = aes(x = x, y = y, label = Label.Terrestrial), size = 2.5)

#can add trophic guild labels, but it gets busy!
#p <- p + geom_text(data = Guild_Label, inherit.aes = FALSE, mapping = aes(x = x, y = y, label = Label.Guild), size = 2)

p

ggsave("Fig 4_combined.tiff", plot = last_plot(), device = "tiff", width = 6.5, height = 7.5, units = "in", dpi = 600)

ggsave("Fig 4_combined.png", plot = last_plot(), device = "png", width = 6.5, height = 7.5, units = "in", dpi = 600)


##
#The rest of the plots are for each separate trophic guild, combined and faceted.
##

#Suspension Feeders
ggplot(data=susp, aes(x=d13C, y=d15N, color=Genus))+
  annotate("rect", xmin=sources.xmin[1], xmax=sources.xmax[1], ymin=sources.ymin[1], ymax=sources.ymax[1], alpha = 0.2, color = "black", fill = "grey")+
  annotate("text", x = -22, y = 12.4, label="Nearshore POM", size = 5)+
  annotate("rect", xmin=sources.xmin[2], xmax=sources.xmax[2], ymin=sources.ymin[2], ymax=sources.ymax[2], alpha = 0.2, color = "black", fill = "grey")+
  annotate("text", x = -15.2, y = 11.6, label="MPB", size = 5)+
  annotate("rect", xmin=sources.xmin[3], xmax=sources.xmax[3], ymin=sources.ymin[3], ymax=sources.ymax[3], alpha = 0.2, color = "black", fill = "grey")+
  annotate("text", x = -26.5, y = 4.4, label="Terrestrial POM", size=5)+
  geom_point(size = 4, alpha = 0.9)+
  ylab(expression(paste(delta^{15},'N (\u2030)'))) +
  xlab(expression(paste(delta^{13},'C (\u2030)')))+
  theme_bw()+
  theme(text = element_text(size=20),
        axis.text = element_text(size=20, color = "black"),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  coord_fixed(ratio = 1)+
  scale_color_manual(values = cbpalette)+
  scale_y_continuous(limits = c(3.5, 16), breaks = seq(5, 15, 5))+
  scale_x_continuous(limits = c(-30, -11), breaks = seq(-30, -15, 5))

ggsave("suspension_biplot.png", plot = last_plot(), device = "png", width = 10, height = 8, units = "in", dpi = 300)

#Separate genera
ggplot(data=susp, aes(x = d13C, y = d15N, shape = Year, color = Location))+
  annotate("rect", xmin=sources.xmin[1], xmax=sources.xmax[1], ymin=sources.ymin[1], ymax=sources.ymax[1], alpha = 0.2, color="black", fill="grey")+
  annotate("text", x = -22, y = 13.3, label="Nearshore POM", size = 4)+
  annotate("rect", xmin=sources.xmin[2], xmax=sources.xmax[2], ymin=sources.ymin[2], ymax=sources.ymax[2], alpha = 0.2, color="black", fill="grey")+
  annotate("text", x = -15.2, y = 11.4, label="MPB", size = 4)+
  annotate("rect", xmin=sources.xmin[3], xmax=sources.xmax[3], ymin=sources.ymin[3], ymax=sources.ymax[3], alpha = 0.2, color="black", fill="grey")+
  annotate("text", x = -26, y = 4.2, label="Terrestrial POM", size = 4)+
  geom_point(size = 4)+
  ylab(expression(paste(delta^{15},'N (\u2030)'))) +
  xlab(expression(paste(delta^{13},'C (\u2030)')))+
  theme_bw()+
  theme(text = element_text(size=20),
        axis.text = element_text(size=20, color = "black"),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        panel.grid.major=element_blank(),#removes panel grids
        panel.grid.minor=element_blank())+
  coord_fixed(ratio=1)+
  scale_y_continuous(limits=c(3.5, 16), breaks = seq(5, 15, 5))+
  scale_x_continuous(limits=c(-30, -11), breaks = seq(-30, -15, 5))+
  facet_wrap(~Genus, ncol = 2)+
  scale_color_manual(values = facet_palette)

ggsave("suspension_biplot_facet.png", plot = last_plot(), device = "png", width = 10, height = 8, units = "in", dpi = 600)


#Deposit Feeders
ggplot(data=depo, aes(x=d13C, y=d15N, color=Genus))+
  annotate("rect", xmin=sources.xmin[1], xmax=sources.xmax[1], ymin=sources.ymin[1], ymax=sources.ymax[1], alpha=0.2, color="black", fill="grey")+
  annotate("text", x=-22, y=12.4, label="Nearshore POM", size = 5)+
  annotate("rect", xmin=sources.xmin[2], xmax=sources.xmax[2], ymin=sources.ymin[2], ymax=sources.ymax[2], alpha=0.2, color="black", fill="grey")+
  annotate("text", x=-15.2, y=11.6, label="MPB", size = 5)+
  annotate("rect", xmin=sources.xmin[3], xmax=sources.xmax[3], ymin=sources.ymin[3], ymax=sources.ymax[3], alpha=0.2, color="black", fill="grey")+
  annotate("text", x=-26.5, y = 5, label="Terrestrial POM", size=5)+
  geom_point(size=4, alpha=.9)+
  ylab(expression(paste(delta^{15},'N (\u2030)'))) +
  xlab(expression(paste(delta^{13},'C (\u2030)')))+
  theme_bw()+
  theme(text = element_text(size=20),
        axis.text = element_text(size=20, color = "black"),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  coord_fixed(ratio=1)+
  scale_color_manual(values = cbpalette)+
  scale_y_continuous(limits=c(3.5, 16), breaks = seq(5, 15, 5))+
  scale_x_continuous(limits=c(-30, -11), breaks = seq(-30, -15, 5))

ggsave("deposit_biplot.png", plot = last_plot(), device = "png", width = 10, height = 8, units = "in", dpi = 300)

#Deposit Feeder facet
ggplot(data=depo, aes(x=d13C, y=d15N, shape = Year, color = Location))+
  annotate("rect", xmin=sources.xmin[1], xmax=sources.xmax[1], ymin=sources.ymin[1], ymax=sources.ymax[1], alpha = 0.2, color="black", fill="grey")+
  annotate("text", x = -22, y = 13.3, label="Nearshore POM", size = 4)+
  annotate("rect", xmin=sources.xmin[2], xmax=sources.xmax[2], ymin=sources.ymin[2], ymax=sources.ymax[2], alpha = 0.2, color="black", fill="grey")+
  annotate("text", x = -15.2, y = 11.4, label="MPB", size = 4)+
  annotate("rect", xmin=sources.xmin[3], xmax=sources.xmax[3], ymin=sources.ymin[3], ymax=sources.ymax[3], alpha = 0.2, color="black", fill="grey")+
  annotate("text", x = -26, y = 4.2, label="Terrestrial POM", size = 4)+
  geom_point(size=4)+
  ylab(expression(paste(delta^{15},'N (\u2030)'))) +
  xlab(expression(paste(delta^{13},'C (\u2030)')))+
  theme_bw()+
  theme(text = element_text(size=20),
        axis.text = element_text(size=16, color = "black"),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  coord_fixed(ratio=1)+
  scale_y_continuous(limits=c(3.5, 16), breaks = seq(5, 15, 5))+
  scale_x_continuous(limits=c(-30, -11), breaks = seq(-30, -15, 5))+
  facet_wrap(~Genus, ncol = 3)+
  scale_color_manual(values = facet_palette)

ggsave("deposit_biplot_facet.png", plot = last_plot(), device = "png", width = 10, height = 8, units = "in", dpi = 600)

#Omnivores
ggplot(data=omni, aes(x=d13C, y=d15N, color=Genus))+
  annotate("rect", xmin=sources.xmin[1], xmax=sources.xmax[1], ymin=sources.ymin[1], ymax=sources.ymax[1], alpha=0.2, color="black", fill="grey")+
  annotate("text", x=-22, y=12.4, label="Shelf POM", size = 5)+
  annotate("rect", xmin=sources.xmin[2], xmax=sources.xmax[2], ymin=sources.ymin[2], ymax=sources.ymax[2], alpha=0.2, color="black", fill="grey")+
  annotate("text", x=-15.2, y=11.6, label="MPB", size = 5)+
  annotate("rect", xmin=sources.xmin[3], xmax=sources.xmax[3], ymin=sources.ymin[3], ymax=sources.ymax[3], alpha=0.2, color="black", fill="grey")+
  annotate("text", x=-26.5, y = 5, label="Terrestrial POM", size=5)+
  geom_point(size=4, alpha=.9)+
  ylab(expression(paste(delta^{15},'N (\u2030)')))+
  xlab(expression(paste(delta^{13},'C (\u2030)')))+
  theme_bw()+
  theme(text = element_text(size=20),
        axis.text = element_text(size=20, color = "black"),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  coord_fixed(ratio=1)+
  scale_color_manual(values = cbpalette)+
  scale_y_continuous(limits=c(3.5, 16), breaks = seq(5, 15, 5))+
  scale_x_continuous(limits=c(-30, -11), breaks = seq(-30, -15, 5))

ggsave("omnivore_biplot.png", plot = last_plot(), device = "png", width = 10, height = 8, units = "in", dpi = 300)

#omnivore Separate by genus
ggplot(data=omni, aes(x=d13C, y=d15N, shape = Year, color = Location))+
  annotate("rect", xmin=sources.xmin[1], xmax=sources.xmax[1], ymin=sources.ymin[1], ymax=sources.ymax[1], alpha = 0.2, color="black", fill="grey")+
  annotate("text", x = -22, y = 13.4, label="Nearshore POM", size = 3)+
  annotate("rect", xmin=sources.xmin[2], xmax=sources.xmax[2], ymin=sources.ymin[2], ymax=sources.ymax[2], alpha = 0.2, color="black", fill="grey")+
  annotate("text", x = -15.2, y = 11.4, label="MPB", size = 3)+
  annotate("rect", xmin=sources.xmin[3], xmax=sources.xmax[3], ymin=sources.ymin[3], ymax=sources.ymax[3], alpha = 0.2, color="black", fill="grey")+
  annotate("text", x = -25, y = 4.1, label="Terrestrial POM", size = 3)+
  geom_point(size=4)+
  ylab(expression(paste(delta^{15},'N (\u2030)'))) +
  xlab(expression(paste(delta^{13},'C (\u2030)')))+
  theme_bw()+
  theme(text = element_text(size=20),
        axis.text = element_text(size=16, color = "black"),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  coord_fixed(ratio=1)+
  scale_y_continuous(limits=c(3.5, 16), breaks = seq(5, 15, 5))+
  scale_x_continuous(limits=c(-30, -11), breaks = seq(-30, -15, 5))+
  facet_wrap(~Genus, ncol = 4)+
  scale_color_manual(values = facet_palette)

ggsave("omnivore_biplot_facet.png", plot = last_plot(), device = "png", width = 10, height = 8, units = "in", dpi = 600)


#Fish
ggplot(data=fish, aes(x=d13C, y=d15N, color=Genus))+
  geom_point(size=4, alpha=.8)+
  ylab(expression(paste(delta^{15},'N (\u2030)'))) +
  xlab(expression(paste(delta^{13},'C (\u2030)')))+
  theme_bw()+
  theme(text = element_text(size=20),
        axis.text = element_text(size=20, color = "black"),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  coord_fixed(ratio=1)+
  scale_color_manual(values = cbpalette)+
  annotate("rect", xmin=sourcesHTL.xmin[1], xmax=sourcesHTL.xmax[1], ymin=sourcesHTL.ymin[1], ymax=sourcesHTL.ymax[1], alpha=0.2, color="black", fill="grey")+
  annotate("text", x=-21, y=17.2, label="Shelf POM", size=4)+
  annotate("rect", xmin=sourcesHTL.xmin[2], xmax=sourcesHTL.xmax[2], ymin=sourcesHTL.ymin[2], ymax=sourcesHTL.ymax[2], alpha=0.2, color="black", fill="grey")+
  annotate("text", x=-13, y = 17.3, label="MPB", size=4)+
  annotate("rect", xmin=sourcesHTL.xmin[3], xmax=sourcesHTL.xmax[3], ymin=sourcesHTL.ymin[3], ymax=sourcesHTL.ymax[3], alpha=0.2, color="black", fill="grey")+
  annotate("text", x=-25, y=8, label="Terrestrial POM", size=4)

ggsave("fish_biplot.png", plot = last_plot(), device = "png", width = 10, height = 8, units = "in", dpi = 300)

#Fish Facet
ggplot(data=fish, aes(x=d13C, y=d15N, shape = Year, color=Location))+
  annotate("rect", xmin=sourcesHTL.xmin[1], xmax=sourcesHTL.xmax[1], ymin=sourcesHTL.ymin[1], ymax=sourcesHTL.ymax[1], alpha=0.2, color="black", fill="grey")+
  annotate("text", x=-21.6, y=17.8, label="Shelf POM", size=4)+
  annotate("rect", xmin=sourcesHTL.xmin[2], xmax=sourcesHTL.xmax[2], ymin=sourcesHTL.ymin[2], ymax=sourcesHTL.ymax[2], alpha=0.2, color="black", fill="grey")+
  annotate("text", x=-13.5, y = 17.5, label="MPB", size=4)+
  annotate("rect", xmin=sourcesHTL.xmin[3], xmax=sourcesHTL.xmax[3], ymin=sourcesHTL.ymin[3], ymax=sourcesHTL.ymax[3], alpha=0.2, color="black", fill="grey")+
  annotate("text", x=-25, y=6.1, label="Terrestrial POM", size=4)+
  geom_point(size=4)+
  ylab(expression(paste(delta^{15},'N (\u2030)'))) +
  xlab(expression(paste(delta^{13},'C (\u2030)')))+
  theme_bw()+
  theme(text = element_text(size=20),
        axis.text = element_text(size=16, color = "black"),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        panel.grid.major=element_blank(),#removes panel grids
        panel.grid.minor=element_blank())+
  coord_fixed(ratio=1)+
  scale_y_continuous(limits=c(5, 22), breaks = seq(5, 20, 5))+
  scale_x_continuous(limits=c(-30, -11), breaks = seq(-30, -15, 5))+
  facet_wrap(~Genus, ncol = 3)+
  scale_color_manual(values = facet_palette)

ggsave("fish_biplot_facet.png", plot = last_plot(), device = "png", width = 10, height = 8, units = "in", dpi = 600)


#####
#Mammals
#####

ggplot(data=dplyr::filter(iso, Feeding.Type == "Mam/Carn"), 
       aes(x=d13C, y=d15N, color=Genus))+  
  annotate("rect", xmin=sources.mamm.xmin[1], xmax=sources.mamm.xmax[1], ymin=sources.mamm.ymin[1], ymax=sources.mamm.ymax[1], alpha=0.2, color="black", fill="grey")+
  annotate("text", x=-17.5, y=25.5, label="Marine POM", size=6)+
  annotate("rect", xmin=sources.mamm.xmin[2], xmax=sources.mamm.xmax[2], ymin=sources.mamm.ymin[2], ymax=sources.mamm.ymax[2], alpha=0.2, color="black", fill="grey")+
  annotate("text", x=-12, y=20, label="Ice Algae", size=6)+
  annotate("rect", xmin=sources.mamm.xmin[3], xmax=sources.mamm.xmax[3], ymin=sources.mamm.ymin[3], ymax=sources.mamm.ymax[3], alpha=0.2, color="black", fill="grey")+
  annotate("text", x=-25.1, y=22.5, label="Terrestrial POM", size=6)+
  geom_point(size = 5)+
  ylab(expression(paste(delta^{15},'N (\u2030)')))+
  xlab(expression(paste(delta^{13},'C (\u2030)')))+
  theme_bw()+
  theme(text = element_text(size=20),
        axis.text = element_text(size=20, color = "black"),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  coord_fixed(ratio=1)+
  scale_color_manual(values = facet_palette)

ggsave("mammal_biplot.png", plot = last_plot(), device = "png", width = 10, height = 8, units = "in", dpi = 300)

ggplot(data=dplyr::filter(iso, Feeding.Type == "Mam/Carn"),  
       aes(x=d13C, y=d15N, shape = Year, fill = Location, color = Location))+  
  annotate("rect", xmin=sources.mamm.xmin[1], xmax=sources.mamm.xmax[1], ymin=sources.mamm.ymin[1], ymax=sources.mamm.ymax[1], alpha=0.2, color="black", fill="grey")+
  annotate("text", x=-17.5, y=22.5, label="Marine POM", size=4)+
  annotate("rect", xmin=sources.mamm.xmin[2], xmax=sources.mamm.xmax[2], ymin=sources.mamm.ymin[2], ymax=sources.mamm.ymax[2], alpha=0.2, color="black", fill="grey")+
  annotate("text", x=-12.8, y=20, label="Ice Algae", size=4)+
  annotate("rect", xmin=sources.mamm.xmin[3], xmax=sources.mamm.xmax[3], ymin=sources.mamm.ymin[3], ymax=sources.mamm.ymax[3], alpha=0.2, color="black", fill="grey")+
  annotate("text", x=-25.5, y=22, label="Terrestrial POM", size=4)+
  geom_point(size=4)+
  ylab(expression(paste(delta^{15},'N (\u2030)'))) +
  xlab(expression(paste(delta^{13},'C (\u2030)')))+
  theme_bw()+
  theme(text = element_text(size=20),
        axis.text = element_text(size=20, color = "black"),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  coord_fixed(ratio=1)+
  facet_wrap(~Genus, ncol = 2)+
  scale_fill_manual(values = facet_palette)+
  scale_color_manual(values = facet_palette)+
  scale_shape_manual(values = shape_years2)

ggsave("mammal_biplot_facet.png", plot = last_plot(), device = "png", width = 10, height = 8, units = "in", dpi = 600)
