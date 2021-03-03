# Author: Bryden Fields
# Email: bryden.fields@york.ac.uk
# Last updated: 02/02/2021

# Initial setup -------------------------------------------------------------------

setwd('/Users/brydenfields/Documents/Publications/2021_Rhizobiuminteractions_paper/scripts/supernatant_scripts')

library(tidyverse)
library(RColorBrewer)
library(plyr); library(dplyr)
library(rcompanion)
library(car)
if(!require(psych)){install.packages("car")}
if(!require(MASS)){install.packages("MASS")}
if(!require(rcompanion)){install.packages("rcompanion")}
library(readxl)
library(reshape2)
library(pheatmap)
library(grid)
library(emmeans)

source('../summarySE.R')


#reading csv containing all data for 24 supernatant experiment
supernatant1 <- read.csv("../../data/raw_data/supernatant_data/bes_project_all_no1.csv")


# Calculate descriptive stats -------------------------------------------------------------------



OD_means <- summarySE(supernatant1, measurevar = 'OD', groupvars = c('inoculant', 'supernatant', 'time'))
write.csv(OD_means, file = '../../data/intermediate_data/supernatant_data/OD_supernatant_descriptstats.csv')


#Creating descriptive statistics summary of OD values for strains, genospecies and farm practise 
#and comp index descriptive stats for strains, genospecies and farmtypes 
#(farmtypes come within genospecies categorisation)
means_OD <- summarySE(supernatant1, measurevar = 'OD', groupvars = c('time', 'supernatant', 'inoculant'))
write.csv(means_OD, file = '../../data/intermediate_data/supernatant_data/OD_inoculant_descriptstats_no1.csv')

means_OD_genospecies <- summarySE(supernatant1, measurevar = 'OD', groupvars = c('time','sup_geno', 'inoc_geno'))
write.csv(means_OD_genospecies, file = '../../data/intermediate_data/supernatant_data/OD_geno_descriptstats_no1.csv')


means_comp <- summarySE(supernatant1, measurevar = 'comp_ind', groupvars = c('time', 'supernatant', 'inoculant'))
write.csv(means_comp, file = '../../data/intermediate_data/supernatant_data/compind_inoculant_descriptstats_no1.csv')

means_comp_genospecies <- summarySE(supernatant1, measurevar = 'comp_ind', groupvars = c('time','sup_geno', 'inoc_geno'))
write.csv(means_comp_genospecies, file = '../../data/intermediate_data/supernatant_data/supernatant1_compind_descriptstats_no1.csv')


# Graphs -------------------------------------------------------------------


#plot: inoculant growth curves grouped by genospecies/environment
(growth_geno <- ggplot(means_OD_genospecies, aes(x = time, y = OD, colour = inoc_geno))+
    geom_point() +
    geom_line(size = 1)+
    geom_errorbar(aes(ymin = OD - sd, ymax = OD + sd, width = 0.3))+
    labs(x = "Time (hrs)", y = "Optical Density (OD600)")+
    facet_wrap(~ sup_geno)+
    scale_colour_manual(values=c("#9ACD32", "#466EA9", "#4D9A7A", "#DA367E"))+
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme_bw()+
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 25),
          text = element_text(size=30),
          axis.text = element_text(size = 25)))
#function to save specific ggplot graph.  additional parameters(width=,height=,units='cm')
ggsave('../../data/intermediate_data/supernatant_data/growthcurves_genospecies.pdf', plot = growth_geno, width = 30, height = 25, units = 'cm' )



# filter for 62 hr timepoint data
supernatant1_62 <- supernatant1 %>% 
  filter(time %in% c('62'))


# calculate descriptive statistics for supernatant and inoculant genospecies groups. 
means_compheat_62 <- summarySE(supernatant1_62, measurevar = 'comp_ind', groupvars = c('time', 'sup_geno',
                                                                                       'supernatant', 'inoc_geno',
                                                                                       'inoculant'))


# add SM to strain names
means_compheat_62$supernatant <- paste0("SM", means_compheat_62$supernatant)
means_compheat_62$inoculant <- paste0("SM", means_compheat_62$inoculant)

# remove SM from the 100% and 50% TY treatments
means_compheat_62$supernatant[means_compheat_62$supernatant == 'SM100% TY'] <- '100% TY'
means_compheat_62$supernatant[means_compheat_62$supernatant == 'SM50% TY'] <- '50% TY'

# reformat supernatant levels and inoculant levels
means_compheat_62$supernatant <- factor(means_compheat_62$supernatant, levels = c("100% TY","50% TY","SM152B","SM137B","SM152A",
                                                                            "SM145B","SM154C", "SM144A", "SM147A", "SM158", 
                                                                            "SM170C", "SM157B", "SM165A", "SM122A", "SM126B",
                                                                            "SM149A", "SM135B", "SM135A", "SM159", "SM168A",
                                                                            "SM41", "SM53", "SM57", "SM74", "SM77", "SM67"))
means_compheat_62$inoculant <- factor(means_compheat_62$inoculant, levels = c("SM152B","SM137B","SM152A",
                                                                              "SM145B","SM154C", "SM144A", "SM147A", "SM158", 
                                                                              "SM170C", "SM157B", "SM165A", "SM122A", "SM126B",
                                                                              "SM149A", "SM135B", "SM135A", "SM159", "SM168A",
                                                                              "SM41", "SM53", "SM57", "SM74", "SM77", "SM67"))


# create horizontal (h) and vertical (v) lines for distinguishing genospecies groups in the heatmaps
h <- c(6.5, 13.5, 18.5)
v <- c(2.5, 8.5, 15.5, 20.5)

#plot: supernatant interactions heatmap
(compind_heat <- ggplot(means_compheat_62, aes(x = supernatant, y = inoculant, fill = comp_ind))+
    geom_tile(colour = "white", size=1) + 
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1)+
    labs(x = "Supernatant strain", y = "Inoculant strain", fill = "RGI")+
    geom_hline(yintercept = h, size = 0.75, alpha = 0.5)+
    geom_vline(xintercept = v, size = 0.75, alpha = 0.5)+
    #geom_vline(yintercept, linetype, color, size)
    scale_y_discrete(expand=c(0,0))+
    scale_x_discrete(expand=c(0,0))+
    coord_fixed()+
    guides(fill = guide_colourbar(barwidth = 2, barheight = 10)) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16),
          axis.text.y = element_text(size = 16), 
          axis.ticks=element_line(size=0.4), 
          # axis.title.x = element_text(size = 16), 
          # axis.title.y = element_text(size = 16),
          #legend.text = element_text(size = 16),
          panel.border=element_blank(), 
          text = element_text(size=25)))
#function to save specific ggplot graph.  additional parameters(width=,height=,units='cm')
ggsave('../../data/intermediate_data/supernatant_data/supernatant_interactions_heatmap.pdf', width = 10, height = 10, plot = compind_heat)








# Create a second heatmap for comparing growth in 50% TY to other strain supernatants 


# load RGI for 50% TY. Convert this into data we can use to compare against RGI of own supernatant.

supernatant_50TY <- read.csv("../../data/raw_data/supernatant_data/bes_project_all_50TY.csv")


#Creating descriptive statistics summary of OD values for strains, genospecies and farm practise 
#and comp index descriptive stats for strains, genospecies and farmtypes 
#(farmtypes come within genospecies categorisation)

means_50TY_comp <- summarySE(supernatant_50TY, measurevar = 'comp_ind_50TY', groupvars = c('time', 'supernatant', 'inoculant'))

means_50TY_comp_genospecies <- summarySE(supernatant_50TY, measurevar = 'comp_ind_50TY', groupvars = c('time','sup_geno', 'inoc_geno'))



#filter for 62 h time point
supernatant_50TY_62 <- supernatant_50TY %>% 
  filter(time %in% c('62'))

# make descriptive stats for inoculant strain and supernatant
means_50TY_compheat_62 <- summarySE(supernatant_50TY_62, measurevar = 'comp_ind_50TY', groupvars = c('time', 'sup_geno',
                                                                                                     'supernatant', 'inoc_geno',
                                                                                                     'inoculant'))

# reformat data
means_50TY_compheat_62$supernatant <- factor(means_50TY_compheat_62$supernatant, levels = c("100_TY","50_TY","152B","137B","152A",
                                                                                            "145B","154C", "144A", "147A", "158", 
                                                                                            "170C", "157B", "165A", "122A", "126B",
                                                                                            "149A", "135B", "135A", "159", "168A",
                                                                                            "41", "53", "57", "77", "74", "67"))
means_50TY_compheat_62$inoculant <- factor(means_50TY_compheat_62$inoculant, levels = c("152B","137B","152A",
                                                                                        "145B","154C", "144A", "147A", "158", 
                                                                                        "170C", "157B", "165A", "122A", "126B",
                                                                                        "149A", "135B", "135A", "159", "168A",
                                                                                        "41", "53", "57", "77", "74", "67"))

means_50TY_compheat_62$supernatant <- as.factor(means_50TY_compheat_62$supernatant)
means_50TY_compheat_62$inoculant <- as.factor(means_50TY_compheat_62$inoculant)

# vectors for splitting genospecies on heatmap horizontally (h) and vertically (v)
h <- c(6.5, 13.5, 18.5)
v <- c(2.5, 8.5, 15.5, 20.5)

#plot: 50% TY RGI heatmap
(compind_heat2 <- ggplot(means_50TY_compheat_62, aes(x = supernatant, y = inoculant, fill = comp_ind_50TY))+
    geom_tile(colour = "white", size=1) + 
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1)+
    labs(x = "Supernatant strain", y = "Inoculant strain")+
    geom_hline(yintercept = h, colour = "grey50", size = 0.75, alpha = 0.5)+
    geom_vline(xintercept = v, colour = "grey50", size = 0.75, alpha = 0.5)+
    #geom_vline(yintercept, linetype, color, size)
    scale_y_discrete(expand=c(0,0))+
    scale_x_discrete(expand=c(0,0))+
    coord_fixed()+
    theme_bw()+
    theme(legend.title = element_blank(), 
          axis.text.x = element_text(size = 16, angle = 90, hjust = 1, vjust = 0.5, colour = "grey50"),
          axis.text.y = element_text(size = 16, colour = "grey50"), 
          axis.ticks=element_line(size=0.4), 
          axis.title.x = element_text(size = 18), 
          axis.title.y = element_text(size = 18),
          legend.text = element_text(size = 16, colour = "grey50"),
          panel.border=element_blank(), text = element_text(size=10)))

#function to save specific ggplot graph.  additional parameters(width=,height=,units='cm')
ggsave('../../data/intermediate_data/supernatant_data/compind_heatmap_50TY.pdf', width = 10, height = 10, plot = compind_heat2)





# Descriptive statistics for genospecies inoculant and supernatant groups -----------------

# data using is: means_compheat_62
# dataset to make compind_heatmap_no1_v2.pdf
# using comp_ind column

# make histogram of the inoculant RGIs
hist(means_compheat_62$comp_ind)

(hist_plot <- ggplot(data = means_compheat_62, aes(x = comp_ind), col="red", fill="green", alpha = .2) + 
    geom_histogram(breaks = seq(0, 3, by = 0.01)) +
    labs(x = "Relative growth index of inoculant strains",
         y = "Count") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20), 
          axis.title.x = element_text(size = 20), 
          axis.title.y = element_text(size = 20))
)
#function to save specific ggplot graph.  additional parameters(width=,height=,units='cm')
ggsave('../../data/intermediate_data/supernatant_data/compind_62hheat_histogramdistribution.pdf', width = 10, height = 10, plot = hist_plot)



# what is the median and mean of the inoculant RGIs
mean(means_compheat_62$comp_ind)
median(means_compheat_62$comp_ind)

# how many RGI's less than 0.75 and RGIs more than 1.25?
nrow(means_compheat_62[means_compheat_62$comp_ind < 0.75, ])
nrow(means_compheat_62[means_compheat_62$comp_ind > 1.25, ])



# correlate ANI to RGI of inoculant
ANIdf <- as.matrix(read.csv('../../data/raw_data/supernatant_data/ani_sorted_by_genospecies_snps_new_data_6kgenes_supernatantsamples.csv', row.names = 1))

# melt and add column names for ANI
ANIdf_melt <- melt(ANIdf)
names(ANIdf_melt) <- c('supernatant', 'inoculant', 'ANI_6K')

# merge ANI with comp_ind RGI
ANI_compind <- merge(means_compheat_62[,c('supernatant', 'inoculant','sup_geno','inoc_geno', 'comp_ind')], 
                     ANIdf_melt, by = c('supernatant', 'inoculant'),all.x = T)

# remove supernatant that are 100% TY and 50% TY
ANI_compind <- ANI_compind[!ANI_compind$supernatant %in% c('100% TY', '50% TY'), ]

# or do it by distance from 1 (neutral interaction)
ANI_compind$dist_from_1 <- abs(ANI_compind$comp_ind - 1)


# correlate
plot(ANI_compind$ANI_6K, ANI_compind$comp_ind)


# correlate again 
plot(ANI_compind$ANI_6K, ANI_compind$dist_from_1,
     xlab = "Average Nucleotide Identity (6,529 genes, 441,287 SNPs)",
     ylab = "RGI distance from 1 (neutral interaction)",
     cex=1.5, cex.axis = 1.5,
     cex.lab = 1.5)

# pearsons correlation statistic
cor.test(ANI_compind$ANI_6K, ANI_compind$dist_from_1)



# remove organic or conventional notation so just genospecies category (gsA,gsC or gsE). 
ANI_compind$inoc_geno <- sub('.', '', ANI_compind$inoc_geno)
ANI_compind$sup_geno <- sub('.', '', ANI_compind$sup_geno)

# adding genospecies group
ANI_compind$geno_interaction <- ifelse(ANI_compind$inoc_geno == ANI_compind$sup_geno, 'within', 'between')
# correct for genospecies OC and CC = make just C so all C are within regardless of farm site. 

# make groups for each genospecies within/between RGIs
ANI_compind$geno_interaction_2 <- paste(ANI_compind$inoc_geno, ANI_compind$geno_interaction)

# mean pairwise strain ANIs for each genospecies group combos
ANI_compind %>%
  group_by(inoc_geno, sup_geno) %>% 
  summarise_at(vars("ANI_6K"), mean)


# model
lmMod <- lm(comp_ind ~ ANI_6K * inoc_geno * geno_interaction, data = ANI_compind)
summary(lmMod)

# Call:
#   lm(formula = comp_ind ~ ANI_6K * inoc_geno * geno_interaction, 
#      data = ANI_compind)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.51033 -0.10700 -0.02236  0.08984  1.10544 
# 
# Coefficients:
#                                             Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)                               -264.16      48.40  -5.458 7.34e-08 ***
#   ANI_6K                                     293.45      53.60   5.475 6.70e-08 ***
#   inoc_genoC                                 138.73      64.92   2.137   0.0330 *  
#   inoc_genoE                                 322.22      68.44   4.708 3.18e-06 ***
#   geno_interactionwithin                     264.49      48.49   5.454 7.51e-08 ***
#   ANI_6K:inoc_genoC                         -153.49      71.88  -2.135   0.0332 *  
#   ANI_6K:inoc_genoE                         -356.50      75.81  -4.702 3.27e-06 ***
#   ANI_6K:geno_interactionwithin             -292.80      53.69  -5.454 7.53e-08 ***
#   inoc_genoC:geno_interactionwithin         -138.56      65.01  -2.131   0.0335 *  
#   inoc_genoE:geno_interactionwithin         -321.67      68.78  -4.677 3.69e-06 ***
#   ANI_6K:inoc_genoC:geno_interactionwithin   153.38      71.97   2.131   0.0335 *  
#   ANI_6K:inoc_genoE:geno_interactionwithin   355.99      76.13   4.676 3.70e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1778 on 540 degrees of freedom
# Multiple R-squared:  0.2961,	Adjusted R-squared:  0.2817 
# F-statistic: 20.65 on 11 and 540 DF,  p-value: < 2.2e-16


# checking assumptions
plot(lmMod)

# anova version of model
lmMod.aov <- aov(lmMod)
summary(lmMod.aov)




# model no interaction for geno_interaction
lmMod_noint <- lm(comp_ind ~ ANI_6K + inoc_geno * geno_interaction, data = ANI_compind) 
summary(lmMod_noint)

# anova to compare model with and without interaction
anova(lmMod_noint, lmMod)

# Analysis of Variance Table
# 
# Model 1: comp_ind ~ ANI_6K + inoc_geno * geno_interaction
# Model 2: comp_ind ~ ANI_6K * inoc_geno * geno_interaction
# Res.Df    RSS Df Sum of Sq      F    Pr(>F)    
# 1    545 18.334                                  
# 2    540 17.078  5    1.2559 7.9424 3.099e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# model no interaction for geno_interaction
lmMod_noint2 <- lm(comp_ind ~ ANI_6K + inoc_geno * geno_interaction, data = ANI_compind) 
summary(lmMod_noint2)

# anova to compare model with and without interaction
anova(lmMod_noint2, lmMod)

# Analysis of Variance Table
# 
# Model 1: comp_ind ~ ANI_6K + inoc_geno * geno_interaction
# Model 2: comp_ind ~ ANI_6K * inoc_geno * geno_interaction
# Res.Df    RSS Df Sum of Sq      F    Pr(>F)    
# 1    545 18.334                                  
# 2    540 17.078  5    1.2559 7.9424 3.099e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# genospecies palette
geno_palette <- c("#466EA9","#4D9A7A", "#DA367E")


# make new labels
ANI_compind$geno_interaction_3 <- ANI_compind$geno_interaction_2

ANI_compind$geno_interaction_3 <- factor(ANI_compind$geno_interaction_3, 
                                         levels = c("A between", "A within", 
                                                    "C between", "C within", 
                                                    "E between", "E within"),
                                         labels = c("gsA - other gs", "gsA - gsA", 
                                                    "gsC - other gs", "gsC - gsC", 
                                                    "gsE - other gs", "gsE - gsE"))

ANI_compind$geno_interaction <- factor(ANI_compind$geno_interaction, 
                                         levels = c("between", "within"),
                                         labels = c("between genospecies", "within genospecies"))


# plot: ANI vs RGI as inoculant. facet graph for between genospecies and within genospecies interactions
(graph_lmMod_facet <- ggplot(ANI_compind, aes(x = ANI_6K, y = comp_ind,
                                                 colour = inoc_geno)) +
    geom_point(alpha = 0.1, size = 3) +
    geom_smooth(method = lm, se = FALSE) +
    scale_colour_manual(values = geno_palette) +
    facet_wrap(. ~ geno_interaction, scales = "free_x") +
    labs(x = "Average Nucleotide Identity (6,529 genes, 441,287 SNPs)", 
         y = "Relative Growth Index (RGI)",
         colour = "Inoculant\ngenospecies") +
    theme_bw() +
    theme(legend.title = element_text(size = 14),
          legend.text=element_text(size=14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          strip.text.x = element_text(size = 16)))
#function to save specific ggplot graph.  additional parameters(width=,height=,units='cm')
ggsave('../../data/intermediate_data/supernatant_data/ANIvsRGI_genospeciesfacet.pdf', plot = graph_lmMod_facet, height = 20, width = 25, units = 'cm')


# replot data without SM168A and SM154C.

# remove rows that have 'SM168A' or 'SM154C' in any of the inoculant or supernatant columns
ANI_compind_168A154Crm <- ANI_compind %>% 
  filter(!supernatant %in% c('SM168A', 'SM154C')) %>%
  filter(!inoculant %in% c('SM168A', 'SM154C'))


# replot: ANI vs RGI as inoculant without SM168A and SM154C. facet graph for between genospecies and within genospecies interactions
(graph_lmMod_facet_168A154Crm <- ggplot(ANI_compind_168A154Crm, aes(x = ANI_6K, y = comp_ind,
                                               colour = inoc_geno)) +
    geom_point(alpha = 0.1, size = 3) +
    geom_smooth(method = lm, se = FALSE) +
    scale_colour_manual(values = geno_palette) +
    facet_wrap(. ~ geno_interaction, scales = "free_x") +
    labs(x = "Average Nucleotide Identity (6,529 genes, 441,287 SNPs)", 
         y = "Relative Growth Index (RGI)",
         colour = "Inoculant\ngenospecies") +
    theme_bw() +
    theme(legend.title = element_text(size = 14),
          legend.text=element_text(size=14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          strip.text.x = element_text(size = 16)))
#function to save specific ggplot graph.  additional parameters(width=,height=,units='cm')
ggsave('../../data/intermediate_data/supernatant_data/ANIvsRGI_genospeciesfacet_168A154Crm.pdf', plot = graph_lmMod_facet_168A154Crm, height = 20, width = 25, units = 'cm')


# emmeans RGI vs ANI, for each genospecies and between/within combinations --------------

#library(emmeans)

# using full model (lmMod)
# full model again is:
# lmMod <- lm(comp_ind ~ ANI_6K * inoc_geno * geno_interaction, data = ANI_compind)
summary(lmMod)

#summary(emtrends(lmMod, pairwise ~ inoc_geno*geno_interaction, var = "ANI_6K"))

# because we are interested its linear effect (slope), we use the emtrends function to estimate the slope in each species individually. In terms of simple slopes, we test whether the between and within interactions slopes are non-zero in each genospecies. The infer argument in the summary of emtrends requests t-tests and p-values for the slopes.
# in the result output ANI_6K.trend refers to the slope value. 
summary(emtrends(lmMod, pairwise ~ inoc_geno*geno_interaction, var = "ANI_6K"), infer=TRUE)

# $emtrends
# inoc_geno geno_interaction ANI_6K.trend    SE  df lower.CL upper.CL t.ratio p.value
# A         between               293.454 53.60 540   188.17   398.74  5.475  <.0001 ***
# C         between               139.966 47.90 540    45.87   234.06  2.922  0.0036 **
# E         between               -63.050 53.62 540  -168.38    42.27 -1.176  0.2401 
# A         within                  0.654  3.16 540    -5.54     6.85  0.207  0.8358 
# C         within                  0.541  1.61 540    -2.63     3.71  0.335  0.7375 
# E         within                  0.138  6.19 540   -12.02    12.29  0.022  0.9822 
# 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast              estimate    SE  df lower.CL upper.CL t.ratio p.value
# A between - C between  153.489 71.88 540   -52.10    359.1  2.135  0.2707 
# A between - E between  356.505 75.81 540   139.68    573.3  4.702  <.0001 ***
# A between - A within   292.800 53.69 540   139.25    446.4  5.454  <.0001 ***
# A between - C within   292.914 53.62 540   139.56    446.3  5.463  <.0001 ***
# A between - E within   293.316 53.95 540   139.01    447.6  5.437  <.0001 ***
# C between - E between  203.016 71.90 540    -2.62    408.6  2.824  0.0553 .
# C between - A within   139.311 48.01 540     2.01    276.6  2.902  0.0444 *
# C between - C within   139.425 47.93 540     2.35    276.5  2.909  0.0436 *
# C between - E within   139.827 48.30 540     1.69    278.0  2.895  0.0453 *
# E between - A within   -63.705 53.71 540  -217.32     89.9 -1.186  0.8435 
# E between - C within   -63.591 53.64 540  -217.01     89.8 -1.185  0.8438 
# E between - E within   -63.189 53.97 540  -217.55     91.2 -1.171  0.8507 
# A within - C within      0.113  3.54 540   -10.02     10.2  0.032  1.0000 
# A within - E within      0.516  6.95 540   -19.35     20.4  0.074  1.0000 
# C within - E within      0.403  6.39 540   -17.89     18.7  0.063  1.0000 
# 
# Confidence level used: 0.95 
# Conf-level adjustment: tukey method for comparing a family of 6 estimates 
# P value adjustment: tukey method for comparing a family of 6 estimates 





# rerun slope analysis on dataset without SM168A and SM154C. 
# ANI_compind_168A154Crm dataset

# we must remake the model but with the SM168A and SM154C removed dataset
lmMod_168A154Crm <- lm(comp_ind ~ ANI_6K * inoc_geno * geno_interaction, data = ANI_compind_168A154Crm)
summary(lmMod_168A154Crm)

#summary(emtrends(lmMod_168A154Crm, pairwise ~ inoc_geno*geno_interaction, var = "ANI_6K"))

# because we are interested its linear effect (slope), we use the emtrends function to estimate the slope in each species individually. In terms of simple slopes, we test whether the between and within interactions slope is non-zero in each genospecies. The infer argument in the summary of emtrends requests t-tests and p-values for the slopes
# in the result output ANI_6K.trend refers to the slope value. 
summary(emtrends(lmMod_168A154Crm, pairwise ~ inoc_geno*geno_interaction, var = "ANI_6K"), infer=TRUE)

# $emtrends
# inoc_geno geno_interaction    ANI_6K.trend    SE  df lower.CL upper.CL t.ratio p.value
# A         between genospecies      191.038 53.66 450    85.58   296.50  3.560  0.0004 ***
# C         between genospecies      138.262 42.74 450    54.27   222.25  3.235  0.0013 **
# E         between genospecies      -31.103 54.30 450  -137.82    75.61 -0.573  0.5671 
# A         within genospecies         3.068  3.20 450    -3.23     9.36  0.958  0.3386 
# C         within genospecies         0.541  1.33 450    -2.07     3.15  0.407  0.6839 
# E         within genospecies        -0.842  6.34 450   -13.31    11.62 -0.133  0.8945 
# 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                                      estimate    SE  df lower.CL upper.CL t.ratio p.value
# A between genospecies - C between genospecies    52.78 68.60 450  -143.56    249.1  0.769  0.9725 
# A between genospecies - E between genospecies   222.14 76.34 450     3.65    440.6  2.910  0.0437 *
# A between genospecies - A within genospecies    187.97 53.76 450    34.11    341.8  3.497  0.0068 **
# A between genospecies - C within genospecies    190.50 53.68 450    36.87    344.1  3.549  0.0057 **
# A between genospecies - E within genospecies    191.88 54.04 450    37.23    346.5  3.551  0.0056 **
# C between genospecies - E between genospecies   169.36 69.10 450   -28.40    367.1  2.451  0.1413 
# C between genospecies - A within genospecies    135.19 42.86 450    12.54    257.9  3.155  0.0211 *
# C between genospecies - C within genospecies    137.72 42.76 450    15.35    260.1  3.221  0.0171 *
# C between genospecies - E within genospecies    139.10 43.21 450    15.45    262.8  3.220  0.0172 *
# E between genospecies - A within genospecies    -34.17 54.39 450  -189.85    121.5 -0.628  0.9889 
# E between genospecies - C within genospecies    -31.64 54.32 450  -187.10    123.8 -0.583  0.9922 
# E between genospecies - E within genospecies    -30.26 54.67 450  -186.73    126.2 -0.554  0.9938 
# A within genospecies - C within genospecies       2.53  3.47 450    -7.39     12.4  0.729  0.9783 
# A within genospecies - E within genospecies       3.91  7.11 450   -16.43     24.2  0.550  0.9940 
# C within genospecies - E within genospecies       1.38  6.48 450   -17.17     19.9  0.213  0.9999 
# 
# Confidence level used: 0.95 
# Conf-level adjustment: tukey method for comparing a family of 6 estimates 
# P value adjustment: tukey method for comparing a family of 6 estimates 




