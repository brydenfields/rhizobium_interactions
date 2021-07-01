# Author: Bryden Fields
# Email: bryden.fields@york.ac.uk
# Last updated: 01/07/2021

# Initial setup -------------------------------------------------------------------

#set the working directory to the directory this script is in.
#setwd("~/what/ever/folder/you/are/working/from") 

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
# [1] 1.005847
median(means_compheat_62$comp_ind)
# [1] 1.003787

# how many RGI's less than 0.75 and RGIs more than 1.25?
nrow(means_compheat_62[means_compheat_62$comp_ind < 0.75, ])
# [1] 55
nrow(means_compheat_62[means_compheat_62$comp_ind > 1.25, ])
# [1] 55


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
# Pearson's product-moment correlation
# 
# data:  ANI_compind$ANI_6K and ANI_compind$dist_from_1
# t = -3.8207, df = 550, p-value = 0.0001482
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.24101588 -0.07839274
# sample estimates:
#        cor 
# -0.1607954



# adding genospecies group
ANI_compind$geno_interaction <- ifelse(ANI_compind$inoc_geno == ANI_compind$sup_geno, 'within', 'between')
# correct for genospecies OC and CC = make just C so all C are within regardless of farm site. 

# make groups for each genospecies within/between RGIs
ANI_compind$geno_interaction_2 <- paste(ANI_compind$inoc_geno, ANI_compind$geno_interaction)

# mean pairwise strain ANIs for each genospecies group combos
ANI_compind %>%
  group_by(inoc_geno, sup_geno) %>% 
  summarise_at(vars("ANI_6K"), mean)

#  relevel factors for model so OA is the intercept
ANI_compind$inoc_geno <- factor(ANI_compind$inoc_geno,
                                levels = c("OA", "OC", "OE", "CC"))
ANI_compind$sup_geno <- factor(ANI_compind$sup_geno,
                                levels = c("OA", "OC", "OE", "CC"))

# model
lmMod <- lm(comp_ind ~ ANI_6K * inoc_geno * geno_interaction, data = ANI_compind)
summary(lmMod)

# Call:
#   lm(formula = comp_ind ~ ANI_6K * inoc_geno * geno_interaction, 
#      data = ANI_compind)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.43215 -0.11171 -0.01716  0.08896  1.10544 
# 
# Coefficients:
#                                             Estimate Std. Error t value Pr(>|t|)    
#  (Intercept)                                -264.16      48.32  -5.467 7.02e-08 ***
#   ANI_6K                                      293.45      53.51   5.484 6.41e-08 ***
#   inoc_genoOC                                 263.83      48.32   5.460 7.28e-08 ***
#   inoc_genoOE                                 322.22      68.33   4.716 3.08e-06 ***
#   inoc_genoCC                                 265.79      48.32   5.501 5.86e-08 ***
#   geno_interactionwithin                      264.49      48.42   5.463 7.19e-08 ***
#   ANI_6K:inoc_genoOC                         -292.01      53.51  -5.457 7.41e-08 ***
#   ANI_6K:inoc_genoOE                         -356.50      75.69  -4.710 3.16e-06 ***
#   ANI_6K:inoc_genoCC                         -294.17      53.51  -5.497 5.97e-08 ***
#   ANI_6K:geno_interactionwithin              -292.80      53.60  -5.462 7.20e-08 ***
#   inoc_genoOC:geno_interactionwithin         -261.34      48.49  -5.390 1.06e-07 ***
#   inoc_genoOE:geno_interactionwithin         -321.67      68.67  -4.684 3.56e-06 ***
#   inoc_genoCC:geno_interactionwithin         -265.43      49.06  -5.410 9.49e-08 ***
#   ANI_6K:inoc_genoOC:geno_interactionwithin   289.53      53.67   5.394 1.03e-07 ***
#   ANI_6K:inoc_genoOE:geno_interactionwithin   355.99      76.01   4.684 3.57e-06 ***
#   ANI_6K:inoc_genoCC:geno_interactionwithin   293.87      54.20   5.422 8.92e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1775 on 536 degrees of freedom
# Multiple R-squared:  0.3036,	Adjusted R-squared:  0.2841 
# F-statistic: 15.58 on 15 and 536 DF,  p-value: < 2.2e-16


# checking assumptions
plot(lmMod)

# anova version of model
lmMod.aov <- aov(lmMod)
summary(lmMod.aov)

#                                      Df Sum Sq Mean Sq F value   Pr(>F)    
#   ANI_6K                              1  0.138  0.1379   4.373   0.0370 *  
#   inoc_geno                           3  5.069  1.6895  53.596  < 2e-16 ***
#   geno_interaction                    1  0.009  0.0089   0.284   0.5944    
#   ANI_6K:inoc_geno                    3  0.919  0.3063   9.715 2.98e-06 ***
#   ANI_6K:geno_interaction             1  0.000  0.0000   0.001   0.9821    
#   inoc_geno:geno_interaction          3  0.254  0.0845   2.681   0.0462 *  
#   ANI_6K:inoc_geno:geno_interaction   3  0.977  0.3256  10.330 1.28e-06 ***
#   Residuals                         536 16.896  0.0315                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# model no interaction for ANI_6K interaction
lmMod_noint <- lm(comp_ind ~ ANI_6K + inoc_geno * geno_interaction, data = ANI_compind) 
summary(lmMod_noint)

# anova to compare model with and without interaction
anova(lmMod_noint, lmMod)

# Analysis of Variance Table
# 
# Model 1: comp_ind ~ ANI_6K + inoc_geno * geno_interaction
# Model 2: comp_ind ~ ANI_6K * inoc_geno * geno_interaction
#   Res.Df    RSS Df Sum of Sq      F    Pr(>F)    
# 1    543 18.194                                  
# 2    536 16.896  7    1.2972 5.8785 1.357e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# model no interaction for geno_interaction
lmMod_noint2 <- lm(comp_ind ~ ANI_6K * inoc_geno + geno_interaction, data = ANI_compind) 
summary(lmMod_noint2)

# anova to compare model with and without interaction
anova(lmMod_noint2, lmMod)

# Analysis of Variance Table
# 
# Model 1: comp_ind ~ ANI_6K * inoc_geno + geno_interaction
# Model 2: comp_ind ~ ANI_6K * inoc_geno * geno_interaction
# Res.Df    RSS Df Sum of Sq      F    Pr(>F)    
# 1    543 18.127                                  
# 2    536 16.896  7    1.2305 5.5763 3.232e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# genospecies palette
geno_palette <- c("#466EA9", "#4D9A7A", "#DA367E", "#9ACD32")


# make new labels
ANI_compind$geno_interaction_3 <- ANI_compind$geno_interaction_2

ANI_compind$geno_interaction_3 <- factor(ANI_compind$geno_interaction_3, 
                                         levels = c("OA between", "OA within", 
                                                    "OC between", "OC within", 
                                                    "OE between", "OE within",
                                                    "CC between", "CC within"),
                                         labels = c("OA - other gs group", "OA - OA", 
                                                    "OC - other gs group", "OC - OC", 
                                                    "OE - other gs group", "OE - OE",
                                                    "CC - other gs group", "CC - CC"))

ANI_compind$geno_interaction <- factor(ANI_compind$geno_interaction, 
                                         levels = c("between", "within"),
                                         labels = c("between genospecies group", "within genospecies group"))


# find minimum, maximum and mean "between genospecies" ANI
min(ANI_compind[ANI_compind$geno_interaction == "within genospecies group", c("ANI_6K")])
# [1] 0.968164
max(ANI_compind[ANI_compind$geno_interaction == "within genospecies group", c("ANI_6K")])
# [1] 0.9999338
mean(ANI_compind[ANI_compind$geno_interaction == "within genospecies group", c("ANI_6K")])
# [1] 0.9829003

# find minimum, maximum and between "within genospecies" ANI
min(ANI_compind[ANI_compind$geno_interaction == "between genospecies group", c("ANI_6K")])
# [1] 0.9017105
max(ANI_compind[ANI_compind$geno_interaction == "between genospecies group", c("ANI_6K")])
# [1] 0.9767628
mean(ANI_compind[ANI_compind$geno_interaction == "between genospecies group", c("ANI_6K")])
# [1] 0.9162912

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


# plot: ANI vs RGI as inoculant. facet graph for between genospecies and within genospecies interactions
# zoomed in on between genospecies
(graph_lmMod_facet2 <- ggplot(ANI_compind[ANI_compind$geno_interaction == "between genospecies group",], aes(x = ANI_6K, y = comp_ind,
                                              colour = inoc_geno)) +
    geom_point(alpha = 0.1, size = 3) +
    geom_smooth(method = lm, se = FALSE) +
    scale_colour_manual(values = geno_palette) +
    labs(x = "Average Nucleotide Identity (6,529 genes, 441,287 SNPs)", 
         y = "Relative Growth Index (RGI)",
         colour = "Inoculant\ngenospecies") +
    coord_cartesian(xlim = c(0.9015, 0.904)) +
    theme_bw() +
    theme(legend.title = element_text(size = 14),
          legend.text=element_text(size=14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          strip.text.x = element_text(size = 16)))
#function to save specific ggplot graph.  additional parameters(width=,height=,units='cm')
ggsave('../../data/intermediate_data/supernatant_data/ANIvsRGI_genospeciesfacet_xlimedit.pdf', plot = graph_lmMod_facet2, height = 20, width = 25, units = 'cm')





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


# replot: ANI vs RGI as inoculant without SM168A and SM154C. facet graph for between genospecies and within genospecies interactions
# zoomed in on between genospecies
(graph_lmMod_facet_168A154Crm2 <- ggplot(ANI_compind_168A154Crm[ANI_compind_168A154Crm$geno_interaction == "between genospecies group",], aes(x = ANI_6K, y = comp_ind, colour = inoc_geno)) +
    geom_point(alpha = 0.1, size = 3) +
    geom_smooth(method = lm, se = FALSE) +
    scale_colour_manual(values = geno_palette) +
    labs(x = "Average Nucleotide Identity (6,529 genes, 441,287 SNPs)", 
         y = "Relative Growth Index (RGI)",
         colour = "Inoculant\ngenospecies") +
    xlim(0.9015, 0.904) +
    theme_bw() +
    theme(legend.title = element_text(size = 14),
          legend.text=element_text(size=14),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          strip.text.x = element_text(size = 16)))
#function to save specific ggplot graph.  additional parameters(width=,height=,units='cm')
ggsave('../../data/intermediate_data/supernatant_data/ANIvsRGI_genospeciesfacet_168A154Crm_xlimedit.pdf', plot = graph_lmMod_facet_168A154Crm2, height = 20, width = 25, units = 'cm')


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
# inoc_geno geno_interaction ANI_6K.trend     SE  df lower.CL upper.CL t.ratio p.value
# OA        between               293.454 53.510 536   188.34  398.569  5.484  <.0001 ***
# OC        between                 1.443  0.500 536     0.46    2.425  2.883  0.0041 **
# OE        between               -63.050 53.530 536  -168.21   42.105 -1.178  0.2394 
# CC        between                -0.718  0.514 536    -1.73    0.292 -1.397  0.1630 
# OA        within                  0.654  3.150 536    -5.53    6.842  0.208  0.8356 
# OC        within                 -1.823  2.694 536    -7.12    3.470 -0.677  0.4990 
# OE        within                  0.138  6.178 536   -12.00   12.273  0.022  0.9822 
# CC        within                  0.352  7.997 536   -15.36   16.062  0.044  0.9649 
# 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                estimate     SE  df  lower.CL upper.CL t.ratio p.value
# OA between - OC between  292.012 53.512 536  129.1784   454.85  5.457  <.0001 ***
# OA between - OE between  356.505 75.689 536  126.1890   586.82  4.710  0.0001 ***
# OA between - CC between  294.173 53.512 536  131.3390   457.01  5.497  <.0001 ***
# OA between - OA within   292.800 53.602 536  129.6920   455.91  5.462  <.0001 ***
# OA between - OC within   295.277 53.577 536  132.2448   458.31  5.511  <.0001 ***
# OA between - OE within   293.316 53.865 536  129.4084   457.22  5.445  <.0001 ***
# OA between - CC within   293.103 54.104 536  128.4678   457.74  5.417  <.0001 ***
# OC between - OE between   64.493 53.533 536  -98.4035   227.39  1.205  0.9304 
# OC between - CC between    2.161  0.718 536   -0.0224     4.34  3.012  0.0547 .
# OC between - OA within     0.788  3.190 536   -8.9176    10.49  0.247  1.0000 
# OC between - OC within     3.265  2.740 536   -5.0733    11.60  1.192  0.9342 
# OC between - OE within     1.304  6.198 536  -17.5551    20.16  0.210  1.0000 
# OC between - CC within     1.091  8.013 536  -23.2919    25.47  0.136  1.0000 
# OE between - CC between  -62.332 53.533 536 -225.2289   100.56 -1.164  0.9416 
# OE between - OA within   -63.705 53.623 536 -226.8759    99.47 -1.188  0.9352 
# OE between - OC within   -61.227 53.598 536 -224.3232   101.87 -1.142  0.9471 
# OE between - OE within   -63.189 53.886 536 -227.1592   100.78 -1.173  0.9394 
# OE between - CC within   -63.402 54.124 536 -228.0995   101.29 -1.171  0.9397 
# CC between - OA within    -1.373  3.192 536  -11.0854     8.34 -0.430  0.9999 
# CC between - OC within     1.105  2.743 536   -7.2423     9.45  0.403  0.9999 
# CC between - OE within    -0.857  6.199 536  -19.7196    18.01 -0.138  1.0000 
# CC between - CC within    -1.070  8.014 536  -25.4556    23.32 -0.134  1.0000 
# OA within - OC within      2.477  4.145 536  -10.1366    15.09  0.598  0.9989 
# OA within - OE within      0.516  6.934 536  -20.5850    21.62  0.074  1.0000 
# OA within - CC within      0.302  8.595 536  -25.8526    26.46  0.035  1.0000 
# OC within - OE within     -1.961  6.740 536  -22.4692    18.55 -0.291  1.0000 
# OC within - CC within     -2.175  8.439 536  -27.8538    23.50 -0.258  1.0000 
# OE within - CC within     -0.214 10.105 536  -30.9635    30.54 -0.021  1.0000 
# 
# Confidence level used: 0.95 
# Conf-level adjustment: tukey method for comparing a family of 8 estimates 
# P value adjustment: tukey method for comparing a family of 8 estimates 




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
# inoc_geno geno_interaction          ANI_6K.trend     SE  df lower.CL upper.CL t.ratio p.value
# OA        between genospecies group      191.038 53.244 446   86.398  295.678  3.588  0.0004 ***
# OC        between genospecies group        1.763  0.425 446    0.928    2.599  4.146  <.0001 ***
# OE        between genospecies group      -31.103 53.877 446 -136.988   74.782 -0.577  0.5640 
# CC        between genospecies group       -0.514  0.439 446   -1.376    0.347 -1.173  0.2414 
# OA        within genospecies group         3.068  3.178 446   -3.177    9.313  0.966  0.3348 
# OC        within genospecies group        -1.823  2.204 446   -6.155    2.510 -0.827  0.4087 
# OE        within genospecies group        -0.842  6.294 446  -13.211   11.527 -0.134  0.8936 
# CC        within genospecies group         0.352  6.543 446  -12.508   13.211  0.054  0.9571 
# 
# Confidence level used: 0.95 
# 
# $contrasts
# contrast                                                    estimate     SE  df lower.CL upper.CL t.ratio p.value
# OA between genospecies group - OC between genospecies group  189.275 53.246 446   27.121   351.43  3.555  0.0099 *
# OA between genospecies group - OE between genospecies group  222.141 75.748 446   -8.540   452.82  2.933  0.0686 .
# OA between genospecies group - CC between genospecies group  191.552 53.246 446   29.399   353.71  3.598  0.0085 *
# OA between genospecies group - OA within genospecies group   187.970 53.339 446   25.533   350.41  3.524  0.0110 *
# OA between genospecies group - OC within genospecies group   192.861 53.290 446   30.574   355.15  3.619  0.0079 *
# OA between genospecies group - OE within genospecies group   191.880 53.615 446   28.603   355.16  3.579  0.0091 *
# OA between genospecies group - CC within genospecies group   190.686 53.645 446   27.318   354.05  3.555  0.0099 *
# OC between genospecies group - OE between genospecies group   32.866 53.879 446 -131.216   196.95  0.610  0.9987 
# OC between genospecies group - CC between genospecies group    2.278  0.611 446    0.417     4.14  3.729  0.0053 *
# OC between genospecies group - OA within genospecies group    -1.305  3.206 446  -11.069     8.46 -0.407  0.9999 
# OC between genospecies group - OC within genospecies group     3.586  2.245 446   -3.251    10.42  1.597  0.7519 
# OC between genospecies group - OE within genospecies group     2.605  6.308 446  -16.606    21.82  0.413  0.9999 
# OC between genospecies group - CC within genospecies group     1.411  6.557 446  -18.557    21.38  0.215  1.0000 
# OE between genospecies group - CC between genospecies group  -30.588 53.879 446 -194.671   133.49 -0.568  0.9992 
# OE between genospecies group - OA within genospecies group   -34.171 53.971 446 -198.533   130.19 -0.633  0.9984 
# OE between genospecies group - OC within genospecies group   -29.280 53.922 446 -193.494   134.93 -0.543  0.9994 
# OE between genospecies group - OE within genospecies group   -30.261 54.244 446 -195.453   134.93 -0.558  0.9993 
# OE between genospecies group - CC within genospecies group   -31.454 54.273 446 -196.737   133.83 -0.580  0.9991 
# CC between genospecies group - OA within genospecies group    -3.583  3.208 446  -13.352     6.19 -1.117  0.9530 
# CC between genospecies group - OC within genospecies group     1.309  2.248 446   -5.537     8.15  0.582  0.9991 
# CC between genospecies group - OE within genospecies group     0.328  6.309 446  -18.886    19.54  0.052  1.0000 
# CC between genospecies group - CC within genospecies group    -0.866  6.558 446  -20.838    19.11 -0.132  1.0000 
# OA within genospecies group - OC within genospecies group      4.891  3.868 446   -6.887    16.67  1.265  0.9113 
# OA within genospecies group - OE within genospecies group      3.910  7.051 446  -17.561    25.38  0.555  0.9993 
# OA within genospecies group - CC within genospecies group      2.716  7.274 446  -19.436    24.87  0.373  1.0000 
# OC within genospecies group - OE within genospecies group     -0.981  6.669 446  -21.290    19.33 -0.147  1.0000 
# OC within genospecies group - CC within genospecies group     -2.175  6.905 446  -23.202    18.85 -0.315  1.0000 
# OE within genospecies group - CC within genospecies group     -1.194  9.079 446  -28.843    26.46 -0.131  1.0000 
# 
# Confidence level used: 0.95 
# Conf-level adjustment: tukey method for comparing a family of 8 estimates 
# P value adjustment: tukey method for comparing a family of 8 estimates 



