# Author: Bryden Fields
# Email: bryden.fields@york.ac.uk
# Last updated: 02/02/2021

# Initial setup -------------------------------------------------------------------

setwd('/Users/brydenfields/Documents/Publications/2021_Rhizobiuminteractions_paper/scripts/supernatant_scripts')

library(tidyverse)
library(RColorBrewer)
library(plyr); library(dplyr)
library(ggrepel)

source('../summarySE.R')



# Descriptive stats -------------------------------------------------------------------


#Creating descriptive statistics summary of OD values for strains, genospecies and environment
#filter to get only 62 hr RGI data.
# data_62 <- data %>% 
#   filter(time %in% c('62')) %>%
#    filter(!supernatant %in% c('100%_TY'))
#    filter(!supernatant %in% c('50%_TY'))
# #group RGI (comp_ind) by supernatant (not averaging replicates of each inoc/sup condition first to maintain variation)
# compind_sup <- summarySE(data_62, measurevar = 'comp_ind', groupvars = c('supernatant', 'sup_geno'))
# #group RGI (comp_ind) by inoculant (not averaging replicates of each inoc/sup condition first to maintain variation)
# compind_inoc <- summarySE(data_62, measurevar = 'comp_ind', groupvars = c('inoculant', 'inoc_geno'))
# # save data. This will be combined and reuploaded as 'combined_compind' dataframe
# # write.csv(compind_sup, file = 'compind_sup.csv')
# # write.csv(compind_inoc, file = 'compind_inoc.csv')


# this was created previously with the above commented out code:
combined_compind <- read.csv('../../data/raw_data/supernatant_data/combined_compind_sup_inoc.csv')


combined_compind$strain <- factor(combined_compind$strain)
combined_compind$geno <- factor(combined_compind$geno)



# Supernatant/Inoculant strain effects correlation graph -------------------------------------------------------------------

# make all strain names have SM on the front
combined_compind$strain <- paste0("SM",combined_compind$strain)


#plot each strain as an inhibitor/facilitator based on inoculum and supernatant
(compind_facil_inhib <- ggplot(combined_compind, aes(x = comp_ind_sup, y = comp_ind_inoc, colour = geno), group = strain)+
      geom_point()+
      geom_smooth(method=lm, se = FALSE, colour = 'grey60')+
      geom_text_repel(aes(label = strain),
                      max.overlaps = Inf,
                      box.padding   = 0.4, 
                      point.padding = 0.6,
                      segment.color = 'grey40', 
                      show.legend = FALSE) +
      geom_errorbar(aes(ymin = comp_ind_inoc - ci_inoc, ymax = comp_ind_inoc + ci_inoc, width = 0))+
      geom_errorbarh(aes(xmin = comp_ind_sup - ci_sup, xmax = comp_ind_sup + ci_sup, width = 0))+
      labs(x = "Supportiveness as supernatant", 
           y = "Growth as inoculant") +
      xlim(0.5, 1.5)+
      #scale_x_reverse(limits = c(1.5, 0.5)) +
      ylim(0.5, 1.5)+
      scale_colour_manual(values=c(CC = "#9ACD32", OA = "#466EA9", OC = "#4D9A7A", OE ="#DA367E"))+
      geom_hline(yintercept = 1, colour = "grey50", size = 0.5, alpha = 0.5)+
      geom_vline(xintercept = 1, colour = "grey50", size = 0.5, alpha = 0.5)+
      #guides(colour = guide_legend(override.aes = list(size=2))) +
      theme_bw()+
      theme(legend.title = element_blank(),
            axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 16),
            legend.text = element_text(size = 16),
            axis.title = element_text(size = 20)))
#function to save specific ggplot graph.  additional parameters(width=,height=,units='cm')
ggsave('../../data/intermediate_data/supernatant_data/RGIinoc_RGIsup_correlation.pdf', plot = compind_facil_inhib, width = 19, height = 19, unit = 'cm')


# Pearson's correlation R statistic
cor.test(combined_compind$comp_ind_sup, combined_compind$comp_ind_inoc, method="pearson")

# Pearson's product-moment correlation
# 
# data:  combined_compind$comp_ind_sup and combined_compind$comp_ind_inoc
# t = -4.6114, df = 22, p-value = 0.0001355
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.8609787 -0.4150665
# sample estimates:
#        cor 
# -0.7010773 


# Simple linear regression
# linear regression looking at the RGI as inoculant by RGI as supernatant in the dataframe called combined_compind
mod <- lm(data = combined_compind, comp_ind_inoc ~ comp_ind_sup)
# summary shows the main results of the model including the coefficients, p-values, r2 values, f-statistic values etc.
summary(mod)

# Call:
#    lm(formula = comp_ind_inoc ~ comp_ind_sup, data = combined_compind)
# 
# Residuals:
#    Min        1Q    Median        3Q       Max 
# -0.193666 -0.068554 -0.005469  0.051245  0.303662 
# 
# Coefficients:
#    Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    1.8186     0.1812  10.036 1.13e-09 ***
#    comp_ind_sup  -0.8387     0.1819  -4.611 0.000136 ***
#    ---
#    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.107 on 22 degrees of freedom
# Multiple R-squared:  0.4915,	Adjusted R-squared:  0.4684 
# F-statistic: 21.27 on 1 and 22 DF,  p-value: 0.0001355

# testing for normality and homogeneity as done previously with the one way anova. 
plot(mod, which=1)
plot(mod, which=2)

shapiro.test(mod$residuals)
# Shapiro-Wilk normality test
# 
# data:  mod$residuals
# W = 0.95503, p-value = 0.3468




