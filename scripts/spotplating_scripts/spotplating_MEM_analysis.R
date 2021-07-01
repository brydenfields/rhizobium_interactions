# Author: Bryden Fields
# Email: bryden.fields@york.ac.uk
# Last updated: 01/07/2021

# Initial set up ---------

#set the working directory to the directory this script is in.
#setwd("~/what/ever/folder/you/are/working/from") 

library(lme4)
library(effects)
library(lmerTest)
library(tidyverse)
library(plyr)
library(dplyr)
library(boot)
library(emmeans)

source('../summarySE.R')



# data formatting for MEM -------------------


# read in averaged pairwise interaction inhibition zone results
all_means <- read.csv("../../data/raw_data/spotplating_data/dataset_IZspot_diameterdiffs_averages.csv", 
                      row.names = 1)

# The model is the same pipeline as for the supernatant assay:
# fit <- lmer(comp_ind ~ inoc_geno * sup_geno + (1|inoculant) + (1|supernatant), data=all62_means)
# with the following alterations
# fit <- lmer(IZ_diff ~ spot_geno * soft_geno + (1|spot_strain) + (1|soft_strain), data=all_means)

#There are 2 sets of random variables (the strain ID of spot and strain ID of lawn) which are related but not nested - we can draw a fully factoral matrix! (crossed random effects). 

# this model fits inhibition zone size with spot and lawn (soft) genotypes as the fixed effects and the strain IDs as random effects. 


# relevel factors so OA is the intercept instead of CC
all_means$spot_geno <- factor(all_means$spot_geno, levels = c('OA', 'OC', 'OE', 'CC'))
all_means$soft_geno <- factor(all_means$soft_geno, levels = c('OA', 'OC', 'OE', 'CC'))



# MEM for spot plating: strain effect ----------------------


# for strain effect to assess the random effects do the same ml, remove the random effects and compare models.

# full model
fit <- lmer(IZ_diff ~ spot_geno * soft_geno + (1|spot_strain) + (1|soft_strain), REML = F, data=all_means)
summary(fit)

# Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: IZ_diff ~ spot_geno * soft_geno + (1 | spot_strain) + (1 | soft_strain)
#    Data: all_means
# 
#      AIC      BIC   logLik deviance df.resid 
#   2274.8   2357.6  -1118.4   2236.8      557 
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -7.7118 -0.1607 -0.0159  0.1220  3.7799 
# 
# Random effects:
#  Groups      Name        Variance Std.Dev.
#  spot_strain (Intercept) 0.1207   0.3474  
#  soft_strain (Intercept) 8.2418   2.8709  
#  Residual                2.2852   1.5117  
# Number of obs: 576, groups:  spot_strain, 24; soft_strain, 24
# 
# Fixed effects:
#                         Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)               3.4994     1.2072  26.3757   2.899  0.00745 ** 
# spot_genoOC               3.2443     0.3940 115.5787   8.234 3.12e-13 ***
# spot_genoOE               4.2026     0.4288 115.5787   9.800  < 2e-16 ***
# spot_genoCC               2.0357     0.4089 115.5787   4.979 2.26e-06 ***
# soft_genoOC              -3.4994     1.6337  25.6712  -2.142  0.04185 *  
# soft_genoOE              -3.4994     1.7781  25.6712  -1.968  0.05995 .  
# soft_genoCC              -3.4994     1.6954  25.6712  -2.064  0.04925 *  
# spot_genoOC:soft_genoOC  -3.2443     0.4679 528.9775  -6.934 1.20e-11 ***
# spot_genoOE:soft_genoOC  -4.2026     0.5093 528.9775  -8.252 1.25e-15 ***
# spot_genoCC:soft_genoOC  -2.0357     0.4856 528.9775  -4.192 3.24e-05 ***
# spot_genoOC:soft_genoOE  -2.9500     0.5093 528.9775  -5.793 1.19e-08 ***
# spot_genoOE:soft_genoOE  -3.9534     0.5543 528.9775  -7.132 3.26e-12 ***
# spot_genoCC:soft_genoOE  -1.6708     0.5285 528.9775  -3.161  0.00166 ** 
# spot_genoOC:soft_genoCC  -3.2443     0.4856 528.9775  -6.681 6.01e-11 ***
# spot_genoOE:soft_genoCC  -3.9154     0.5285 528.9775  -7.408 5.08e-13 ***
# spot_genoCC:soft_genoCC  -2.0357     0.5039 528.9775  -4.040 6.14e-05 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# removed inoculant random effect
fit_noinoc <- lmer(IZ_diff ~ spot_geno * soft_geno + (1|soft_strain), REML = F, data=all_means)

# removed supernatant random effect
fit_nosup <- lmer(IZ_diff ~ spot_geno * soft_geno + (1|spot_strain), REML = F, data=all_means)


# LR test comparing full and no inoculant random effect model
anova(fit_noinoc, fit)

# Data: all_means
# Models:
# fit_noinoc: IZ_diff ~ spot_geno * soft_geno + (1 | soft_strain)
# fit: IZ_diff ~ spot_geno * soft_geno + (1 | spot_strain) + (1 | soft_strain)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)   
# fit_noinoc   18 2282.4 2360.8 -1123.2   2246.4                        
# fit          19 2274.8 2357.6 -1118.4   2236.8 9.5939  1   0.001952 **
#   ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# LR test comparing full and no supernatant random effect model
anova(fit_nosup, fit)

# Data: all_means
# Models:
# fit_nosup: IZ_diff ~ spot_geno * soft_geno + (1 | spot_strain)
# fit: IZ_diff ~ spot_geno * soft_geno + (1 | spot_strain) + (1 | soft_strain)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
# fit_nosup   18 3032.5 3110.9 -1498.3   2996.5                         
# fit         19 2274.8 2357.6 -1118.4   2236.8 759.68  1  < 2.2e-16 ***
#   ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# MEM for spot plating: genospecies effect ----------------------


# rerun full model
fit2 <- lmer(IZ_diff ~ spot_geno * soft_geno + (1|spot_strain) + (1|soft_strain), REML = F, data=all_means)
plot(allEffects(fit2))
summary(fit2)

# Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: IZ_diff ~ spot_geno * soft_geno + (1 | spot_strain) + (1 | soft_strain)
#    Data: all_means
# 
#      AIC      BIC   logLik deviance df.resid 
#   2274.8   2357.6  -1118.4   2236.8      557 
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -7.7118 -0.1607 -0.0159  0.1220  3.7799 
# 
# Random effects:
#  Groups      Name        Variance Std.Dev.
#  spot_strain (Intercept) 0.1207   0.3474  
#  soft_strain (Intercept) 8.2418   2.8709  
#  Residual                2.2852   1.5117  
# Number of obs: 576, groups:  spot_strain, 24; soft_strain, 24
# 
# Fixed effects:
#                         Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)               3.4994     1.2072  26.3757   2.899  0.00745 ** 
# spot_genoOC               3.2443     0.3940 115.5787   8.234 3.12e-13 ***
# spot_genoOE               4.2026     0.4288 115.5787   9.800  < 2e-16 ***
# spot_genoCC               2.0357     0.4089 115.5787   4.979 2.26e-06 ***
# soft_genoOC              -3.4994     1.6337  25.6712  -2.142  0.04185 *  
# soft_genoOE              -3.4994     1.7781  25.6712  -1.968  0.05995 .  
# soft_genoCC              -3.4994     1.6954  25.6712  -2.064  0.04925 *  
# spot_genoOC:soft_genoOC  -3.2443     0.4679 528.9775  -6.934 1.20e-11 ***
# spot_genoOE:soft_genoOC  -4.2026     0.5093 528.9775  -8.252 1.25e-15 ***
# spot_genoCC:soft_genoOC  -2.0357     0.4856 528.9775  -4.192 3.24e-05 ***
# spot_genoOC:soft_genoOE  -2.9500     0.5093 528.9775  -5.793 1.19e-08 ***
# spot_genoOE:soft_genoOE  -3.9534     0.5543 528.9775  -7.132 3.26e-12 ***
# spot_genoCC:soft_genoOE  -1.6708     0.5285 528.9775  -3.161  0.00166 ** 
# spot_genoOC:soft_genoCC  -3.2443     0.4856 528.9775  -6.681 6.01e-11 ***
# spot_genoOE:soft_genoCC  -3.9154     0.5285 528.9775  -7.408 5.08e-13 ***
# spot_genoCC:soft_genoCC  -2.0357     0.5039 528.9775  -4.040 6.14e-05 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# make the plot all effects as a ggplot
fit2_graph_data <- summarySE(data = all_means, measurevar = 'IZ_diff', groupvars = c('spot_geno', 'soft_geno'))

# genospecies palette
genoPalette <- c("#466EA9", "#4D9A7A", "#DA367E","#9ACD32")

# rename to include Soft or Inoc
all_means$spot_geno <- paste0('Spot-\n', all_means$spot_geno)
all_means$soft_geno <- paste0('Lawn-\n', all_means$soft_geno)

fit2_graph_data$spot_geno <- paste0('Spot-\n', fit2_graph_data$spot_geno)
fit2_graph_data$soft_geno <- paste0('Lawn-\n', fit2_graph_data$soft_geno)

fit2_graph_data$spot_geno <- factor(fit2_graph_data$spot_geno, levels = c('Spot-\nOA', 'Spot-\nOC', 'Spot-\nOE', 'Spot-\nCC'))
fit2_graph_data$soft_geno <- factor(fit2_graph_data$soft_geno, levels = c('Lawn-\nOA', 'Lawn-\nOC', 'Lawn-\nOE', 'Lawn-\nCC'))

all_means$spot_geno <- factor(all_means$spot_geno, levels = c('Spot-\nOA', 'Spot-\nOC', 'Spot-\nOE', 'Spot-\nCC'))
all_means$soft_geno <- factor(all_means$soft_geno, levels = c('Lawn-\nOA', 'Lawn-\nOC', 'Lawn-\nOE', 'Lawn-\nCC'))


(fit2_graph <- ggplot(data = fit2_graph_data, aes(x = spot_geno, y = IZ_diff,
                                                  fill = spot_geno, col = spot_geno)) +
    geom_jitter(data = all_means, alpha = 0.2, col = "black", width = 0.1) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = IZ_diff - ci, ymax = IZ_diff + ci), width=0.6, size=1.5) +
    facet_wrap(. ~ soft_geno, nrow = 1) +
    labs(x = 'Spotted genospecies group',
         y = 'Mean inhibition zone diameter (mm)') +
    scale_y_continuous(breaks = seq(0, 21, by = 5)) +
    scale_colour_manual(values = genoPalette) +
    scale_fill_manual(values = genoPalette) +
    theme_bw()+
    theme(legend.position = "none",
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          strip.text = element_text(size = 16),
          plot.title = element_text(size = 16, hjust = 0.5))
)
ggsave(plot = fit2_graph, '../../data/intermediate_data/spotplating_data/fullmodel_MEM_spotplating_ci_spotlawnfacet.pdf', width = 25, height = 15, units = 'cm')




# now run without interaction and do LR test
fit2_noint <- lmer(IZ_diff ~ spot_geno + soft_geno + (1|spot_strain) + (1|soft_strain), REML = F, data=all_means)
plot(allEffects(fit2_noint))
summary(fit2_noint)


# fit2_noint alleffects as a ggplot
fit2_graph_data_inoc <- summarySE(data = all_means, measurevar = 'IZ_diff', groupvars = c('spot_geno'))
fit2_graph_data_sup <- summarySE(data = all_means, measurevar = 'IZ_diff', groupvars = c('soft_geno'))

(fit2_graph_inoc <- ggplot(data = fit2_graph_data_inoc, aes(x = spot_geno, y = IZ_diff,
                                                            col = spot_geno, fill = spot_geno)) +
    geom_jitter(data = all_means, alpha = 0.2, col = "black", width = 0.1) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = IZ_diff - ci, ymax = IZ_diff + ci), width=0.6, size=1.5) +
    labs(x = 'Spotted genospecies group',
         y = 'Mean inhibition zone diameter (mm)') +
    scale_y_continuous(breaks = seq(0, 21, by = 5)) +
    scale_colour_manual(values = genoPalette) +
    scale_fill_manual(values = genoPalette) +
    theme_bw()+
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          strip.text = element_text(size = 16),
          legend.position = "none")
)
ggsave(plot = fit2_graph_inoc, '../../data/intermediate_data/spotplating_data/fullmodelnoint_MEM_spotplating_ci_spot.pdf', width = 12, height = 15, units = 'cm')

(fit2_graph_sup <- ggplot(data = fit2_graph_data_sup, aes(x = soft_geno, y = IZ_diff,
                                                          col = soft_geno, fill = soft_geno)) +
    geom_jitter(data = all_means, alpha = 0.2, col = "black", width = 0.1) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = IZ_diff - ci, ymax = IZ_diff + ci), width=0.6, size=1.5) +
    labs(x = 'Soft agar lawn genospecies group',
         y = 'Mean inhibition zone diameter (mm)') +
    scale_y_continuous(breaks = seq(0, 21, by = 5)) +
    scale_colour_manual(values = genoPalette) +
    scale_fill_manual(values = genoPalette) +
    theme_bw()+
    theme(axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          strip.text = element_text(size = 16),
          legend.position = "none")
)
ggsave(plot = fit2_graph_sup, '../../data/intermediate_data/spotplating_data/fullmodelnoint_MEM_spotplating_ci_soft.pdf', width = 12, height = 15, units = 'cm')


# comparison of models with and without interaction (likelihood ratio tests via anova of the full model with the effect in question against the model without the effect in question)
anova(fit2_noint, fit2)

# refitting model(s) with ML (instead of REML)
# Data: all_means
# Models:
#   fit2_noint: IZ_diff ~ spot_geno + soft_geno + (1 | spot_strain) + (1 | soft_strain)
# fit2: IZ_diff ~ spot_geno * soft_geno + (1 | spot_strain) + (1 | soft_strain)
#           Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# fit2_noint 10 2352.8 2396.3 -1166.4   2332.8                             
# fit2       19 2274.8 2357.6 -1118.4   2236.8 95.933      9  < 2.2e-16 ***

  

# bootstrapping of model comparison with interaction
# https://jofrhwld.github.io/teaching/courses/2017_lsa/lectures/Session_8.nb.html
# The function you want to apply to each re-fit model. This will almost always be just fixef
# nsim = The number of bootstrap simulations you want to do
fit2_bootstraps <- bootMer(fit2, 
                          FUN = fixef, 
                          nsim = 1000,
                          verbose = T)

# In the bootstrap summary, the original column corresponds to the original coefficient estimate from the model.
fit2_bootstraps

# PARAMETRIC BOOTSTRAP
# Call:bootMer(x = fit2, FUN = fixef, nsim = 1000, verbose = T)
# Bootstrap Statistics :
#   original       bias    std. error
# t1*   3.499417 -0.014229620   1.1567458
# t2*   3.244258  0.001411378   0.3804049
# t3*   4.202583 -0.006498104   0.4445434
# t4*   2.035685 -0.004532050   0.4098585
# t5*  -3.499417  0.045358717   1.6232129
# t6*  -3.499417  0.015176738   1.6701986
# t7*  -3.499417  0.021828128   1.6226244
# t8*  -3.244258 -0.016376631   0.4586267
# t9*  -4.202583  0.001433537   0.5260449
# t10* -2.035685  0.007291640   0.4927942
# t11* -2.950020  0.010069478   0.4945079
# t12* -3.953370  0.020309408   0.5691480
# t13* -1.670752  0.033324578   0.5369333
# t14* -3.244258  0.011330616   0.4602029
# t15* -3.915361  0.026262104   0.5240933
# t16* -2.035685  0.009934989   0.5033019
# 
# 24 message(s): boundary (singular) fit: see ?isSingular
# 52 warning(s): Model failed to converge with max|grad| = 0.00203828 (tol = 0.002, component 1) (and others)


# To get a confidence interval for these parameters, we need to use the boot.ci() function.
# If the confidence interval excludes 0 we can call this a reliable effect.
# these are the order of the parameters with their beta estimates.
fixef(fit2)
# confidence intervals I am particularly interested in is the spot_genoOA:soft_genoOE interaction
# and inocOE which were particularly good at growing as inoculants in all treatments spot_genoOE
# they were also seen as significant in the initial interaction
# I've also kept any other parameters that are significant 
boot.ci(fit2_bootstraps, index = 1, type = "perc") # intercept - inoc_OA: sup_OA
boot.ci(fit2_bootstraps, index = 3, type = "perc") # spot_genoOE
# BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
# Based on 1000 bootstrap replicates
# CALL :  boot.ci(boot.out = fit2_bootstraps, type = "perc", index = 3)
# Intervals : 
#   Level     Percentile     
# 95%   ( 3.355,  5.062 )  
# Calculations and Intervals on Original Scale



# rerun MEM excluding extreme strains -----------------------


# remove the extreme OA strains 144A, 154C, 145B
all_means2 <- all_means[all_means$spot_strain != "144A" & all_means$soft_strain != "144A" ,]
all_means3 <- all_means2[all_means2$spot_strain != "154C" & all_means2$soft_strain != "154C" ,]
all_means4 <- all_means3[all_means3$spot_strain != "145B" & all_means3$soft_strain != "145B" ,]


# fit without strains 144A, 154C, 145B
fit2_remstrains <- lmer(IZ_diff ~ spot_geno * soft_geno + (1|spot_strain) + (1|soft_strain), REML = F, data=all_means4)
summary(fit2_remstrains)

# Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: IZ_diff ~ spot_geno * soft_geno + (1 | spot_strain) + (1 | soft_strain)
#    Data: all_means4
# 
#      AIC      BIC   logLik deviance df.resid 
#    483.5    561.2   -222.8    445.5      422 
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -4.5293 -0.1878  0.0000  0.1658  9.8538 
# 
# Random effects:
#  Groups      Name        Variance Std.Dev.
#  spot_strain (Intercept) 0.00000  0.0000  
#  soft_strain (Intercept) 0.07447  0.2729  
#  Residual                0.14289  0.3780  
# Number of obs: 441, groups:  spot_strain, 21; soft_strain, 21
# 
# Fixed effects:
#                           Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)             -2.741e-14  2.017e-01  4.680e+01   0.000   1.0000    
# spot_genoOC              1.310e-14  1.506e-01  4.200e+02   0.000   1.0000    
# spot_genoOE              1.501e+00  1.594e-01  4.200e+02   9.420  < 2e-16 ***
# spot_genoCC              1.273e-14  1.543e-01  4.200e+02   0.000   1.0000    
# soft_genoOC              2.754e-14  2.411e-01  4.680e+01   0.000   1.0000    
# soft_genoOE              2.746e-14  2.552e-01  4.680e+01   0.000   1.0000    
# soft_genoCC              2.700e-14  2.471e-01  4.680e+01   0.000   1.0000    
# spot_genoOC:soft_genoOC -1.337e-14  1.800e-01  4.200e+02   0.000   1.0000    
# spot_genoOE:soft_genoOC -1.501e+00  1.905e-01  4.200e+02  -7.881 2.82e-14 ***
# spot_genoCC:soft_genoOC -1.317e-14  1.844e-01  4.200e+02   0.000   1.0000    
# spot_genoOC:soft_genoOE  2.942e-01  1.905e-01  4.200e+02   1.545   0.1232    
# spot_genoOE:soft_genoOE -1.252e+00  2.016e-01  4.200e+02  -6.211 1.27e-09 ***
# spot_genoCC:soft_genoOE  3.649e-01  1.952e-01  4.200e+02   1.870   0.0622 .  
# spot_genoOC:soft_genoCC -1.272e-14  1.844e-01  4.200e+02   0.000   1.0000    
# spot_genoOE:soft_genoCC -1.214e+00  1.952e-01  4.200e+02  -6.220 1.20e-09 ***
# spot_genoCC:soft_genoCC -1.243e-14  1.890e-01  4.200e+02   0.000   1.0000    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



# fit without strains and interaction
fit2_remstrain_noint <- lmer(IZ_diff ~ spot_geno + soft_geno + (1|spot_strain) + (1|soft_strain), REML = F, data=all_means4)

# fit without strains and soft_geno parameter
fit2_remstrains_nosupgeno <- lmer(IZ_diff ~ spot_geno + (1|spot_strain) + (1|soft_strain), REML = F, data=all_means4)

# LR test to compare if new models are significantly different from each other
anova(fit2_remstrains_nosupgeno, fit2_remstrain_noint, fit2_remstrains)


# bootstrapping of model comparison with interaction but without 3 strains strains
# https://jofrhwld.github.io/teaching/courses/2017_lsa/lectures/Session_8.nb.html
# The function you want to apply to each re-fit model. This will almost always be just fixef
# nsim = The number of bootstrap simulations you want to do
fit2_remstrains_bootstraps <- bootMer(fit2_remstrains, 
                                     FUN = fixef, 
                                     nsim = 1000,
                                     verbose = T)

# In the bootstrap summary, the original column corresponds to the original coefficient estimate from the model.
fixef(fit2_remstrains)
fit2_remstrains_bootstraps

# PARAMETRIC BOOTSTRAP
# 
# 
# Call:
#   bootMer(x = fit2_remstrains, FUN = fixef, nsim = 1000, verbose = T)
# 
# 
# Bootstrap Statistics :
#   original        bias    std. error
# t1*  -2.740775e-14 -0.0045082818   0.1936732
# t2*   1.309906e-14  0.0044858000   0.1515905
# t3*   1.501333e+00  0.0046241766   0.1553949
# t4*   1.272842e-14  0.0058127994   0.1549163
# t5*   2.753896e-14  0.0117738320   0.2283037
# t6*   2.745655e-14  0.0045563732   0.2488514
# t7*   2.699863e-14  0.0015714290   0.2386301
# t8*  -1.337026e-14 -0.0059780560   0.1802470
# t9*  -1.501333e+00 -0.0076847258   0.1780591
# t10* -1.317346e-14 -0.0079753518   0.1833082
# t11*  2.942381e-01  0.0027316598   0.1966744
# t12* -1.252120e+00  0.0013393128   0.2030266
# t13*  3.649333e-01 -0.0010121292   0.1933290
# t14* -1.271920e-14 -0.0032533239   0.1860973
# t15* -1.214111e+00 -0.0008824865   0.1914126
# t16* -1.243450e-14 -0.0052132781   0.1912803
# 
# 698 message(s): boundary (singular) fit: see ?isSingular
# 1 warning(s): Model failed to converge with max|grad| = 0.00252574 (tol = 0.002, component 1)


# To get a confidence interval for these parameters, we need to use the boot.ci() function.
# If the confidence interval excludes 0 we can call this a reliable effect.
# confidence intervals I am particularly interested in is the spot_genoOA:soft_genoOE interaction
# and inocOE which were particularly good at growing as inoculants in all treatments spot_genoOE
# they were also seen as significant in the initial interaction
# I've also kept any other parameters that are significant (spoiler only 2 are significant and one is the intercept)
fixef(fit2_remstrains)
boot.ci(fit2_remstrains_bootstraps, index = 1, type = "perc") # Intercept
boot.ci(fit2_remstrains_bootstraps, index = 3, type = "perc") # spot_genoOE 
# BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
# Based on 1000 bootstrap replicates
# 
# CALL : 
#   boot.ci(boot.out = fit2_remstrains_bootstraps, type = "perc", 
#           index = 3)
# 
# Intervals : 
#   Level     Percentile     
# 95%   ( 1.196,  1.820 )  
# Calculations and Intervals on Original Scale



# emmeans for OC / CC comparison on OA-reflevel model with all strains included ----------


# using full model (fit2)
# comparing pairwise the genospecies groups spotted inoculant on different agar lawns

emmip(fit2,  ~ spot_geno | soft_geno)

fit2_emmeans_spotgeno <- emmeans(fit2, specs = pairwise ~ spot_geno | soft_geno)
fit2_emmeans_spotgeno$emmeans
fit2_emmeans_spotgeno$contrasts

fit2_emmeans_spotgeno$contrasts %>%
  summary(infer = TRUE)


# soft_geno = Soft-OA:
#   contrast    estimate    SE  df lower.CL upper.CL t.ratio p.value
# I-OA - I-OC  -3.2443 0.409 132   -4.307   -2.181 -7.942  <.0001 
# I-OA - I-OE  -4.2026 0.445 132   -5.360   -3.046 -9.452  <.0001 
# I-OA - I-CC  -2.0357 0.424 132   -3.139   -0.933 -4.802  <.0001 
# I-OC - I-OE  -0.9583 0.430 132   -2.077    0.160 -2.229  0.1208 
# I-OC - I-CC   1.2086 0.409 132    0.146    2.272  2.958  0.0190 *
# I-OE - I-CC   2.1669 0.445 132    1.010    3.324  4.874  <.0001 
# 
# soft_geno = Soft-OC:
#   contrast    estimate    SE  df lower.CL upper.CL t.ratio p.value
# I-OA - I-OC   0.0000 0.387 108   -1.010    1.010  0.000  1.0000 
# I-OA - I-OE   0.0000 0.421 108   -1.099    1.099  0.000  1.0000 
# I-OA - I-CC   0.0000 0.402 108   -1.048    1.048  0.000  1.0000 
# I-OC - I-OE   0.0000 0.407 108   -1.063    1.063  0.000  1.0000 
# I-OC - I-CC   0.0000 0.387 108   -1.010    1.010  0.000  1.0000 
# I-OE - I-CC   0.0000 0.421 108   -1.099    1.099  0.000  1.0000 
# 
# soft_geno = Soft-OE:
#   contrast    estimate    SE  df lower.CL upper.CL t.ratio p.value
# I-OA - I-OC  -0.2942 0.437 167   -1.428    0.840 -0.673  0.9070 
# I-OA - I-OE  -0.2492 0.476 167   -1.483    0.985 -0.524  0.9532 
# I-OA - I-CC  -0.3649 0.453 167   -1.541    0.812 -0.805  0.8520 
# I-OC - I-OE   0.0450 0.460 167   -1.148    1.238  0.098  0.9997 
# I-OC - I-CC  -0.0707 0.437 167   -1.204    1.063 -0.162  0.9985 
# I-OE - I-CC  -0.1157 0.476 167   -1.350    1.118 -0.243  0.9949 
# 
# soft_geno = Soft-CC:
#   contrast    estimate    SE  df lower.CL upper.CL t.ratio p.value
# I-OA - I-OC   0.0000 0.409 132   -1.063    1.063  0.000  1.0000 
# I-OA - I-OE  -0.2872 0.445 132   -1.444    0.870 -0.646  0.9168 
# I-OA - I-CC   0.0000 0.424 132   -1.103    1.103  0.000  1.0000 
# I-OC - I-OE  -0.2872 0.430 132   -1.406    0.832 -0.668  0.9089 
# I-OC - I-CC   0.0000 0.409 132   -1.063    1.063  0.000  1.0000 
# I-OE - I-CC   0.2872 0.445 132   -0.870    1.444  0.646  0.9168 
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# Conf-level adjustment: tukey method for comparing a family of 4 estimates 
# P value adjustment: tukey method for comparing a family of 4 estimates 



# comparing pairwise the genospecies groups agar lawns with different spot inoculants

emmip(fit2,  ~ soft_geno | spot_geno)

fit2_emmeans_softgeno <- emmeans(fit2, specs = pairwise ~ soft_geno | spot_geno)
fit2_emmeans_softgeno$emmeans
fit2_emmeans_softgeno$contrasts

fit2_emmeans_softgeno$contrasts %>%
  summary(infer = TRUE)

# spot_geno = I-OA:
#   contrast          estimate   SE   df lower.CL upper.CL t.ratio p.value
# Soft-OA - Soft-OC    3.499 1.78 30.8   -1.346     8.34  1.961  0.2247 
# Soft-OA - Soft-OE    3.499 1.94 30.8   -1.774     8.77  1.802  0.2919 
# Soft-OA - Soft-CC    3.499 1.85 30.8   -1.529     8.53  1.890  0.2533 
# Soft-OC - Soft-OE    0.000 1.88 30.8   -5.100     5.10  0.000  1.0000 
# Soft-OC - Soft-CC    0.000 1.78 30.8   -4.846     4.85  0.000  1.0000 
# Soft-OE - Soft-CC    0.000 1.94 30.8   -5.274     5.27  0.000  1.0000 
# 
# spot_geno = I-OC:
#   contrast          estimate   SE   df lower.CL upper.CL t.ratio p.value
# Soft-OA - Soft-OC    6.744 1.78 30.4    1.908    11.58  3.789  0.0036 *
# Soft-OA - Soft-OE    6.449 1.94 30.4    1.186    11.71  3.329  0.0116 *
# Soft-OA - Soft-CC    6.744 1.85 30.4    1.725    11.76  3.651  0.0051 *
# Soft-OC - Soft-OE   -0.294 1.87 30.4   -5.384     4.80 -0.157  0.9986 
# Soft-OC - Soft-CC    0.000 1.78 30.4   -4.836     4.84  0.000  1.0000 
# Soft-OE - Soft-CC    0.294 1.94 30.4   -4.969     5.56  0.152  0.9987 
# 
# spot_geno = I-OE:
#   contrast          estimate   SE   df lower.CL upper.CL t.ratio p.value
# Soft-OA - Soft-OC    7.702 1.79 31.3    2.843    12.56  4.300  0.0009 *
# Soft-OA - Soft-OE    7.453 1.95 31.3    2.164    12.74  3.823  0.0031 *
# Soft-OA - Soft-CC    7.415 1.86 31.3    2.372    12.46  3.989  0.0020 *
# Soft-OC - Soft-OE   -0.249 1.89 31.3   -5.363     4.86 -0.132  0.9992 
# Soft-OC - Soft-CC   -0.287 1.79 31.3   -5.146     4.57 -0.160  0.9985 
# Soft-OE - Soft-CC   -0.038 1.95 31.3   -5.327     5.25 -0.019  1.0000 
# 
# spot_geno = I-CC:
#   contrast          estimate   SE   df lower.CL upper.CL t.ratio p.value
# Soft-OA - Soft-OC    5.535 1.78 30.8    0.690    10.38  3.102  0.0203 *
# Soft-OA - Soft-OE    5.170 1.94 30.8   -0.104    10.44  2.662  0.0563 *
# Soft-OA - Soft-CC    5.535 1.85 30.8    0.507    10.56  2.989  0.0266 *
# Soft-OC - Soft-OE   -0.365 1.88 30.8   -5.465     4.73 -0.194  0.9973 
# Soft-OC - Soft-CC    0.000 1.78 30.8   -4.846     4.85  0.000  1.0000 
# Soft-OE - Soft-CC    0.365 1.94 30.8   -4.909     5.64  0.188  0.9976 
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# Conf-level adjustment: tukey method for comparing a family of 4 estimates 
# P value adjustment: tukey method for comparing a family of 4 estimates 



