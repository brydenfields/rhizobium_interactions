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

# Read in supernatant results but for AUC data - the growth_rate_comparisons_analysis.r script should be run first to generate this data in the 'intermediate_data' output directory.
all <- read.csv("../../data/intermediate_data/growth_metric_comparison_data/AUC_growth_data.csv")

# removing the TY controls and where inoculant strain is the same as the supernatant strain
all62 <- all[all$supernatant != "100_TY" & all$supernatant != "50_TY" & all$supernatant != all$inoculant,]

# averaging over the means
all62_means <- na.omit(as.data.frame.table(tapply(all62$AUC_index, list(all62$inoc_geno, all62$sup_geno, all62$supernatant, all62$inoculant), mean)))
names(all62_means) <- c("inoc_geno", "sup_geno", "supernatant", "inoculant", "AUC_index")

# relevel categorical levels MEM where OA is the intercept
all62_means$inoc_geno <- factor(all62_means$inoc_geno, levels = c('OA', 'OC', 'OE', 'CC'))
all62_means$sup_geno <- factor(all62_means$sup_geno, levels = c('OA', 'OC', 'OE', 'CC'))



# MEM for supernatant: strain effect ----------------------

# for strain effect to assess the random effects do the same ml, remove the random effects and compare models.

# full model
fit <- lmer(AUC_index ~ inoc_geno * sup_geno + (1|inoculant) + (1|supernatant), REML = F, data=all62_means)
summary(fit)

# removed inoculant random effect
fit_noinoc <- lmer(AUC_index ~ inoc_geno * sup_geno + (1|supernatant), REML = F, data=all62_means)

# removed supernatant random effect
fit_nosup <- lmer(AUC_index ~ inoc_geno * sup_geno + (1|inoculant), REML = F, data=all62_means)


# LR test comparing full and no inoculant random effect model
anova(fit_noinoc, fit)

# Data: all62_means
# Models:
#     fit_noinoc: comp_ind ~ inoc_geno * sup_geno + (1 | supernatant)
# fit: comp_ind ~ inoc_geno * sup_geno + (1 | inoculant) + (1 | supernatant)
#            npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)    
# fit_noinoc   18 -505.44 -427.79 270.72  -541.44                         
# fit          19 -858.25 -776.29 448.13  -896.25 354.82  1  < 2.2e-16 ***

# LR test comparing full and no supernatant random effect model
anova(fit_nosup, fit)

# Data: all62_means
# Models:
#     fit_nosup: comp_ind ~ inoc_geno * sup_geno + (1 | inoculant)
# fit: comp_ind ~ inoc_geno * sup_geno + (1 | inoculant) + (1 | supernatant)
# npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)    
# fit_nosup   18 -614.24 -536.60 325.12  -650.24                         
# fit         19 -858.25 -776.29 448.13  -896.25 246.01  1  < 2.2e-16 ***





# MEM for supernatant: genospecies effect-------------------

# run full model
#There are 2 sets of random variables (the strain IDs for the innocula and the supernatants) which are related but not nested - therefore we can draw a fully factoral matrix! (crossed random effects).
fit2 <- lmer(AUC_index ~ inoc_geno * sup_geno + (1|inoculant) + (1|supernatant), REML = F, data=all62_means)
plot(allEffects(fit2))
summary(fit2)

# Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: comp_ind ~ inoc_geno * sup_geno + (1 | inoculant) + (1 | supernatant)
#    Data: all62_means
# 
#      AIC      BIC   logLik deviance df.resid 
#   -858.3   -776.3    448.1   -896.3      533 
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -3.7393 -0.4899 -0.0425  0.4248  6.9848 
# 
# Random effects:
#  Groups      Name        Variance Std.Dev.
#  inoculant   (Intercept) 0.010893 0.10437 
#  supernatant (Intercept) 0.007039 0.08390 
#  Residual                0.008800 0.09381 
# Number of obs: 552, groups:  inoculant, 24; supernatant, 24
# 
# Fixed effects:
#                          Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)              0.968512   0.057290  50.202794  16.905  < 2e-16 ***
# inoc_genoOC              0.070689   0.062247  29.375171   1.136 0.265290    
# inoc_genoOE              0.208752   0.067683  29.260421   3.084 0.004426 ** 
# inoc_genoCC              0.074302   0.064567  29.322551   1.151 0.259127    
# sup_genoOC              -0.111758   0.051786  32.171426  -2.158 0.038488 *  
# sup_genoOE              -0.382209   0.056283  31.990536  -6.791 1.13e-07 ***
# sup_genoCC              -0.006766   0.053705  32.088460  -0.126 0.900528    
# inoc_genoOC:sup_genoOC   0.103648   0.030363 505.043471   3.414 0.000693 ***
# inoc_genoOE:sup_genoOC   0.039788   0.032367 505.043471   1.229 0.219540    
# inoc_genoCC:sup_genoOC   0.005379   0.030933 505.043471   0.174 0.862029    
# inoc_genoOC:sup_genoOE   0.247622   0.032367 505.043471   7.650 1.02e-13 ***
# inoc_genoOE:sup_genoOE   0.219678   0.036332 505.043471   6.046 2.88e-09 ***
# inoc_genoCC:sup_genoOE   0.253220   0.033533 505.043471   7.551 2.03e-13 ***
# inoc_genoOC:sup_genoCC   0.040557   0.030933 505.043471   1.311 0.190414    
# inoc_genoOE:sup_genoCC   0.013217   0.033533 505.043471   0.394 0.693638    
# inoc_genoCC:sup_genoCC  -0.003805   0.032796 505.043471  -0.116 0.907690    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# make the plot all effects as a ggplot
fit2_graph_data <- summarySE(data = all62_means, measurevar = 'AUC_index', groupvars = c('inoc_geno', 'sup_geno'))

# genospecies palette
genoPalette <- c("#466EA9", "#4D9A7A", "#DA367E","#9ACD32")

# add S- and I- to the genospecies names
all62_means$inoc_geno <- paste0('I-', all62_means$inoc_geno)
all62_means$sup_geno <- paste0('Sup-', all62_means$sup_geno)

fit2_graph_data$inoc_geno <- paste0('I-', fit2_graph_data$inoc_geno)
fit2_graph_data$sup_geno <- paste0('Sup-', fit2_graph_data$sup_geno)

fit2_graph_data$inoc_geno <- factor(fit2_graph_data$inoc_geno, levels = c('I-OA', 'I-OC', 'I-OE', 'I-CC'))
fit2_graph_data$sup_geno <- factor(fit2_graph_data$sup_geno, levels = c('Sup-OA', 'Sup-OC', 'Sup-OE', 'Sup-CC'))

all62_means$inoc_geno <- factor(all62_means$inoc_geno, levels = c('I-OA', 'I-OC', 'I-OE', 'I-CC'))
all62_means$sup_geno <- factor(all62_means$sup_geno, levels = c('Sup-OA', 'Sup-OC', 'Sup-OE', 'Sup-CC'))

(fit2_graph <- ggplot(data = fit2_graph_data, aes(x = inoc_geno, y = AUC_index,
                      fill = inoc_geno, col = inoc_geno)) +
    geom_jitter(data = all62_means, alpha = 0.2, col = "black", width = 0.1) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = AUC_index - ci, ymax = AUC_index + ci), width=0.6, size=1.5) +
    geom_hline(yintercept = 1) +
    facet_wrap(. ~ sup_geno, nrow = 1) +
    labs(x = 'Inoculant genospecies group',
         y = 'AUC Index') +
    scale_y_continuous(breaks = seq(0, 2.5, by = 0.5)) +
    scale_colour_manual(values = genoPalette) +
    scale_fill_manual(values = genoPalette) +
    theme_bw()+
    theme(legend.position = "none",
          axis.title = element_text(size = 16),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          strip.text = element_text(size = 16),
          plot.title = element_text(size = 16, hjust = 0.5))
)
ggsave(plot = fit2_graph, '../../data/intermediate_data/growth_metric_comparison_data/MEMsupernatant_ci_inocsupfacets_AUCindex.pdf', width = 25, height = 15, units = 'cm')



# now run without interaction and do LR test
fit2_noint <- lmer(AUC_index ~ inoc_geno + sup_geno + (1|inoculant) + (1|supernatant), REML = F, data=all62_means)
plot(allEffects(fit2_noint))
summary(fit2_noint)

# Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: comp_ind ~ inoc_geno + sup_geno + (1 | inoculant) + (1 | supernatant)
#    Data: all62_means
# 
#      AIC      BIC   logLik deviance df.resid 
#   -773.7   -730.6    396.9   -793.7      542 
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -3.7557 -0.4435  0.0064  0.4405  6.3025 
# 
# Random effects:
#  Groups      Name        Variance Std.Dev.
#  inoculant   (Intercept) 0.010811 0.1040  
#  supernatant (Intercept) 0.006956 0.0834  
#  Residual                0.010780 0.1038  
# Number of obs: 552, groups:  inoculant, 24; supernatant, 24
# 
# Fixed effects:
#             Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)  0.91196    0.05571 44.87296  16.371  < 2e-16 ***
# inoc_genoOC  0.16393    0.05909 23.85308   2.774 0.010576 *  
# inoc_genoOE  0.27044    0.06431 23.85308   4.205 0.000317 ***
# inoc_genoCC  0.13089    0.06132 23.85308   2.135 0.043278 *  
# sup_genoOC  -0.07061    0.04794 23.64710  -1.473 0.153960    
# sup_genoOE  -0.19990    0.05218 23.64710  -3.831 0.000823 ***
# sup_genoCC   0.01008    0.04975 23.64710   0.203 0.841114    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr) inc_OC inc_OE inc_CC sp_gOC sp_gOE
# inoc_genoOC -0.572                                   
# inoc_genoOE -0.526  0.495                            
# inoc_genoCC -0.551  0.519  0.477                     
# sup_genoOC  -0.465  0.002  0.001  0.001              
# sup_genoOE  -0.427  0.001  0.002  0.001  0.495       
# sup_genoCC  -0.448  0.001  0.001  0.002  0.519  0.477




# comparison of models with and without interaction (likelihood ratio tests via anova of the full model with the effect in question against the model without the effect in question)
anova(fit2_noint, fit2)

# refitting model(s) with ML (instead of REML)
# Data: all62_means
# Models:
#     fit2_noint: comp_ind ~ inoc_geno + sup_geno + (1 | inoculant) + (1 | supernatant)
# fit2: comp_ind ~ inoc_geno * sup_geno + (1 | inoculant) + (1 | supernatant)
# Df     AIC     BIC logLik deviance Chisq Chi Df Pr(>Chisq)    
# fit2_noint 10 -773.75 -730.61 396.87  -793.75                            
# fit2       19 -858.25 -776.29 448.13  -896.25 102.5      9  < 2.2e-16 ***

# remove the supernatant genotype parameter to see if it significantly affects the model fit in lr test.
fit2_nosupgeno <- lmer(AUC_index ~ inoc_geno + (1|inoculant) + (1|supernatant), REML = F, data=all62_means)
anova(fit2_nosupgeno, fit2_noint, fit2)

# refitting model(s) with ML (instead of REML)
# Data: all62_means
# Models:
#     fit2_nosupgeno: comp_ind ~ inoc_geno + (1 | inoculant) + (1 | supernatant)
# fit2_noint: comp_ind ~ inoc_geno + sup_geno + (1 | inoculant) + (1 | supernatant)
# fit2: comp_ind ~ inoc_geno * sup_geno + (1 | inoculant) + (1 | supernatant)
# Df     AIC     BIC logLik deviance   Chisq Chi Df Pr(>Chisq)    
# fit2_nosupgeno  7 -765.16 -734.96 389.58  -779.16                              
# fit2_noint     10 -773.75 -730.61 396.87  -793.75  14.591      3   0.002202 ** 
# fit2           19 -858.25 -776.29 448.13  -896.25 102.504      9  < 2.2e-16 ***


# bootstrapping of model comparison with interaction
# https://jofrhwld.github.io/teaching/courses/2017_lsa/lectures/Session_8.nb.html
# The function you want to apply to each re-fit model. This will almost always be just fixef
# nsim = The number of bootstrap simulations you want to do
fit2_bootstraps <- bootMer(fit2, 
                          FUN = fixef, 
                          nsim = 1000,
                          verbose = T)

# In the bootstrap summary, the original column corresponds to the original coefficient estimate from the model.
fixef(fit2)
fit2_bootstraps

# PARAMETRIC BOOTSTRAP
# Call:bootMer(x = fit2, FUN = fixef, nsim = 1000, verbose = T)
# Bootstrap Statistics :
#     original        bias    std. error
# t1*   0.968511559 -0.0015269315  0.05500357
# t2*   0.070689440  0.0030224823  0.06321701
# t3*   0.208752170 -0.0032234528  0.06676182
# t4*   0.074301727 -0.0019718697  0.06355803
# t5*  -0.111758261  0.0003665257  0.05299713
# t6*  -0.382209186 -0.0012399888  0.05466065
# t7*  -0.006766203  0.0010016949  0.05277809
# t8*   0.103648321 -0.0003208637  0.02986778
# t9*   0.039788422  0.0007326463  0.03154622
# t10*  0.005378629  0.0009654303  0.03153673
# t11*  0.247621794 -0.0008289471  0.03270996
# t12*  0.219678118  0.0008163374  0.03698619
# t13*  0.253220118  0.0009823197  0.03394518
# t14*  0.040556639  0.0010734713  0.03073006
# t15*  0.013216990  0.0001018318  0.03403038
# t16* -0.003804701  0.0013608433  0.03261010
# 
# 30 warning(s): Model failed to converge with max|grad| = 0.00201935 (tol = 0.002, component 1) (and others)

# To get a confidence interval for these parameters, we need to use the boot.ci() function.
# If the confidence interval excludes 0 we can call this a reliable effect.
# these are the order of the parameters with their beta estimates.
fixef(fit2)

# confidence intervals I am particularly interested in are:
# inoc_genoOA:sup_genoOE interaction
# and inocOE which were particularly good at growing as inoculants in all treatments inoc_genoOE
# they were also seen as significant in the initial interaction
# I've also kept any other parameters that are significant (spoiler: only 2 are significant and one is the intercept)
boot.ci(fit2_bootstraps, index = 1, type = "perc") # intercept - inoc_OA: sup_OA

boot.ci(fit2_bootstraps, index = 3, type = "perc") # inoc_genoOE

# BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
# Based on 1000 bootstrap replicates
# CALL :  boot.ci(boot.out = fit2_bootstraps, type = "perc", index = 3)
# Intervals : 
#     Level     Percentile     
# 95%   ( 0.0721,  0.3275 )  
# Calculations and Intervals on Original Scale


boot.ci(fit2_bootstraps, index = 6, type = "perc") # sup_genoOE

# BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
# Based on 1000 bootstrap replicates
# CALL : boot.ci(boot.out = fit2_bootstraps, type = "perc", index = 6)
# Intervals : 
#     Level     Percentile     
# 95%   (-0.4915, -0.2797 )  
# Calculations and Intervals on Original Scale




# emmeans for OC / CC comparison on OA-reflevel model with all strains included ----------

# using full model (fit2)
# comparing pairwise the genospecies groups inoculant in different supernatants

emmip(fit2,  ~ inoc_geno | sup_geno)

fit2_emmeans_inocgeno <- emmeans(fit2, specs = pairwise ~ inoc_geno | sup_geno)
fit2_emmeans_inocgeno$emmeans
fit2_emmeans_inocgeno$contrasts

fit2_emmeans_inocgeno$contrasts %>%
    summary(infer = TRUE)


# sup_geno = Sup-OA:
#     contrast    estimate     SE   df lower.CL upper.CL t.ratio p.value
# I-OA - I-OC -0.07069 0.0672 34.3  -0.2520  0.11064 -1.052  0.7202
# I-OA - I-OE -0.20875 0.0730 34.2  -0.4060 -0.01153 -2.858  0.0347 *
# I-OA - I-CC -0.07430 0.0697 34.3  -0.2624  0.11381 -1.066  0.7119
# I-OC - I-OE -0.13806 0.0703 33.5  -0.3281  0.05196 -1.964  0.2219
# I-OC - I-CC -0.00361 0.0668 33.5  -0.1842  0.17694 -0.054  0.9999
# I-OE - I-CC  0.13445 0.0727 33.5  -0.0621  0.33097  1.849  0.2690
# 
# sup_geno = Sup-OC:
#     contrast    estimate     SE   df lower.CL upper.CL t.ratio p.value
# I-OA - I-OC -0.17434 0.0665 32.9  -0.3543  0.00567 -2.620  0.0605
# I-OA - I-OE -0.24854 0.0722 32.4  -0.4439 -0.05314 -3.444  0.0083 **
# I-OA - I-CC -0.07968 0.0688 32.4  -0.2660  0.10663 -1.158  0.6570
# I-OC - I-OE -0.07420 0.0700 32.8  -0.2636  0.11520 -1.060  0.7157
# I-OC - I-CC  0.09466 0.0665 32.9  -0.0853  0.27466  1.423  0.4946
# I-OE - I-CC  0.16886 0.0722 32.4  -0.0265  0.36426  2.340  0.1099
# 
# sup_geno = Sup-OE:
#     contrast    estimate     SE   df lower.CL upper.CL t.ratio p.value
# I-OA - I-OC -0.31831 0.0675 35.1  -0.5003 -0.13632 -4.717  0.0002 ***
# I-OA - I-OE -0.42843 0.0741 36.4  -0.6278 -0.22907 -5.785  <.0001 ***
# I-OA - I-CC -0.32752 0.0700 35.1  -0.5164 -0.13866 -4.677  0.0002 ***
# I-OC - I-OE -0.11012 0.0717 36.5  -0.3030  0.08274 -1.537  0.4267
# I-OC - I-CC -0.00921 0.0675 35.1  -0.1912  0.17278 -0.136  0.9991
# I-OE - I-CC  0.10091 0.0741 36.4  -0.0984  0.30026  1.363  0.5303
# 
# sup_geno = Sup-CC:
#     contrast    estimate     SE   df lower.CL upper.CL t.ratio p.value
# I-OA - I-OC -0.11125 0.0668 33.5  -0.2918  0.06931 -1.665  0.3574
# I-OA - I-OE -0.22197 0.0727 33.5  -0.4185 -0.02545 -3.053  0.0219 *
# I-OA - I-CC -0.07050 0.0697 34.3  -0.2586  0.11762 -1.012  0.7437
# I-OC - I-OE -0.11072 0.0703 33.5  -0.3008  0.07930 -1.575  0.4062
# I-OC - I-CC  0.04075 0.0672 34.3  -0.1406  0.22207  0.607  0.9293
# I-OE - I-CC  0.15147 0.0730 34.2  -0.0458  0.34870  2.074  0.1821
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# Conf-level adjustment: tukey method for comparing a family of 4 estimates 
# P value adjustment: tukey method for comparing a family of 4 estimates 



# comparing pairwise the genospecies groups supernatants with different inoculants

emmip(fit2,  ~ sup_geno | inoc_geno)

fit2_emmeans_supgeno <- emmeans(fit2, specs = pairwise ~ sup_geno | inoc_geno)
fit2_emmeans_supgeno$emmeans
fit2_emmeans_supgeno$contrasts

fit2_emmeans_supgeno$contrasts %>%
    summary(infer = TRUE)


# inoc_geno = I-OA:
#     contrast        estimate     SE   df lower.CL upper.CL t.ratio p.value
# Sup-OA - Sup-OC  0.11176 0.0555 37.2 -0.03745   0.2610  2.014  0.2012 
# Sup-OA - Sup-OE  0.38221 0.0603 37.0  0.21997   0.5444  6.336  <.0001 ***
# Sup-OA - Sup-CC  0.00677 0.0576 37.1 -0.14801   0.1615  0.118  0.9994 
# Sup-OC - Sup-OE  0.27045 0.0579 35.9  0.11442   0.4265  4.669  0.0002 ***
# Sup-OC - Sup-CC -0.10499 0.0550 35.9 -0.25325   0.0433 -1.908  0.2430 
# Sup-OE - Sup-CC -0.37544 0.0599 35.9 -0.53680  -0.2141 -6.267  <.0001 ***
# 
# inoc_geno = I-OC:
#     contrast        estimate     SE   df lower.CL upper.CL t.ratio p.value
# Sup-OA - Sup-OC  0.00811 0.0547 35.0 -0.13946   0.1557  0.148  0.9988 
# Sup-OA - Sup-OE  0.13459 0.0592 34.2 -0.02539   0.2946  2.272  0.1248 
# Sup-OA - Sup-CC -0.03379 0.0565 34.2 -0.18632   0.1187 -0.598  0.9319 
# Sup-OC - Sup-OE  0.12648 0.0576 34.9 -0.02878   0.2817  2.197  0.1439 
# Sup-OC - Sup-CC -0.04190 0.0547 35.0 -0.18947   0.1057 -0.766  0.8693 
# Sup-OE - Sup-CC -0.16838 0.0592 34.2 -0.32836  -0.0084 -2.842  0.0360 *
# 
# inoc_geno = I-OE:
#     contrast        estimate     SE   df lower.CL upper.CL t.ratio p.value
# Sup-OA - Sup-OC  0.07197 0.0559 38.4 -0.07806   0.2220  1.288  0.5759 
# Sup-OA - Sup-OE  0.16253 0.0615 40.5 -0.00235   0.3274  2.641  0.0546 
# Sup-OA - Sup-CC -0.00645 0.0580 38.4 -0.16215   0.1492 -0.111  0.9995 
# Sup-OC - Sup-OE  0.09056 0.0596 40.6 -0.06898   0.2501  1.520  0.4350 
# Sup-OC - Sup-CC -0.07842 0.0559 38.4 -0.22845   0.0716 -1.404  0.5049 
# Sup-OE - Sup-CC -0.16898 0.0615 40.5 -0.33386  -0.0041 -2.746  0.0427 
# 
# inoc_geno = I-CC:
#     contrast        estimate     SE   df lower.CL upper.CL t.ratio p.value
# Sup-OA - Sup-OC  0.10638 0.0550 35.9 -0.04187   0.2546  1.933  0.2327 
# Sup-OA - Sup-OE  0.12899 0.0599 35.9 -0.03237   0.2903  2.153  0.1561 
# Sup-OA - Sup-CC  0.01057 0.0576 37.1 -0.14420   0.1653  0.184  0.9978 
# Sup-OC - Sup-OE  0.02261 0.0579 35.9 -0.13342   0.1786  0.390  0.9795 
# Sup-OC - Sup-CC -0.09581 0.0555 37.2 -0.24502   0.0534 -1.727  0.3247 
# Sup-OE - Sup-CC -0.11842 0.0603 37.0 -0.28066   0.0438 -1.963  0.2203 
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# Conf-level adjustment: tukey method for comparing a family of 4 estimates 
# P value adjustment: tukey method for comparing a family of 4 estimates 




