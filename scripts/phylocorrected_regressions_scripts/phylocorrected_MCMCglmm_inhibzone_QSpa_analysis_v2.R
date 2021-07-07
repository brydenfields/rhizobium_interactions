# Author: Bryden Fields
# Email: bryden.fields@york.ac.uk
# Last updated: 01/07/2021

# Initial set-up -----------

library(tidyverse)
library(plyr); library(dplyr)
library(reshape2)
library(MCMCglmm)
library(ape)
library(phytools)

source('../../scripts/summarySE.R')

#set the working directory to the directory this script is in.
#setwd("~/what/ever/folder/you/are/working/from") 


# useful links -----

# correcting for phylogeny with MCMCglmm useful links:

# MCMCglmm overview and course notes:
# https://cran.r-project.org/web/packages/MCMCglmm/vignettes/Overview.pdf
# https://cran.r-project.org/web/packages/MCMCglmm/vignettes/CourseNotes.pdf

# Simple MCMCglmm tutorial without phylo-correction
# https://ourcodingclub.github.io/tutorials/mcmcglmm/

# Simple MCMCglmm with phylo-correction example but not for multiple observations per species
# https://github.com/TGuillerme/mulTree/blob/master/doc/Vanilla_flavoured_phylogenetic_analyses.Rmd

# MCMCglmm with phylo-correction example and multiple observations per species
# https://github.com/MPCMEvolution/MPCMArchive/wiki/11.-Chapter-11

# MCMCglmm with phylo-correction example for binary data at the end (not multiple observations)
# https://devillemereuil.legtux.org/wp-content/uploads/2012/12/tuto_en.pdf

# more examples
# http://ms.mcmaster.ca/~bolker/R/misc/foxchapter/bolker_chap.html
# http://www.wildanimalmodels.org/tiki-download_wiki_attachment.php?attId=24
# https://rstudio-pubs-static.s3.amazonaws.com/159004_213879bf716a4d44833689c03ccf17ab.html
# https://stackoverflow.com/questions/67941391/how-to-properly-code-a-scaled-inverse-wishart-prior-for-a-mcmcglmm-model


# Formatting the data for MCMCglmm ------


# read in QS gene data. This is a binary matrix of the presence/absence of 3 "complete" QS systems in the strains.
# 1 = QS system present, 0 = absent.
# strains are the column names - note they do not have SM at the start.
# for some reason error message appears with "incomplete final line...." but ignore there is nothing wrong with the file.
QS_data <- read.csv('quorumsensingsystem_presenceabsence.csv', row.names = 1)

# reverse the column names
QS_data <- QS_data[,order(ncol(QS_data):1)]

# transpose the QS_data dataframe
# make back into dataframe
QS_data <- as.data.frame(t(QS_data))

# make new column called strain so that we can merge dataframes by strain
QS_data$strain <- row.names(QS_data)

# remove X from strain names
QS_data$strain <- gsub("X", "", QS_data$strain)


# read in means summary inhibition zone diameter data 
# IZ_diff is the mean size of the inhibition zone for the row specified strain combinations.
diameter_data_means <- read.csv('../../data/raw_data/spotplating_data/dataset_IZspot_diameterdiffs_averages.csv', row.names = 1)
# keep only spot and lawn strain IDs and the mean inhibition zone diameter
diameter_data_means <- diameter_data_means[c("spot_strain", "soft_strain", "spot_geno", "soft_geno","IZ_diff")]

# in the diameter_data_means dataframe make copies of both the spot_strain and soft_strain columns. 
# These will be used to associate with the phylogeny for correction in the model, while still allowing for the 
# soft_strain and spot_strain columns to be included in the model as separate random effects which will account
# for variance produced as a result of the multiple observations per strain. and the new spot_strain_phylogeny and soft_strain_phylogeny will be used to account for phylogeny of the strains in the model. 
diameter_data_means$spot_strain_phylogeny <- diameter_data_means$spot_strain
diameter_data_means$soft_strain_phylogeny <- diameter_data_means$soft_strain

# combine QS_data and diameter_data_means dataframes 
# combine via spotted strain column
# make new dataframe with gene presence/absence for SPOTTED strain
IZ_means_spotpa <- merge(diameter_data_means, QS_data, 
                      by.x = c('spot_strain'),
                      by.y = c('strain'),
                      all.x = T)

# combine with diameter_data_means 
# combine via lawn strain column
# AGAIN make new dataframe with gene presence/absence for SOFT AGAR LAWN strain
IZ_means_lawnpa <- merge(diameter_data_means, QS_data, 
                         by.x = c('soft_strain'),
                         by.y = c('strain'),
                         all.x = T)



# read in phylogenetic tree of 24 strains
tree <- read.tree('305_core_genes_24subset.nw')

# view the tree
plot(tree)

# check the following:
# is it binary
is.binary.tree(tree)
# [1] TRUE
# is it ultrametric
is.ultrametric(tree)
# [1] FALSE
# is it rooted
is.rooted(tree)
# [1] TRUE - Checking if our tree is rooted, and that the root is non-zero length:
# Now we’ll set any zero-length branches to one-ten-thousandth of the tree size
tree$edge.length[tree$edge.length==0]<-max(nodeHeights(tree))*1e-4


# We can use the branch lengths of a non-ultrametric tree and transform them to be ultrametric. The first step is to see how long the tree is from root to tip. Then we can use that information to generate an ultrametric tree that is the same length. 
# gives you distances of each node from the base of the tree (the base is at zero).
# The first column of the output is the node Height from the base of the tree to the ancestral side of a given branch, and the second is the height from the base to the descendent side of each branch. If you take the maximum of this object (the greatest number), you have the height of the entire tree.
nodeHeights(tree)
max(nodeHeights(tree))

# use the chronopl function (part of the ape package) to transform the tree to a chronogram, which is another term for an ultrametric tree. 
ultra_tree <- chronopl(tree, lambda = 1, age.min=max(nodeHeights(tree)))
# now do the tree checks again:
is.ultrametric(ultra_tree)
# [1] TRUE
is.binary(ultra_tree)
# [1] TRUE
is.rooted(ultra_tree)
# [1] TRUE
# plot the tree
plot(ultra_tree)

# remove the first part of the node name before and up to the '-' so then it matches the strain names in the IZ_means_lawnpa and IZ_means_spotpa dataframes
ultra_tree$tip.label
ultra_tree$tip.label <- gsub(".*-","",ultra_tree$tip.label)
plot(ultra_tree)

# take the inverse of the phylogenic correlations which will calculate the branch lengths to be used in the model
inv.phylo <- inverseA(ultra_tree, nodes="TIPS",scale=TRUE)

# the data looks to be zero-inflated so we will fit a zero-inflated model
hist(IZ_means_lawnpa$IZ_diff, 
     #ylim = c(0,20),
     breaks = seq(0,20,0.1)
     )

# how much of the data is 0?
100*sum(IZ_means_lawnpa$IZ_diff == 0)/nrow(IZ_means_lawnpa)
# [1] 84.02778



# MCMCglmm set up -----------


# produce 8 regression models. 
# tra + vicibactin geneset x whether the gene presence absence is based on the spotted strain or the lawn strain x phylogenetic correction or not
# use the same settings for all regressions


# First set up our prior for the random (G) and residual variances (R) as a list within a list.
# We don't need to set up the priors for fixed effects (B) because the MCMCglmm does a good job of this for us.
# Because we don’t want our estimates in the model heavily influenced by our prior we use weakly informative prior values such as `V = 1` and `nu = 0.002`, as suggested in the links attached at the start of this section.
# Because we are using binary data, caution is advised, and it is suggested that we fix the residual priors to 1
prior1 <- list(R = list(V = 1, fix = 1),
               G = list(G1 = list(V = 1, nu = 1,  alpha.mu = 0, alpha.V = 1000),
                        G2 = list(V = 1, nu = 1,  alpha.mu = 0, alpha.V = 1000)
                        ))

prior_nophy <- list(R = list(V = 1, fix = 1),
               G = list(G1 = list(V = 1, nu = 1,  alpha.mu = 0, alpha.V = 1000)
               ))

prior2 <- list(R = list(V = 1, fix = 1),
               G = list(G1 = list(V = 1, nu = 1,  alpha.mu = 0, alpha.V = 1000),
                        G2 = list(V = 1, nu = 1,  alpha.mu = 0, alpha.V = 1000),
                        G3 = list(V = 1, nu = 1,  alpha.mu = 0, alpha.V = 1000)
               ))

# number of iterations, burning and thinning parameters
nsamp <- 1000
THIN <- 5000
BURNIN <- 50000
NITT <- BURNIN + THIN*nsamp


### model for Vicibactin presence/absence in spot strain - corrected for phylogeny in spot strain. ----

# Make a MCMCglmm model 
# define the covariance structure in correlation argument
model_vicspot <- MCMCglmm(fixed = vicibactin ~ IZ_diff, 
                          random = ~ spot_strain_phylogeny + spot_strain, 
                          family = "categorical", 
                          prior = prior1,
                          ginverse=list(spot_strain_phylogeny = inv.phylo$Ainv),
                          data = IZ_means_spotpa,
                          nitt = NITT,
                          burnin = BURNIN,
                          thin = THIN) 


# obtain the 'trace' of the sampling (to check for convergence and auto-correlation) and posterior density of each parameter.
# the graphs should look like horizontally stable "fluffy caterpillars".
# plot fixed effects
plot(model_vicspot$Sol)
#plot random effects
plot(model_vicspot$VCV)

# check for autocorrelation in traces. autocorr. diag calculates the autocorrelation function for the Markov chain mcmc. obj at the lags given by lags
# fixed and residual effects
autocorr.diag(model_vicspot$Sol)
# (Intercept)   IZ_diff
# Lag 0       1.00000000 1.0000000
# Lag 5000    0.16341068 0.9694068
# Lag 25000   0.14129506 0.8732675
# Lag 50000   0.15130624 0.7860252
# Lag 250000  0.05473707 0.2767877

#acf plot for the first fixed estimate in our model (the intercept)
acf(model_vicspot$Sol[,1], lag.max = 20)
# inhibition zone size fixed effect
acf(model_vicspot$Sol[,2], lag.max = 20)

# random effects
autocorr.diag(model_vicspot$VCV)
# spot_strain_phylogeny spot_strain units
# Lag 0                  1.0000000   1.0000000   NaN
# Lag 5000               0.7374403   0.5833175   NaN
# Lag 25000              0.7192795   0.4924022   NaN
# Lag 50000              0.7117516   0.3992556   NaN
# Lag 250000             0.5480141   0.2663834   NaN

# acf plot for the first random term in our model (the phylogeny term)
acf(model_vicspot$VCV[,1], lag.max = 20)


# view results of the model
summary(model_vicspot)

# Iterations = 50001:5045001
# Thinning interval  = 5000
# Sample size  = 1000 
# 
# DIC: 0.2798471 
# 
# G-structure:  ~spot_strain_phylogeny
# 
#                       post.mean l-95% CI u-95% CI eff.samp
# spot_strain_phylogeny     17984    331.2    49717     3.95
# 
# ~spot_strain
# 
#              post.mean  l-95% CI u-95% CI eff.samp
# spot_strain     436.2 0.0002601     1958    38.51
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units         1        1        1        0
# 
# Location effects: vicibactin ~ IZ_diff 
# 
#             post.mean  l-95% CI  u-95% CI eff.samp pMCMC
# (Intercept)  -72.8446 -234.0560   76.7614   118.94 0.260
# IZ_diff        0.2945   -3.9912    4.8072    13.34 0.946


# calculate the posterior probability of the phylogenetic signal lambda
lambda_vicspot <- model_vicspot$VCV[,'spot_strain_phylogeny']/
  (model_vicspot$VCV[,'spot_strain_phylogeny'] + model_vicspot$VCV[,'units'])

# calculate the posterior mean (mean of the posterior distribution)
# 1 is following a brownian motion model - so phylogenetic signal is high.
# 0 is no phylogenetic signal
mean(lambda_vicspot)
# [1] 0.9998028
# posterior mode
posterior.mode(lambda_vicspot)
# 0.9999653
# 95% credible interval of lambda
HPDinterval(lambda_vicspot)
# lower     upper
# var1 0.9991302 0.9999906
# attr(,"Probability")
# [1] 0.95





### model for vicibactin presence/absence in spot strain - NOT corrected for phylogeny in spot strain.----


# Make a MCMCglmm model 
# define the covariance structure in correlation argument
model_vicspot_nophy <- MCMCglmm(fixed = vicibactin ~ IZ_diff, 
                                random = ~ spot_strain, 
                                family = "categorical", 
                                prior = prior_nophy,
                                data = IZ_means_spotpa,
                                nitt = NITT,
                                burnin = BURNIN,
                                thin = THIN)


# obtain the 'trace' of the sampling (to check for convergence and auto-correlation) and posterior density of each parameter.
# plot fixed effects
plot(model_vicspot_nophy$Sol)
#plot random effects
plot(model_vicspot_nophy$VCV)

# check for autocorrelation in traces. autocorr. diag calculates the autocorrelation function for the Markov chain mcmc. obj at the lags given by lags
# fixed and residual effects
autocorr.diag(model_vicspot_nophy$Sol)
# (Intercept)   IZ_diff
# Lag 0        1.0000000 1.0000000
# Lag 5000     0.4663682 0.7868087
# Lag 25000    0.4247404 0.4590207
# Lag 50000    0.3957113 0.2465722
# Lag 250000   0.2120729 0.1357849

#acf plot for the first fixed estimate in our model (the intercept)
acf(model_vicspot_nophy$Sol[,1], lag.max = 20)
# inhibition zone size fixed effect
acf(model_vicspot_nophy$Sol[,2], lag.max = 20)

# random effects
autocorr.diag(model_vicspot_nophy$VCV)
# spot_strain units
# Lag 0        1.0000000   NaN
# Lag 5000     0.6261420   NaN
# Lag 25000    0.5762481   NaN
# Lag 50000    0.5838761   NaN
# Lag 250000   0.3714414   NaN

# acf plot for the first random term in our model (the phylogeny term)
acf(model_vicspot_nophy$VCV[,1], lag.max = 20)

# view results
summary(model_vicspot_nophy)

# printed result

# Iterations = 50001:5045001
# Thinning interval  = 5000
# Sample size  = 1000 
# 
# DIC: 1.259176 
# 
# G-structure:  ~spot_strain
# 
#             post.mean l-95% CI u-95% CI eff.samp
# spot_strain      4455    302.9     9564    9.916
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units         1        1        1        0
# 
# Location effects: vicibactin ~ IZ_diff 
# 
#             post.mean  l-95% CI  u-95% CI eff.samp  pMCMC    
# (Intercept) -49.87829 -87.01003 -10.23022    22.41 <0.001 ***
# IZ_diff       0.03326  -1.33018   2.13654    75.86  0.968    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# compare models by DIC values
model_vicspot$DIC
# [1] 0.2798471
model_vicspot_nophy$DIC
# [1] 1.259176



### model for vicibactin presence/absence in lawn strain - corrected for phylogeny in lawn strain.---- 

# Make a MCMCglmm model .
# define the covariance structure in correlation argument
model_viclawn <- MCMCglmm(fixed = vicibactin ~ IZ_diff, 
                          random = ~ soft_strain_phylogeny + soft_strain, 
                          family = "categorical", 
                          prior = prior1,
                          ginverse=list(soft_strain_phylogeny = inv.phylo$Ainv),
                          data = IZ_means_lawnpa,
                          nitt = NITT,
                          burnin = BURNIN,
                          thin = THIN) 



# obtain the 'trace' of the sampling (to check for convergence and auto-correlation) and posterior density of each parameter.
# the graphs should look like horizontally stable "fluffy caterpillars".
# plot fixed effects
plot(model_viclawn$Sol)
#plot random effects
plot(model_viclawn$VCV)

# check for autocorrelation in traces. autocorr. diag calculates the autocorrelation function for the Markov chain mcmc. obj at the lags given by lags
# fixed and residual effects
autocorr.diag(model_viclawn$Sol)
# (Intercept)   IZ_diff
# Lag 0        1.0000000 1.0000000
# Lag 5000     0.1787369 0.9737656
# Lag 25000    0.1522805 0.9051713
# Lag 50000    0.1904282 0.8360595
# Lag 250000   0.1409063 0.5494185

#acf plot for the first fixed estimate in our model (the intercept)
acf(model_viclawn$Sol[,1], lag.max = 20)
# inhibition zone size fixed effect
acf(model_viclawn$Sol[,2], lag.max = 20)

# random effects
autocorr.diag(model_viclawn$VCV)
# soft_strain_phylogeny soft_strain units
# Lag 0                  1.0000000   1.0000000   NaN
# Lag 5000               0.6787479   0.5865283   NaN
# Lag 25000              0.6339458   0.4303658   NaN
# Lag 50000              0.5472347   0.2852276   NaN
# Lag 250000             0.3389404   0.0240114   NaN

# acf plot for the first random term in our model (the phylogeny term)
acf(model_viclawn$VCV[,1], lag.max = 20)


# view results of the model
summary(model_viclawn)

# Iterations = 50001:5045001
# Thinning interval  = 5000
# Sample size  = 1000 
# 
# DIC: 0.4486576 
# 
# G-structure:  ~soft_strain_phylogeny
# 
#                       post.mean l-95% CI u-95% CI eff.samp
# soft_strain_phylogeny      6118    110.7    16751    10.54
# 
# ~soft_strain
# 
#             post.mean l-95% CI u-95% CI eff.samp
# soft_strain       192 0.001262    901.9    66.86
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units         1        1        1        0
# 
# Location effects: vicibactin ~ IZ_diff 
# 
#             post.mean l-95% CI u-95% CI eff.samp pMCMC
# (Intercept)   -43.852 -155.514   59.780    41.74 0.310
# IZ_diff        -6.196  -16.564    3.318    10.07 0.294


# calculate the posterior probability of the phylogenetic signal lambda
lambda_viclawn <- model_viclawn$VCV[,'soft_strain_phylogeny']/
  (model_viclawn$VCV[,'soft_strain_phylogeny'] + model_viclawn$VCV[,'units'])

# calculate the posterior mean (mean of the posterior distribution)
# 1 is following a brownian motion model - so phylogenetic signal is high.
# 0 is no phylogenetic signal
mean(lambda_viclawn)
# [1] 0.9993464

# posterior mode
posterior.mode(lambda_viclawn)
# 0.9998933

# 95% credible interval of lambda
HPDinterval(lambda_viclawn)
# lower     upper
# var1 0.9971746 0.9999723
# attr(,"Probability")
# [1] 0.95





### model for vicibactin presence/absence in lawn strain - NOT corrected for phylogeny in lawn strain.----


# Make a MCMCglmm model
# define the covariance structure in correlation argument
model_viclawn_nophy <- MCMCglmm(fixed = vicibactin ~ IZ_diff, 
                                random = ~ soft_strain, 
                                family = "categorical", 
                                prior = prior_nophy,
                                data = IZ_means_lawnpa,
                                nitt = NITT,
                                burnin = BURNIN,
                                thin = THIN)


# obtain the 'trace' of the sampling (to check for convergence and auto-correlation) and posterior density of each parameter.
# plot fixed effects
plot(model_viclawn_nophy$Sol)
#plot random effects
plot(model_viclawn_nophy$VCV)

# check for autocorrelation in traces. autocorr. diag calculates the autocorrelation function for the Markov chain mcmc. obj at the lags given by lags
# fixed and residual effects
autocorr.diag(model_viclawn_nophy$Sol)
# (Intercept)   IZ_diff
# Lag 0        1.0000000 1.0000000
# Lag 5000     0.3683229 0.9276720
# Lag 25000    0.4147777 0.7730613
# Lag 50000    0.3273075 0.6826259
# Lag 250000   0.2116437 0.2879652

#acf plot for the first fixed estimate in our model (the intercept)
acf(model_viclawn_nophy$Sol[,1], lag.max = 20)
# inhibition zone size fixed effect
acf(model_viclawn_nophy$Sol[,2], lag.max = 20)

# random effects
autocorr.diag(model_viclawn_nophy$VCV)
# soft_strain units
# Lag 0        1.0000000   NaN
# Lag 5000     0.7725660   NaN
# Lag 25000    0.7273274   NaN
# Lag 50000    0.7030394   NaN
# Lag 250000   0.4076166   NaN

# acf plot for the first random term in our model (the phylogeny term)
acf(model_viclawn_nophy$VCV[,1], lag.max = 20)

# view results
summary(model_viclawn_nophy)


# printed result
# Iterations = 50001:5045001
# Thinning interval  = 5000
# Sample size  = 1000 
# 
# DIC: 1.455214 
# 
# G-structure:  ~soft_strain
# 
#             post.mean l-95% CI u-95% CI eff.samp
# soft_strain      3452    180.5     9346     7.55
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units         1        1        1        0
# 
# Location effects: vicibactin ~ IZ_diff 
# 
#             post.mean l-95% CI u-95% CI eff.samp pMCMC   
# (Intercept)   -38.240  -65.823   -9.795    33.77 0.006 **
# IZ_diff        -2.383   -9.677    2.921    16.13 0.574   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# compare models by DIC values
model_viclawn$DIC
# [1] 0.4486576
model_viclawn_nophy$DIC
# [1] 1.455214


# model for tra genes presence/absence in lawn strain - corrected for phylogeny in spot strain ---------------

# in this next model, we correct for phylogeny in the spotted strain, as we technically account for phylogeny in the lawn strain by only doing the test on OA lawn strains.

# filter lawn_pa data for only OA strains. 
IZ_means_lawnpa2 <- IZ_means_lawnpa[IZ_means_lawnpa$soft_geno == "OA", ]

# calculate the mean IZ_diff for each strain of OA across all spotted strains
IZ_means_lawnpa3 <- summarySE(IZ_means_lawnpa2, measurevar = 'IZ_diff', 
         groupvars = c("soft_strain", "soft_geno", 'tra'))


nsamp2 <- 1000 
THIN2 <- 5000
BURNIN2 <- 50000 
NITT2 <- BURNIN2 + THIN2*nsamp2

model_trasoftOA <- MCMCglmm(fixed = tra ~ IZ_diff, 
                          random = ~ spot_strain_phylogeny + spot_strain, # spot strain phylo tho
                          family = "categorical", 
                          prior = prior1,
                          ginverse=list(spot_strain_phylogeny = inv.phylo$Ainv),
                          data = IZ_means_lawnpa2, # soft OA strain tra pres/abs
                          nitt = NITT2,
                          burnin = BURNIN2,
                          thin = THIN2) 


# obtain the 'trace' of the sampling (to check for convergence and auto-correlation) and posterior density of each parameter.
# the graphs should look like horizontally stable "fluffy caterpillars".
# plot fixed effects
plot(model_trasoftOA$Sol)
#plot random effects
plot(model_trasoftOA$VCV)

# check for autocorrelation in traces. autocorr. diag calculates the autocorrelation function for the Markov chain mcmc. obj at the lags given by lags
# fixed and residual effects
autocorr.diag(model_trasoftOA$Sol)
# (Intercept)   IZ_diff
# Lag 0       1.0000000000 1.0000000
# Lag 5000    0.3334194930 0.9619851
# Lag 25000   0.2320114407 0.8261159
# Lag 50000   0.1856722109 0.7151219
# Lag 250000 -0.0001084637 0.1482906

#acf plot for the first fixed estimate in our model (the intercept)
acf(model_trasoftOA$Sol[,1], lag.max = 20)
# inhibition zone size fixed effect
acf(model_trasoftOA$Sol[,2], lag.max = 20)

# random effects
autocorr.diag(model_trasoftOA$VCV)
# spot_strain_phylogeny  spot_strain units
# Lag 0                 1.00000000  1.000000000   NaN
# Lag 5000              0.39372353  0.604971260   NaN
# Lag 25000             0.23220694  0.257340446   NaN
# Lag 50000             0.16825089  0.045037014   NaN
# Lag 250000            0.02036452 -0.007566601   NaN

# acf plot for the first random term in our model (the phylogeny term)
acf(model_trasoftOA$VCV[,1], lag.max = 20)


# view results of the model
summary(model_trasoftOA)

# Iterations = 50001:5045001
# Thinning interval  = 5000
# Sample size  = 1000 
# 
# DIC: 32.08732 
# 
# G-structure:  ~spot_strain_phylogeny
# 
#                       post.mean l-95% CI u-95% CI eff.samp
# spot_strain_phylogeny     570.5   0.2346     1699    170.7
# 
# ~spot_strain
# 
#              post.mean  l-95% CI u-95% CI eff.samp
# spot_strain     80.53 5.349e-06    361.7    149.3
# 
# R-structure:  ~units
# 
#        post.mean l-95% CI u-95% CI eff.samp
# units         1        1        1        0
# 
# Location effects: tra ~ IZ_diff 
# 
#             post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)  -29.1480 -62.8012  -0.3066    70.91  0.032 *  
# IZ_diff        8.6868   3.1231  13.3649    20.20 <0.001 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# calculate the posterior probability of the phylogenetic signal lambda
lambda_trasoftOA <- model_trasoftOA$VCV[,'spot_strain_phylogeny']/
  (model_trasoftOA$VCV[,'spot_strain_phylogeny'] + model_trasoftOA$VCV[,'units'])

# calculate the posterior mean (mean of the posterior distribution)
# 1 is following a brownian motion model - so phylogenetic signal is high.
# 0 is no phylogenetic signal
mean(lambda_trasoftOA )
# [1] 0.992775

# posterior mode
posterior.mode(lambda_trasoftOA )
# 0.9984571

# 95% credible interval of lambda
HPDinterval(lambda_trasoftOA )
# lower     upper
# var1 0.9834796 0.9998075
# attr(,"Probability")
# [1] 0.95


