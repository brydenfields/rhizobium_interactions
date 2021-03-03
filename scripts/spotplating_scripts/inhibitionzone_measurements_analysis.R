# Author: Bryden Fields
# Email: bryden.fields@york.ac.uk
# Last updated: 03/02/2021

# Initial set-up -----------

setwd('/Users/brydenfields/Documents/Publications/2021_Rhizobiuminteractions_paper/scripts/spotplating_scripts')

library(tidyverse)
library(plyr); library(dplyr)
library(reshape2)
library(devtools)
library(ggpubr)
library(FSA)

source('../summarySE.R')


# For reference only: Raw files were read, combined and formatted --------------

# # read in results file names for the first replicate
# 
# # # rep1
# # file_list <- Sys.glob("/Users/brydenfields/Documents/PhD/spotplating_images/0619_exp/rep1/72h/rep1_72h_BW/*_results.csv")
# # # rep2
# # file_list <- append(file_list, Sys.glob("/Users/brydenfields/Documents/PhD/spotplating_images/0619_exp/rep2/72h/rep2_72h_BW/*_results.csv"), after = length(file_list))
# # # rep3
# # file_list <- append(file_list, Sys.glob("/Users/brydenfields/Documents/PhD/spotplating_images/0619_exp/rep3/72h/rep3_72h_BW/*_results.csv"), after = length(file_list))
# # # rep4
# # file_list <- append(file_list, Sys.glob("/Users/brydenfields/Documents/PhD/spotplating_images/0619_exp/rep4/72h/rep4_72h_BW/*_results.csv"), after = length(file_list))
# 
# 
# # rep1
# file_list <- Sys.glob("/Users/brydenfields/Documents/PhD/supernatant/spotplating_images/0619_exp/0320_analysis/rep1_72h_BW/*_results.csv")
# # rep2
# file_list <- append(file_list, Sys.glob("/Users/brydenfields/Documents/PhD/supernatant/spotplating_images/0619_exp/0320_analysis/rep2_72h_BW/*_results.csv"), after = length(file_list))
# # rep3
# file_list <- append(file_list, Sys.glob("/Users/brydenfields/Documents/PhD/supernatant/spotplating_images/0619_exp/0320_analysis/rep3_72h_BW/*_results.csv"), after = length(file_list))
# 
# # code taken from:
# # https://psychwire.wordpress.com/2011/06/03/merge-all-files-in-a-directory-using-r-into-a-single-dataframe/
# # The final step is to iterate through the list of files in the current working directory and put them together to form a dataframe.
# # When the script encounters the first file in the file_list, it creates the main dataframe to merge everything into (called dataset here). This is done using the !exists conditional:
# # If dataset already exists, then a temporary dataframe called temp_dataset is created and added to dataset. The temporary dataframe is removed when we’re done with it using the rm(temp_dataset) command. 
# # If dataset doesn’t exist (!exists is true), then we create it.
# 
# for (file in file_list){
#   
#   # if the merged dataset doesn't exist, create it
#   if (!exists("dataset")){
#     dataset <- read.csv(file)
#   }
#   
#   # if the merged dataset does exist, append to it
#   if (exists("dataset")){
#     temp_dataset <-read.csv(file)
#     dataset<-rbind(dataset, temp_dataset)
#     rm(temp_dataset)
#   }
#   
# }
# 
# # the only poor thing about this above script that I copied is that it copies the first dataframe twice.
# # therefore remove the top 24 rows
# dataset <- dataset[25:nrow(dataset),]
# 
# dataset$imagejID <- NULL
# 
# 
# # checking number of plates is correct when called in noinhib_files at this point (now moved elsewhere)
# # 24 strains x 2 plates = 48 ≠ 2 TY control plates = 50
# # rep 1 = 48 (lost 152B because accidentally did rep for 152A), rep 2 = 50, rep3 = 50
# # 116 + 32 (remember no 152B rep1) = 148 
# 
# # take datasets with IZ
# # filter out both the IZ and the spot measurements
# dataset_IZ <- dataset %>%
#   filter(measurement %in% 'IZ')
# 
# dataset_spot <- dataset %>%
#   filter(measurement %in% 'spot')
# 
# # combine the columns together
# dataset_IZspot <- merge(dataset_IZ, dataset_spot, by = c("soft_strain", "spot_strain", "plate", "rep"), 
#                         all.x = T)
# 
# rm(dataset)
# rm(dataset_IZ)
# rm(dataset_spot)
# 
# write.csv(dataset_IZspot, '/Users/brydenfields/Documents/PhD/supernatant/spotplating_images/0619_exp/0320_analysis/dataset_IZspot_0320.csv')



# Format data to correlate IZ against spot diameter ------------

# read in data
dataset_IZspot <- read.csv('../../data/raw_data/spotplating_data/dataset_IZspot_0320.csv', row.names = 1)

# filter for only the data where we calculated the inhibition zone and culture spot diameters
# thats rows that are full with data rather than my manually added rows for pairwise comparisons
# with no inhibition zones. 
dataset_IZspot2 <- dataset_IZspot[complete.cases(dataset_IZspot), ]

# then correlate the IZ ferret against spot Ferret
# calculate: IZ ferret - spot ferret
cor.test(dataset_IZspot2$Feret.x ,dataset_IZspot2$Feret.y)


# Pearson's product-moment correlation
# 
# data:  dataset_IZspot2$Feret.x and dataset_IZspot2$Feret.y
# t = 5.768, df = 242, p-value = 2.439e-08
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.2322066 0.4534425
# sample estimates:
#       cor 
# 0.3476542 



# regression analysis
# linear regression looking at the size of inhibition zone by liquid culture spot size in the dataframe called dataset_IZspot2 (which only looks at the combinations that produce inhibition zones).
mod <- lm(data = dataset_IZspot2, Feret.x ~ Feret.y)
# summary shows the main results of the model including the coefficients, p-values, r2 values, f-statistic values etc.
summary(mod)
# testing for normality and homogeneity as done previously with the one way anova. 
plot(mod, which=1)
plot(mod, which=2)
shapiro.test(mod$residuals)
hist(dataset_IZspot2$Feret.x)
hist(dataset_IZspot2$Feret.y)

# Call:
#   lm(formula = Feret.x ~ Feret.y, data = dataset_IZspot2)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -11.781  -5.921   1.532   4.365  11.399 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   5.6402     1.6486   3.421 0.000732 ***
# Feret.y       2.1873     0.3792   5.768 2.44e-08 ***
#   ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 5.892 on 242 degrees of freedom
# Multiple R-squared:  0.1209,	Adjusted R-squared:  0.1172 
# F-statistic: 33.27 on 1 and 242 DF,  p-value: 2.439e-08

# for a 1 unit increase in spot diameter then there is a 5 unit increase in the inhibition zone diameter.


# graph: correlation plot
(graph_cor <- ggplot(dataset_IZspot2, aes(Feret.x, Feret.y)) +
  geom_point(size=3, alpha = 0.2) +
    coord_fixed() +
  geom_smooth(method=lm, se=FALSE, col = "black") +
  #stat_cor(method = "pearson", label.x = 3, label.y = 10) +
  xlab("Diameter of Inhibition Zone (mm)") +
  ylab("Diameter of bacterial\nculture spot (mm)") + 
  theme_bw() +
  theme(legend.title = element_blank(), 
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title=element_text(size = 20),
        legend.text = element_text(size = 16))
)
#function to save specific ggplot graph.  additional parameters(width=,height=,units='cm')
ggsave('../../data/intermediate_data/spotplating_data/correlation_IZvsspot_lmline_v1.pdf', plot = graph_cor, width = 25, height = 20, unit = 'cm')


# remove results below 1 and see if this changes the correlation
dataset_IZspot3 <- dataset_IZspot2 %>%
  filter(Feret.y > 1)

# redo correlation
cor.test(dataset_IZspot3$Feret.x ,dataset_IZspot3$Feret.y)
# Pearson's product-moment correlation
# 
# data:  dataset_IZspot3$Feret.x and dataset_IZspot3$Feret.y
# t = 2.7295, df = 236, p-value = 0.006821
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.04886139 0.29552410
# sample estimates:
#       cor 
# 0.1749363 

# redo regression for chapter stats
# linear regression looking at the size of inhibition zone by liquid culture spot size
mod3 <- lm(data = dataset_IZspot3, Feret.x ~ Feret.y)
# summary shows the main results of the model including the coefficients, p-values, r2 values, f-statistic values etc.
summary(mod3)
# testing for normality and homogeneity
plot(mod3, which=1)
plot(mod3, which=2)
shapiro.test(mod3$residuals)

# Call:
#   lm(formula = Feret.x ~ Feret.y, data = dataset_IZspot3)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -10.432  -5.733   1.095   4.301  10.853 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   9.6807     2.0770   4.661 5.27e-06 ***
# Feret.y       1.2881     0.4719   2.730  0.00682 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 5.844 on 236 degrees of freedom
# Multiple R-squared:  0.0306,	Adjusted R-squared:  0.0265 
# F-statistic:  7.45 on 1 and 236 DF,  p-value: 0.006821


# replot with subsetted data
(graph_cor2 <- ggplot(dataset_IZspot3, aes(Feret.x, Feret.y)) +
    geom_point(size=3, alpha = 0.2) +
    coord_fixed() +
    geom_smooth(method=lm, se=FALSE, col = "black") +
    #stat_cor(method = "pearson", label.x = 3, label.y = 10) +
    xlab("Diameter of Inhibition Zone (mm)") +
    ylab("Diameter of bacterial\nculture spot (mm)") + 
    theme_bw() +
    theme(legend.title = element_blank(), 
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title=element_text(size = 20),
          legend.text = element_text(size = 16))
)
#function to save specific ggplot graph.  additional parameters(width=,height=,units='cm')
ggsave('../../data/intermediate_data/spotplating_data/correlation_IZvsspot_lmline_v2.pdf', plot = graph_cor2, width = 25, height = 20, unit = 'cm')





# Graphs for average Ferret diameter (IZ - spot diameter = IZ_diff) ---------

# make new dataset to work with
all_data <- dataset_IZspot

# must add columns for soft_strain and spot_strain groups (OA,OC,OE,CC)
all_data$soft_geno <- all_data$soft_strain
all_data$soft_geno <- factor(all_data$soft_geno, levels = c('152B', '137B','152A','145B','154C','144A',
                                                                        '147A', '158', '170C', '157B', '165A','122A',
                                                                        '126B','149A','135B','135A','159','168A',
                                                                        '41','53','57','74','77','67'),
                                   labels = c('OA','OA','OA','OA','OA','OA',
                                              'OC','OC','OC','OC','OC','OC',
                                              'OC','OE','OE','OE','OE','OE',
                                              'CC','CC','CC','CC','CC','CC'))


all_data$spot_geno <- all_data$spot_strain 
all_data$spot_geno <- factor(all_data$spot_geno, levels = c('152B', '137B','152A','145B','154C','144A',
                                                                        '147A', '158', '170C', '157B', '165A','122A',
                                                                        '126B','149A','135B','135A','159','168A',
                                                                        '41','53','57','74','77','67'),
                                   labels = c('OA','OA','OA','OA','OA','OA',
                                              'OC','OC','OC','OC','OC','OC',
                                              'OC','OE','OE','OE','OE','OE',
                                              'CC','CC','CC','CC','CC','CC'))


# make a new column which will contain the diameter of the inhibition zone - the bacterial culture spot diameter.
# Ferret.x = IZ diameter
# Ferret.y = spot diameter
all_data$IZ_diff <- all_data$Feret.x - all_data$Feret.y
# if there is no IZ then it should be 0?
# so convert NA to 0 for IZ_diff
# this will make a difference when averaging the IZ_diff
# all_data$IZ_diff[is.na(all_data$IZ_diff)] <- 0

# write csv for data
write.csv(all_data, '../../data/raw_data/spotplating_data/dataset_IZspot_diameterdiffs.csv')

# get an average IZ_diff diameter for each combination
IZdiff_summary <- summarySE(all_data, measurevar = 'IZ_diff', 
                            groupvars = c("soft_strain", "spot_strain", "soft_geno", "spot_geno"))


# work out the number of combinations that produce inhibition zones
nrow(IZdiff_summary[IZdiff_summary$IZ_diff > 0,])
# answer = 92
nrow(IZdiff_summary[IZdiff_summary$IZ_diff == 0,])
# answer 484
# 92 + 484 = 576 (which is 24 strains x 24 strains combinations)

# write csv for averages for combinations
write.csv(IZdiff_summary, '../../data/raw_data/spotplating_data/dataset_IZspot_diameterdiffs_averages.csv')

# add SM to strain names
IZdiff_summary$spot_strain <- paste0("SM", IZdiff_summary$spot_strain)
IZdiff_summary$soft_strain <- paste0("SM", IZdiff_summary$soft_strain)

#write out the strain order that should appear on the graphs (should be the same as the supernatant heatmaps)
IZdiff_summary$soft_strain<- factor(IZdiff_summary$soft_strain, levels = c('SM152B', 'SM137B','SM152A','SM145B','SM154C','SM144A',
                                                            'SM147A', 'SM158', 'SM170C', 'SM157B', 'SM165A','SM122A',
                                                            'SM126B','SM149A','SM135B','SM135A','SM159','SM168A',
                                                            'SM41','SM53','SM57','SM74','SM77','SM67'))
IZdiff_summary$spot_strain <- factor(IZdiff_summary$spot_strain, levels = c('SM152B', 'SM137B','SM152A','SM145B','SM154C','SM144A',
                                                            'SM147A', 'SM158', 'SM170C', 'SM157B', 'SM165A','SM122A',
                                                            'SM126B','SM149A','SM135B','SM135A','SM159','SM168A',
                                                            'SM41','SM53','SM57','SM74','SM77','SM67'))


# circles get larger with diameter matrix
spot.theme <- list(
  theme_classic(),
  theme(axis.ticks.x=element_blank(), axis.text.x=element_text(size = 19, angle = 90, hjust = 0)),
  theme(axis.ticks.y=element_blank(), axis.text.y=element_text(size = 19)),
  theme(axis.line=element_blank()),
  theme(text = element_text(size = 32)),
  #theme(legend.position = "none"),
  theme(plot.margin = unit(c(10,10,10,10), "mm")),
  scale_size_continuous(range = c(-0.3, 15))
  #scale_x_discrete(position = "top")
  )

h <- c(6.5, 13.5, 18.5)
v <- c(6.5, 13.5, 18.5)
(graph1 <- ggplot(data=IZdiff_summary, aes(spot_strain, soft_strain)) + 
  geom_point(aes(size = IZ_diff), colour = "grey") +
  geom_hline(yintercept = h, size = 0.75, alpha = 0.5) +
  geom_vline(xintercept = v, size = 0.75, alpha = 0.5) +
  labs(x = 'Spotted strain',
       y = 'Soft agar lawn strain',
       size = 'Inhibition\nZone\nDiameter\n(mm)') +
  spot.theme
)
#function to save specific ggplot graph.  additional parameters(width=,height=,units='cm')
ggsave(plot = graph1, '../../data/intermediate_data/spotplating_data/inhibitionzone_spotmap.pdf', width = 35, height = 35, unit = 'cm')



# remake inhibition zone diameters but by genospecies/environment groups instead
IZdiff_geno_summary <- summarySE(all_data, measurevar = 'IZ_diff', 
                            groupvars = c("soft_geno", "spot_geno"))

# write genospecies summary into a csv file
write.csv(IZdiff_geno_summary, '../../data/intermediate_data/spotplating_data/dataset_IZspot_diameterdiffs_averages_genospecies.csv')



# Calculate min and max inhibition zone diameter for SM145B, SM154C and SM144A soft agar strains -------

IZdiff_summary_3s <- IZdiff_summary %>%
  filter(soft_strain %in% c('SM145B', 'SM154C', 'SM144A'))

OA_sizes <- IZdiff_summary_3s[IZdiff_summary_3s$IZ_diff != 0, ]
min(OA_sizes$IZ_diff)
# [1] 7.107333
max(OA_sizes$IZ_diff)
# [1] 19.84667



# Correlate inhibition zone diameter with the initial liquid culture spot ODs ------------
  
# read in initial OD data
OD_data <- read.csv('../../data/raw_data/spotplating_data/inoculum_inital_ODs_allreps.csv')
# remove genospecies column
OD_data$genospecies <- NULL
# rename columns to match diameter_data for merging and consistency
names(OD_data) <- c('spot_strain', 'rep','inoc_OD_blanked')
# remove SM from the strain names in OD_data to match with diameter_data
OD_data$spot_strain <- gsub('SM', '', OD_data$spot_strain)

# read in inhibition zone diameter data 
diameter_data <- read.csv('../../data/raw_data/spotplating_data/dataset_IZspot_diameterdiffs.csv', row.names = 1)

# merge diameter data with initial OD data
# must make sure it matches up to the correct replicate for the spot_strain
di_OD_data <- merge(diameter_data[c(1,2,4,26,27)], OD_data, by = c('spot_strain', 'rep'), all.x = T)


# make ggplot scatter plot to show correlation between inhib diameter and initial OD
di_OD_data$spot_geno <- factor(di_OD_data$spot_geno, levels = c("OA", "OC", "OE", "CC"))
genoPalette <- c("#466EA9", "#4D9A7A", "#DA367E", "#9ACD32")

(graph_a <- ggplot(data = di_OD_data, aes(x = inoc_OD_blanked, y = IZ_diff, 
                                          col = spot_geno, fill = spot_geno)) +
    geom_point(size = 3, alpha = 0.6) +
    labs(x = expression(paste("Optical Density of inoculum for liquid culture spot (",OD[600],")",sep="")),
         y = "Inhibition zone diameter (mm)",
         colour = 'Genospecies group', fill = 'Genospecies group') +
    scale_color_manual(values = genoPalette) +
    scale_fill_manual(values = genoPalette) +
    guides(colour = guide_legend(override.aes = list(size=10))) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title=element_text(size = 20),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18))
)
#function to save specific ggplot graph.  additional parameters(width=,height=,units='cm')
ggsave('../../data/intermediate_data/spotplating_data/IZdiameterdiff_inocinitialOD_scatterplot_genospeciescol.pdf', plot = graph_a, width = 25, height = 20, units = 'cm')


# linear regression looking at the size of inhibition zone by liquid culture spot OD.
mod_iz <- lm(data = di_OD_data, IZ_diff ~ inoc_OD_blanked)
#mod_iz <- lm(data = di_OD_data, IZ_diff ~ inoc_OD_blanked * spot_geno)
# summary shows the main results of the model including the coefficients, p-values, r2 values, f-statistic values etc.
summary(mod_iz)
# testing for normality and homogeneity as done previously with the one way anova. 
plot(mod_iz, which=1)
plot(mod_iz, which=2)
shapiro.test(mod_iz$residuals)
hist(di_OD_data$IZ_diff)
hist(di_OD_data$inoc_OD_blanked)

# Call:
#   lm(formula = IZ_diff ~ inoc_OD_blanked, data = di_OD_data)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -1.790 -1.629 -1.583 -1.506 19.978 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)       2.0687     0.7016   2.949  0.00324 **
# inoc_OD_blanked  -6.4703     9.4545  -0.684  0.49384   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 4.454 on 1630 degrees of freedom
# Multiple R-squared:  0.0002872,	Adjusted R-squared:  -0.0003261 
# F-statistic: 0.4683 on 1 and 1630 DF,  p-value: 0.4938



# Even remove the 0s just to see what happens
di_OD_data2 <- di_OD_data[di_OD_data$IZ_diff != 0, ]

# simple linear model 
mod_iz2 <- lm(data = di_OD_data2, IZ_diff ~ inoc_OD_blanked)
summary(mod_iz2)

# Call:
#   lm(formula = IZ_diff ~ inoc_OD_blanked, data = di_OD_data2)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -10.502  -7.364   1.058   4.280  11.191 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        9.954      2.775   3.587 0.000405 ***
# inoc_OD_blanked    9.783     37.841   0.259 0.796227    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6.009 on 242 degrees of freedom
# Multiple R-squared:  0.0002761,	Adjusted R-squared:  -0.003855 
# F-statistic: 0.06683 on 1 and 242 DF,  p-value: 0.7962



# ANOVA to see if optical densities of spotted liquid cultures differed between genotypes and genospecies ---------

data_spotcultureOD <- read.csv('/Users/brydenfields/Documents/PhD/supernatant/spotplating_images/0619_exp/inoculum_inital_ODs_allreps.csv')

# remove repeat 4
data_spotcultureOD <- data_spotcultureOD %>%
  filter(!replicate %in% '4')

# one way anova 
# first we can test homogeneity of variance with bartlett.test
bartlett.test(data_spotcultureOD$OD_blanked, data_spotcultureOD$genospecies)

# Bartlett test of homogeneity of variances
# 
# data:  data_spotcultureOD$OD_blanked and data_spotcultureOD$genospecies
# Bartlett's K-squared = 12.573, df = 3, p-value = 0.005658

# we can also test normality of the myoglobin result distribution within each of the species groups in the 'species' variable. 
tapply(data_spotcultureOD$OD_blanked, data_spotcultureOD$genospecies, shapiro.test)
# $CC
# Shapiro-Wilk normality test
# data:  X[[i]]
# W = 0.92841, p-value = 0.1821
# 
# $OA
# Shapiro-Wilk normality test
# data:  X[[i]]
# W = 0.72622, p-value = 0.0001638
# 
# $OC
# Shapiro-Wilk normality test
# data:  X[[i]]
# W = 0.92221, p-value = 0.09606
# 
# $OE
# Shapiro-Wilk normality test
# data:  X[[i]]
# W = 0.95005, p-value = 0.5252

# non-parametric one-way anova: Kruskal-wallis
kruskal.test(data_spotcultureOD$OD_blanked ~ data_spotcultureOD$genospecies)
# Kruskal-Wallis rank sum test
# data:  data_spotcultureOD$OD_blanked by data_spotcultureOD$genospecies
# Kruskal-Wallis chi-squared = 20.036, df = 3, p-value = 0.0001669

# post hoc with Dunns test which will give adjusted p-values
DT_geno = dunnTest(OD_blanked ~ genospecies,
                   data = data_spotcultureOD,
                   method="bh")
DT_geno
# Dunn (1964) Kruskal-Wallis multiple comparison
# p-values adjusted with the Benjamini-Hochberg method.
# 
# Comparison          Z      P.unadj        P.adj
# 1    CC - OA -4.3892099 1.137632e-05 0.00006825793 ***
# 2    CC - OC -1.5271646 1.267201e-01 0.1520641 
# 3    OA - OC  3.0277338 2.463950e-03 0.007391851 **
# 4    CC - OE -1.8442657 6.514443e-02 0.09771665
# 5    OA - OE  2.3406818 1.924857e-02 0.03849713 *
# 6    OC - OE -0.4561983 6.482474e-01 0.6482474


# kruskal wallis for individual strains as well
kruskal.test(data_spotcultureOD$OD_blanked ~ data_spotcultureOD$strain)
# Kruskal-Wallis rank sum test
# data:  data_spotcultureOD$OD_blanked by data_spotcultureOD$strain
# Kruskal-Wallis chi-squared = 58.244, df = 23, p-value = 6.825e-05

# post hoc 
DT_strain = dunnTest(OD_blanked ~ strain,
                     data = data_spotcultureOD,
                     method="bh") 
DT_strain
# Dunn Kruskal-Wallis multiple comparison
# p-values adjusted with the Benjamini-Hochberg method.
# I have only printed the significant interactions to save space
# Comparison                      Z      P.unadj      P.adj
# 5   SM126B OC - SM135B OA     -2.958811419 0.0030882803 0.04735363
# 8   SM126B OC - SM137B OA     -3.427534020 0.0006090900 0.03362177
# 12  SM126B OC - SM144A OA     -3.271293153 0.0010705687 0.04221100
# 17  SM126B OC - SM145B OA     -3.300588316 0.0009648234 0.04438188
# 38  SM126B OC - SM152A OA     -2.929516257 0.0033949005 0.04684963
# 60  SM137B OA - SM154C OA     3.251763045 0.0011469159 0.03956860
# 61  SM144A OA - SM154C OA     3.095522178 0.0019646668 0.04518734
# 62  SM145B OA - SM154C OA     3.124817340 0.0017791550 0.04464062
# 107 SM126B OC - SM165A OC     -2.929516257 0.0033949005 0.04931540
# 158   SM137B OA - SM41 CC     3.154112503 0.0016098701 0.04443242
# 159   SM144A OA - SM41 CC     2.997871636 0.0027187216 0.04413924
# 160   SM145B OA - SM41 CC     3.027166798 0.0024685770 0.04258295
# 164   SM152B OA - SM41 CC     3.427534020 0.0006090900 0.04202721
# 182   SM152B OA - SM53 CC     3.027166798 0.0024685770 0.04542182
# 258   SM137B OA - SM77 CC     3.202937774 0.0013603337 0.04171690
# 259   SM144A OA - SM77 CC     3.046696907 0.0023137088 0.04561312
# 260   SM145B OA - SM77 CC     3.075992069 0.0020980343 0.04454288
# 264   SM152B OA - SM77 CC     3.476359291 0.0005082709 0.04676092

# box plot of inoculum strain density. 
(graph_c <- ggplot(data = data_spotcultureOD, aes(x = genospecies, y = OD_blanked, group = genospecies)) +
    geom_boxplot() +
    labs(x = 'Genospecies group',
         y = expression(paste("Optical Density of inoculum for liquid culture spot (",OD[600],")",sep=""))
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title=element_text(size = 15))
)
# save plot with ggsave
ggsave(plot = graph_c, '../../data/intermediate_data/spotplating_data/initialODs_spotcultures_reps13.pdf')




# correlate inhibition zone size to RGI as inoculant ----------------

# read in inhibition zone diameter data 
diameter_data_means <- read.csv('../../data/raw_data/spotplating_data/dataset_IZspot_diameterdiffs_averages.csv', row.names = 1)

# read in RGI supernatant data 
sup_data <- read.csv('../../data/raw_data/supernatant_data/bes_project_all_no1.csv')

# average all the RGIs reps for each pairwise treatment
sup_data <- sup_data[sup_data$time == 62 & sup_data$supernatant != "100% TY" & sup_data$supernatant != "50% TY",]

# averaging over the means as with MEM.
sup_data_means <- summarySE(sup_data, 
                            groupvars = c('inoc_geno','sup_geno', 'supernatant', 'inoculant'),
                            measurevar = c('comp_ind'))

# merge mean RGI of inoculant with the mean diameter data
#### IMPORTANT: SUPERNATANT STRAIN IS MATCHED WITH SOFT_STRAIN 
# if you grow better and the inoculant then you are more likely to produce larger inhibition zones as the spot
di_sup_means <- merge(sup_data_means[c(1:4, 6)], diameter_data_means[c(1:4, 6)], 
                      # by.x = c('supernatant', 'inoculant'),
                      # by.y = c('spot_strain', 'soft_strain'),
                      by.x = c('supernatant', 'inoculant'),
                      by.y = c('soft_strain', 'spot_strain'),
                      all.x = T)

# make IZ_diffs above 5 mm one category and those below another category
di_sup_means$IZ_group <- NA
di_sup_means[di_sup_means$IZ_diff > 5, c('IZ_group')] <- 'more than 5 mm'
di_sup_means[di_sup_means$IZ_diff <= 5, c('IZ_group')] <- 'less than 5 mm'
di_sup_means[di_sup_means$IZ_diff == 0, c('IZ_group')] <- '0 mm'

# colour palette for different IZ_diff groups
cbPalette <- c("#D55E00","#0072B2","#E69F00")

# make a graph correlating the mean RGI as inoculant to mean inhibition zone diameter
(graph_d <- ggplot(data = di_sup_means, aes(y = IZ_diff, x = comp_ind, )) +
    geom_point(aes(col = IZ_group, fill = IZ_group), size = 5, alpha = 0.6) +
    geom_smooth(method= lm, aes(group = IZ_group), colour="black", se = FALSE) +
    #geom_smooth(method= lm, colour="blue") +
    labs(y = 'Inhibition zone diameter produced by spotted strain (mm)', 
         x = 'Relative Growth Index of inoculant',
         colour = 'Inhibition zone \ndiameter',
         fill = 'Inhibition zone \ndiameter') +
    scale_colour_manual(values = cbPalette) +
    scale_fill_manual(values = cbPalette) +
    guides(colour = guide_legend(override.aes = list(size=7))) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title=element_text(size = 20),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18))
  
)
# ggsave(plot = graph_c, '/Users/brydenfields/Documents/PhD/supernatant/spotplating_images/0619_exp/0320_analysis/RGIinoculant_IZdiff_lmcorrelation_v3.pdf', width = 25, height = 20, units = 'cm')
ggsave(plot = graph_d, '../../data/intermediate_data/spotplating_data/RGIinoculant_IZdiff_lmcorrelation_v4.pdf', width = 25, height = 20, units = 'cm')


# regression for data with inhibition zones

# filter for inhibition zones above 5mm
disupmeans_m5 <- di_sup_means[di_sup_means$IZ_group == 'more than 5 mm', ]

# carry out regression
mod_IZsupmeans_m5 <- lm(IZ_diff ~ comp_ind, data = disupmeans_m5)
summary(mod_IZsupmeans_m5)

# Call:
#   lm(formula = IZ_diff ~ comp_ind, data = disupmeans_m5)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -6.538 -1.259 -0.171  2.214  5.613 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    3.834      1.649   2.325   0.0234 *  
# comp_ind       8.541      1.511   5.652 4.45e-07 ***
#   ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 2.615 on 61 degrees of freedom
# Multiple R-squared:  0.3437,	Adjusted R-squared:  0.333 
# F-statistic: 31.95 on 1 and 61 DF,  p-value: 4.449e-07

# checks 
plot(mod_IZsupmeans_m5, which=1)
plot(mod_IZsupmeans_m5, which=2)
shapiro.test(mod_IZsupmeans_m5$residuals)
hist(disupmeans_m5$IZ_diff)
hist(disupmeans_m5$comp_ind)


# filter for inihibition zones less than 5 mm but more than 0
disupmeans_l5 <- di_sup_means[di_sup_means$IZ_group == 'less than 5 mm', ]

# carry out regression
mod_IZsupmeans_l5 <- lm(IZ_diff ~ comp_ind, data = disupmeans_l5)
summary(mod_IZsupmeans_l5)

# Call:
#   lm(formula = IZ_diff ~ comp_ind, data = disupmeans_l5)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.57724 -0.46843 -0.01654  0.43830  1.77202 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)   -0.817      1.131  -0.722   0.4763  
# comp_ind       2.751      1.086   2.532   0.0175 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.8161 on 27 degrees of freedom
# Multiple R-squared:  0.1919,	Adjusted R-squared:  0.162 
# F-statistic: 6.411 on 1 and 27 DF,  p-value: 0.01747


# checks 
plot(mod_IZsupmeans_l5, which=1)
plot(mod_IZsupmeans_l5, which=2)
shapiro.test(mod_IZsupmeans_l5$residuals)
hist(disupmeans_l5$IZ_diff)
hist(disupmeans_l5$comp_ind)
