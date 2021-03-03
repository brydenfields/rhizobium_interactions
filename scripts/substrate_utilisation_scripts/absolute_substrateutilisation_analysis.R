# Author: Bryden Fields
# Email: bryden.fields@york.ac.uk
# Last updated: 02/02/2021

# Initial setup -------------------------------------------------------------------

setwd('/Users/brydenfields/Documents/Publications/2021_Rhizobiuminteractions_paper/scripts/supernatant_scripts')

library(tidyverse)
library(RColorBrewer)
library(plyr); library(dplyr)
library(rcompanion)
library(readxl)
library(vegan)
library(openxlsx)
library(ggpubr)
library(factoextra)
library(ggfortify)
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
library(pheatmap)
library(FSA)
library(reshape2)

source('../summarySE.R')


# Reading and formatting data ---------------

#reading excel and formatting containing emmas growth and ecoplate data for 24 supernatant experiment strains
data <- read_excel("../../data/raw_data/substrate_utilisation_data/24sup_growth_phenotype_formatted.xlsx")
data$geno <- factor(data$geno, levels = c("CC","OA","OC","OE"))

data <- data %>% 
  filter(resource %in% c("Tween 40", "Tween 80", "a-Cyclodextrin", "Glycogen", "D-Cellobiose",
                          "a-D-Lactose", "B-Methyl-D-Glucoside", "D-Xylose", "i-Erythrithol", "D-Mannitol",
                          "N-Acetyl-D-Glucosamine", "Glucose-1-Phosphate", "D,L-a-Glycerol Phosphate", "Pyruvic Acid Methl Ester",
                          "D-Glucosaminic Acid", "D-Galactonic Acid y-Lactone", "D-Galacturonic Acid", "2-Hydroxy-Benzoic Acid",
                          "4-Hydroxy-Benzoic Acid", "y-Hydroxybutyric Acid", "Itaconic Acid", "a-Ketobutyric Acid", "D-Malic Acid",
                          "L-Arginine", "L-Asparagine", "L-Phenylalanine", "L-Serine", "L-Threonine", "Glycyl-L-Glytamic Acid", 
                          "Phenylethyl-amine", "Putrescine"))


#data for correlating AWCD against inoculant/supernatant RGI
RGI_AWCD <- read.csv('../../data/raw_data/substrate_utilisation_data/compind_sup_inoc_AWCD.csv')




# Calculate descriptive statistics -------------

# calculate mean OD for each strain by resource group
data_strainmean <-  summarySE(data, measurevar = 'OD', groupvars = c('strain', 'geno', 'resource_group'))
write.csv(data_strainmean, file = '../../data/intermediate_data/substrate_utilisation_data/resourcegroup_strainmeanOD_absolute.csv')

# calculate mean OD for each strain by resource characteristics
data_strainmean2 <-  summarySE(data, measurevar = 'OD', groupvars = c('strain', 'geno', 'resource_characteristics'))
write.csv(data_strainmean2, file = '../../data/intermediate_data/substrate_utilisation_data/resourcecharacteristics_strainmeanOD_absolute.csv')

# calculate mean OD for each genospecies by resource
data_means <- summarySE(data, measurevar = 'OD', groupvars = c('geno', 'resource_group', 'resource'))

# calculate mean OD for each genospecies by resource group
data_means1 <- summarySE(data, measurevar = 'OD', groupvars = c('geno', 'resource_group'))

# format resource name labels so they can include their special characters
data$resource <- factor(data$resource, levels = c("Tween 40", "Tween 80", "a-Cyclodextrin", "Glycogen", "D-Cellobiose",
                                                 "a-D-Lactose", "B-Methyl-D-Glucoside", "D-Xylose", "i-Erythrithol", "D-Mannitol",
                                                 "N-Acetyl-D-Glucosamine", "Glucose-1-Phosphate", "D,L-a-Glycerol Phosphate", "Pyruvic Acid Methl Ester",
                                                 "D-Glucosaminic Acid", "D-Galactonic Acid y-Lactone", "D-Galacturonic Acid", "2-Hydroxy-Benzoic Acid",
                                                 "4-Hydroxy-Benzoic Acid", "y-Hydroxybutyric Acid", "Itaconic Acid", "a-Ketobutyric Acid", "D-Malic Acid",
                                                 "L-Arginine", "L-Asparagine", "L-Phenylalanine", "L-Serine", "L-Threonine", "Glycyl-L-Glytamic Acid", 
                                                 "Phenylethyl-amine", "Putrescine"),
                                        labels = c("Tween 40", "Tween 80", "\U03B1-Cyclodextrin", "Glycogen", "D-Cellobiose",
                                                   "\U03B1-D-Lactose", "\U03B2-Methyl-D-Glucoside", "D-Xylose", "i-Erythrithol", "D-Mannitol",
                                                   "N-Acetyl-D-Glucosamine", "Glucose-1-Phosphate", "D,L-\U03B1-Glycerol Phosphate", "Pyruvic Acid Methl Ester",
                                                   "D-Glucosaminic Acid", "D-Galactonic Acid \U03B3-Lactone", "D-Galacturonic Acid", "2-Hydroxy-Benzoic Acid",
                                                   "4-Hydroxy-Benzoic Acid", "\U03B3-Hydroxybutyric Acid", "Itaconic Acid", "\U03B1-Ketobutyric Acid", "D-Malic Acid",
                                                   "L-Arginine", "L-Asparagine", "L-Phenylalanine", "L-Serine", "L-Threonine", "Glycyl-L-Glutamic Acid", 
                                                   "Phenylethyl-amine", "Putrescine"))




# Graphs: genospecies resource utilisation and comparison to RGI -------------------------------------------------------------------

#creating cut off
cutoff <- data.frame(yintercept=0, cutoff=factor(0))

# genospecies colours
genoPalette <- c(CC = "#9ACD32", OA = "#466EA9", OC = "#4D9A7A", OE ="#DA367E")


# graph all genospecies resources comparison.

(allgeno_comparison <- ggplot(data, aes(x = resource, y = OD, colour = geno), group = geno)+
    geom_boxplot( stat = "boxplot",
                  position = position_dodge(1), 
                  size = .5)+
    ggtitle("Genospecies Resource utilisation")+
    labs(x = "Resource", y = "Optical Density (OD600)")+
    scale_y_continuous(breaks=seq(0,4.5,by=0.5))+
    scale_colour_manual(values=genoPalette)+
    geom_hline(aes(yintercept=yintercept, linetype=cutoff), data=cutoff, show.legend=FALSE)+
    theme_bw()+
    theme(legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust=1), 
          text = element_text(size=6), legend.text = element_text(size=6)))
#function to save specific ggplot graph.  additional parameters(width=,height=,units='cm')
ggsave('../../data/intermediate_data/substrate_utilisation_data/allgeno_resourcecomparison.pdf', plot = allgeno_comparison, width = 25, height = 10, units = 'cm')


#graph for AWCD - identify strains that are more generalist
awcd_no159 <- data %>%
  filter(strain %in% c("41", "53", "57", "67", "74", "77", 
                       "137B","144A","145B","152A", "152B", "154C",
                       "122A","126A","147A","157B","158", "165A", "170C",
                       "135A", "135B", "149A", "168A"))

# relevel strain order for graphs
awcd_no159$strain <- factor(awcd_no159$strain, levels = c("41", "53", "57", "67",
                                                    "74", "77", "137B", "144A", "145B", "152A",
                                                    "152B", "154C", "122A", "126A", 
                                                    "147A", "157B", "158", "165A",
                                                    "170C", "135A", "135B", "149A", 
                                                    "159", "168A"))

# add SM to the strain names
awcd_no159$strain <- paste0("SM", awcd_no159$strain)


#graph correlating AWCD against inoculant/supernatant RGI
RGI_AWCD <- RGI_AWCD %>%
  filter(strain %in% c("41", "53", "57", "67", "74", "77", 
                       "137B","144A","145B","152A", "152B", "154C",
                       "122A","126A","147A","157B","158", "165A", "170C",
                       "135A", "135B", "149A", "168A"))


(AWCDgraph <- ggplot(RGI_AWCD, aes(x = comp_ind_inoc, y = AWCD, colour = geno, size = 0.25), group = geno)+
    geom_point()+
    labs(x = "Mean RGI of inoculant strain", 
         #y = expression(paste("Metabolic capacity (",OD[590],")",sep=""))
         y = "Metabolic capacity")+
    scale_colour_manual(values=genoPalette)+
    scale_size(guide = "none")+
    geom_smooth(method='lm', col = "grey60", se=FALSE, aes(size = 0.2 ))+
    xlim(0.5,1.5)+
    guides(colour = guide_legend(override.aes = list(size=8))) +
    theme_bw()+
    theme(legend.title = element_blank(), 
          #axis.text.x = element_text(angle = 45, hjust = 1, vjust=1), 
          text = element_text(size=20)))
#function to save specific ggplot graph.  additional parameters(width=,height=,units='cm')
ggsave('../../data/intermediate_data/substrate_utilisation_data/AWCD_RGI_correlation.pdf', plot = AWCDgraph, width = 15, height = 15, unit = 'cm')

# correlation coefficient
cor.test(RGI_AWCD$comp_ind_inoc, RGI_AWCD$AWCD)

# data:  RGI_AWCD$comp_ind_inoc and RGI_AWCD$AWCD
# t = 2.9644, df = 20, p-value = 0.007665
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.1706506 0.7900745
# sample estimates:
#   cor 
# 0.5525074 




# PCA individual resources -------------------------------------------------------------------

# using the following protocol: https://www.r-bloggers.com/computing-and-visualizing-pca-in-r/


#upload data for PCA
#individual resources
pca_data <- read_excel("../../data/raw_data/substrate_utilisation_data/24sup_growth_phenotype_AWCDgroups_PCAdata.xlsx")

# remove 159 because outlier - ecoplate results are likely contaminated.
pca_data <- pca_data[!pca_data$strain == "159", ]

#add 1 to all the AWCD scores because otherwise will not be able to log some of the negative scores
pca_data[ ,6:37] <-  pca_data[ ,6:37] + 1

#set strain names to be row names
rownames(pca_data) <- pca_data$strain

#Since skewness and the magnitude of the variables influence the resulting PCs, it is good practice to apply skewness transformation,
#center and scale the variables prior to the application of PCA.

# log transform continuous variables (AWCD data for carbon groups). will use categorical data later to visualise tables.
pca_datalog <- log(pca_data[, 7:37])
rownames(pca_datalog) <- pca_data$strain
#pca_datalog <- log(pca_data[, 17:30])
#rownames(pca_datalog) <- pca_data$strain


# apply PCA - scale. = TRUE is highly advisable, but default is FALSE. 
eco_pca <- prcomp(pca_datalog,
                 center = TRUE,
                 scale. = TRUE)

#The prcomp function returns an object of class prcomp, which have some methods available. 
#The print method returns the standard deviation of each of the four PCs, and their rotation (or loadings), 
#which are the coefficients of the linear combinations of the continuous variables.

# print method
print(eco_pca)

# The plot method returns a plot of the variances (y-axis) associated with the PCs (x-axis). 
# The Figure below is useful to decide how many PCs to retain for further analysis. 
# See which PCs explain most of the variability in the data.
# In this case, it seems like the first 5 PCs explain most of the variability within the data.
# plot method
plot(eco_pca, type = "l")

# The summary method describe the importance of the PCs. 
# The first row describe again the standard deviation associated with each PC. 
# The second row shows the proportion of the variance in the data explained by each component.
# The third row describe the cumulative proportion of explained variance. 
# summary method
summary(eco_pca)


#access to PCA results
# Eigenvalues
eig.val <- get_eigenvalue(eco_pca)
eig.val

# Results for Variables
res.var <- get_pca_var(eco_pca)
res.var$coord          # Coordinates
var_contrib <- res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 
# Results for individuals
res.ind <- get_pca_ind(eco_pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 




#Graph of variables. 
#Positive correlated variables point to the same side of the plot. Negative correlated variables point to opposite sides of the graph.
(p <- fviz_pca_var(eco_pca,
                   col.var = "contrib", # Color by contributions to the PC
                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                   legend.title = "Contribution\n of variable",
                   repel = TRUE     # Avoid text overlapping
))
#change dim axes titles to pc titles with correct % value.
p + xlab("PC1 variance (38.9%)") +
  ylab("PC2 variance (23.6%)") +
  theme_bw()
ggsave('../../data/intermediate_data/substrate_utilisation_data/individualsubstrates_pca_PC1P2_159rem_varcontrib.pdf', plot = p, width = 25, height = 20, unit = 'cm')


# genospecies colours palette
genoPalette <- c(CC = "#9ACD32", OA = "#466EA9", OC = "#4D9A7A", OE ="#DA367E")



#get co-ordinates from prcomp
#add geno variables etc.
pca_coord <- as.data.frame(eco_pca$x)
pca_coord



# calculate the contribution to the total variance for each component
percentVar <- eco_pca$sdev^2 / sum( eco_pca$sdev^2 )
percentVar <- round(100 * percentVar)


# pca graph
(pcaplot <- ggplot(eco_pca, aes(PC1, PC2, color = pca_data$geno, fill = pca_data$geno)) +
    coord_fixed() +
    geom_point(size=3) +
    #must manually apply the percentage variance in this case
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    stat_chull(geom = "polygon", alpha = 0.2) +
    stat_mean(size = 5) +
    scale_color_manual(values = genoPalette) +
    scale_fill_manual(values = genoPalette) +
    theme_bw() +
    theme(legend.title = element_blank(), 
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title=element_text(size = 20),
          legend.text = element_text(size = 16))
)
#function to save specific ggplot graph.  additional parameters(width=,height=,units='cm')
ggsave('../../data/intermediate_data/substrate_utilisation_data/individualsubstrates_pca_PC1P2_159rem.pdf', plot = pcaplot, width = 25, height = 20, unit = 'cm')






# PCA grouped resources -------------------------------------------------------------------

# using the following protocol: https://www.r-bloggers.com/computing-and-visualizing-pca-in-r/


#upload data for PCA
#grouped resources
#Use AWCD data rather than relative OD values averaged (mean).
#I have also removed some data such as soil type and biofilm data. 
pca_data_groupsub <- read_excel("../../data/raw_data/substrate_utilisation_data/24sup_growth_phenotype_AWCDgroups_PCAdata.xlsx")

# remove 159 because outlier?
pca_data_groupsub <- pca_data_groupsub[!pca_data_groupsub$strain == "159", ]

#add 1 to all the AWCD scores because otherwise will not be able to log some of the negative scores
pca_data_groupsub[ ,39:44] <-  pca_data_groupsub[ ,39:44] + 1

#set strain names to be row names
rownames(pca_data_groupsub) <- pca_data_groupsub$strain

# log transform continuous variables (AWCD data). will use categorical data later to visualise tables.
pca_data_groupsub_log <- log(pca_data_groupsub[, 39:44])
rownames(pca_data_groupsub_log) <- pca_data_groupsub$strain

# rename columns
names(pca_data_groupsub_log)[c(2,4:6)] <- c('Amino acids', 'Carboxylic acids', 
                                  'Complex carbon sources', 'Phosphate carbon')

#taking the genospecies/descriptive data for use later
pca_data_groupsub_geno <- pca_data_groupsub[, 2]
rownames(pca_data_groupsub_geno) <- pca_data_groupsub$strain

pca_data_groupsub_descipt <- pca_data_groupsub[, 2:5]
rownames(pca_data_groupsub_descipt) <- pca_data_groupsub$strain

# apply PCA - scale. = TRUE is highly advisable, but default is FALSE. 
eco_pca_groupsub <- prcomp(pca_data_groupsub_log,
                  center = TRUE,
                  scale. = TRUE)

#The prcomp function returns an object of class prcomp, which have some methods available. 
#The print method returns the standard deviation of each of the four PCs, and their rotation (or loadings), 
#which are the coefficients of the linear combinations of the continuous variables.

# print method
print(eco_pca_groupsub)

#The plot method returns a plot of the variances (y-axis) associated with the PCs (x-axis). 
#The Figure below is useful to decide how many PCs to retain for further analysis. 
#See which PCs explain most of the variability in the data.
# plot method
plot(eco_pca_groupsub, type = "l")

#The summary method describe the importance of the PCs. 
#The first row describe again the standard deviation associated with each PC. 
#The second row shows the proportion of the variance in the data explained by each component.
#The third row describe the cumulative proportion of explained variance. 
# summary method
summary(eco_pca_groupsub)

#access to PCA results
# Eigenvalues
eig.val_groupsub <- get_eigenvalue(eco_pca_groupsub)
eig.val_groupsub

# Results for Variables
res.var_groupsub <- get_pca_var(eco_pca_groupsub)
res.var_groupsub$coord          # Coordinates
res.var_groupsub$contrib        # Contributions to the PCs
res.var_groupsub$cos2           # Quality of representation 
# Results for individuals
res.ind_groupsub <- get_pca_ind(eco_pca_groupsub)
res.ind_groupsub$coord          # Coordinates
res.ind_groupsub$contrib        # Contributions to the PCs
res.ind_groupsub$cos2           # Quality of representation 




#Graph of variables. 
#Positive correlated variables point to the same side of the plot. Negative correlated variables point to opposite sides of the graph.
(p_groupsub <- fviz_pca_var(eco_pca_groupsub,
                   col.var = "contrib", # Color by contributions to the PC
                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                   legend.title = "Contribution\nof\nvariable",
                   labelsize = 7,
                   repel = TRUE     # Avoid text overlapping
))
#change dim axes titles to pc titles with correct % value.
p_groupsub + xlab("PC1 variance (52.1%)") +
  ylab("PC2 variance (23.9%)") +
  theme_bw() +
  theme(legend.title = element_text(size = 20),
        text = element_text(size = 25),
        axis.title = element_text(size = 25),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
ggsave('../../data/intermediate_data/substrate_utilisation_data/groupedsubstrates_pca_PC1P2_159rem_varcontrib.pdf', plot = p_groupsub, width = 20, height = 15, unit = 'cm')


# genospecies colour palette
genoPalette <- c(CC = "#9ACD32", OA = "#466EA9", OC = "#4D9A7A", OE ="#DA367E")


#get co-ordinates from prcomp
#add geno variables etc.
pca_coord_groupsub <- as.data.frame(eco_pca_groupsub$x)


# the contribution to the total variance for each component
percentVar_groupsub <- eco_pca_groupsub$sdev^2 / sum( eco_pca_groupsub$sdev^2 )
percentVar_groupsub <- round(100 * percentVar_groupsub)


# pca graph
(pcaplot_groupsub <- ggplot(eco_pca_groupsub, aes(PC1, PC2, color = pca_data_groupsub$geno, fill = pca_data_groupsub$geno)) +
    coord_fixed() +
    geom_point(size=3) +
    #must manually apply the percentage variance in this case
    xlab(paste0("PC1: ",percentVar_groupsub[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar_groupsub[2],"% variance")) + 
    stat_chull(geom = "polygon", alpha = 0.2) +
    stat_mean(size = 5) +
    scale_color_manual(values = genoPalette) +
    scale_fill_manual(values = genoPalette) +
    theme_bw() +
    theme(legend.title = element_blank(), 
          text = element_text(size = 30))
    )
#function to save specific ggplot graph.  additional parameters(width=,height=,units='cm')
ggsave('../../data/intermediate_data/substrate_utilisation_data/groupedsubstrates_pca_PC1P2_159rem.pdf', plot = pcaplot_groupsub, width = 20, height = 15, unit = 'cm')


# PERMANOVA analysis for confirmation of PCA groups
# turn data into euclidean distance with vegan package vegdist
permanova_data <- vegdist(pca_data_groupsub_log, method = "euclidean")
# now do permanova with adonis
adonis(permanova_data ~ geno, data = pca_data_groupsub[2])

# Call:adonis(formula = permanova_data ~ geno, data = pca_data[2]) 
# Permutation: free
# Number of permutations: 999
# Terms added sequentially (first to last)
#             Df SumsOfSqs  MeanSqs F.Model     R2 Pr(>F)   
# geno       3   0.12932 0.043108  3.8293 0.3768  0.006 **
# Residuals 19   0.21389 0.011257         0.6232          
# Total     22   0.34321                  1.0000          


#POST-HOC TEST for Adonis
#see which plots are significant for adonis?
pairwise.adonis(permanova_data, pca_data_groupsub$geno, p.adjust.m = 'bonferroni')

#     pairs  Df   SumsOfSqs   F.Model         R2 p.value p.adjusted sig
# 1 CC vs OC  1 0.015661763 1.6382824 0.12962856   0.175      1.000    
# 2 CC vs OE  1 0.108988651 9.6210099 0.54599651   0.005      0.030   # significant
# 3 CC vs OA  1 0.004962939 0.5651670 0.05349343   0.521      1.000    
# 4 OC vs OE  1 0.058457361 4.1729968 0.31678417   0.021      0.126    
# 5 OC vs OA  1 0.007094862 0.6331383 0.05442541   0.644      1.000    
# 6 OE vs OA  1 0.088615836 6.5199759 0.44903490   0.015      0.090    




# Does metabolic similarity of grouped substrates correlate with ANI? -----------------

# using data from above section: 'pca_data_groupsub_log'

# turn data into euclidean distance with vegan package vegdist but makes it symmetrical so can do mantel test. 
permanova_data2 <- as.matrix(dist(pca_data_groupsub_log, method = "euclidean"))

# read in 6K ANI data 
ANI_data <- read.csv('../../data/raw_data/supernatant_data/ani_sorted_by_genospecies_snps_new_data_6kgenes_supernatantsamples.csv', row.names = 1)

# remove SM from strain names in rows and column names
row.names(ANI_data) <- gsub('SM', '', row.names(ANI_data))
names(ANI_data) <- gsub('SM', '', names(ANI_data))


# remove 159 because outlier?
ANI_data <- ANI_data[!row.names(ANI_data) == "159", !names(ANI_data) == "159"]
# reorder matrix columns and rows to be the same as metabolic similarlity matrix
ANI_data <- ANI_data[row.names(pca_data_groupsub), row.names(pca_data_groupsub)]
# make into a dist object
ANI_data <- as.matrix(ANI_data)

# do a mantel test to correlate metabolic similarity (euclidean distance) to 6K ANI
# all written arguments are default i.e. 9999 permutations

mantel(permanova_data2, ANI_data, method = "pearson")
# Mantel statistic based on Pearson's product-moment correlation 
# Call:mantel(xdis = permanova_data, ydis = ANI_data, method = "pearson") 
# Mantel statistic r: -0.231 
#       Significance: 0.987 
# Upper quantiles of permutations (null model):
#   90%   95% 97.5%   99% 
# 0.107 0.130 0.139 0.150 
# Permutation: free
# Number of permutations: 999


# the two graphs
plot(permanova_data2, ANI_data) # AWCD with log



# heatmap for presence/absence of metabolic capacity ------------------

# read csv
data_heat <- read_excel("../../data/raw_data/substrate_utilisation_data/24sup_growth_phenotype_AWCDgroups_PCAdata.xlsx")

# remove 159 because outlier?
data_heat <- data_heat[!data_heat$strain == "159", ]

#set strain names to be row names
rownames(data_heat) <- data_heat$strain

# Order by genospecies
data_heat <- data_heat[order(data_heat$geno),]

# split data into two dataframes for metadata and OD data
data_heat_OD <- data_heat[, 7:37]
rownames(data_heat_OD) <- data_heat$strain


# make a dataframe for the genospecies label colours
datTraits_heat <- data_heat[1:2]
row.names(datTraits_heat) <- datTraits_heat$strain
datTraits_heat$strain <- NULL
datTraits_heat <- as.data.frame(datTraits_heat)
datTraits_heat$geno <- factor(datTraits_heat$geno)
names(datTraits_heat)[1] <- "Genospecies"

# reorder columns by most metabolised on average
mns <- colSums(data_heat_OD <= 0)
order(mns)
data_heat_OD <- data_heat_OD[,order(mns)]
row.names(data_heat_OD) <- row.names(datTraits_heat)

# shorten OD readings to 2 decimal places (instead of 4)
data_heat_OD_round <- data_heat_OD %>% 
  mutate_if(is.numeric, round, digits = 2)
row.names(data_heat_OD_round) <- row.names(data_heat_OD)

# column names as labels 
sublabels = c("Glycogen", "D-Xylose", "D-Mannitol", "Pyruvic Acid Methl Ester", 
              "i-Erythrithol", "\U03B3-Hydroxybutyric Acid", "D-Malic Acid", "D-Cellobiose", 
              "N-Acetyl-D-Glucosamine", "Glucose-1-Phosphate", "D-Galactonic Acid \U03B3-Lactone", "L-Serine",
              "L-Threonine", "\U03B1-D-Lactose", "D-Glucosaminic Acid", "D,L-\U03B1-Glycerol Phosphate", 
              "Glycyl-L-Glutamic Acid", "\U03B2-Methyl-D-Glucoside", "\U03B1-Cyclodextrin", "D-Galacturonic Acid",
              "L-Asparagine", "Itaconic Acid", "Phenylethyl-amine", "L-Phenylalanine",
              "Tween 80", "L-Arginine", "Putrescine", "Tween 40",
              "\U03B1-Ketobutyric Acid", "4-Hydroxy-Benzoic Acid", "2-Hydroxy-Benzoic Acid"
              )
names(data_heat_OD_round) <- sublabels

# make colour palette for the gene expression
#breaks
bk <- c(-0.5,0.00,seq(0.001,1.10,by=0.01))
# #colors (one less than breaks) colours below 0 get red, colours above 0 get a gradient blue
mycols <- c("red", colorRampPalette(colors = c("white","blue"))(length(bk)-2))


#colour palette for genospecies labels
AnnColour = list(
  Genospecies = c(CC = "#9ACD32", OA = "#466EA9", OC = "#4D9A7A", OE ="#DA367E")
)


# make heatmap with pretty heatmap
pheatmap(data_heat_OD_round,
         display_numbers = T, number_format = "%.2f",
         cluster_rows = F, cluster_cols = F,
         gaps_row = c(6,12,19),
         color = mycols,
         breaks = bk,
         annotation_row = datTraits_heat,
         annotation_colors = AnnColour,
         cellwidth = 20, cellheight = 20,
         width = 14, height = 14, 
         fontsize = 7,
         fontsize_col = 10, fontsize_row = 10,
         legend=T,
         filename = "../../data/intermediate_data/substrate_utilisation_data/substrateutilisation_heatmap.pdf")




# Statistical analyses -------------------------------------------------------------------



# 1) Differences in metabolic capacity (AWCD) between genospecies

# normal quantile plot (normal Q–Q plot).
# without SM159
AWCD_1 <- awcd_no159$AWCD + 1
plotNormalHistogram(AWCD_1, breaks = 20)
qqnorm(AWCD_1, ylab="Sample Quantiles for data$AWCD")
qqline(AWCD_1, col="red")

# normality test with SM159
shapiro.test(data$AWCD)


# non-parametric one way anova
kruskal.test(AWCD~geno, data[1:24,])

# Kruskal-Wallis rank sum test
# 
# data:  AWCD by geno
# Kruskal-Wallis chi-squared = 11.152, df = 3, p-value = 0.01093

# Dunns post hoc test.
PT2 = dunnTest(AWCD ~ geno,
              data=data[1:24,],
              method="bonferroni")
PT2
# Comparison          Z     P.unadj      P.adj
# 1    CC - OA -0.4490731 0.653378911 1.00000000
# 2    CC - OC -1.2407164 0.214710518 1.00000000
# 3    OA - OC -0.7746912 0.438522127 1.00000000
# 4    CC - OE -3.1217808 0.001797607 0.01078564 # sig
# 5    OA - OE -2.6936063 0.007068359 0.04241016 # sig
# 6    OC - OE -2.0494947 0.040413764 0.24248259





# 2) Number of substrates metabolised vs genospecies groups

# Kruskal-Wallis non-parametric one way anova
num_substrates_data <- read.csv('../../data/raw_data/substrate_utilisation_data/num_substrates_used_ecoplate.csv')

# looking at data distribution
num_substrates_data2 <- num_substrates_data$num_substrates
plotNormalHistogram(num_substrates_data2, breaks = 10)
qqnorm(num_substrates_data2, ylab="Sample Quantiles for num_substrates_data$num_substrates")
qqline(num_substrates_data2, col="red")


# or non-parametric one way anova
kruskal.test(num_substrates ~ genospecies, num_substrates_data)

# Kruskal-Wallis rank sum test
# 
# data:  num_substrates by genospecies
# Kruskal-Wallis chi-squared = 7.9402, df = 3, p-value = 0.04726

# Dunns test for kruskal wallis
PT3 = dunnTest(num_substrates ~ genospecies,
              data = num_substrates_data,
              method = "bonferroni")
PT3

# Dunn (1964) Kruskal-Wallis multiple comparison
# p-values adjusted with the Bonferroni method.
# 
# Comparison          Z     P.unadj      P.adj
# 1    CC - OA  1.4514452 0.146655931 0.87993558
# 2    CC - OC -0.2436558 0.807497406 1.00000000
# 3    OA - OC -1.7498917 0.080137008 0.48082205
# 4    CC - OE -1.4700342 0.141552481 0.84931489
# 5    OA - OE -2.7682463 0.005635885 0.03381531 # sig
# 6    OC - OE -1.2976504 0.194407487 1.00000000



# Correlate no. substrates metabolised and RGI? Graph: metabolic capacity and no. substrates metabolised -----------------

# Read in file with RGI as inoculant and AWCD
RGI_AWCD <- read.csv('../../data/raw_data/substrate_utilisation_data/compind_sup_inoc_AWCD.csv')
RGI_AWCD <- RGI_AWCD %>%
  filter(strain %in% c("41", "53", "57", "67", "74", "77", 
                       "137B","144A","145B","152A", "152B", "154C",
                       "122A","126A","147A","157B","158", "165A", "170C",
                       "135A", "135B", "149A", "168A"))

# add SM to the beginning of strain names
RGI_AWCD$strain <- paste0("SM", RGI_AWCD$strain)

# merge with num_substrate_data
RGI_AWCD_numsub <- merge(RGI_AWCD, num_substrates_data, by = "strain", all.x = T)


# genospecies colour palette
genoPalette <- c(CC = "#9ACD32", OA = "#466EA9", OC = "#4D9A7A", OE ="#DA367E")

# graph correlating number of substrates able to use and RGI
(numsubgraph <- ggplot(RGI_AWCD_numsub, aes(x = comp_ind_inoc, y = num_substrates, colour = geno, size = 0.25), group = geno)+
   geom_point()+
   labs(x = "Mean RGI of inoculant strain", y = "Number of substrates can metabolised") +
   scale_colour_manual(values=genoPalette)+
   scale_size(guide = "none")+
   geom_smooth(method='lm', col = "grey60", se=FALSE, aes(size = 0.2 ))+
   xlim(0.5,1.5)+
   guides(colour = guide_legend(override.aes = list(size=10))) +
   theme_bw()+
   theme(legend.title = element_blank(), 
         #axis.text.x = element_text(angle = 45, hjust = 1, vjust=1), 
         text = element_text(size=18)))
#function to save specific ggplot graph.  additional parameters(width=,height=,units='cm')
ggsave('../../data/intermediate_data/substrate_utilisation_data/numsubstrates_RGI_correlation.pdf', plot = numsubgraph, width = 20, height = 20, unit = 'cm')

cor.test(RGI_AWCD_numsub$comp_ind_inoc, RGI_AWCD_numsub$num_substrates)

# Pearson's product-moment correlation
# 
# data:  RGI_AWCD_numsub$comp_ind_inoc and RGI_AWCD_numsub$num_substrates
# t = 2.3683, df = 20, p-value = 0.02806
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.05778116 0.74299788
# sample estimates:
#       cor 
# 0.4679889 


# cast to make num_sub and AWCD into one column
RGI_AWCD_numsub_melt <- melt(RGI_AWCD_numsub[, c("strain", "geno", "AWCD", "num_substrates")])
names(RGI_AWCD_numsub_melt)[3] <- c("AWCD_substrate")

# relevel AWCD and num substrates
RGI_AWCD_numsub_melt$AWCD_substrate <- factor(RGI_AWCD_numsub_melt$AWCD_substrate,
                                              level = c('num_substrates', 'AWCD'))

# facet to make plot together with number of substrates metabolised.
(AWCDgenograph4 <- ggplot(RGI_AWCD_numsub_melt, aes(x = geno, y = value, colour = geno, group = geno))+
    geom_boxplot(aes(fill = geno), alpha = 0.5)+
    geom_jitter(col = "black", width = 0.1, size = 2) +
    labs(x = "genospecies group")+
    facet_grid(AWCD_substrate ~ ., scales = "free",
               labeller = as_labeller(c("AWCD" = "Metabolic capacity", 
                                        "num_substrates" = "Number of substrates\nmetabolised"))) +
    scale_colour_manual(values=genoPalette)+
    scale_fill_manual(values=genoPalette)+
    #scale_y_continuous(breaks=seq(0,1.20,by=0.1))+
    theme_bw()+
    theme(axis.title.y = element_blank(),
          legend.title = element_blank(), 
          legend.position="none", 
          text = element_text(size=22)))
#function to save specific ggplot graph.  additional parameters(width=,height=,units='cm')
ggsave('../../data/intermediate_data/substrate_utilisation_data/numsubstrates_AWCD_geno_facet.pdf', plot = AWCDgenograph4, width = 14.5, height = 15, unit = 'cm')
