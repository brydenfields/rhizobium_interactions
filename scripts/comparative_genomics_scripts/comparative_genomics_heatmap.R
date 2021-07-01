# Author: Bryden Fields
# Email: bryden.fields@york.ac.uk
# Last updated: 01/07/2021


# initial set-up -------------------

library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(rcompanion)

#set the working directory to the directory this script is in.
#setwd("~/what/ever/folder/you/are/working/from") 

# percentage identity heatmap --------------

# read in percentage identity table and fill blanks with NA
data <- read.csv('../../data/raw_data/comparative_genomics_data/quorumsensing_percentageidentity_blasttable.csv', na.strings=c("","NA"))

# remove the % signs
data <- as.data.frame(sapply(data, gsub, pattern='%', replacement=''))

# make gene names as row names 
row.names(data) <- data$X
data$X <- NULL

# remove X from column names
names(data) <- sub("X", "SM", names(data))

# reverse the column names
data=data[,order(ncol(data):1)]

# make a dataframe for the genospecies label colours
datTraits <- data.frame( Strain = names(data), 
                         Genospecies = c("OA", "OA", "OA", "OA", "OA", "OA",
                                         "OC", "OC", "OC", "OC", "OC", "OC", "OC", 
                                         "OE", "OE", "OE", "OE", "OE",
                                         "CC", "CC", "CC", "CC", "CC", "CC"),
                         row.names = names(data))

datTraits$Strain <- NULL

# make heatmap numeric
data[1:24] <- lapply(data[1:24], function(x) as.numeric(as.character(x)))

# rename rows
row.names(data)
row.names(data) <- c("tfxG", "medium bacteriocin", "cinI", "cinR", "cinS",              
                     "expR", "bisR", "traI", "traR", "raiI",              
                     "raiR", "rhiI", "rhiR", "rhiA", "rhiB", "rhiC")


# make colour palette for the gene percentage identity
my_palette <- colorRampPalette(brewer.pal(9,"RdPu"))(256)

#colour palette for genospecies labels
AnnColour = list(
  Genospecies = c(CC = "#9ACD32", OA = "#466EA9", OC = "#4D9A7A", OE ="#DA367E"))


# make heatmap with pretty heatmap
pheatmap(data,
         cluster_rows = F, cluster_cols = F,
         gaps_col = c(6,13,18),
         col = my_palette,
         annotation_col = datTraits,
         annotation_colors = AnnColour,
         cellwidth = 10, cellheight = 10,
         width = 7, height = 4, 
         fontsize = 7,
         fontsize_col = 10, fontsize_row = 10,
         legend=T,
         filename = "../../data/intermediate_data/comparative_genomics_data/QS_blasthits_percentageidentity_heatmap.pdf")




# Heatmap for secondary metabolites --------------------

### MUST RUN % IDENTITY HEATMAP FIRST

data_sm <- read.csv('../../data/raw_data/comparative_genomics_data/antiSMASH_totalclusters_summary.csv')

# make gene names as row names 
row.names(data_sm) <- data_sm$X
data_sm$X <- NULL

# remove X from column names and add SM
names(data_sm) <- sub("X", "SM", names(data_sm))

# make col names same order as in % identity heatmap
data_sm <- data_sm[,names(data)]


# make colour palette for the gene expression
my_palette_sm <- colorRampPalette(brewer.pal(9,"OrRd"))(6)


# make heatmap with pretty heatmap
pheatmap(data_sm,
         cluster_rows = F, cluster_cols = F,
         gaps_col = c(6,13,18),
         col = my_palette_sm,
         annotation_col = datTraits,
         annotation_colors = AnnColour,
         cellwidth = 10, cellheight = 10,
         width = 7, height = 4, 
         fontsize = 7,
         fontsize_col = 10, fontsize_row = 10,
         legend=T,
         filename = "../../data/intermediate_data/comparative_genomics_data/secondarymetabolite_count_heatmap.pdf")




# Heatmap for prophages ---------------------------------


### MUST RUN % IDENTITY HEATMAP FIRST

data_pha <- read.csv('../../data/raw_data/comparative_genomics_data/phagecount_summary.csv')

# add SM to strain names
data_pha$X <- paste0("SM", data_pha$X)

# make gene names as row names 
row.names(data_pha) <- data_pha$X
data_pha$X <- NULL

# transpose dataframe
data_pha <- data.frame(t(data_pha))

# make col names same order as in % identity heatmap
data_pha <- data_pha[,names(data)]

# make colour palette for the gene expression
my_palette_pha <- colorRampPalette(brewer.pal(9,"GnBu"))(6)


# make heatmap with pretty heatmap
pheatmap(data_pha,
         cluster_rows = F, cluster_cols = F,
         gaps_col = c(6,13,18),
         col = my_palette_pha,
         annotation_col = datTraits,
         annotation_colors = AnnColour,
         cellwidth = 10, cellheight = 10,
         width = 7, height = 4, 
         fontsize = 7,
         fontsize_col = 10, fontsize_row = 10,
         legend=T,
         filename = "../../data/intermediate_data/comparative_genomics_data/phages_count_heatmap.pdf")


# Heatmap for 196 Rlt strain presence/absence percentage analysis -------------


# read in data 
data_196 <- read.csv('../../data/raw_data/comparative_genomics_data/comparative_genomics_data/genospecies_196strain_presabs_percentage.csv')

# make gene names as row names 
row.names(data_196) <- data_196$gene
data_196$gene <- NULL

# change 0 to NA
data_196[data_196 == 0] <- NA

# colour palette
# make colour palette for the gene presence/absence and colour
my_palette_196 <- colorRampPalette(brewer.pal(9,"YlGn"))(256)


# make heatmap with pretty heatmap
pheatmap(data_196,
         cluster_rows = F, cluster_cols = F,
         col = my_palette_196,
         cellwidth = 10, cellheight = 10,
         width = 7, height = 4, 
         fontsize = 7,
         fontsize_col = 10, fontsize_row = 10,
         legend=T,
         filename = "../../data/intermediate_data/comparative_genomics_data/genospecies196_presabs_percentage_heatmap.pdf")





# chi-square for significance of non-random quorum sensing gene distributions ------------------


# presence and absence of whole gene qs system for each genospecies.
# do individual chi-square tests for each QS pathway (tra, rhi, rai)
whole_contindata <- read.csv('/../../data/raw_data/comparative_genomics_data/genospecies_196strain_presabs_absolutevalues_v2.csv')


# use fisher's exact test for each gene instead.
#http://rcompanion.org/rcompanion/b_07.html


for (i in unique(whole_contindata$gene)){
  # make a contingency table for each gene
  subset_data <- whole_contindata[whole_contindata$gene == i, ]
  row.names(subset_data) <- subset_data$genospecies
  subset_data[c("genospecies", "gene")] <- NULL
  # make into a matrix for the test
  subset_data <- as.matrix(subset_data, header=TRUE)
  # run the fisher exact test
  print(i)
  print(fisher.test(subset_data,
              alternative="two.sided"))
  
}


# [1] "bisR"
# 
# Fisher's Exact Test for Count Data
# 
# data:  subset_data
# p-value = 0.0004895
# alternative hypothesis: two.sided
# 
# [1] "traI"
# 
# 	Fisher's Exact Test for Count Data
# 
# data:  subset_data
# p-value = 0.001102
# alternative hypothesis: two.sided
# 
# [1] "traR"
# 
# Fisher's Exact Test for Count Data
# 
# data:  subset_data
# p-value = 0.0004895
# alternative hypothesis: two.sided
# 
# [1] "raiI"
# 
# 	Fisher's Exact Test for Count Data
# 
# data:  subset_data
# p-value < 2.2e-16
# alternative hypothesis: two.sided
# 
# [1] "raiR"
# 
# Fisher's Exact Test for Count Data
# 
# data:  subset_data
# p-value < 2.2e-16
# alternative hypothesis: two.sided
# 
# [1] "rhiI"
# 
# 	Fisher's Exact Test for Count Data
# 
# data:  subset_data
# p-value < 2.2e-16
# alternative hypothesis: two.sided
# 
# [1] "rhiR"
# 
# Fisher's Exact Test for Count Data
# 
# data:  subset_data
# p-value < 2.2e-16
# alternative hypothesis: two.sided
# 
# [1] "rhiA"
# 
# 	Fisher's Exact Test for Count Data
# 
# data:  subset_data
# p-value < 2.2e-16
# alternative hypothesis: two.sided
# 
# [1] "rhiB"
# 
# Fisher's Exact Test for Count Data
# 
# data:  subset_data
# p-value < 2.2e-16
# alternative hypothesis: two.sided
# 
# [1] "rhiC"
# 
# 	Fisher's Exact Test for Count Data
# 
# data:  subset_data
# p-value < 2.2e-16
# alternative hypothesis: two.sided
# 
# [1] "vbsS"
# 
# Fisher's Exact Test for Count Data
# 
# data:  subset_data
# p-value < 2.2e-16
# alternative hypothesis: two.sided






