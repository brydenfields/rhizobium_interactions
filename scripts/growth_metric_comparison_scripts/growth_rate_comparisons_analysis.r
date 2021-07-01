# Author: Bryden Fields
# Email: bryden.fields@york.ac.uk
# Last updated: 01/07/2021


# Initial set-up ---------------------

library(ggplot2)
library(flux)
library(Rmisc)
library(patchwork)
library(reshape2)
library(readxl)

#set the working directory to the directory this script is in.
#setwd("~/what/ever/folder/you/are/working/from") 

source('../summarySE.R')


# Calculating AUC and correlating to 62h OD growth -----------------

# Read in growth data. 
GROdf <- read_excel("../../data/raw_data/supernatant_data/bes_project_all.xlsx")

# remove comp_ind and comp_ind_50TY columns
GROdf$comp_ind <- NULL
GROdf$comp_ind_50TY <- NULL

# make a new column which is the inoc_strain and sup_strain combined. inoc strain comes first and then sup strain.
GROdf$inoc_sup_strain <- paste(GROdf$inoculant, GROdf$supernatant, sep = "_")


# find how many replicates are for each treatment (should be 5)
rowreps <- ddply(GROdf,.(time, supernatant, sup_geno, inoculant, inoc_geno),nrow)
# there are 5 reps in all treatments so:
# Create a last colum for technical reps.
repseq <- rep(c(1,2,3,4,5), times = nrow(rowreps))
GROdf$rep <- repseq

# make a new column which is the inoc_strain and sup_strain combined and the rep - so each individual well. inoc strain comes first and then sup strain.
GROdf$inoc_sup_strain_rep <- paste(GROdf$inoculant, GROdf$supernatant, GROdf$rep, sep = "_")

# make sure GROdf is a dataframe
GROdf <- data.frame(GROdf)
str(GROdf)


# Create summary stats

# make a blank data frame with the columns to keep with a row for each unique combination
growth<-unique(GROdf[,c("inoc_sup_strain_rep", "supernatant", "sup_geno", "inoculant", "inoc_geno", "inoc_sup_strain", "rep")])

#calculate summary metrics - taken from Rosanna Wright's code. Explanation of code by Ellie Harrison.

# this uses a loop. Each interation of the loop will make a subset of the data for 1 curve, then calculates the different metrics from it.
# the loop term (for i in 1:NROW(growth)) means for each number between 1 and 'the number of rows in the dataframe GROdf do the following actions. On the first round of the loop i=1, then i=2 etc etc
# to follow what is happening try specifying i, e.g. i <- 1, then going line by line within the loop to see the data being added to the growth dataframe
for(i in 1:NROW(growth)){
  asub <- GROdf[GROdf$inoc_sup_strain_rep == growth[i,'inoc_sup_strain_rep'], ] #subset of data and make a mini dataframe called asub. Because I kept the 'inoc_sup_strain_rep' code I can just pull out the rows from the GROdf data frame where the inoc_sup_strain_rep code matches the inoc_sup_strain_rep code in the specified row (i) of the new summary dataframe
  growth[i, "AUC"] <- auc(asub[,"time"],asub[,'OD']) # Integral of a curve = area under the curve and is a good metric for growth. 
  growth[i,"MaxOD"] <-max(asub[,'OD']) # maximum absorbance reached (just picks out the maximum OD value in the asub dataframe)
  growth[i, "OD_62h"] <- asub[asub$time == 62, "OD"] # OD at 62 h time point. This should be the same as MaxOD but just to make sure.
  
  #growth rates = change in absorbance per hour. this is a litle tricky as you cant just do it across the curve (because the growth slows down). So you need to calculate rate for a buch of time windows through the growth cycle and then you can pick the highest (i.e. maximum growth rate)
  #1. data frame containing growth rates through time for each curve
  time<-c(0,24,48,62) # just makes a data frame of time - by hours
  time<-time[order(time)] # ensures it is in order
  rates <- data.frame(Rate=as.numeric(), Time=as.numeric(), Abs=as.numeric()) #make a new little dataframe to collect the rate info
  for(j in c(2:3)){  # note that j starts at 2 (24h) and ends 1 time points before the end (48h), hopefully you can work out why
    a0<-asub[asub$time == time[j-1], 'OD'] # takes the OD at time point j-24 (one hr before j)
    a1<-asub[asub$time == time[j], 'OD'] # takes the OD at time point j
    a2<-asub[asub$time == time[j+1], 'OD'] # takes the OD at time point j+24 (one hr after j)
    rates[NROW(rates)+1, 'Rate']<-(a2-a0)/48 # calculate the rate around point J (the demonimator is 48 becuase its across 48 hours)
    rates[NROW(rates),'Abs']<-a1 #  add the OD for row j
    rates[NROW(rates), 'Time']<-time[j] # add the time for row j
  }
  #maximum growth rate just taken as the highest value in the list of rates
  max <-rates[rates$Rate == max(rates$Rate),] #make a tiny dataframe containing the maximum values
  growth[i,'MaxGrowthRate']<-max$Rate[1] #if more than 1 row, we want the first time this occurs
  growth[i,'TimeMaxGrowthRate']<-max$Time[1]
  growth[i,'AbsMaxGrowthRate']<-max$Abs[1]
}


#decide which metric is best for your data. this is a function copied from Crawley's R book
# you specify the function first and then run it using pairs()
# this can help us tell what values might be beneficial to examine. 
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "grey", ...)
}
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}


pairs(~ OD_62h + AUC + MaxGrowthRate , upper.panel=panel.cor, diag.panel=panel.hist, data=growth)




# recreating heatmap for mean AUC ------------


# calculate descriptive statistics for supernatant and inoculant genospecies groups. 
growth_AUC_mean <- summarySE(growth, measurevar = 'AUC', groupvars = c( 'sup_geno',
                                                                        'supernatant', 'inoc_geno',
                                                                                       'inoculant'))
# add SM to strain names
growth_AUC_mean$supernatant <- paste0("SM", growth_AUC_mean$supernatant)
growth_AUC_mean$inoculant <- paste0("SM", growth_AUC_mean$inoculant)

# remove SM from the 100% and 50% TY treatments
growth_AUC_mean$supernatant[growth_AUC_mean$supernatant == 'SM100_TY'] <- '100% TY'
growth_AUC_mean$supernatant[growth_AUC_mean$supernatant == 'SM50_TY'] <- '50% TY'

# reformat supernatant levels and inoculant levels
growth_AUC_mean$supernatant <- factor(growth_AUC_mean$supernatant, levels = c("100% TY","50% TY","SM152B","SM137B","SM152A",
                                                                                  "SM145B","SM154C", "SM144A", "SM147A", "SM158", 
                                                                                  "SM170C", "SM157B", "SM165A", "SM122A", "SM126B",
                                                                                  "SM149A", "SM135B", "SM135A", "SM159", "SM168A",
                                                                                  "SM41", "SM53", "SM57", "SM74", "SM77", "SM67"))
growth_AUC_mean$inoculant <- factor(growth_AUC_mean$inoculant, levels = c("SM152B","SM137B","SM152A",
                                                                              "SM145B","SM154C", "SM144A", "SM147A", "SM158", 
                                                                              "SM170C", "SM157B", "SM165A", "SM122A", "SM126B",
                                                                              "SM149A", "SM135B", "SM135A", "SM159", "SM168A",
                                                                              "SM41", "SM53", "SM57", "SM74", "SM77", "SM67"))


# create horizontal (h) and vertical (v) lines for distinguishing genospecies groups in the heatmaps
h <- c(6.5, 13.5, 18.5)
v <- c(2.5, 8.5, 15.5, 20.5)

#plot: supernatant interactions heatmap
(compind_heat <- ggplot(growth_AUC_mean, aes(x = supernatant, y = inoculant, fill = AUC))+
    geom_tile(colour = "white", size=1) + 
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 5)+
    labs(x = "Supernatant strain", y = "Inoculant strain", fill = "Mean\nAUC")+
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
ggsave('../../data/intermediate_data/growth_metric_comparison_data/supernatant_interactions_heatmap_AUC.pdf', width = 10, height = 10, plot = compind_heat)








# make heatmap again where AUC in other strain supernatant is divided by AUC own supernatant (like RGI)

# make a dataframe which contains the OD values where a sample is grown in its own supernatant
own_ODs <- growth[growth$supernatant == growth$inoculant, ]
# make new column showing only inoculant and the rep number (removing all before and upto the first _ in the string)
own_ODs$inoc_rep <- gsub("[^_]*_(.*)", "\\1", own_ODs$inoc_sup_strain_rep)

# make loop to add new column onto dataframe

for(i in 1:NROW(growth)){
  # paste inoculant strain id and the replicate number. 
  id_inoc <- paste(growth[i, "inoculant"], growth[i, "rep"], sep = "_")
  # calculate auc_index. divide auc of sample by own auc. (remember to include rep number)
  growth[i, "AUC_index"] <- growth[i, "AUC"] / own_ODs[own_ODs$inoc_rep == id_inoc, "AUC"]
}

# do same for OD_62h to get RGI
for(i in 1:NROW(growth)){
  # paste inoculant strain id and the replicate number. 
  id_inoc <- paste(growth[i, "inoculant"], growth[i, "rep"], sep = "_")
  # calculate auc_index. divide auc of sample by own auc. (remember to include rep number)
  growth[i, "RGI"] <- growth[i, "OD_62h"] / own_ODs[own_ODs$inoc_rep == id_inoc, "OD_62h"]
}


# then calculate the mean with summaryse again and plot the data
growth_AUCindex_mean <- summarySE(growth, measurevar = 'AUC_index', groupvars = c( 'sup_geno',
                                                                        'supernatant', 'inoc_geno',
                                                                        'inoculant'))
# add SM to strain names
growth_AUCindex_mean$supernatant <- paste0("SM", growth_AUCindex_mean$supernatant)
growth_AUCindex_mean$inoculant <- paste0("SM", growth_AUCindex_mean$inoculant)

# remove SM from the 100% and 50% TY treatments
growth_AUCindex_mean$supernatant[growth_AUCindex_mean$supernatant == 'SM100_TY'] <- '100% TY'
growth_AUCindex_mean$supernatant[growth_AUCindex_mean$supernatant == 'SM50_TY'] <- '50% TY'

# reformat supernatant levels and inoculant levels
growth_AUCindex_mean$supernatant <- factor(growth_AUCindex_mean$supernatant, levels = c("100% TY","50% TY","SM152B","SM137B","SM152A",
                                                                              "SM145B","SM154C", "SM144A", "SM147A", "SM158", 
                                                                              "SM170C", "SM157B", "SM165A", "SM122A", "SM126B",
                                                                              "SM149A", "SM135B", "SM135A", "SM159", "SM168A",
                                                                              "SM41", "SM53", "SM57", "SM74", "SM77", "SM67"))
growth_AUCindex_mean$inoculant <- factor(growth_AUCindex_mean$inoculant, levels = c("SM152B","SM137B","SM152A",
                                                                          "SM145B","SM154C", "SM144A", "SM147A", "SM158", 
                                                                          "SM170C", "SM157B", "SM165A", "SM122A", "SM126B",
                                                                          "SM149A", "SM135B", "SM135A", "SM159", "SM168A",
                                                                          "SM41", "SM53", "SM57", "SM74", "SM77", "SM67"))


# create horizontal (h) and vertical (v) lines for distinguishing genospecies groups in the heatmaps
h <- c(6.5, 13.5, 18.5)
v <- c(2.5, 8.5, 15.5, 20.5)

#plot: supernatant interactions heatmap
(compind_heat2 <- ggplot(growth_AUCindex_mean, aes(x = supernatant, y = inoculant, fill = AUC_index))+
    geom_tile(colour = "white", size=1) + 
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1)+
    labs(x = "Supernatant strain", y = "Inoculant strain", fill = "Mean\nAUC\nindex")+
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
ggsave('../../data/intermediate_data/growth_metric_comparison_data/supernatant_interactions_heatmap_AUCindex.pdf', width = 10, height = 10, plot = compind_heat2)


# correlate mean AUC_index with mean RGI

pairs(~ RGI + AUC_index , upper.panel=panel.cor, diag.panel=panel.hist, data=growth)


# save growth data 
write.csv(growth, '../../data/intermediate_data/growth_metric_comparison_data/AUC_growth_data.csv', row.names = F)





# recreating heatmap for Max Growth Rate ------------


# calculate descriptive statistics for supernatant and inoculant genospecies groups. 
growth_MGR_mean <- summarySE(growth, measurevar = 'MaxGrowthRate', groupvars = c( 'sup_geno',
                                                                        'supernatant', 'inoc_geno',
                                                                        'inoculant'))
# add SM to strain names
growth_MGR_mean$supernatant <- paste0("SM", growth_MGR_mean$supernatant)
growth_MGR_mean$inoculant <- paste0("SM", growth_MGR_mean$inoculant)

# remove SM from the 100% and 50% TY treatments
growth_MGR_mean$supernatant[growth_MGR_mean$supernatant == 'SM100_TY'] <- '100% TY'
growth_MGR_mean$supernatant[growth_MGR_mean$supernatant == 'SM50_TY'] <- '50% TY'

# reformat supernatant levels and inoculant levels
growth_MGR_mean$supernatant <- factor(growth_MGR_mean$supernatant, levels = c("100% TY","50% TY","SM152B","SM137B","SM152A",
                                                                              "SM145B","SM154C", "SM144A", "SM147A", "SM158", 
                                                                              "SM170C", "SM157B", "SM165A", "SM122A", "SM126B",
                                                                              "SM149A", "SM135B", "SM135A", "SM159", "SM168A",
                                                                              "SM41", "SM53", "SM57", "SM74", "SM77", "SM67"))
growth_MGR_mean$inoculant <- factor(growth_MGR_mean$inoculant, levels = c("SM152B","SM137B","SM152A",
                                                                          "SM145B","SM154C", "SM144A", "SM147A", "SM158", 
                                                                          "SM170C", "SM157B", "SM165A", "SM122A", "SM126B",
                                                                          "SM149A", "SM135B", "SM135A", "SM159", "SM168A",
                                                                          "SM41", "SM53", "SM57", "SM74", "SM77", "SM67"))


# create horizontal (h) and vertical (v) lines for distinguishing genospecies groups in the heatmaps
h <- c(6.5, 13.5, 18.5)
v <- c(2.5, 8.5, 15.5, 20.5)

#plot: supernatant interactions heatmap
(MGR_heat <- ggplot(growth_MGR_mean, aes(x = supernatant, y = inoculant, fill = MaxGrowthRate))+
    geom_tile(colour = "white", size=1) + 
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.004)+
    labs(x = "Supernatant strain", y = "Inoculant strain", fill = "Mean\nMax Growth Rate")+
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
ggsave('../../data/intermediate_data/growth_metric_comparison_data/supernatant_interactions_heatmap_MaxGrowthRate.pdf', width = 10, height = 10, plot = MGR_heat)








# make heatmap again where AUC in other strain supernatant is divided by AUC own supernatant (like RGI)


# make loop to add new column onto dataframe: MaxGrowthRate_index

for(i in 1:NROW(growth)){
  # paste inoculant strain id and the replicate number. 
  id_inoc <- paste(growth[i, "inoculant"], growth[i, "rep"], sep = "_")
  # calculate auc_index. divide auc of sample by own auc. (remember to include rep number)
  growth[i, "MaxGrowthRate_index"] <- growth[i, "MaxGrowthRate"] / own_ODs[own_ODs$inoc_rep == id_inoc, "MaxGrowthRate"]
}



# then calculate the mean with summaryse again and plot the data
growth_MGRindex_mean <- summarySE(growth, measurevar = 'MaxGrowthRate_index', groupvars = c( 'sup_geno',
                                                                                   'supernatant', 'inoc_geno',
                                                                                   'inoculant'))
# add SM to strain names
growth_MGRindex_mean$supernatant <- paste0("SM", growth_MGRindex_mean$supernatant)
growth_MGRindex_mean$inoculant <- paste0("SM", growth_MGRindex_mean$inoculant)

# remove SM from the 100% and 50% TY treatments
growth_MGRindex_mean$supernatant[growth_MGRindex_mean$supernatant == 'SM100_TY'] <- '100% TY'
growth_MGRindex_mean$supernatant[growth_MGRindex_mean$supernatant == 'SM50_TY'] <- '50% TY'

# reformat supernatant levels and inoculant levels
growth_MGRindex_mean$supernatant <- factor(growth_MGRindex_mean$supernatant, levels = c("100% TY","50% TY","SM152B","SM137B","SM152A",
                                                                                        "SM145B","SM154C", "SM144A", "SM147A", "SM158", 
                                                                                        "SM170C", "SM157B", "SM165A", "SM122A", "SM126B",
                                                                                        "SM149A", "SM135B", "SM135A", "SM159", "SM168A",
                                                                                        "SM41", "SM53", "SM57", "SM74", "SM77", "SM67"))
growth_MGRindex_mean$inoculant <- factor(growth_MGRindex_mean$inoculant, levels = c("SM152B","SM137B","SM152A",
                                                                                    "SM145B","SM154C", "SM144A", "SM147A", "SM158", 
                                                                                    "SM170C", "SM157B", "SM165A", "SM122A", "SM126B",
                                                                                    "SM149A", "SM135B", "SM135A", "SM159", "SM168A",
                                                                                    "SM41", "SM53", "SM57", "SM74", "SM77", "SM67"))


# create horizontal (h) and vertical (v) lines for distinguishing genospecies groups in the heatmaps
h <- c(6.5, 13.5, 18.5)
v <- c(2.5, 8.5, 15.5, 20.5)

#plot: supernatant interactions heatmap
(MGR_heat2 <- ggplot(growth_MGRindex_mean, aes(x = supernatant, y = inoculant, fill = MaxGrowthRate_index))+
    geom_tile(colour = "white", size=1) + 
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1)+
    labs(x = "Supernatant strain", y = "Inoculant strain", fill = "Mean\nMaxGrowthRate\nindex")+
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
ggsave('../../data/intermediate_data/growth_metric_comparison_data/supernatant_interactions_heatmap_MaxGrowthRateindex.pdf', width = 10, height = 10, plot = MGR_heat2)


# correlate mean AUC_index with mean RGI

pairs(~ RGI + AUC_index + MaxGrowthRate_index , upper.panel=panel.cor, diag.panel=panel.hist, data=growth)


# save growth data 
write.csv(growth, '../../data/intermediate_data/growth_metric_comparison_data/AUC_MaxGrowthRate_growth_data.csv', row.names = F)




# correlate the growth_AUCindex_mean values in the heatmap with the original RGI_mean values in the original heatmap -----

# read in the growth dataset "fresh" - called "AUC_MaxGrowthRate_growth_data.csv"
growth <- read.csv("../../data/intermediate_data/growth_metric_comparison_data/AUC_MaxGrowthRate_growth_data.csv")

# calculate RGI descriptive statistics for supernatant and inoculant genospecies groups. 
RGI_mean <- summarySE(growth, measurevar = 'RGI', groupvars = c( 'sup_geno',
                                                                        'supernatant', 'inoc_geno',
                                                                        'inoculant'))

# calculate AUC descriptive statistics for supernatant and inoculant genospecies groups. 
AUCindex_mean <- summarySE(growth, measurevar = 'AUC_index', groupvars = c( 'sup_geno',
                                                                      'supernatant', 'inoc_geno',
                                                                      'inoculant'))


# merge datasets together
merged_data <- merge(RGI_mean[1:6], AUCindex_mean[1:6], by = c('sup_geno',
                                      'supernatant', 'inoc_geno',
                                      'inoculant', 'N'))

# make a correlation plot

(plot_indexcomp <- ggplot(data = merged_data, aes(x = RGI, y = AUC_index)) +
    geom_point() +
    geom_smooth(method = "lm", se = F) + 
    labs(x = "Mean RGI heatmap values", y = "Mean AUC index heatmap values")+
    theme_bw()+
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16), 
          axis.ticks=element_line(size=0.4), 
          text = element_text(size=25))
)
#function to save specific ggplot graph.  additional parameters(width=,height=,units='cm')
ggsave('../../data/intermediate_data/growth_metric_comparison_data/correlation_RGI_AUCindex_heatmapvalues.pdf', width = 7, height = 7, plot = plot_indexcomp)



# simple linear regression between the RGI and AUC_index
model1 <- lm(RGI ~ AUC_index, data = merged_data)
summary(model1)

# Call:
#   lm(formula = RGI ~ AUC_index, data = merged_data)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.33286 -0.05988 -0.00553  0.05541  0.35639 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.29117    0.01392   20.92   <2e-16 ***
#   AUC_index    0.72368    0.01358   53.29   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.09347 on 622 degrees of freedom
# Multiple R-squared:  0.8203,	Adjusted R-squared:   0.82 
# F-statistic:  2840 on 1 and 622 DF,  p-value: < 2.2e-16


