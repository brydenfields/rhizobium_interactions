rhizobium_interactions: scripts for figure generation and statistical analyses
===

This repository provides the R scripts and raw data used to generate the figures and statistical analyses in the manuscript by Fields et al. "Genetic variation is associated with differences in facilitative and competitive interactions in the Rhizobium leguminosarum species complex", currently in prep for publication.

The R scripts and raw data are provided in this repository to enable regeneration of the data analysis undertaken in the study. The figures and intermediate data files generated within the scripts are also saved here. The scripts have been edited so that all folders can be downloaded and the files run locally on a user's computer. The working directory is hard-coded and will need to be altered by the user for the scripts to run - user's should set the working directory for each script to the local directory path of the specific script being run.

The following directories are included:

## scripts/
-----

The following scripts are organised into 4 separate folders depending on which experimental data is being used. 
Scripts are run with R unless otherwise stated. 
All scripts include comments detailing the code function at each step.

## comparative_genomics_scripts

**BLASTn_quorumsensing_genes.txt** used for processing BLASTn output of quorum sensing genes in linux.

**comparative_genomics_heatmap.R** used for ***Figure 5***


## spotplating_scripts

**inhibitionzone_measurements_analysis.R** used for ***Figures 2b, S1, S5, S9, S10*** 	

**spotplating_MEM_analysis.R** used for ***Figure S2b and S3c-d***	


## substrate_utilisation_scripts			

**absolute_substrateutilisation_analysis.R** used for ***Figures 4, S6, S7***	


##supernatant_scripts

**inoculant_supernatant_effects_correlation.R** used for ***Figure 3***  
			
**supernatant_MEM_analysis.R** used for ***Figure S2a, S3a-b***
		
**supernatant_RGI_analysis.R** used for ***Figure 2a, S4, S9***


## data/raw_data/
-----

The following raw data files are organised into 4 separate folders depending on which experimental data is being used - as with the script directory layout above. 

## comparative_genomics_data

**antiSMASH_totalclusters_summary.csv** contains a count matrix of antiSMASH clusters identified for each strain.

**genospecies_196strain_presabs_percentage.csv** contains a matrix of the percentage of strains containing a gene in each genospecies. Data generated from Cavassim et al., 2020. https://doi.org/10.1099/mgen.0.000351

**phagecount_summary.csv** contains a count matrix of PHASTER prophage regions identified for each strain.

**quorumsensing_percentageidentity_blasttable.csv** contains a percentage identity matrix of quorum sensing genes identified for each strain.


## spotplating_data (direct interactions experiment)

**dataset_IZspot_0320.csv** contains inhibition zone Feret diameter measurements produced by ImageJ for each strain combination. 

**dataset_IZspot_diameterdiffs_averages.csv** contains inhibition zone diameter summary descriptive statistics (average = "IZ_diff") for each strain combination.

**dataset_IZspot_diameterdiffs.csv** contains inhibition zone Feret diameter measurements produced by ImageJ for each strain combination, along with a difference calculation between the inhibition zone diameter and the culture spot diameter.

**inoculum_inital_ODs_allreps.csv** contains the initial culture spot OD600 values of strains spotted onto soft agar lawns. Strain combinations were repeated up to 4 times as a result of some soft agar lawn strains failing to grow efficiently. 


## substrate_utilisation_data

**compind_sup_inoc_AWCD.csv** contains growth summary descriptive statistics for strains when acting as the inoculant or the supernatant in the indirect interaction experiment, including the average well colour development values for growth in 31 ecoplate substrates.

**num_substrates_used_ecoplate.csv** contains ecoplate substrate utilisation data for the number of substrates used by each strain.

**24sup_growth_phenotype_AWCDgroups_PCAdata.xlsx** contains ecoplate substrate utilisation data for PCA analysis.

**24sup_growth_phenotype_formatted.xlsx** contains raw OD values of ecoplate substrate utilisation data for each strain.


## supernatant_data (indirect interactions experiment)

**ani_sorted_by_genospecies_snps_new_data_6kgenes_supernatantsamples.csv** contains an average nucleotide identity matrix for strains used in this study. Original data from Cavassim et al., 2020. https://doi.org/10.1099/mgen.0.000351

**bes_project_all_50TY.csv** contains raw OD values and relative growth index values for growth compared to in 50% TY. 

**bes_project_all_no1.csv** contains raw OD values and relative growth index values for growth compared to in a strain's own supernatant (indirect interactions experiment).

**combined_compind_sup_inoc.csv** contains growth summary descriptive statistics for strains when acting as the inoculant or the supernatant in the indirect interaction experiment.


## data/intermediate_data
-----

An empty set of folder directories where the output data files and figure will be saved when the scripts are run. The intermediate_data folders are organised into 4 separate folders depending on which experimental data is being used - as with the script/data directory layouts above. 


## data/intermediate_data_expected_output
-----

A folder containing the expected output data files and figures generated by the scripts. (For user reference).



---------
---
Author: Bryden Fields - bryden.fields@york.ac.uk