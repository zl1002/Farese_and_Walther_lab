# Script: FWL_lipidomics_3.1.R
# Author: Wenting 
# Notes:  To start, typing command in the console-----> source("FWL_lipidomics_3.1.R") 
#         or press the source butthon. 
#         Please make sure Mac users installed XQuartz.
#         This script is designed as a reference for lipidomics experiment. 
#         It is based on Niklas and Kenny's previous work (Their work files can be found in folders 
#         quality_control and statistics_quantification). Acknowledge to Weng's technical 
#         guidance, Laura's fatty acid saturation analysis project and Sebstian's shadow experiment
#
# Usage:  
#      
# Status:  In progress             

#####################################################################################
rm(list=ls())  ##### empty the environment first. If you need to keep former variables in workspace, just comment this part
setSessionTimeLimit(cpu = Inf, elapsed = Inf)
# source function from FWL_lipidomics.**.functions.R script
source("FWL_lipidomics_3.2.2.FUNCTIONS.R")
# Check directory existence, if not then make one
# all the plots will stored in plot directory and data in the data directory
dirs <- c("plot", "data", "plot/classes", "plot/QC", "plot/Quantification", "plot/Saturation", "plot/Ether", "plot/Length",
          "plot/Volc", "data/Volc","data/QC", "data/Quantification", "data/Saturation", "data/Ether", "data/Length")
mkdirs(dirs)

message("Are you using PC or MAC?")
robot <- retype_choice("PC/MAC")
if(robot == "mac"){
  message("\n\nPlease make sure you installed Quartz before running pipeline.
          \nDo you want to continue?\n")
  process <- retype_choice("Y/N")
}


# read the file from csv directory
csv_files <- list.files(path = "converted", pattern = "\\.csv$")
csv_list <- c()
for(i in seq_along(csv_files)) csv_list[i] <- addquotes(!!as.character(i), " ", !!csv_files[i], "\n")
message("\nThe following files had been generated. \nSelect ONE for subsequent the list of file names:\n", csv_list)
file_option <- readline("Please input the index of the file: ") %>% as.numeric()
# the file path and name
target_file <- paste("converted/", csv_files[file_option], sep = "")
print(target_file)
lipidomics <- read_csv(target_file, col_types = cols())


##########################################################################################
# QC part I
##########################################################################################
source("FWL_lipidomics_QC_3.2.2.R")




##########################################################################################
# Input group information
##########################################################################################
# Make groups
# input the group information and reterieve the data from the csv file.
samples <- Input(filtered_lipidomics2)
# sample info
sample_info <- samples[[1]]
# # group name info
# group_names <- samples[[2]]
# # group number info
# ngroups <- samples[[3]]


##########################################################################################                                   
# PCA and correlation visualization
##########################################################################################
# pca and correlation plots
label <- "initial"
info_list <- PCA_pairs_Plot(sample_info, filtered_lipidomics2, label)
# check if deleting samples needed and plot new PCA
message("\nPlots can visualized under 'plot' directory or r studio plots panel.\nDo you want to edit sample information for subsequent analyses?")
pca_check <- retype_choice("Y/N")
info_list <- PCAcheck(pca_check, filtered_lipidomics2)

# making group repeats according to its position for making groups of later PCA
sample_raw_list <- info_list[[1]]
group_repeats <- info_list[[2]]
# make a data frame contains sample information and group information
group_info <- data.frame(samples=sample_raw_list, 
                         groups=group_repeats, 
                         stringsAsFactors = FALSE) %>% 
  group_by(groups) 
write_csv(group_info, "data/group_information.csv")
group_names <- unique(group_repeats)
ngroups <- length(group_names)
# # delete unchoosed samples from sample list
# total_samples <- filtered_lipidomics %>% select(contains("MainArea[s")) %>% colnames()
# deleted_samples <- total_samples %>% subset(!total_samples %in% sample_raw_list) 
# filtered_lipidomics <- filtered_lipidomics %>% select(-all_of(deleted_samples))
# deleted <- deleted_samples %>% str_remove_all(., "MainArea\\[") %>% str_remove_all(., "\\]") %>% paste0(., sep = ",", collapse = "") %>% remove(.,  pos = -1)
# if(length(deleted)>0){
#   message("Samples ", deleted, " will be removed from analysis")
# }





###########################################################################################
# Background subtraction, filter potential invalid lipids
###########################################################################################
message("\nsubtract sample area from background/solvent run?" )
background_option <- retype_choice("Y/N")
if(background_option == "y"){
    subtracted_data <- subtract_background(filtered_lipidomics2, sample_raw_list, "MainArea[c]")
    write_csv(subtracted_data, "data/subtracted_lipids.csv")
    filtered_lipidomics_copy <- filtered_lipidomics2   ######## filtered_lipidomics_copy will store the raw filtered data (filtered_lipidomics)
    filtered_lipidomics3 <- subtracted_data   ######### if doing background subtraction, the filtered lipidomics data will be replaced by subtracted data
    # detect potential invalid lipid molecules (empty or negative value in the samples)
    invalid_lipids <- filtered_lipidomics3 %>% filter_at(sample_raw_list, any_vars(. <= 0)) 
    # add SFAE and UNSAFE flag for indicating potential invalid lipid molecules
    filtered_lipidomics3 <- filtered_lipidomics3  %>% 
      rowwise() %>% 
      mutate(FLAG_invalid = ifelse(LipidMolec %in% invalid_lipids$LipidMolec, "UN_SAFE", "SAFE")) 
    # filter potential invalid lipid molecules
    filtered_lipidomics <- filter_invalid(filtered_lipidomics3, group_info, invalid_lipids)
    subtracted_data <- subtracted_data %>% filter(LipidMolec %in% filtered_lipidomics$LipidMolec)
}else{
  filtered_lipidomics <- filtered_lipidomics2
}






###########################################################################################
# Quantification of lipid class and individual lipid molecules
###########################################################################################
source("FWL_lipidomics_QUANTIFICATION_3.2.2.R")






##########################################################################################
# saturation analysis
##########################################################################################
# calculate the saturation for different lipid class
message("\nThe filtered data will be used for saturation analysis.")
source("FWL_FattyAcids_Saturation_3.3.3.R")

message("\nThe SFA, MUFA, PUFA information will be stored in the count_lipid.csv and aggregated.csv")

##########################################################################################
# fatty acids length analysis
##########################################################################################
source("FWL_FattyAcids_Length_3.2.2.R")




##########################################################################################
# ether lipid analysis
##########################################################################################
message("\nThe filtered data will be used for ether lipid analysis.\n")
source("FWL_EtherLipids_3.2.2.R")










##########################################################################################
# test random sample distribution
##########################################################################################
# test sample distribution
# log transformation and visualize approximately normal distribution of the transformed data
log2_lipids <- filtered_lipidomics %>% mutate_at(sample_raw_list, log2trans)
# log2_lipids %>% filter_at(sample_raw_list, any_vars(.<= 0)) 
# randomly choose a sample to see its distribution approximately normal distribution
i <- sample(length(unique(sample_raw_list)), 1)
message(paste("MainArea[s", i, "]", sep=""), " is chosen for plotting distribution")
# the data is approximately normal distribution
plot_all(log2_lipids, sample_raw_list[i]) +
  geom_density() +
  labs(x = "sample distribution (log transformed AUC value)", y = "lipid molecule count") 


##########################################################################################
# log transformation of subtracted/filtered data
##########################################################################################
write.csv(log2_lipids, "data/log.molec.csv")






##########################################################################################
# volcano plot
##########################################################################################

# Create a design matrix 
samples <- factor(group_repeats)
design <- model.matrix(~0+samples)
colnames(design) <- levels(samples)
# message("Begin to plot volcano \nPlease note that you NEED make contrast groups manually to Compare the difference between/among groups.\ne.g. compare B+C against A: A-(B+C)/2; A against B:  A-B; A-B against C-D: (A-B)-(C-D), etc.")

ncomparisons <- readline("How many volcano plots to generate: \n") %>% as.numeric()
impute_dt <- impute_not("y", filtered_lipidomics, sample_raw_list)
imputated_lipids <- impute_dt[[2]]
imputated_lipids <- imputated_lipids %>% ungroup() %>% as.data.frame(., stringsAsFactors=TRUE)
rownames(imputated_lipids) <- imputated_lipids %>% select(LipidMolec) %>% unlist()
processed_lipids <- imputated_lipids %>% select(-Class, -LipidMolec)

# Fit model and extract contrasts 
fit <- lmFit(processed_lipids, design)
# plot volcano graph
VolPlot(ncomparisons, fit, processed_lipids)



