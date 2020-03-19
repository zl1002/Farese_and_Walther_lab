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
source("FWL_lipidomics_3.2.0.FUNCTIONS.R")
# Check directory existence, if not then make one
# all the plots will stored in plot directory and data in the data directory
dirs <- c("plot", "data", "plot/classes")
mkdirs(dirs)

# read the file from csv directory
csv_files <- list.files(path = "converted", pattern = "\\.csv$")
csv_list <- c()
for(i in seq_along(csv_files)) csv_list[i] <- addquotes(!!as.character(i), " ", !!csv_files[i], "\n")
message("The following files had been generated. Select ONE for subsequent the list of file names:\n", csv_list)
file_option <- readline("Please input the index of the file: ") %>% as.numeric()
# the file path and name
target_file <- paste("converted/", csv_files[file_option], sep = "")
print(target_file)
lipidomics <- read_csv(target_file, col_types = cols())


##########################################################################################
# check experiment samples
##########################################################################################
# check if there is internal control samples or samples don't use for analysis
# and extract information of different grades and p value 
# which are less or equal than 0.001 for each molecule and store this info in the variable lipid_count
message("\nDo you have experimental controls like internal standards or extraction control?")
delete_option <- readline("Please type Y/N: ") %>% str_to_lower(.)
lipid_check <- delete_samples(delete_option, lipidomics)
# counter of Grades A and B, filter P value less than 0.001 for LipidSearch
lipid_count <- lipid_check[[1]]
# remove samples
lipid_remove <- lipid_check[[2]]
# add this count into original strucutes
#lipid_select <- lipidomics %>% bind_cols(lipid_count) #### uncomment it, if reserve the experiment controls.
lipid_select <- lipid_remove %>% bind_cols(lipid_count)  ##### comment it, if reserve the experiment controls.


##########################################################################################
# set filter parameters
##########################################################################################
# Filter the dataset based on your criteria  
# flexible parameter k and j which depends on the experiment for total number of Grade A and B and APvalue
message("\nData are filtered using 3 criteria.  
\n 1. Not rejected by LipidSearch (Rej n=0); 
        \n 2. minimum number of Grade A+B required; 
        \n 3. minimum number significant identification (p-value p<=0.001) for LipidSearch standard.
        \n Filtered data is stored in the filtered.raw.data.csv.")
k <- readline("Minimum number of identified molecules (Grade A+B) required in all samples, n>=3:  ") %>% as.numeric(.)
j <- readline("Minumum number of significantly identified lipids p<=0.001 in all samples, n>=3: ") %>% as.numeric(.)
filtered_lipidomics <- lipid_select %>% 
  rowwise() %>% 
  filter(Rej == 0 & sum(A, B) >= k & APValue.001 >= j) %>% 
  as.data.frame(.)
# the filtered raw data is stored in data/filtered.raw.data.csv
write.csv(filtered_lipidomics, "data/filtered.raw.data.csv")

##########################################################################################
# check background information
##########################################################################################
# plot the background information of the blank sample c
blank_sample <-filtered_lipidomics %>% 
  select(Class, contains("MainArea[c]")) %>% 
  group_by(Class) %>% 
  summarise_all(.funs=sum)
parameters_bw <- colnames(blank_sample) <- c("class", "value")

# plot the blank sample
plot_all(blank_sample, parameters_bw) +
  geom_bar(position="dodge", stat="identity") +
  labs(x = "Lipid Classes", 
       y = (bquote("Fold change: " ~ log[10]~"("~AUC~")")), 
       title = "Area under Curve for blank sample (background) in each class") +
  # scale_y_continuous(trans = pseudo_log_trans(base = 10),
  #                    breaks = c(0, 10^3, 10^6, 10^9, 10^12),
  #                    labels = c(0, bquote(10^3), bquote(10^6), bquote(10^9), bquote(10^12)),
  #                    expand = c(0, 0, 0.1, 0)) +
  add_scales() +
  coord_flip() 
# the blank sample plot is stored in background.png
ggsave(filename = "background.png", path="plot/", device = "png")


##########################################################################################
# lipid class summary
##########################################################################################
message("\nThe info below and summary plot will show the summary information of classes after filtering the data")
# two methods check how many lipids passed filtering
print(describe(filtered_lipidomics$Class))
filtered_lipidomics %>% group_by(Class) %>% summarise(lipid_class_num = n()) %>% select(Class, lipid_class_num) %>% arrange(lipid_class_num)  %>% formattable(.)

# get lipid class proportion information
prop_data <- filtered_lipidomics
class(prop_data$Class) <- factor(prop_data$Class)
prop_data <- prop_data %>% 
  group_by(Class) %>% 
  summarise(count = n()) %>% 
  mutate(prop = count/sum(count))

# plot lipid class proportion information
plot_all(prop_data, c("Class", "prop")) +
#  ggplot(prop_data, aes(x = reorder(Class, prop), y = prop)) +
  theme_bw() + 
  set_theme() +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels=scales::percent, expand = c(0, 0, 0.2, 0)) +
  ggtitle("Class Summary plot: \n Proportion of different class numbers among all the samples")+
  xlab("classes")+ ylab("Relative frequencies (%)")  +
  coord_flip() + 
  scale_fill_grey()
ggsave(filename = "prop.bw.summary.png", path = 'plot/', device = "png")

# save lipid class summary information into file proportion_classes.csv in data directory
prop_data %>%
  mutate(prop = percent(prop*100/100)) %>%
  rename(number_of_lipid_molec = count, proportion = prop)%>% 
  write_csv(., "data/proportion_classes.csv")


##########################################################################################
# retention time analysis
##########################################################################################
###  Abundance vs. retention time for all samples
retention_data <- filtered_lipidomics %>% 
  select(contains("MainArea[s"), BaseRt, Class) %>% 
  gather(sample, MainArea, -c(BaseRt, Class))
n_classes <- retention_data$Class %>% unique() %>% length()
pars <- c("MainArea", "BaseRt", n_classes)
rt_plot <- plot_rt_allsamples(retention_data, pars)
print(rt_plot)
ggsave(filename = "all.retention.png", path = 'plot/', device = "png", 
       width = 10, height = 8, dpi = 150, units = "in")



##########################################################################################
# mark odd chains
##########################################################################################
# mark classes which contain the odd chains classes 
odd_index <- filtered_lipidomics$LipidMolec %>%
  str_locate(., "(\\d[13579]:)|([13579]:)")
odd_chains <- filtered_lipidomics[unique(which(!is.na(odd_index), arr.ind=TRUE)[,1]), ]
# store odd chains in file odd_chains
write_csv(odd_chains, "data/odd_chains.csv")
percent_odd <- nrow(odd_chains)/nrow(filtered_lipidomics) 
percent_odd <- percent(percent_odd*100/100)
# get odd chain percentage in the raw data
message("\nThere are ", nrow(odd_chains), " lipid molecules contain odd chains. ",
        "\nThe odd chain of fatty acids percent is ", percent_odd, " in total.")
# mark the odd chains in TG class
message("The odd chain of TG information is stored in TG_odd.csv.")

##########################################################################################
# test standard TG(17:1/17:1/17:1) abundance in all samples
# Please note that this part could be editted for other funture standards
###########################################################################################
#Please have a look and decide if you wanna delete any lipid molecules. If you do, just delete the corresponding lipid in TG.odd.csv. The pipeline will take care of the rest
TGs <- odd_chains %>% filter(Class=="TG") %>% arrange(FA)
write.csv(TGs, "data/TG_odd.csv")
# you can replace the other trusted standards instead of TG(17:1/17:1/17:1)
TG17 <- TGs %>% filter(LipidMolec == "TG(17:1/17:1/17:1)")  %>% select(LipidMolec, contains("MainArea"))
# If used TG(17:1/17:1/17:1) as internal standards in the experiment
if(nrow(TG17) != 0){
  colnames(TG17) <- colnames(TG17) %>% str_replace_all(., "MainArea\\[", "") %>% str_remove_all(., "\\]")
  TG17_trick <- TG17[,-1] + 1
  TG17_trick <- t(TG17_trick) %>% data.frame()
  TG17_trick <- cbind(rownames(TG17_trick), TG17_trick)
  paras2 <- colnames(TG17_trick) <- c("samples", "value")
  tp17 <- plot_all(TG17_trick, paras2)+ 
    geom_bar(position="dodge", stat="identity") +
    labs(title = "TG(17:1/17:1/17:1) lipid molecules in all samples", x= "samples",
         y = "aggregated AUC (Main Area under curve)") +
    # add_scales(scale.params = list(trans = tn,
    #                                breaks = c(-10^3, 0, 10^3, 10^6, 10^9, 10^12, 10^15), 
    #                                labels = c(bquote(-10^3), 0, bquote(10^3), bquote(10^6),
    #                                           bquote(10^9), bquote(10^12), bquote(10^15)),
    #                                expand = c(0, 0, 0.2, 0))) +
   add_scales() +
    coord_flip()
  print(tp17)
  ggsave(filename = "TG17_all.png", path="plot/", device="png")
}


##########################################################################################
# Fix duplicated lipid molecules
##########################################################################################
# detect same lipid molecules with different retention time
# Filter the lipid molecule contains same name but different retention time based on your criteria  
# all the dupliated molecules will be stored in dupliates.csv
duplicate_molecs <- detect_duplicates(filtered_lipidomics)
# if move on, fix method for duplicated molecules
paras <- c("Class", "LipidMolec", "BaseRt", "MainIon")
if(nrow(duplicate_molecs) != 0){
  message("\n!!!Attention: Identical lipid molecules with multiples retention time. Please note that the duplicate lipid molecules are stored in reserved_duplicates.csv \n !!!!!! Potential sample contamination. \n To PROCEED, pick one: ")
  write_csv(duplicate_molecs, "data/duplicated.molecules.csv")
  # diff_baseRT will store the retention time differences for the same lipid molecule
  diff_baseRT <- duplicate_molecs %>% 
    select(LipidMolec, BaseRt, MainIon) %>% 
    group_by(LipidMolec, MainIon) %>% 
    summarise(gap=max(BaseRt)-min(BaseRt)) 
  message("\nDifferences in retention time for identical lipid molecule are stored under diff_RT.csv")
  write_csv(diff_baseRT, "data/diff_RT.csv")
  # two methods to fix duplicate lipid molecules
  message("\n A: use only ONE lipid molecule with largest main area under curve, OR \n B: Summation of main area under curve of ALL identical lipid molecule.")
  fix_method <- readline("ENTER 'A' or 'B': ")
  fixed_dt<- fix_duplicate(filtered_lipidomics, duplicate_molecs, fix_method, paras)
  # transform duplicate lipid molecules
  duplicate_output <- fixed_dt[[1]]
  # store reserved duplicate lipid molecules
  write_csv(duplicate_output, "data/reserved_duplicates.csv")
  # get filtered lipid information within processed duplicate lipid molecules
  filtered_lipidomics <- fixed_dt[[2]]
  message("Filtered lipid molecules sans duplicates are stored under removeduplicates.csv")
  write_csv(filtered_lipidomics, "data/removeduplicates.csv")
} else{
  message("There is no duplicate molecules in your data")
}


##########################################################################################
# Input group information
##########################################################################################
# Make groups
# input the group information and reterieve the data from the csv file.
samples <- Input(filtered_lipidomics)
# sample info
sample_info <- samples[[1]]
# group name info
group_names <- samples[[2]]
# group number info
ngroups <- samples[[3]]


##########################################################################################                                   
# PCA and correlation visualization
##########################################################################################
# pca and correlation plots
label <- "initial"
info_list <- PCA_pairs_Plot(sample_info, group_names, filtered_lipidomics, label)
# check if deleting outlier samples needed and plot new PCA
message("\nPlots can visualized under 'plot' directory or r studio plots panel.\nDo you want to edit sample information for subsequent analyses?")
pca_check <- readline("Please type Y/N: ")
info_list <- PCAcheck(pca_check, filtered_lipidomics)

# making group repeats according to its position for making groups of later PCA
sample_raw_list <- info_list[[1]]
group_repeats <- info_list[[2]]
# make a data frame contains sample information and group information
group_info <- data.frame(samples=sample_raw_list, 
                         groups=group_repeats, 
                         stringsAsFactors = FALSE) %>% 
  group_by(groups) 
write_csv(group_info, "data/group_information.csv")

# delete unchoosed samples from sample list
total_samples <- filtered_lipidomics %>% select(contains("MainArea[s")) %>% colnames()
deleted_samples <- total_samples %>% subset(!total_samples %in% sample_raw_list) 
filtered_lipidomics <- filtered_lipidomics %>% select(-all_of(deleted_samples))
deleted <- deleted_samples %>% str_remove_all(., "MainArea\\[") %>% str_remove_all(., "\\]") %>% paste0(., sep = ",", collapse = "") %>% remove(.,  pos = -1)
if(length(deleted)>0){
  message("Samples ", deleted, " will be removed from analysis")
}


##########################################################################################
# saturation analysis
##########################################################################################
# calculate the saturation for different lipid class
message("\nThe filtered data will be used for saturation analysis.")
source("FWL_FattyAcids_3.2.0.R")

message("\nThe SFA, MUFA, PUFA information will be stored in the count_lipid.csv and aggregated.csv")



##########################################################################################
# ether lipid analysis
##########################################################################################
message("\nThe filtered data will be used for ether lipid analysis.")
source("FWL_EtherLipids_3.2.0.R")







###########################################################################################
# Background subtraction or not for quantification analysis
###########################################################################################
subtracted_data <- subtract_background(filtered_lipidomics, sample_raw_list, "MainArea[c]")
message("\nsubtract sample area from background/solvent run?" )
background_option <- readline("Please type Y/N: ") %>% str_to_lower()
while(background_option !="n"){
  if(background_option == "y"){
    write_csv(subtracted_data, "data/subtracted_lipids.csv")
    filtered_lipidomics_copy <- filtered_lipidomics   ######## filtered_lipidomics_copy will store the raw filtered data (filtered_lipidomics)
    
    filtered_lipidomics <- subtracted_data   ######### if doing background subtraction, the filtered lipidomics data will be replaced by subtracted data
    # message("Background removed lipids are stored in the file subtracted_lipids.csv")
    break
  }else{
    print("Invalid command, please type: Y/N")
    background_option <- readline("Please type Y/N: ") %>% str_to_lower()
  }
}


##########################################################################################
# Check and process potential invalid lipids (0 or negative values)
# manually or automatically
##########################################################################################
# test how many 0 value or negative value after background subtraction
# detect potential invalid lipid molecules (empty or negative value in the samples)
invalid_lipids <- filtered_lipidomics %>% filter_at(sample_raw_list, any_vars(. <= 0)) 
# add SFAE and UNSAFE flag for indicating potential invalid lipid molecules
filtered_lipidomics <- filtered_lipidomics  %>% 
  rowwise() %>% 
  mutate(FLAG_invalid = ifelse(LipidMolec %in% invalid_lipids$LipidMolec, "UN_SAFE", "SAFE")) 
# automatically dealing lipid molecules which negative percentage is over 50% in all groups,
# or manually processing potential invalid lipid molecules
if(nrow(invalid_lipids) != 0){
  # data negative and empty value percent information 
  neg_percent_info <- detect_invalid(invalid_lipids, group_info) 
  neg_percent <- neg_percent_info[[1]]
  # transform all negative value with NA
  neg_info <- neg_percent_info[[2]] %>% 
    select(-contains("APValue"), -contains("ARatio"),
           -contains("Grade"), -c(A, B, C, D, FA))
  # replace all values for a group into NA if the negative percentage is over 50%
  negs_all <- neg_percent_info[[3]] %>% 
    select(-contains("APValue"), -contains("ARatio"),
           -contains("Grade"), -c(A, B, C, D, FA))
  # the invalid data after background subtraction.
  message("\n
For lipid molecules that contain zero values or negative values (background subtracted), 
These values are subsequently replaced as non-valid values (NA). 
Fold change analyses is performed using only samples containing valid values")
  # write negative and 0 percentage inforamtion
  write_csv(neg_percent, "data/neg.percent.csv")
  # write potential invalid lipid molecules information
  write_csv(neg_info, "data/checkInvalid.csv")
  # write copy of potential lipid molecules information
  write_csv(neg_info, "data/invalid.csv")     ## copy data for invalid lipids information
  message("Please view file imputeNA.csv for all the data contains negative values after background subtraction.")
  # write data which transform negative into NA
  write_csv(negs_all, "data/imputeNA.csv")
  message("\nType 1 if you would like the pipleline to proceed with this function \nType 2 if you prefer to exlcude certain lipid molecules for fold change analysis ")
  option <- readline("Please type 1/2: ") %>% str_to_upper()
  filtered_negs <- negs_all %>% filter_at(sample_raw_list, any_vars(!is.na(.)))
  if(option == "2"){
    # this step need manually editting the invalid lipid molecules on your computer for advanced users
    message("Select 'checkInvalid.csv' to manually exclude specific lipid molecules and click SAVE.")
    # here stops 10 seconds
    Sys.sleep(10)
    # deleting the lipid molecules you select in the checkInvalid.csv
    #message("Now we need to open the changed file checkInvalid.csv after you deleting the invalid lipid molecules.")
    continues <- readline("If you finished preprocess the data, please continue and press Y: ")
    # advanced users
    manual_data <- read_csv("data/checkInvalid.csv", col_types = cols())
    # deleting the invalid lipid molecules in raw data based on your standards and saved it in pre_filtered.lipidomics.csv
    filtered_lipidomics <- fix_invalid_by_choice(filtered_lipidomics, manual_data, filtered_negs)
    write_csv(filtered_lipidomics, "data/manual_filtered.lipidomics.csv")
    subtracted_data <- fix_invalid_by_choice(subtracted_data, manual_data, filtered_negs)
    write_csv(subtracted_data, "data/post_subtracted.csv")
  }
  while(option != "2"){
    if(option ==  "1"){
      message("The pipeline will first transform all the negative value into NA.")
      message("If negative percentage is over 50% in a group, all the values in the group for the molecule will be transformed into NA.")
      message("If a molecule which negative percentage is over 50% for all groups, it will then be deleted.")
      # delete the molecule which negative values of replicates for all group are over 50% (all NA.)
      deleted_neg_molec <- negs_all %>% filter_at(sample_raw_list, all_vars(is.na(.))) 
      deleted_molec <- negs_all %>% filter_at(sample_raw_list, all_vars(.==0)) %>% bind_rows(., deleted_neg_molec)
      if(nrow(deleted_molec) != 0){
        deleted <- deleted_molec$LipidMolec %>% unlist() %>% paste0(., ", ", collapse = "") %>% substr(., 1, nchar(.)-2)
        message("\n
Since the lipid molecule ", deleted, " is invalid (all negative values) after background subtraction. 
It will be deleted in the filtered data.")
        # delete corresponding lipid molecules in total data
        filtered_lipidomics <- anti_join(filtered_lipidomics, deleted_molec, by = "LipidMolec") #%>% select(LipidMolec, contains("MainArea"))
        filtered_lipidomics_copy <- anti_join(filtered_lipidomics_copy, deleted_molec, by = "LipidMolec")
        subtracted_data <- filtered_lipidomics
      }
      write_csv(filtered_lipidomics, "data/filtered_to_NA.csv")
      write_csv(filtered_lipidomics_copy, "data/post_filtered.lipids.csv")
      write_csv(subtracted_data, "data/post_subtracted.csv")
      break
    } else{
      option = readline("You typed wrong, please type again, 1/2: ")
    }
  }
}



##########################################################################################
# Quantification of total lipid classes (mean, sd)
##########################################################################################
# making barplot graphs
message("\n Mean value for each group is used for quantification visualization.")
# sum lipid molecues by class
aggregated_class <- filtered_lipidomics %>% 
  select(Class, sample_raw_list) %>% 
  group_by(Class) %>% 
  summarise_at(sample_raw_list, list(~sum(.)))
# get the mean and sd of each group for each class
total_group_mean <- cal_lipid_statistics(aggregated_class, group_info, "mean", "Class")
#class_mean_wide <- total_group_mean[[1]]
class_mean_long <- total_group_mean[[2]]
total_group_sd <- cal_lipid_statistics(aggregated_class, group_info, "sd", "Class")
class_sd_long <- total_group_sd[[2]]
# bind mean and sd for group data
total_data <- left_join(class_mean_long, class_sd_long, by = c("Class", "Groups")) %>% 
  select(Class, mean, sd, Groups) 
#axis_set <- total_data %>% filter_at(vars(mean), any_vars(.<0)) %>% nrow()
# total class lipids of unprocessed raw data visualization

#paras1 <- list("Class", "mean", "Groups", "sd", scale.params)
paras1 <- c("Class", "mean", "Groups", "sd")
total_data$Groups <- factor(total_data$Groups, levels = rev(group_names))
total_plot <- plot_all(total_data, paras1, se=TRUE) + 
  geom_bar(position="dodge", stat="identity") + 
  labs(title = "Total lipid classes", x = "Lipid Classes", y= "AUC (Area under curve)",
       fill = "Experiment Groups") +
   coord_flip()  + 
  guides(fill = guide_legend(reverse = TRUE)) 

print(total_plot)
# save filtered total class lipid abundance in plot total.class.png
ggsave(filename="total.class.png", path = 'plot/', device = "png")



##########################################################################################
# Quantification of individual lipid class
##########################################################################################
# get mean data of each group for each lipid molecule
# you can also use median
each_group_mean <- cal_lipid_statistics(filtered_lipidomics, group_info, "mean", "LipidMolec")
lipid_mean_wide <- each_group_mean[[1]]
lipid_mean_long <- each_group_mean[[2]]
each_group_sd <- cal_lipid_statistics(filtered_lipidomics, group_info,  "sd", "LipidMolec")
lipid_sd_wide <- each_group_sd[[1]]
lipid_sd_long <- each_group_sd[[2]]
# merge mean, sd value and reformat data
each_class <- left_join(lipid_mean_long, lipid_sd_long, by = c("Class", "Groups")) %>% 
  select(Class, mean, sd, Groups) 
message("\nQuantification analysis for individule lipid class")
par_eachclass <- c("LipidMolec", "mean", "Groups", "sd")
# par_eachclass <- list("LipidMolec", "mean", "Groups", "sd", scale_params)
lipidmolecNO_max <- filtered_lipidomics %>% group_by(Class) %>% tally() %>% select(n) %>% unlist() %>% max()
# setting plot labs titles
labs1 <-  labs(x="Acyl composition", y="Main Area", 
               caption="Error bar is the standard deviation for each class in each group", fill = NULL)

#scale.params <- alist(trans = "identity",breaks = waiver(),labels = waiver(), expand = c( 0, 0, 0.2, 0))

# setting the plot limits when the bar numbers exceed the threshold nbar
nbar <- 70    # estimation of threshold which can be modified and at least bigger than group number
post_name <- ""
pars <- list(nbar, ngroups, par_eachclass, plot_all, post_name, labs1)
message("\nEach plot is split no more than ", nbar, " bars for display")
EachClassPlot(each_class, pars)

# overview of each class plot 
nbar <- lipidmolecNO_max
post_name <- "all"
pars <- list(nbar, ngroups, par_eachclass, plot_all, post_name, labs1)
message("\nAlternative display quantification of individule lipid class (all lipids in a class in the same png)")
EachClassPlot(each_class, pars)


##########################################################################################
dev.off()
options(device = "RStudioGD")                                                                                                                                       
##########################################################################################


##########################################################################################
# visualization of lipid class data by pick median value
##########################################################################################
message("Please note that all the lipid molecules contain negative or 0 value in control group
will be deleted for fold change analysis since the data below is not imputed, ")
# visualize for samples
# visualization repliates distribution among all lipid classes
vs_dt <- filtered_lipidomics %>% 
  select(LipidMolec, Class, sample_raw_list) 
vs_sumdt <- vs_dt %>% group_by(Class) %>% summarise_at(vars(sample_raw_list), list(~sum(.)))
sum_dt <- cal_lipid_statistics(vs_sumdt, group_info, "median", "Class")
sum_dt_wide <- sum_dt[[1]]
# remove lipid classes which contains negative values
# sum_wide <- sum_dt_wide %>% group_by(Class) %>% filter_all(all_vars(.>=0))
message("Please type control group name among all experiment groups for total visualization")
control <- check_group(group_names, "control")
control_name <- addquotes(!!control,"_median" )

# remove lipid classes which control group value is 0 
sum_wide <- sum_dt_wide
sum_wide <- sum_wide %>% filter_at(vars(!!control_name), any_vars(.>0))
conserved_class <- vs_sumdt %>% filter(Class %in% sum_wide$Class) 
sum_ctr <- sum_wide %>% select(Class, contains(!!control))
dt <- left_join(conserved_class, sum_ctr, by = "Class")
deleted_classes <- vs_sumdt %>% filter(!Class %in% conserved_class$Class) %>% select(Class) %>% unlist()

# show the deleted lipid class name if it is negative or 0 in control group
if(length(deleted_classes)>=1){
  deleted_classes <- paste0(deleted_classes, " ,", collapse = "") %>% substr(., 1, nchar(.)-1)
  message("Lipid Class ", deleted_classes, "is/are deleted for fold change analysis.")
}

# normalize each value of replicates by median value in control group 
class_nm <- norm_by_mean2(dt, "Class", sample_raw_list, control_name, "n")

# reformat the data
class_long <- class_nm %>%  gather("SAMPLES", "VALs", -Class)
class_long <- class_long %>% 
  rowwise() %>%  
  mutate(GROUPS = case_when(SAMPLES %in% group_info$samples ~ 
                              group_info[samples==SAMPLES, ]$groups, TRUE~"NA"))
class_long$SAMPLES <- class_long$SAMPLES %>% str_remove_all(., ".*\\[") %>% str_remove_all(., "\\]")
class_long$GROUPS <- factor(class_long$GROUPS, levels = group_names) 
# plot relative fold change of lipid class (normalized by median), dotplot
pd <- plot_all(class_long, c("GROUPS", "VALs")) +
  geom_dotplot(aes(fill = GROUPS), binaxis = "y",         # which axis to bin along
               stackdir='center', dotsize = 1, position = "dodge") +
  stat_summary(fun.y = median, geom = "point",  size = 5, color = "red", shape = 95, alpha = 0.8, stroke = 3) +
  facet_wrap(~Class, scales = "free") +
  #scale_fill_manual(values = wes_palette("GrandBudapest2")) +
  scale_fill_startrek() +
  theme_bw() +
  set_theme(theme_params = list(axis.text.x =  element_text(angle=45, hjust=1))) +
  scale_y_continuous(expand = c(0, 0, 0.1, 0)) +
  expand_limits(y = 0) +
  geom_text_repel(aes(label = SAMPLES), size = 2) +
  #scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  labs(title = "Normalized by median", caption = "each dot represent the replicate of one lipid class", 
       x = "Experiment Groups", y = "Relative Fold Change")
print(pd)
ggsave(filename = "class_median_dot.png", path="plot/", device = "png")
ggsave(filename = "class_median_dot.svg", path="plot/")

# plot relative fold change of lipid class (normalized by median), boxplot
pb <- plot_all(class_long, c("Class", "VALs", "GROUPS")) +
  geom_boxplot(outlier.shape = NA, fatten = 1) +
  facet_wrap(~Class, scales = "free") +
  labs(title = "Lipid class data (normalized by median)", x = "Lipid Class", 
       y = "Relative Foldchange", fill = "Experiment Groups") +
  scale_fill_manual(values = wes_palette("GrandBudapest2")) +
  theme_bw() +
  set_theme() +
  expand_limits(y = 0) +
  scale_y_continuous(expand = c(0, 0, 0.1, 0))  +
  scale_fill_d3()
print(pb)
ggsave(filename = "class_median.png", path="plot/", device = "png")




##########################################################################################
# visualization of lipid molecule data, normalized by mean value
##########################################################################################
message("\nFold change analysis for lipid molecules.")


filtered_dt <- filtered_lipidomics%>% mutate_at(sample_raw_list, list(~ifelse(.>= 0, ., NA)))
dt <- filtered_dt  %>%  cal_lipid_statistics(., group_info, "mean", "LipidMolec")
dts <- filtered_dt %>% ungroup() %>% select(LipidMolec, Class, sample_raw_list) 
message("\nPlease type the name of the control group.")
control <- check_group(group_names, "control")
control_name <- paste0(control, "_mean")
dt_control <- dt[[1]] %>% select(Class, contains(control)) %>% rename(LipidMolec = Class)
dt_merged <- left_join(dts, dt_control)
# filter 0 or negative in control group
dt_all <- dt_merged %>% filter_at(vars(all_of(control_name)), any_vars(.>0))
# # transform negative data into NA
# dt_NA <- dt_all  %>% mutate_at(sample_raw_list, list(~ifelse(.>= 0, ., NA)))

dt_nm <- norm_by_mean2(dt_all, "LipidMolec", sample_raw_list, control_name, "n") %>% ungroup()
#dt_nm2 <- dt_nm %>% filter_all(all_vars(!is.na(.)))

molec <-  cal_lipid_statistics(dt_nm, group_info, method, "LipidMolec")

molec_long <- molec[[2]] %>% filter(!is.na(mean)) %>% select(-TYPE)
#molec_long$Groups <- factor(molec_long$Groups, levels = group_names) 
# customize times of comparison  
message("\nHow many fold change analysis for groups you want to do?")
fc_times <- readline("Please input the number: ")
for(i in 1:fc_times){
  message("\nPlease type the name of the control group.")
  control <- check_group(group_names, "control")
  message("\nPlease type the name of the contrast group(s).")
  contrasts <- check_group(group_names, "contrast")
  selected_dt <- molec_long %>% filter(Groups %in% c(control, contrasts))
    par_eachclass <- c("LipidMolec", "mean", "Groups", "sd")  
    labs<- labs(x="Acyl composition", y="Main Area", title = "normalized fold change",
                caption="Error bar is the standard deviation for each class in each group", fill = NULL)
  post_name <- paste0("fc.", i, control, "_", contrasts)
  #par_eachclass <- list("LipidMolec", "mean", "Groups", "sd", scales)  
  pars <- list(nbar, ngroups, par_eachclass, plot_fc, post_name, labs)
  EachClassPlot(selected_dt, pars)

}


##########################################################################################
dev.off()
options(device = "RStudioGD")                                                                                                                                       
##########################################################################################




molecules <- dt_nm  %>%  gather("SAMPLES", "VALs", -LipidMolec)
molecules <- molecules%>% 
  rowwise() %>%  
  mutate(GROUPS = case_when(SAMPLES %in% group_info$samples ~ 
                              group_info[samples==SAMPLES, ]$groups, TRUE~"NA"))
molecules$SAMPLES <- molecules$SAMPLES %>% str_remove_all(., ".*\\[") %>% str_remove_all(., "\\]")
molecules <- molecules %>% separate(., LipidMolec, into=c("Class", "LipidMolec"), sep = "\\(")
molecules$LipidMolec <- molecules$LipidMolec %>% str_remove_all(., "\\)")
molecules$GROUPS <- factor(molecules$GROUPS, levels = group_names)

# violin plot
pv <- ggplot(molecules, aes(x = GROUPS, y = VALs)) +
  geom_violin(aes(fill = GROUPS), trim = FALSE)+
  #geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  # geom_dotplot(aes(fill = GROUPS), binaxis = "y",         # which axis to bin along
  #   #              stackdir='center', position = "dodge") +
  geom_jitter(size = 0.1,  alpha = 0.4, width = 0.1) +
  # stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
  #              geom = "crossbar", width = 0.3, alpha = 0.4, position = "dodge") +
  #stat_summary(fun.y = mean,  geom = "line",  alpha = 0.4, position = "dodge") +
  stat_summary(fun.y=eval(expr(method)), geom="point", na.rm = TRUE,
               size=5, color = "red", shape = 95, alpha = 0.5) +
  facet_wrap(~Class, scales = "free") +
  theme_bw()+
  set_theme(theme_params = list(axis.text.x =  element_text(angle=45, hjust=1))) +
  # theme(#panel.grid  = element_blank(),
  #   panel.border = element_blank(),
  #   axis.line = element_line(colour = "black",
  #                            size = .2,
  #                            #arrow = arrow(length = unit(0.06, "cm"), type = "closed")
  #   ),
  #   title = element_text(size = 10),
  #   strip.background = element_blank(),
  #   strip.text = element_text(size = 6),
  #   axis.title = element_text(size = 8),
  #   axis.text.y = element_text(size = 6),
  #   axis.text.x = element_text(angle = 45, size = 6, hjust = 1),
  #   legend.title = element_text(size = 5),
  #   legend.text = element_text(size = 6)) +
  expand_limits(y = 0) +
  scale_y_continuous(expand = c(0, 0, 0.2, 0)) +
  labs(title = paste0("Individual Lipidmolecules (normalized by ", method, ")"),
       caption = paste0("each dot represent a replicate of one lipid molecule\nThe red line represent ", method),
       x = "Experiment Groups", y = "Relative Foldchange") +
  #scale_fill_manual(values = wes_palette("GrandBudapest2"))
scale_fill_d3()
print(pv)
ggsave(filename = "molec_log_mean.png", path="plot/", device = "png")
ggsave(filename = "molec_log_mean.eps", path="plot/")



##########################################################################################
# break chains
##########################################################################################
fas <- filtered_dt %>% mutate(patterns = str_remove_all(LipidMolec, ".*\\(") %>% 
                                        str_remove_all(., "\\)")) %>% 
  separate(patterns, c("FA1", "FA2", "FA3"), sep = "\\/") %>% 
  select(Rej, LipidMolec, FA1, FA2, FA3, colnames(filtered_lipidomics)[3:ncol(filtered_lipidomics)])
write_csv(fas, "data/fa_chains.csv")
# extract auc information
fas_chains <- fas %>% 
  select(Class, FA1, FA2, FA3, contains("MainArea")) %>% 
  gather(TYPE, FA, -c(Class, contains("MainArea"))) %>% 
  select(-TYPE, Class, FA, contains("MainArea")) %>% 
  filter(!is.na(FA)) %>% 
  group_by(Class, FA) 
write_csv(fas_chains, "data/fa_chains_AUC.csv")
fa_count <- fas_chains %>%   count(FA) 
fas_class <- fas_chains %>% summarise_all(sum) %>% left_join(fa_count, .)
write_csv(fas_class, "data/class_chain_auc.csv")

method <- readline("Please type median/mean of normalization method for total lipid class fold change analysis: ") %>% str_to_lower()

# find control
control <- check_group(group_names, "control")
control_name <- paste0(control, "_", method)
# cal mean 
dt <- fas_class %>% ungroup() 
dt_mean <- cal_lipid_statistics(dt, group_info, method, "Class")
dt_control <- dt_mean[[1]] %>% select(contains(control))
dt_merged <- cbind(dt, dt_control) 
# filter NA in control group
dt_merged <- dt_merged %>% filter(!is.na(dt_control))

# normalized by control group
dt_nm <- norm_by_mean2(dt_merged, "Class", sample_raw_list, control_name, "n") %>% ungroup()
molec_mean <-  cal_lipid_statistics(dt_nm, group_info, method, "Class")
molecules_long <- molec_mean[[2]]
if(method == "mean"){
  molec_sd <- cal_lipid_statistics(dt_nm, group_info, "sd", "Class")
  molecules_long <- left_join(molec_long, molec_sd[[2]])
}


# reformat the data
ds <- dt_merged %>% 
  select(FA) %>% 
  cbind(., dt_nm)  %>% 
  gather("SAMPLES", !!method, -c("Class", "FA"))

ds <- ds %>% filter(!is.na(mean))
ds<- ds %>% 
  ungroup() %>% 
  rowwise() %>%  
  mutate(Groups = case_when(SAMPLES %in% group_info$samples ~ 
                              group_info[group_info$samples==SAMPLES, ]$groups, TRUE~"NA"))
ds$SAMPLES <- ds$SAMPLES %>% str_remove_all(., ".*\\[") %>% str_remove_all(., "\\]")
ds$Groups <- factor(ds$Groups, levels = group_names)

### statistics 
# ymax <- molecules_long %>% filter(!is.na(mean)) %>% mutate(sum = mean + se) %>% group_by(Class) %>% summarise_at(vars(sum), list(~max(.)))
# ypos <- ymax %>% mutate(y1 = sum + 0.05, y2 = sum + 0.05) %>% unite(labels, y1, y2, sep = ",")
# ylabels <- ypos$labels %>% str_split(., ",") %>% lapply(., as.numeric) %>% unlist() %>% round(., digits = 2)
# 
# stat_test <-  compare_means(mean~Groups,  data = ds, ref.group = "control",
#                             method = "wilcox.test", paired = FALSE, group.by = "Class",
# )
# 
# stat_test <- stat_test %>% mutate(y.position = ylabels)
# my_comparisons = rep(list(c("het", "control" ), c("ko","control")), length(unique(molec_long$Class)))
# ggbarplot(ds, x = "FA", y = "mean",
#           fill = "Groups",
#           add = "mean_se", add.params = list(group = "Groups"), 
#           position = position_dodge(0.6),error.plot = "upper_errorbar") +
  # stat_pvalue_manual(
  #   stat_test, x = "Class",
  #   label = "p.signif",
  #   xmin = "group2", xmax = "group1",
  #   position = position_dodge(0.6),
  #   size = 3, 
  # ) +
  # labs(title = "PSs", y="Fold Change", x = "Fatty Acids") +
  # theme(axis.text.y = element_text(angle = 30),
  #       legend.position = "right") +
  # scale_fill_manual(values = wes_palette("GrandBudapest2")) +
  # scale_y_continuous(expand = c(0, 0, 0.2, 0)) +
  # #stat_compare_means(comparisons = my_comparisons) +
  # coord_flip() 

classes <- ds$Class %>% unique()
ds$Groups <- factor(ds$Groups, levels = rev(group_names))
if(method == "mean"){
  for(i in classes){
    fa_chain <- ds %>% filter(Class == i)
    name <- paste0(i, " Relative fold change (Normaized by ", method, ")")
    file1 <- paste0("plot/", i, ".fc.png")
    file2 <- paste0("plot/", i, ".fc.svg")
    p <-  ggbarplot(fa_chain, x = "FA", y = "mean",
                    fill = "Groups",
                    add = "mean_se", add.params = list(group = "Groups"),
                    position = position_dodge(0.6),error.plot = "upper_errorbar") +
      labs(title = name, y="Fold Change", x = "Fatty Acids") +
      theme(axis.text.y = element_text(angle = 30),
            legend.position = "right") +
      #scale_fill_manual(values = wes_palette("GrandBudapest2")) +
      scale_fill_simpsons() +
      # add_scales()
      scale_y_continuous(expand = c(0, 0, 0.2, 0)) +
      guides(fill = guide_legend(reverse=FALSE)) +
      coord_flip()
    print(p)
    ggsave(file1, device = "png")
    ggsave(file2, device = "svg")
  }
}else{
  for(i in classes){
    fa_chain <- ds %>% filter(Class == i)
    name <- paste0(i, " Relative fold change (Normaized by ", method, ")")
    file1 <- paste0("plot/", i, ".fc.png")
    file2 <- paste0("plot/", i, ".fc.svg")
    p <-  ggbarplot(fa_chain, x = "FA", y = "median",
                    fill = "Groups",
                    add = "median", add.params = list(group = "Groups"),
                    position = position_dodge(0.6),error.plot = "upper_errorbar") +
      labs(title = name, y="Fold Change", x = "Fatty Acids") +
      theme(axis.text.y = element_text(angle = 30),
            legend.position = "right") +
      #scale_fill_manual(values = wes_palette("GrandBudapest2")) +
      scale_fill_simpsons() +
      # add_scales()
      scale_y_continuous(expand = c(0, 0, 0.2, 0)) +
      guides(fill = guide_legend(reverse=FALSE)) +
      coord_flip()
    print(p)
    ggsave(file1, device = "png")
    ggsave(file2, device = "svg")
  }
}




# heat map   

































































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
message("\n")
ncomparisons <- readline("How many volcano plots to generate: ") %>% as.numeric()
impute_dt <- impute_not("y", filtered_lipidomics, sample_raw_list)
imputated_lipids <- impute_dt[[2]]
rownames(imputated_lipids) <- imputated_lipids$LipidMolec
processed_lipids <- imputated_lipids %>% select(-Class, -LipidMolec)

# Fit model and extract contrasts 
fit <- lmFit(processed_lipids, design)
# plot volcano graph
VolPlot(ncomparisons, fit, processed_lipids)

