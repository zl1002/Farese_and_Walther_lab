####################################################################################
# Script: FWL_FattyAcids_3.1.R
# Author: Wenting Lyu
# Notes: This script assist executing for main script FWL_lipidomics_3.1.R which
#         helps generating the graph and data for the workflow of lipidomics.
#         it also calculate the SFA, MUFA, PUFA and the aggregated main area 
#         for each sample.
#         First, Make sure that your R and Rstudio are newly enough for installing the 
#         packages needed in the script. Otherwise the script will pop out warnings 
#         for packages and won't run.
#         Second, typing command in the console-----> source("FWL_lipidomics_3.1.R")
#         or press the source button.
#         Third, users can independently implement this analysis by running 
#         "fattyAcidsaturation_analysis2.1.r" in directory fattyAcids_saturation_analysis.
#         This script is derived from Laura's project
#####################################################################################
data <- cal_sample_saturation(filtered_lipidomics,  group_info)

data <- data %>% filter_all(any_vars(!is.na(.)))

# find all groups which mean or median value is 0
empty_group <- data %>% 
  select(-contains("MainArea"), -contains("se")) %>% 
  group_by(Class) %>% 
  filter_all(any_vars(.<=0))
write_csv(empty_group, "data/empty_groups.csv")


# normalized the data for visualization
############
control <- check_group(group_names, "control") 
control_info <- addquotes(!!control, "_mean")
control_data <- pick_control(data, "mean", group_info, control)
post_fix <- c("_SFA", "_MUFA", "_PUFA")
control_names <- c()
samples <- group_info$samples %>% unlist()
sample_list <- list()
for(j in post_fix) sample_list[[j]] <- addlists(samples, j) 

# find non detected lipid class which control group mean is 0 
non_detected_lipid <- control_data %>% 
  select(Class, contains(!!control_info)) %>% 
  mutate(control_group_mean = select(., contains(control_info)) %>% 
           rowSums()) %>% 
  filter(control_group_mean == 0)
if(nrow(non_detected_lipid) != 0){
  message("\nNon detected lipid classes (0 value) are stored under non_detected_class.csv.")
  write_csv(non_detected_lipid, "data/non_detected_class_control.csv")
}

# find fatty acids of control group which mean is 0, and store it as non_detected_FA_control.csv
for(i in post_fix) control_names[i] <- addquotes(!!control_info, !!i)
non_detected <- control_data %>% 
  select(Class, contains(!!control_info))  %>% 
  filter_all(any_vars(.==0))
message("\nNon detected fatty acid types (0 value) are stored under non_detected_class.csv.")
write_csv(non_detected, "data/non_detected_FA_control.csv")

# get data which class are detected
detected_data <- anti_join(control_data, non_detected, by="Class")

# pick classes you don't want to display in the plot
message("\n")
delete_not <- readline("Any classes of lipid to be excluded from analysis? (Y/N): ") %>% str_to_lower()
if(delete_not != "n"){
  delete_class <- readline("Which class of lipid to be excluded? e.g. TG Cer: ")
  data_copy <- data # data_copy stores the unfiltered data
  data <- data %>% filter(!Class %in% delete_class)
  detected_data <- detected_data %>% filter(!Class %in% delete_class)
} 

group_level <- group_names
legend_level <- c("SFA", "MUFA", "PUFA")

# reformat the data for stack plot
dtFA_mean <- FormatData2(data, "mean")
dtFA_se <- FormatData2(data, "se")
dtFA <- left_join(dtFA_mean, dtFA_se, by = c("Class", "TYPE", "Groups"))
dtFA_mean$Groups <- factor(dtFA_mean$Groups, levels = group_level)
dtFA$Groups <- factor(dtFA$Groups, levels = group_level)
dtFA$TYPE <- factor(dtFA$TYPE, levels= legend_level)
write_csv(dtFA, "data/fa_mean.csv")

dtFA_median <- FormatData2(data, "median")
dtFA_median$Groups <- factor(dtFA_median$Groups, levels = group_level)
dtFA_median$TYPE <- factor(dtFA_median$TYPE, levels= legend_level)
write_csv(dtFA_median, "data/fa_median.csv")


# normalized all the samples by the mean of control group
samples_SFA <- norm_by_mean(detected_data, "_SFA", sample_list[[1]], control_names[1])
samples_MUFA <- norm_by_mean(detected_data, "_MUFA", sample_list[[2]], control_names[2])
samples_PUFA <- norm_by_mean(detected_data, "_PUFA", sample_list[[3]], control_names[3])



# store samples which are normalized by control mean
norm_samples <- left_join(samples_SFA, samples_MUFA, by = "Class") %>% 
  left_join(., samples_PUFA, by = "Class")
write_csv(norm_samples, "data/normalized_samples.csv")

# get detected data 
#detected <- data[which(data$Class %in% detected_data$Class), ]

# get fold change information
data_mean <- FC_fun(group_info, "mean", norm_samples)
data_se <- FC_fun(group_info, "se", norm_samples)
data_mean_se <- cbind(Class = norm_samples$Class, data_mean, data_se)
write_csv(data_mean_se, "data/fc_data.csv")


# reformat data for ploting
dtFA_mean_long <- FormatData2(data_mean_se, "mean")
dtFA_se_long <- FormatData2(data_mean_se, "se")
dtFA_long <- left_join(dtFA_mean_long, dtFA_se_long, by = c("Class", "TYPE", "Groups"))
dtFA_long$Groups <- factor(dtFA_long$Groups, levels = group_level)
dtFA_long$TYPE <- factor(dtFA_long$TYPE, levels= legend_level)
write_csv(dtFA_long, "data/fa_mean_long.csv")


plot_fas <- function(dtFA_mean, dtFA_median, dtFA_long, dtFA){
  options(warn=-1)
  # stack plots based on mean
  pars1 <- c("Groups", "mean", "TYPE")
  p1 <-  plot_all(dtFA_mean, pars1) +
    geom_bar(stat = "identity") + 
    facet_wrap(~Class, scales = "free") +
    scale_y_continuous(labels = label_scientific(digits = 2), expand = c(0, 0, 0.2, 0)) +
    labs(x = "Fatty Acids types", y = "Main Area value (AUC)", 
         title ="mean based stackplot", fill = NULL,
         caption = "error bar is standard error") +
    scale_fill_jco() +
    theme(
   axis.text.x  = element_text(angle = 30, hjust = 1, size = 6)
    )
  print(p1)    
  ggsave(filename = "meanBased.stackplots.png", path="plot/", device="png")
  ggsave(filename = "meanBased.stackplots.svg", path="plot/")
  
  # stack plots based on median
  pars2 <- c("Groups", "median", "TYPE")
  p2 <- plot_all(dtFA_median, pars2) +
    geom_bar(stat = "identity") + 
    facet_wrap(~Class, scales = "free") +
    scale_y_continuous(labels = label_scientific(digits = 2), expand = c(0, 0, 0.2, 0)) +
    labs(x = "Fatty Acids types", y = "Main Area value (AUC)", 
         title ="median based stackplot", fill = NULL,
         caption = "error bar is standard error") +
    scale_fill_simpsons() +
    theme(
      axis.text.x  = element_text(angle = 30, hjust = 1, size = 6)
    )
  print(p2)
  ggsave(filename = "medianBased.stackplots.png", path="plot/", device="png")
  ggsave(filename = "medianBased.stackplots.svg", path="plot/")
  
  # bar plots split by lipid class
  pars3 <- c("Groups", "mean", "TYPE", "se")
  p3 <- plot_all(dtFA, pars3, se = TRUE) +
    geom_bar(stat = "identity",  position=position_dodge()) +
    facet_wrap(~Class, scales = "free") +
    scale_fill_jama() +
   scale_y_continuous(labels = label_scientific(digits = 2), 
                      expand = c(0, 0, 0.2, 0)) +
    labs(x = "Fatty Acids types", y = "Main Area value (AUC)", 
         title = "Mean based data", fill = NULL,
         caption = "error bar is standard error") +
    theme(
      axis.text.x  = element_text(angle = 30, hjust = 1, size = 6)
    ) 
    
  print(p3)
                                                                                                                                                                      ggsave(filename = "meanBased.fattyAcids.png", path="plot/", device="png")
  ggsave(filename = "meanBased.fattyAcids.svg", path="plot/")
  
  # fold change for mean based data
  p4 <- plot_all(dtFA_long, pars3, se = TRUE)  +  
    geom_bar(stat = "identity",  position=position_dodge()) +
    facet_wrap(~Class, scales = "free") +
   scale_y_continuous(expand = c(0, 0, 0.2, 0)) +
    labs(x = "Fatty Acids types", y = "fold change (normalized)", 
         title = "Fold Change plots", fill = NULL,
         caption = "error bar is standard error") +
   scale_fill_nejm()+
    theme(
      axis.text.x  = element_text(angle = 30, hjust = 1, size = 6)
    )
  print(p4)
  ggsave(filename = "meanBased.fc.png", path="plot/", device="png")
  ggsave(filename = "meanBased.fc.svg", path="plot/")
  
  # percentage plots
  p5 <- plot_all(dtFA, pars1) +
    labs(title = "mean based stackplot", x ="Fatty Acids types", y = "Percentage") + 
    geom_bar(stat = "identity", position = "fill") +
    facet_wrap(~Class, scales = "free") +
    scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0, 0, 0)) +
    labs(title = "mean based stackplot", fill = NULL,
         x ="Fatty Acids types", y = "Percentage") +
   # scale_fill_igv()
   scale_fill_manual(values = wes_palette("GrandBudapest2")) +
    theme(
      axis.text.x  = element_text(angle = 30, hjust = 1, size = 6)
    )
 # scale_fill_simpsons()
  print(p5)
  ggsave(filename = "percentage.png", path="plot/", device = "png")
  ggsave(filename = "percentage.svg", path="plot/")
  
}


plot_fas(dtFA_mean, dtFA_median, dtFA_long, dtFA)

# setting order for groups and legends
#ms <- message("Do you need reorder the group names for display? Please enter Y/N: ")
message("\n")
options <- readline("Do you need reorder the group names for display? Please enter Y/N: ") %>% str_to_lower()
if(options == "y"){
  message("Please list the order of group names for diplay from left to right: ")
  group_level <- readline(" e.g. Untreated Palmitate :  ") %>% str_split(., "\\s") %>% unlist()
  dtFA$Groups <- factor(dtFA$Groups, levels = group_level)
  #dtFA$TYPE <- factor(dtFA$TYPE, levels= legend_level)
  dtFA_mean$Groups <- factor(dtFA_mean$Groups, levels = group_level)
  dtFA_median$Groups <- factor(dtFA_median$Groups, levels = group_level)
  plot_fas(dtFA_mean, dtFA_median, dtFA_long, dtFA)
}
#return(dtFA)
#}
