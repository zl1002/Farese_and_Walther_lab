# Script: FWL_lipidomics_QC_3.2.2.R
# Author: Wenting 
# Notes:  To start, typing command in the console-----> source("FWL_lipidomics_3.2.2.R") 
#         or press the source butthon. 
#         Please make sure Mac users installed XQuartz.
#         This script is designed as a reference for lipidomics experiment. 
#         It is based on Niklas and Kenny's previous work (Their work files can be found in folders 
#         quality_control and statistics_quantification). Acknowledge to Weng's technical 
#         guidance, Laura's fatty acid saturation analysis project and Sebstian's shadow experiment
###########################################################################################
#

#########################################################################################
###### read data
###
# if you want to run the QC independently, please unquote the part below
##########################################################################################
# # source function from FWL_lipidomics.**.functions.R script
# source("FWL_lipidomics_3.2.2.FUNCTIONS.R")
# # Check directory existence, if not then make one
# # all the plots will stored in plot directory and data in the data directory
# dirs <- c("plot", "data", "plot/classes")
# mkdirs(dirs)
# 
# # read the file from csv directory
# csv_files <- list.files(path = "converted", pattern = "\\.csv$")
# csv_list <- c()
# for(i in seq_along(csv_files)) csv_list[i] <- addquotes(!!as.character(i), " ", !!csv_files[i], "\n")
# message("The following files had been generated. Select ONE for subsequent the list of file names:\n", csv_list)
# file_option <- readline("Please input the index of the file: ") %>% as.numeric()
# # the file path and name
# target_file <- paste("converted/", csv_files[file_option], sep = "")
# print(target_file)
# lipidomics <- read_csv(target_file, col_types = cols())


##########################################################################################
# check experiment samples
##########################################################################################
# check if there is internal control samples or samples don't use for analysis
# and extract information of different grades and p value 
# which are less or equal than 0.001 for each molecule and store this info in the variable lipid_count
message("\nDo you have samples excluded for analysis which including experiment controls like internal standards")
delete_option <- retype_choice("Y/N")

#colnames(lipidomics) <- colnames(lipidomics) %>% ifelse(str_detect(., "GroupArea"), str_replace_all(., "GroupArea", "MainArea"), .)
#colnames(lipidomics) <- colnames(lipidomics) %>% ifelse(str_detect(., "Grade"), str_replace_all(., "Grade", "MainGrade"), .)



lipid_check <- delete_samples(delete_option, lipidomics)
# counter of Grades A and B, filter P value less than 0.001 for LipidSearch
lipid_count <- lipid_check[[1]]
# remove samples
lipid_remove <- lipid_check[[2]]
# add this count into original strucutes
lipid_select <- lipid_remove %>% bind_cols(lipid_count)  

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
k <- readline("Minimum number of identified molecules (Grade A+B) required in all samples, n>=3: ") %>% as.numeric(.)
j <- readline("Minumum number of significantly identified lipids p<=0.001 in all samples, n>=3: ") %>% as.numeric(.)


filtered_lipidomics1 <- lipid_select %>% 
  rowwise() %>% 
  filter(Rej == 0 & sum(A, B) >= k & APValue.001 >= j) %>% 
  #filter(Rej. == 0 & sum(A, B) >= k ) %>% 
  as.data.frame(.)
# the filtered raw data is stored in data/filtered.raw.data.csv
write.csv(filtered_lipidomics1, "data/QC/filtered.raw.data.csv")



##########################################################################################
# Fix duplicated lipid molecules
##########################################################################################
# detect same lipid molecules with different retention time
# Filter the lipid molecule contains same name but different retention time based on your criteria  
# all the dupliated molecules will be stored in dupliates.csv
duplicate_molecs <- detect_duplicates(filtered_lipidomics1)
# if move on, fix method for duplicated molecules
paras <- c("Class", "LipidMolec", "BaseRt", "MainIon")

filtered_lipidomics2  <- filter_duplicate(duplicate_molecs, filtered_lipidomics1, paras)




#filtered_lipidomics1 <- filtered_lipidomics1 %>% separate(., LipidIon, into=c("LipidMolec", "MainIon"), sep = "\\)")

##########################################################################################
# check background information
##########################################################################################
# plot the background information of the blank sample c

blank_sample <-filtered_lipidomics2 %>% 
  select(Class, contains("MainArea[c]")) %>% 
  group_by(Class) %>% 
  summarise_all(.funs=sum)
parameters_bw <- colnames(blank_sample) <- c("class", "value")

# plot the blank sample
p1 <- plot_all(blank_sample, parameters_bw) +
  geom_bar(position="dodge", stat="identity") +
  labs(x = "Lipid Classes", 
       y = (bquote("Fold change: " ~ log[10]~"("~AUC~")")), 
       title = "Area under Curve for blank sample (background) in each class") +
  # scale_y_continuous(trans = pseudo_log_trans(base = 10),
  #                    breaks = c(0, 10^3, 10^6, 10^9, 10^12),
  #                    labels = c(0, bquote(10^3), bquote(10^6), bquote(10^9), bquote(10^12)),
  #                    expand = c(0, 0, 0.1, 0)) +
  add_scales() +
 # scale_y_continuous(labels = scientific_format(), expand = c(0,0, 0.1, 0)) +
  coord_flip() 
print(p1)
# the blank sample plot is stored in background.png
ggsave(filename = "background.png", path="plot/QC/", device = "png")


##########################################################################################
# lipid class summary
##########################################################################################
message("\nThe info below and summary plot will show the summary information of classes after filtering the data")
# two methods check how many lipids passed filtering
print(describe(filtered_lipidomics2$Class))
filtered_lipidomics2 %>% group_by(Class) %>% summarise(lipid_class_num = n()) %>% select(Class, lipid_class_num) %>% arrange(lipid_class_num)  %>% formattable(.)

# get lipid class proportion information
prop_data <- filtered_lipidomics2
class(prop_data$Class) <- factor(prop_data$Class)
prop_data <- prop_data %>% 
  group_by(Class) %>% 
  summarise(count = n()) %>% 
  mutate(prop = count/sum(count))

# plot lipid class proportion information
p2 <- plot_all(prop_data, c("Class", "prop")) +
  #  ggplot(prop_data, aes(x = reorder(Class, prop), y = prop)) +
  theme_bw() + 
  set_theme() +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels=scales::percent, expand = c(0, 0, 0.2, 0)) +
  ggtitle("Class Summary plot: \n Proportion of different class numbers among all the samples")+
  xlab("classes")+ ylab("Relative frequencies (%)")  +
  coord_flip() + 
  scale_fill_grey()
print(p2)
ggsave(filename = "prop_summary.png", path = 'plot/QC/', device = "png")

# save lipid class summary information into file proportion_classes.csv in data directory
prop_data %>%
  mutate(prop = percent(prop*100/100)) %>%
  rename(number_of_lipid_molec = count, proportion = prop)%>% 
  write_csv(., "data/QC/proportion_classes.csv")


##########################################################################################
# retention time analysis
##########################################################################################
###  Abundance vs. retention time for all samples
retention_data <- filtered_lipidomics2 %>% 
  select(contains("MainArea[s"), BaseRt, Class) %>% 
  gather(sample, MainArea, -c(BaseRt, Class))
n_classes <- retention_data$Class %>% unique() %>% length()
pars <- c("MainArea", "BaseRt", n_classes)
rt_plot <- plot_rt_allsamples(retention_data, pars)
print(rt_plot)
ggsave(filename = "all_retention.png", path = 'plot/', device = "png", 
       width = 10, height = 8, dpi = 150, units = "in")



##########################################################################################
# mark odd chains
##########################################################################################
# mark classes which contain the odd chains classes 
odd_index <- filtered_lipidomics2$LipidMolec %>%
  str_locate(., "(\\d[13579]:)|([13579]:)")
odd_chains <- filtered_lipidomics2[unique(which(!is.na(odd_index), arr.ind=TRUE)[,1]), ]
# store odd chains in file odd_chains
write_csv(odd_chains, "data/QC/odd_chains.csv")
percent_odd <- nrow(odd_chains)/nrow(filtered_lipidomics2) 
percent_odd <- percent(percent_odd*100/100)
# get odd chain percentage in the raw data
message("\nThere are ", nrow(odd_chains), " lipid molecules contain odd chains. ",
        "\nThe odd chain of fatty acids percent is ", percent_odd, " in total.")

message("The odd chain information is stored in odd_chains.csv.")

##########################################################################################
# test standard TG(17:1/17:1/17:1) abundance in all samples
# Please note that this part could be editted for other funture standards
###########################################################################################
#Please have a look and decide if you wanna delete any lipid molecules. If you do, just delete the corresponding lipid in TG.odd.csv. The pipeline will take care of the rest
TGs <- odd_chains %>% filter(Class=="TG") %>% arrange(FA)
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
    add_scales() +
    coord_flip()
  print(tp17)
  ggsave(filename = "TG17_all.png", path="plot/", device="png")
}















