####################################################################################
# Script: FWL_EtherLipids_3.2.0.R
# Author: Wenting Lyu
# Notes: This script assist executing for main script FWL_lipidomics_3.1.R which
#         helps generating the graph and data for the workflow of lipidomics.
#         it also calculate help visualize the ether saturation information
#         First, Make sure that your R and Rstudio are newly enough for installing the 
#         packages needed in the script. Otherwise the script will pop out warnings 
#         for packages and won't run.
#         Second, typing command in the console-----> source("FWL_lipidomics_3.1.R")
#         or press the source button.
#         Third, users can independently implement this analysis by running 
#         "fattyAcidsaturation_analysis2.1.r" in directory fattyAcids_saturation_analysis.
#         This script is derived from Laura's project
#####################################################################################
ether_detection <- filtered_lipidomics %>% rowwise() %>% mutate(ether = ifelse(str_detect(LipidMolec, "\\(.*[ep]"), "YES", "NO"))

# transform negative to NA
ether_detection <- ether_detection %>% mutate_at(vars(all_of(sample_raw_list)), list(~ifelse(.<0, NA, .)))

ether_yes <- ether_detection %>% filter(ether == "YES")
ether_all <- ether_yes %>% 
  select(Class, contains("MainArea"), -"MainArea[c]") %>% 
  group_by(Class)%>% 
  summarise_at(vars(sample_raw_list), list(~sum(., na.rm = TRUE)))%>% 
  gather(SAMPLES, ether_AUC, -Class) 
  

all_lipids <- ether_detection %>% 
  filter(Class %in% unique(ether_all$Class)) %>% 
  select(Class, contains("MainArea"), -"MainArea[c]") %>% 
  group_by(Class)%>% 
  summarise_at(vars(sample_raw_list), list(~sum(., na.rm = TRUE))) %>% 
  gather(SAMPLES, all_AUC, -Class)


ether_percent <- left_join(all_lipids, ether_all)
 
ether_percent <- ether_percent %>% 
  rowwise() %>%
   mutate(rest_AUC = all_AUC - ether_AUC) %>% 
   select(-all_AUC) %>% 
  mutate(GROUPS = ifelse(SAMPLES %in% group_info$samples, 
                         unlist(group_info[group_info$samples==SAMPLES, 2]), "NA")) %>% 

  mutate(SAMPLES = str_remove_all(SAMPLES, "MainArea\\[") %>% str_remove_all(., "\\]")) %>% 
  gather(type, value, -c("Class", "SAMPLES", "GROUPS"))

ether_percent <- ether_percent %>% mutate_if(is.character, as.factor)
ordered_samples<- unique(ether_percent$SAMPLES) %>% 
  str_remove_all(., "s") %>% 
  as.numeric() %>% 
  sort() %>% 
  paste0("s", .)
ether_percent$SAMPLES <- factor(ether_percent$SAMPLES, levels = ordered_samples)
ether_percent$GROUPS <- factor(ether_percent$GROUPS, levels = group_names)








params1 <- c("SAMPLES", "value", "type")
p1 <- plot_all(data = ether_percent, params1) +
  geom_bar(stat = "identity") +
  scale_fill_simpsons(labels = c("ether", "rest_lipid")) +
  facet_wrap(~Class, scales = "free") +
  theme(axis.text.x = element_text(angle = 45, size = 8, hjust = 1),
        axis.line = element_line(size = 0.2)) +
  scale_y_continuous(labels = scientific_format(), expand = c(0, 0, 0.2, 0)) +
  labs(x = "experiment samples", y = "value", title = "ether in lipid class for each sample", fill = "", caption = "This visualization is only for ether lipids here") 
print(p1)
ggsave("plot/Ether/ether.png", device = "png")
ggsave("plot/Ether/ether.svg", device = "svg")

p2 <- plot_all(data = ether_percent, params1) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_simpsons(labels = c("ether", "rest_lipid")) +
  facet_wrap(~Class, scales = "free") +
  set_theme(theme_params = list(axis.text.x = element_text(angle = 45, size = 8, hjust = 1))) +
  scale_y_continuous( expand = c(0, 0, 0.1, 0), labels = scales::percent_format()) +
  labs(x = "experiment samples", y = "value", title = "ether in lipid class for each sample", fill = "", caption = "This visualization is only for ether lipids here") 
print(p2)
ggsave("plot/Ether/ether_percentage.png", device = "png")
ggsave("plot/Ether/ther_percentage.svg", device = "svg")






ether_percent_group <- ether_percent %>% group_by(Class, GROUPS, type)  %>% summarise_at(vars(value), list(~mean(., na.rm = TRUE)))
params2 <- c("GROUPS", "value", "type")
p3 <- plot_all(data = ether_percent_group, params2) +
  geom_bar(stat = "identity") +
  scale_fill_simpsons(labels = c("ether", "rest_lipid")) +
  facet_wrap(~Class, scales = "free") +
  set_theme(theme_params = list(axis.text.x = element_text(angle = 45, size = 8, hjust = 1))) +
  scale_y_continuous(expand = c(0, 0, 0.2, 0), labels = scientific_format()) +
  labs(x = "experiment Groups", y = "value", title = "ether in lipid class for each group", fill = "", caption = "This visualization is only for ether lipids here") 
print(p3)
ggsave("plot/Ether/ether_group.png", device = "png")
ggsave("plot/Ether/ether_group.svg", device = "svg")

p4 <- plot_all(data = ether_percent_group, params2) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_simpsons(labels = c("ether", "rest_lipid")) +
  facet_wrap(~Class, scales = "free") +
  set_theme(theme_params = list(axis.text.x = element_text(angle = 45, size = 8, hjust = 1))) +
  scale_y_continuous(expand = c(0, 0, 0.1, 0), labels = scales::percent_format()) +
  labs(x = "experiment samples", y = "value", title = "ether in lipid class for each group", fill = "", caption = "This visualization is only for ether lipids here") 
print(p4)
ggsave("plot/Ether/ether_group_percent.png", device = "png")
ggsave("plot/Ether/ether_group_percent.svg", device = "svg")


sd_group <- ether_percent %>% group_by(Class, GROUPS, type) %>% select(-SAMPLES) %>% filter(type == "ether_AUC") %>% summarise_at(vars(value), list(~sd(., na.rm=TRUE)))  %>% rename(sd=value)
data <- ether_percent_group %>% filter(type == "ether_AUC") %>% left_join(sd_group) %>% ungroup()
data$GROUPS <- factor(data$GROUPS, levels= rev(group_names))
p5 <- plot_all(data, c("Class", "value", "GROUPS", "sd"), se=FALSE) +
  geom_bar(position=position_dodge(preserve = "single", width = 0.6), stat="identity", width = 0.6 ) + 
  geom_errorbar(aes(x=Class, ymin=value, ymax = value + sd), 
                position=position_dodge(preserve = "single", width = 0.6), 
                stat="identity", width= 0.1, size= 0.2) +
  labs(title = "Total lipid classes", x = "Experiment Groups", y= "AUC (Area under curve)",
       fill = "Experiment Groups") +
  coord_flip() + 
  scale_fill_simpsons() +
  guides(fill = guide_legend(reverse = TRUE)) 
print(p5)
ggsave("plot/Ether/ether_abundance.png", device = "png")  

group_par <- c("LipidMolec", "Class")
data_mean <- cal_lipid_statistics(ether_yes, group_info, "mean", group_par)
ether_mean_long <- data_mean[[2]]
data_sd <- cal_lipid_statistics(ether_yes, group_info,  "sd", group_par)
ether_sd_long <- data_sd[[2]]
# merge mean, sd value and reformat data
each_class <- left_join(ether_mean_long, ether_sd_long) 
each_class$Groups <- factor(each_class$Groups, levels = rev(group_names))
each_class <- each_class %>% ungroup() %>% rowwise() %>% mutate(LipidMolec = str_remove_all(LipidMolec, ".*\\(") %>% str_remove_all(., "\\)")) %>% ungroup()

p6 <- plot_all(each_class, c("LipidMolec", "mean", "Groups", "sd"), se=FALSE) +
  geom_bar(position=position_dodge(preserve = "single", width = 0.6), stat="identity", width = 0.6 ) + 
  geom_errorbar(aes(x=LipidMolec, ymin=mean, ymax = mean + sd), 
                position=position_dodge(preserve = "single", width = 0.6), 
                stat="identity", width= 0.1, size= 0.2) +
  labs(title = "Ether lipid molecules", x = "Lipid Classes", y= "AUC (Area under curve)",
       fill = "Experiment Groups") +
  facet_wrap(~Class, scales="free") +
  coord_flip() + 
  scale_fill_jco() +
  guides(fill = guide_legend(reverse = TRUE)) 
print(p6)
ggsave("plot/Ether/ether_molec_abundance.png", device = "png")  


message("Please type control group name among all experiment groups for total visualization")
control <- check_group(group_names, "control")
control_name <- addquotes(!!control,"_mean" )
control_group <- data_mean[[1]] %>% select(LipidMolec, all_of(control_name)) 
ether_each <- left_join(ether_yes, control_group)
normalized_ether <- norm_by_mean2(ether_each, group_par, sample_raw_list, control_name, "n")
nm_ether_mean <- cal_lipid_statistics(normalized_ether, group_info, "mean", group_par)
nm_ether_sd <- cal_lipid_statistics(normalized_ether, group_info, "sd", group_par)
nm_ether <- left_join(nm_ether_mean[[2]], nm_ether_sd[[2]])
nm_ether <- nm_ether %>% ungroup() %>% mutate(LipidMolec = str_remove_all(LipidMolec, ".*\\(") %>% str_remove_all(., "\\)"))
nm_ether$Groups <- factor(nm_ether$Groups, levels = rev(group_names))

p7 <- plot_all(nm_ether, c("LipidMolec", "mean", "Groups", "sd"), se=FALSE) +
  geom_bar(position=position_dodge(preserve = "single", width = 0.6), stat="identity", width = 0.6 ) + 
  geom_errorbar(aes(x=LipidMolec, ymin=mean, ymax = mean + sd), 
                position=position_dodge(preserve = "single", width = 0.6), 
                stat="identity", width= 0.1, size= 0.2) +
  labs(title = "Ether lipid molecules (Normalized by mean)", x = "Lipid Classes", y= "AUC (Area under curve)",
       fill = "Experiment Groups") +
  facet_wrap(~Class, scales="free") +
  scale_y_continuous(expand = c(0, 0, 0.2, 0)) +
  theme(axis.line = element_line(size = .2)) + 
  coord_flip() + 
  scale_fill_jco() +
  scale_y_continuous(labels = scientific_format(), expand= c(0, 0, 0.2, 0)) +
  guides(fill = guide_legend(reverse = TRUE)) 
print(p7)
ggsave("plot/Ether/normalized_ether_molec_abundance.png", device = "png")  



# althernative method
# test <- cal_lipid_statistics(normalized_ether, group_info, "mean", "LipidMolec")
# test1 <- test[[2]]
# test_sd <- cal_lipid_statistics(normalized_ether, group_info, "sd", "LipidMolec")
# test2 <- test_sd[[2]]
# test_all <- left_join(test1, test2)
# 
# test_all$Groups <- factor(test_all$Groups, levels = rev(group_names))
# plot_all(test_all, c("LipidMolec", "mean", "Groups", "sd"), se=TRUE) +
#   geom_bar(stat = "identity",  position=position_dodge()) +
#   scale_fill_simpsons() +
#   scale_y_continuous(expand = c(0,0, 0.2, 0)) +
#   coord_flip()
# 



# ether lipid
saturation_data <- read_csv("data/Saturation/count_lipid_filtered_lipidomics.csv", col_types = cols())

ether_lipids <- saturation_data %>% rowwise() %>% mutate(ether = ifelse(str_detect(LipidMolec, "\\(.*[ep]"), "YES", "NO"))

ether_dt <- ether_lipids %>% filter(ether == "YES")
write_csv(ether_dt,"data/Ether/ether_lipids_saturation.csv")


all_lipids <- ether_lipids %>% 
  filter(Class %in% unique(ether_dt$Class)) %>% 
  select(Class, contains("MainArea"), -"MainArea[c]") %>% 
  group_by(Class)%>% 
  summarise_all(sum) %>% 
  gather(SAMPLES, all_AUC, -Class)


ether_saturation <- ether_dt %>% filter(!FA_types == "None")

ether_long <- ether_saturation %>% 
  group_by(Class) %>% 
  summarise_at(vars(all_of(sample_raw_list)), list(~sum(., na.rm = TRUE))) %>% 
  gather(SAMPLES, total_value, -Class)




pufa_percent <- ether_saturation %>% 
  mutate("PUFA%" = PUFA/(SFA + MUFA + PUFA))
write_csv(pufa_percent, "data/Ether/pufa_individuleMolec.csv")

sfa_percent <-  ether_saturation %>% 
  mutate("SFA%" = SFA/(SFA + MUFA + PUFA))

mufa_percent <-  ether_saturation %>% 
  mutate("MUFA%" = PUFA/(SFA + MUFA + PUFA))

pufa_value <- pufa_percent %>% transmute_at(vars(all_of(sample_raw_list)), list(~.*`PUFA%`)) %>% bind_cols(Class=pufa_percent$Class, .)
pufa_samples <- bind_cols(LipidMolec = pufa_percent$LipidMolec, LipidIon = pufa_percent$LipidIon, pufa_value)
write_csv(pufa_samples, "data/Ether/pufa_in_each_sample.csv")

pufa_ether <- pufa_value %>% group_by(Class) %>% summarise_at(vars(all_of(sample_raw_list)), list(~sum(.))) %>% arrange(Class)
pufa_long <- pufa_ether %>% gather(SAMPLES, PUFA_value, -Class)

sfa_value <- sfa_percent %>% 
  transmute_at(vars(all_of(sample_raw_list)), list(~.*`SFA%`)) %>% 
  bind_cols(Class=sfa_percent$Class, .)

sfa_long <- sfa_value %>% 
  group_by(Class) %>% 
  summarise_at(vars(all_of(sample_raw_list)), list(~sum(.))) %>% 
  arrange(Class) %>% 
  gather(SAMPLES, SFA_value, -Class)

mufa_value <- mufa_percent %>% 
  transmute_at(vars(all_of(sample_raw_list)), list(~.*`MUFA%`)) %>% 
  bind_cols(Class=mufa_percent$Class, .) 
mufa_long <- mufa_value %>% 
  group_by(Class) %>% 
  summarise_at(vars(all_of(sample_raw_list)), list(~sum(.))) %>% 
  arrange(Class) %>% 
  gather(SAMPLES, MUFA_value, -Class)

#ether_fa <- left_join( pufa_long, sfa_long, by = c("Class", "SAMPLES")) %>% left_join(., mufa_long, by = c("Class", "SAMPLES"))
ether_fa <- left_join(pufa_long, ether_long)
# ether_info <- ether_fa %>%
#   rowwise()  %>% 
#   mutate(GROUPS = ifelse(SAMPLES %in% group_info$samples, 
#                          unlist(group_info[group_info$samples==SAMPLES, 2]), "NA")) %>% 
#   mutate("SFA_MUFA" = SFA_value + MUFA_value) %>% 
#   mutate("total_value" = PUFA_value + SFA_MUFA) %>% 
#   mutate(SAMPLES = str_remove_all(SAMPLES, "MainArea\\[") %>% str_remove_all(., "\\]")) %>% 
#   select(c("Class", "SAMPLES", "PUFA_value",  "SFA_value", "MUFA_value", "SFA_MUFA", "total_value", "GROUPS" ))
# write_csv(ether_info, "ether_pufa.csv")
ether_info <- ether_fa %>% 
  rowwise() %>%
  mutate(GROUPS = ifelse(SAMPLES %in% group_info$samples, 
                         unlist(group_info[group_info$samples==SAMPLES, 2]), "NA")) %>% 
  mutate("SFA_MUFA" = total_value - PUFA_value) %>% 
  mutate(SAMPLES = str_remove_all(SAMPLES, "MainArea\\[") %>% str_remove_all(., "\\]")) %>% 
  select(Class, SAMPLES, PUFA_value, SFA_MUFA, total_value, GROUPS)



ether1 <- ether_info %>% 
  select(c("Class", "SAMPLES",  "PUFA_value", "SFA_MUFA", "GROUPS")) %>% 
  gather(type, value, -c("Class", "SAMPLES", "GROUPS"))
name_factor <- sample_raw_list %>% str_remove_all(., "MainArea\\[") %>% str_remove_all(., "\\]")
ether1$SAMPLES <- factor(ether1$SAMPLES, levels = name_factor)


params1 <- c("SAMPLES", "value", "type")
p8 <- plot_all(data = ether1, params1) +
  geom_bar(stat = "identity") +
  scale_fill_simpsons(labels = c("PUFA", "SFA_MUFA")) +
  facet_wrap(~Class, scales = "free") +
  set_theme(theme_params = list(axis.text.x = element_text(angle = 45, size = 8, hjust = 1))) +
  add_scales() +
  scale_y_continuous(labels = scientific_format(), expand = c(0, 0, 0.2, 0)) +
  labs(x = "experiment samples", y = "value", title = "PUFA in ehter lipids for each sample", fill = "", caption = "This visualization is only for ether lipids here") 
print(p8)
ggsave("plot/Ether/pufa_ether.png", device = "png")


p9 <- ggplot(data = ether1, aes(x=SAMPLES, y = value, fill = type)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous( expand = c(0, 0, 0.1, 0), labels = scales::percent_format()) +
  scale_fill_jco(labels = c("PUFA", "SFA_MUFA")) +
  facet_wrap(~Class, scales = "free") +
  theme_bw() +
  set_theme(theme_params = list(axis.text.x = element_text(angle = 45, size = 8, hjust = 1))) +
  labs(x = "experiment samples", y = "value", title = "PUFA in ehter lipids for each sample", color = "", caption = "This visualization is only for ether lipids here") 
print(p9)
ggsave("plot/Ether/pufa_ether_percentage.png", device = "png")



aggregated_ether <- ether_saturation %>% group_by(Class) %>% summarise_at(vars(all_of(sample_raw_list)), list(~sum(.))) %>% arrange(Class)
ether_mean <- cal_lipid_statistics(aggregated_ether, group_info, "mean", "Class")
ether_mean_long <- ether_mean[[2]] 
pufa_mean <- cal_lipid_statistics(pufa_ether, group_info, "mean", "Class")
pufa_mean_long <- pufa_mean[[2]]
colnames(pufa_mean_long)[3] <- "PUFA_mean"
ether_pufa_mean <- left_join(ether_mean_long, pufa_mean_long, by = c("Class", "Groups")) %>% mutate("SFA/MUFA_mean" = mean - PUFA_mean)
ether2 <- ether_pufa_mean %>% select(-mean) %>% gather(type, value, -c("Class", "Groups"))

ether2 <- ether_info %>% 
  group_by(Class, GROUPS) %>% 
  summarise_at(vars(c(PUFA_value, SFA_MUFA)), list(~mean(.))) %>% 
  gather(type, value, -c(Class, GROUPS))
ether2$GROUPS <- factor(ether2$GROUPS, levels = group_names)

p10 <- plot_all(ether2, c("GROUPS", "value", "type")) +
  geom_bar(stat = "identity") +
  scale_fill_d3(labels = c("PUFA", "SFA_MUFA")) +
  facet_wrap(~Class, scales = "free") +
  set_theme() +
  scale_y_continuous(labels = scientific_format(), expand = c(0, 0, 0.2, 0)) +
  labs(x = "Experiment Groups", y = "value", title = "PUFA (mean) in ether lipids for each group in each lipid class", color = "", caption = "This visualization is only for ether lipids here") 
print(p10)
ggsave("plot/Ether/pufa_ether_group.png", device = "png")
ggsave("plot/Ether/pufa_ether_group.svg", device = "svg")





