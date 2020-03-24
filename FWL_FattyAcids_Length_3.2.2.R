####################################################################################
# Script: FWL_FattyAcids_Length_3.1.R
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
#         "FWL_FattyAcids_Length_3.1.R" in directory fattyAcids_saturation_analysis.
#         This script is derived from Laura's project
#####################################################################################
transformed_dt <- filtered_lipidomics %>% mutate_at(vars(sample_raw_list), list(~ifelse(.>= 0, ., NA_real_)))
split_dt <- transformed_dt %>% 
  # mutate(patterns = str_remove_all(LipidMolec, ".*\\(") %>% 
  #                                  str_remove_all(., "\\)")) %>% 
  # select(patterns, colnames(transformed_dt))
  select(FA, all_of(colnames(transformed_dt))) %>% 
#  select(contains("FA"), all_of(colnames(transformed_dt))) %>%
  rename(patterns = FA)
# filter lipid molecules which don't have chains, e.g. QA
filtered_dt <- split_dt %>% filter(str_detect(patterns, "\\:"))
if(nrow(split_dt) != nrow(filtered_dt)){
  deleted_dt <- setdiff(split_dt$LipidMolec, filtered_dt$LipidMolec)
  message("Lipid Molecules: ", paste(deleted_dt, sep = ", ", collapse = ", "), " are deleted for analysis.")
}

# get percentage for fatty acids length types for each lipid molecules
fas <- filtered_dt %>% 
   separate("patterns", c("FA1", "FA2", "FA3", "FA4"), sep = "\\/", fill = "right") %>% 
  rowwise() %>% 
  mutate("FA1%" = sum(!is.na(FA1))/sum(!is.na(c(FA1, FA2, FA3, FA4))),
         "FA2%" = sum(!is.na(FA2))/sum(!is.na(c(FA1, FA2, FA3, FA4))),
         "FA3%" = sum(!is.na(FA3))/sum(!is.na(c(FA1, FA2, FA3, FA4))),
         "FA4%" = sum(!is.na(FA4))/sum(!is.na(c(FA1, FA2, FA3, FA4)))) %>% 
   select(Rej, LipidMolec, FA1, "FA1%", FA2, "FA2%", FA3, "FA3%", FA4, "FA4%", 
          all_of(colnames(filtered_dt)[-c(1:3)]))
write_csv(fas, "data/fa_chains.csv")

# unite 2 columns into 1 column information
fas_merged <- fas %>% 
  unite(FA1, c("FA1", "FA1%"), sep = ", ") %>% 
  unite(FA2, c("FA2", "FA2%"), sep = ", ") %>% 
  unite(FA3, c("FA3", "FA3%"), sep = ", ") %>% 
  unite(FA4, c("FA4", "FA4%"), sep = ", ")

# extract AUC information
# calculate each lipid molecules fatty acids length type percentage AUC
fas_separate <- fas_merged %>% 
  select(Class, FA1, FA2, FA3, FA4, contains("MainArea")) %>% 
  gather(TYPE, FA, -c(Class, contains("MainArea"))) %>% 
  select(-TYPE, Class, FA, contains("MainArea")) %>% 
  separate(FA, c("FA", "percentage"), sep = ", ") %>% 
  mutate(percentage = as.numeric(percentage)) %>% 
  mutate_at(vars(sample_raw_list), list(~.*percentage)) %>% 
  filter(!str_detect(FA, "NA")) %>% 
  group_by(Class, FA)

write_csv(fas_separate, "data/Length/fa_chains_AUC.csv") 

# count observations of lipid molecules
fas_count <- fas_separate  %>% count(FA) 

fas_class <- fas_separate %>% summarise_at(vars(sample_raw_list), list(~sum(., na.rm = TRUE))) %>% left_join(fas_count, .)
write_csv(fas_class, "data/Length/class_chain_auc.csv")

# get pattern
fas_length <- fas_class %>% 
  mutate(length = str_remove_all(FA, "^.*\\:") %>% 
          str_remove_all(., "\\:|[A-z]") %>%
           as.numeric())
#  mutate(length = str_remove_all(FA, "\\:.*") %>% str_remove_all(., "\\(|[A-z]") %>% as.numeric())


fas_length_pattern <- fas_length %>% 
  mutate(length_type = case_when(length <= 5 ~"SCFA",
                                 length < 13~"MCFA",
                                 length < 22 ~"LCFA",
                                 length > 21 ~ "VLCFA")) %>% 
  select(Class, FA, length_type, length, all_of(sample_raw_list)) %>% 
  write_csv(., "data/Length/fas_length.csv")

fas_sum <- fas_length_pattern %>% 
  ungroup() %>% 
  select(Class, length_type, all_of(sample_raw_list)) %>% 
  mutate_if(is.character, as.factor) %>% 
  group_by(Class, length_type) %>% 
  summarise_all(list(~sum(., na.rm = TRUE))) %>% 
  ungroup() %>% 
  write_csv(., "data/Length/fas_length_sum.csv")

 fas_sum <- fas_sum %>% filter(!is.na(length_type))


length_long <- fas_sum %>% 
  gather(SAMPLES, value, -c("Class", "length_type")) %>% 
  rowwise() %>% 
  mutate(GROUPS = ifelse(SAMPLES %in% group_info$samples, 
                         unlist(group_info[group_info$samples==SAMPLES, 2]), "NA")) %>% 
  mutate(SAMPLES = str_remove_all(SAMPLES, "MainArea\\[") %>% str_remove_all(., "\\]")) %>% 
  ungroup() %>% 
  mutate_if(is.character, as.factor) 

length_type_levels <- c("SCFA", "MCFA", "LCFA", "VLCFA")
length_long$length_type <- factor(length_long$length_type, 
                                  levels = length_type_levels)
ordered_samples<- unique(length_long$SAMPLES) %>% 
  str_remove_all(., "s") %>% 
  as.numeric() %>% 
  sort() %>% 
  paste0("s", .)
length_long$SAMPLES <- factor(length_long$SAMPLES, levels = ordered_samples)
write_csv(length_long, "data/Length/fas_length_long.csv")

length_pars1 <- c("SAMPLES", "value", "length_type")
p1 <- plot_all(data = length_long, c("SAMPLES", "value", "length_type")) +
  geom_bar(stat = "identity", position = "stack") + 
  facet_wrap(~Class, scales = "free") +
  scale_y_continuous(labels = label_scientific(digits = 2), expand = c(0, 0, 0.2, 0)) +
  theme(axis.text.x  = element_text(angle = 30, hjust = 1, size = 6),
        axis.line = element_line(size = .2)) +
  labs(x = "experiment samples", y = "AUC", 
       title = "Fatty acids length of lipid classfor each sample", fill = "") +
  scale_fill_jama() 

print(p1)
ggsave("plot/Length/fa_length.png", device = "png")
ggsave("plot/Length/fa_length.svg", device = "svg")

p2 <- plot_all(data = length_long, length_pars1) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_startrek() +
  facet_wrap(~Class, scales = "free") +
  theme(axis.text.x = element_text(angle = 45, size = 8, hjust = 1),
        axis.line = element_line(size = .2)) +
  scale_y_continuous( expand = c(0, 0, 0.1, 0), labels = scales::percent_format()) +
  labs(x = "experiment samples", y = "AUC", 
       title = "Fatty acids length of lipid classfor each sample", fill = "") 
print(p2)
ggsave("plot/Length/fa_length_percentage.png", device = "png")
ggsave("plot/Length/fa_length_percentage.svg", device = "svg")


length_groups <- length_long %>% 
  group_by(Class, GROUPS, length_type) %>% 
  summarise(mean =  mean(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  write_csv(., "data/Length/length_groups.csv")
length_groups$GROUPS <- factor(length_groups$GROUPS, levels = group_names)
length_groups$length_type <- factor(length_groups$length_type, levels = length_type_levels)

length_pars2 <- c("GROUPS", "mean", "length_type")
p3 <- plot_all(data = length_groups, length_pars2) +
  geom_bar(stat = "identity") + 
  scale_y_continuous(labels = label_scientific(digits = 2), 
                     expand = c(0.01, 0, 0.2, 0)) +
  facet_wrap(~Class, scales = "free") +
  theme(axis.text.x  = element_text(angle = 30, hjust = 1, size = 6),
        axis.line = element_line(size = .2)) +
  labs(x = "Experiment groups", y = "AUC (mean value)", 
       title = "Fatty acids length of lipid class for each group", fill = "") +
  scale_fill_nejm() 
print(p3)
ggsave("plot/Length/fa_length_group.png", device = "png")
ggsave("plot/Length/fa_length_group.svg", device = "svg")

p4 <- plot_all(data = length_groups, length_pars2) +
  geom_bar(stat = "identity", position = "fill") + 
  facet_wrap(~Class, scales = "free") +
  theme(axis.text.x  = element_text(angle = 30, hjust = 1, size = 6),
        axis.line = element_line(size = .2)) +
  labs(x = "Experiment groups", y = "AUC (mean value)", 
       title = "Fatty acids length of lipid class for each group", fill = "") +
  scale_fill_jama() + 
  scale_y_continuous(expand = c(0, 0, 0.1, 0), labels = scales::percent_format())
print(p4)  
ggsave("plot/Length/fa_length_gr_percentage.png", device = "png")
ggsave("plot/Length/fa_length_gr_percentage.svg", device = "svg")


length_pars3 <- c("GROUPS", "mean", "length_type", "sd")
p5 <- plot_all(length_groups, length_pars3, se = TRUE)  +  
  geom_bar(stat = "identity",  position=position_dodge(preserve="single")) +
  facet_wrap(~Class, scales = "free") +
  scale_y_continuous(labels = label_scientific(digits = 2),
                     expand = c(0.01, 0, 0.2, 0)) +
  labs(x = "Experiment groups", y = "AUC (mean value)", 
       title = "Fatty acids length of lipid class for each group", fill = "",
       caption = "error bar is standard deviation") +
  theme(axis.text.x  = element_text(angle = 30, hjust = 1, size = 6),
        axis.line = element_line(size = .2)) + 
  scale_fill_nejm()
print(p5)  
ggsave("plot/Length/fa_length_gr.png", device = "png")
ggsave("plot/Length/fa_length_gr.svg", device = "svg")


control <- check_group(group_names, "control")
message("Please choose normalization method from mean or median for experiment group.")
method <- retype_choice("MEAN/MEDIAN")
control_names <- paste0(control, "_", method)
by_group <- c("Class", "length_type")
fa_grs <- cal_lipid_statistics(fas_sum, group_info, method, by_group)
fa_grs_control <- fa_grs[[1]] %>% select(contains(control))
fas_info <- cbind(fas_sum, fa_grs_control)

# filter lipid if control group value is 0
fas_info <- fas_info %>%  filter(eval(sym(control_names)) != 0 )
normalized_length <- norm_by_mean2(fas_info, by_group, sample_raw_list, control_names, "n")
length_info <- cal_lipid_statistics(normalized_length, group_info, method, by_group)
length_info_long <- length_info[[2]] 
if(method == "mean"){
  length_sd <- cal_lipid_statistics(normalized_length, group_info, "sd", by_group)
  length_sd_long <- length_sd[[2]] 
  length_info_long <- left_join(length_info_long, length_sd_long)
  captions <- "error bar is standard deviation"
  length_pars4 <- length_pars3
}else{
  captions = ""
  length_pars4 <- c("GROUPS", method, "length_type")
}
length_info_long <- length_info_long %>% rename(GROUPS=Groups)
write_csv(length_info_long, "data/Length/fa_length_normalized_mean.csv")
titles <- paste0("Fatty acids length of lipid class for each group (normalized by ", method, " )")
p6 <- plot_all(length_info_long, length_pars4, se = TRUE)  +  
  geom_bar(stat = "identity",  position=position_dodge(preserve="single")) +
  scale_y_continuous(expand = c(0, 0, 0.2, 0)) +
  facet_wrap(~Class, scales = "free") +
  labs(x = "Experiment groups", y = "AUC (mean value)", 
       title = titles, fill = "",
       caption = captions) +
  theme(axis.text.x  = element_text(angle = 30, hjust = 1, size = 6)) + 
  scale_fill_nejm()
print(p6)
ggsave("plot/Length/fa_length_normalized.png", device = "png")
ggsave("plot/Length/fa_length_normalized.svg", device = "svg")


# for individuale chains
message("Please choose normalization method from mean/median for chains in individule lipid molecules.") 
method <- retype_choice("MEAN/MEDIAN")
control <- check_group(group_names, "control")
control_nm <- paste0(control, "_", method)
par_group <- c("Class", "FA")
fas_ctr <- cal_lipid_statistics(fas_class, group_info, method, par_group)
fas_ctr_gr <- fas_ctr[[1]] %>% select(all_of(par_group), all_of(control_nm))
fas_all <- left_join(fas_class, fas_ctr_gr)
fas_nm <- norm_by_mean2(fas_all, par_group, sample_raw_list, control_nm, "n")


fas_class_long <- fas_nm %>% 
  ungroup() %>% 
  gather(SAMPLES, VALs, -c(Class, FA)) %>% 
  rowwise() %>% 
  mutate(GROUPS = ifelse(SAMPLES %in% group_info$samples, 
                         unlist(group_info[group_info$samples==SAMPLES, 2]), "NA")) %>% 
  mutate(SAMPLES = str_remove_all(SAMPLES, "MainArea\\[") %>% str_remove_all(., "\\]")) %>% 
  ungroup() %>% 
  filter(!is.na(VALs))
fas_class_long$SAMPLES <- factor(fas_class_long$SAMPLES, levels = ordered_samples)
fas_class_long$GROUPS <- factor(fas_class_long$GROUPS, levels = rev(group_names))
classes <- fas_class_long$Class %>% unique() 
write_csv(fas_class_long, "data/Length/normalized_chain.csv")

if(method == "mean"){
  for(i in classes){
    fa_chain <- fas_class_long %>% filter(Class == i)
    name <- paste0(i, " Relative fold change (Normaized by ", method, ")")
    file1 <- paste0("plot/Length/", i, ".fc.png")
    file2 <- paste0("plot/Length/", i, ".fc.svg")
    p <-  ggbarplot(fa_chain, x = "FA", y = "VALs",
                    fill = "GROUPS",
                    add = "mean_se", add.params = list(group = "GROUPS"),
                    position = position_dodge(0.6, preserve = "single"),error.plot = "upper_errorbar") +
      labs(title = name, y="Fold Change", x = "Fatty Acids") +
      theme(axis.text.y = element_text(angle = 30),
            legend.position = "right") +
      #scale_fill_manual(values = wes_palette("GrandBudapest2")) +
      scale_fill_simpsons() +
      # add_scales()
      scale_y_continuous(expand = c(0, 0, 0.2, 0)) +
      guides(fill = guide_legend(reverse=TRUE)) +
      coord_flip()
    print(p)
    ggsave(file1, device = "png")
    ggsave(file2, device = "svg")
  }
}else{
  for(i in classes){
    fa_chain <- fas_class_long %>% filter(Class == i)
    name <- paste0(i, " Relative fold change (Normaized by ", method, ")")
    file1 <- paste0("plot/Length/", i, ".fc.png")
    file2 <- paste0("plot/Length/", i, ".fc.svg")
    p <-  ggbarplot(fa_chain, x = "FA", y = "VALs",
                    fill = "GROUPS",
                    add = "median", add.params = list(group = "GROUPS"),
                    position = position_dodge(0.6, preserve = "single"),error.plot = "upper_errorbar") +
      labs(title = name, y="Fold Change", x = "Fatty Acids") +
      theme(axis.text.y = element_text(angle = 30),
            legend.position = "right") +
      #scale_fill_manual(values = wes_palette("GrandBudapest2")) +
      scale_fill_simpsons() +
      # add_scales()
      scale_y_continuous(expand = c(0, 0, 0.2, 0)) +
      guides(fill = guide_legend(reverse=TRUE)) +
      coord_flip()
    print(p)
    ggsave(file1, device = "png")
    ggsave(file2, device = "svg")
  }
}











