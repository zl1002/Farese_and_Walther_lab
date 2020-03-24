##########################################################################################
# Quantification of total lipid classes (mean, sd)
##########################################################################################
# making barplot graphs
message("\n Mean value for each group is used for quantification visualization.")
# sum lipid molecues by class
aggregated_class <- filtered_lipidomics %>% 
  ungroup() %>% 
  select(Class, all_of(sample_raw_list)) %>% 
  group_by(Class) %>% 
  summarise_at(sample_raw_list, list(~sum(.)))
write_csv(aggregated_class, "data/Quantification/aggregated_class.csv")

# get the mean and sd of each group for each class
total_group_mean <- cal_lipid_statistics(aggregated_class, group_info, "mean", "Class")
#class_mean_wide <- total_group_mean[[1]]
class_mean_long <- total_group_mean[[2]]
total_group_sd <- cal_lipid_statistics(aggregated_class, group_info, "sd", "Class")
class_sd_long <- total_group_sd[[2]]
# bind mean and sd for group data
total_data <- left_join(class_mean_long, class_sd_long, by = c("Class", "Groups"))
write_csv(total_data, "data/Quantification/total_class.csv")

#paras1 <- list("Class", "mean", "Groups", "sd", scale.params)
paras1 <- c("Class", "mean", "Groups", "sd")
total_data$Groups <- factor(total_data$Groups, levels = rev(group_names))
total_plot <- plot_all(total_data, paras1, se=TRUE) + 
  geom_bar(position=position_dodge(preserve="single"), stat="identity") + 
  labs(title = "Total lipid classes", x = "Lipid Classes", y= "AUC (Area under curve)",
       fill = "Experiment Groups") +
  coord_flip() + 
  guides(fill = guide_legend(reverse = TRUE))  

  #scale_y_continuous(labels = scientific_format(), expand = c(0, 0, 0.1,0))

print(total_plot)
# save filtered total class lipid abundance in plot total.class.png
ggsave(filename="total.class.png", path = 'plot/Quantification/', device = "png")


##########################################################################################
# Quantification of individual lipid class
##########################################################################################
# get mean data of each group for each lipid molecule
# you can also use median
group_par <- c("LipidMolec", "Class")
each_group_mean <- cal_lipid_statistics(filtered_lipidomics, group_info, "mean", group_par)
lipid_mean_wide <- each_group_mean[[1]]
lipid_mean_long <- each_group_mean[[2]]
each_group_sd <- cal_lipid_statistics(filtered_lipidomics, group_info,  "sd", group_par)
lipid_sd_wide <- each_group_sd[[1]]
lipid_sd_long <- each_group_sd[[2]]
# merge mean, sd value and reformat data
each_class <- left_join(lipid_mean_long, lipid_sd_long) %>% ungroup()
write_csv(each_class, "data/Quantification/all_lipidmolec.csv")
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
EachClassPlot(each_class, pars, robot)

# # overview of each class plot 
# nbar <- lipidmolecNO_max
# post_name <- "all"
# pars <- list(nbar, ngroups, par_eachclass, plot_all, post_name, labs1)
# message("\nAlternative display quantification of individule lipid class (all lipids in a class in the same png)")
# EachClassPlot(each_class, pars)
# 

##########################################################################################
# dev.off()
# options(device = "RStudioGD")
##########################################################################################


##########################################################################################
# visualization of normalized lipid class and individule molecules

##########################################################################################
# visualization of lipid class data, normalized by median value
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

# remove lipid classes which control group value is 0 or negative
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
write_csv(dt, "data/Quantification/raw_class_median.csv")
write_csv(class_nm, "data/Quantification/normalized_class_median.csv")

# reformat the data
class_long <- class_nm %>%  gather("SAMPLES", "VALs", -Class) %>% ungroup()
class_long <- class_long %>% 
  rowwise() %>%  
  mutate(GROUPS = ifelse(SAMPLES %in% group_info$samples, 
                         unlist(group_info[group_info$samples==SAMPLES, 2]), "NA")) %>% 
  mutate(SAMPLES = str_remove_all(SAMPLES, "MainArea\\[") %>% str_remove_all(., "\\]")) %>% 
  ungroup()
  
class_long$GROUPS <- factor(class_long$GROUPS, levels = group_names) 
write_csv(class_long, "data/normalized_class_median_long.csv")
# plot relative fold change of lipid class (normalized by median), dotplot
paras2 <- c("GROUPS", "VALs", "GROUPS")
pd <- plot_all(class_long, paras2) +
  geom_dotplot(binaxis = "y",         # which axis to bin along
               stackdir='center', dotsize = 1, position = "dodge") +
  stat_summary(fun.y = median, geom = "point",  size = 5, color = "red", shape = 95, alpha = 0.8, stroke = 3) +
  facet_wrap(~Class, scales = "free") +
  #scale_fill_manual(values = wes_palette("GrandBudapest2")) +
  scale_fill_startrek() +
  theme(axis.text.x =  element_text(angle=45, hjust=1, size =6),
        axis.line = element_line(size = 0.2)) +
  scale_y_continuous(expand = c(0, 0, 0.1, 0)) +
  expand_limits(y = 0) +
  geom_text_repel(aes(label = SAMPLES), size = 2) +
  #scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  labs(title = "Normalized by median", caption = "each dot represent the replicate of one lipid class", 
       x = "Experiment Groups", y = "Relative Fold Change", fill="Experiment Groups")
print(pd)
ggsave(filename = "class_median_dot.png", path="plot/Quantification", device = "png")
ggsave(filename = "class_median_dot.svg", path="plot/Quantification", device = "svg")

# plot relative fold change of lipid class (normalized by median), boxplot
pb <- plot_all(class_long, paras2) +
  geom_boxplot(outlier.shape = NA, fatten = 1) +
  facet_wrap(~Class, scales = "free") +
  labs(title = "Lipid class data (normalized by median)", x = "Lipid Class", 
       y = "Relative Foldchange", fill = "Experiment Groups") +
  scale_y_continuous(expand = c(0, 0, 0.1, 0)) +
  expand_limits(y = 0) +
  theme(axis.text.x =  element_text(angle=45, hjust=1, size =6),
        axis.line = element_line(size = 0.2)) +
  #scale_fill_d3()
  scale_fill_manual(values = wes_palette("GrandBudapest2"))
print(pb)
ggsave(filename = "class_median_box.png", path="plot/Quantification", device = "png")



##########################################################################################
# visualization of lipid molecule data, normalized by mean/median value
##########################################################################################
message("\nFold change analysis for lipid molecules.")
# mutate negative value into NA 
filtered_dt <- filtered_lipidomics %>% mutate_at(sample_raw_list, list(~ifelse(.>= 0, ., NA)))
group_bar <- c("LipidMolec", "Class")
message("\nPlease choose a normalization method from Mean or Median.")
method <- retype_choice("mean/median")
molec_group_mean <- filtered_dt  %>%  cal_lipid_statistics(., group_info, "mean", group_bar)
molec_group_median <- filtered_dt  %>%  cal_lipid_statistics(., group_info, "median", group_bar)
selected_lipid <- filtered_dt %>% ungroup() %>% select(all_of(group_bar), all_of(sample_raw_list)) 
molec_all <- left_join(selected_lipid, molec_group_mean[[1]]) %>% 
  left_join(., molec_group_median[[1]])
write_csv(molec_all, "data/Quantification/molecules_group_statistics.csv")
# customize times of comparison  
message("\nHow many times of fold change analysis for groups you want to do?")
fc_times <- readline("Please input the number: ")
for(i in 1:fc_times){
  message("\nPlease type the name of the control group.")
  control <- check_group(group_names, "control")
  control_name <- paste0(control, "_", method)
  message("\nPlease type the name of the contrast group(s).")
  contrasts <- check_group(group_names, "contrast")
  selected_group <- molec_all %>% select(LipidMolec, Class, all_of(sample_raw_list), !!control_name)
  # filter 0 or negative in control group
  filtered_molec <- selected_group %>% filter_at(vars(contains(control)), any_vars(.>0))
  write_csv(filtered_molec, "data/Quantification/raw_molec.csv")
  # normaliz lipid molecule
  molec_nm <- norm_by_mean2(filtered_molec, group_bar, sample_raw_list, control_name, "n") 
  write_csv(molec_nm, "data/Quantification/normalized_molec.csv")
  molec_info <- cal_lipid_statistics(molec_nm, group_info, method, group_bar)
  molec_info_long <- molec_info[[2]] 
  gtitle <- paste0("Normalized by (", method, ")")
  if(method == "mean"){
    molec_sd <- cal_lipid_statistics(molec_nm, group_info, "sd", group_bar)
    molec_sd_long <- molec_sd[[2]] 
    molec_info_long <- left_join(molec_info_long, molec_sd_long)
    captions <- "error bar is standard deviation"
    par_eachclass <- c("LipidMolec", "mean", "Groups", "sd")  
    labs<- labs(x="Acyl composition", y="Relative fold change", title = gtitle,
                caption="Error bar is the standard deviation for each class in each group", fill = NULL)
  }else{
    captions = ""
    par_eachclass <- c("LipidMolec", "median", "Groups") 
    labs<- labs(x="Acyl composition", y="Main Area", title = gtitle, fill = NULL)
  }
  file_name <- paste0("data/Quantification/normalized_molec_", method, ".csv")
  write_csv(molec_info_long, file_name)
  selected_molec <- molec_info_long %>% 
    filter(Groups %in% c(control, contrasts)) %>% 
    ungroup()
  post_name <- paste0("fc/", control, "_", contrasts, "_", method, "_", i, "_")
  #par_eachclass <- list("LipidMolec", "mean", "Groups", "sd", scales)  
  n <- length(c(control, contrasts))
 # nbar <- selected_molec %>% group_by(Class) %>% tally() %>% select(n) %>% unlist() %>% max()
  nbar <- 80
  pars <- list(nbar, n, par_eachclass, plot_fc, post_name, labs)
  EachClassPlot(selected_molec, pars, robot)
 }
  
  

##########################################################################################
# dev.off()
# options(device = "RStudioGD")
##########################################################################################

# violin plot
message("\nPlease choose a value (mean/median) for visualize lipid molecules in violin plot.
        \n Please choose a normalization method from MEAN or MEDIAN.")
method <- retype_choice("mean/median")
control <- check_group(group_names, "control")
control_gr <- paste0(control, "_", method)
selected_molec_all <- molec_all %>% select(LipidMolec, Class, contains("MainArea"),!!control_gr)
# filter lipid molecules which control group value is 0 or negative
filtered_molec_all <- selected_molec_all %>% filter(eval(!!control_gr)>0)
molec_all_nm <- norm_by_mean2(filtered_molec_all, c("LipidMolec", "Class"), sample_raw_list, control_gr, "n") 
molecules <- molec_all_nm  %>%  gather("SAMPLES", "VALs", -c("LipidMolec", "Class"))
molecules <- molecules %>% 
  ungroup() %>%  
  rowwise() %>%  
  mutate(GROUPS = ifelse(SAMPLES %in% group_info$samples, 
                         unlist(group_info[group_info$samples==SAMPLES, 2]), "NA")) %>% 
  mutate(SAMPLES = str_remove_all(SAMPLES, "MainArea\\[") %>% str_remove_all(., "\\]")) %>% 
  ungroup()
# set information for tooltip
molecules <- molecules %>% 
  mutate(Molec = str_remove_all(LipidMolec, ".*\\(") %>% str_remove_all(., "\\)")) %>% 
  mutate(Molecule = paste0(Molec , "\nVALs: ", formatC(VALs,  digits = 2), "\n"))

molecules$GROUPS <- factor(molecules$GROUPS, levels = group_names)
name <- paste0("data/Quantification/molecules_", method, "_violin.csv")
write_csv(molecules, name)

# violin plot
pv <- ggplot(molecules, aes(GROUPS, VALs, label=Molecule)) +
  geom_violin(aes(fill = GROUPS), trim = TRUE, na.rm = FALSE)+
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
  theme_bw() +
  set_theme() +
  theme(axis.text.x =  element_text(angle=45, hjust=1, size = 6),
        axis.line = element_line(size = 0.3, colour = "black")) +
  expand_limits(y = 0) +
  scale_y_continuous(expand = c(0, 0, 0.2, 0)) +
  labs(title = paste0("Individual Lipidmolecules (normalized by ", method, ")"),
       caption = paste0("each dot represent a replicate of one lipid molecule\nThe red line represent ", method),
      y = "Relative Foldchange") +
  scale_fill_manual(values = wes_palette("GrandBudapest2")) 
  
print(pv)
ggsave(filename = "molec_violin.png", path="plot/Quantification", device = "png")
ggsave(filename = "molec_violin.svg", path="plot/Quantification", device = "svg")

# set path
dir_path <- paste0(getwd(),"/plot/Quantification") 
# save interavtive verstion of violin plat
pv %>% ggplotly(tooltip="Molecule") %>% 
  layout(dragmode = "select") %>% 
  saveWidget(., file.path(dir_path, "molec_violin_all.html"))

# customized interactive violin plot
lipid <- molecules$Class %>% 
  unique()
classes <-  lipid %>% 
  paste0(., sep=", ", collapse = "") %>% 
  substr(., 1, nchar(.)-2) 
molecules$VALs <- round(molecules$VALs, 2)

message("To save running time, the violin plot will just display the lipid class that you choose.
        \nPlease choose a lipid class name from the list below: \n", classes, ".")

choice <- readline("Please type the name of lipid class: ") %>% 
  str_split(., "\\s+") %>% 
  unlist()

if(!choice %in% lipid){
  message("Please type exact lipid class name from the list: ", classes, ".")
  choice <- readline("Please type the name of lipid class: ") %>% 
    str_split(., "\\s+") %>% 
    unlist()
}

for(i in choice){
  dt_molec <- molecules %>% filter(Class == i)
  heads <- paste0("Lipid Molecules distribution of Lipid Class ", i, "\n(normalized by ", method, ")")
  number <- dt_molec$GROUPS %>% unique() %>% length()
  position <- 1:number
  annotates_list <- lapply(group_names, function(x) list(showarrow = FALSE, x = 1, y = -0.5, text = x))
  annotates <- Map(modifyList, annotates_list, val = lapply(position, function(x) list("x" = x)))
  
  p <- dt_molec %>%
    plot_ly() %>% 
    add_trace(x = ~as.numeric(GROUPS), y = ~VALs, color = ~GROUPS, 
              type = "violin", box = list(visible = F), 
              hoverinfo = 'name+y' 
              ) %>%
    add_markers(x = ~jitter(as.numeric(GROUPS)), y = ~VALs, color = ~GROUPS,
                marker = list(size = 6),
                hoverinfo = "text",
                text = ~paste0("LipidMolec: ", LipidMolec,
                               "<br>VALs: ", VALs,
                               "<br>Sample: ", SAMPLES),
                showlegend = FALSE) %>%
    layout(title = heads,
           legend = list(orientation = "h",
                         x = 0.1, xanchor = "center",
                         y = 1, yanchor = "bottom"
    ),
    xaxis = list(title = "Experiment Groups",
                 tickvals = position,
                 showticklabels = FALSE,
                 zeroline = TRUE,
                 zerolinecolor = toRGB("black"),
                 zerolinewidth = 5,
                 showline = FALSE),
    yaxis = list(title = "Relative Foldchange",
                 showline = TRUE,
                 zeroline = TRUE,
                 zerolinecolor = toRGB("black"),
                 zerolinewidth = 4),
    annotations = annotates)
  print(p)
  name <- paste0(i, "_violin.html")
  saveWidget(as_widget(p), file.path(dir_path, name))
  
}











