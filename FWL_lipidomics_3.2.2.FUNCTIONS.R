
####################################################################################
# Script: FWL_lipidomics_3.2.0.FUNCTIONS.R
# Author: Wenting, Niklas, Kenny
# Notes:  This script assist executing for main script FWL_lipidomics_3.2.0.R which
#         helps generating the graph and data for the workflow of lipidomics.
#         To start, typing command in the console-----> source("FWL_lipidomics_3.2.0.R")
#         or press the source button.
#         This pipeline is based on Niklas's and Kenny's previous work. Their original 
#         source scripts are in the quality_control and statistics_quantification directories
#
# Warning: For mac users, please make sure XQuartz is installed. Otherwise the pipeline 
#           might crash when generating summary plots for each class
#####################################################################################

##### Please install the lacking packages ------> install.packages("package")
# source("https://install-github.me/tidyverse/rlang")
initial.package <- c("BiocManager", "rlang", "devtools", "svglite")
need.package <- initial.package[!(initial.package %in% installed.packages()[, "Package"])]   
if(length(need.package) > 0) install.packages(need.package)


list.of.packages <- c( "FactoMineR", "factoextra", "scales", "magrittr", "ggrepel", "reshape2",
                       "stargazer", "Hmisc", "limma", "factoextra", "scales", "RColorBrewer",
                       "stringr", "readxl", "RCy3", "igraph", "tidyverse", "dplyr","viridis", 
                       "igraph", "network", "visNetwork", "extrafont", "ggforce", "kableExtra", 
                       "ggpubr", "wesanderson", "formattable", "ggsci", "plotly", "htmlwidgets", 
                       "gridExtra", "grid", "DT")
need.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if(length(need.packages) > 0) BiocManager::install(need.packages)


# Suppress the package messages
lapply(list.of.packages, function(x) 
  suppressMessages(require(x, character.only = TRUE, quietly = TRUE)))




# color which are kind for color blinded people, source from Internet
clPalette1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
clPalette2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
clPalette3 <- rep("#bababa", 30)


##########################################################################################
### Functions
#########################################################################################
# function name: set_theme
# utility: setting customized theme which also can be modified
##########################################################################################
set_theme <- function(..., theme_params = list()){
  params <- list(...)
  theme_params <- modifyList(params, theme_params)
  themes <- do.call("theme", modifyList(
    list(panel.grid  = element_blank(),
         panel.border = element_blank(),
         axis.line = element_line(colour = "black", 
                                  size = 1,
                                  #arrow = arrow(length = unit(0.06, "cm"), type = "closed")
                                  arrow = NULL
         ),
         axis.ticks.length = unit(1.5, "mm"),
         axis.ticks = element_line(size = 1.2),
         axis.text = element_text(size = 10,  colour = "black"),
         axis.title =  element_text(size = 14),
         # axis.text = element_text(size = 10),
         #axis.ticks = element_line(size = .2),
         #axis.ticks.length = unit(.5, "mm"),
         strip.text = element_text(size = 10),
         title = element_text(size = 16),
         strip.background = element_blank()),
    theme_params
  ))
  themes
}



##########################################################################################
### Functions
#########################################################################################
# function name: tn
# utility: setting customized log transformed scales. Please note that it is not friendly 
# for fatty acids saturation/length stack plots display.
##########################################################################################
tn <- trans_new("custom_log_scale",
                transform = function(x)ifelse(x!= 0, sign(x)*log10(abs(x)), 0),
                inverse = function(x) ifelse(x != 0, sign(x)*10^abs(x), 0),
                domain = c(-Inf, Inf)
)



##########################################################################################
### Functions
#########################################################################################
# function name: mkdirs
# utility: Check directory existence, if not then make one
##########################################################################################
mkdirs <- function(x){
  for(i in 1:length(x))
    if(!file.exists(x[i])){
      mkdirs(dirname(x[i]))
      dir.create(x[i])
    } 
}




#########################################################################################
# function name: addquotes
# utility: add quotes for strings. It could be replaced with paste0 function
# source: Advanced R
#########################################################################################
addquotes <- function(...){
  args <- ensyms(...)
  paste(purrr::map(args, rlang::as_string), collapse = "")
}


#########################################################################################
# function name: add_all
# utility: add quotes for each string separately
#         eg. add_all(Apple, Grapes) will generate: "Apple", "Grapes"
#########################################################################################
add_all <- function(...){
  args <- ensyms(...)
  paste(purrr::map(args, rlang::as_string), sep = ",")
}


#########################################################################################
# function name: addlists
# parameters: x, y, e.g. f(c("x", "y", "z"), "_style")
# utility: add y symbol as post fix for each item in list of x
#         eg. addlists(c("x", "y", "z"), "_style") will generate: "
#         x_style", "y_style", "z_style"
#########################################################################################
addlists <- function(x, y){
  quotes <- sapply(x, function(x) addquotes(!!x, !!y))
  return(quotes)
}


#########################################################################################
# function name: reverse_addlists
# parameters: x, y, e.g. f(c("x", "y", "z"), "pre")
# utility: add y symbol as post fix for each item in list of x
#         eg. addlists("pre_", c("x", "y", "z")) will generate: "pre_x", "", "z_style"
#########################################################################################
reverse_addlists <- function(x, y){
  quotes <- sapply(x, function(x) addquotes(!!y, !!x))
  return(quotes)
}





retype_choice <- function(choices){
  choice <- choices %>% str_to_lower() %>% str_split(., "/") %>% unlist
  command <- paste0("Please type ", choices, ": ")
  option <-  readline(prompt = command) %>% str_to_lower()
  if(option %in% choice){
    return(option)
  }else{
    message("Type wrong.")
    option <- retype_choice(choices)
  }
}

#########################################################################################
# function name: select_expr_out
# utility: build expression for list of variables which can be used by data 
# frame to delete the corresponding variables
#########################################################################################
select_expr_out <- function(options){
  option_list <- sapply(options, function(x) expr(-contains(!!!x)))
  return(option_list)
}






#########################################################################################
# function name: delete_samples
# parameter: options (y/n)
# utility: This function will check if there are option control samples exist. 
#          It will delete the option controls for QC and later PCA and correlations analysis. 
#         It will also make new columns which store the information of sums of Grade A, B, 
#         C and D, plus the p values which are less or equal than 0.001 for each row. It will return the 
#         the extra information
#########################################################################################
delete_samples <- function(options, lipid){
  if(options=="y"){
    message("\nIndicate which samples used as controls or will be deleted for analysis\nSample ID , eg. s22 s23")
    # message("Sample ID for extracted and unextracted option control, eg. s22 s23")
    option_standards <- as.character(readline("option standards -----> ")) %>% 
      str_to_lower(.) %>% 
      str_split(., "\\s") %>% 
      unlist(.)
    # build formats for samples, e.g. option_Grades: MainGrade[s3], MainGrade[s4]
    option_Grades <- paste0("MainGrade[", option_standards)
    option_P.Nos <- paste0("APValue[", option_standards)
    option_MainArea <- paste0("MainArea[", option_standards)
    # build a list contains the expression variables for removal
    list <- select_expr_out(c(option_P.Nos, option_MainArea, option_Grades))
    # updated lipid by subtracting the options
    for(i in list){
      lipid <- lipid %>% select(!!i)
    }
    selected_lipidomics <- lipid
  }else if(options=="n"){
    selected_lipidomics <- lipid %>%
      select(contains("Grade"), contains("APValue"))
  }else{
    print("You are typing wrong information, please type again") 
    internal <- readline("Please type Y/N: ") %>% str_to_lower(.)
    return(delete_samples(internal))
  }
  lipid_count1 <- selected_lipidomics %>% 
    select(contains("MainGrade")) %>%
    transmute(A = rowSums(.=="A", na.rm = TRUE), B = rowSums(. == "B", na.rm = TRUE),
              C = rowSums(. == "C", na.rm = TRUE), D = rowSums(. == "D", na.rm = TRUE), No.grade = rowSums(.=="-", na.rm = TRUE))
  lipid_count <- selected_lipidomics %>% 
    select(contains("APValue")) %>% 
    transmute(APValue.001 = rowSums(.<=0.001)) %>% 
    bind_cols(., lipid_count1)
  #selected_lipidomics <- bind_cols(selected_lipidomics, lipid_count)
  return(list(lipid_count, lipid))
}



# #########################################################################################
# function name: geom_se
# parameters: ..., se.params
# utility: add error bar for plots. Parameters ... are 2 lists. Its first list is made of 4 
#         parameters while the last one element is standard deviation.
#         The second list is logical element, e.g. se=FALSE or se=TRUE.
#         parameter se.params are the list user can pass into function for modification.
# #########################################################################################
# geom_se <- function(...){
#   parameters <- syms(...)
#   p2 <- parameters[[2]]
#   p4 <- parameters[[4]]
#   se_condition <- paste0(syms(...)[[5]]) %>% as.logical()
#   if(se_condition){
#     se <- geom_errorbar(aes(ymin=sign(eval(p2))*abs(eval(p2)), 
#                       ymax=sign(eval(p2))* (abs(eval(p2)) + abs(eval(p4)))),
#                   position=position_dodge(0.9),
#                   width = 0.1,
#                   size=.2)
#     return(se)
#   }
# }
geom_se <- function(parameters, se){
  parameters <- syms(parameters)
  p2 <- parameters[[2]]
  p4 <- parameters[[4]]
  list(
    if(se)
      geom_errorbar(aes(ymin=sign(eval(p2))*abs(eval(p2)), 
                        ymax=sign(eval(p2))* (abs(eval(p2)) + eval(p4))),
                    position=position_dodge(preserve = "single", .9),
                    width = 0.1,
                    size=.2) 
  )
}




se_position <- function(..., position.params = list()){
  params <- list(...)
  position.params <- modifyList(params, bar.params)
  se.position <- do.call("geom_errobar", modifyList(
    list(position = position_dodge(width = 0.6), stat="identity", width = 0.4),
    position.params
  ))
  se.position
}



plot_bars <- function(..., bar.params = list()){
  params <- list(...)
  bar.params <- modifyList(params, bar.params)
  bar <- do.call("geom_bar", modifyList(
    list(position = position_dodge(width = 0.6), stat="identity", width = 0.4),
    bar.params
  ))
  bar
}



# add_scales <- function(..., scale.params = list()){
#   params <- list(...)
#   scale.params <- modifyList(params, scale.params)
#   scale_fun <- do.call("scale_y_continuous", modifyList(
#     list(trans = tn,
#          breaks = c( -10^3, 0, 10^3, 10^6, 10^9, 10^12, 10^15),
#          labels = c( bquote(-10^3), 0, bquote(10^3), bquote(10^6), bquote(10^9), bquote(10^12), bquote(10^15)),
#          expand = c(0, 0, 0.1, 0)),
#     scale.params
#   ))
#   list(scale_fun)
# }



# scale_params <- list(trans = tn,
#                      breaks = c(-10^3, 0, 10^3, 10^6, 10^9, 10^12, 10^15, 10^18),
#                      # labels = c(bquote(-10^3), 0, bquote(10^3), bquote(10^6),
#                      #            bquote(10^9), bquote(10^12), bquote(10^15), bquote(10^18)),
#                      labels = c(expression(paste("10"^"-3")), 0, expression(paste("10"^"3")),
#                                 expression(paste("10"^"6")),expression(paste("10"^"9")),
#                                 expression(paste("10"^"-3")),expression(paste("10"^"-3")),
#                                 expression(paste("10"^"-3"))),
#                      expand = c(0, 0, 0.1, 0))
# 
# 
# 
# 



#########################################################################################
# function name: add_scales
# parameters: ..., scale.params 
# utility: add customized log transformed scale layer for axis display. The parameters trans,
#         breaks, labels and expand could be modified.
##########################################################################################
add_scales <- function(..., scale.params = list()){
  params <- list(...)
  scale.params <- modifyList(params, scale.params)
  scale_fun <- do.call("scale_y_continuous", modifyList(
    # list(trans = "identity",
    #      breaks = waiver(),
    #      labels=scales::percent,
    #      expand = c(0, 0, 0.2, 0)),
    list(trans = tn,
         breaks = c(-10^3, 0, 10^3, 10^6, 10^9, 10^12, 10^15, 10^18),
         labels = c(bquote(-10^3), 0, bquote(10^3), bquote(10^6),
                    bquote(10^9), bquote(10^12), bquote(10^15), bquote(10^18)),
         expand = c(0, 0, 0.2, 0)),
    # scale_params,
    scale.params
  ))
  scale_fun
}






########################################################################
# function name: plot_all
# parameters: data, parameters (list of paramters), ...
# utility: passing data and a list of parameters for plotting
########################################################################
plot_all <- function(data, parameters, se){
  #paras <- syms(parameters)
  if(length(parameters) == 1){
    paras <- syms(parameters)
    p <- ggplot(data, aes(eval(paras[[1]]))) +
      theme_bw() +
      set_theme() 
    return(p)
  } else if(length(parameters) == 2){
    paras <- syms(parameters)
    p <- ggplot(data, aes(x = reorder(eval(paras[[1]]), eval(paras[[2]])), y = eval(paras[[2]]))) +
      # p <- ggplot(data, aes(x = eval(paras[[1]]), y = eval(paras[[2]]))) +
      theme_bw() + 
      set_theme() 
    # add_scales()
    return(p) 
  } else if(length(parameters) == 3){
    paras <- syms(parameters)
    p <- ggplot(data, aes(x = eval(paras[[1]]), y = eval(paras[[2]]), fill = eval(paras[[3]]))) +
      theme_bw() +
      set_theme() 
    axis_st <- data %>% filter_at(vars(!!paras[[2]]), any_vars(.<0)) %>% nrow()
    if(axis_st >0){
      p <- p + 
        add_scales(scale.params = list(expand = c(0.02, 0, 0.2, 0))) +
        geom_hline(yintercept = 0, color = "black",size = 1, linetype = "dashed") +
        scale_fill_d3()
    }else{
      p <-  p + add_scales()
    }
    return(p)
  } else{
    
    paras <- syms(parameters[1:4])
    # scale_info <- parameters[[5]]
    ranges <- data %>% 
      mutate(bounds=sign(eval(paras[[2]])) * (abs(eval(paras[[2]])) + abs(eval(paras[[4]])))) %>% 
      summarise(max(bounds), min(bounds))
    limits <- sapply(ranges, function(x)sign(x)*10^(ceil(log10(abs(x)), 0)))
    # limits[[2]] <- limits[[2]] - 10^5
    axis_st <- data %>% filter_at(vars(!!paras[[2]]), any_vars(.<0)) %>% nrow()
    
    p <- ggplot(data, aes(x = eval(paras[[1]]),
                          y = eval(paras[[2]]), fill = eval(paras[[3]]))) +
      #   expand_limits(y = limits) +
      theme_bw() +
      set_theme() +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.y = element_text(angle = 30, hjust=1)) + 
      # scale_y_continuous(trans = tn,
      #                    breaks = c(-10^3, 0, 10^3, 10^6, 10^9, 10^12, 10^15),
      #                    labels = c(bquote(-10^3), 0, bquote(10^3), bquote(10^6), bquote(10^9), bquote(10^12), bquote(10^15)),
      #                    expand = c(0, 0, 0.1, 0)) +
      #add_scales(scale.params = scale_info)+ 
      #add_scales() +
      # ylim(limits) +
      scale_fill_manual(values = clPalette1) +
      geom_se(parameters[1:4], se) 
    
    if(axis_st >0){
      p <- p + 
        add_scales(scale.params = list(expand = c(0.02, 0, 0.2, 0))) +
        geom_hline(yintercept = 0, color = "black",size = 1, linetype = "dashed") +
        scale_fill_d3()
    }else{
      p <-  p + add_scales()
    }
    return(p)
  }
}




# plot_all <- function(data, parameters, se, order_ax){
#   #paras <- syms(parameters)
#   if(length(parameters) == 1){
#     paras <- syms(parameters)
#     
#     p <- ggplot(data, aes(eval(paras[[1]]))) +
#       theme_bw() +
#       set_theme() 
#     return(p)
#   } else if(length(parameters) == 2){
#     paras <- syms(parameters)
#     
#     print(head(data))
#     p <- ggplot(data, aes(x = reorder(eval(paras[[1]]), eval(paras[[2]])), y = eval(paras[[2]]))) +
#       theme_bw() + 
#       set_theme() 
#     # add_scales()
#     return(p) 
#   } else if(length(parameters) == 3){
#     paras <- syms(parameters)
#     
#     p <- ggplot(data, aes(x = eval(paras[[1]]), y = eval(paras[[2]]), fill = eval(paras[[3]]))) +
#       theme_bw() +
#       set_theme() 
#     return(p)
#   } else{
#     
#     paras <- syms(parameters[1:4])
#     # scale_info <- parameters[[5]]
#     ranges <- data %>% 
#       mutate(bounds=sign(eval(paras[[2]])) * (abs(eval(paras[[2]])) + abs(eval(paras[[4]])))) %>% 
#       summarise(max(bounds), min(bounds))
#     limits <- sapply(ranges, function(x)sign(x)*10^(ceil(log10(abs(x)), 0)))
#     # limits[[2]] <- limits[[2]] - 10^5
#     axis_st <- data %>% filter_at(vars(!!paras[[2]]), any_vars(.<0)) %>% nrow()
#     
#     if(order_ax){
#       p <- ggplot(data, aes(x = reorder(eval(paras[[1]]), eval(paras[[2]])), 
#                             y = eval(paras[[2]]), fill = eval(paras[[3]])))
#     }else{
#       p <- ggplot(data, aes(x = eval(paras[[1]]), 
#                        y = eval(paras[[2]]), fill = eval(paras[[3]])))
#     }
#     p <- p +
#       #   expand_limits(y = limits) +
#       theme_bw() +
#       set_theme() +
#       theme(plot.title = element_text(hjust = 0.5),
#             axis.text.y = element_text(angle = 30, hjust=1)) + 
#       # scale_y_continuous(trans = tn,
#       #                    breaks = c(-10^3, 0, 10^3, 10^6, 10^9, 10^12, 10^15),
#       #                    labels = c(bquote(-10^3), 0, bquote(10^3), bquote(10^6), bquote(10^9), bquote(10^12), bquote(10^15)),
#       #                    expand = c(0, 0, 0.1, 0)) +
#       #add_scales(scale.params = scale_info)+ 
#       #add_scales() +
#       # ylim(limits) +
#       scale_fill_manual(values = clPalette1) +
#       geom_se(parameters[1:4], se)  
#     
#     if(axis_st >0){
#       p <- p + 
#         add_scales(scale.params = list(expand = c(0.02, 0, 0.2, 0))) +
#         geom_hline(yintercept = 0, color = "black",size = 1, linetype = "dashed") +
#         scale_fill_d3()
#     }else{
#       p <-  p + add_scales()
#     }
#     return(p)
#   }
# }
# 
















subtract_background <- function(data, sample_list, subtraction){
  # total class lipids of raw data by subtracting background (blank sample c) from all samples
  data <- data %>% 
    mutate_at(vars(all_of(sample_list)), list(~.-eval(sym(!!subtraction)))) %>%  
    mutate_at(vars(!!subtraction), list(~ifelse(.==0, ., 0)))
  return(data)
}




########################################################################
# function name: plot_rt_allsamples
# parameter: data (retention_data), parameters (MainArea, BaseRt)
# utility: plot the retention time for all sample abundances (log transformed) 
#         and return the plot object
########################################################################
plot_rt_allsamples <- function(data, parameters){
  paras <- syms(parameters[1:2])
  n <- parameters[3] %>% as.numeric()
  retention.all.plot <- ggplot(data=data, aes(x=log10(eval(paras[[1]])), y = eval(paras[[2]]))) +
    geom_point() +
    theme_bw() +
    # facet_grid(.~Class) +
    labs("Abundance VS. Retention time",
         x="All samples (log10(AUC))", y="Retention time (mins)") +
    theme(panel.border = element_rect(size = .2, colour = "black"),
          strip.text = element_text(size = 6),
          axis.text = element_text(size = 6),
          axis.ticks = element_line(size = .2),
          axis.ticks.length = unit(.5, "mm"),
          strip.background = element_rect(size = .2, colour = "black", fill = "white"),
          axis.title = element_text(size = 8)) 
  if(n<11){
    retention.all.plot <- retention.all.plot +
      facet_grid(.~Class)
  }else{
    retention.all.plot <- retention.all.plot +
      facet_wrap(.~Class, nrow = 2)
  }
  return(retention.all.plot)
  #ggsave(filename = "all.retention.pdf", path = 'plot/', units = "in", height=4, width=5, dpi=300, dev="pdf")
}


########################################################################
# function name: IndividuleRetentionPlot
# parameter: lipid_data, e.g. filtered_lipidomics; sample, e.g. s1
# utility: This function is NOT used. If you want to use it, please uncomment the 
#          correponding part in the main script. it will 
#           plot retention time for individule sample and save the plot 
########################################################################
IndividuleRetentionPlot <- function(lipid_data, sample){
  retention.plot <-   ggplot(data=lipid_data,
                             aes_string(x = sprintf("log10(%s)", sample), y = "BaseRt")) +
    geom_point() +
    theme_bw() +
    facet_grid(.~Class) +
    labs("Abundance VS. Retention time",
         x="lipid class (log10(MainArea))", y="Retention time")
  print(retention.plot)
  message("Please input the name you want to store for the graph. e.g. retention.pdf")
  plot_name <- readline("QC plot name: ")
  ggsave(filename = plot_name, path = 'plot/', device = "pdf", dpi=300)
}


########################################################################
# function name: detect_duplicates
# parameter: lipid_data
# utility: This function detect same lipid molecules with different retention time, 
#           return the duplicated lipid molecules
########################################################################
detect_duplicates <- function(lipid_data){
  # get same molecule which contains 0 value and filter the lipids by standard
  options(scipen = 999)
  # group lipid molecs by class names
  class_lipids <- lipid_data %>% 
    arrange(Class, LipidMolec) %>% 
    select(Class, LipidMolec, BaseRt, MainIon, contains("MainArea")) 
  # gain size information for each class
  class_size <-   class_lipids%>% 
    group_by(Class, LipidMolec) %>% 
    group_size()
  # get unique names of lipid moleculed
  molec_names <- unique(class_lipids$LipidMolec)
  # retrieve duplicated molecule names
  unique_molec_names <- molec_names[which(class_size>1)]
  # get duplicated lipid molecules
  duplicates <- class_lipids[which(class_lipids$LipidMolec %in% unique_molec_names), ] %>% 
    arrange(Class, LipidMolec, BaseRt)
  return(duplicates)
}


# detect same lipid molecules with different retention time
# Filter the lipid molecule contains same name but different retention time based on your criteria  
# all the dupliated molecules will be stored in dupliates.csv
filter_duplicate <- function(duplicate_molecs, data, selections){
  Class <- syms(selections)[[1]]
  LipidMolec <- syms(selections)[[2]]
  BaseRt <- syms(selections)[[3]]
  MainIon <- syms(selections)[[4]]
  if(nrow(duplicate_molecs) != 0){
    message("\n!!!Attention: Identical lipid molecules with multiples retention time. Please note that the duplicate lipid molecules are stored in reserved_duplicates.csv \n !!!!!! Potential sample contamination. \n To PROCEED, pick one: ")
    write_csv(duplicate_molecs, "data/QC/duplicated.molecules.csv")
    # diff_baseRT will store the retention time differences for the same lipid molecule
    diff_baseRT <- duplicate_molecs %>% 
      select(LipidMolec, BaseRt, MainIon) %>% 
      group_by(LipidMolec, MainIon) %>% 
      summarise(gap=round(max(BaseRt)-min(BaseRt), 2))
    message("\nDifferences in retention time for identical lipid molecule are stored under diff_RT.csv")
    write_csv(diff_baseRT, "data/QC/diff_RT.csv")
    # two methods to fix duplicate lipid molecules
    message("\n A: use only ONE lipid molecule with largest main area under curve, OR \n B: Summation of main area under curve of ALL identical lipid molecule.")
    fix_method <- readline("ENTER 'A' or 'B': ")
    fixed_dt<- fix_duplicate(data, duplicate_molecs, fix_method, selections)
    # transform duplicate lipid molecules
    duplicate_output <- fixed_dt[[1]]
    # store reserved duplicate lipid molecules
    write_csv(duplicate_output, "data/reserved_duplicates.csv")
    # get filtered lipid information within processed duplicate lipid molecules
    data <- fixed_dt[[2]]
    message("Filtered lipid molecules sans duplicates are stored under removeduplicates.csv")
    write_csv(data, "data/QC/rm_duplicates.csv")
    return(data)
  } else{
    message("There is no duplicate molecules in your data")
  }
}
########################################################################
# function name: fix_duplicate
# parameter: fix_method, prelipids
# utility: fixing the same lipid molecule with different retention time
#           1) pick the lipid molecule which summation of main area is the 
#             biggest one;
#           2) aggregate the two different retention time lipid molecules
########################################################################
fix_duplicate <- function(filtered_lipids, duplicates, fixmethod, selections){
  Class <- syms(selections)[[1]]
  LipidMolec <- syms(selections)[[2]]
  BaseRt <- syms(selections)[[3]]
  MainIon <- syms(selections)[[4]]
  # get all names of sample
  all_names <- colnames(duplicates)[-c(1:4)]
  if(str_to_upper(fixmethod) == "A"){
    # get sum main areas of all samples for each duplicate lipid molecule
    sum_dt <- duplicates %>%  
      select(!!!syms(selections), all_of(all_names))  %>% 
      mutate(sumROC=rowSums(.[5:ncol(.)]))
    # get the lipid molecule with maximum value in duplicates
    sum_dt <- sum_dt %>% group_by(LipidMolec) %>% arrange(LipidMolec, sumROC) %>% top_n(1, sumROC)
    # get the lipid molecules which only keep bigger main Area value for duplicated molecules
    duplicate_dt <- semi_join(filtered_lipids, sum_dt, by = c("BaseRt", "LipidMolec", "MainIon")) 
    unique_dt <- anti_join(filtered_lipids, duplicates, by = c("BaseRt", "LipidMolec", "MainIon"))
    dt <- bind_rows(duplicate_dt, unique_dt) %>% arrange(LipidMolec)
  }else if(str_to_upper(fixmethod) == "B"){
    # aggregate all the duplicated molecules by type
    sum_dt <- duplicates %>%  
      select(Class, LipidMolec, all_of(all_names)) %>% 
      group_by(Class, LipidMolec) %>% 
      summarise_at(all_names, list(~sum(.)))
    # retrieve all column information for duplicates
    dt1 <- merge(sum_dt, filtered_lipids, all.x = TRUE) %>% select(colnames(filtered_lipids))
    # keep unchanged information like Rej, LipidMolec, Formula, etc. in the duplicated molecules
    dt_replace <- filtered_lipids %>% 
      select(Rej, LipidMolec, Class, contains("FA"), `Calc Mass`, Formula) %>% 
      group_by(LipidMolec)%>% 
      distinct()  %>% 
      filter(LipidMolec %in% sum_dt$LipidMolec)
    dt2 <- dt1 %>% 
      select(-!!colnames(dt_replace)) %>% 
      bind_cols(., dt_replace) %>% 
      select(colnames(filtered_lipids))
    # retrieve raw unduplicated molecule 
    dt3 <- filtered_lipids %>% 
      filter(! LipidMolec %in% unique(duplicates$LipidMolec))
    # bind processed duplicated with unduplicated lipid molecules
    dt <- bind_rows(dt2, dt3) %>% arrange(LipidMolec)
  } else{
    fixmethod <- readline("You type wrong, please type again, A/B: ")
    return(fix_duplicate(filtered_lipids, duplicates, fixmethod))
  }
  return(list(sum_dt, dt))
}



########################################################################
# function name: Input
# parameter: data (filtered_lipids)
# utility: input group information of samples from the console and reterieve 
#           corresponding sample information from the csv file
########################################################################
Input <- function(data){
  ###  making group
  message("\nProvide infomation of experimental groups")
  # info about group size
  ngroups <- readline("How many experimental groups: ") %>% as.numeric(.)
  ngroups <- check_group_number(ngroups)
  # Type the group info and store it into a variable
  total_info <- InputGroups(ngroups)
  # store the sample and index in a list
  info <- GrepIndex(total_info, data)
  message("Take a look at the sample info and its column position information in the file below")
  glimpse(info)
  # retrieve the group names
  group_names <- dimnames(info)[[2]]
  
  return(list(info, group_names, ngroups))
}


########################################################################
# function name: check_group
# parameter: x (ngroups from Input function/group number)
# utility: check if the group number is numeric and store the correct form
#         of group number information
########################################################################
check_group_number <- function(x){
  if((x%%1 != 0) | (is.na(x))){
    print("Group number must be numeric!")
    x <- readline("Group number: ") %>% as.numeric(.)
    return(check_group_number(x))
  }else{
    if(x<=0){
      print("You input wrong group number, please retype")
      x <- readline("Group number: ") %>% as.numeric(.)
      return(check_group_number(x))
    }else{
      return(x)
    }
  }
}




########################################################################
# function name: InputGroups
# parameter: n (ngroups from input function/group number)
# utility: read standard input from console and check if it's correct and
#           store group info into list
########################################################################
InputGroups <- function(n){
  group_info  <- c()
  info <- c()
  total_info <- list()
  for(i in 1:n){
    ms1                   <- paste("Description for Group ", i, " (name): ")
    group_info[i]         <- readline(prompt = ms1)
    ms2                   <- paste("Which samples assigned to Group ", i, "(sample number, e.g. s1 s2 s3 ): ")
    info[i]               <- readline(prompt = ms2)
    total_info[[group_info[i]]] <- info[i]
  }
  message("CONFIRM the group information below")
  glimpse(total_info)
  message("Do you want to edit group infomation? ")
  group_condition <- readline("Y/N: ")
  if(str_to_lower(group_condition) !="n"){
    if(str_to_lower(group_condition) =="y"){
      return(InputGroups(n))
    }else{
      print("ERROR, please enter Y/N.")
      return(InputGroups(n))
    }
  }
  return(total_info)
}




# Preparation for pair-wise correlations
# source from stackoverflow 
inf2NA      <- function(x) { x[is.infinite(x)] <-  NA; x }

panel.cor   <- function(x, y, digits = 2, cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient
  r   <- cor(inf2NA(x), inf2NA(y), use = "pairwise.complete.obs", method = "pearson")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste("r= ", txt, sep = "")
  text(0.5, 0.6, txt)
  
  
  #############################################################################################################
  # p-value calculation
  p                       <- cor.test(x, y)$p.value
  txt2                    <- format(c(p, 0.123456789), digits = digits)[1]
  txt2                    <- paste("p= ", txt2, sep = "")
  
  # edited by WL. Need more work for this part
  
  if((is.na(p<0.01) | p < 0.01)) txt2  <- paste("p= ", "<0.01", sep = "")
  text(0.5, 0.4, txt2)
  #############################################################################################################
  
}

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h   <- hist(x, plot = FALSE)
  breaks  <- h$breaks; nB <- length(breaks)
  y       <- h$counts; y  <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "#4393C3", ...)
}

# this part reformats the axis label positions on same sides. Source: stackflow
pairs2 <- function (x, labels, panel = points, ..., 
                    lower.panel = panel.smooth, diag.panel=panel.hist, 
                    upper.panel = panel.cor, text.panel = textPanel, 
                    label.pos = 0.5 + has.diag/3, cex.labels = NULL, font.labels = 1, 
                    row1attop = TRUE, gap = 1) {
  textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x, 
                                                               y, txt, cex = cex, font = font)
  localAxis <- function(side, x, y, xpd, bg, col = NULL, main, 
                        oma, ...) {
    if (side%%2 == 1) 
      Axis(x, side = side, xpd = NA, ...)
    else Axis(y, side = side, xpd = NA, ...)
  }
  localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
  localLowerPanel <- function(..., main, oma, font.main, cex.main) lower.panel(...)
  localUpperPanel <- function(..., main, oma, font.main, cex.main) upper.panel(...)
  localDiagPanel <- function(..., main, oma, font.main, cex.main) diag.panel(...)
  dots <- list(...)
  nmdots <- names(dots)
  if (!is.matrix(x)) {
    x <- as.data.frame(x)
    for (i in seq_along(names(x))) {
      if (is.factor(x[[i]]) || is.logical(x[[i]]))
        x[[i]] <- as.numeric(x[[i]])
      if (!is.numeric(unclass(x[[i]])))
        stop("non-numeric argument to 'pairs'")
    }
  }else if(!is.numeric(x)){
    stop("non-numeric argument to 'pairs'")
  }
  
  panel <- match.fun(panel)
  if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel))
    lower.panel <- match.fun(lower.panel)
  if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel))
    upper.panel <- match.fun(upper.panel)
  if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel))
    diag.panel <- match.fun(diag.panel)
  if (row1attop) {
    tmp <- lower.panel
    lower.panel <- upper.panel
    upper.panel <- tmp
    tmp <- has.lower
    has.lower <- has.upper
    has.upper <- tmp
  }
  nc <- ncol(x)
  if (nc < 2)
    stop("only one column in the argument to 'pairs'")
  has.labs <- TRUE
  if (missing(labels)) {
    labels <- colnames(x)
    if (is.null(labels))
      labels <- paste("var", 1L:nc)
  }else if (is.null(labels))
    has.labs <- FALSE
  oma <- if ("oma" %in% nmdots){
    dots$oma
  }else NULL
  main <- if ("main" %in% nmdots){
    dots$main
  }else NULL
  if (is.null(oma)) {
    oma <- c(4, 4, 4, 4)
    if (!is.null(main))
      oma[3L] <- 6
  }
  opar <- par(mfrow = c(nc, nc), mar = rep.int(gap/2, 4), oma = oma)
  on.exit(par(opar))
  dev.hold()
  on.exit(dev.flush(), add = TRUE)
  for (i in if (row1attop)
    1L:nc
    else nc:1L) for (j in 1L:nc) {
      localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE,
                type = "n", ...)
      if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
        box()
        # edited here...
        #           if (i == 1 && (!(j%%2) || !has.upper || !has.lower))
        #           localAxis(1 + 2 * row1attop, x[, j], x[, i],
        #                       ...)
        # draw x-axis
        if (i == nc & j != nc)
          localAxis(1, x[, j], x[, i],
                    ...)
        # draw y-axis
        if (j == 1 & i != 1)
          localAxis(2, x[, j], x[, i], ...)
        #           if (j == nc && (i%%2 || !has.upper || !has.lower))
        #             localAxis(4, x[, j], x[, i], ...)
        mfg <- par("mfg")
        if (i == j) {
          if (has.diag)
            localDiagPanel(as.vector(x[, i]), ...)
          if (has.labs) {
            par(usr = c(0, 1, 0, 1))
            if (is.null(cex.labels)) {
              l.wid <- strwidth(labels, "user")
              cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
            }
            text.panel(0.5, label.pos, labels[i], cex = cex.labels,
                       font = font.labels)
          }
        }
        else if (i < j)
          localLowerPanel(as.vector(x[, j]), as.vector(x[, i]), ...)
        else localUpperPanel(as.vector(x[, j]), as.vector(x[, i]), ...)
        if (any(par("mfg") != mfg))
          stop("the 'panel' function made a new plot")
      }
      else par(new = FALSE)
    }
  if (!is.null(main)) {
    font.main <- if ("font.main" %in% nmdots)
      dots$font.main
    else par("font.main")
    cex.main <- if ("cex.main" %in% nmdots)
      dots$cex.main
    else par("cex.main")
    mtext(main, 3, 3, TRUE, 0.5, cex = cex.main, font = font.main)
  }
  invisible(NULL)
}




######################################################################

########################################################################
# function name: GrepIndex
# parameter: sample, data
# utility: get sample names and its column position information from filtered data csv file
########################################################################
GrepIndex <- function(sample, data){
  sample.names  <- list()
  group_names   <- c()
  col.index     <- list()
  for( i in 1:length(names(sample))){
    group_names[i]                  <- names(sample)[i]
    group_sample_names              <- sample[[i]] %>%
      strsplit(., "\\s+") %>%
      unlist() %>%
      str_to_lower()
    sample.names[[group_names[i]]]  <- paste("MainArea[", group_sample_names, "]", sep="")
    col.index[[group_names[i]]]     <- which(colnames(data) %in% sample.names[[i]])
  }
  info      <- rbind(sample.names, col.index)
  return(info)
}


########################################################################
# function name: PCA_pairs_Plot
# parameter: info, group_names, filtered_lipids, mark
# utility: this function will plot PCA and pair-wise correlations among the samples.
#           This will also produce two kinds of correlation plots based on different
#           axis display style. 
#           Please notice that this function is really long format. 
########################################################################
PCA_pairs_Plot <- function(info, filtered_lipids, mark){
  # retrieve the sample info position 
  index         <- 1:length(info)
  sample_index  <- index[index %% 2 ==0]
  
  ######################################################################
  # QC PLOT 2 - Pair-wise correlation between replicates, reformatted axis
  for(i in 1:length(sample_index)){
    range     <- sample_index[i]
    plot_name <- paste("pairs.plot.", i, ".", mark, ".png",sep="")
    path      <- file.path("plot/", plot_name)
    pairs2(log10(filtered_lipids[, info[[range]]]), 
           lower.panel = panel.smooth, diag.panel=panel.hist, 
           upper.panel = panel.cor, xlim=c(4,12), ylim=c(4,12))
    dev.copy(png, path)
    dev.off()
  }
  
  #################################################################
  ### old axis for PCA
  # 
  # # Storing pairwise plot in the plot directory
  # for(i in 1:length(sample_index)){
  #   range     <- sample_index[i]
  #   plot_name <- paste("pairs.plot.", i, ".", mark, ".oldversion.pdf",sep="")
  #   path      <- file.path("plot/", plot_name)
  #   pairs(log10(filtered_lipids[, info[[range]]]), 
  #         lower.panel = panel.smooth, diag.panel=panel.hist, 
  #         upper.panel = panel.cor, xlim=c(4,12), ylim=c(4,12))
  #   dev.copy(pdf, path)
  #   dev.off()
  # }
  # 
  
  #################################################################
  # making group repeats according to its position for making groups of later PCA
  information <- retrieve_info(info)
  group_repeat <- information[[1]]
  sample_list <- information[[2]]
  sample_index <- information[[3]]
  
  
  
  # Formatting the table for PCA
  filtered_lipids_PCA <-  filtered_lipids %>% 
    select(all_of(sample_list)) %>% 
    t() %>% 
    as.data.frame()
  
  colnames(filtered_lipids_PCA) <- filtered_lipids$LipidMolec    
  
  # if just show the sample names by removing the MainArea prefix
  rownames(filtered_lipids_PCA) <- filtered_lipids_PCA %>% 
    rownames(.) %>% 
    str_remove_all(., "MainArea\\[") %>% str_remove_all(., "\\]")
  
  log2_filtered_lipids_PCA      <- log((filtered_lipids_PCA+1), 2)
  
  # making group info for PCA 
  filtered_lipids_PCA$Group       <- group_repeat
  log2_filtered_lipids_PCA$Group  <- group_repeat
  
  # Perform PCA 
  res.pca         <-  PCA(log2_filtered_lipids_PCA, scale.unit=TRUE, 
                          quali.sup=ncol(filtered_lipids_PCA), graph=FALSE)
  
  var <- get_pca_var(res.pca)
  pvar <- fviz_pca_var(res.pca, col.var = "contrib",
                       gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                       geom.var = c("point", "text"))
  print(pvar)
  
  p1 <- fviz_screeplot(res.pca, ncp=10)
  concat          <-  cbind.data.frame(log2_filtered_lipids_PCA[, ncol(log2_filtered_lipids_PCA)], res.pca$ind$coord)
  ellipse.coord   <-  coord.ellipse(concat, bary=TRUE)
  print(p1)
  
  p2 <- plot.PCA(res.pca, habillage=ncol(log2_filtered_lipids_PCA), 
                 ellipse=ellipse.coord, cex=0.8, label="all")
  print(p2)
  pca_name <- addquotes("plot/sample.pca.", !!mark, ".png")
  dev.copy(png, pca_name)
  dev.off()
  ###################################### 
  return(list(sample_list, group_repeat))
}



########################################################################
# function name: PCAcheck
# parameter: pca_check
# utility: check if deleting samples needed and replot PCA
########################################################################
PCAcheck <- function(pca_check, data){
  if(str_to_lower(pca_check)!="n"){
    if(str_to_lower(pca_check)=="y"){
      samples <- Input(data)
      info <- samples[[1]]
      #group_name <- samples[[2]]
      # ngroups <- samples[[3]]
      label <- "new"
      info_list <-  PCA_pairs_Plot(info, data, label)
    }else{
      print("You typed wrong, please retype.")
      pca_check <- readline("Please type Y/N: ")
      return(PCAcheck(pca_check, data))
    }
  }
  return(info_list)
}  


########################################################################
# function name: retrieve_info
# parameter: info
# utility: retrieve group and its sample information
########################################################################
retrieve_info <- function(info){
  # retrieve the sample info position 
  index         <- 1:length(info)
  sample_index  <- index[index %% 2 ==0]
  # making group repeats according to its position for making groups of later PCA
  index_list  <- index[!index %% 2 ==0]
  sample_list <- c()
  for(i in length(index_list):1){
    j           <- index_list[i]
    sample_list <- c(info[[j]], sample_list)
  }
  group_repeats <- c()
  group_name <- dimnames(info)[[2]]
  for(i in length(group_name):1){
    k             <- index_list[i]
    len           <- length(info[[k]])
    repeats       <- rep(group_name[i], len)
    group_repeats <- c(repeats, group_repeats)
  }
  return(list(group_repeats, sample_list, sample_index))
}


#########################################################################
# function name: count_pathway
# parameter: x (a list or vector of lipid molecules)
# utility: find the saturation pattern and its counts
########################################################################
count_pathway <- function(x){
  pattern_name <- i <- j <- k <- c()
  i <- str_count(x, ":0")
  j <- str_count(x, ":1")
  k <- str_count(x, ":[2-9]")
  SFA <- paste(rep("SFA", i), collapse = "/")
  MUFA <- paste(rep("MUFA", j), collapse = "/")
  PUFA <- paste(rep("PUFA", k), collapse = "/")
  pattern_name <- paste( SFA, MUFA, PUFA, sep = "/")
  patterns <- data.frame(pattern= pattern_name, SFA=i, MUFA=j, PUFA=k, stringsAsFactors = FALSE) %>% as.matrix(.)
  return(patterns)
}


###############################################################
# transform data attributes to numeric type
###############################################################
transform_to_numeric <- function(x){class(x) <- as.numeric(x)}


########################################################################
# function name: calc_group
# parameter: fa_percent, groups
# utility: find mean or median of samples in each group,
#         return saturantion percentage information, mean and median in each group
########################################################################
calc_group <- function(percent_info, groups){
  #  group_list <- unique(groups[, 2]) %>% unlist()
  all_samples <- unique(groups[, 1]) %>% unlist()
  names(all_samples) <- all_samples
  dt <- percent_info %>%
    mutate_at(all_samples, list(SFA = ~.*`%SFA`, MUFA = ~.*`%MUFA`, PUFA = ~.*`%PUFA`))
  
  dd <- dt %>%
    select(Class, contains("MainArea"), -ends_with("]")) %>%
    group_by(Class) %>%
    summarise_all(sum)
  
  dd_mean <- cal_class(dd, groups, mean, "mean")
  dd_median <- cal_class(dd, groups, median, "median")
  dd_sd <- cal_class(dd, groups, sd, "sd")
  dd_se <- cal_class(dd, groups, se, "se")
  data <- cbind(dd, dd_median, dd_mean, dd_se)
  return(list(dt, dd, data))
}


########################################################################
# function name: calc_class
# parameter: data, info, fun, pick
# utility: auxiliary function for calc_group, which passing percentage pattern 
#         informaion, group information, methods (fun, pick) into the function to 
#         get corresponding statistics for methods
########################################################################
cal_class <- function(data, info, fun, pick){  
  group_list <- unique(info[, 2]) %>% unlist()
  dt1 <- data.frame(row.names = rownames(data))
  for (i in group_list){
    sample_list <- filter(info, groups == i) %>% ungroup() %>% select(samples) %>% unlist()
    # bulid names
    name_group <- addquotes(!!i, "_", !!pick)
    post_fix <- c("_SFA", "_MUFA", "_PUFA")
    names <- c()
    for(i in post_fix) names[i] <- addquotes(!!name_group, !!i)
    samples <- list()
    for(j in post_fix) samples[[j]] <- addlists(sample_list, j)
    dt <- Fun_statistics(data, names, fun, samples)
    dt1 <- cbind(dt1, dt)
  }
  return(dt1)
}


########################################################################
# function name: Fun_statisitcs
# parameter: data, names, fun, samples
# utility: independent or auxiliary function for calculation for selected samples by row
# ########################################################################
# Fun_statistics <- function(data, names, fun, samples){
#   dts <- data.frame(row.names = rownames(data))
#   # get median, mean and sd for samples in each group 
#   for(i in seq_along(names)){
#     dt <-  data %>% 
#       rowwise() %>% 
#       transmute(!!(names[i]) := !!fun(c(!!!syms(samples[[i]])), na.rm = TRUE))
#     dts <- cbind(dts, dt)
#   }
#   return(dts)
# }








########################################################################
# function name: Fun_statisitcs
# parameter: data, names, fun, samples
# utility: independent or auxiliary function for calculation for selected samples by row
#########################################################################
Fun_statistics <- function(data, names, fun, samples){
  dts <- data.frame(row.names = rownames(data))
  # get median, mean and sd for samples in each group 
  for(i in seq_along(names)){
    dt <-  data %>% 
      rowwise() %>% 
      transmute(!!(names[i]) := eval(expr((!!fun)(c(!!!syms(samples[[i]])), na.rm = TRUE))))
    
    dts <- cbind(dts, dt)
  }
  return(dts)
}

########################################################################
# function name: cal_lipid_statisitics
# parameter: data, group_information, method, type
# utility: independent or auxiliary function for calculation for selected samples by row
#         return the output and its formated style
#########################################################################
cal_lipid_statistics <- function(data, group_information, method,  type){
  selections <- syms(type)
  type.name <- data %>% select(!!!selections) %>% ungroup()
 # type.name <- data %>% select(!!type)
  group_names <- group_information$groups %>% unique() %>% unlist()
  group_samples <- list()
  post_name <- addquotes("_", !!method)
  
  # get group information for samples
  for(i in group_names) group_samples[[i]] <- subset(group_information, groups == i) %>% ungroup() %>% select(samples) %>% unlist()
  # calculate the mean for each group of each class
  lipid_group <- Fun_statistics(data, group_names, method, group_samples) %>% 
    cbind(Class = type.name, .)
  colnames(lipid_group) <- c(type, group_names)
  # add post name, e.g. _mean
  colnames(lipid_group)[-c(1:length(type))] <- colnames(lipid_group)[-c(1:length(type))] %>% addlists(., post_name)
  # reformat the data
  data_group <- FormatData(lipid_group, method, type)
  # # remove the post name, e.g. _mean
  # data_group <- data_group %>% mutate_at(vars(TYPE, Groups), list(~str_remove_all(., post_name)))
  data_group$Groups <- factor(data_group$Groups, levels = group_names)
  return(list(lipid_group, data_group))
}





########################################################################
# function name: pick_control
# parameter: data, pick, info, control
# utility: get control group information for later calculation
########################################################################
pick_control <- function(data, pick, info, control){
  groups <- unique(info[, 2]) %>% unlist()
  control_info <- addquotes(!!control, "_", !!pick)
  selected_data <- data %>% 
    select(Class, contains(!!control_info), contains("MainArea"))
  return(selected_data)
}



# se <- function(x)sd(x,na.rm=TRUE)/sqrt(length(x[complete.cases(x)]))







#######################################################################
# function name: norm_by_mean
# parameter: data, pattern, sample_list, control_names
# utility: get mean/median of the control group, and normalized every sample by selected value
########################################################################
norm_by_mean <- function(data, pattern, sample_list, control_names){
  dt <- data %>% 
    select(Class, contains(!!pattern)) %>% 
    group_by(Class) %>% 
    transmute_at(sample_list , list(~./!!sym(control_names)))
  colnames(dt)[-1] <- addlists(colnames(dt)[-1], pattern)
  return(dt)
} 

#######################################################################
# function name: norm_by_mean2
# parameter: data, ll, sample_list, control_names, data_type
# utility: ll is the group variable, e.g Class, LipidMolec.
#         get mean/median of the control group, and normalized every sample by 
#         the selected value. The normalization method depends on if the data 
#         type is log transformed or not.
########################################################################
norm_by_mean2 <- function(data, ll, sample_list, control_names, data_type){
  if(data_type != "y"){
    dt <- data %>% 
      group_by(!!!syms(ll)) %>% 
      transmute_at(sample_list , list(~./!!sym(control_names)))
    return(dt)
  }else{
    dt <- data %>% 
      group_by(!!!syms(ll)) %>% 
      transmute_at(sample_list , list(~.-!!sym(control_names)))
    return(dt)
  }
  
} 




#######################################################################
# function name: FC_fun
# parameter: group_info, pick, norm_samples, method
# utility: calculate the fold change for each normalized group
########################################################################
FC_fun <- function(group_info, pick, norm_samples){
  dt <- data.frame(row.names = rownames(norm_samples))
  groups <- unique(group_info$groups)
  for(i in groups){
    sample_list <- filter(group_info, groups == i) %>% ungroup() %>% select(samples) %>% unlist()
    name_group <- addquotes(!!i, "_", !!pick)
    post_fix <- c("_SFA", "_MUFA", "_PUFA")
    names <- c()
    for(j in post_fix) names[j] <- addquotes(!!name_group, !!j)
    samples <- list()
    for(k in post_fix) samples[[k]] <- addlists(sample_list, k)
    dt1 <- Fun_statistics(norm_samples, names, pick, samples)
    dt <- cbind(dt, dt1)
  }
  return(dt)
} 


########################################################################
# function name: check_group
# parameter: groups
# utility: check if the user type correct group name
########################################################################
check_group <- function(groups, name){
  gr <- ifelse(name =="", 
               readline(paste("Enter the name of group: ")), 
               readline(paste("Enter the name(s) of", name, "group: ")) )%>%  
    str_split(., "\\s") %>% 
    unlist()
  if(length(gr) > length(groups)){
    message("Exceed the length of groups!")
    gr <- check_group(groups, name)
  }
  if(!all(gr %in% groups)){
    message("Name cannot be identified, please type again.")
    gr <- check_group(groups, name)
  } 
  return(gr)
}


########################################################################
# function name: FormatData
# parameter: data, pick
# utility: reformat the data to long formats for plotting
########################################################################
FormatData <- function(data, pick, fixed_cols){
  dt <- data %>% select(!!!syms(fixed_cols), contains(pick))
  pattern <- addquotes("_", !!pick)
  dt <- dt %>% gather("Groups", !!pick, -all_of(fixed_cols))
  dt$Groups <- str_remove_all(dt$Groups, pattern)
  return(dt)
}






########################################################################
# function name: FormatData2
# parameter: data, pick
# utility: reformat the data to long formats for plotting
########################################################################
FormatData2 <- function(data, pick){
  dt <- data %>% select(Class, contains(pick))
  pattern1 <- addquotes("_", !!pick, "_.*")
  pattern2 <- addquotes(".*", !!pick, "_")
  dt <- dt %>% gather("TYPE", !!pick, -Class)
  dt$Groups <- str_replace_all(dt$TYPE, pattern1, "")
  dt$TYPE <- str_replace_all(dt$TYPE, pattern2, "")
  return(dt)
}
##
########################################################################
# function name: PlotTypes
# parameter: data, m (x, y and color parameters for plots)
# utility: passing data and 3 plot parameters for plot
########################################################################
PlotTypes <- function(data, m){
  m <- syms(m)
  p1 <- m[[1]]
  p2 <- m[[2]]
  p3 <- m[[3]]
  min_y <- data %>% select(addquotes(!!p2)) %>% min()
  max_y <- data %>% select(addquotes(!!p2)) %>% max()
  plot_types <-   ggplot(data, aes(x=eval(p1), y=eval(p2), fill = eval(p3))) +
    # facet_wrap(~Class, scales = "free") +
    scale_fill_manual(values = clPalette1) +
    theme_bw()+
    set_theme()
  return(plot_types) 
} 


########################################################################
# function name: subtract_background
# parameter: data (filtered_lipidomics)
# utility: total class lipids of raw data by subtracting background (blank sample c) from all samples,
#         return subtracted information and its mean and sd statistics data
########################################################################
subtract_background <- function(data, sample_list, subtraction){
  # total class lipids of raw data by subtracting background (blank sample c) from all samples
  data <- data %>% 
    mutate_at(vars(all_of(sample_list)), list(~.-eval(sym(!!subtraction)))) %>%  
    mutate_at(vars(!!subtraction), list(~ifelse(.==0, ., 0)))
  return(data)
}


########################################################################
# function name: cal_lipid_statistics
# parameter: data, group_information, method, pick, type
# utility: calculate picked method (e.g. mean, sd, median, sum..) of each group for each class type.
#           The type could be lipid Class or lipid molecule colname name
########################################################################
# cal_lipid_statistics <- function(data, group_information, method, pick, type){
#   type.name <- data %>% select(!!type)
#   group_names <- group_information$groups %>% unique() %>% unlist()
#   group_samples <- list()
#   post_name <- addquotes("_", !!pick)
#   
#   # get group information for samples
#   for(i in group_names) group_samples[[i]] <- subset(group_information, groups == i) %>% ungroup() %>% select(samples) %>% unlist()
#   # calculate the mean for each group of each class
#   lipid_group <- Fun_statistics(data, group_names, method, group_samples) %>% 
#     cbind(Class = type.name, .)
#   colnames(lipid_group)[1] <- "Class"
#   # add post name, e.g. _mean
#   colnames(lipid_group)[-1] <- colnames(lipid_group)[-1] %>% addlists(., post_name)
#   # reformat the data
#   data_group <- FormatData(lipid_group, pick)
#   # remove the post name, e.g. _mean
#   data_group <- data_group %>% mutate_at(vars(TYPE, Groups), list(~str_remove_all(., post_name)))
#   return(list(lipid_group, data_group))
# }














########################################################################
# function name: ClassPlot
# parameter: data, list of parameters (e.g. x, y, color)
# utility: make color bar plots for the data
########################################################################
ClassPlot <- function(data, parameters){
  parameters <- syms(parameters)
  p1 <- parameters[[1]]
  p2 <- parameters[[2]]
  p3 <- parameters[[3]]
  p4 <- parameters[[4]]
  ranges <- data %>% 
    mutate(limits=sign(eval(p2))* (abs(eval(p2))+eval(p4))) %>% 
    summarise(max(limits), min(limits))
  limits <- sapply(ranges, function(x)sign(x)*10^(ceil(log10(abs(x)), 0)))
  classplot <- ggplot(data, aes(x=reorder(eval(p1), eval(p2)), fill=eval(p3), y=eval(p2))) +
    geom_errorbar(aes(ymin=sign(eval(p2))*abs(eval(p2)), ymax=sign(eval(p2))* (abs(eval(p2))+eval(p4))),
                  position=position_dodge(width = .9, preserve = "single"), 
                  width = 0.1,
                  size=.2) +
    theme_bw() +
    set_theme()+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.y = element_text(angle = 30, hjust=1)
    ) + 
    # scale_y_log10(
    #    breaks = scales::trans_breaks("log10", function(x) 10^(x)),
    #    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    #   ) +
    # scale_y_continuous(trans = 'log10',
    #                    breaks = trans_breaks('log10', function(x) 10^x),
    #                    labels = trans_format('log10', math_format(10^.x))) +
    
    scale_y_continuous(trans = tn,
                       breaks = c(0, 10^3, 10^6, 10^9, 10^12),
                       labels = c(0, bquote(10^3), bquote(10^6), bquote(10^9), bquote(10^12)),
                       expand = c(0, 0, 0.1, 0)
    ) +
    expand_limits(y = limits) +
    scale_fill_manual(values = clPalette1) 
  return(classplot)
}




















geom_mean <- function(){
  list(
    stat_summary(fun.y = "mean", geom = "bar", fill = "grey70"),
    stat_summary(fun.data = "mean_cl_normal", geom = "errorbar", width = 0.4)
  )
}

ggplot(mpg, aes(class, cty)) + geom_mean()




ClassPlot2 <- function(data, parameters, se){
  parameters <- syms(parameters)
  p1 <- parameters[[1]]
  p2 <- parameters[[2]]
  p3 <- parameters[[3]]
  p4 <- parameters[[4]]
  ranges <- data %>% 
    mutate(limits=sign(eval(p2))* (abs(eval(p2))+eval(p4))) %>% 
    summarise(max(limits), min(limits))
  limits <- sapply(ranges, function(x)sign(x)*10^(ceil(log10(abs(x)), 0)))
  classplot <-  ggplot(data, aes(x=reorder(eval(p1), eval(p2)), fill=eval(p3), y=eval(p2))) +
    theme_bw() +
    set_theme()+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.y = element_text(angle = 30, hjust=1)) + 
    
    scale_y_continuous(trans = tn,
                       breaks = c(0, 10^3, 10^6, 10^9, 10^12),
                       labels = c(0, bquote(10^3), bquote(10^6), bquote(10^9), bquote(10^12)),
                       expand = c(0, 0, 0.1, 0)) +
    expand_limits(y = limits) +
    scale_fill_manual(values = clPalette1) +
    geom_se(parameters, se)
  return(classplot)
}












#########################################################################
# function name: detect_invalid
# parameter: data, group information
# This function will count 0 and negative values and their percentage in each 
# group based on he group information
#########################################################################
detect_invalid <- function(data, group_information){
  data <- as.data.frame(data)
  sample_list <- group_info$samples %>% unique() %>% unlist()
  dt1 <- data.frame(row.names = 1:nrow(data))
  for(i in unique(group_information$groups)){
    samples <- subset(group_information, groups == i, samples) %>% unlist()
    names <- c()
    names[1] <- paste0(i, "_empty")
    names[2] <- paste0(i, "_empty_percent")
    names[3] <- paste0(i, "_negative")
    names[4] <- paste0(i, "_negative_percent")
    data <- data %>% 
      mutate(!!names[1] := rowSums(.[which(colnames(.) %in% samples)] == 0), 
             !!names[3] := rowSums(.[which(colnames(.) %in% samples)] < 0)) %>%
      mutate(!!names[2] := eval(sym(!!names[1]))/length(samples),
             !!names[4] := eval(sym(!!names[3]))/length(samples)) %>% 
      select(colnames(data), all_of(names))
    
    # transform all values in a groun into NA if the negative percentage is over 50%
    # transform all negative values into NA
    dt2 <- data %>% transmute_at(vars(samples), list(~ifelse(eval(sym(!!names[4])) > 0.5 | . < 0 , NA, .)))
    dt1 <- cbind(dt1, dt2)
  }
  # transform all negative values into NA and format percentage value
  dt3 <- data %>%
    rowwise()%>%
    mutate_at(vars(sample_list), list(~ifelse(. < 0, NA, .))) %>%
    mutate_at(vars(contains("percent")), list(~scales::percent(.)))
  
  # get information for all values into NA in a group
  dt4 <- data %>% select(-all_of(sample_list)) %>% cbind(., dt1) %>% select(colnames(data))
  return(list(data, dt3, dt4))
}






plot_fc <- function(data, parameters, se){
  if(length(parameters) == 3){
    print("HAHA")
    paras <- syms(parameters)
    p <- ggplot(data, aes(x = eval(paras[[1]]), y = eval(paras[[2]]), fill = eval(paras[[3]]))) +
      theme_bw() +
      set_theme() +
      theme(axis.line = element_line(colour = "black", 
                                     size = .15),
            plot.title = element_text(hjust = 0.5),
            axis.text.y = element_text(angle = 30, hjust=1, size = 6),
            axis.text.x = element_text(size = 6),
            axis.ticks = element_line(size = .8),
            axis.ticks.length = unit(1, "mm")) +
      # scale_y_continuous(expand = c(0, 0, 0.2, 0)) +
      # ylim(limits) +
      scale_fill_manual(values = clPalette1)
    }else{
      print("yoyo")
    paras <- syms(parameters)
    ranges <- data %>% 
      mutate(bounds=sign(eval(paras[[2]])) * (abs(eval(paras[[2]])) + abs(eval(paras[[4]])))) %>% 
      summarise(max(bounds), min(bounds))
    limits <- ranges
    # limits[[2]] <- limits[[2]] - 10^5
    p <- ggplot(data, aes(x = reorder(eval(paras[[1]]), eval(paras[[2]])), 
                          y = eval(paras[[2]]), fill = eval(paras[[3]]))) +
      #   expand_limits(y = limits) +
      theme_bw() +
      set_theme() +
      theme(axis.line = element_line(colour = "black", 
                                     size = .2),
            plot.title = element_text(hjust = 0.5),
            axis.text.y = element_text(angle = 30, hjust=1, size = 6),
            axis.text.x = element_text(size = 6),
            axis.ticks = element_line(size = .8),
            axis.ticks.length = unit(1, "mm")) +
      # scale_y_continuous(expand = c(0, 0, 0.2, 0)) +
      # ylim(limits) +
      scale_fill_manual(values = clPalette1) +
      geom_se(parameters, se)
    } 
    axis_st <-  data %>% filter(eval(paras[[2]]) < 0) %>% nrow()
    if(axis_st >0){
      p <- p + 
        scale_y_continuous(expand = c(0.02, 0, 0.2, 0)) +
        geom_hline(yintercept = 0, color = "black",size = 1, linetype = "dashed") +
        scale_fill_d3()
    }else{
      p <-  p + scale_y_continuous(expand = c(0, 0, 0.2, 0))
    }
}







############## migrate the part here or not?







#########################################################################
# function name: fix_invalid_by_choice
# parameter: data, invalid_data
# utility: open the csv (e.g. checkInvalid_raw.csv or checkInvalid.csv) which 
#         already mannually deleted the suspcious invalid lipid molecules. This 
#         function will deleted same invalid lipid molecules in the 
#         original big data filtered_lipidomics based on your standards.
#########################################################################
fix_invalid_by_choice <- function(data, filtered_data, invalid_data){
  deleted_molecules <- setdiff(invalid_data$LipidMolec, filtered_data$LipidMolec)
  data <- data %>% filter(!LipidMolec %in% deleted_molecules)
  return(data)
}


########################################################################
# function name: EachClassPlot
# parameter: data, n_control
# utility: calculate the fold changes among median samples in groups
#         if you are a mac user, please install quartz via 
########################################################################
EachClassPlot <- function(long_data, paras, computer){
  #quartz(type = "png")
  if(computer == "mac"){
    quartz() 
  }
  long_data <- long_data %>% separate(., LipidMolec, into=c("Class", "LipidMolec"), sep = "\\(")
  long_data$LipidMolec <- long_data$LipidMolec %>% str_remove_all(., "\\)")
  n_bar <- paras[[1]]
  n_groups <- paras[[2]]
  symbols <- syms(paras[[3]])
  molec <- symbols[[1]]
  plot_func <- paras[[4]]
  pars <- paras[[3]]
  post_name <- paras[[5]]
  labs_info <- paras[[6]]
  
  #list1 <- paras[[3]][[5]]
  #list2 <- list1
  #list2$expand <- as.call(quote(c(0.02, 0, 0.2, 0)))
  #print(list2)
  # select lipid class
  classes <- long_data %>% select(Class) %>% unlist() %>% unique()
  counts <- long_data %>% group_by(Class) %>% tally()
  
  #if(paras[[3]][2] == "mean"){}

  
  
 # p4 <- parameters[[4]]
  # list2 <- paras[[7]]
  print(length(symbols))
  #print("\n\n")
  for(i in 1:length(classes)){
    options(warn=-1)
    pick_class <- classes[i]
    print(pick_class)
    observations <- subset(counts, Class == classes[i], n) %>% unlist()
    class_data <- subset(long_data, Class == classes[i]) %>% arrange(., !!molec)
    if(observations <= n_bar){
     # print("1")
      axis_st <- class_data %>% filter(eval(symbols[[2]])<0) %>% nrow()
     # p1 <- plot_func(class_data, pars, se = FALSE) 
      #p1 <- plot_func(class_data, pars, se = TRUE) 
      if(length(symbols) < 4){
        p1 <- plot_func(class_data, pars) +
          plot_bars(bar.params = list(position=position_dodge(width = 0.6, preserve = "single"), width = 0.6)) 
      }else{
        p1 <- plot_func(class_data, pars, se = FALSE) +
          plot_bars(bar.params = list(position=position_dodge(width = 0.6, preserve = "single"), width = 0.6)) +
          geom_errorbar(aes(ymin=sign(eval(symbols[[2]]))*abs(eval(symbols[[2]])), ymax=sign(eval(symbols[[2]]))* (abs(eval(symbols[[2]]))+abs(eval(symbols[[4]])))),
                        position=position_dodge(width = .6, preserve = "single"),
                        width = 0.1,
                        size=.2)
        
      }
      p1 <- p1 +
        # plot_bars(bar.params = list(position=position_dodge(width = 0.6, preserve = "single"), width = 0.7)) +
        
        # plot_bars() +
        labs_info +
        #   add_scales(scale.params = list1) +
        ggtitle(pick_class) +
        # geom_se(se.params = list(position=position_dodge(width = .6, preserve = "single"),
        #                          width = 0.1,
        #                          size=.2)) +
        # geom_errorbar(aes(ymin=sign(eval(symbols[[2]]))*abs(eval(symbols[[2]])), ymax=sign(eval(symbols[[2]]))* (abs(eval(symbols[[2]]))+abs(eval(symbols[[4]])))),
        #             position=position_dodge(width = .6, preserve = "single"),
        #             width = 0.1,
        #              size=.2) +
        coord_flip() 
     
      print(p1)
      #ggsave(filename = paste(pick_class, ".png", sep=""), path = 'plot/classes', device = "png", width=15, height=15, dpi=300)
      ggsave(filename = paste(post_name, pick_class,".png", sep=""), path = 'plot/classes', device = "png")
      
    }else{
     # print("3")
      if(n_bar %% n_groups != 0){
        n_bar <- (n_bar%/% n_groups)*n_groups
      }
      nfacet <- observations %/% n_bar 
      for(k in 0:nfacet){
        if(k < nfacet){
       #   print("4")
          # bar ranges
          ranges <- (n_bar*k+1):(n_bar*(k+1))
          print(ranges)

          data <- class_data %>% slice(ranges)
          axis_st  <- data %>% filter(eval(symbols[[2]])<0) %>% nrow()
          #p2 <- plot_func(data, pars, se=TRUE)
          if(length(symbols) < 4){
            p2 <- plot_func(data, pars) +
              plot_bars(bar.params = list(position=position_dodge(width = 0.9, preserve = "single"), width = 0.9)) 
              
            
          }else{
            p2 <- plot_func(data, pars, se = TRUE) +
              plot_bars(bar.params = list(position=position_dodge(width = 0.9, preserve = "single"), width = 0.9)) 
            
          }
          p2 <- p2 + 
            # plot_bars(bar.params = list(position=position_dodge(width = 0.9, preserve = "single"), width = 0.9)) +
            labs_info +
            # add_scales(scale.params = list1) +
            
            ggtitle(pick_class) +
            coord_flip() 
        
        
          print(p2)
          ggsave(filename = paste(pick_class, ".", k+1, post_name,".png", sep=""), path = 'plot/classes', device = "png")
        }else{
        #  print("6")
          ranges <- (n_bar*k+1):observations
          print(ranges)
          data <- class_data %>% slice(ranges)
          axis_st <- data %>% filter(eval(symbols[[2]])<0) %>% nrow()
          #p2 <- plot_func(data, pars, se = FALSE)
          #p2 <- plot_func(data, pars, se = TRUE)
          #p2 <- #p2 + 
            # plot_bars(bar.params = list(position=position_dodge(width = 0.6, preserve = "single"), width = 0.4)) +
            # geom_se(se.params = list(position=position_dodge(width = .6, preserve = "single"),
            #                          width = 0.1,
            #                          size=.2)) +
            # plot_bars() +
            #  add_scales(scale.params = list1) +
            # geom_errorbar(aes(ymin=sign(eval(symbols[[2]]))*abs(eval(symbols[[2]])), ymax=sign(eval(symbols[[2]]))* (abs(eval(symbols[[2]]))+abs(eval(symbols[[4]])))),
            #               position=position_dodge(width = .6, preserve = "single"),
            #               width = 0.1,
            #               size=.2) +
          
          if(length(symbols) < 4){
            p2 <- plot_func(data, pars) +
              plot_bars(bar.params = list(position=position_dodge(width = 0.6, preserve = "single"), width = 0.6)) 
            
         
            
          }else{
            p2 <- plot_func(data, pars, se = FALSE) +
              plot_bars(bar.params = list(position=position_dodge(width = 0.6, preserve = "single"), width = 0.6)) +
              geom_errorbar(aes(ymin=sign(eval(symbols[[2]]))*abs(eval(symbols[[2]])), ymax=sign(eval(symbols[[2]]))* (abs(eval(symbols[[2]]))+abs(eval(symbols[[4]])))),
                            position=position_dodge(width = .6, preserve = "single"),
                            width = 0.1,
                            size=.2)
          }
          p2 <- p2 + 
            labs_info +
            ggtitle(pick_class) +
            coord_flip() 
          print(p2)
          ggsave(filename = paste(pick_class, ".", k+1, post_name,".png", sep=""), path = 'plot/classes', device = "png")
        }
      }
    } 
  }
   dev.off()
   options(device = "RStudioGD")
}

















########################################################################
# function name: cal_foldchange
# parameter: data, control, pick
# utility: calculate fold change for each class as auxiliary or independent function
# Please note that the lipid which value is 0, negative or NA for control group is deleted for analysis
########################################################################
cal_foldchange <-function(data, control, pick, data_type){
  # filter the control group which value is 0, negative, or NA
  ko_class <- data %>% filter(eval(sym(control)) <= 0 | is.na(eval(sym(control)))) %>% select(Class) %>% unlist() %>% addlists(., ", ")
  if(length(ko_class) != 0){
    ko_class[length(ko_class)] <- substr(last(ko_class), 1, nchar(last(ko_class)) - 2)
    message("\n\nPlease note that the control group of class: ", ko_class, " are 0 or invalid for its ", pick,  " value which are deleted for analysis.\n\n")
  }
  data <- data %>% filter(eval(sym(control)) > 0 & !is.na(eval(sym(control))))
  control_data <- data[, colnames(data) %in% control]
  if(data_type == "y"){
    all_fc <- data[, -1] - control_data    ##### for imputed data
    all_fc <- bind_cols(Class = data[, 1], all_fc)
    return(all_fc)
  }else{
    all_fc <- data[, -1]/control_data #### for un-imputed data
    all_fc <- bind_cols(Class = data[, 1], all_fc)
    # log_fc <- log2(data[, -1]) - log2(control_data) 
    # log_fc <- bind_cols(Class = data[, 1], log_fc)
    # return(list(all_fc, log_fc))
    return(all_fc)
  }
}

#######################################################################
# function name: cal_multiple_fc
# parameter: times, data, groups, pick
# utility: get fold change and log fold change plot and data for each contrast
########################################################################
# cal_multiple_fc <- function(times, data, groups, pick, data_type){
#   #for(i in 1:times){
#   control <- check_group(groups, "control")
#   control_name <- addquotes(!!control, "_", !!pick)
#   
#   fc_data <- cal_foldchange(data, control_name, pick, data_type )
#   fc <- fc_data
#   fc_long <- FormatData(fc, pick)
#   post_name <- addquotes("_", !!pick)
#   fc_long <- fc_long %>% mutate_at(vars(TYPE, Groups), list(~str_remove_all(., post_name)))
#   
#   contrast <- check_group(groups, "comparison")
#   # get fold change info
#   fc_data_long <- fc_long %>% filter(Groups %in% c(control, contrast)) 
#   fc_data_wide <- fc_data_long %>% select(-TYPE) %>% spread(., Groups, pick)
#   file1 <- addquotes("data/fc_", !!control_name, "_",  !!(as.character(i)), ".csv")
#   write_csv(fc_data_wide, path=file1)
#   
#   # transform cntrol group value to 1 for displaying on plot
#   #fc_data_long <- fc_data_long %>% mutate_at(vars(!!pick), function(x)ifelse(eval(expr(`$`(., !!pick))) == 0, 1, x))
#   return(fc_data_long)
# }





########################################################################
# function name: log2trans
# parameter: x
# utility: take log transformation of the data
########################################################################
log2trans <-  function(x, na.rm=FALSE){log(x)}

########################################################################
# function name: ReplaceInf
# parameter: x
# utility: replace all the zero values to negative values, such that the 
#         log transformed data of zero value will produce NaN other than -Inf
########################################################################
ReplaceInf <- function(x, na.rm=FALSE){x=replace(x, x == -Inf, NaN)}
ReplaceNa <- function(x, na.rm=FALSE){x=replace(x, is.na(x), 1)}
se <- function(x, na.rm = TRUE) sqrt(var(x)/length(x))

########################################################################
# function name: ImputeMinProb
# parameter: data, q, tune.sigma
# utility: this function applys the imputation of left-censored missing data 
#           by random draws from a Gaussian distribution centered in a minimal 
#           value. The minimal value observed is estimated as being the q-th 
#           quantile (e.g. q = 0.01) of the observed values in that sample. 
#           q is a scalar used to determine a low expression value to be used 
#           for missing data imputation. 0 < q < 1, in this case q should be 
#           set to a low value. The default value is q = 0.01. tune.sigma is 
#           a scalar used to control the standard deviation of the Gaussian 
#           distribution used for random draws. If the sd is overestimated, 
#           than 0 < sigma.coef < 1. The default value is tune.sigma = 1.
#           The function is obtained from imputeLCMD package
########################################################################
ImputeMinProb <- function (data, q = 0.01, tune.sigma = 1) 
{
  n_samples <-  dim(data)[2]
  n_observations <-  dim(data)[1]
  imputated_data <-  data
  min_samples <-  apply(imputated_data, 2, quantile, prob = q, na.rm = T)
  count_NAs <-  apply(!is.na(data), 1, sum)
  count_NAs <-  count_NAs/n_samples
  filtered_data <-  data[which(count_NAs > 0.5), ]
  data_sd <-  apply(filtered_data, 1, sd)
  sd_temp <-  median(data_sd, na.rm = T) * tune.sigma
  print(sd_temp)
  for (i in 1:(n_samples)) {
    input <-  rnorm(n_observations, 
                    mean = min_samples[i], 
                    sd = sd_temp)
    imputated_data[which(is.na(data[, i])), i] <- input[which(is.na(data[,i]))]
  }
  return(imputated_data)
}


########################################################################
# function name: VolPlot
# parameter: n_comparisons, fit, filtered_data
# utility: passing parameters to BuildContrast function to plot volcano graph
#         for number of n_comparisons
########################################################################
VolPlot <- function(n_comparisons, fit, filtered_data){
  replicate( n_comparisons, BuildContrast(fit, filtered_data))
}



check_lipid <- function(lipid){
  customized_lipid <- readline("Input the name of lipid class(es), e.g. Cer TG: ") %>% 
    str_split(., "\\s+") %>% 
    unlist()
  if(!all(customized_lipid %in% lipid)){
    message("Typed wrong, the lipid class should be exactly from lipid class list above.")
    customized_lipid <- check_lipid(lipid)
  }else{
    return(customized_lipid)
  }
}
########################################################################
# function name: BuildContrast
# parameter: fit, filtered_data
# utility: Build contrast matrix for limma 
########################################################################
BuildContrast <- function(fit, filtered_data){
  # contmatrix <- makeContrasts(readline("Enter groups names for comparison, spaced by 'vs', e.g. KOvsWT: "),
  #                             levels = design
  # )
  option <- readline("Enter groups names for comparison, spaced by 'vs', e.g. KO vs WT: ") %>% 
    strsplit(., paste0("(?i) ", "vs", "( |[!\",.:;?})\\]])"),perl=TRUE) %>% 
    unlist()
  
  if(!all(option  %in% group_names)){
    message("You typed wrong.")
    option <- readline("Enter groups names for comparison, spaced by 'vs', e.g. KO vs WT: ") %>% 
      strsplit(., paste0("(?i) ", "vs", "( |[!\",.:;?})\\]])"),perl=TRUE) %>% 
      unlist()
  }
  
  contrast_group <- paste(option[1], "-", option[2], sep="") 
  cmd <- paste("tmp <- makeContrasts(", contrast_group, ", levels =
fit$design)", sep = '"')
  
  # contmatrix <- makeContrasts(eval(parse(text = cmd)), levels = fit$design)
  #eval(parse(text = cmd))
  fit2 <- contrasts.fit(fit, eval(parse(text = cmd)))
  fit2 <- eBayes(fit2)
  # print(summary(fit2))
  fold_change <- as.numeric(readline("ENTER Fc-values (log2) threshold required, recommended values '1' or '2' : "))
  comparison <- output1 <- topTable(fit2, coef=1, adjust.method = 'fdr',
                                    lfc= 0, number=nrow(filtered_data)) %>% as.data.frame()
  # store the result into csv
  
  #name1 <- readline("Enter file name to export fold change values (csv): ")
  name1 <- paste(option[1], "vs.", option[2], ".csv", sep = "")
  write_csv(comparison, file.path("data/", name1))
  # store significant result into csv
  output2 <- topTable(fit2, coef=1, adjust.method = 'fdr',
                      p.value =0.05, lfc=log2(1), number=nrow(filtered_data)) %>% as.data.frame()
  #name2 <- readline("Enter file name to export significant fold change values (p<0.05) (csv): ")
  
  name2 <- paste(option[1], "vs.", option[2], ".sig.csv", sep = "")
  write_csv(x=output2, file.path("data/Volc/", name2))
  # volcano plot
  # volcano input info.
  input <- comparison  # input of comparison group info which can be changed
  # significant data
  input$sig <- factor(input$adj.P.Val< 0.05 & abs(input$logFC) > fold_change)
  # size of significant data
  significant_lipids <- sum(input$adj.P.Val< 0.05 & abs(input$logFC) > fold_change)
  message("When the tests' q value treshold is 0.05 and the fold change threshold is ",fold_change,". The number of lipids which are statistically significant are: ", significant_lipids)
  
  
  # Define significant data for volcano plot graphing
  points <- input %>%
    rownames_to_column('lipid') %>%
    filter(abs(logFC) > fold_change & adj.P.Val< 0.05 ) %>%
    column_to_rownames('lipid')
  # subset the foldchange larger than 5 for text information                    
  fc_points <- input %>% 
    rownames_to_column('lipid') %>% 
    filter(abs(logFC) > 1 & adj.P.Val>= 0.05 ) %>% 
    column_to_rownames('lipid')
 
  
  # making volcano plot
  volc1 <- PlotVolc(input, points, fold_change) 
  volc1 <- volc1 + 
    geom_point(aes(col=sig, size=AveExpr)) + 
    geom_point(data = fc_points, colour = "black") +
    scale_size_continuous(range = c(0.25,2.5), breaks=c(10, 20, 30), guide = FALSE)+
    scale_color_manual(values=c("#bdbdbd", "#de2d26"), labels = c("Non-significant", "Significant")) +
    #geom_point(aes(fc_points, fill= "black")) +
    guides(color = guide_legend(title = "Fold Change")) +
    xlab(bquote("Fold change, " ~ log[2]~"("~textstyle(frac({.(option[1])}, {.(option[2])}))~")")) 
  print(volc1)
  name3 <- paste(option[1], "vs.", option[2], ".png", sep = "")
  # plot_name <- readline("Please input the volcano plot name: ")
  ggsave(filename = name3, path = 'plot/Volc/', device = "png")
  
  
  # build mutated data frame
  class_names <- rownames(input) %>% str_extract_all(., "(.+)\\(") %>% str_remove_all(., "\\(")
  #Add a new column to specify lipids
  input_lipids <- input
  input_lipids$lipclass <- class_names
  input_lipids$lipclass <- ifelse(input_lipids$sig=="TRUE", 
                                  case_when( input_lipids$lipclass %in% c("PC", "PE", "PG", "PI", "PS", "LPC", "LPE",
                                                                          "LPI", "dMePE", "CL")~ "Glycerophospholipids",
                                             input_lipids$lipclass %in% c("TG", "DG", "MG") ~ "Neutral lipids",
                                             input_lipids$lipclass %in% c("SM", "So", "SoG1", "Cer", "CerG1","CerG2",
                                                                          "GM2","GM3", "CerG2GNAc1", "CerG3GNAc1", "CerG3NAc2") ~ "Sphingolipids",
                                             input_lipids$lipclass %in% c("ChE","Cholestoral") ~ "Sterols",
                                             TRUE ~ "Other lipids"), 
                                  "n.s")
  # making levels for the class category
  class_levels <- c("Glycerophospholipids", "Neutral lipids", "Sphingolipids", "Sterols", "Other lipids", "n.s")
  # make the class category as factors
  input_lipids$lipclass <- factor(input_lipids$lipclass, levels=class_levels)
  # rearrange the data by its levels 
  input_lipids <-  input_lipids %>% 
    rownames_to_column('lipid') %>% 
    arrange(lipclass) %>% 
    column_to_rownames('lipid')
  # subset the significant points for text information
  significant_points <- input_lipids %>% 
    rownames_to_column('lipid') %>% 
    filter(!lipclass=="n.s") %>% 
    column_to_rownames('lipid')
  
  
  # plot the color version 
  volc2 <- PlotVolc(input_lipids, significant_points, fold_change)
  volc2 <- volc2 + 
    #   geom_text_repel(data=significant_points, aes(label=rownames(significant_point)), size=3) + 
    geom_point(aes(col=lipclass, size=AveExpr))  +
    scale_size_continuous(range = c(0.25,2.5), breaks=c(10, 20, 30), guide = FALSE) +
    scale_color_manual(values=c( "darkorange3","dodgerblue3", 'chartreuse4', 
                                 '#c51b8a', '#756bb1',"#bdbdbd"), drop=FALSE) +
    guides(color = guide_legend(title = "Lipid Class")) +
    #xlab(bquote("Fold change, " ~ log[2]~"("~expression(frac({.(option[1])}, {.(option[2])}))~scriptstyle(, AUC)")")) 
    xlab(bquote("Fold change, " ~ log[2]~"("~textstyle(frac({.(option[1])}, {.(option[2])}))~")")) 
  print(volc2)
  #plot_name <- readline("Please input the volcano plot name: ")
  name4 <- paste(option[1], "vs.", option[2], ".color.png", sep = "")
  ggsave(filename = name4, path = 'plot/Volc/', device = "png") #  width=15, height=15, dpi=300
  
  lipid_class <- rownames(input_lipids) %>% 
    str_remove_all(., "\\(.*\\)") %>%
    unique()
  lipids <- lipid_class %>% 
    paste0(., sep=", ", collapse = "") %>% 
    substr(., 1, nchar(.)-2) 
  
  message("\nPlease input the lipid name from the list below for displaying.\n", lipids, 
          "\nPlease note that the input is caps sensitive!")
  
  # customized_class <- readline("Input the name of lipid class(es), e.g. Cer TG: ") %>% 
  #   str_split(., "\\s+") %>% 
  #   unlist()
  # 
    
    customized_class <- check_lipid(lipid_class)

  
  if(length(customized_class)>1){
    customized_class <- paste(customized_class, collapse = "\\(|")
  }
  
  customized_points <- input_lipids %>% 
    rownames_to_column('lipid') %>% 
    filter(grepl(customized_class, lipid) & !lipclass=="n.s") %>% 
    column_to_rownames('lipid')
  
  title1 <- paste0("Volcanot plot for showing lipid class(es): ", customized_class)
  volc3 <- PlotVolc(input_lipids, customized_points, fold_change)
  volc3 <- volc3 +
    geom_point(aes(col=sig, size=AveExpr)) + 
    geom_point(data = fc_points, colour = "black") +
    scale_size_continuous(range = c(0.25,2.5), breaks=c(10, 20, 30), guide = FALSE)+
    scale_color_manual(values=c("#bdbdbd", "#de2d26"), labels = c("Non-significant", "Significant")) +
    #geom_point(aes(fc_points, fill= "black")) +
    guides(color = guide_legend(title = "Fold Change")) +
    xlab(bquote("Fold change, " ~ log[2]~"("~textstyle(frac({.(option[1])}, {.(option[2])}))~")")) +
    labs(title = title1)
  print(volc3)
  name5 <- paste(option[1], "vs.", option[2], ".customized.png", sep = "")
  # plot_name <- readline("Please input the volcano plot name: ")
  ggsave(filename = name5, path = 'plot/Volc/', device = "png")
  
  
  
  ether <- input_lipids 
  ether <- ether %>% 
    mutate(LipidMolec = rownames(ether), 
           ether = ifelse(str_detect(rownames(.), "\\(.*[ep]"), "YES", "NO"),
           Class = str_remove_all(rownames(.), "\\(.*")) 
  rownames(ether) <- ether$LipidMolec
  
  ether_class <-  ether %>% filter(sig == TRUE & ether == "YES") %>% select(Class) %>% unlist() %>% unique()
  extra_color_size <- length(ether_class)
  # Define significant data for volcano plot graphing
  ether_points <- ether %>%
    rownames_to_column('lipid') %>%
    filter(sig == TRUE & ether == "YES" ) %>%
    column_to_rownames('lipid')
  ether_points <- ether_points %>% 
    # filter(Class %in% c( "PE", "PC")) %>% 
    mutate(
      text = ifelse(sig== TRUE & ether == "YES", paste0("Ether Lipid: ", LipidMolec , "\nAdj.Pvalue: ", formatC(adj.P.Val, format = "e", digits = 2), "\n"), NA))
  #text = ifelse(adj.P.Val < 0.05 & ether == "YES", paste0("Ether Lipid: ", LipidMolec , "\nAdj.Pvalue: ", adj.P.Val, "\n"), NA))
  rownames(ether_points) <- ether_points$LipidMolec
  # excluded_ethers <- ether_class[!ether_class %in% c("PE", "PC")]
  #other_class <- ether %>% filter(!Class %in% c("PE", "PC" )) %>% select(Class) %>% unlist() %>% unique()
  #input$Class <- factor(ether$Class, levels = c("DG", "LPC", "PC", "PE", "PI", "PS", "TG", other_class))
  #ether$Class <- factor(ether$Class, levels = c( "PE", "PC", other_class))
  #ether_points$Class <- factor(ether$Class, levels = c("PE", "PC"))
  title2 <- paste0("Volcano plot for showing significant point from Ether lipid")
  volc4 <- PlotVolc(ether, ether_points, fold_change)
  volc4 <- volc4 +
    #geom_text_repel(data=ether_points, aes(label=rownames(ether_points)), size=3) + 
    geom_point(aes(col=lipclass, size=AveExpr))  +
    scale_size_continuous(range = c(0.25,2.5), breaks=c(10, 20, 30), guide = FALSE) +
    scale_color_manual(values=c( "darkorange3","dodgerblue3", 'chartreuse4', 
                                 '#c51b8a', '#756bb1',"#bdbdbd"), drop=FALSE) +
    guides(color = guide_legend(title = "Lipid Class")) +
    #xlab(bquote("Fold change, " ~ log[2]~"("~expression(frac({.(option[1])}, {.(option[2])}))~scriptstyle(, AUC)")")) 
    xlab(bquote("Fold change, " ~ log[2]~"("~textstyle(frac({.(option[1])}, {.(option[2])}))~")")) +
    labs(title = title2)
  print(volc4)
  
  name6 <- paste(option[1], "vs.", option[2], ".ether.png", sep = "")
  # plot_name <- readline("Please input the volcano plot name: ")
  ggsave(filename = name6, path = 'plot/Volc/', device = "png")
  
  
}






########################################################################
# function name: PlotVolc
# parameter: input, points, fold_change
# utility: plot volcano graph
########################################################################
PlotVolc <- function(input, points, fold_change){
  plot <- ggplot(input, aes(logFC, -log10(adj.P.Val))) + 
    geom_point(alpha = 0.8) +
    geom_vline(xintercept = c(-fold_change, fold_change), linetype="dashed") +
    geom_hline(yintercept= -log10(0.05), linetype="dashed", colour="black", size=0.2) +
    ylab("-Log10(q value)") +
    ggtitle("volcano plot") +
    theme_bw()+ 
    set_theme()+
    # theme(line=element_blank(),
    #       axis.line = element_line(colour = "black", size = 1),
    #       # panel.border = element_blank(),
    #       # panel.background = element_blank(),
    #       # legend.title = element_blank(),
    # ) +  
    theme(#line=element_blank(),
      axis.line = element_line(colour = "black", size = 1),
      legend.position = "right",
      axis.ticks.length = unit(1.5, "mm"),
      # width of tick marks in mm
      axis.ticks = element_line(size = 1.2),
      axis.text = element_text(size = 16, face = "bold", colour = "black"),
      axis.title = element_text(size = 18),
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 16)
      # panel.border = element_blank(),
      # panel.background = element_blank(),
      # legend.title = element_blank(),
    ) +
    geom_label_repel(data = points, 
                    aes(label = rownames(points)), 
                    size = 3,
                    direction    = "both",
                    hjust        = 0,
                    nudge_x = 0.2,
                    segment.size = 0.2,
                    box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.3, "lines"),
    ) +
    scale_y_continuous(breaks = c(0, -log10(0.05), 2, 4, 6, 8, 10), 
                       labels = c(0, "-log10(q value)", 2, 4, 6, 8, 10), 
                       expand = c(0, 0, 0.2, 0)) #+
  return(plot)
}


########################################################################
# function name: build_igraph
# parameter: data, group
# utility: build and plot pathway 
########################################################################
build_igraph <- function(data, group){
  info <- data  %>% filter(Groups==!!group) 
  # build node
  node <- data.frame(id = info$Class, label=info$Class, 
                     value = info$fc, shape = "dot" , 
                     saturation = info$saturation, color_saturation = info$color_saturation)  %>% 
    arrange(color_saturation)
  
  # build edges
  link <- data.frame(from=c("fatty acids", "LCBs", "Cer", "fatty acids", "G3P", "LPA", "PA", 
                            "PA", "DAG", "DAG", "DAG","PC", "PE", "LPE","PS","CDP-DAG", "CDP-DAG",
                            "CDP-DAG","PI", "PG", "PS", "LPC","LPS", "LPG", "LPI", "PC"),
                     to=c("LCBs", "Cer", "SM", "LPA","LPA", "PA", "DAG", "CDP-DAG", "TAG", 
                          "PC", "PE", "LPC", "LPE", "PE","PE", "PI", "PG", "PS", "LPI", "LPG", 
                          "LPS", "PC", "PS", "PG", "PI", "PE"),
                     label = c(rep("", 4), "GPAT", "AGPAT", "Lipin", "", "DGAT", rep("", 17)),
                     font.color = "green") %>% 
    mutate(type = "mention")
  
  link <- link %>% arrange(from)
  visNetwork(node, link, height="600px", width="100%", main="Network!")
  
  node$shadow <- TRUE # Nodes will drop shadow
  node$size = node$value # set fold change of group mean as node size
  node$color.background <- brewer.pal(nrow(data), "Blues")[node$color_saturation] # set node color
  node$color.border <- "black" 
  graph <- visNetwork(node, link, main = addquotes(!!group, " pathway"), layout = "layout_nicely") %>% 
    visConfigure(enabled = TRUE) %>%
    visLayout(randomSeed = 25)
  return(graph)
  
}






cal_sample_saturation <- function(data_fa, group_info){
  data_name <- deparse(substitute(data_fa))
  sample_raw_list <- group_info$samples
  # make a data store the pattern and count information
  lipid_list <- sapply(data_fa$LipidMolec, function(x) count_pathway(x)) 
  lipid_list <- t(lipid_list) %>% data.frame(., stringsAsFactors = FALSE)
  # name the pattern data
  colnames(lipid_list) <- c("FA_types", "SFA", "MUFA", "PUFA")
  
  # reformat the FA_types style by deleting extra speration mark "/"
  message("\nLipid molecules are classified as SFA, MUFA, PUFA or None, e.g Q10")
  lipid_list$FA_types <- apply(lipid_list, 1, function(x) x[1] <- x[1] %>% 
                                 str_replace(., "//", "/") %>% 
                                 str_remove(., "^/") %>% 
                                 str_remove(., "/$"))
  
  # assign the other pattern as None, e.g. Q10
  lipid_list$FA_types[lipid_list$FA_types==""] <- "None"
  
  # reorder the column position and put the counting patterns in between FA and FA Group key columns
  count_lipid <- data_fa %>% select(1:4) %>% cbind(., lipid_list)
  count_lipid <- count_lipid %>% cbind(., data_fa[5:ncol(data_fa)]) 
  rownames(count_lipid) <- NULL
  
  
  
  # write the count of pattern into count_lipid.csv file. 
  message("\nClassification of SFA, MUFA and PUFA are stored under count_lipid.csv and aggregated.csv")
  file_name1 <- paste("data/Saturation/count_lipid_", data_name, ".csv", sep = "")
  write.csv(count_lipid, file_name1)
  
  
  # delete lipid which can't do saturation analysis
  dis_lipid <- count_lipid %>% filter(FA_types=="None") %>% select(LipidMolec) %>% unlist()
  if(length(dis_lipid)>=1){
    lip <- paste(dis_lipid, sep = ", ", collapse = " ")
    message("\n\n", lip, " can't be used for saturation analysis and will delete in saturation analysis.\n\n") 
  }
  
  count_lipid <- count_lipid %>% filter(!LipidMolec %in% dis_lipid)
  
  
  # select columns of information for calculating the aggregation of MainArea of samples
  selected_lipids <- count_lipid %>% select(LipidMolec, Class, FA_types, SFA, MUFA, PUFA, contains("MainArea"))
  
  # transform attributes of SFA, MUFA, PUFA for calculations
  names <- c("SFA", "MUFA", "PUFA")
  
  # transformed_lipid <- selected_lipids %>% mutate_at(names, transform_to_numeric)
  transformed_lipid <-  selected_lipids %>% mutate_at(names, as.numeric) 
  
  # # filter all negative or 0 values
  # deleted_lipid <- transformed_lipid %>% filter_at(sample_raw_list, all_vars(.<=0))
  # # reserve last lipids
  # filtered_lipid <- anti_join(transformed_lipid, deleted_lipid)
  # delete lipid molecules contains negative values
  # reserved_lipid <- filtered_lipid %>% filter_at(sample_raw_list, all_vars(.>0))
  # if(nrow(filtered_lipid) > nrow(reserved_lipid)){
  #   deleted_molec <- anti_join(transformed_lipid, reserved_lipid) %>% 
  #     select(LipidMolec) %>% 
  #     unlist() %>% 
  #     paste0(., ", ", collapse = " ")
  #   message("Lipid Molecule(s): ", deleted_molec, "are deleted for saturantion and later ether composition analysis.\n\n")
  # }
  # 
  # reserved_lipid <- reserved_lipid %>% select(-LipidMolec)
  # selected_lipids <- selected_lipids
  
  # transform negative into NA
  reserved_lipid <- transformed_lipid %>% 
    mutate_at(sample_raw_list, list(~ifelse(.<0, NA, .))) 
  
  #calculate the aggregation of lipid by same class and FA_types
  aggregate_lipids <-  reserved_lipid %>% select(-LipidMolec) %>% group_by(Class, FA_types) %>% summarise_all(list(~sum(., na.rm = TRUE)))
  observation_count <-  reserved_lipid %>% group_by(Class, FA_types) %>% tally()
  aggregate_lipids <- left_join(aggregate_lipids, observation_count, by = c("Class", "FA_types")) 
  aggregate_lipids <- aggregate_lipids %>% select(Class, n, FA_types, SFA, MUFA, PUFA, all_of(sample_raw_list))
  
  # write the aggregation information into aggregated.csv file
  file_name2 <- paste("data/Saturation/aggregated_", data_name, ".csv", sep = "")
  aggregate_lipids %>% arrange(Class, FA_types) %>% write.csv(., file_name2)
  
  # select three variables Class, n, FA_types for analysis
  group_lipids <- aggregate_lipids %>% select(Class, n, FA_types)
  
  
  # calculate the SFA, MUFA and PUFA's percentages 
  fa_percent <- group_lipids %>% rowwise() %>% 
    # calculate how many SFA, MUFA and PUFA by row (by different SFA, MUFA and PUFA combination patterns)
    mutate(SFA = str_count(FA_types, "SFA"),
           MUFA= str_count(FA_types, "MUFA"), 
           PUFA= str_count(FA_types, "PUFA"), ) %>%   
    # calcutate SFA, MUFA and PUFA's percentage by row
    mutate("%SFA" = SFA/(SFA+ MUFA + PUFA), 
           "%MUFA" = MUFA/(SFA+ MUFA + PUFA), 
           "%PUFA"= PUFA/(SFA+ MUFA + PUFA)) 
  
  # combine the percentage information with aggregated lipid info
  fa_percent <- aggregate_lipids %>% ungroup() %>% select(all_of(sample_raw_list)) %>% cbind(fa_percent, .)
  file_name3 <- paste("data/Saturation/fa_percentage_", data_name, ".csv", sep="")
  write_csv(fa_percent, file_name3)
  
  # find median and mean for each sample in each group and calculate its value
  data_list <- calc_group(fa_percent, group_info)
  dt <- data_list[[1]]
  dd <- data_list[[2]]
  data <- data_list[[3]]
  
  dt <- dt %>% mutate_if(is.numeric, list(~formatC(., format = "f", digits = 3)))
  
  
  file_name4 <- paste("data/Saturation/individual_saturations_", data_name, ".csv", sep = "")
  file_name5 <- paste("data/Saturation/class_saturations_", data_name, ".csv", sep = "")
  file_name6 <- paste("data/Saturation/mean_median_", data_name, ".csv", sep = "")
  write_csv(dt, file_name4)
  write_csv(dd, file_name5)
  write.csv(data, file_name6)
  return(data)
}




########################################################################
# function name: impute_not
# parameter: condition, data, sample_list
# utility: impute 0 or negative value for log transformed data
########################################################################
impute_not <- function(condition, data, sample_list){
  if(condition == "y"){
    message("\n\nWarnings!!!!!!\nYou are now using imputated data for analysis.\n\n")
    log2_lipids <- data %>% mutate_at(sample_list, log2trans) 
    replace_inf_lipids <- log2_lipids %>% mutate_at(sample_list, ReplaceInf)
    preImputated_lipids <- replace_inf_lipids %>% select(sample_list) 
    imputated_lipids <- ImputeMinProb(preImputated_lipids, 0.01, 1)
    imputated_lipids$LipidMolec <- log2_lipids$LipidMolec
    imputated_lipids$Class <- log2_lipids$Class
    imputated_lipids <- imputated_lipids %>% select(LipidMolec, Class, sample_list)
    write_csv(imputated_lipids, "data/Volc/imputeMolec.csv")
    message("\nLog 2 transformed data are stored under log.molec.csv")
    message("Imputed data are stored under imputeMolec.csv")
    # impute same value in raw data
    log2_filtered_lipids <- data
    log2_filtered_lipids <- log2_filtered_lipids %>% arrange(LipidMolec)
    imputated_lipids <- imputated_lipids %>% arrange(LipidMolec)
    log2_filtered_lipids[colnames(log2_filtered_lipids) %in% c( "LipidMolec", "Class", sample_raw_list)] <- imputated_lipids
    # store log2 transformed raw data into log2_filtered_data.csv
    write_csv(log2_filtered_lipids, "data/Volc/log2_filtered_data.csv")
    return(list(condition, imputated_lipids))
  }else if(condition == "n"){
    # return data with selected columns
    data2 <- data %>% select(Class, LipidMolec, sample_list)
    write_csv(data2, "data/Volc/non_impute_molecs.csv")
    return(list(condition, data2))
  } else{
    condition <- readline("Typed wrong. Please Type Y/N for impute or not: ") %>% str_to_lower()
    return(impute_not(condition, data, sample_list))
  }
}


filter_invalid <- function(data, group_info, invalid_data){
  if(nrow(invalid_data) != 0){
    sample_list <- group_info$samples
    # data negative and empty value percent information 
    neg_percent_info <- detect_invalid(invalid_data, group_info) 
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
    option <- retype_choice("1/2")
    filtered_negs <- negs_all %>% filter_at(sample_list, any_vars(!is.na(.)))
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
      ref_data <- read_csv("data/invalid.csv", col_types = cols())
      # deleting the invalid lipid molecules in raw data based on your standards and saved it in pre_filtered.lipidomics.csv
      data <- fix_invalid_by_choice(data, manual_data, ref_data)
      write_csv(data, "data/manual_filtered.lipidomics.csv")
    }else{
      message("The pipeline will first transform all the negative value into NA.")
      message("If negative percentage is over 50% in a group, all the values in the group for the molecule will be transformed into NA.")
      message("If a molecule which negative percentage is over 50% for all groups, it will then be deleted.")
      # delete the molecule which negative values of replicates for all group are over 50% (all NA.)
      deleted_neg_molec <- negs_all %>% filter_at(sample_list, all_vars(is.na(.))) 
      deleted_molec <- negs_all %>% filter_at(sample_list, all_vars(.==0)) %>% bind_rows(., deleted_neg_molec)
      if(nrow(deleted_molec) != 0){
        deleted <- deleted_molec$LipidMolec %>% unlist() %>% sort() %>% paste0(., ", ", collapse = "") %>% substr(., 1, nchar(.)-2) 
        message("\n
  Since the lipid molecule ", deleted, " is invalid (all negative or all 0 values) after background subtraction. 
  \nIt will be deleted in the filtered data.")
        # delete corresponding lipid molecules in total data
        data <- anti_join(data, deleted_molec, by = "LipidMolec") #%>% select(LipidMolec, contains("MainArea"))
      }
      write_csv(data, "data/auto_filtered.lipidomics.csv")
    }
  return(data)
  }
}


# 
# filter_invalid <- function(data, data_copy, group_info, invalid_data){
#   if(nrow(invalid_data) != 0){
#     sample_list <- group_info$samples
#     # data negative and empty value percent information 
#     neg_percent_info <- detect_invalid(invalid_data, group_info) 
#     neg_percent <- neg_percent_info[[1]]
#     # transform all negative value with NA
#     neg_info <- neg_percent_info[[2]] %>% 
#       select(-contains("APValue"), -contains("ARatio"),
#              -contains("Grade"), -c(A, B, C, D, FA))
#     # replace all values for a group into NA if the negative percentage is over 50%
#     negs_all <- neg_percent_info[[3]] %>% 
#       select(-contains("APValue"), -contains("ARatio"),
#              -contains("Grade"), -c(A, B, C, D, FA))
#     # the invalid data after background subtraction.
#     message("\n
#   For lipid molecules that contain zero values or negative values (background subtracted), 
#   These values are subsequently replaced as non-valid values (NA). 
#   Fold change analyses is performed using only samples containing valid values")
#     # write negative and 0 percentage inforamtion
#     write_csv(neg_percent, "data/neg.percent.csv")
#     # write potential invalid lipid molecules information
#     write_csv(neg_info, "data/checkInvalid.csv")
#     # write copy of potential lipid molecules information
#     write_csv(neg_info, "data/invalid.csv")     ## copy data for invalid lipids information
#     message("Please view file imputeNA.csv for all the data contains negative values after background subtraction.")
#     # write data which transform negative into NA
#     write_csv(negs_all, "data/imputeNA.csv")
#     message("\nType 1 if you would like the pipleline to proceed with this function \nType 2 if you prefer to exlcude certain lipid molecules for fold change analysis ")
#     option <- readline("Please type 1/2: ") %>% str_to_upper()
#     filtered_negs <- negs_all %>% filter_at(sample_list, any_vars(!is.na(.)))
#     if(option == "2"){
#       # this step need manually editting the invalid lipid molecules on your computer for advanced users
#       message("Select 'checkInvalid.csv' to manually exclude specific lipid molecules and click SAVE.")
#       # here stops 10 seconds
#       Sys.sleep(10)
#       # deleting the lipid molecules you select in the checkInvalid.csv
#       #message("Now we need to open the changed file checkInvalid.csv after you deleting the invalid lipid molecules.")
#       continues <- readline("If you finished preprocess the data, please continue and press Y: ")
#       # advanced users
#       manual_data <- read_csv("data/checkInvalid.csv", col_types = cols())
#       # deleting the invalid lipid molecules in raw data based on your standards and saved it in pre_filtered.lipidomics.csv
#       data <- fix_invalid_by_choice(data, manual_data, filtered_negs)
#       write_csv(data, "data/manual_filtered.lipidomics.csv")
#       
#     }
#     while(option != "2"){
#       if(option ==  "1"){
#         message("The pipeline will first transform all the negative value into NA.")
#         message("If negative percentage is over 50% in a group, all the values in the group for the molecule will be transformed into NA.")
#         message("If a molecule which negative percentage is over 50% for all groups, it will then be deleted.")
#         # delete the molecule which negative values of replicates for all group are over 50% (all NA.)
#         deleted_neg_molec <- negs_all %>% filter_at(sample_list, all_vars(is.na(.))) 
#         deleted_molec <- negs_all %>% filter_at(sample_list, all_vars(.==0)) %>% bind_rows(., deleted_neg_molec)
#         if(nrow(deleted_molec) != 0){
#           deleted <- deleted_molec$LipidMolec %>% unlist() %>% paste0(., ", ", collapse = "") %>% substr(., 1, nchar(.)-2)
#           message("\n
#   Since the lipid molecule ", deleted, " is invalid (all negative or all 0 values) after background subtraction. 
#   \nIt will be deleted in the filtered data.")
#           # delete corresponding lipid molecules in total data
#           data <- anti_join(data, deleted_molec, by = "LipidMolec") #%>% select(LipidMolec, contains("MainArea"))
#           data_copy <- anti_join(data_copy, deleted_molec, by = "LipidMolec")
#         }
#         write_csv(data, "data/filtered_to_NA.csv")
#         write_csv(data_copy, "data/post_filtered.lipids.csv")
#         break
#       } else{
#         option = readline("You typed wrong, please type again, 1/2: ")
#       }
#     }
#   }
#   
#   return(data)
#   
# }








