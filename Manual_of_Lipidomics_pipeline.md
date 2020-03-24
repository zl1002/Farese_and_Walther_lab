Manual of Lipidomics Pipeline
================
Wenting
2/21/2020
3/34/2020 ZWL edit





### Motivation

This R-script is developed by the joint laboratory of Drs. Robert Farese, Jr. 
and Tobias Walther at Department of Molecular Metabolism,
Harvard T.H. Chan School of Pulic Health. This analysis pipeline facilitates
intepretation of mass-spretrometry derived data for lipidomics experiments
Lipid identification and quantification was performed using LipidSearch
software (Thermo Fischer Scientific). Key features of this R-script pipeline
includes quality control, statistical and normalization analyses, as well as
data visualization

### Installation and Set up

  - Requires R version 3.6.0 or later.  
  - For Mac users, XQuartz must to be installed manually.  
  - Any additional packages/libraries will be installed and uploaded
    automatically.
  - Pipeline will automatically generate store all the output under invidually generated
    directories, converted data, and plots.

<!-- ![an image caption Source: screenshot.](button_source.png) -->

### Additional Information

![](Manual_of_Lipidomics_pipeline_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

### Step-by-step Instructions

Please note that pipeline will generate 2 folders: `plot` and `data`.
And within each folder, there will be additional subfolders, including
`QC`, `Quantification`, `Saturation`, `Length`, `Ether`, `Volc`.  
*Please note that for plot display, the data are collected from
different experiments.*  
Data stored within the folder `data` can be used to generate different
plots.  
Users commands input are case insensitive except lipid classes and
group names.  


##### File Conversion (txt -> csv):

  - Objective: Convert txt file to csv file for subsequent analysis 
  - Utility: To run the pipeline,
      - 1. Download all the scripts in same directory with R
        project.  
      - 2. Run script **txt\_to\_csv\_version.R** for converting raw txt file by
        typing *`source`* button, or type command
        *`source("txt_to_csv_version.R")`* in the console. Please note
        that, the converted csv files are stored in directory
        `converted`.  
      - 4. Before running the main script, make sure that XQuartz is
        installed for mac user.
      - 3. Run main script **FWL\_lipidomics\_version.R**. Press the *`source`*
        button, or type *`source("FWL_lipidomics_VERSION.R")`* in console.  

<!-- end list -->

``` r
# commands
source("txt_to_csv_version.R")

source("FWL_lipidomics_version.R")
```

    Are you using PC or MAC?
    Please type PC/MAC: MAC

``` 


Please make sure you installed Quartz before running pipeline.
          
Do you want to continue?
Please type Y/N: y
```

  - Useful short cuts: `cmd (ctr) + Enter` (run selected command), `cmd
    (ctr) + shift + c` (comment/uncomment).

##### Preprocessing and QC:

**(1) Read data**

  - Objective: Console will list the files in the converted directory
    and it will pop out a command line to ask user to input the index of
    the file for analysis.  
    <!-- - Code display: -->

  - Utility: Please enter the number of file index. After input the file
    index, console will diplay the file you choose.  

  - Example display:

<!-- end list -->

    The following files had been generated. 
    Select ONE for subsequent the list of file names:
    1 030420_cGSL_mouse-brains_publication.csv
    2 10252019_PGRN_ctl-het-ko.csv
    3 190911_IDX_WH_pB_Job1702_pos_MC_01.raw.csv
    4 20171127_AH_MEFpilot.csv
    5 20171127_JC_Seipin.csv
    6 20180823_LB_FIT2.csv
    7 20190130_Aditi_DGATGPAT.csv
    8 20191220_CJ_s2.csv
    9 20200210_JSCJ_s2.csv
    10 20200317_CJ_s2.csv
    11 20200317_CJ.csv
    12 Basti's work.csv
    13 BRAIIIN_single.csv
    14 cells_single.csv
    15 col_test.csv
    16 FGL_Ld_Lipidomics_0525.csv
    17 RNF213_LB_07292019.csv
    
    Please input the index of the file: 8

    [1] "[1] converted/20171127_JC_Seipin.csv"

**(2) Confirmation of Experimental Samples**

  - Objective: This part will check if any samples need to be
    excluded from analysis.  
    <!-- - Code Display: -->

  - Utility: Please enter Y or N (case insensitive).  

  - Example diplay:
    
      - When typing N, no samples are deleted for analysis.  
      - When typing Y, console will pop out message below, user then
        input the sample name. Please note that, the sample name has
        fixed format, sample 1 would be s1 or S1. And the samples are
        separated by space.  
      - If your data are merged from multiple sets of experiments, it is strongly advised to delete samples from different
        experiments from your data. Otherwise data processing, statical analysis and normalization will be affected.

<!-- end list -->

``` 

Do you have experimental controls like internal standards or extraction control?

Please type Y/N: N
```

``` 

Do you have samples excluded for analysis which including experiment controls (internal standards)?


Please type Y/N: Y

Indicate which samples used as controls or will be deleted for analysis
Sample ID , eg. s22 s23

option standards -----> s17 s18
```

##### **Quality Check**

Related data and plots will be under `QC` subfolders of `plot` and
`data`.

**(1) Set Filtering Parameters**

  - Objective: Data filtering based on user’s preference. Set flexible
    parameters k and j, depending on the experimental set up. It sets the minimum number
    of Grade A and B for lipid idenfication and required p-value.  
    Please note that Filter standard Rej = 0 will be automatically be
    applied. <!-- - Code display:  -->

  - Utility: Please input two corresponding numbers.  

  - Example display:

<!-- end list -->

``` 

Data are filtered using 3 criteria.  

 1. Not rejected by LipidSearch (Rej n=0); 
        
 2. minimum number of Grade A+B required; 
        
 3. minimum number significant identification (p-value p<=0.001) for LipidSearch standard.
        
 Filtered data is stored in the filtered.raw.data.csv.

Minimum number of identified molecules (Grade A+B) required in all samples, n>=3:

Minumum number of significantly identified lipids p<=0.001 in all samples, n>=3: 
```

**(2) Check Background Noise from Solvent-only Run**

  - Objective: Plot the back ground information (sample blank)
    abundance in each detected lipid class. The plot will be saved as
    `background.png` in `plot` directory.

  - Output
example:

<img src="display/background.png" width="40%" style="display: block; margin: auto;" />

**(3) Lipid Class Summary**

  - Objective: Generate tables and plots for showing detected
    lipid class summary information. The plot and data will be saved as
    `prop_summary.png` and `proportion_classes.csv`.  
  - Output
example:

<center>

<img align="middle" src="display/prop_summary.png"  width="50%" height="50%"/>

</center>

**(4) Rention Time Analysis**

  - Objective: Diplay retention time (min) vs. log transformed AUC (area under the curve) 
  value of all samples for each lipid class. And the plot is saved as `all_retention.png`.  
  - Output
example:

<img src="display/all.retention.png" width="100%" style="display: block; margin: auto;" />

**(5) Marking Odd Chains and External Standard e.g. TG(17:1/17:1/17:1) Abundance in
all samples**

  - Objective:
      - Console will display the number of lipid molecules containing
        odd chains and its percentage. The odd chain lipid molecules are
        stored in `odd_chains.csv`.  
      - And for external standards TG(17:171:17:1), if it is used, pipeline could
        detect its abundance in all samples and make a plot
        `TG17_all.png`. The standard could be replaced with other future
        standard.  
      - Please note that the bars are ordered by AUC value for each
        sample.
  - Example display:

<!-- end list -->

``` 

There are 157 lipid molecules contain odd chains. 
          
The odd chain of fatty acids percent is 24.04% in total.

The odd chain information is stored in odd_chains.csv.
```

<img src="display/TG17_all.png" width="40%" style="display: block; margin: auto;" />

**(6) Resolving Duplicates of Same Lipid Molecules**

  - Objective: Pipeline will dectect same lipid molecules with
    different retention time. And it will filter the duplicates based on
    2 criteria. User will need to input criteria A or B (case
    insensitive). Duplicated lipid molecules will be under
    `duplicated.molecules.csv` and `diff_RT.csv`, reserved lipid
    molecules will be under `reserved_duplicates.csv`, and filtered data
    will be stored under `rm_duplicates.csv`.
      - Criteria A: use only ONE lipid molecule with largest AUC (Main
        Area Under Curve).  
      - Criteria B: Sum AUC for all duplicates.  
        `Please note that, method B will produce non-valid values in some columns
        since the non-valid values cannot be summed.`
  - Example display ( when exist identical lipid molecules)

<!-- end list -->

``` 

!!!Attention: Identical lipid molecules with multiples retention time. Please note that the duplicate lipid molecules are stored in reserved_duplicates.csv 
 !!!!!! Potential sample contamination. 
 To PROCEED, pick one: 

Differences in retention time for identical lipid molecule are stored under diff_RT.csv

 A: Use only ONE lipid molecule with largest main area under curve, OR 
 B: Summation of main area under curve of ALL identical lipid molecule.

Filtered lipid molecules sans duplicates are stored under removeduplicates.csv
Enter 'A' or 'B': A
Filterted lipid molecules sans duplicates are stored under removeduplicates.csv
```

**(7) Enter Biological Groupping Information**

  - Objective: Enter the different group information (how the experiment was designed)
    i.e. experiment group names and its samples. Pipeline will extract sample AUC’s
    corresponding columns information. Group information will be stored
    in `group_information.csv`.
      - Please note that the group names will be asked for input in
        later analysis and it will be *CASE SENSITIVE*.  
      - Group number must be numerical, the sample naming format is
        consistent with previous example, e.g. s1 s2.  
      - Pipeline will ask if user want to edit group information for
        correcting group information later.  
      - User can check group information under `group_information.csv`
  - Code display: The variable `label` in the code will be suffix in the
    genrated plot names and could be modified. And sample information
    for subsequent analyses could be re-edit after user typing Y in
    command line.

<!-- end list -->

``` r
# pca and correlation plots
label <- "initial"
info_list <- PCA_pairs_Plot(sample_info, group_names, filtered_lipidomics2, label)
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
# variable store group names
group_names <- unique(group_repeats)
# variable store group number
ngroups <- length(group_names)
```

  - Example display

<!-- end list -->

``` 

Provide infomation of experimental groups
        
How many experimental groups: 4
        
Description for Group  1  (name): Control
        
Which samples assigned to Group  1 (sample number, e.g. s1 s2 s3 ): s1 s2 s3 s4
        
Description for Group  2  (name): KO
        
Which samples assigned to Group  2 (sample number, e.g. s1 s2 s3 ): s5 s7 s8
        
Description for Group  3  (name): OE1
        
Which samples assigned to Group  3 (sample number, e.g. s1 s2 s3 ): s9 s10 s11 s12
        
Description for Group  4  (name): OE2
        
Which samples assigned to Group  4 (sample number, e.g. s1 s2 s3 ): s13 s14 s15 s16
        
CONFIRM the group information below
        
List of 4
        
$ Control: chr "s1 s2 s3 s4"
        
$ KO     : chr "s5 s7 s8"
        
$ OE1    : chr "s9 s10 s11 s12"
        
$ OE2    : chr "s13 s14 s15 s16"

Do you want to edit group infomation? 

Y/N: n

Take a look at the sample info and its column position information in the file below
List of 8
 $ : chr [1:4] "MainArea[s1]" "MainArea[s2]" "MainArea[s3]" "MainArea[s4]"
 $ : int [1:4] 48 49 50 51
 $ : chr [1:4] "MainArea[s5]"  "MainArea[s7]" "MainArea[s8]"
 $ : int [1:4] 52 54 55
 $ : chr [1:3] "MainArea[s9]" "MainArea[s10]" "MainArea[s11]" "MainArea[s12]"
 $ : int [1:3] 56 57 58 59
 $ : chr [1:4] "MainArea[s13]" "MainArea[s14]" "MainArea[s15]" "MainArea[s16]"
 $ : int [1:4] 59 60 61 62
 - attr(*, "dim")= int [1:2] 2 4
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:2] "sample.names" "col.index"
  ..$ : chr [1:4] "Control" "KO" "OE1" "OE2"
```

<img src="display/pairs.plot.1.initial.png" width="40%" /><img src="display/pairs.plot.2.initial.png" width="40%" /><img src="display/pairs.plot.3.initial.png" width="40%" /><img src="display/pairs.plot.4.initial.png" width="40%" />
<img align="middle" src="display/sample.pca.initial.png"  width="50%" height="50%"/>

##### **Background Subtraction**

  - Objective: Provides the option for background subtraction (if preferred). 
    If user choose to remove background noise, the pipeline will subtract 
    individual AUC with background AUC for subsequent analyses. 
    Non-background substracted data will be stored as
    *filtered\_lipidomics\_copy* variable and `rm_duplicates.csv` in
    previous step. 
    Background substrated data will be saved under
    `subtracted_lipids.csv`. After background subtraction, pipeline will
    provide 2 optional methods for fixing potential invalid lipid
    molecules. 
    Method 1 will automatically delete lipid molecules which
    AUC are all negative or 0. Method 2 will ask user to delete
    potential invalid lipids in file `checkInvalid.csv` manually. Please
    note that `invalid.csv` is its copy and can be used as reference for
    user. Console then will pop out
  - Utility:
      - Background substrated data could contain 0 or
        negative values. Pipeline will detect potential non-valid value
        and calculate its sample size percentage information. Please
        check files `neg.percent.csv`, `checkInvalid.csv` or
        `invalid.csv`, `imputeNA.csv`.  
      - type 1 or 2  
  - Example display

<!-- end list -->

``` 

For lipid molecules that contain zero values or negative values (background subtracted), 
        
These values are subsequently replaced as non-valid values (NA). 
        
Fold change analyses is performed using only samples containing valid values
        
Please view file imputeNA.csv for all the data contains negative values after background subtraction.
        
Type 1 if you would like the pipleline to proceed with this function
        
Type 2 if you prefer to exlcude certain lipid molecules for fold change analysis
```

    - if choose method 1

``` 

Please type 1/2: 1
        
The pipeline will first transform all the negative value into NA.
        
If negative percentage is over 50% in a group, all the values in the group for the molecule will be transformed into NA.
        
If a molecule which negative percentage is over 50% for all groups, it will then be deleted.
```

    - if choose method 2

``` 

Please type 1/2: 2

Select 'checkInvalid.csv' to manually exclude specific lipid molecules and click SAVE.

If you finished preprocess the data, please continue and press Y: y
```

##### **Quantification Analysis**

Related data and plots will be under `Quantification` subfolders of
`plot` and `data`.

**(1) Quantification of Total Lipid Classes (mean and stdev)**

  - Objective: Summation of all lipid molecules under same lipid class
    (`aggregated_class.csv, total_class.csv`). Display its mean and
    standard deviation for each experimental group. Note:
    data will include negative values if background subtraction was performed.
  - Example
display:

<img src="display/total.class.png" width="40%" style="display: block; margin: auto;" />

**(2) Quantification Individual Lipid Molecules of Same Class (mean and stdev)**

  - Objective: Visualize mean and standard deviation in each
    experiment group for individual lipid class (`all_lipidmolec.csv`).
    And all the plots generated by `EachClassPlot` function
    are under folder `classes` in plot.  
  - Utility:
      - Check the code below, the `nbar` variable which could be
        modified is set for approximate max number of bars displayed in
        one plot. If you want display molecules of each lipid class in
        one page, please uncommented the part by choosing nbar as
        variable `lipidNO_max`.  
      - `post_name` variable is common suffix name for plot and also
        could be modified.  
      - Please note that, the function `EachClassPlot` for plotting
        could also be used for plotting other value like median of
        group. Correspondingly, the function `cal_lipid_statistics`
        could be used for calculating different demands, like mean, sd
        or other statistics. And please make sure these methods process
        the negative value properly, e.g. na.rm = TRUE.  
  - Code display:

<!-- end list -->

``` r
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
each_class <- left_join(lipid_mean_long, lipid_sd_long) 
write_csv(each_class, "data/Quantification/all_lipidmolec.csv")
message("\nQuantification analysis for individule lipid class")
par_eachclass <- c("LipidMolec", "mean", "Groups", "sd")
# maximum bar number limits
lipidmolecNO_max <- filtered_lipidomics %>% group_by(Class) %>% tally() %>% select(n) %>% unlist() %>% max()
# setting plot labs titles
labs1 <-  labs(x="Acyl composition", y="Main Area", 
               caption="Error bar is the standard deviation for each class in each group", fill = NULL)
# setting the plot limits when the bar numbers exceed the threshold nbar
nbar <- 70    # estimation of threshold which can be modified and at least bigger than group number
post_name <- ""
pars <- list(nbar, ngroups, par_eachclass, plot_all, post_name, labs1)
message("\nEach plot is split no more than ", nbar, " bars for display")
EachClassPlot(each_class, pars)

# # if uncommented part below, 
## overview of each class plot by its largest bar setting 
# nbar <- lipidmolecNO_max
# post_name <- "all"
# pars <- list(nbar, ngroups, par_eachclass, plot_all, post_name, labs1)
# message("\nAlternative display quantification of individule lipid class (all lipids in a class in the same png)")
# EachClassPlot(each_class, pars)


################################################################################ turn off Quartz for mac
dev.off()
options(device = "RStudioGD")                                                                                                                                 
###############################################################################
```

  - Example
display:

<img src="display/CL.1.png" width="40%" /><img src="display/Cer.png" width="40%" />

**(3) Visualization of Lipid Class, normalized by median**

  - Objective: Samples of each experiment group for lipid class are
    normalized (`raw_class_median.csv, normalized_class_median.csv,
    normalized_class_median_long.csv`) by median value of control group,
    and then median for each group is displayed as dot plot and box
    plot. If median value of lipid in control group is 0 or negative,
    the molecule will be removed in fold-change analysis. User will need
    to enter the control group name as comparison for fold change
    analysis.

<img src="display/class_median_dot.png" width="45%" /><img src="display/class_median_box.png" width="45%" />

**(4) Visualization of Lipid Molecules, normalized by mean/median
value**

  - Objective: Samples of lipid molecules are normalized by either mean or
    median value. And corresponding mean or median value will be plotted.
    This visualization also implies `EachClassPlot` function,
    thus the plots are under classes, data are stored as
    `molecules_group_statics.cs, raw_molec.csv, normalized_molec.csv,
    normalized_molec_mean( or normalized_molec_median)` under
    Quantification folder. And there are two types of violin plots for
    display normalized data by choosing from mean or median as well. One
    violin (`molec_violin.png/svg`) will display all lipid class in one
    plot. One will only display the interactive customized lipid
    class(es) violin plot (`molec_violin_all.html`) to save time.  
    Please note that the interactive plots will display in the `Viewer` console,
    click panel `Plots` to view other plots.  
  - Example display: The first plot is normalized by mean, while second
    by median. The third one is violin plot, and only static voilin plot
    is showed in the manual
book.

<img src="display/PE.3_mean_in_CONTROL_HET_1.png" width="30%" /><img src="display/DG.1_median_in_control_het_1.png" width="30%" /><img src="display/molec_violin.png" width="30%" />

##### **Fatty Acids Saturation Analysis**

  - Objective: This part will analyze saturation of fatty acids
    (`SFA`), Monounsaturated fatty acids (`MUFA`) and Polyunsaturated
    fatty acids (`PUFA`). All data from analysis will be stored under subfolder
    `Saturation`. And the data used for plotting are: `fa_mean.csv,
    fa_median.csv, fa_normalized_data.csv, fa_normalized_mean_long.csv`.
    Note: If control group contains mean value 0 for fold-change, the lipid 
    will be excluded from the analysis.Negative values are transformed into
    non-valid values after background subtraction and will also be excluded from the analysis.
  - Utility: User will need to enter the control group name for
    data normalization.  
  - Example display: The first stack-plot uses either mean or median value from
    each group. And the third plot uses the mean value of each group. The
    fourth plot uses data which is normalized by the mean value of
    control group. And final plot shows the percentage for fatty acids
    saturation.

<img src="display/meanBased.stackplots.png" width="45%" /><img src="display/medianBased.stackplots.png" width="45%" /><img src="display/meanBased.fattyAcids.png" width="45%" /><img src="display/meanBased.fc.png" width="45%" /><img src="display/percentage.png" width="45%" />

##### **Fatty Acids Length Analysis**

  - Objective: Fatty Acids Length analysis will analyze the
    abundance, percentage of Short-chain fatty acids (SCFA),
    Medium-chain fatty acids (MCFA), Long-chain fatty acids (LCFA) and
    Very long chain fatty acids (VLCFA) in each sample and group. All
    data and plots will under subfolder Length of folder data and
    plot.  
    There will be many plots which name suffix 'fc.png.' User will need to choose
    control group and normalization by either mean or median values. Please
    note that negative are transformed into non-valid values after background
    subtraction and will be excluded for analysis.
  - Example display: The first plot displays the different fatty acids
    length abundance in each sample. The second plot displays the
    percentage of different length in each sample. The third stack-plot
    displays different length abundance in each group. The fourth plot
    displays the percentage of different length in each group. The fifth
    plot display the different length abundance in each group. The final
    plot shows all the chains fold-change in PE lipid
class.

<img src="display/fa_length.png" width="45%" /><img src="display/fa_length_percentage.png" width="45%" /><img src="display/fa_length_group.png" width="45%" /><img src="display/fa_length_gr_percentage.png" width="45%" /><img src="display/fa_length_gr.png" width="45%" /><img src="display/PE.fc.png" width="45%" />

##### **Ether lipid analysis**

  - Objective: This part will analyze the abundance of ether lipids.
    The first 7 plots are ether lipids abundance analysis. The
    rest is combined with saturation analysis to identify PUFA percentage in
    ether lipids. And all the data and plots are under subfolder 'Ether'
    of data and plot.

  - Example
display:

<img src="display/ether.png" width="40%" /><img src="display/ether_percentage.png" width="40%" /><img src="display/ether_group.png" width="40%" /><img src="display/ether_group_percent.png" width="40%" /><img src="display/ether_abundance.png" width="40%" /><img src="display/ether_molec_abundance.png" width="40%" /><img src="display/normalized_ether_molec_abundance.png" width="40%" /><img src="display/pufa_ether.png" width="40%" /><img src="display/pufa_ether_percentage.png" width="40%" /><img src="display/pufa_ether_group.png" width="40%" />

##### **Data Imputation for Differential Expression Analysis**

**(1) Test Random Sample Distribution**

  - Objective: One random sample will be log transformed and
    visualized its distribution. The note generated is just a reference
    and will stored. The log transformed data will be under
    `log.molec.csv`.  
  - Example display:

<img src="display/sample_distr.png" width="40%" />

**(2) Volcano Plots**

  - Objective: Volcano plots will generate 4 formats of same plot.
    The first plot will only display significant lipid molecule. The
    second plot will classify significant lipid molecules into
    Glycerophospholipids, Neutral lipids, Sphingolipids, Sterols and
    other lipids. The third plot is customized style, and it will only
    display the lipid class user input from the lipid class list. The
    fourth plot will only mark significant ether lipids if exist. Data
    and plots will be stored under subfolder 'Volc' of data and plot.  
  - Utility: User will need to enter number of comparisons, group names for
    comparison, fold-change threshold and customized lipid class(es).
    Please note that the group names need to be consistent with previous
    and separate by `vs` (or `VS`).  
  - Example display:

<!-- end list -->

``` 


Warnings!!!!!!

You are now using imputated data for analysis.


[1] 0.8253833

Log 2 transformed data are stored under log.molec.csv

Imputed data are stored under imputeMolec.csv


How many volcano plots to generate: 1
        
Enter groups names for comparison, spaced by 'vs', e.g. KO vs WT: LINOLEATE VS OLEATE
        ENTER Fc-values (log2) threshold required, recommended values '1' or '2' : 1
        

Please input the lipid name from the list below for displaying.
PE, PC, CL, LPC, TG, DG, PA, PS, LPE, PI, dMePE, LPS, So, PG, LPI, DGDG, SM, MGMG, MGDG, ChE, Cer, CerG3, CerG1, Co

Please note that the input is caps sensitive!

Input the name of lipid class(es), e.g. Cer TG: TG
```

<img src="display/LINOLEATEvs.OLEATE.png" width="40%" /><img src="display/LINOLEATEvs.OLEATE.color.png" width="40%" /><img src="display/LINOLEATEvs.OLEATE.customized.png" width="40%" /><img src="display/LINOLEATEvs.OLEATE.ether.png" width="40%" />
