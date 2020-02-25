Manual
================
Wenting
2/21/2020

### Motivation

The laboratory of Drs. Robert Farese, Jr. and Tobias Walther at the
Harvard T.H. Chan School of Pulic Health has been using LipidSearch
software (version 4.1) for some time, currently we build a pipeline
based on Mass Spectrometry data processed by this software. We filtered
the data with flexible standards, did QC, statistics analysis and
visualization with extracted information from original txt file.

### Set up

Required R version 3.6.0 or later.

For Mac users, XQaurtz need to be installed manually.

Any additional packages/libraries will be installed and uploaded
automatically.

Pipeline will store all the output under 4 automatically generated
directories, converted, data, plot, plot/classes.

To run the pipeline,

First, please download all the scripts in same directory with R project.

Second, after click R project, please run script
**txt\_to\_csv\_VERSION.R** for converting raw txt file by typing
*`source`* button, or type *`source("txt_to_csv_VERSION.R")`*.

Third, please go to the main script **FWL\_lipidomics\_VERSION.R**.
Press the *`source`* button, or type
*`source("FWL_lipidomics_VERSION.R")`* in
console.

<!-- ![an image caption Source: screenshot.](button_source.png) -->

### Overflow

![](Manual_of_Lipidomics_pipeline_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

### Each Step Detail

##### Initialization

1.  Convert txt file to csv file for consistent analysis.

<!-- end list -->

``` r
source("txt_to_csv_version.R")
```

2.  Run main script **FWL\_lipidomics\_Version.R**

Console will list the files in the converted directory and it will pop
out a command line to ask user to input the index of the
    file.}

    The following files had been generated. Select ONE for subsequent the list of file names:
    1 20171127_JC_Seipin.csv
    2 20191220_CJ_s2.csv
    3 20200210_JSCJ_s2.csv
    4 BRAIIIN_single.csv
    
    Please input the index of the file:
