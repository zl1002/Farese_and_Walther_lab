---
title: "Manual"
author: "Wenting"
date: "2/21/2020"
output:
  html_document:
    keep_md: yes
    self_contained: no
  rmarkdown::github_document: default
---

```
{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(DiagrammeR)
```



### Motivation
The laboratory of Drs. Robert Farese, Jr. and Tobias Walther at the Harvard T.H. Chan School of Pulic Health has been using LipidSearch software (version 4.1) for some time, currently we build a pipeline based on Mass Spectrometry data processed by this software. We filtered the data with flexible standards, did QC, statistics analysis and visualization with extracted information from original txt file.

### Set up
Required R version 3.6.0 or later. 

For Mac users, XQaurtz need to be installed manually.

Any additional packages/libraries will be installed and uploaded automatically.

Pipeline will store all the output under 4 automatically generated directories, converted, data, plot, plot/classes.

To run the pipeline, 

First, please download all the scripts in same directory with R project. 

Second, after click R project, please run script **txt_to_csv_VERSION.R** for converting raw txt file by typing *`source`* button, or type *`source("txt_to_csv_VERSION.R")`*.

Third, please go to the main script **FWL_lipidomics_VERSION.R**. Press the *`source`* button, or type *`source("FWL_lipidomics_VERSION.R")`* in console.

<!-- ![an image caption Source: screenshot.](button_source.png) -->

### Overflow
<!--html_preserve--><div id="htmlwidget-31396468ef2f3aa18277" style="width:672px;height:480px;" class="grViz html-widget"></div>
<script type="application/json" data-for="htmlwidget-31396468ef2f3aa18277">{"x":{"diagram":"digraph flowchart {\n      # node definitions with substituted label text\n      node [fontname = Helvetica, shape = rectangle, style = filled]\n      tab1 [label = \"Raw data (txt file)\", fillcolor = Lavender]\n      node [fontname = Helvetica, shape = rectangle, style = filled]\n      tab2 [label = \"Initialization\", fillcolor = SteelBlue]\n      tab3 [label = \"Preprocess and Quality Check\", fillcolor = SteelBlue]\n      tab4 [label = \"Fatty Acids Saturation Analysis \nEther lipids abundance Analysis \nQuantification Analysis\", fillcolor = SteelBlue]\n      tab5 [label = \"Impute Lipid Molecules \nfor Differential Expression analysis\", fillcolor = SteelBlue]\n      node [fontname = Helvetica, shape = oval, style = filled]\n      tab6 [label = \"Background Subtration\", fillcolor = Plum]\n\n      # edge definitions with the node IDs\n      tab1 -> tab2-> tab3 -> tab6 -> tab4 -> tab5;\n      }","config":{"engine":"dot","options":null}},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->

### Each Step Detail

##### Initialization


1. Convert txt file to csv file for consistent analysis.

```r
source("txt_to_csv_version.R")
```
2. Run main script **FWL_lipidomics_Version.R**

  Console will list the files in the converted directory and it will pop out a command line to ask user to input the index of the file.}


```
The following files had been generated. Select ONE for subsequent the list of file names:
1 20171127_JC_Seipin.csv
2 20191220_CJ_s2.csv
3 20200210_JSCJ_s2.csv
4 BRAIIIN_single.csv

Please input the index of the file:
```




