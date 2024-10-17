This repository contains the R script used for sub-clustering analyses, Supplementary Table 2, and a ShinyApp script for the scRNA-seq analysis used for the manuscript **Barcoding Notch signaling in the developing brain**, Abigail M. Siniscalco, Roshan Priyarangana Perera, Jessie E. Greenslade, Hemagowri Veeravenkatasubramanian,  Aiden Masters, Hannah M. Doll, and Bushra Raj.

All raw sequencing files and processed objects can be downloaded at GEO accession GSE268356.

# Overview
We provide an extensive scRNA-seq dataset collected from 21-23 dpf juvenile zebrafish brains.  We sequenced 156,353 cells and retained 148,853 cells for downstream analysis after quality filtering. Using a combination of initial coarse-grained computational clustering and subsequent rounds of iterative clustering, we classified 137 neuronal, non-neuronal and progenitor cell subtypes and cell states.  

### Master script
This R script details the iterative sub-clustering steps used for clustering analyses and plot generation.  To replicate this analysis with our sub-cluster annotations, the processed objects are available for download at GEO accession GSE268356.

### Supplementary Table 2
This table contains cell type annotations for the full dataset, as well as a list of DEGs for each sub-cluster.  (R script used to create DEG lists for each sub-cluster provided within the master script).

### Exploring the gene expression using R Shiny App
This Shiny App is created to visualize gene expression data across all the Seurat objects used for the study **Barcoding Notch signaling in the developing brain**. 
### Use the following steps to exploring the gene expression using R Shiny App.

1. Install the following R Packages: 
* install.packages("shiny")
* install.packages("shinyFiles")
* install.packages("Seurat") # this app created using Seurat v5.1.0
* install.packages("ggplot2")
2.	Download the Seurat Objects.
3.	Download the **app.R** Script.
4.	Open the app.R Script in RStudio.
5.	Click on the Run App button in RStudio to execute the app.R script.
6.	Use the file browsing field to upload your Seurat object.
7.	Upon successful loading, the Cluster UMAP will be displayed by default.
8.	Enter the gene name you are interested in the "Enter Gene Name" search field.
9.	The app will display the FeaturePlot and Violin Plot for the specified gene.
10.	Modify the "plot width" and "plot height" parameters to adjust the plot size as needed.
