# Exploring the gene expression using R Shiny App
This Shiny App is created to visualize gene expression data across all the Seurat objects used for the study **Barcoding Notch signaling in the developing brain** by Abigail Siniscalco, Roshan Priyarangana Perera, Jessie E. Greenslade, Aiden Masters, Hannah Doll, and Bushra Raj.
## Use the following steps to exploring the gene expression using R Shiny App.

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
