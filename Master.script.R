library(Seurat)
#NOTE:  Run on Seurat version 5.0.1 & Seurat Object 5.0.1
library(SeuratObject)
library(SeuratWrappers)
library(patchwork)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scCustomize)
#update rlang >1.1.1 to run scCustomize package
#install.packages("rlang")

palette(brewer.pal(n = 9, name = "YlGnBu"))

#Due to the stochastic nature of certain processing steps, clustering outcomes may have slight variation from published objects/figures.
#The master script shows the processing and annotation steps for each sub-clustered object.
#Processed objects available on GEO.

# BROAD CELLS
w3.merge.res.1.0 <- readRDS("w3.merge.res.1.0.rds")
#25857 features across 148853 samples within 1 assay

ids.w3 <- c(    '0'='Midbrain',
                '1'='OT',
                '2'='Midbrain',
                '3'='Neuron',
                '4'='GC/TL',
                '5'='Prog_A',
                '6'='Pallium',
                '7'='Hyp',
                '8'='Prog_B',
                '9'='Prog_C',
                '10'='GC/TL',
                '11'='Pallium',
                '12'='Prog_D',
                '13'='GC/TL',
                '14'='Prog_E',
                '15'='RG',
                '16'='URL Progenitor',
                '17'='Prog_F',
                '18'='Microglia',
                '19'='Oligodendrocyte_A',
                '20'='Diencephalon',
                '21'='Ependymal_2',
                '22'='Unk_1',
                '23'='Pericytes',
                '24'='Ependymal_1',
                '25'='OPC',
                '26'='Endothelial_1',
                '27'='Oligodendrocyte_B',
                '28'='Neuron',
                '29'='Endothelial_2',
                '30'='Erythrocytes',
                '31'='Epithelial cells')

w3.merge.res.1.0.rename <- RenameIdents(w3.merge.res.1.0, ids.w3)

DimPlot(w3.merge.res.1.0.rename, label = T) + NoLegend()

#create metadata column with coarse clustering labels
w3.merge.res.1.0.rename[["coarse.ident"]] <- Idents(object = w3.merge.res.1.0.rename)

saveRDS(w3.merge.res.1.0.rename, "w3.merge.res.1.0.rename.rds")

#Optional: create table of cluster sizes
w3.table <- table(w3.merge.res.1.0.rename@meta.data$seurat_clusters)
write.csv(w3.table, "w3_cluster.size.csv")          


#############################################################################
# NEURO SUB

neur.sub <- subset(w3.merge.res.1.0.rename, idents = c("Midbrain", "OT", "Neuron", "Hyp", "Pallium", "Diencephalon"))

neur.sub <- NormalizeData(neur.sub)
neur.sub <- FindVariableFeatures(neur.sub, 
                                 nfeatures = 2000, 
                                 selection.method = "vst")

neur.sub <- ScaleData(neur.sub,vars.to.regress = "percent.mt")
neur.sub <- RunPCA(neur.sub, npcs = 50)
neur.sub <- RunUMAP(neur.sub, dims = 1:15)
neur.sub <- FindNeighbors(neur.sub)
neur.sub.res.3.0 <- FindClusters(neur.sub, resolution = 3.0)

neur.sub.res.3.0 <- JoinLayers(neur.sub.res.3.0)
#Find DEGs
neur.sub.res.3.0.deg <- FindAllMarkers(neur.sub.res.3.0, 
                                       only.pos = T,
                                       min.pct = 0.1, 
                                       min.diff.pct = 0.1, 
                                       logfc.threshold = 0.1)

#Order DEGs
neur.sub.res.3.0.deg  <- neur.sub.res.3.0.deg %>%
                         mutate(pct.diff=(pct.1-pct.2)) %>% 
                         group_by(cluster)  %>% 
                         arrange(desc(avg_log2FC), 
                         (pct.diff),.by_group = TRUE)

#Export DEGs
write.csv(neur.sub.res.3.0.deg,"neur.sub.res.3.0.deg.csv")
          
DimPlot(neur.sub.res.3.0)

saveRDS(neur.sub.res.3.0, "neur.sub.res.3.0_from.master.rds")

neur.sub.res.3.0
#25857 features across 73690 samples within 1 assay
#Note processed object has 74,566 cells compared to 73,690 due to change in some sub-cluster labels after processed neuron object was made.  
#Current labels are accurate and will not affect the sub-clustering results.


#############################################################################
#Neuron object 
#Processed object neur.sub.res.3.0_starter.rds

neuro.ids  <- c(  '0'='Neu(diff)',
                  '1'='OT',
                  '2'='M',
                  '3'='MTh',
                  '4'='Hyp',
                  '5'='Pa',
                  '6'='Hyp',
                  '7'='Hyp',
                  '8'='SPa',
                  '9'='MTh',
                  '10'='OT',
                  '11'='MTh',
                  '12'='Pa',
                  '13'='TL',
                  '14'='M',
                  '15'='Pa',
                  '16'='Hyp',
                  '17'='M',
                  '18'='SPa',
                  '19'='M',
                  '20'='M',
                  '21'='Pa',
                  '22'='OT',
                  '23'='Hb',
                  '24'='OT',
                  '25'='MTh',
                  '26'='Di',
                  '27'='Hyp',
                  '28'='Neuron_2',
                  '29'='MTh',
                  '30'='M',
                  '31'='Hb',
                  '32'='Pa',
                  '33'='DHab',
                  '34'='TL',
                  '35'='Di',
                  '36'='Hyp',
                  '37'='VHab',
                  '38'='PurN_1',
                  '39'='POA_1',
                  '40'='Pa',
                  '41'='Tel',
                  '42'='Pa',
                  '43'='OT',
                  '44'='Hyp',
                  '45'='Hyp',
                  '46'='Neu(diff)',
                  '47'='Tel',
                  '48'='Oligodendrocyte_2',
                  '49'='M',
                  '50'='Hyp',
                  '51'='Ganglia',
                  '52'='PurN_2',
                  '53'='Pa',
                  '54'='PurN_3',
                  '55'='PurN_4',
                  '56'='Oligodendrocyte_3',
                  '57'='PurN_5',
                  '58'='OB')


#Rename neur.sub.res.3.0 idents
neur.sub.res.3.0.rename <- RenameIdents(neur.sub.res.3.0, neuro.ids)

DimPlot(neur.sub.res.3.0.rename, label = T)

#create metadata column with coarse clustering labels
neur.sub.res.3.0.rename[["coarse.ident"]] <- Idents(object = neur.sub.res.3.0.rename)

saveRDS(neur.sub.res.3.0.rename, "neur.sub.res.3.0.rename.rds")


###########################################################################
#GC/TL
#Processed object GC.res.0.45.rds
#NOTE: GC.TL obj run on Seurat version 5.0.3 & Seurat Object 5.0.1

GC_TL <- subset(w3.merge.res.1.0.rename, idents = "GC/TL")

#subset cluster 13 and 34 (TL) from neurons

TL <- subset(neur.sub.res.3.0.rename, idents = "TL")

#merge
GC.TL <- merge(GC_TL, TL)

GC.TL <- NormalizeData(GC.TL)
GC.TL <- FindVariableFeatures(GC.TL, 
                               nfeatures = 2000, 
                               selection.method = "vst")

GC.TL <- ScaleData(GC.TL, vars.to.regress = "percent.mt")
GC.TL <- RunPCA(GC.TL, npcs = 50)
GC.TL <- RunUMAP(GC.TL, dims = 1:15)
GC.TL <- FindNeighbors(GC.TL)
GC.TL.res.0.45<- FindClusters(GC.TL, resolution = 0.45)
DimPlot(GC.TL.res.0.45, label = T)

#saveRDS(GC.TL.res.0.45, "GC.res.0.45.rds")

GC.TL.res.0.45 <- JoinLayers(GC.TL.res.0.45)

#Find DEGs
GC.TL.res.0.45.deg <- FindAllMarkers(GC.TL.res.0.45, 
                                     only.pos = T,
                                     min.pct = 0.1, 
                                     min.Telff.pct = 0.1, 
                                     logfc.threshold = 0.1)

#Reorder DEGs
GC.TL.res.0.45.deg <- GC.TL.res.0.45.deg %>% 
          mutate(pct.diff=(pct.1-pct.2)) %>% 
          group_by(cluster)  %>% 
          arrange(desc(avg_log2FC), 
          (pct.diff),.by_group = TRUE)

#export DEGs
write.csv(GC.TL.res.0.45.deg,"GC.TL.res.0.45.deg_master.csv")


GC.Mgenes <- c("neurod1", "cspg5b", "fat2", "lmo4b","otx1",
               "galn","pvalb7","unc5a",
               "tenm3","ebf3a",
               "vsnl1a","ntng2a",
               "insm1a", "oprd1b",
               "kcnd3","camkvb",
               "atp1a3b","vgf",
               "en1b","ptprz1b", 
               "kcnip3b","trpc3", "tac1", 
               "dscamb", "grin2aa",
               "lmx1al", "kif26ab")

#processed object cluster annotations
GC.ids <-c('4'='TL_1',
           '8'='TL_2',
           '0'='GC_1',
           '1'='GC_2',
           '2'='GC_3',
           '3'='GC_4',
           '5'='Neuron_cdh4+_roraa+', #minimal expression of GC/TL markers.
           '6'='GC_5',# markers such as tubb5, insm1a indicate an immature GC cell type.
           '7'='GC_6',
           '9'='GC_7',
           '10'='Neuron_neto1l+_cntn4+', #minimal expression of GC/TL markers.
           '11'='GC_8',
           '12'='GC_9',
           '13'='Unk_2') #minimal expression of GC/TL marker genes. Ambiguous mix of progenitor/neuron markers. 

GC.TL.res.0.45.rename <- RenameIdents(GC.TL.res.0.45,GC.ids)

idents.to.keep <- c("TL_1",
                    "TL_2",
                    "GC_1",
                    "GC_2",
                    "GC_3",
                    "GC_4",
                    "GC_5",
                    "GC_6",
                    "GC_7",
                    "GC_8",
                    "GC_9")

                                       
DotPlot_scCustom(GC.TL.res.0.45.rename,
                 GC.Mgenes,
                 idents= idents.to.keep,
                 dot.scale=8, 
                 dot.min = 0,
                 x_lab_rotate = T, 
                 scale.min = 0,
                 scale.max = 100,
                 colors_use = palette())   


plot.cells = WhichCells(GC.TL.res.0.45.rename, idents = idents.to.keep)
DimPlot(GC.TL.res.0.45.rename,cells = plot.cells, label = F, label.size = 5)+
  scale_color_manual(values =c('#F8766D',"#ABA300",'#28CECA','#00C19A',"#CD9600",'#FF61CC', '#31C51F', '#25aff5',
                               "#8494FF", '#C77CFF' , '#ff9a36')) 


saveRDS(GC.TL.res.0.45.rename, "GC.res.0.45.rds")


###########################################################################
# HYPOTHALAMUS
#Processed object Hyp.res.0.2.rds

hyp.sub <- subset(neur.sub.res.3.0.rename, idents = "Hyp")

hyp.sub <- NormalizeData(hyp.sub)
hyp.sub <- FindVariableFeatures(hyp.sub, nfeatures = 1500, 
                                    selection.method = "vst")
          
hyp.sub <- ScaleData(hyp.sub, vars.to.regress = "percent.mt")
hyp.sub <- RunPCA(hyp.sub, npcs = 50)
hyp.sub <- RunUMAP(hyp.sub, dims = 1:15)
hyp.sub <- FindNeighbors(hyp.sub)
hyp.sub.res.0.2 <- FindClusters(hyp.sub, resolution = 0.2)
DimPlot(hyp.sub.res.0.2)

#saveRDS(hyp.sub.res.0.2,"hyp.sub.res.0.2.rds")
          

hyp.sub.res.0.2 <- JoinLayers(hyp.sub.res.0.2)
hyp.sub.res.0.2.deg <- FindAllMarkers(hyp.sub.res.0.2, 
                                      only.pos = T,
                                      min.pct = 0.1, 
                                      min.Diff.pct = 0.1, 
                                      logfc.threshold = 0.1)

hyp.sub.res.0.2.deg <- hyp.sub.res.0.2.deg  %>% 
                       mutate(pct.Diff=(pct.1-pct.2)) %>% 
                       group_by(cluster)  %>% 
                       arrange(desc(avg_log2FC), 
                       (pct.Diff),.by_group = TRUE)

write.csv(hyp.sub.res.0.2.deg, "hyp.sub.res.0.2.deg.csv")


hyp.Mgenes <-c( "kctd4", "gad2", "slc6a1b", "slc32a1", "nkx2.4a",
                "scg2b", "tmsb",
                "tafa5a","crema", "nrn1a",
                "zfhx3",
                "ptgdsb.2","gad1b",
                "lhx5","nr5a1b", "htr2cl1",
                "prdx1","scg3","chgb",
                "gsx1",  "isl1", "dlx2a", "dlx5a")


#processed object cluster annotations
hyp.ids <- c("0"= "Hyp_1",
             "1"= "Hyp_2",
             "2"= "Hyp_3",
             "3"= "Hyp_4",
             "4"= "SubPal_7", 
             "5"= "Hyp_5",
             "6"= "Hyp_6",
             "7"= "SubPal_8")

# Further sub-clustering of the Hypothalamus resolved two SubPallium sub-clusters.
# SubPal_7 expression of synpr, sp8a aligns with known SubPallium cluster.
# SubPal_8 expression of lhx6, nkx2.1, sst1.1, & sox6 aligns with markers of interneurons derived from the SubPallium.

hyp.sub.res.0.2.rename <- RenameIdents(hyp.sub.res.0.2, hyp.ids)

idents.to.keep <- c("Hyp_1", 
                    "Hyp_2",
                    "Hyp_3",
                    "Hyp_4",
                    "Hyp_5",
                    "Hyp_6")


DotPlot_scCustom(hyp.res.0.2.rename,
                 hyp.Mgenes,
                 idents= idents.to.keep,
                 dot.scale=8, 
                 dot.min = 0,  
                 x_lab_rotate = T, 
                 scale.min = 0,
                 scale.max = 100,
                 colors_use = palette())  


plot.cells = WhichCells(hyp.res.0.2.rename, idents = idents.to.keep)
DimPlot(hyp.res.0.2.rename, cells = plot.cells, label = T, label.size = 5)

saveRDS(hyp.sub.res.0.2.rename, "hyp.sub.res.0.2.rename.rds")          


###########################################################################
# OPTIC TECTUM
#processed object OT.res.0.2.rds

opt.sub <- subset(neur.sub.res.3.0.rename, idents = "OT")

opt.sub <- NormalizeData(opt.sub)
opt.sub <- FindVariableFeatures(opt.sub,nfeatures = 2000,
                                selection.method = "vst")

opt.sub <- ScaleData(opt.sub,vars.to.regress = "percent.mt")
opt.sub <- RunPCA(opt.sub, npcs = 50)
opt.sub <- RunUMAP(opt.sub, dims = 1:15)
opt.sub <- FindNeighbors(opt.sub)
opt.sub.res.0.2 <- FindClusters(opt.sub, resolution = 0.2)
opt.sub.res.0.5 <- FindClusters(opt.sub, resolution = 0.5)

DimPlot(opt.sub.res.0.2, label = T)

#saveRDS(opt.sub.res.0.2, "opt.sub.res.0.2.rds")
#saveRDS(opt.sub.res.0.5, "opt.sub.res.0.5.rds")


opt.sub.res.0.2 <- JoinLayers(opt.sub.res.0.2)
opt.sub.res.0.2.deg <- FindAllMarkers(opt.sub.res.0.2, 
                                      only.pos = T,
                                      min.pct = 0.1, 
                                      min.diff.pct = 0.1,
                                      logfc.threshold = 0.1)


opt.sub.res.0.2.deg <- opt.sub.res.0.2.deg %>% 
                       mutate(pct.diff=(pct.1-pct.2)) %>% 
                       group_by(cluster)  %>% 
                       arrange(desc(avg_log2FC), (pct.diff),.by_group = TRUE)

write.csv(opt.sub.res.0.2.deg, "opt.sub.res.0.2.deg.csv")

OT.Mgenes <- c("pax7a", "sox14", "gata3", "gata2a", 
               "gad1b", "emx2",
               "tal1", "otx2b", "pax7b","slc6a1b","gad2","tubb5",
               "grm8a", "il1rapl2",
               "pbx4","id2a","irx3a",
               "nlgn1","slc6a1a","si:dkey-7j14.5",  
               "znf536","crabp1a", "tns1b",
               "bcam","robo1",
               "pcbp3","col14a1a", "gjd2b")

#processed object cluster annotations
OT.ids <-c("0"= "OT_1",
           "1"= "OT_2",
           "2"= "OT_3",
           "3"= "OT_4",
           "4"= "OT_5",
           "5"= "OT_6")

opt.sub.res.0.2.rename <- RenameIdents(opt.sub.res.0.2,OT.ids)

DotPlot_scCustom(opt.sub.res.0.2.rename,
                 OT.Mgenes,
                 dot.scale=8, 
                 dot.min = 0,  
                 x_lab_rotate = T, 
                 scale.min = 0,
                 scale.max = 100,
                 colors_use = palette())         


DimPlot(opt.sub.res.0.2.rename,label = T, label.size = 5)

saveRDS(opt.sub.res.0.2.rename, "OT.res.0.2.rds") 

#Cluster 9 from opt.sub.res.0.5 showed HB markers.  Labeled as sub-cluster HB_2 in final object.
Cluster9.hb <- WhichCells(opt.sub.res.0.5, idents = "9")

saveRDS(Cluster9.hb, "Cluster9.hb.rds")


###########################################################################
# MIDBRAIN
#processed obj MB.res.0.1.rds

mb.sub <- subset(neur.sub.res.3.0.rename, idents = c("M","MTh"))

mb.sub<- NormalizeData(mb.sub)
mb.sub<- FindVariableFeatures(mb.sub, nfeatures = 2000,
                              selection.method = "vst")

mb.sub <- ScaleData(mb.sub, vars.to.regress = "percent.mt")
mb.sub <- RunPCA(mb.sub, npcs = 50)
mb.sub <- RunUMAP(mb.sub, dims = 1:15)
mb.sub <- FindNeighbors(mb.sub)
mb.sub.res.0.1 <- FindClusters(mb.sub, resolution = 0.1)
DimPlot(mb.sub.res.0.1, label = T)

#saveRDS(mb.sub.res.0.1, "mb.sub.res.0.1.rds")

mb.sub.res.0.1 <- JoinLayers(mb.sub.res.0.1)
mb.sub.res.0.1.deg <- FindAllMarkers(mb.sub.res.0.1,
                                     only.pos = T,
                                     min.pct = 0.1,
                                     min.diff.pct = 0.1,
                                     logfc.threshold = 0.1)

mb.sub.res.0.1.deg  <- mb.sub.res.0.1.deg  %>% 
                       mutate(pct.diff=(pct.1-pct.2)) %>% 
                       group_by(cluster)  %>% 
                       arrange(desc(avg_log2FC), 
                       (pct.diff),.by_group = TRUE)

write.csv(mb.sub.res.0.1.deg, "mb.sub.res.0.1.deg.csv")


MB.Mgenes <- c( "robo1", "islr2",  "tcf7l2",
                "lhx9", "barhl2","pou3f1", 
                "slc17a6b", "ddit3",
                "tfap2a", "shox2", "onecut1",
                "nefma", "dkk3a", "pvalb6",
                "calb2a", "lef1")

#processed object cluster annotations
mb.ids <- c("0"= "MB_1",
            "1"= "MB_2",
            "2"= "MB_3",
            "3"= "MB_4",
            "4"= "MB_5"
)

mb.sub.res.0.1.rename <- RenameIdents(mb.sub.res.0.1,mb.ids)


DotPlot_scCustom(mb.sub.res.0.1.rename,
                 MB.Mgenes,
                 dot.scale=8, 
                 dot.min = 0,  
                 x_lab_rotate = T, 
                 scale.min = 0,
                 scale.max = 100,
                 colors_use = palette())

DimPlot(mb.sub.res.0.1.rename,label = T, label.size = 5)

saveRDS(mb.sub.res.0.1.rename,"mb.sub.res.0.1.rename.rds")


###########################################################################
# PALLIUM
#processed object Pallium.res.0.5.rds

pallium.sub <- subset(neur.sub.res.3.0.rename, idents = c("Pa"))

pallium.sub <- NormalizeData(pallium.sub)
pallium.sub <- FindVariableFeatures(pallium.sub, 
                                    nfeatures = 2000, 
                                    selection.method = "vst")

pallium.sub <- ScaleData(pallium.sub, vars.to.regress = "percent.mt")

pallium.sub <- RunPCA(pallium.sub, npcs = 50)
pallium.sub <- RunUMAP(pallium.sub, dims = 1:15)
pallium.sub <- FindNeighbors(pallium.sub)
pallium.sub.res.0.5 <- FindClusters(pallium.sub, resolution = 0.5)
DimPlot(pallium.sub.res.0.5, label = T)

#saveRDS(pallium.sub.res.0.5 ,"pallium.sub.res.0.5.rds")


pallium.sub.res.0.5 <- JoinLayers(pallium.sub.res.0.5)
pallium.sub.res.0.5.deg <- FindAllMarkers(pallium.sub.res.0.5, 
                                          only.pos = T,
                                          min.pct = 0.1, 
                                          min.diff.pct = 0.1, 
                                          logfc.threshold = 0.1)


pallium.sub.res.0.5.deg  <- pallium.sub.res.0.5.deg %>% 
                            mutate(pct.diff=(pct.1-pct.2)) %>% 
                            group_by(cluster)  %>% 
                            arrange(desc(avg_log2FC), 
                            (pct.diff),.by_group = TRUE)

write.csv(pallium.sub.res.0.5.deg, "pallium.sub.res.0.5.deg.csv")


Pal.Mgenes <- c("neurod1", "tbr1b","eomesa",
                "bhlhe22", "zbtb18", "emx3","tubb5",
                "rcc2", "nhlh2",
                "fndc5b","c1ql3b.1", 
                "lhx9",  "tbr1a", 
                "myca", "rhbdf1b",
                "erbb4b", "cntnap5b",
                "limch1b", "slc24a3", 
                "ebf1b", "nfia",
                "mcama","eno3",
                "nrgna", 
                "pvalb7", "chrm2a", 
                "foxp4","hpcal4",
                "grm2a","msi2b","ebf3a.1", "lhx1a")

# removed Cluster 0 from further Pallium analysis due to extremely limited expression of any Pallium markers genes including neurod1, tbr1b, eomesa, bhlhe22, zbtb18, emx3.
# Cluster 0 re-labeled as "miscellaneous" sub-cluster below

pallium.sub.res.0.5.final <- subset(x = pallium.sub.res.0.5, 
                                  idents = c("1",
                                             "2",
                                             "3",
                                             "4",
                                             "5",
                                             "6",
                                             "7",
                                             "8",
                                             "9",
                                             "10",
                                             "11",
                                             "12"))

DimPlot(pallium.sub.res.0.5.final, label = T)

#processed object cluster annotations
Pal.ids <-c("1"= "Pal_1",
            "2"= "Pal_2",
            "3"= "Pal_3",
            "4"= "Pal_4",
            "5"= "Pal_5",
            "6"= "Pal_6",
            "7"= "Pal_7",
            "8"= "Pal_8",
            "9"= "Neuron_grm2a+",
            "10"= "Pal_9",
            "11"= "Pal_10",
            "12"= "Pal_11")

#cluster 5 expressing markers associated with newborn neurons (NBN) (msi2b,ebf3a.1, rtn4r, etv5a), as well as OB marker grm2a.
#labeled Neuron_grm2a+ due to likelihood of both NBN & OB cell types present but unable to be resolved.

Pal.res.0.5.rename <- RenameIdents(pallium.sub.res.0.5.final,Pal.ids)

levels(Pal.res.0.5.rename)

levels(x = Pal.res.0.5.rename) <- c("Pal_1", 
                                    "Pal_2", 
                                    "Pal_3", 
                                    "Pal_4", 
                                    "Pal_5", 
                                    "Pal_6", 
                                    "Pal_7", 
                                    "Pal_8", 
                                    "Pal_9", 
                                    "Pal_10", 
                                    "Pal_11", 
                                    "Neuron_grm2a+")


DotPlot_scCustom(Pal.res.0.5.rename,
                 Pal.Mgenes,
                 dot.scale=8, 
                 dot.min = 0,  
                 x_lab_rotate = T, 
                 scale.min = 0,
                 scale.max = 100,
                 colors_use = palette()) 

DimPlot(Pal.res.0.5.rename, label = T, label.size = 5)

saveRDS(Pal.res.0.5.rename, "Pallium.res.0.5.rds")


###########################################################################
# MISCELLANEOUS
#processed object Misc.res.0.3.rds
#subset the cluster 0 cells that were removed from the Pallium analysis.

misc.sub <- subset(pallium.sub.res.0.5, idents = c("0"))

misc.sub <- NormalizeData(misc.sub)
misc.sub <- FindVariableFeatures(misc.sub, 
                               nfeatures = 1500, 
                               selection.method = "vst")

misc.sub <- ScaleData(misc.sub, vars.to.regress = "percent.mt")
misc.sub <- RunPCA(misc.sub, npcs = 50)
misc.sub <- RunUMAP(misc.sub, dims = 1:15)
misc.sub <- FindNeighbors(misc.sub)
misc.sub.res.0.3 <- FindClusters(misc.sub, resolution = 0.3)
DimPlot(misc.sub.res.0.3, label = T)

#saveRDS(misc.sub.res.0.3, "misc.sub.res.0.3.rds")


misc.sub <- JoinLayers(misc.sub)
misc.sub.res.0.3.deg <- FindAllMarkers(misc.sub.res.0.3, 
                                      only.pos = T,
                                      min.pct = 0.1, 
                                      min.diff.pct = 0.1, 
                                      logfc.threshold = 0.1)


misc.sub.res.0.3.deg <- misc.sub.res.0.3.deg %>% 
                        mutate(pct.diff=(pct.1-pct.2)) %>% 
                        group_by(cluster)  %>% 
                        arrange(desc(avg_log2FC), 
                        (pct.diff),.by_group = TRUE)

write.csv(misc.sub.res.0.3.deg,"misc.sub.res.0.3.deg.csv")

#processed object cluster annotations
misc.ids <-c('0'='MB_het',
           '1'='MB_het',
           '5'='Neuron_eomesa+_zic2a+',
           '3'='Di_7',
           '4'='HB_4',
           '2'='Neuron_zbtb18+_neurod1+')

misc.sub.res.0.3.rename <- RenameIdents(misc.sub.res.0.3, misc.ids)

DimPlot(misc.sub.res.0.3.rename, label = T)

saveRDS(misc.sub.res.0.3.rename, "Misc.res.0.3.rds")


###########################################################################
# SUB-PALLIUM
#processed object SubPal.res.0.4.rds

SPa.sub <- subset(neur.sub.res.3.0.rename, idents = "SPa")

SPa.sub<- NormalizeData(SPa.sub)
SPa.sub <- FindVariableFeatures(SPa.sub, 
                                nfeatures = 2000, 
                                selection.method = "vst")

SPa.sub <- ScaleData(SPa.sub,vars.to.regress = "percent.mt")
SPa.sub <- RunPCA(SPa.sub, npcs = 50)
SPa.sub <- RunUMAP(SPa.sub, dims = 1:15)
SPa.sub <- FindNeighbors(SPa.sub)
SPa.sub.res.0.4 <- FindClusters(SPa.sub, resolution = 0.4)
DimPlot(SPa.sub.res.0.4, label = T)

#saveRDS(SPa.sub.res.0.4 ,"SPa.sub.res.0.4.rds")


SPa.sub.res.0.4 <- JoinLayers(SPa.sub.res.0.4)
SPa.sub.res.0.4.deg <- FindAllMarkers(SPa.sub.res.0.4, 
                                      only.pos = T,
                                      min.pct = 0.1, 
                                      min.diff.pct = 0.1, 
                                      logfc.threshold = 0.1)

SPa.sub.res.0.4.deg  <- SPa.sub.res.0.4.deg  %>% 
                        mutate(pct.diff=(pct.1-pct.2)) %>% 
                        group_by(cluster)  %>% 
                        arrange(desc(avg_log2FC), (pct.diff),.by_group = TRUE)

write.csv(SPa.sub.res.0.4.deg, "SPa.sub.res.0.4.deg.csv")


SubPal.Mgenes <- c( "dlx2a","dlx5a","dlx6a", "gad2",
                    "nrp2b","grid1b",
                    "arrdc3b", "atf4b",
                    "mibp2", "dlx1a", 
                    "lhx6","sox6","lhx8a", 
                    "nr2f2", "nfia", 
                    "erbb4b","mpped1",
                    "six3a","six3b")

#processed object cluster annotations
SPa.ids <-c("0"= "Neuron_1",#does not express SubPallium markers such as, dlx2a, dlx5a, dlx6a
             "1"= "SubPal_1",
             "2"= "SubPal_2",
             "3"= "SubPal_3",
             "4"= "SubPal_4",
             "5"= "SubPal_5",
             "6"= "SubPal_6")

SPa.sub.res.0.4.rename <- RenameIdents(SPa.sub.res.0.4, SPa.ids)


DotPlot_scCustom(SPa.sub.res.0.4.rename,
                 SubPal.Mgenes,
                 dot.scale=8, 
                 dot.min = 0, 
                 x_lab_rotate = T, 
                 scale.min = 0,
                 scale.max = 100,
                 colors_use = palette())


DimPlot(SPa.sub.res.0.4.rename, label = T, label.size = 5)

saveRDS(SPa.sub.res.0.4.rename, "SPa.sub.res.0.4.rename.rds")


###########################################################################
# VENTRAL HABENULA
#processed object VHab.res.0.7.rds
#neur.sub.res.3.0.rename <-JoinLayers(neur.sub.res.3.0.rename)

vH.sub <- subset(neur.sub.res.3.0.rename, idents = c("VHab"))

vH.sub<- NormalizeData(vH.sub)
vH.sub <- FindVariableFeatures(vH.sub,
                               nfeatures = 2000, 
                               selection.method = "vst")

vH.sub <- ScaleData(vH.sub, vars.to.regress = "percent.mt")
vH.sub <- RunPCA(vH.sub, npcs = 50)
vH.sub <- RunUMAP(vH.sub, dims = 1:15)
vH.sub <- FindNeighbors(vH.sub)
vH.sub.res.0.7 <- FindClusters(vH.sub, resolution = 0.7)
DimPlot(vH.sub.res.0.7, label = T)

#saveRDS(vH.sub.res.0.7, "vH.sub.res.0.7.rds")

vH.sub.res.0.7 <- JoinLayers(vH.sub.res.0.7)
vH.sub.res.0.7.deg <- FindAllMarkers(vH.sub.res.0.7, 
                                     only.pos = T,
                                     min.pct = 0.1, 
                                     min.diff.pct = 0.1, 
                                     logfc.threshold = 0.1)

vH.sub.res.0.7.deg <- vH.sub.res.0.7.deg %>% 
                      mutate(pct.diff=(pct.1-pct.2)) %>% 
                      group_by(cluster)  %>% 
                      arrange(desc(avg_log2FC), 
                      (pct.diff),.by_group = TRUE)

write.csv(vH.sub.res.0.7.deg, "vH.sub.res.0.7.deg.csv")

vH.Mgenes <- c("kiss1", "prkcq", "ppp1r14ab",
               "pcdh17", "synpr",
               "hist2h2l", "arrdc3b",
               "egr2b",
               "rbfox1", "nfixb", "zfpm2a",
               "robo1", "cntn2",
               "agrn",
               "atp2b1a","tubb5")

#processed object cluster annotations
vH.ids <-c("0"= "VHab_1",
           "1"= "VHab_2",
           "2"= "VHab_3",
           "3"= "VHab_4",
           "4"= "VHab_5",
           "5"= "VHab_6")

vH.sub.res.0.7.rename <- RenameIdents(vH.sub.res.0.7,vH.ids)


DotPlot_scCustom(vH.sub.res.0.7.rename, 
                 features = vH.Mgenes,
                 dot.scale=8, 
                 dot.min = 0, 
                 x_lab_rotate = T, 
                 scale.min = 0,
                 scale.max = 100,
                 colors_use = palette())


DimPlot(vH.sub.res.0.7.rename, label = T, label.size = 5)

saveRDS(vH.sub.res.0.7.rename,"VHab.res.0.7.rds")


###########################################################################
# DORSAL HABENULA
#processed object DHab.res.0.4.rds
#neur.sub.res.3.0.rename <-JoinLayers(neur.sub.res.3.0.rename)

dH.sub <- subset(neur.sub.res.3.0.rename, idents = "DHab")

dH.sub<- NormalizeData(dH.sub)
dH.sub <- FindVariableFeatures(dH.sub, 
                               nfeatures = 2000, 
                               selection.method = "vst")

dH.sub <- ScaleData(dH.sub, vars.to.regress = "percent.mt")
dH.sub <- RunPCA(dH.sub, npcs = 50)
#DimHeatmap(dH.sub, dims = 1:10 , cells = 500, balanced = TRUE)
#used dims 1:8 based in the PC heatmaps
dH.sub <- RunUMAP(dH.sub, dims = 1:8)
dH.sub <- FindNeighbors(dH.sub)
dH.sub.res.0.4 <- FindClusters(dH.sub, resolution = 0.4)
DimPlot(dH.sub.res.0.4, label = T)

#saveRDS(dH.sub.res.0.4, "dH.sub.res.0.4.rds")

dH.sub.res.0.4 <- JoinLayers(dH.sub.res.0.4)
dH.sub.res.0.4.deg <- FindAllMarkers(dH.sub.res.0.4, 
                                     only.pos = T,
                                     min.pct = 0.1, 
                                     min.diff.pct = 0.1, 
                                     logfc.threshold = 0.1)


dH.sub.res.0.4.deg <- dH.sub.res.0.4.deg %>% 
                      mutate(pct.diff=(pct.1-pct.2)) %>% 
                      group_by(cluster)  %>% 
                      arrange(desc(avg_log2FC), (pct.diff),.by_group = TRUE)

write.csv(dH.sub.res.0.4.deg, "dH.sub.res.0.4.deg.csv")


dH.Mgenes <- c("gng8","synpr","tac3a", "sst1.1","tubb5",
               "fxyd1", "g0s2", 
               "zfhx4", "tuba8l3", 
               "vav3b","mpped1", 
               "luzp2","wnt7aa", 
               "grid2","adcyap1a", 
               "lrrtm1", "foxa1",
               "fabp7a", "her4.1",
               "id4","cbln2b", "tnr", 
               "sdc4", "c1qtnf4", "atf3" )

#processed object cluster annotations
dH.ids <-c("0"= "DHab_1",
           "1"= "DHab_2",
           "2"= "DHab_3",
           "3"= "DHab_4",
           "4"= "DHab_5",
           "5"= "DHab_6",
           "6"= "DHab_7",
           "7"= "DHab_8",
           "8"= "DHab_9")

dH.sub.res.0.4.rename <- RenameIdents(dH.sub.res.0.4,dH.ids)


DotPlot_scCustom(dH.sub.res.0.4.rename,
                 dH.Mgenes, 
                 dot.scale=8, 
                 dot.min = 0, 
                 x_lab_rotate = T, 
                 scale.min = 0,
                 scale.max = 100,
                 colors_use = palette())

DimPlot(dH.sub.res.0.4.rename, label = T, label.size = 5)

saveRDS(dH.sub.res.0.4.rename, "DHab.res.0.4.rds")


###########################################################################
# TELENCEPHALON
#processed object Tel.res.0.4.rds

Tel.sub <- subset(neur.sub.res.3.0.rename, idents= c("Tel"))

Tel.sub <- NormalizeData(Tel.sub)
Tel.sub <- FindVariableFeatures(Tel.sub, 
                                nfeatures = 2000, 
                                selection.method = "vst")

Tel.sub <- ScaleData(Tel.sub, vars.to.regress = "percent.mt")
Tel.sub <- RunPCA(Tel.sub, npcs = 50)
#DimHeatmap(Tel.sub, dims = 1:12 , cells = 500, balanced = TRUE)
Tel.sub <- RunUMAP(Tel.sub, dims = 1:9)
Tel.sub <- FindNeighbors(Tel.sub)
Tel.sub.res.0.4 <- FindClusters(Tel.sub, resolution = 0.4)
DimPlot(Tel.sub.res.0.4, label = T)

#saveRDS(Tel.sub.res.0.4, "Tel.sub.res.0.4.rds")


Tel.sub.res.0.4 <- JoinLayers(Tel.sub.res.0.4)
Tel.sub.res.0.4.deg <- FindAllMarkers(Tel.sub.res.0.4, 
                                      only.pos = T,
                                      min.pct = 0.1, 
                                      min.Telff.pct = 0.1, 
                                      logfc.threshold = 0.1)


Tel.sub.res.0.4.deg <- Tel.sub.res.0.4.deg %>% 
                       mutate(pct.Diff=(pct.1-pct.2)) %>% 
                       group_by(cluster)  %>% 
                       arrange(desc(avg_log2FC), 
                       (pct.Diff),.by_group = TRUE)

write.csv(Tel.sub.res.0.4.deg,"Tel.sub.res.0.4.deg.csv")


Tel.genes <- c( "tubb5", "elavl4", "elavl3",
                "tbr1b", "lhx5", "lhx1a",
                "uncx", "uncx4.1", "slc17a6b",
                "tac1", "pou6f2","zeb2a", "mef2cb")

#processed object cluster annotations
Tel.ids <-c('0'='Tel_IN',
            '1'='Tel',
            '2'='Neuron_zeb2a+',
            '3'='Neuron_mef2cb+',
            '4'='Hyp_7',#Hyp markers: nkx2.4, gad1b, gad2, tac1
            '5'='POA_2',
            '6'='Neuron_tac1+_pou6f2+',
            '7'='HB_3')

#most telencephalon neurons represented as pallium or subpallium.
#only cluster 0 & 1 showed some Telencephalon markers, as well as have markers of immature/new neurons.
#re-labeled clusters 2-7 based on DEG markers. 

Tel.res.0.4.rename <- RenameIdents(Tel.res.0.4.rename, Tel.ids)


DotPlot_scCustom(Tel.res.0.4.rename, 
                 features = Tel.genes,
                 dot.scale=8, 
                 dot.min = 0,  
                 x_lab_rotate = T, 
                 scale.min = 0,
                 scale.max = 100,
                 colors_use = palette())

 
DimPlot(Tel.res.0.4.rename, label = T, label.size = 5)

saveRDS(Tel.res.0.4.rename, "Tel.res.0.4.rds")


###########################################################################
# DIENCEPHALON
#processed object Di.res.0.2.rds

Di.sub <- subset(neur.sub.res.3.0.rename, idents= c("Di"))

Di.sub <- NormalizeData(Di.sub)
Di.sub <- FindVariableFeatures(Di.sub, 
                               nfeatures = 2000, 
                               selection.method = "vst")

Di.sub <- ScaleData(Di.sub, vars.to.regress = "percent.mt")
Di.sub <- RunPCA(Di.sub, npcs = 50)
Di.sub <- RunUMAP(Di.sub, dims = 1:15)
Di.sub <- FindNeighbors(Di.sub)
Di.sub.res.0.2 <- FindClusters(Di.sub, resolution = 0.2)
DimPlot(Di.sub.res.0.2, label = T)

#saveRDS(Di.sub.res.0.2, "Di.sub.res.0.2.rds")

Di.sub.res.0.2 <- JoinLayers(Di.sub.res.0.2)
Di.sub.res.0.2.deg <- FindAllMarkers(Di.sub.res.0.2, 
                                        only.pos = T,
                                        min.pct = 0.1, 
                                        min.diff.pct = 0.1, 
                                        logfc.threshold = 0.1)


Di.sub.res.0.2.deg <- Di.sub.res.0.2.deg %>% 
                      mutate(pct.diff=(pct.1-pct.2)) %>% 
                      group_by(cluster)  %>% 
                      arrange(desc(avg_log2FC), 
                      (pct.diff),.by_group = TRUE)

write.csv(Di.sub.res.0.2.deg,"Di.sub.res.0.2.deg.csv")


Di.genes <- c( "pitx2", "nme2a", "CR936442.1", 
               "calb2a", "neurod6a","slc17a6b",
               "tubb5", "elavl3", "foxa1",  
               "lmo4a","irx2a", 
               "bhlhe22", 
               "rps6", "rpl3", "rps17","rgmb", 
               "epha6", "pvalb7", "mgll", "rbfox3b","tpma", 
               "bcam", "zeb2a","kidins220a", "slit3","atf4b",
               "scrt1b", "dusp19b", "arrdc3b",
               "rasd4", "nrgna", "calb1", "tmsb1","atp2b4","uchl1",
               "prkcda",
               "col14a1a", "gfra4a","rergla")

#processed object cluster annotations
Di.ids <-c("0"= "Di_1_IN",
           "1"= "Di_2",
           "2"= "Di_3",
           "3"= "Di_4",
           "4"= "Di_5",
           "5"= "Di_6")

Di.sub.res.0.2.rename <- RenameIdents(Di.sub.res.0.2, Di.ids)


DotPlot_scCustom(Di.sub.res.0.2.rename, 
                 features = Di.genes,
                 dot.scale=8, 
                 dot.min = 0,  
                 x_lab_rotate = T, 
                 scale.min = 0,
                 scale.max = 100,
                 colors_use = palette())


DimPlot(Di.sub.res.0.2.rename,label = T, label.size = 5)

saveRDS(Di.sub.res.0.2.rename, "Di.sub.res.0.2.rename.rds")


###########################################################################
# HINDBRAIN
#processed object HB.res.0.2.rds

HB.sub<- subset(neur.sub.res.3.0.rename, idents= c("Hb"))

HB.sub <- NormalizeData(HB.sub)
HB.sub <- FindVariableFeatures(HB.sub, 
                               nfeatures = 2000, 
                               selection.method = "vst")

HB.sub <- ScaleData(HB.sub,vars.to.regress = "percent.mt")
HB.sub <- RunPCA(HB.sub, npcs = 50)
HB.sub <- RunUMAP(HB.sub, dims = 1:15)
HB.sub <- FindNeighbors(HB.sub)
HB.sub.res.0.2 <- FindClusters(HB.sub, resolution = 0.2)
DimPlot(HB.sub.res.0.2, label = T)

#saveRDS(HB.sub.res.0.2, "HB.sub.res.0.2.rds")


HB.sub.res.0.2 <- JoinLayers(HB.sub.res.0.2)
HB.sub.res.0.2.deg <- FindAllMarkers(HB.sub.res.0.2, 
                                        only.pos = T,
                                        min.pct = 0.1, 
                                        min.HBff.pct = 0.1, 
                                        logfc.threshold = 0.1)


HB.sub.res.0.2.deg <- HB.sub.res.0.2.deg %>% 
                      mutate(pct.HBff=(pct.1-pct.2)) %>% 
                      group_by(cluster)  %>% 
                      arrange(desc(avg_log2FC), 
                      (pct.HBff),.by_group = TRUE)

write.csv(HB.sub.res.0.2.deg,"HB.sub.res.0.2.deg.csv")


HB.genes <- c( "lmx1ba", "lmx1bb", "phox2bb", 
               "hoxb3a", "hoxb5a", "hoxb5b",  "lmo4a", 
               "gad2", "gad1b", "slc32a1", "neurod1", "neurod6b", "lhx5",
               "shox2", "tlx2", "hoxb2a",
               "pcdh11", "pcdh7b", "nlgn1", "ptprua",
               "nfixb", "her4.2.1", "zbtb18", "nfixa",
               "robo3", "foxp4", "mmp17a", "dbx1a",
               "phox2a", "tlx3b",
               "adarb2", "tafa5a","grm8a", "pax2b",
               "gabrz", "galnt18a", "adcy2a", "foxb1a")

#processed object cluster annotations
HB.ids <-c("0"= "Neuron_fstl5+",
           "1"= "Prog_4",
           "2"= "HB_5",
           "3"= "HB_1",
           "4"= "Neuron_adarb2+",
           "5"= "HB_6")

#Non HB clusters lacked HB markers

HB.sub.res.0.2.rename <- RenameIdents(HB.sub.res.0.2, HB.ids)


DotPlot_scCustom(HB.sub.res.0.2.rename ,
                 features = HB.genes,
                 dot.scale=8, 
                 dot.min = 0,  
                 x_lab_rotate = T, 
                 scale.min = 0,
                 scale.max = 100,
                 colors_use = palette())


DimPlot(HB.sub.res.0.2.rename, label = T, label.size = 5)

saveRDS(HB.sub.res.0.2.rename, "HB.res.0.2.rds")


###########################################################################
# NEURONS DIFFERENTIATING
#processed object Neur.Diff.res.0.3.rds

####Subset
Neur_D.sub<- subset(neur.sub.res.3.0.rename, idents= c("Neu(diff)"))

Neur_D.sub <- NormalizeData(Neur_D.sub)
Neur_D.sub <- FindVariableFeatures(Neur_D.sub, 
                                   nfeatures = 2000, 
                                   selection.method = "vst")

Neur_D.sub <- ScaleData(Neur_D.sub,vars.to.regress = "percent.mt")
Neur_D.sub <- RunPCA(Neur_D.sub, npcs = 50)
#DimHeatmap(Neur_D.sub, dims = 10:15 , cells = 500, balanced = TRUE)
Neur_D.sub <- RunUMAP(Neur_D.sub, dims = 1:12)
Neur_D.sub <- FindNeighbors(Neur_D.sub)
Neur_D.sub.res.0.3 <- FindClusters(Neur_D.sub, resolution = 0.3)
DimPlot(Neur_D.sub.res.0.3, label = T)

#saveRDS(Neur_D.sub.res.0.3, "Neur_D.sub.res.0.3.rds")


Neur_D.sub.res.0.3 <- JoinLayers(Neur_D.sub.res.0.3)
Neur_D.sub.res.0.3.deg <- FindAllMarkers(Neur_D.sub.res.0.3, 
                                            only.pos = T,
                                            min.pct = 0.1, 
                                            min.Neur_Dff.pct = 0.1, 
                                            logfc.threshold = 0.1)

Neur_D.sub.res.0.3.deg <- Neur_D.sub.res.0.3.deg %>% 
                          mutate(pct.Neur_Dff=(pct.1-pct.2)) %>% 
                          group_by(cluster)  %>% 
                          arrange(desc(avg_log2FC), 
                          (pct.Neur_Dff),.by_group = TRUE)


write.csv(Neur_D.sub.res.0.3.deg,"Neur_D.sub.res.0.3.deg.csv")


ND.genes <- c( "tubb5", "elavl4", "elavl3")

#processed object cluster annotations
Neur_Diff.ids <-c('0'='IN_1',
                  '1'='Precursor_1',
                  '2'='Precursor_2',
                  '3'='Precursor_3',
                  '4'='IN_2',
                  '5'='Prog_5',#no differentiating genes (tubb5, elavl3/4)
                  '6'='Pineal')

Neur_Diff.res.0.3.rename <- RenameIdents(Neur_D.sub.res.0.3, Neur_Diff.ids)

DotPlot_scCustom(Neur_Diff.res.0.3.rename, 
                 features = ND.genes,
                 dot.scale=8, 
                 dot.min = 0,  
                 x_lab_rotate = T, 
                 scale.min = 0,
                 scale.max = 100,
                 colors_use = palette()) 


DimPlot(Neur_Diff.res.0.3.rename, label = T, label.size = 5)

saveRDS(Neur_Diff.res.0.3.rename, "Neur.Diff.res.0.3.rds")


###########################################################################
# PROGENITORS
#processed object Prog.res.0.13.rds
#make sure to subset from full w3 object
prog.sub <- subset(w3.merge.res.1.0.rename, idents = c("Prog_A", 
                                                       "Prog_B", 
                                                       "Prog_C", 
                                                       "Prog_D", 
                                                       "Prog_E", 
                                                       "Prog_F",
                                                       "URL Progenitor",
                                                       "RG"))

prog.sub <- NormalizeData(prog.sub)
prog.sub <- FindVariableFeatures(prog.sub, 
                                 nfeatures = 1500, 
                                 selection.method = "vst")

prog.sub <- ScaleData(prog.sub, vars.to.regress = "percent.mt")
prog.sub <- RunPCA(prog.sub, npcs = 50)
prog.sub <- RunUMAP(prog.sub, dims = 1:15)
prog.sub <- FindNeighbors(prog.sub)
prog.sub.res.0.13 <- FindClusters(prog.sub, resolution = 0.13)
DimPlot(prog.sub.res.0.13, label = T)

#saveRDS(prog.sub.res.0.13, "prog.sub.res.0.13.rds")

prog.sub.res.0.13 <- JoinLayers(prog.sub.res.0.13)
prog.sub.res.0.13.deg <- FindAllMarkers(prog.sub.res.0.13, 
                                       only.pos = T,
                                       min.pct = 0.1, 
                                       min.diff.pct = 0.1, 
                                       logfc.threshold = 0.1)


prog.sub.res.0.13.deg <- prog.sub.res.0.13.deg %>% 
  mutate(pct.diff=(pct.1-pct.2)) %>% 
  group_by(cluster)  %>% 
  arrange(desc(avg_log2FC), 
          (pct.diff),.by_group = TRUE)

write.csv(prog.sub.res.0.13.deg,"prog.sub.res.0.13.deg.csv")


#DEGs indicated the cluster 6 has marker genes associated with purkinje cells (aldoca, pvalb7,itpr1b)
#Removed cluster 6 from further progenitor analysis.  Cluster is re-labeled PurN_6 in full dataset.
prog.sub.res.0.13.new <- subset(prog.sub.res.0.13, idents= c(0,1,2,3,4,5))


prog.Mgenes <-c("her4.2","her4.3","pcna","mki67","nusap1","tuba8l",
                "fgfrl1a","fgfr3",
                "ccnd1","pax7a","mab21l2","plp1a",
                "sox11b","dlb","neurod4",
                "ccnd2a","atoh1c","neurod1","oprd1b",
                "cx43",
                "ptgdsb.2","ptn",
                "snap25a","sncb")

#ids based off of the processed object used for paper analysis
prog.ids <- c("0"="Prog_1",
              "1"="Prog_2",
              "2"="Prog_3",
              "3"="URL_Prog",
              "4"="RG",
              "5"="Prog/Diff")

prog.sub.res.0.13.rename <- RenameIdents(prog.sub.res.0.13.new, prog.ids)


DotPlot_scCustom(prog.sub.res.0.13.rename,
                 features = prog.Mgenes,
                 dot.scale=8, 
                 dot.min = 0, 
                 x_lab_rotate = T, 
                 scale.min = 0,
                 scale.max = 100,
                 colors_use = palette())


DimPlot(prog.sub.res.0.13.rename, label = T, label.size = 5)

saveRDS(prog.sub.res.0.13.rename, "Prog.res.0.13.rds")


###########################################################################
# RGC and URL
#processed object RGC.URL.res.0.2.rds

RGC.URL.subset <- subset(prog.sub.res.0.13.rename, idents= c("URL_Prog", "RG"))

RGC.URL.subset <- NormalizeData(RGC.URL.subset)
RGC.URL.subset <- FindVariableFeatures(RGC.URL.subset, 
                                       nfeatures = 1500, 
                                       selection.method = "vst")

RGC.URL.subset <- ScaleData(RGC.URL.subset, vars.to.regress = "percent.mt")
RGC.URL.subset <- RunPCA(RGC.URL.subset, npcs = 50)
RGC.URL.subset <- RunUMAP(RGC.URL.subset, dims = 1:15)
RGC.URL.subset <- FindNeighbors(RGC.URL.subset)
RGC.URL.subset.res.0.2 <- FindClusters(RGC.URL.subset, resolution = 0.2)
DimPlot(RGC.URL.subset.res.0.2, label = T)

#saveRDS(RGC.URL.subset.res.0.2, "RGC.URL.subset.res.0.2.rds")


RGC.URL.subset.res.0.2 <- JoinLayers(RGC.URL.subset.res.0.2)
RGC.URL.subset.res.0.2.deg <- FindAllMarkers(RGC.URL.subset.res.0.2, 
                                             only.pos = T,
                                             min.pct = 0.1, 
                                             min.diff.pct = 0.1, 
                                             logfc.threshold = 0.1)


RGC.URL.subset.res.0.2.deg <- RGC.URL.subset.res.0.2.deg %>% 
  mutate(pct.diff=(pct.1-pct.2)) %>% 
  group_by(cluster)  %>% 
  arrange(desc(avg_log2FC), 
          (pct.diff),.by_group = TRUE)

write.csv(RGC.URL.subset.res.0.2.deg,"RGC.URL.subset.res.0.2.deg.csv")

RG.URL.genes <- c("fabp7a", "glula", "slc1a2b", "cx43", 
                  "her4.1","her4.3","her15.1",
                  "insm1a",  "nusap1",
                  "neurod1", "golga7ba", "oprd1b",
                  "stmn1a","pcna", "atoh1c", "atoh1b",
                  "crabp1b", "myo18aa", "gria2b",
                  "ntrk2a", "grm2b", "wnt7bb",
                  "atf3", "cxcl18b","her12",
                  "robo4", "pax7a", "mab21l2",
                  "msrb2","pcp4l1","lgals2a",
                  "spata18","capsla","dcdc2b",
                  "ccl19a.1","cx31.7","foxd3",
                  "eya1","kcnj10a","bmp2a")

#processed object cluster annotations
ids <- c("0"="URL_1",
         "2"="URL_2",
         "4"="URL_3",
         "1"="RG_1",
         "3"="RG_2",
         "5"="RG_3",
         "6"="RG_4",
         "7"="RG_5",
         "8"="RG_6",
         "9"="RG_7")


RGC.URL.subset.res.0.2.rename <- RenameIdents(RGC.URL.subset.res.0.2,ids)


DotPlot_scCustom(RGC.URL.subset.res.0.2.rename,
                 features = RG.URL.genes,
                 dot.scale=8, 
                 dot.min = 0,  
                 x_lab_rotate = T, 
                 scale.min = 0,
                 scale.max = 100,
                 colors_use = palette())


DimPlot(RGC.URL.subset.res.0.2.rename, label=T, label.size = 5) 

saveRDS(RGC.URL.subset.res.0.2.rename, "RGC.URL.res.0.2.rds")


###########################################################################
#Purkinje Neurons
#processed object PurN_6
#cluster 6 from progenitor sub-cluster expresses purkinje markers

PurN <- subset(prog.sub.res.0.13, idents= c(6))

PurN.id <- c("6"="PurN_6")
#PurN_1-5 found in the coarse clustering of the full neur.sub.res.0.3 object.  Clusters 38, 52, 54, 55, 57.

PurN_6 <- RenameIdents(PurN,PurN.id)

saveRDS(PurN_6, "PurN_6.rds")


###########################################################################
#Update full neuron object with sub-cluster annotations.

# Generate a new column called sub_cluster in the metadata
neur.sub.res.3.0.rename$sub_cluster <- as.character(Idents(neur.sub.res.3.0.rename ))

# Add sub-cluster annotations into the new meta data column.

neur.sub.res.3.0.rename$sub_cluster[Cells(hyp.sub.res.0.2.rename)] <- paste(Idents(hyp.sub.res.0.2.rename))

neur.sub.res.3.0.rename$sub_cluster[Cells(opt.sub.res.0.2.rename)] <- paste(Idents(opt.sub.res.0.2.rename))

neur.sub.res.3.0.rename$sub_cluster[Cells(mb.sub.res.0.1.rename)] <- paste(Idents(mb.sub.res.0.1.rename))

neur.sub.res.3.0.rename$sub_cluster[Cells(Pal.res.0.5.rename)] <- paste(Idents(Pal.res.0.5.rename))

neur.sub.res.3.0.rename$sub_cluster[Cells(misc.sub.res.0.3.rename)] <- paste(Idents(misc.sub.res.0.3.rename))

neur.sub.res.3.0.rename$sub_cluster[Cells(SPa.sub.res.0.4.rename)] <- paste(Idents(SPa.sub.res.0.4.rename))

neur.sub.res.3.0.rename$sub_cluster[Cells(vH.sub.res.0.7.rename)] <- paste(Idents(vH.sub.res.0.7.rename))

neur.sub.res.3.0.rename$sub_cluster[Cells(dH.sub.res.0.4.rename)] <- paste(Idents(dH.sub.res.0.4.rename))

neur.sub.res.3.0.rename$sub_cluster[Cells(Tel.res.0.4.rename)] <- paste(Idents(Tel.res.0.4.rename))

neur.sub.res.3.0.rename$sub_cluster[Cells(Di.sub.res.0.2.rename)] <- paste(Idents(Di.sub.res.0.2.rename))

neur.sub.res.3.0.rename$sub_cluster[Cells(HB.sub.res.0.2.rename)] <- paste(Idents(HB.sub.res.0.2.rename))

neur.sub.res.3.0.rename$sub_cluster[Cells(Neur_Diff.res.0.3.rename)] <- paste(Idents(Neur_Diff.res.0.3.rename))

neur.sub.res.3.0.rename$sub_cluster[Cells(GC.TL.res.0.45.rename)] <- paste(Idents(GC.TL.res.0.45.rename))

#subset from OT.sub (cluster9), run after OT
neur.sub.res.3.0.rename $sub_cluster[Cluster9.hb] <- paste("HB_2")

#rename Idents
Idents(object = neur.sub.res.3.0.rename) <- "sub_cluster"

saveRDS(neur.sub.res.3.0.rename, "final.neur.sub.res.3.0.rds")


###########################################################################
# Update full week 3 object with sub-cluster annotations.

# Generate a new column called sub_cluster in the metadata
w3.merge.res.1.0.rename$sub_cluster <- as.character(Idents(w3.merge.res.1.0.rename))

# Add sub-cluster annotations into the new meta data column.

w3.merge.res.1.0.rename$sub_cluster[Cells(neur.sub.res.3.0.rename)] <- paste(Idents(neur.sub.res.3.0.rename))

#GC/TL subset from w3 object
w3.merge.res.1.0.rename$sub_cluster[Cells(GC.TL.res.0.45.rename)] <- paste(Idents(GC.TL.res.0.45.rename))

w3.merge.res.1.0.rename$sub_cluster[Cells(prog.sub.res.0.13.rename)] <- paste(Idents(prog.sub.res.0.13.rename))

w3.merge.res.1.0.rename$sub_cluster[Cells(RGC.URL.subset.res.0.2.rename)] <- paste(Idents(RGC.URL.subset.res.0.2.rename))

w3.merge.res.1.0.rename$sub_cluster[Cells(PurN_6)] <- paste(Idents(PurN_6))


#rename Idents
Idents(object = w3.merge.res.1.0.rename) <- "sub_cluster"

w3.table <- table(w3.merge.res.1.0.rename@meta.data$sub_cluster)
write.csv(w3.table, "upd.w3_cluster.size.csv")          


saveRDS(w3.merge.res.1.0.rename, "final.w3.merge.res.1.0.rds")
