#R scripts for processing the Xenium data from HSV-1 infected skin organoids

#Install Voltron from github
#devtools::install_github("Artur-man/VoltRon", force=TRUE)
#install other packages from CRAN or Bioconductor

options(java.parameters = "-Xmx28g")
options(rgl.useNULL=TRUE)
library(VoltRon)

#Packages rJava and RBioFormats are not needed are only needed for initial import
library(rJava)
library(RBioFormats)


library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(igraph)
library(reshape2)
library(DESeq2)
library(apeglm)
library(dendextend)
library(stringr)
library(gridExtra)
library(cowplot)
library(sf)
library(concaveman)
library(data.table)
library(pals)
library(presto)
library(purrr)
library(ggpubr)
library(RColorBrewer)



#The summarySE.R function was taken from http://www.cookbook-r.com/Manipulating_data/Summarizing_data/
#Note that it uses plyr, which can create confusion when calling dplyr functions without dplyr:: 
source("./summarySE.R")


#From https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#Based on function above and https://jackw01.github.io/HCLPicker/
gg_color_lum <- function(h, c, n) {
  lums = seq(from=40, 100, length = n + 1)
  rev(hcl(h = h, l = lums, c = c)[1:n])
}


#Colors for cell types
celltype_cols <- c(`Adipocytes/Sebaceous Glands` = "#F8766D", `ciliated epithelial cells` = "#F27C56", `Arrector pili muscle` = "#EC823A",
                   `Endothelial cells` = "#DC8D00", `Melanocytes` = "#D89000", `MerkelCells` = "#CE9500", `Neuronal cells_1` = "#C39A00", `Neuronal cells_2` = "#B79F00",
                   `SchwannCells` = "#AAA300", `Chondrocytes` = "#93AA00",
                   `Fibroblasts: Dermal Papilla` = "#007700", `Fibroblasts: Myofibroblasts_1` = "#00BA38", `Fibroblasts: Mesenchymal_1` = "#55FF8A",
                   `Fibroblasts: Mesenchymal_2` = "#00807E", `Fibroblasts: Pro-inflammatory` = "#00C0BD", `Fibroblasts: Dermal Sheath (cyc)` = "#00FFE4", 
                   `Fibroblasts: Papillary_1` = "#0078AF", `Fibroblasts: Papillary_2` = "#00B5EE", `Fibroblasts: Papillary_3` = "#00F4FF",
                   `Fibroblasts: Dermal Sheath` = "#0957CE", `Fibroblasts: undefined` = "#BFDDFF",
                   `Keratinocytes: Stratum Basale` = "#A600BB", `Keratinocytes: Stratum corneum` = "#E16FF7", `Keratinocytes: IRS` = "#FFC7FF",
                   `Keratinocytes: Matrix` = "#BB0092", `Keratinocytes: Bulge region` = "#FD61D1", `Keratinocytes: ORS/companion layer` = "#FFB5CF",
                   `Keratinocytes: Stratum spinosum/Granulosum` = "#B42049", 
                   Undefined1 ="gray90", 
                   Undefined2 ="gray88", Undefined3 ="gray86", Undefined4 ="gray84", Undefined5 ="gray82", Undefined6 ="gray80", 
                   Undefined7 ="gray78", Undefined8 ="gray76", Undefined9 ="gray74", Undefined10 ="gray72", Undefined11 ="gray70", UndefinedCorePap="gray30")

time_levels <- c("0dpi", "2dpi", "3dpi", "4dpi")

#This function is occasionally needed if errors occur with VoltRon objects
fixVoltRon <- function(object){
  
  # sample.metadata
  sample.metadata <- SampleMetadata(object)
  
  # fix samples and layers
  for(samp in unique(sample.metadata$Sample)){
    
    # sample
    object_sample <- object[[samp]]
    
    # correct
    catch_connect <- try(slot(object_sample, name = "zlocation"), silent = TRUE)
    if(methods::is(catch_connect, 'try-error') || methods::is(catch_connect,'error')){
      object_sample@zlocation <- numeric(0)
    }
    catch_connect <- try(slot(object_sample, name = "adjacency"), silent = TRUE)
    if(methods::is(catch_connect, 'try-error') || methods::is(catch_connect,'error')){
      object_sample@adjacency <- matrix()
    }
    
    # fix layers
    for(lyr in unique(sample.metadata$Layer[sample.metadata$Sample == samp])){
      object_layer <- object_sample[[lyr]]
      
      # correct
      catch_connect <- try(slot(object_layer, name = "connectivity"), silent = TRUE)
      if(methods::is(catch_connect, 'try-error') || methods::is(catch_connect,'error')){
        object_layer@connectivity <- igraph::make_empty_graph()
      }
      
      object_sample[[lyr]] <- object_layer
    }
    object[[samp]] <- object_sample
  }
  object
}


#Import Xenium data folder (can be obtained from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313919)
#Note that this step is memory intense. Consider increasing the resolution_level if it does not work
skin_all <- importXenium("/Volumes/Storage/MoreApplications/stuff/Xenium/Xenium_Stockholm_Second_Run/TC_MDC/catalyst_release_TC_MDC_Dec23/Skin_organoids/", resolution_level = 2, overwrite_resolution = TRUE)
saveRDS(skin_all, "skin_all.rds")

#Filter out low quality cells
skin_all$HSV1Count <- colSums(vrData(skin_all, norm = FALSE)[c("HSV1_LAT", "HSV1_UL27", "HSV1_UL29", "HSV1_UL54", "HSV1_US1"),])
skin_all$NonVirusCount <- skin_all$Count - skin_all$HSV1Count
skin_all <- subset(skin_all, NonVirusCount > 10)

#Annotate the 12 individual organoids. In the interface, use simple labels (R1, R2, R3, R21 etc.) that are then expanded
#The data can also be transferred from the fully annotated object at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313919 , simply copy the annotation column from the @metadata@cell dataframe
skin_all_ann <- annotateSpatialData(skin_all, label="annotation", pt.size=0.1)
skin_all_ann@metadata@cell <- skin_all_ann@metadata@cell %>% mutate(annotation = case_when(annotation == "R1" ~ "region_0dpi_rep1", 
                                                                                            annotation == "R2" ~ "region_0dpi_rep2", 
                                                                                            annotation == "R3" ~ "region_0dpi_rep3", 
                                                                                            annotation == "R21" ~ "region_2dpi_rep1", 
                                                                                            annotation == "R22" ~ "region_2dpi_rep2", 
                                                                                            annotation == "R23" ~ "region_2dpi_rep3",
                                                                                            annotation == "R31" ~ "region_3dpi_rep1", 
                                                                                            annotation == "R32" ~ "region_3dpi_rep2", 
                                                                                            annotation == "R33" ~ "region_3dpi_rep3",
                                                                                            annotation == "R41" ~ "region_4dpi_rep1", 
                                                                                            annotation == "R42" ~ "region_4dpi_rep2", 
                                                                                            annotation == "R43" ~ "region_4dpi_rep3",
                                                                                            annotation == "undefined" ~ "undefined"))

#Add timepoint column
skin_all_ann@metadata@cell$timepoint <- gsub("^[^\\_]*_([^\\_]*)_[^\\_]*$", "\\1", skin_all_ann@metadata@cell$annotation, perl=TRUE)

#Plot and check
vrSpatialPlot(skin_all_ann, group.by = "annotation", alpha = 1, plot.segments = TRUE, background.color = "black", assay="Assay1")
vrSpatialPlot(skin_all_ann, group.by = "timepoint", alpha = 1, plot.segments = TRUE, background.color = "black", assay="Assay1")

#To strengthen the annotation, include the three uninfected organoids from the Merkel cell polyomavirus experiment (Albertini et al.)
pyv_0dpi <- importXenium("/Volumes/Storage/MoreApplications/stuff/Xenium/20240418__105901__ST004_X0057_X0058_Landthaler/output-XETG00046__0021832__Region_4__20240418__110821", resolution_level = 2, overwrite_resolution = TRUE)


#Temporary
pyv_0dpi <- readRDS("pyv_0dpi.rds")
saveRDS(skin_all_ann, "skin_all_ann_temp.rds")
