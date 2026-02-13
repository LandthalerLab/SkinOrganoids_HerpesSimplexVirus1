options(rgl.useNULL=TRUE)

#Install Voltron from github
#devtools::install_github("Artur-man/VoltRon", force=TRUE)
#install other packages from CRAN or Bioconductor
#Packages rJava and RBioFormats are not needed are only needed for initial import

library(Seurat)
library(dplyr)
library(ggplot2)
library(VoltRon)
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
#source("~/temp/transfer/notes and scripts/summarySE.R")
source("~/Documents/Largescale-data/notes and scripts/summarySE.R")
library(sf)
library(concaveman)
library(data.table)
library(pals)
library(presto)
library(purrr)
library(ggpubr)
library(RColorBrewer)

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
