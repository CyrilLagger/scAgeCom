library(ontologyIndex)
library(ontoProc)
library(ontologyPlot)

setwd("C:/Users/Cyril/Desktop/A faire pour Cyril/scDiffCom/Cell_types")

# -------- load cell ontology index and extract list of names
cl = getCellOnto()
cl_names <- cl$name

# -------- load document with cell types
our_cell_types <- read.csv("../data_scAgeCom/analysis/scDiffCom_cell_types.csv")
all_tissues <- unique(our_cell_types$Tissue)
all_cellfamily <- unique(our_cell_types$Family...broad)
all_cellfamily <- all_cellfamily[!(all_cellfamily == "" | grepl("\\+|\\?| or |\\,",all_cellfamily))]
setwd("../data_scAgeCom/analysis/")

# -------- function to plot cell types per tissue and family
plot_celltypes <- function(tissue=all_tissues, family=all_cellfamily, action="save"){
  our_tissue_cell_types <- subset(our_cell_types, Tissue %in% tissue & Family...broad %in% family)
  #our_tissue_cell_types <- unique(c(our_tissue_cell_types$scDiffCom.cell.type,
  #                                  our_tissue_cell_types$TMS.Calico.cell.type))
  # to use only scDiffCom cell types, comment line above and uncomment line below
  our_tissue_cell_types <- unique(our_tissue_cell_types$scDiffCom.cell.type)
  
  tissue_names <- cl_names[cl_names %in% our_tissue_cell_types]
  
  if(action == "show"){
    onto_plot2(cl, names(tissue_names), cex = 0.8)
  } else {
    if(length(tissue) == 1 & length(family) == 1){
      png(file=paste0("Plots/Mixed/",tissue,"_",family,".png"),width=2000, height=1000)
      onto_plot2(cl, names(tissue_names), cex = 0.8)
      dev.off()
    } else if(length(family) == 1 & identical(tissue,all_tissues)) {
      png(file=paste0("Plots/Families/",family,".png"),width=2000, height=1000)
      onto_plot2(cl, names(tissue_names), cex = 0.8)
      dev.off()
    } else if(length(tissue) == 1 & identical(family,all_cellfamily)) {
      png(file=paste0("Plots/Tissues/",tissue,".png"),width=2000, height=1000)
      onto_plot2(cl, names(tissue_names), cex = 0.8)
      dev.off()
    } else {
      onto_plot2(cl, names(tissue_names), cex = 0.8)
      message("Please save plot manually to give it a proper name")
    }
  }
}

# -------- Using function

# Plot and show all cell types for one tissue
plot_celltypes(tissue = "Liver", action = "show")

# Plot and save one cell type for one tissue
plot_celltypes(tissue = "Spleen", family = "leukocyte")

# Plot and save all cell types for each tissue
sapply(all_tissues, function(x){
  plot_celltypes(tissue = x)
}
)
# Plot and save each cell type families for all tissues
sapply(all_cellfamily, function(x){
  plot_celltypes(family = x)
}
)
