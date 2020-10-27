library(ontoProc)
library(ontologyIndex)
library(ontologyPlot)
library(data.table)

cl = getCellOnto()

cl_names <- cl$name
head(cl_names)

our_cell_types <- read.csv("../../../../../cell_types_with_counts.csv")


celltype <- function(tissue){
our_tissue_cell_types <- subset(our_cell_types, Tissue == tissue, select=Cell.type)
our_tissue_cell_types <- our_tissue_cell_types$Cell.type
our_tissue_cell_types <- sub(".*_","", our_tissue_cell_types)
our_tissue_cell_types <- unique(our_tissue_cell_types)
tissue_names <- cl_names[cl_names %in% our_tissue_cell_types]
onto_plot2(cl, names(tissue_names), cex = 0.8)
}

celltype_bis <- function(tissue){
  our_tissue_cell_types <- subset(our_cell_types, Tissue == tissue, select=Cell.type)
  our_tissue_cell_types <- our_tissue_cell_types$Cell.type
  our_tissue_cell_types <- sub(".*_","", our_tissue_cell_types)
  our_tissue_cell_types <- unique(our_tissue_cell_types)
  tissue_names <- cl_names[cl_names %in% our_tissue_cell_types]
}

celltype("Thymus")

celltype("Marrow")
onto_plot2(cl, c(names(immune_names), "CL:0000232"), cex = 0.8)

test <- lapply(c("Marrow", "Lung"), celltype)


data(hpo)
remove_uninformative_terms(hpo, list(Patient1=c("HP:0001873","HP:0000118")))
remove_terms_with_less_than_n_occurrences(hpo, 
                                          term_sets=list("HP:0001873", "HP:0001902"), n=2)
cl

?minimal_set
minimal_set(hpo, c("HP:0001873", "HP:0001872"))

subset(our_cell_types, Tissue == "Aorta", select=Cell.type)$Cell.type


onto_plot(cl, terms=remove_links(cl,names(celltype_bis("Kidney"))))

onto_plot2(cl, c(names(celltype_bis("Kidney")), "CL:0000653", "CL:0002188") )

onto_plot2(cl, c(names(celltype_bis("Lung")), "CL:0000158") )

onto_plot2(cl, c(names(celltype_bis("Marrow")), "CL:0000038") )



celltype("GAT")



###############################




##
onto_plot2(cl, names(cl_names[cl_names %in% our_cell_types[Tissue == "Aorta"]$cell_type_short]))


# epithelial cells
our_epithelial_cell_types <- subset(our_cell_types, Cell.type.subfamily == "epithelial cell" | Cell.type.subfamily == "endothelial cells", select=Cell.type)
our_epithelial_cell_types <- our_epithelial_cell_types$Cell.type
our_epithelial_cell_types <- sub(".*_","", our_epithelial_cell_types)
our_epithelial_cell_types <- unique(our_epithelial_cell_types)

epithelial_names <- cl_names[cl_names %in% our_epithelial_cell_types]

onto_plot2(cl, names(epithelial_names), cex = 0.8)

# immune cells
our_immune_cell_types <- subset(our_cell_types, Cell.type.subfamily == "leukocyte" | Cell.type.subfamily == "lymphocyte" | Cell.type.subfamily == "myeloid cell" | Cell.type.subfamily == "professional antigen presenting cell", select=Cell.type)
our_immune_cell_types <- our_immune_cell_types$Cell.type
our_immune_cell_types <- sub(".*_","", our_immune_cell_types)
our_immune_cell_types <- unique(our_immune_cell_types)

immune_names <- c(cl_names[cl_names %in% our_immune_cell_types])



# connective tissue cells
our_connective_cell_types <- subset(our_cell_types, Cell.type.subfamily == "connective tissue cell", select=Cell.type)
our_connective_cell_types <- our_connective_cell_types$Cell.type
our_connective_cell_types <- sub(".*_","", our_connective_cell_types)
our_connective_cell_types <- unique(our_connective_cell_types)

connective_names <- cl_names[cl_names %in% our_connective_cell_types]

onto_plot2(cl, names(connective_names), cex = 0.8)

# endothelial cells
our_endothelial_cell_types <- subset(our_cell_types, Cell.type.subfamily == "endothelial cell", select=Cell.type)
our_endothelial_cell_types <- our_endothelial_cell_types$Cell.type
our_endothelial_cell_types <- sub(".*_","", our_endothelial_cell_types)
our_endothelial_cell_types <- unique(our_endothelial_cell_types)

endothelial_names <- cl_names[cl_names %in% our_endothelial_cell_types]

onto_plot2(cl, names(endothelial_names), cex = 0.8)
#onto_plot(cl, names(endothelial_names))
# onto_plot(cl, terms=remove_links(cl, get_ancestors(cl, names(endothelial_names))))


# neural cells
our_neural_cell_types <- subset(our_cell_types, Cell.type.subfamily == "neuron" | Cell.type.subfamily == "macroglial cell" | Cell.type.subfamily == "glial cell", select=Cell.type)
our_neural_cell_types <- our_neural_cell_types$Cell.type
our_neural_cell_types <- sub(".*_","", our_neural_cell_types)
our_neural_cell_types <- unique(our_neural_cell_types)

neural_names <- cl_names[cl_names %in% our_neural_cell_types]

onto_plot2(cl, names(neural_names), cex = 0.8)

# stem cells and precursor cells