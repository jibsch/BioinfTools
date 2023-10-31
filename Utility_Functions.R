
source("~/tools/JASMINE/JASMINE_V1_11October2021.r")
library(Seurat)
library(tidyverse)
library(escape)


AddSignature <- function(seurat, genes, name, method = "AddModuleScore"){
  if(method == "AddModuleScore") {
    seurat = AddModuleScore(seurat, features = list(x = genes))
    names(seurat@meta.data) = recode(names(seurat@meta.data),
                                     "Cluster1" = name)
  }
  else if(method == "JASMINE_odds") {
    s = JASMINE(seurat@assays[[seu@active.assay]]@data, genes, method = "oddsratio")
    seurat = AddMetaData(seurat, s$JAS_Scores, col.name = name)
  }
  else if(method == "JASMINE_like") {
    s = JASMINE(seurat@assays[[seu@active.assay]]@data, genes, method = "likelihood")
    seurat = AddMetaData(seurat, s$JAS_Scores, col.name = name)
  }
  else if(method == "enrichIt_gsea"){
    s = enrichIt(obj = seurat, method = "ssGSEA", gene.sets = list(x = genes), min.size = NULL)
    seurat = AddMetaData(seurat, metadata = s$x, col.name = name)
  }
  else if(method == "enrichIt_ucell"){
    s = enrichIt(obj = seurat, method = "UCell", gene.sets = list(x = genes), min.size = NULL)
    seurat = AddMetaData(seurat, metadata = s$x, col.name = name)
  }
  seurat
}

AddAllSignatures = function(seurat, genes, name) {
  nop = genes[!genes %in% row.names(seurat)]
  print(paste("The following genes are not in the data and won't be used:", nop))

  seurat = AddSignature(seurat, genes, paste0("Seurat_",name), method = "AddModuleScore")
  seurat = AddSignature(seurat, genes, paste0("JASMINE_oddsratio_",name), method = "JASMINE_odds")
  seurat = AddSignature(seurat, genes, paste0("JASMINE_likelihood_",name), method = "JASMINE_like")
  seurat = AddSignature(seurat, genes, paste0("enrichtIt_ssGSEA_",name), method = "enrichIt_gsea")
  seurat = AddSignature(seurat, genes, paste0("enrichtIt_UCell_",name), method = "enrichIt_ucell")
}

# PlotAll
# ggarrange(plots = lapply(as.list(names(seu@meta.data)[endsWith(names(seu@meta.data), "Mem")]), function(x) plotto_signature_scoring_plot(x, seu, size = 1)))
