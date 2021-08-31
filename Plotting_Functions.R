

## Heatmap Colour Gradient Scale
library(RColorBrewer)
plotto_color_gradient_blue_red = rev(brewer.pal(n = 10, name = "RdYlBu"))


## Load Arial
library(extrafont)
font_import(paths = "/home/jsch0032/projects/Ethan_iBlastoids/data", pattern = "arial.ttf")

loadfonts()

## Size parameters
stroke = 0.3
size = 0.15
lines = 0.25

## Custom theme for panels with small text
plotto_theme_panel<- function(){
  font <- "Arial"   #assign font family up front

  theme_classic() %+replace%    #replace elements we want to change

    theme(

      #text elements
      plot.title = element_text(             #title
        family = font,            #set font family
        size = 7,                #set font size
        # face = 'bold',            #bold typeface
        hjust = 0),               #raise slightly

      axis.title = element_text(             #axis titles
        family = font,            #font family
        size = 7),               #font size

      axis.text = element_text(              #axis text
        family = font,            #axis famuly
        size = 6),                #font size

      legend.text = element_text(             #legend items
        family = font,            #font family
        size = 7),                #font size

      legend.title = element_text(family = font, size=8),

      line = element_line(size=lines)
    )
}

## Single Cell Plotting Functions
library(Seurat)
library(ggplot2)
plotto_marker_plot = function(m, seurat, reduction = "umap", size=0.15) {
  Embeddings(seurat, reduction = reduction) %>% as_tibble() %>%
    mutate(g = as.matrix(seurat@assays[[seurat@active.assay]][m,])[1,]) -> pdat
  components = names(pdat)[1:2]
  cnames = gsub("_","",components)

  ggplot() +
    geom_point(data = pdat[pdat$g<0,], aes_string(components[1], components[2]), stroke=stroke, size=size, colour="lightgrey") +
    geom_point(data = pdat[pdat$g>=0,], aes_string(components[1], components[2], color="g"),  stroke=stroke, size=size) +
    #    scale_color_gradientn(colours = rev(brewer.pal(n = 7, name =
    # "RdYlBu"))) +
    scale_color_gradientn(colours = c("lightgrey", rev(brewer.pal(n = 11, name =
                                                                    "Spectral")[1:5]))) +
    plotto_theme_panel() +
    theme(legend.position = "bottom", legend.key.height = unit(0.5,"line"),
          legend.spacing.x = unit(0.2, 'cm'),
          legend.box.margin = margin(t=-0.4, unit = "cm")) +
    labs(x=cnames[1], y=cnames[2], colour=m) -> p

  p
}

plotto_signature_scoring_plot = function(sig, seurat, reduction = "umap", size=0.15) {
  # components = paste(ifelse(reduction=="umap", "UMAP", "PC"), 1:2, sep="_")
  # cnames = paste(ifelse(reduction=="umap", "UMAP", "PC"), 1:2, sep="")
  components = colnames(Embeddings(seurat, reduction = reduction))[1:2]
  cnames = gsub("_","",components)
  Embeddings(seurat, reduction = reduction) %>% as_tibble() %>%
    mutate(s = seurat@meta.data[,sig]) %>%
    ggplot(., aes_string(components[1], components[2])) +
    geom_point(aes(colour=s), stroke=stroke, size=size) +
    scale_color_gradientn(colours = rev(brewer.pal(n = 10, name =
                                                     "RdYlBu"))) +
    plotto_theme_panel() +
    theme(legend.position = "bottom", legend.key.height = unit(0.5,"line"),
          legend.spacing.x = unit(0.2, 'cm'),
          line = element_line(size=lines),
          legend.box.margin = margin(t=-0.4, unit = "cm")) +
    labs(x=cnames[1], y=cnames[2], colour=sig) -> p
  p
}

## Arrange Multiple Plots in a grid, setting the panel size
#library(gridExtra)
library(egg)
plotto_panel_it = function(plot_list, width = 2.5, height = 2.5, nrow=1) {
  grid.arrange(grobs = lapply(
    plot_list,
    set_panel_size,
    width = unit(width, "cm"),
    height = unit(height, "cm")
  ),nrow=nrow) -> p
  p
}

