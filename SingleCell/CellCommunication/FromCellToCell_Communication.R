#' Generate chemokine/cytokine DotPlots for epithelial/stroma and immune cells
#'
#' @param seurat_obj A Seurat object containing your single-cell data.
#' @param CellchatInteteracome The detailed cell-cell communications among all the cell types within different timepoints
#' @param Cellchatpathways The detailed cell-cell pathways among all the cell types within different timepoints note this is a wraped up prob; three columns: pathway_name, total_prob and date
#' @param requested_pathway The requested pathway
#' @param requested_pathway_tag The requested pathway tag
#' 
#' @param ligand_celltypes_tag Character value for the epithelial/stroma metadata tag.
#' @param ligand_celltypes Character vector of identities to keep in epithelial plot.
#' 
#' @param receptor_celltypes_tag Character vector of immune cell metadata tag
#' @param receptor_celltypes Character vector of receptor identities.
#' @param output_dir Folder to save plots (default is current).
#'
#' @return Saves two dot plots to PDF.
#' 
Interestedpathways_celltocell_dotplots <- function(seurat_obj,
                                                   CellchatInteteracome,
                                                   Cellchatpathways,
                                                   requested_pathway_tag,
                                                   requested_pathway,
                                                   ligand_celltypes_tag,
                                                   ligand_celltypes,
                                                   receptor_celltypes_tag,
                                                   receptor_celltypes,
                                                   output_dir) {
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))  # set up the color
  Requestedpathwaystrength <-
    Cellchatpathways %>% filter(pathway_name %in% requested_pathway)
  Requestedpathwaystrength$date <-
    factor(Requestedpathwaystrength$date,
           levels = c("D0", "D3", "D7", "D28"))
  Requestedpathwaystrength$pathway_name <-
    factor(Requestedpathwaystrength$pathway_name, levels = requested_pathway)
  # Generate a heatmap for the pathway strength
  RequestedStrengthOverTimeFigure <-
    ggplot(Requestedpathwaystrength,
           aes(
             x = date,
             y = pathway_name,
             fill = log2(total_prob)
           )) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(name = "Total Prob", option = "D") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      #strip.text = element_text(face = "bold"),
      #panel.grid = element_blank()) +
      labs(title = "Pathway Activity Over Time", x = "Date", y = requested_pathway_tag)
    )
  StreHeight = max(5, length(requested_pathway)/2)
  ggsave(
    filename = paste0(output_dir, requested_pathway_tag, "_stengthovertime.pdf"),
    plot = RequestedStrengthOverTimeFigure,
    height = StreHeight,
    width = 5
  )
  # Subset ligands and receptor for the requested pathways from requested ligand and receptor cell types
  CellCommunicationSubtypes <-
    CellchatInteteracomedetails_named_combined %>% filter(
      pathway_name %in% requested_pathway &
        source %in% ligand_celltypes &
        target %in% receptor_celltypes
    ) %>% mutate (SenderReceiverCell = paste0(source, " -> ", target)) %>% arrange(desc(prob))
  # set up the levels of date
  CellCommunicationSubtypes$date <-
    factor(CellCommunicationSubtypes$date,
           levels = c("D0", "D3", "D7", "D28"))
  # generate the heatmaps based on the ligand and receptors
  RequestedInteractionBetweenLigandReceptor <-
    ggplot(
      CellCommunicationSubtypes,
      aes(
        x = SenderReceiverCell,
        y = interaction_name_2,
        fill = log2(prob)
      )
    ) +
    geom_dotplot(dotsize = 0.3) +
    #geom_tile(color = "white") +
    ggplot2::scale_fill_gradientn(colours = myPalette(4)) +
    facet_grid(rows = vars(date),
               scales = "free_y",
               space = "free_y") +   # separate by date
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          #strip.text = element_text(face = "bold"),
          #panel.grid = element_blank()) +
          labs(
            x = paste0("From ", ligand_celltypes_tag, " to ", receptor_celltypes_tag),
            y = paste0(requested_pathway_tag, " Interaction")
          ))
  # set up the heights and width
  InteractNumber<- CellCommunicationSubtypes %>% mutate(interactomewithintime= paste0(interaction_name,"date")) %>% pull(interactomewithintime) %>%  length()
  Interheight <- min(30,max(6, InteractNumber/2))
  Interwidth <- min(40,max(12, nrow(CellCommunicationSubtypes)/2))
  ggsave(
    paste0(
      output_dir,
      requested_pathway_tag,
      "From",
      ligand_celltypes_tag,
      "To",
      receptor_celltypes_tag ,
      "_Interactome_Dotplot_Withintime.pdf"
    ),
    plot = RequestedInteractionBetweenLigandReceptor,
    #height = Interheight,
    #width = Interwidth
    height = 7,
    width = 20
    
  )
  # get the ligand features
  ligand_features <-
    CellCommunicationSubtypes$ligand %>% unique()
  # get the receptor features
  receptor_features <-
    CellCommunicationSubtypes$receptor %>% unique()
  # generate the dotplots for these ligand and receptors for the requested ligand cell types and receptor cell types
  ligand_obj <-
    subset(seurat_obj, idents = ligand_celltypes)
  Idents(ligand_obj) <-
    factor(Idents(ligand_obj), levels = ligand_celltypes)
  # Subset immune
  receptor_obj <-
    subset(seurat_obj, idents = receptor_celltypes)
  Idents(receptor_obj) <-
    factor(Idents(receptor_obj), levels = receptor_celltypes)
  # Epithelial ligand DotPlot
  ligand_plot <-
    DotPlot(
      ligand_obj,
      features = rev(ligand_features),
      cols = c("lightgrey", "blue")
    ) +
    coord_flip() +
    scale_x_discrete(position = "top") +
    theme_classic() +
    theme(
      panel.border = element_rect(
        color = "black",
        fill = NA,
        size = 1
      ),
      axis.line = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.4, "cm"),
      axis.text.x = element_text(
        angle = 60,
        hjust = 1,
        size = 8
      )
    ) + guides(color = guide_colorbar(title = "AverExp"),
               size = guide_legend(title = "PercExp"))
  # receptor DotPlot
  receptor_plot <-
    DotPlot(receptor_obj,features = rev(receptor_features),cols = c("lightgrey", "red")) +
    coord_flip() +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.line = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.4, "cm"),
      axis.text.x = element_text(
        angle = 60,
        hjust = 1,
        size = 8
      )
    ) +  guides(color = guide_colorbar(title = "AverExp"), size = guide_legend(title = "PercExp"))
  combined_plot <- ligand_plot + receptor_plot + plot_layout(ncol = 2,widths = c(1, 3))
  #combined_plot <- ligand_plot + receptor_plot + plot_layout(ncol = 2)
  # set up the heights and width
  LRheight <- max(4, length(ligand_features)/2)
  LRwidth <- max(8, (length(ligand_celltypes)+ length(receptor_celltypes))/2)
  ggsave(
    filename = paste0(
      output_dir,
      ligand_celltypes_tag,
      "_",
      receptor_celltypes_tag,
      requested_pathway_tag,
      "CellComminication_DotPlot.pdf"
    ),
    plot = combined_plot,
    #height = LRheight,
    height = 6,
    #width = LRwidth
    width = 14
  )
}