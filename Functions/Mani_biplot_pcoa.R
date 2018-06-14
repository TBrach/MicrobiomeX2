biplot_pcoa <- function (physeq, color = "Group", shape = NULL, axis1 = 1, axis2 = 2, 
          show.taxa = TRUE, label = "Genus", repel = FALSE, show.legend = TRUE, 
          custom_palette = NULL) 
{
        dist <- phyloseq::distance(physeq = physeq, method = "bray")
        ordination <- phyloseq::ordinate(physeq, method = "PCoA", 
                                         distance = dist)
        axes = c(axis1, axis2)
        DF <- phyloseq::plot_ordination(physeq, ordination, axes = axes, 
                                        type = "biplot", justDF = TRUE)
        x = colnames(DF)[1]
        y = colnames(DF)[2]
        ord_map = aes_string(x = x, y = y, color = color, shape = shape, 
                             na.rm = TRUE)
        label_map <- aes_string(x = x, y = y, label = label, na.rm = TRUE)
        if (length(extract_eigenvalue(ordination)[axes]) > 0) {
                eigvec = extract_eigenvalue(ordination)
                fracvar = eigvec[axes]/sum(eigvec)
                percvar = round(100 * fracvar, 1)
        }
        DF_Taxa <- DF[DF$id.type == "Taxa", ]
        DF_Samples <- DF[DF$id.type == "Samples", ]
        score <- fracvar[1] * (DF_Taxa[x]^2) + fracvar[2] * (DF_Taxa[y]^2)
        colnames(score) <- NULL
        DF_Taxa <- cbind(DF_Taxa, Score = score)
        rm(score)
        DF_Taxa <- dplyr::arrange(DF_Taxa, desc(Score))
        DF_Taxa <- head(DF_Taxa, 50)
        Tr <- ggplot(DF_Samples, ord_map) + geom_point(na.rm = TRUE, 
                                                       size = 2, show.legend = show.legend) + theme(aspect.ratio = 1)
        color_values <- get_my_palette_colors(custom_palette, color, 
                                              levels(DF_Samples[[color]]), offset = 0)
        Tr <- Tr + scale_color_manual(values = color_values)
        if (show.taxa) {
                if (repel) {
                        Tr <- Tr + geom_text_repel(label_map, data = rm.na.phyloseq(DF_Taxa, 
                                                                                    label), size = 3, color = "black", vjust = "center", 
                                                   hjust = "middle", show.legend = FALSE, na.rm = TRUE)
                }
                else {
                        Tr <- Tr + geom_text(label_map, data = rm.na.phyloseq(DF_Taxa, 
                                                                              label), check_overlap = TRUE, size = 3, color = "black", 
                                             vjust = "center", hjust = "middle", show.legend = FALSE, 
                                             na.rm = TRUE)
                }
        }
        strivar = paste0("PC", axes, "  [", percvar, "%]")
        Tr = Tr + xlab(strivar[1]) + ylab(strivar[2])
        return(Tr)
}