# - load packages -
library(phyloseq); # packageVersion("phyloseq")
library(dada2); # packageVersion("dada2")
library(vegan); # packageVersion("vegan")
library(ggplot2); # packageVersion("ggplot2")
library(dplyr); # packageVersion("dplyr")
library(tidyr); # packageVersion("tidyr")
library(gridExtra); # packageVersion("gridExtra")
library(xtable); # packageVersion("xtable")
library(RVAideMemoire); #packageVersion("RVAideMemoire")
library(viridis); # packageVersion("viridis")
library(scales); # packageVersion("scales") # for the oob = squish option in gradient plots
library(ggthemes); # packageVersion("ggthemes")
library(DESeq2); # packageVersion("DESeq2")
library(ggpubr); # packageVersion("ggpubr")
library(RColorBrewer)
library(pheatmap)
#library(grid)
#library(breakaway)
#library(bookdown)
# --


# - source functions -
functionpath <- "."
source(file.path(functionpath, "_n_000_helper_functions.R"))

source(file.path(functionpath, "_n_050_diff_abundance_functions.R"))
# --



# - load your phyloseq object - 
name_phyloseq_rds <- "physeq_microdiab_ngt_Men.rds"

datapath <- "."


ps <- readRDS(file.path(datapath, name_phyloseq_rds))
# --


# - Define the group variable for sample comparisons -
group_var <- "Country" # MANDATORY: a variable in sample_data(ps) based on which samples will be compared

group_var_levels <- c("IN", "DK") # defines the order of the groups in all plots. If set to NULL:

color_levels <- c(cbPalette[2], cbPalette[4]) # choose your preferred colors for each group in your group_var. If set to NULL:

names(color_levels) <- group_var_levels

shape <- "Gender" 
# -- 



# - do a simple prevalence filtering to speed up the example -
prevalence <- 20
min_obs <- 0L
ps_filt <- phyloseq::filter_taxa(ps, function(x){(sum(x > min_obs) > (prevalence/100)*length(x))}, prune = TRUE)
# --


# - determine the phylum colors -
PhylaForColor <- check_phyla_distribution(ps_filt)

if (nrow(PhylaForColor) < 16) {
        phylum_colors <- make_color_vector(as.character(PhylaForColor$Phylum), QuantColors15)
} else {
        phylum_colors <- make_color_vector(as.character(PhylaForColor$Phylum), c(QuantColors15, viridis(nrow(PhylaForColor)-15)))
}
# --




# - do a fisher prevalence test -
physeq_to_test <- ps_filt

diff_ab_df <- test_diffs_in_prevalence_single(physeq = physeq_to_test, group_var = group_var, compare = group_var_levels, p.adj.method = "fdr", minCount = 0L, symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))
# --


# - catch the hits -
hit_list <- format_hit_table(diff_ab_df, p.adjust.threshold = 0.05, p.adjust.method = NULL)
taxa_hit_df <- hit_list[["hit_table"]]
# --


# - define colors for heatmap -
significance_colors <- brewer.pal(4, "Reds")
significance_colors <- c(rev(significance_colors), "gray", "violet")
names(significance_colors) = c("****", "***", "**", "*", "ns", "?")
taxa_colors <- list("signi_adj" = significance_colors, "Phylum" = phylum_colors)

if (!is.null(shape)){
        sample_colors <- list(color_levels, NA)
        names(sample_colors) <- c(group_var, shape)
} else {
        sample_colors <- list(color_levels)
        names(sample_colors) <- group_var
}
# --

# - get a more informative taxa annotation for the heatmap -
taxa_annotation <- taxa_hit_df$Annotation
# --


# - do the heatmap -
p <- plot_heatmap_physeq(physeq_to_test, sample_colors = sample_colors, taxa_info_df = head(taxa_hit_df, 40), taxa_colors = taxa_colors, 
                         taxa_annotation = head(taxa_annotation, 40), max_abundance_for_color = NULL, gradient_steps = c(0.15, 0.3, 0.45, 1), 
                         zero_color = "gray", color_function = viridis, color_steps_bw_markers = 10, log_transform = FALSE, drop_color_levels = TRUE,
                         border_color = NA, 
                         cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = TRUE, show_colnames = FALSE, annotation_names_row = FALSE, 
                         annotation_names_col = FALSE, annotation_legend = TRUE, legend = TRUE, font_size = 14, 
                         fontsize_row = 12, fontsize_col = 12, fontsize_number = 12)


# IF YOU WANT TO SAVE IT:

# pdf(file = "heatmapprev.pdf", width = 9, height = 6)
# grid::grid.newpage()
# grid::grid.draw(p$gtable)
# dev.off()
# --


