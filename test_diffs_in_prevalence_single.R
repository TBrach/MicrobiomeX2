# --
#######################################
### FUNCTION: test_diffs_in_prevalence_single
#######################################
# Function performs fisher exact test on prevalence (absence presence) between the levels
# in a grouping factor
# INPUT:
# physeq: phyloseq
# group_var: name of the column in sample_data(physeq) that defines the groups
# p.adj.method, used in p.adjust
# minCount: present are taxa in species with more counts than minCount
# OUTPUT:
# list of data.frames, one data frame for each combi of levels in your grouping factor
# The data frames are ordered by p_value, and the tax_table has been cbound:)

test_diffs_in_prevalence_single <- function(physeq, group_var, compare = NULL, p.adj.method = "fdr", minCount = 0L, symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) {
        
        
        if (taxa_are_rows(physeq)) { 
                physeq <- t(physeq)
        }
        
        if(! group_var %in% colnames(sample_data(physeq))) {
                stop("The given group_var is not a variable in the sample data of the loaded phyloseq object.")
        }
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        if (!is.null(compare)) {
                group_var_levels <- compare
        } else {
                group_var_levels <- levels(group_fac)
        }
        
        if (length(group_var_levels) != 2) {
                stop(paste0("compare (group_var_levels) must consist of two groups - you asked for ", 
                            paste(group_var_levels, collapse = ", ")))
        }
        
        if (!all(group_var_levels %in% levels(group_fac))) {
                stop("Not all given compare (group_var_levels) are actually levels in group_var column.")
        }
        
        
        
        CT <- as(otu_table(physeq), "matrix")
        CT <- CT > minCount
        
        
        prev_list <- lapply(group_var_levels, function(level){
                data.frame(Present = colSums(CT[group_fac == level, ]),
                           Absent = colSums(!CT[group_fac == level, ]))
        })
        
        pr_ab_gr1 <- prev_list[[1]] #pr_ab = presence absence
        pr_ab_gr2 <- prev_list[[2]]
        
        rowwise_compare_matrix <- cbind(pr_ab_gr1, pr_ab_gr2)
        
        FisherTests <- lapply(1:nrow(rowwise_compare_matrix), function(e){
                mat_fisher <- matrix(c(rowwise_compare_matrix[e, 1],
                                       rowwise_compare_matrix[e, 3],
                                       rowwise_compare_matrix[e, 2],
                                       rowwise_compare_matrix[e, 4]), ncol = 2)
                fisher.test(mat_fisher, conf.int = TRUE, conf.level = 0.95)
                # fisher.test(mat_fisher, conf.int = FALSE)
        })
        
        # - take out oddsRatios -
        oddRatios <- sapply(FisherTests, function(TestResult){TestResult$estimate})
        oddRatios_lb <- sapply(FisherTests, function(TestResult){TestResult$conf.int[1]})
        oddRatios_ub <- sapply(FisherTests, function(TestResult){TestResult$conf.int[2]})

        direction <- rep(group_var_levels[1], length(oddRatios))
        direction[oddRatios < 1] <- group_var_levels[2]
        
        oddRatios_lb[oddRatios < 1] <- 1/oddRatios_lb[oddRatios < 1]
        oddRatios_ub[oddRatios < 1] <- 1/oddRatios_ub[oddRatios < 1]
        oddRatios_lb_final <- pmin(oddRatios_lb, oddRatios_ub)
        oddRatios_ub_final <- pmax(oddRatios_lb, oddRatios_ub)
        oddRatios[oddRatios < 1] <- 1/oddRatios[oddRatios < 1]
        # --
        
        
        # - take out p_vals and assign significance levels -
        p_vals <- sapply(FisherTests, function(TestResult){TestResult$p.value})
        p_vals_adj <- p.adjust(p_vals, p.adj.method)
        
        symnum.args$x <- p_vals
        significance <- do.call(stats::symnum, symnum.args) %>% as.character()
        symnum.args$x <- p_vals_adj
        significance_adj <- do.call(stats::symnum, symnum.args) %>% as.character()
        # --
        
        # - take out prevalence percentages -
        prev_PC_gr1 <- round(100*(pr_ab_gr1[, "Present"]/sum(group_fac == group_var_levels[1])), 1)
        prev_PC_gr2 <- round(100*(pr_ab_gr2[, "Present"]/sum(group_fac == group_var_levels[2])), 1)
        # --
        
        df <- data.frame(p_val = p_vals, p_val_adj = p_vals_adj, signi = significance, signi_adj = significance_adj, oddsRatio = round(oddRatios, 2), oddsRatio_lb = round(oddRatios_lb_final, 2),
                         oddsRatio_ub = round(oddRatios_ub_final, 2), direction = direction, prev_PC_gr1 = prev_PC_gr1,  prev_PC_gr2 = prev_PC_gr2)

        df <- cbind(as.data.frame(df), tax_table(physeq))
        df$Taxon <- colnames(CT)
        df <- arrange(df, p_val) %>% select(Taxon, 1:(ncol(df)-1))
        df
        
}
# --

format_hit_table <- function (result_df, p.adjust.threshold = 0.1, p.adjust.method = NULL) {
        
        if (!all(c("p_val", "p_val_adj", "Taxon", "direction", "signi", "signi_adj") %in% colnames(result_df))) {
                stop("result_df should be a data.frame generated by one of the differential abundance tests and contain all corresponding columns.")
        }

        if (!is.null(p.adjust.method)) {
                result_df$p_val_adj = p.adjust(result_df$p_val, method = p.adjust.method)
        }
        
        
        result_df <- arrange(result_df, p_val_adj, p_val)
        
        no_hits <- sum(result_df$p_val_adj <= p.adjust.threshold)
        
        keepTaxa <- no_hits
        
        if (keepTaxa < 5 && nrow(result_df) >= 5) {
                keepTaxa <- 5
        } else if (keepTaxa < 5 && nrow(result_df < 5)){
                keepTaxa <- nrow(result_df)
        }
        
        
        df <- result_df[1:keepTaxa,]
        
        df$Annotation <- get_taxon_names(df)
        
        df <- select(df, Taxon, Annotation, p_val, p_val_adj, signi, signi_adj, direction, colnames(df)[!(colnames(df) %in% c("Taxon", "Annotation", "p_val", "p_val_adj", "signi", "signi_adj", "direction"))])
        
        rownames(df) <- df$Taxon
        
        list(hit_table = df, no_hits = no_hits)
}



# --
#######################################
### FUNCTION: get_taxon_names
#######################################


get_taxon_names <- function(df) {
        
        df1 <- df[, colnames(df) %in% c("Kingdom", "Phylum", "Class", "Order", "Family", 
                                       "Genus", "Species")]
        
        if (ncol(df1) == 0) {stop("the provided data frame did not contain the expected taxonomy columns such as Phylum, Class etc.")}
        
        df1[] <- lapply(df1, as.character)
        
        Last_NotNA_Position <- apply(df1, 1, function(x){length(which(!is.na(x)))})
        Last_NotNA_Position[Last_NotNA_Position == 0] <- 1
        
        Names <- vector(mode = "character", length = nrow(df1))
        
        for (i in 1:nrow(df1)) {
                if (Last_NotNA_Position[i] == 7){
                        Names[i] <- paste(df1[i, "Genus"], df1[i, "Species"], sep = " ")
                } else {
                        Names[i] <- df1[i, Last_NotNA_Position[i]] 
                }
        }
        
        Names[is.na(Names)] <- "NA"
        Names
        
}
#--







#########################
## plot_heatmap_physeq ##
#########################
# function that uses pheatmap to draw a heatmap of the physeq object, or a pruned version of it determined by taxa_info_df
# INPUTs
# physeq = physeq object
# sample_colors = named list of named color factors. The sample_data(physeq) data.frame is used as sample_info_df to color the samples (columns) in the plot.
#    Only those variables/columns in sample_info_data that are names in sample_colors are considered for coloring the samples. I.e. if sample_colors = NULL no coloring occurs.
     # NB: if a variable is a name in sample_colors but the color_factor is NULL, default colors will be used
# taxa_info_df = a data frame similar to sample_info_df (see sample_colors) that together with taxa_colors will be used to label the taxa.
     # NB: MORE importantly, it is also used to restrict the included taxa and to order the taxa!! The rownames(taxa_info_df) must be the taxa_names in physeq 
     # of the taxa you want to include in the heatmap in the order you want the taxa to be shown!
# taxa_colors: see sample_colors

hm.parameters <- list(DF_CT, color = myColors, breaks = myBreaks, border_color = NA, cluster_cols = FALSE, cluster_rows = FALSE,
                      show_rownames = (show_rownames && nrow(DF_CT) < 70), show_colnames = FALSE, annotation_col = sam_info,
                      annotation_row = taxa_data, annotation_colors = annotation_colors, labels_row = taxa_annotation, annotation_names_row = FALSE,
                      annotation_names_col = FALSE, annotation_legend = TRUE, legend = TRUE, font_size = 10, fontsize_row = 8, fontsize_col = 8, fontsize_num = 6)

plot_heatmap_physeq <- function (physeq, sample_colors = NULL, taxa_info_data, taxa_colors, taxa_annotation = NULL, 
          max_abundance_for_color = NULL, gradient_steps = c(0.15, 0.3, 0.45, 1), taxa_colors = NULL, block = NULL, block_colors = NULL,
          log_transform = FALSE, border_color = "grey60", filename = NA, show_rownames = TRUE, cellwidth = NA, cellheight = NA, gaps_col = NULL, 
          custom_palette = NULL) {
        
        # - keep only taxa_to_include in physeq -
        if (!is.null(taxa_to_include)) {
                
                if (!all(taxa_to_include %in% taxa_names(physeq))) {
                        stop("not all taxa in taxa_to_include are taxa in physeq!")
                }
                pruned_ps <- prune_taxa(taxa_to_include, physeq)
        } else {
                taxa_to_include <- taxa_names(physeq)
                pruned_ps <- physeq
        }
        # --
        
        
        
        # - test or set group_var, compare, and sample_colors input -
        sam_info <- as(sample_data(pruned_ps), "data.frame")
        
        if(! group_var %in% colnames(sam_info)) {
                stop("The given group_var is not a variable in the sample data of physeq.")
        }
        
        
        if (is.null(compare)) { 
                compare <- unique(sam_info[[group_var]])
        }
        
        if (!all(compare %in% unique(sam_info[[group_var]]))) {
                stop("Not all given compare groups are actually levels in group_var column.")
        }
        
        if (is.null(sample_colors)) { 
                if (length(compare) < 9) {
                        sample_colors <- cbPalette[1:length(compare)]
                } else if (length(compare) < 16) {
                        sample_colors <- QuantColors15[1:length(compare)]
                } else {
                        sample_colors <- viridis(length(compare))
                }
        }
        
        if (length(sample_colors) != length(compare)) {
                stop("The number of colors in the given sample_colors did not fit to the number of groups in compare.")
        }
        
        if (!all(areColors(sample_colors))) {
                stop("Stopped because not all entries in sample_colors were valid R colors.")
        }
        
        names(sample_colors) <- compare
        # --
        
        # - prune samples to only keep those covered by compare -
        keepSamples <- sample_names(pruned_ps)[sam_info[[group_var]] %in% compare]
        
        pruned_ps <- prune_samples(keepSamples, pruned_ps)
        # --
        
        
        # - check or set taxa_annotation -
        if (is.null(taxa_annotation)){
                taxa_annotation <- taxa_names(pruned_ps) # or you could use 
                # taxa_annotation <- get_taxon_names(as.data.frame(tax_table(pruned_ps)))
        }
        
        if (length(taxa_annotation) != ntaxa(pruned_ps)) {
                stop("taxa_annotation did not fit in length to taxa_to_include")
        }
        
        taxa_annotation <- as.character(taxa_annotation)

        taxa_annotation[is.na(taxa_annotation)] <- "NA"
        
        taxa_annotation <- make.unique(taxa_annotation)
        # --
        
        
        # - generate count data frame in which taxa are rows in the order determined by taxa_to_include -
        if (taxa_are_rows(pruned_ps)) {
                pruned_ps <- t(pruned_ps)
        }
        
        DF_CT <- as(otu_table(pruned_ps), "matrix")
        DF_CT <- as.data.frame(t(DF_CT))
        DF_CT <- DF_CT[taxa_to_include, ]
        # --
        
        # - order the samples in DF_CT based on compare and block -
        sam_info <- as(sample_data(pruned_ps), "data.frame")
        sam_info$IDSaver <- rownames(sam_info)
        sam_info[[group_var]] <- factor(sam_info[[group_var]], levels = compare, ordered = TRUE)
        
        if (!is.null(block)) {
                if (!(block %in% colnames(sam_info))) {
                        stop("You entered a block entry that was not a variable in sample_data(physeq)")
                }
                sam_info[[block]] <- factor(sam_info[[block]])
                sam_info <- dplyr::arrange_(sam_info, group_var, block)
        }
        else {
                sam_info <- dplyr::arrange_(sam_info, group_var)
        }
        
        DF_CT <- DF_CT[, sam_info$IDSaver]
        # --
        
        
        # - set the annotation colors (i.e. block and taxa colors #NB: you already set sample_colors) -
        if (!is.null(block)) {
                if (is.null(block_colors)){
                        block_colors <- brewer.pal(max(3, length(levels(sam_info[[block]]))), "Paired") 
                        block_colors <- block_colors[1:length(levels(sam_info[[block]]))]
                        names(block_colors) <- levels(sam_info[[block]])       
                } else {
                        
                        if (!all(levels(sam_info[[block]]) %in% names(block_colors))) {
                                stop("Not all levels in the block factor are assigned to a color in block_colors")
                        }
                        
                        if (!all(areColors(block_colors))) {
                                stop("Stopped because not all entries in block_colors were valid R colors.")
                        }
                        
                }
        } else {
                block_colors <- NULL
        }

        
        # option to color the taxa for example based on significance
        if (!is.null(taxa_info_data)) {
                
                if (nrow(taxa_info_data) != length(taxa_to_include)) {
                        stop("provided taxa_info_data is in length mismatch to taxa_to_include")
                }
                
                if (is.null(taxa_color_factor) || !(taxa_color_factor %in% colnames(taxa_info_data))) {
                        stop("taxa_color_factor is not a column in the provided taxa_info_data")
                }
                
                if (!is.factor(taxa_info_data[[taxa_color_factor]])) {
                        taxa_info_data[[taxa_color_factor]] <- factor(taxa_info_data[[taxa_color_factor]])
                }
                
                if (is.null(taxa_colors)){
                        taxa_colors <- brewer.pal(max(3, length(levels(taxa_info_data[[taxa_color_factor]]))), name = "Set2") 
                        taxa_colors <- taxa_colors[1:length(levels(taxa_info_data[[taxa_color_factor]]))]
                        names(taxa_colors) <- levels(taxa_info_data[[taxa_color_factor]])       
                } else {
                        if (!all(levels(taxa_info_data[[taxa_color_factor]]) %in% names(taxa_colors))) {
                                stop("Not all levels in the taxa_color_factor are assigned to a color in taxa_colors")
                        }
                        
                        if (!all(areColors(block_colors))) {
                                stop("Stopped because not all entries in taxa_colors were valid R colors.")
                        }
                        
                }
                
                taxa_info_data <- dplyr::select_(taxa_info_data, taxa_color_factor) 
                rownames(taxa_info_data) <- taxa_to_include
                taxa_colors <- list(taxa_colors)
                names(taxa_colors) <- colnames(taxa_info_data)
                
        } else {
                taxa_colors <- list()
                
        }
        
        annotation_colors <- list(sample_colors, block_colors)
        names(annotation_colors) <- c(group_var, block)
        annotation_colors <- c(annotation_colors, taxa_colors)
        # -- 
        
        # - test and adjust max_abundance_for_color -
        if (is.null(max_abundance_for_color)) {
                max_abundance_for_color <- max(DF_CT)
        }
        if (max_abundance_for_color < min(DF_CT) && max_abundance_for_color > max(DF_CT)) {
                max_abundance_for_color <- max(DF_CT)
        }
        # --
        
        
        # - do a possible log transform using min_in_data/5 as pseudocounts -
        if (log_transform) {
                min_in_data <- min(DF_CT[DF_CT > 0])
                min_in_data <- min_in_data/5
                DF_CT[DF_CT == 0] <- min_in_data
                DF_CT <- log10(DF_CT)
                ZeroValues <- log10(c(min_in_data, min_in_data + 1e-14))
                gradient_steps <- log10(gradient_steps)
                max_abundance_for_color <- log10(max_abundance_for_color)
        } else {
                ZeroValues <- c(0, 1e-14)
        }
        # --
        
        
        # - set breaks and colors for pheatmap -
        min_in_data <- min(DF_CT[DF_CT > 0]) # the lowest non-zero value
        # you want that the final color gradient covers the values min to max_abundance_for_color
        # 0 values will be set to a different color (usually red or white), values above max_abundance_for_color should be all max_color
        # normalise gradient steps with max_abundance_for_color
        myBreaks <- gradient_steps * max_abundance_for_color
        # add Zero and minimum breaks
        myBreaks <- c(min_in_data, myBreaks) 
        # now myBreaks has Zero values and goes up to max_abundance_for_color (provided the last gradient_steps was 1)
        myColors = viridis(length(myBreaks)) # see help pheatmap, breaks should be 1 element longer than color, now it is same legnth
        
        # myColors contains now the viridis colors that represent the gradient_steps values (= markers).
        # now we want to introduce breaks between these markers and make linear color gradients between the markers
        color_steps_bw_markers <- 10
        myBreaks1 <- lapply(1:(length(myBreaks)-1), function(i) {
                seq(from = myBreaks[i], to = myBreaks[i + 1], length.out = color_steps_bw_markers + 1)[1:color_steps_bw_markers] # in each step the right side marker is not in
        })
        myBreaks <- c(unlist(myBreaks1), myBreaks[length(myBreaks)]) #length(myBreaks) is now length(myBreaks) * color_steps_bw_markers + 1, the markers are at positions 1, 1+color_steps_bw_markers, 1+2*color_steps_bw_markers 
        
        myColors <- unlist(lapply(1:(length(myColors)-1), function(i) {
                colorRampPalette(colors = c(myColors[i], myColors[i+1]))(color_steps_bw_markers)
        }))
        
        # ZeroValues should be white
        myBreaks <- c(ZeroValue, myBreaks)
        myColors <- c("white", myColors) #
        # --
        
        # make sure sam_info and taxa_data that will be used as annotation_row and annotation_col data_frames have only the relevant columns
        sam_info <- dplyr::select(sam_info, which(colnames(sam_info) %in% names(annotation_colors)))
        rownames(sam_info) <- colnames(DF_CT)
        taxa_data <- dplyr::select(taxa_data, which(colnames(taxa_data) %in% names(annotation_colors)))
        rownames(taxa_data) <- rownames(DF_CT)
        # --

        
        # --
        hm.parameters <- list(DF_CT, color = myColors, breaks = myBreaks, border_color = NA, cluster_cols = FALSE, cluster_rows = FALSE,
                              show_rownames = (show_rownames && nrow(DF_CT) < 70), show_colnames = FALSE, annotation_col = sam_info,
                              annotation_row = taxa_data, annotation_colors = annotation_colors, labels_row = taxa_annotation, annotation_names_row = FALSE,
                              annotation_names_col = FALSE, annotation_legend = TRUE, legend = TRUE, font_size = 10, fontsize_row = 8, fontsize_col = 8, fontsize_num = 6)
        do.call("pheatmap", hm.parameters)
        
        
}

