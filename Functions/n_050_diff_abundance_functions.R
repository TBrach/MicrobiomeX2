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
# --



# --
#######################################
### FUNCTION: format_hit_table
#######################################

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


# --
#########################
## plot_heatmap_physeq ##
#########################
# function that uses pheatmap to draw a heatmap of the physeq object, or a pruned version of it determined by taxa_info_df
# INPUTs
# physeq = physeq object
# sample_colors = named list of named color character vectors. The sample_data(physeq) data.frame is used as sample_info_df to color the samples (columns) in the heatmap.
    # Only those variables/columns in sample_info_df that are names in sample_colors are considered for coloring the samples. I.e. if sample_colors = NULL no coloring occurs.
    # NB: if a variable is a name in sample_colors but the color_character is NA, default colors will be used.
    # NB!!: The named character vectors in sample_colors are also used to prune the samples!! This works as follows:
    # The character vectors in sample_colors are checked in order. If the first character vector (sample_colors[[1]]) is for example "Group" and it contains only colors for "A" and "B", 
    # but sample_info_df[["Group"]] contains also "C", the C samples will be removed! The procedure continues for all entries in sample_colors
    # finally the order of names(sample_colors) defines the order of the color bars, the first entry being drawn closest to the heatmap and so on. 
# taxa_info_df = a data frame similar to sample_info_df (see sample_colors) that together with taxa_colors will be used to label the taxa.
    # NB: MORE importantly, it is also used to restrict the included taxa and to order the taxa!! The rownames(taxa_info_df) must be the taxa_names in physeq 
    # of the taxa you want to include in the heatmap in the order you want the taxa to be shown! IF NULL all taxa will be plotted in the order given by physeq
# taxa_colors: see sample_colors
# taxa_annotation: if not NULL will be used as labels_row in pheatmap, i.e. a different labelling of the taxa instead of taxa_names(physeq) or rownames(taxa_info_df)
# max_abundance_for_color: the maximum abundance value to which the maximum color is assigned, all abundances above get also this color, if NULL max(otu_table(physeq))
    # is used
# gradient_steps:  gradient_steps: vector of numbers between 0 and 1, default c(0.15, 0.3, 0.45, 1). Defines how the main colors in the gradient will be positioned between min_in_data and
    # max_abundance_for_color. min_in_data (= lowest non-zero count/abundance) will be added to first position after gradient_steps have been normalized with max_abundance_for_color
# zero_color: defines the color of the 0 values (minimum Value in otu_table(physeq).
# color_function: a function such as the default viridis function that generates color strings when given a number.
# color_steps_bw_markers: numeric value defining how many equally distributed breaks will be introduced between the markers defined by gradient_steps.
# log_transform: option to log10 transform the counts
# drop_color_levels: if TRUE, color levels in sample_colors and taxa_colors are removed (not shown in legend so) when not present in the data
# rest pheatmap arguments

plot_heatmap_physeq <- function (physeq, sample_colors = NULL, taxa_info_df = NULL, taxa_colors = NULL, taxa_annotation = NULL, 
                                 max_abundance_for_color = NULL, gradient_steps = c(0.15, 0.3, 0.45, 1), zero_color = "white",
                                 color_function = viridis,
                                 color_steps_bw_markers = 10, log_transform = FALSE, drop_color_levels = FALSE,
                                 border_color = NA, cluster_cols = FALSE, cluster_rows = FALSE,
                                 show_rownames = TRUE, show_colnames = FALSE, annotation_names_row = TRUE, annotation_names_col = TRUE, annotation_legend = TRUE,
                                 legend = TRUE, font_size = 10, fontsize_row = 8, fontsize_col = 8, fontsize_number = 6, filename = NA, ...) {
        
        # - keep only taxa in taxa_info_df in physeq -
        if (!is.null(taxa_info_df)) {
                
                if (!all(rownames(taxa_info_df) %in% taxa_names(physeq))) {
                        stop("not all taxa in taxa_info_df are taxa in physeq!")
                }
                pruned_ps <- prune_taxa(rownames(taxa_info_df), physeq)
        } else {
                pruned_ps <- physeq
                taxa_colors <- NULL
        }
        # --
        
        # - test sample_colors input and prepare sample_info_df -
        if (!is.null(sample_colors)) {
                
                if (!is.list(sample_colors)) {
                        stop("sample_colors must be a named list.")
                }
                
                
                if (!all(names(sample_colors) %in% colnames(sample_data(pruned_ps)))) {
                        stop("all names of the sample_colors list must be variables/columns in sample_data(physeq).")
                }
                
                sample_info_df <- as(sample_data(pruned_ps), "data.frame")
                sample_info_df <- sample_info_df[, match(names(sample_colors), colnames(sample_info_df)), drop = FALSE] # the match construct makes sure that sample_info_df columns are ordered after names(sample_colors)
                # which makes sure that the column coloring occurs in the whished for order
                
                # -- test all entries in sample_colors, prune pruned_ps and order sample_info_df, and set default color where necessary --
                for (i in 1:length(sample_colors)){
                        variable_name <- names(sample_colors)[i]
                        color_character <- sample_colors[[i]]
                        
                        if (!all(is.na(color_character))){
                                
                                # to not remove all samples from pruned_ps:
                                if (!(sum(names(color_character) %in% unique(sample_info_df[[variable_name]])) > 0)) {
                                        stop(paste("Not a single name of a color for ", variable_name, " matched to an entry in the corresponding column in sample_data(physeq)", sep = ""))
                                }
                                
                                # remove all samples that no color was assigned to in sample_colors, i.e.option to restrict the comparison to defined samples
                                if (!all(unique(sample_info_df[[variable_name]]) %in% names(color_character))) {
                                        
                                        keepSamples <- sample_names(pruned_ps)[sample_info_df[[variable_name]] %in% names(color_character)]
                                        pruned_ps <- prune_samples(keepSamples, pruned_ps)
                                        sample_info_df <- as(sample_data(pruned_ps), "data.frame")
                                        sample_info_df <- sample_info_df[, colnames(sample_info_df) %in% names(sample_colors), drop = FALSE]
                                }
                                
                                # remove unused color levels if asked for it
                                if (drop_color_levels && !all(names(color_character) %in% unique(sample_info_df[[variable_name]]))) {
                                        color_character <- color_character[names(color_character) %in% unique(sample_info_df[[variable_name]])]
                                        sample_colors[[i]] <- color_character
                                }
                                
                                if (!all(areColors(color_character))) {
                                        warning(paste("Not all entries in sample_colors entry ", variable_name, " were colors. Assigned default colors.", sep = ""))
                                        color_character <- assign_default_colors(sample_info_df, variable_name)
                                        sample_colors[[i]] <- color_character
                                }
                                
                                # define the order of the samples based on the sample_color vectors (actual ordering below with dplyr::arrange)
                                sample_info_df[[variable_name]] <- factor(sample_info_df[[variable_name]], levels = names(color_character), ordered = TRUE)
                                
                        } else {
                                color_character <- assign_default_colors(sample_info_df, variable_name)
                                sample_colors[[i]] <- color_character
                        }
                }
                sample_info_df$IDSaver <- rownames(sample_info_df)
                sample_info_df <- dplyr::arrange_(sample_info_df, names(sample_colors))
                rownames(sample_info_df) <- sample_info_df$IDSaver
                sample_info_df <- sample_info_df[, -ncol(sample_info_df), drop = FALSE] # to remove IDSaver
                # ----
        } else {
                sample_info_df <- NULL
        }
        # --
        
        # - set the taxa colors (and pool sample_colors and taxa_colors into annotation_colors) -
        if (!is.null(taxa_colors)) { # see above, guarantees also that taxa_info_id is not NULL
                
                if (!is.list(taxa_colors)) {
                        stop("taxa_colors must be a named list.")
                }
                
                if (!all(names(taxa_colors) %in% colnames(taxa_info_df))) {
                        stop("all names of the taxa_colors list must be variables/columns in taxa_info_df.")
                }
                
                taxa_info_df <- taxa_info_df[, match(names(taxa_colors), colnames(taxa_info_df)), drop = FALSE] # again match to get coloring ordered after taxa_colors
                
                # -- test all entries in taxa_colors --
                for (i in 1:length(taxa_colors)){
                        variable_name <- names(taxa_colors)[i]
                        color_character <- taxa_colors[[i]]
                        
                        if (!all(is.na(color_character))){
                                
                                
                                if (!all(unique(taxa_info_df[[variable_name]]) %in% names(color_character))) {
                                        warning(paste("There were levels in taxa_info_df at ", variable_name, " for which no color was assigned in taxa_colors. Assigned default colors.", sep = ""))
                                        color_character <- assign_default_colors(taxa_info_df, variable_name)
                                        taxa_colors[[i]] <- color_character
                                }
                                
                                # remove unused color levels if asked for it
                                if (drop_color_levels && !all(names(color_character) %in% unique(taxa_info_df[[variable_name]]))) {
                                        color_character <- color_character[names(color_character) %in% unique(taxa_info_df[[variable_name]])]
                                        taxa_colors[[i]] <- color_character
                                }
                                
                                
                                if (!all(areColors(color_character))) {
                                        warning(paste("Not all entries in taxa_colors entry ", variable_name, " were colors. Assigned default colors.", sep = ""))
                                        color_character <- assign_default_colors(taxa_info_df, variable_name)
                                        taxa_colors[[i]] <- color_character
                                }
                                
                                # define the order of the taxa coloring based on the taxa_color vectors
                                # taxa_info_df[[variable_name]] <- factor(taxa_info_df[[variable_name]], levels = names(color_character), ordered = TRUE) # was unnecessary, pheatmap orders based on annotation_colors
                                
                        } else {
                                color_character <- assign_default_colors(taxa_info_df, variable_name)
                                taxa_colors[[i]] <- color_character
                        }
                }
                # ----
        } 
        
        annotation_colors <- c(sample_colors, taxa_colors)
        # -- 
        
        # - check or set taxa_annotation -
        if (is.null(taxa_annotation)){
                taxa_annotation <- taxa_names(pruned_ps) # or you could use 
                # taxa_annotation <- get_taxon_names(as.data.frame(tax_table(pruned_ps)))
        }
        
        if (length(taxa_annotation) != ntaxa(pruned_ps)) {
                warning("taxa_annotation did not fit in length to nrow(taxa_info_df) or ntaxa(physeq), used taxa_names")
                taxa_annotation <- taxa_names(pruned_ps)
        }
        
        taxa_annotation <- as.character(taxa_annotation)
        
        taxa_annotation[is.na(taxa_annotation)] <- "NA"
        
        taxa_annotation <- make.unique(taxa_annotation)
        # --
        
        # - generate count data frame in which taxa are rows in the order determined by taxa_info_df -
        if (taxa_are_rows(pruned_ps)) {
                pruned_ps <- t(pruned_ps)
        }
        
        DF_CT <- as(otu_table(pruned_ps), "matrix")
        DF_CT <- as.data.frame(t(DF_CT))
        
        if (!is.null(taxa_info_df)){
                DF_CT <- DF_CT[rownames(taxa_info_df), ]
        }
        # --
        
        # - order the samples in DF_CT based on sample_colors list (see above for sample_info_df) -
        if (!is.null(sample_info_df)){
                DF_CT <- DF_CT[, rownames(sample_info_df)]
        }
        # --
        
        # - test and adjust max_abundance_for_color -
        if (is.null(max_abundance_for_color)) {
                max_abundance_for_color <- max(DF_CT)
        }
        if (max_abundance_for_color < min(DF_CT) && max_abundance_for_color > max(DF_CT)) {
                max_abundance_for_color <- max(DF_CT)
        }
        # --
        
        # - set breaks and colors for pheatmap and do a possible log transform -
        ZeroValue <- min(DF_CT) # should be 0 in most cases of microbiome data, unless pseudocount had been added. In rare cases of high taxonomy it might be a non-zero count so actually min_in_data 
        min_in_data <- min(DF_CT[DF_CT > ZeroValue]) # the lowest non-zero value
        max_in_data <- max(DF_CT)
        
        # -- make sure the last entry in gradient steps = 1 --
        if (!all(gradient_steps >= 0 & gradient_steps <= 1)) {
                gradient_steps <- c(0.15, 0.3, 0.45, 1)
        }
        
        if (gradient_steps[length(gradient_steps)] != 1) {
                gradient_steps <- c(gradient_steps, 1)
        }
        # ----
        
        # you want that the final color gradient covers the values min to max_abundance_for_color
        # 0 values will be set to a different color (usually red or white), values above max_abundance_for_color should be all max_color
        # normalise gradient steps with max_abundance_for_color
        myBreaks <- gradient_steps * max_abundance_for_color
        myBreaks <- myBreaks[myBreaks > min_in_data]
        # add break at min_in_data
        myBreaks <- c(min_in_data, myBreaks)
        
        # now myBreaks goes from min_in_data up to max_abundance_for_color (provided the last gradient_steps was 1)
        myColors = color_function(length(myBreaks)) # see help pheatmap, breaks should be 1 element longer than color, now it is same legnth
        
        # myColors contains now the colors that represent the gradient_steps values (= markers).
        # now we want to introduce breaks between these markers and make linear color gradients between the markers. The number of breaks between the markers is defined
        # by color_steps_bw_markers
        
        myColors <- unlist(lapply(1:(length(myColors)-1), function(i) {
                colorRampPalette(colors = c(myColors[i], myColors[i+1]))(color_steps_bw_markers)
        }))
        
        
        # -- do a possible log transform using min_in_data/5 as pseudocounts --
        if (log_transform) {
                if (ZeroValue > 0){ # > 0 if pseudocount had been added or simply no zero counts in data, < 0 if data h
                        pseudocount <- ZeroValue
                } else if (ZeroValue < 0){
                        stop("you asked for log_transform but the lowest count in your data is already below 0!")
                } else {
                        pseudocount <- min_in_data/2
                }
                DF_CT[DF_CT == ZeroValue] <- pseudocount
                DF_CT <- log10(DF_CT)
                
                myBreaks <- log10(myBreaks)
                myBreaks1 <- lapply(1:(length(myBreaks)-1), function(i) {
                        seq(from = myBreaks[i], to = myBreaks[i + 1], length.out = color_steps_bw_markers + 1)[1:color_steps_bw_markers] # in each step the right side marker is not in
                })
                myBreaks <- c(unlist(myBreaks1), myBreaks[length(myBreaks)]) #length(myBreaks) is now length(myBreaks) * color_steps_bw_markers + 1, the markers are at positions 1, 1+color_steps_bw_markers, 1+2*color_steps_bw_markers 
                
                if (log10(max_in_data) > myBreaks[length(myBreaks)]) {
                        myBreaks <- c(log10(pseudocount), myBreaks, log10(max_in_data))
                        myColors <- c(zero_color, myColors, myColors[length(myColors)])
                } else {
                        myBreaks <- c(log10(pseudocount), myBreaks)
                        myColors <- c(zero_color, myColors)
                }
                
        } else {
                myBreaks1 <- lapply(1:(length(myBreaks)-1), function(i) {
                        seq(from = myBreaks[i], to = myBreaks[i + 1], length.out = color_steps_bw_markers + 1)[1:color_steps_bw_markers] # in each step the right side marker is not in
                })
                myBreaks <- c(unlist(myBreaks1), myBreaks[length(myBreaks)])
                
                # add max_in_data values to myBreaks, otherwise all values above max_abundance_for_color will be white
                if (max_in_data > myBreaks[length(myBreaks)]) {
                        myBreaks <- c(ZeroValue, myBreaks, max_in_data)
                        myColors <- c(zero_color, myColors, myColors[length(myColors)])
                } else {
                        myBreaks <- c(ZeroValue, myBreaks)
                        myColors <- c(zero_color, myColors)
                }
        }
        # ----
        # --
        
        # - not necessary, just for clarity -
        if (is.null(annotation_colors)) {annotation_colors <- NA}
        if (is.null(sample_info_df)) {sample_info_df <- NA}
        if (is.null(taxa_info_df)) {taxa_info_df <- NA}
        # --
        
        # --
        hm.parameters <- list(DF_CT, color = myColors, breaks = myBreaks, border_color = border_color, cluster_cols = cluster_cols, cluster_rows = cluster_rows,
                              show_rownames = show_rownames, show_colnames = show_colnames, annotation_col = sample_info_df,
                              annotation_row = taxa_info_df, annotation_colors = annotation_colors, labels_row = taxa_annotation, annotation_names_row = annotation_names_row,
                              annotation_names_col = annotation_names_col, annotation_legend = annotation_legend, legend = legend, font_size = font_size, 
                              fontsize_row = fontsize_row, fontsize_col = fontsize_col, fontsize_number = fontsize_number)
        do.call("pheatmap", hm.parameters)
        
        
}
# --



# --
#######################################
### FUNCTION: assign_default_colors
#######################################

assign_default_colors <- function(info_df, variable_name){
        
        if (!(variable_name %in% colnames(info_df))) {
                stop("variable_name must be a variable/column in info_df")   
        }
        
        if (is.factor(info_df[[variable_name]])) {
                uniqueEntries <- levels(info_df[[variable_name]])
        } else {
                uniqueEntries <- unique(info_df[[variable_name]])
        }
        
        
        if (length(uniqueEntries) < 9) {
                color_char <- cbPalette[1:length(uniqueEntries)]
        } else if (length(uniqueEntries) < 16) {
                color_char <- QuantColors15[1:length(uniqueEntries)]
        } else {
                color_char <- viridis(length(uniqueEntries))
        }
        
        names(color_char) <- uniqueEntries
        
        color_char
}
# --




















# OLDER Versions #################


# plot_heatmap_physeqOld <- function (physeq, sample_colors = NULL, taxa_info_df = NULL, taxa_colors = NULL, taxa_annotation = NULL, 
#                                     max_abundance_for_color = NULL, gradient_steps = c(0.15, 0.3, 0.45, 1), color_function = viridis, zero_color = "white",
#                                     color_steps_bw_markers = 10, log_transform = FALSE, border_color = NA, cluster_cols = FALSE, cluster_rows = FALSE,
#                                     show_rownames = TRUE, show_colnames = FALSE, annotation_names_row = TRUE, annotation_names_col = TRUE, annotation_legend = TRUE,
#                                     legend = TRUE, font_size = 10, fontsize_row = 8, fontsize_col = 8, fontsize_number = 6, filename = NA, ...) {
#         
#         # - keep only taxa in taxa_info_df in physeq -
#         if (!is.null(taxa_info_df)) {
#                 
#                 if (!all(rownames(taxa_info_df) %in% taxa_names(physeq))) {
#                         stop("not all taxa in taxa_info_df are taxa in physeq!")
#                 }
#                 pruned_ps <- prune_taxa(rownames(taxa_info_df), physeq)
#         } else {
#                 pruned_ps <- physeq
#                 taxa_colors <- NULL
#         }
#         # --
#         
#         # - test sample_colors input and prepare sample_info_df -
#         if (!is.null(sample_colors)) {
#                 
#                 if (!is.list(sample_colors)) {
#                         stop("sample_colors must be a named list.")
#                 }
#                 
#                 
#                 if (!all(names(sample_colors) %in% colnames(sample_data(pruned_ps)))) {
#                         stop("all names of the sample_colors list must be variables/columns in sample_data(physeq).")
#                 }
#                 
#                 sample_info_df <- as(sample_data(pruned_ps), "data.frame")
#                 sample_info_df <- sample_info_df[, colnames(sample_info_df) %in% names(sample_colors), drop = FALSE]
#                 
#                 # -- test all entries in sample_colors, prune pruned_ps and order sample_info_df, and set default color where necessary --
#                 for (i in 1:length(sample_colors)){
#                         variable_name <- names(sample_colors)[i]
#                         color_character <- sample_colors[[i]]
#                         
#                         if (!all(is.na(color_character))){
#                                 
#                                 # to not remove all samples from pruned_ps:
#                                 if (!(sum(names(color_character) %in% unique(sample_info_df[[variable_name]])) > 0)) {
#                                         stop(paste("Not a single name of a color for ", variable_name, " matched to an entry in the corresponding column in sample_data(physeq)", sep = ""))
#                                 }
#                                 
#                                 # remove all samples that no color was assigned to in sample_colors, i.e.option to restrict the comparison to defined samples
#                                 if (!all(unique(sample_info_df[[variable_name]]) %in% names(color_character))) {
#                                         
#                                         keepSamples <- sample_names(pruned_ps)[sample_info_df[[variable_name]] %in% names(color_character)]
#                                         pruned_ps <- prune_samples(keepSamples, pruned_ps)
#                                         sample_info_df <- as(sample_data(pruned_ps), "data.frame")
#                                         sample_info_df <- sample_info_df[, colnames(sample_info_df) %in% names(sample_colors), drop = FALSE]
#                                 }
#                                 
#                                 if (!all(areColors(color_character))) {
#                                         warning(paste("Not all entries in sample_colors entry ", variable_name, " were colors. Assigned default colors.", sep = ""))
#                                         color_character <- assign_default_colors(sample_info_df, variable_name)
#                                         sample_colors[[i]] <- color_character
#                                 }
#                                 
#                                 # define the order of the samples based on the sample_color vectors (actual ordering below with dplyr::arrange)
#                                 sample_info_df[[variable_name]] <- factor(sample_info_df[[variable_name]], levels = names(color_character), ordered = TRUE)
#                                 
#                         } else {
#                                 color_character <- assign_default_colors(sample_info_df, variable_name)
#                                 sample_colors[[i]] <- color_character
#                         }
#                 }
#                 sample_info_df$IDSaver <- rownames(sample_info_df)
#                 sample_info_df <- dplyr::arrange_(sample_info_df, names(sample_colors))
#                 rownames(sample_info_df) <- sample_info_df$IDSaver
#                 sample_info_df <- sample_info_df[, -ncol(sample_info_df), drop = FALSE] # to remove IDSaver
#                 # ----
#         } else {
#                 sample_info_df <- NULL
#         }
#         # --
#         
#         # - set the taxa colors (and pool sample_colors and taxa_colors into annotation_colors) -
#         if (!is.null(taxa_colors)) { # see above, guarantees also that taxa_info_id is not NULL
#                 
#                 if (!is.list(taxa_colors)) {
#                         stop("taxa_colors must be a named list.")
#                 }
#                 
#                 if (!all(names(taxa_colors) %in% colnames(taxa_info_df))) {
#                         stop("all names of the taxa_colors list must be variables/columns in taxa_info_df.")
#                 }
#                 
#                 taxa_info_df <- taxa_info_df[, colnames(taxa_info_df) %in% names(taxa_colors), drop = FALSE]
#                 
#                 # -- test all entries in taxa_colors --
#                 for (i in 1:length(taxa_colors)){
#                         variable_name <- names(taxa_colors)[i]
#                         color_character <- taxa_colors[[i]]
#                         
#                         if (!all(is.na(color_character))){
#                                 
#                                 
#                                 if (!all(unique(taxa_info_df[[variable_name]]) %in% names(color_character))) {
#                                         warning(paste("There were levels in taxa_info_df at ", variable_name, " for which no color was assigned in taxa_colors. Assigned default colors.", sep = ""))
#                                         color_character <- assign_default_colors(taxa_info_df, variable_name)
#                                         taxa_colors[[i]] <- color_character
#                                 }
#                                 
#                                 
#                                 if (!all(areColors(color_character))) {
#                                         warning(paste("Not all entries in taxa_colors entry ", variable_name, " were colors. Assigned default colors.", sep = ""))
#                                         color_character <- assign_default_colors(taxa_info_df, variable_name)
#                                         taxa_colors[[i]] <- color_character
#                                 }
#                                 
#                                 # define the order of the taxa coloring based on the taxa_color vectors
#                                 # taxa_info_df[[variable_name]] <- factor(taxa_info_df[[variable_name]], levels = names(color_character), ordered = TRUE) # was unnecessary, pheatmap orders based on annotation_colors
#                                 
#                         } else {
#                                 color_character <- assign_default_colors(taxa_info_df, variable_name)
#                                 taxa_colors[[i]] <- color_character
#                         }
#                 }
#                 # ----
#         } 
#         
#         annotation_colors <- c(sample_colors, taxa_colors)
#         # -- 
#         
#         # - check or set taxa_annotation -
#         if (is.null(taxa_annotation)){
#                 taxa_annotation <- taxa_names(pruned_ps) # or you could use 
#                 # taxa_annotation <- get_taxon_names(as.data.frame(tax_table(pruned_ps)))
#         }
#         
#         if (length(taxa_annotation) != ntaxa(pruned_ps)) {
#                 warning("taxa_annotation did not fit in length to nrow(taxa_info_df) or ntaxa(physeq), used taxa_names")
#                 taxa_annotation <- taxa_names(pruned_ps)
#         }
#         
#         taxa_annotation <- as.character(taxa_annotation)
#         
#         taxa_annotation[is.na(taxa_annotation)] <- "NA"
#         
#         taxa_annotation <- make.unique(taxa_annotation)
#         # --
#         
#         # - generate count data frame in which taxa are rows in the order determined by taxa_info_df -
#         if (taxa_are_rows(pruned_ps)) {
#                 pruned_ps <- t(pruned_ps)
#         }
#         
#         DF_CT <- as(otu_table(pruned_ps), "matrix")
#         DF_CT <- as.data.frame(t(DF_CT))
#         
#         if (!is.null(taxa_info_df)){
#                 DF_CT <- DF_CT[rownames(taxa_info_df), ]
#         }
#         # --
#         
#         # - order the samples in DF_CT based on sample_colors list (see above for sample_info_df) -
#         if (!is.null(sample_info_df)){
#                 DF_CT <- DF_CT[, rownames(sample_info_df)]
#         }
#         # --
#         
#         # - test and adjust max_abundance_for_color -
#         if (is.null(max_abundance_for_color)) {
#                 max_abundance_for_color <- max(DF_CT)
#         }
#         if (max_abundance_for_color < min(DF_CT) && max_abundance_for_color > max(DF_CT)) {
#                 max_abundance_for_color <- max(DF_CT)
#         }
#         # --
#         
#         
#         # - set breaks and colors for pheatmap and do a possible log transform -
#         ZeroValue <- 0
#         min_in_data <- min(DF_CT[DF_CT > 0]) # the lowest non-zero value
#         max_in_data <- max(DF_CT)
#         
#         # -- make sure the last entry in gradient steps = 1 --
#         if (!all(gradient_steps >= 0 & gradient_steps <= 1)) {
#                 gradient_steps <- c(0.15, 0.3, 0.45, 1)
#         }
#         
#         if (gradient_steps[length(gradient_steps)] != 1) {
#                 gradient_steps <- c(gradient_steps, 1)
#         }
#         # ----
#         
#         # you want that the final color gradient covers the values min to max_abundance_for_color
#         # 0 values will be set to a different color (usually red or white), values above max_abundance_for_color should be all max_color
#         # normalise gradient steps with max_abundance_for_color
#         myBreaks <- gradient_steps * max_abundance_for_color
#         # add break at min_in_data
#         myBreaks <- c(min_in_data, myBreaks)
#         
#         # now myBreaks goes from min_in_data up to max_abundance_for_color (provided the last gradient_steps was 1)
#         myColors = color_function(length(myBreaks)) # see help pheatmap, breaks should be 1 element longer than color, now it is same legnth
#         
#         # myColors contains now the viridis colors that represent the gradient_steps values (= markers).
#         # now we want to introduce breaks between these markers and make linear color gradients between the markers
#         myBreaks1 <- lapply(1:(length(myBreaks)-1), function(i) {
#                 seq(from = myBreaks[i], to = myBreaks[i + 1], length.out = color_steps_bw_markers + 1)[1:color_steps_bw_markers] # in each step the right side marker is not in
#         })
#         myBreaks <- c(unlist(myBreaks1), myBreaks[length(myBreaks)]) #length(myBreaks) is now length(myBreaks) * color_steps_bw_markers + 1, the markers are at positions 1, 1+color_steps_bw_markers, 1+2*color_steps_bw_markers 
#         
#         myColors <- unlist(lapply(1:(length(myColors)-1), function(i) {
#                 colorRampPalette(colors = c(myColors[i], myColors[i+1]))(color_steps_bw_markers)
#         }))
#         
#         # add max_in_data values to myBreaks, otherwise all values above max_abundance_for_color will be white
#         if (max_in_data > myBreaks[length(myBreaks)]) {
#                 myBreaks <- c(ZeroValue, myBreaks, max_in_data)
#                 myColors <- c(zero_color, myColors, myColors[length(myColors)])
#         } else {
#                 myBreaks <- c(ZeroValue, myBreaks)
#                 myColors <- c(zero_color, myColors)
#         } 
#         
#         # -- do a possible log transform using min_in_data/5 as pseudocounts --
#         if (log_transform) {
#                 pseudocount <- min_in_data/5
#                 DF_CT[DF_CT == ZeroValue] <- pseudocount
#                 DF_CT <- log10(DF_CT)
#                 myBreaks[1] <- pseudocount
#                 myBreaks <- log10(myBreaks)
#         } 
#         # ----
#         
#         
#         
#         # --
#         
#         # - not necessary, just for clarity -
#         if (is.null(annotation_colors)) {annotation_colors <- NA}
#         if (is.null(sample_info_df)) {sample_info_df <- NA}
#         if (is.null(taxa_info_df)) {taxa_info_df <- NA}
#         # --
#         
#         # --
#         hm.parameters <- list(DF_CT, color = myColors, breaks = myBreaks, border_color = border_color, cluster_cols = cluster_cols, cluster_rows = cluster_rows,
#                               show_rownames = show_rownames, show_colnames = show_colnames, annotation_col = sample_info_df,
#                               annotation_row = taxa_info_df, annotation_colors = annotation_colors, labels_row = taxa_annotation, annotation_names_row = annotation_names_row,
#                               annotation_names_col = annotation_names_col, annotation_legend = annotation_legend, legend = legend, font_size = font_size, 
#                               fontsize_row = fontsize_row, fontsize_col = fontsize_col, fontsize_number = fontsize_number)
#         do.call("pheatmap", hm.parameters)
#         
#         
# }




