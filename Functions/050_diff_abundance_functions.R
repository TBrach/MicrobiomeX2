# --
#######################################
### FUNCTION: test_diffs_in_prevalence
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

test_diffs_in_prevalence <- function(physeq, group_var, p.adj.method = "fdr", minCount = 0L, symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) {
        
        
        ### The Fisher test of pooled successes
        # MatrixFisher <- matrix(c(200, 80, 200, 320), nrow = 2)
        # rownames(MatrixFisher) <- c("Strain A", "Strain B")
        # colnames(MatrixFisher) <- c("Success", "Fail")
        # fisher.test(MatrixFisher, alternative = "two")
        # PValueFisher <- fisher.test(MatrixFisher, alternative = "two")$p.value
        
        if (taxa_are_rows(physeq)) { physeq <- t(physeq)}
        
        CT <- as(otu_table(physeq), "matrix")
        CT <- CT > minCount
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        fac_levels <- levels(group_fac)
        
        prev_list <- lapply(fac_levels, function(level){
                data.frame(Present = colSums(CT[group_fac == level, ]),
                           Absent = colSums(!CT[group_fac == level, ]))
        })
        
        names(prev_list) <- fac_levels
        
        # - get the level combis -
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
        i_sAndj_s <- get_unique_facLevel_combis(fac_levels)
        i_s <- i_sAndj_s[["i_s"]]
        j_s <- i_sAndj_s[["j_s"]]
        # --
        
        p_val_list <- vector("list", length = length(i_s))
        for (k in seq_along(i_s)){
                i <- i_s[k]
                j <- j_s[k]
                prev_PC_1 <- round(100*(prev_list[[i]][, 1]/(prev_list[[i]][1,1] + prev_list[[i]][1,2])), 1)
                prev_PC_2 <- round(100*(prev_list[[j]][, 1]/(prev_list[[j]][1,1] + prev_list[[j]][1,2])), 1)
                direction <- rep("down", length(prev_PC_1))
                direction[prev_PC_1 > prev_PC_2] <- "up"
                rowwise_compare_matrix <- cbind(prev_list[[i]], prev_list[[j]])
                p_vals <- sapply(1:nrow(rowwise_compare_matrix), function(e){
                        mat_fisher <- matrix(c(rowwise_compare_matrix[e, 1],
                                               rowwise_compare_matrix[e, 3],
                                               rowwise_compare_matrix[e, 2],
                                               rowwise_compare_matrix[e, 4]), ncol = 2)
                        fisher.test(mat_fisher)$p.value
                })
                p_vals_adj <- p.adjust(p_vals, p.adj.method)
                symnum.args$x <- p_vals
                significance <- do.call(stats::symnum, symnum.args) %>% as.character()
                symnum.args$x <- p_vals_adj
                significance_adj <- do.call(stats::symnum, symnum.args) %>% as.character()
                p_val_df <- data.frame(p_vals, p_vals_adj, significance, significance_adj, prev_PC_1, prev_PC_2, direction)
                # colnames(p_val_df)[1:2] <- paste(c("p_val_", "p_val_adj_"), fac_levels[i], "_vs_", fac_levels[j], sep = "")
                # colnames(p_val_df)[5:6] <- paste(c("prev_PC_"), c(fac_levels[i], fac_levels[j]), sep = "")
                colnames(p_val_df)[1:6] <- c("p_val", "p_val_adj", "signi.", "signi_adj", "prev_PC_grp1", "prev_PC_grp2")
                p_val_df <- cbind(as.data.frame(p_val_df), tax_table(physeq))
                p_val_df$Taxon <- colnames(CT)
                p_val_df <- arrange(p_val_df, p_val)
                p_val_df <- select(p_val_df, Taxon, 1:(ncol(p_val_df)-1))
                p_val_list[[k]] <- p_val_df
                names(p_val_list)[[k]] <- paste(fac_levels[i], "_vs_", fac_levels[j], sep = "")
        }
        
        p_val_list
        
}
# --



# --
#######################################
### prepare_diff_abundance_results_for_plotting
#######################################
# function exclusively to save some typing, and for forcing me to make the outcomes of the diff_abundance test functions more uniform
# TbT is just because the TbT method has currently no p_val_adj yet! Has to be changed

prepare_diff_abundance_results_for_plotting <- function(res_list, physeq, taxonomic_level, TbT = "no") {
        
        if(TbT == "yes"){
                head_values <- rep(10, times = length(res_list))
                names(head_values) <- names(res_list)
        } else {
                suppressWarnings(head_values <- sapply(res_list, function(df){
                        max(which(df[, "p_val_adj"] < 0.05))
                }))
        }
        
        original_head_values <- data.frame(Comparison = c(names(head_values), "Total"), NoSignificant = c(head_values, ntaxa(physeq)),
                                           PC_Significant = 100*c(head_values, ntaxa(physeq))/ntaxa(physeq))
        
        # Show at least 10 taxa even if less are significant and show max 25 even if far more are significant
        head_values[head_values < 10] <- 10
        head_values[head_values > 25] <- 25
        
        head_values[head_values > ntaxa(physeq)] <- ntaxa(physeq) # for "Phylum" analysis it is possible that less than 10 taxa are in physeq
        
        res_table_list <- lapply(1:length(res_list), function(i){
                df <- head(res_list[[i]], head_values[i])
                df <- dplyr::select(df, which(colnames(df) == taxonomic_level), 2:ncol(df)) #simply put taxonomic level in front and remove taxon column (= column 1) 
                df
        })
        names(res_table_list) <- names(res_list)
        
        
        row_names_for_heat_maps <- lapply(res_table_list, function(df){
                the_names <- as.character(df[, taxonomic_level])
                the_names <- sapply(strsplit(the_names, split = "/"), `[`, 1)
                the_names <- paste(the_names, df$signi., df$signi_adj, sep = "_")
        }) # in case of ambiguous species assignment keep only first one
        
        tax_orders <- lapply(1:length(res_list), function(i){
                df <- res_list[[i]]
                as.character(df$Taxon[1:head_values[i]])
        })
        
        pruned_physeqs_to_test <- lapply(1:length(tax_orders), function(i){
                prune_taxa(tax_orders[[i]], physeq)
        })
        
        list(original_head_values = original_head_values, head_values = head_values, res_table_list = res_table_list, row_names_for_heat_maps = row_names_for_heat_maps,
             tax_orders = tax_orders, pruned_physeqs_to_test = pruned_physeqs_to_test)
        
}
# --



# --
#######################################
### make_heat_map_physeq_levels##
#################
# same heat map plots as from make_heat_map_physeq but for each combination of levels in your physeq
## Input:
# physeq object
# group_var: the name of the group_fac column in sample_data used to order the samples in the heat map
# max_abundance_for_color: if null the 90% percentile count/relative abundance in the data is used. all counts above this value will be
# shown yellow in the heat map
# tax_order: character vector of the original taxon names in the order you want them shown. if NULL >> tax_order = taxa_names(physeq)
# tax_names: the names that will be used for the taxons, if Null Taxon_1, Taxon_2 and so on will be used. NB: renaming of course after
# ordering. 
# color_sample_names: if TRUE and if you have less than 7 levels, the sample names will be colored using cbPalette
# gradient_steps: the steps the blue to green to yellow will be distributed in the viridis gradient: 4 numbers ranging from 1e-14 to 1,
# see default, you might wanna try c(0.25, 0.5, 0.75, 1) as well
# Output: list of heat maps for each level combination

make_heat_map_physeq_levels <- function(physeq, group_var, color_levels, max_abundance_for_color = NULL, tax_order = NULL,
                                        tax_names = NULL, color_sample_names = TRUE, gradient_steps = c(0.15, 0.3, 0.45, 1)){
        
        gradient_steps <- c(0, 1e-14, gradient_steps)
        
        if(taxa_are_rows(physeq)) { physeq <- t(physeq) }
        
        if (!is.factor(sample_data(physeq)[[group_var]])) {sample_data(physeq)[[group_var]] <- as.factor(sample_data(physeq)[[group_var]])}
        
        CT <- as(otu_table(physeq), "matrix") # taxa are columns
        DF_CT <- as.data.frame(t(CT))
        
        if (is.null(tax_order)){tax_order <- taxa_names(physeq)}
        
        if (length(tax_order) != nrow(DF_CT) || !all(tax_order %in% rownames(DF_CT))) {stop("the given tax_order must contain all taxa_names(physeq). If you want to subset,
                                                                                            do it before with prune_taxa")}
        DF_CT <- DF_CT[tax_order, ]
        
        if (is.null(tax_names)) {tax_names <- paste("Taxon_", 1:nrow(DF_CT), sep = "")}
        
        if (length(tax_names) != nrow(DF_CT) ) {stop("the given tax_names must be a character vector of length = ntaxa(physeq)")}
        
        tax_names[is.na(tax_names)] <- "NA"
        rownames(DF_CT) <- make.unique(tax_names)
        DF_CT$Taxa <- rownames(DF_CT)
        DF_CT <- tidyr::gather(DF_CT, key = Sample , value = Count, -Taxa)
        DF_CT$Taxa <- factor(DF_CT$Taxa, levels = rev(unique(DF_CT$Taxa)), ordered = TRUE)
        LookUpDF <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
        LookUpDF <- LookUpDF[order(match(LookUpDF$Group, levels(LookUpDF$Group))), ]
        DF_CT$Sample <- factor(DF_CT$Sample, levels = LookUpDF$Sample, ordered = TRUE)
        DF_CT$Group <- LookUpDF$Group[match(DF_CT$Sample, LookUpDF$Sample)]
        
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        fac_levels <- levels(group_fac)
        
        # - get the level combis -
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
        i_sAndj_s <- get_unique_facLevel_combis(fac_levels)
        i_s <- i_sAndj_s[["i_s"]]
        j_s <- i_sAndj_s[["j_s"]]
        # --
        
        plot_list <- vector("list", length = length(i_s))
        
        
        
        for (k in seq_along(i_s)){
                i <- i_s[k]
                j <- j_s[k]
                
                group_fac_current <- droplevels(group_fac[as.numeric(group_fac) %in% c(i, j)])
                DF_CT_current <- filter(DF_CT, Group %in% group_fac_current)
                DF_CT_current$Group <- factor(DF_CT_current$Group, levels = c(fac_levels[i], fac_levels[j]), ordered = TRUE)
                LookUpDF_current <- LookUpDF[LookUpDF$Sample %in% DF_CT_current$Sample, ]
                DF_CT_current$Sample <- factor(DF_CT_current$Sample, levels = LookUpDF_current$Sample, ordered = TRUE)
                
                if (is.null(max_abundance_for_color)) {
                        max_abundance_for_color_current <- quantile(DF_CT_current$Count, .9)
                } else {
                        max_abundance_for_color_current <- max_abundance_for_color
                }
                
                if (max_abundance_for_color_current == 0) {max_abundance_for_color_current = max(DF_CT_current$Count)}
                
                # Color the sample names based on color_levels
                colxaxis <- color_levels[LookUpDF_current$Group]
                
                
                hmTr <- ggplot(DF_CT_current, aes(x = Sample, y = Taxa, fill = Count))
                hmTr <- hmTr + 
                        geom_raster() + 
                        scale_fill_gradientn("", limits = c(0, max_abundance_for_color_current), colors = c("red", viridis(5)), values = gradient_steps, oob = squish) +
                        scale_x_discrete(position = "top") +
                        #coord_equal() +
                        labs(x=NULL, y=NULL) +
                        theme_tufte(base_family = "Helvetica") +
                        theme(axis.ticks=element_blank(),
                              axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0,
                                                         colour = colxaxis))
                
                plot_list[[k]] <- hmTr
                names(plot_list)[k] <- paste(fac_levels[i], "_vs_", fac_levels[j], sep = "")
        }
        plot_list
}
# --



# --
#######################################
### DESeq2Apply_physeq
#######################################
## Inputs
# physeq: phyloseq object
# group_var: name of column that defines group fac in sample_data
# SFs: often you might want to give the SizeFactors already because you wanted to calculate them on non-filtered data,
# when SFs are not NULL, type is ignored
# type: type in estimateSizeFactors, ignored when Size factors given
## OUTPUT:
# list of two lists: the first: List of DESeq2 results plus of fisher.exact test, for each level combination in group factor one data_frame = list entry
# in the second list is just the size factor adjusted physeq object


DESeq2Apply_physeq <- function(physeq, group_var, SFs = NULL, type = "ratio", p.adjust.method = "fdr", symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))){
        
        if (taxa_are_rows(physeq)) {physeq <- t(physeq)}
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        fac_levels <- levels(group_fac)
        
        # NB: DESeq2 can not deal with ordered factors, sees them somehow as just one level, therefore
        sample_data(physeq)[[group_var]] <- factor(sample_data(physeq)[[group_var]], levels = levels(sample_data(physeq)[[group_var]]), ordered = FALSE)
        
        DES = phyloseq::phyloseq_to_deseq2(physeq, formula(paste("~", group_var)))
        
        
        if (is.null(SFs)){
                if(type == "ratio"){
                        GM <- apply(otu_table(physeq), 2, gm_own, zeros.count = FALSE)
                }
                
                dds <- estimateSizeFactors(DES, type = type, geoMeans = GM)
                # NB: geoMeans is ignored when type = "iterate"
                # NB2: "iterate" takes much longer than "ratio", and when using GM via gm_own with "ratio" the size factors 
                # correlate with above 99% with size factors from "iterate"
                # SFs2 <- sizeFactors(dds)
                
        } else {
                dds <- DES
                sizeFactors(dds) = SFs
                # identical(sizeFactors(dds), SFs) 
        }
        
        
        dds <- estimateDispersions(dds, quiet = TRUE) 
        dds <- nbinomWaldTest(dds)
        
        # to get the size factor adjusted physeq object
        physeq_out <- physeq
        otu_table(physeq_out) <- otu_table(t(counts(dds, normalized = TRUE)), taxa_are_rows = FALSE)
        
        # - get the level combis -
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
        i_sAndj_s <- get_unique_facLevel_combis(fac_levels)
        i_s <- i_sAndj_s[["i_s"]]
        j_s <- i_sAndj_s[["j_s"]]
        # --
        
        result_list <- vector("list", length = length(i_s))
        
        for (k in seq_along(i_s)) {
                i <- i_s[k]
                j <- j_s[k]
                
                res <- as.data.frame(results(dds, contrast = c(group_var, fac_levels[i], fac_levels[j])))
                
                res$p_val_adj <- p.adjust(res$pvalue, method = p.adjust.method) # NB: in case of "fdr" same as default DESeq2
                
                CT <- counts(dds, normalized = TRUE)
                n1 <- sum(group_fac == fac_levels[i])
                n2 <- sum(group_fac == fac_levels[j])
                res$Median_grp1 <- apply(CT[, group_fac == fac_levels[i]], 1, median)
                res$Median_grp2 <- apply(CT[, group_fac == fac_levels[j]], 1, median)
                res$Mean_grp1 <- apply(CT[, group_fac == fac_levels[i]], 1, mean)
                res$Mean_grp2 <- apply(CT[, group_fac == fac_levels[j]], 1, mean)
                # res$baseMeanSelf <- apply(CT, 1, mean) # exactly the same as baseMean!
                res$Zeros_grp1 <- apply(CT[, group_fac == fac_levels[i]], 1, function(cnts){sum(cnts == 0)})
                res$Zeros_grp2 <- apply(CT[, group_fac == fac_levels[j]], 1, function(cnts){sum(cnts == 0)})
                res$Present_grp1 <- apply(CT[, group_fac == fac_levels[i]], 1, function(cnts){sum(cnts != 0)})
                res$Present_grp2 <- apply(CT[, group_fac == fac_levels[j]], 1, function(cnts){sum(cnts != 0)})
                res$prev_PC_grp1 <- round(100*(res$Present_grp1/n1),1)
                res$prev_PC_grp2 <- round(100*(res$Present_grp2/n2), 1)
                #res$Sparsity_grp1 <- 100*(res$Zeros_grp1/n1)
                #res$Sparsity_grp2 <- 100*(res$Zeros_grp2/n2)
                res$n1 <- n1
                res$n2 <- n2
                
                # - add fisher exact test of sparsity/prevalence again -
                Fisher <- t(sapply(1:nrow(res), FUN = function(i){
                        fisherMat <- matrix(c(res$Present_grp1[i], res$Zeros_grp1[i], res$Present_grp2[i],
                                              res$Zeros_grp2[i]), ncol = 2, dimnames = list(c("Present", "Absent"), c("grp1", "grp2")))
                        Test <- fisher.test(fisherMat)
                        cbind(Test$p.value, Test$estimate)
                }))
                
                p_val_adj_Fisher <- p.adjust(Fisher[,1], method = p.adjust.method)
                
                res$p_val_Fisher <- Fisher[,1]
                res$p_val_Fisher_adj <- p_val_adj_Fisher
                res$oddsRatioFisher <- Fisher[,2]
                # --
                
                # - add sginificance and direction -
                symnum.args$x <- res$pvalue
                res$signi. <- do.call(stats::symnum, symnum.args) %>% as.character()
                symnum.args$x <- res$p_val_adj
                res$signi_adj <- do.call(stats::symnum, symnum.args) %>% as.character()
                symnum.args$x <- res$p_val_Fisher
                res$signif_fisher <- do.call(stats::symnum, symnum.args) %>% as.character()
                symnum.args$x <- res$p_val_Fisher_adj
                res$signif_fisher_adj <- do.call(stats::symnum, symnum.args) %>% as.character()
                
                res$direction <- rep("down", nrow(res))
                res$direction[res$Median_grp1 > res$Median_grp2] <- "up"
                res$direction_fisher <- rep("down", nrow(res))
                res$direction_fisher[res$prev_PC_grp1 > res$prev_PC_grp2] <- "up"
                # --
                
                res$Taxon <- rownames(res)
                
                res <- dplyr::select(res, Taxon, teststat = stat, p_val = pvalue, p_val_adj,
                                     signi., signi_adj, direction, p_val_Fisher,
                                     p_val_Fisher_adj, signif_fisher, signif_fisher_adj,
                                     direction_fisher, Median_grp1, Median_grp2, Mean_grp1,
                                     Mean_grp2, Present_grp1, Present_grp2, Zeros_grp1, Zeros_grp2, prev_PC_grp1, prev_PC_grp2, 
                                     n1, n2, baseMean, log2FoldChange, lfcSE, oddsRatioFisher)
                # NB: I dropped here padj from DESeq since same as p_val_adj in case of p.adjust.method = "fdr"
                res <- cbind(res, tax_table(physeq))
                res <- dplyr::arrange(res, desc(abs(teststat)))
                # res <- res[order(res$p_val),]
                result_list[[k]] <- res
                names(result_list)[[k]] <- paste(fac_levels[i], "_vs_", fac_levels[j], sep = "")
        }
        
        out <- list(result_list, physeq_out)
        
}
# --



# --
#######################################
### make_heat_map_physeq##
#################
## Input:
# physeq object
# group_var: the name of the group_fac column in sample_data used to order the samples in the heat map
# max_abundance_for_color: if null the 90% quantile of count/relative abundance in the data is used. all counts above this value will be
# shown yellow in the heat map
# tax_order: character vector of the original taxon names in the order you want them shown. if NULL >> tax_order = taxa_names(physeq)
# tax_names: the names that will be used for the taxons, if Null Taxon_1, Taxon_2 and so on will be used. NB: renaming of course after
# ordering. 
# color_sample_names: if TRUE and if you have less than 7 levels, the sample names will be colored using cbPalette
# gradient_steps: the steps the blue to green to yellow will be distributed in the viridis gradient: 4 numbers ranging from 1e-14 to 1,
# see default, you might wanna try c(0.25, 0.5, 0.75, 1) as well
## Output:
# the heat map trellis

make_heat_map_physeq <- function(physeq, group_var, color_levels, max_abundance_for_color = NULL, tax_order = NULL,
                                 tax_names = NULL, color_sample_names = TRUE, gradient_steps = c(0.15, 0.3, 0.45, 1)){
        
        gradient_steps <- c(0, 1e-14, gradient_steps)
        
        if(taxa_are_rows(physeq)) { physeq <- t(physeq) }
        
        if (!is.factor(sample_data(physeq)[[group_var]])) {sample_data(physeq)[[group_var]] <- as.factor(sample_data(physeq)[[group_var]])}
        
        
        CT <- as(otu_table(physeq), "matrix") # taxa are columns
        DF_CT <- as.data.frame(t(CT))
        
        if (is.null(tax_order)){tax_order <- taxa_names(physeq)}
        
        if (length(tax_order) != nrow(DF_CT) || !all(tax_order %in% rownames(DF_CT))) {stop("the given tax_order must contain all taxa_names(physeq). If you want to subset,
                                                                                            do it before with prune_taxa")}
        DF_CT <- DF_CT[tax_order, ]
        
        if (is.null(tax_names)) {tax_names <- paste("Taxon_", 1:nrow(DF_CT), sep = "")}
        
        if (length(tax_names) != nrow(DF_CT) ) {stop("the given tax_names must be a character vector of length = ntaxa(physeq)")}
        
        
        rownames(DF_CT) <- make.unique(tax_names)
        DF_CT$Taxa <- rownames(DF_CT)
        DF_CT <- tidyr::gather(DF_CT, key = Sample , value = Count, -Taxa)
        DF_CT$Taxa <- factor(DF_CT$Taxa, levels = rev(unique(DF_CT$Taxa)), ordered = TRUE)
        LookUpDF <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
        LookUpDF <- LookUpDF[order(match(LookUpDF$Group, levels(LookUpDF$Group))), ]
        DF_CT$Sample <- factor(DF_CT$Sample, levels = LookUpDF$Sample, ordered = TRUE)
        
        if (is.null(max_abundance_for_color)) {
                max_abundance_for_color <- quantile(DF_CT$Count, .9)
        }
        
        if (max_abundance_for_color == 0) {max_abundance_for_color <- max(DF_CT$Count)}
        
        colxaxis <- color_levels[LookUpDF$Group]
        
        
        hmTr <- ggplot(DF_CT, aes(x = Sample, y = Taxa, fill = Count))
        hmTr <- hmTr + 
                geom_raster() + 
                scale_fill_gradientn("", limits = c(0, max_abundance_for_color), colors = c("red", viridis(5)), values = gradient_steps, oob = squish) +
                scale_x_discrete(position = "top") +
                #coord_equal() +
                labs(x=NULL, y=NULL) +
                theme_tufte(base_family = "Helvetica") +
                theme(axis.ticks=element_blank(),
                      axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0,
                                                 colour = colxaxis))
        hmTr
}
# --



# --
#######################################
### plot_toptaxa_boxAndviolin##
#################
# generates box and violin plots of the taxons in physeq in the order given by tax_order and named by tax_names
# if tax_order and tax_names are NULL the order and names of physeq will be used
# facet_cols determines ncol in facet_wrap
# color_levels must be a named vector of colors for all levels found in group_var
# if ttestp = "yes" ggpubr is used to add significance asterisk from simple t.tests

# OUTPUT:
# generates for each level combination of the grouping factor (defined by group_var) seven plots, so for each combi a list of 8 plots,
# specifically: boxplot, boxplot faceted, violin plot, violin plot faceted, and the same again for y axis log10 (NB: need to make new plots
# there because that needs pseydocounts for looking good)

plot_toptaxa_boxAndviolin <- function(physeq, group_var, tax_order = NULL, tax_names = NULL, facet_cols = 5, color_levels, ttestp = "yes"){
        
        if(taxa_are_rows(physeq)) { physeq <- t(physeq) }
        # I prefer taxa_are_rows = FALSE so rows (= observations = samples), and colums = variables = taxa
        
        if (!is.factor(sample_data(physeq)[[group_var]])) {sample_data(physeq)[[group_var]] <- as.factor(sample_data(physeq)[[group_var]])}
        
        CT <- as(otu_table(physeq), "matrix") # taxa are columns
        DF_CT <- as.data.frame(t(CT))
        
        if (is.null(tax_order)){tax_order <- taxa_names(physeq)}
        
        if (length(tax_order) != nrow(DF_CT) || !all(tax_order %in% rownames(DF_CT))) {stop("the given tax_order must contain all taxa_names(physeq). If you want to subset,
                                                                                            do it before with prune_taxa")}
        DF_CT <- DF_CT[tax_order, ]
        
        if (is.null(tax_names)) {tax_names <- paste("Taxon_", 1:nrow(DF_CT), sep = "")}
        
        if (length(tax_names) != nrow(DF_CT) ) {stop("the given tax_names must be a character vector of length = ntaxa(physeq)")}
        
        rownames(DF_CT) <- make.unique(tax_names)
        DF_CT$Taxa <- factor(rownames(DF_CT), levels = rownames(DF_CT), ordered = TRUE)
        
        DF_CT <- tidyr::gather(DF_CT, key = Sample , value = Count, -Taxa)
        
        LookUpDF <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
        DF_CT$Group <- LookUpDF$Group[match(DF_CT$Sample, LookUpDF$Sample)]
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        fac_levels <- levels(group_fac)
        
        # - get the level combis -
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
        i_sAndj_s <- get_unique_facLevel_combis(fac_levels)
        i_s <- i_sAndj_s[["i_s"]]
        j_s <- i_sAndj_s[["j_s"]]
        # --
        
        
        plot_list <- vector("list", length = length(i_s))
        
        for (k in seq_along(i_s)){
                i <- i_s[k]
                j <- j_s[k]
                
                group_fac_current <- droplevels(group_fac[as.numeric(group_fac) %in% c(i, j)])
                DF_CT_current <- filter(DF_CT, Group %in% group_fac_current)
                DF_CT_current$Group <- factor(DF_CT_current$Group, levels = c(fac_levels[i], fac_levels[j]), ordered = TRUE)
                
                Tr <- ggplot(DF_CT_current, aes(x = Taxa, y = Count, col = Group))
                Tr <- Tr +
                        geom_boxplot(outlier.color = NA) +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge()) +
                        scale_color_manual("", values = color_levels) +
                        xlab("") +
                        ylab("abundance") +
                        theme_bw() +
                        theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0),
                              legend.position = "none")
                
                if (ttestp == "yes"){
                        Tr <- Tr + ggpubr::stat_compare_means(label = "p.signif", method = "t.test") # t.test significance
                }
                
                
                Tr1 <- ggplot(DF_CT_current, aes(x = Group, y = Count, col = Group))
                Tr1 <- Tr1 +
                        geom_boxplot(outlier.color = NA) +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge()) +
                        facet_wrap(~ Taxa, ncol = facet_cols, scales = "free_y") +
                        scale_color_manual("", values = color_levels) +
                        xlab("") +
                        ylab("abundance") +
                        theme_bw() +
                        theme(legend.position = "none")
                
                if (ttestp == "yes"){
                        Tr1 <- Tr1 + ggpubr::stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.5) # t.test significance
                }
                
                
                Tr2 <- ggplot(DF_CT_current, aes(x = Taxa, y = Count, col = Group))
                Tr2 <- Tr2 +
                        geom_violin(fill = NA) +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge()) +
                        scale_color_manual("", values = color_levels) +
                        xlab("") +
                        ylab("abundance") +
                        theme_bw() +
                        theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0),
                              legend.position = "none")
                
                if (ttestp == "yes"){
                        Tr2 <- Tr2 + ggpubr::stat_compare_means(label = "p.signif", method = "t.test") # t.test significance
                }
                
                Tr3 <- ggplot(DF_CT_current, aes(x = Group, y = Count, col = Group))
                Tr3 <- Tr3 +
                        geom_violin(fill = NA) +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge()) +
                        facet_wrap(~ Taxa, ncol = facet_cols, scales = "free_y") +
                        scale_color_manual("", values = color_levels) +
                        xlab("") +
                        ylab("abundance") +
                        theme_bw() +
                        theme(legend.position = "none")
                
                if (ttestp == "yes"){
                        Tr3 <- Tr3 + ggpubr::stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.5) # t.test significance
                }
                
                
                
                # - add a pseudocount for log10 plots - # necessary for good boxplots
                DF_CT_current$Count[DF_CT_current$Count == 0] <- min(DF_CT_current$Count[DF_CT_current$Count > 0])
                # --
                
                Tr4 <- ggplot(DF_CT_current, aes(x = Taxa, y = Count, col = Group))
                Tr4 <- Tr4 +
                        geom_boxplot(outlier.color = NA) +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge()) +
                        scale_color_manual("", values = color_levels) +
                        xlab("") +
                        ylab("abundance") +
                        theme_bw() +
                        theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0),
                              legend.position = "none") +
                        scale_y_log10()
                
                if (ttestp == "yes"){
                        Tr4 <- Tr4 + ggpubr::stat_compare_means(label = "p.signif", method = "t.test") # t.test significance
                }
                
                
                Tr5 <- ggplot(DF_CT_current, aes(x = Group, y = Count, col = Group))
                Tr5 <- Tr5 +
                        geom_boxplot(outlier.color = NA) +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge()) +
                        facet_wrap(~ Taxa, ncol = facet_cols, scales = "free_y") +
                        scale_color_manual("", values = color_levels) +
                        xlab("") +
                        ylab("abundance") +
                        theme_bw() +
                        scale_y_log10() +
                        theme(legend.position = "none")
                
                if (ttestp == "yes"){
                        Tr5 <- Tr5 + ggpubr::stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.5) # t.test significance
                }
                
                
                
                Tr6 <- ggplot(DF_CT_current, aes(x = Taxa, y = Count, col = Group))
                Tr6 <- Tr6 +
                        geom_violin(fill = NA) +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge()) +
                        scale_color_manual("", values = color_levels) +
                        xlab("") +
                        ylab("abundance") +
                        theme_bw() +
                        theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0),
                              legend.position = "none") +
                        scale_y_log10()
                
                if (ttestp == "yes"){
                        Tr6 <- Tr6 + ggpubr::stat_compare_means(label = "p.signif", method = "t.test") # t.test significance
                }
                
                Tr7 <- ggplot(DF_CT_current, aes(x = Group, y = Count, col = Group))
                Tr7 <- Tr7 +
                        geom_violin(fill = NA) +
                        geom_point(size = 1, alpha = 0.6, position = position_jitterdodge()) +
                        facet_wrap(~ Taxa, ncol = facet_cols, scales = "free_y") +
                        scale_color_manual("", values = color_levels) +
                        xlab("") +
                        ylab("abundance") +
                        theme_bw() +
                        scale_y_log10() +
                        theme(legend.position = "none")
                
                if (ttestp == "yes"){
                        Tr7 <- Tr7 + ggpubr::stat_compare_means(label = "p.signif", method = "t.test", label.x = 1.5) # t.test significance
                }
                
                plot_list[[k]] <- list(Tr, Tr1, Tr2, Tr3, Tr4, Tr5, Tr6, Tr7)
                names(plot_list)[k] <- paste(fac_levels[i], "_vs_", fac_levels[j], sep = "")
        }
        
        plot_list
}
# --



# --
#######################################
### wilcoxTestApply_physeq
#######################################
# Accepts several groups in group_var and results in a list of DF for the pairwise comparisons
# directly related to wilcoxTestApply from SimulationFunctions.R just not looping through a list of phyloseq objects
# but doing the job on a single phyloseq object
# does wilcoxon test and adds further info
# in addition it performs a fisher exact test on the sparsity proportions
# NB: Tested: 
# with excludeZeros you can decide on whether 0 counts should be excluded for the wilcox.test, the fisher sparsity test is
# of course not affected. 
# NB: in case in one of the two groups all counts are 0 and excludeZeros = T, then NA is given for all wilcoxon 
# statistics!
# The teststatistic is based on the standardized teststatistic, equation provided by multtest::mt.minP (compare with mtApply)
# (see equation for standStat2 in the code, results in exactly the same as standStat when uncomment)
## Input
# physeq object
# group_var: refers to a factor in sample_data(phyloseq) that defines the groups
# excludeZeros: decides on whether 0s should be considered when comparing the groups in a wilcox.test
# p.adjust.method, used as method in p.adjust
## Output
# list of dataframes with the results for each pairwise group comparison in group_var associated group_fac

wilcoxTestApply_physeq <- function(physeq, group_var, excludeZeros = FALSE, p.adjust.method = "fdr", symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) {
        
        
        if (taxa_are_rows(physeq)) { physeq <- t(physeq)}
        
        CT <- as(otu_table(physeq), "matrix")
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        fac_levels <- levels(group_fac)
        
        # - get the level combis -
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels) 
        i_sAndj_s <- get_unique_facLevel_combis(fac_levels)
        i_s <- i_sAndj_s[["i_s"]]
        j_s <- i_sAndj_s[["j_s"]]
        # --
        
        result_list <- vector("list", length = length(i_s))
        
        for (k in seq_along(i_s)){
                i <- i_s[k]
                j <- j_s[k]
                
                res_mat <- apply(CT, 2, function(taxon_counts){
                        x <- taxon_counts[as.numeric(group_fac) == i]
                        Zeros_grp1 <- sum(x == 0)
                        # Sparsity_grp1 <- 100*(Zeros_grp1/length(x))
                        Present_grp1 <- length(x)-Zeros_grp1
                        prev_PC_grp1 <- 100*(Present_grp1/length(x))
                        if(excludeZeros){
                                x <- x[x != 0]
                        }
                        Median_grp1 <- median(x, na.rm = T) # NA in case all 0
                        Mean_grp1 <- mean(x, na.rm = T) # NaN in case all 0
                        if (is.na(Mean_grp1)){ Mean_grp1 = NA }
                        y <- taxon_counts[as.numeric(group_fac) == j]
                        Zeros_grp2 <- sum(y == 0)
                        Present_grp2 <- length(y)-Zeros_grp2
                        # Sparsity_grp2 <- 100*(Zeros_grp2/length(y))
                        prev_PC_grp2 <- 100*(Present_grp2/length(y))
                        if(excludeZeros){
                                #if(all(y == 0)){y[1] <- ceiling(mean(taxon_counts))+1}
                                y <- y[y != 0]
                        }
                        Median_grp2 <- median(y, na.rm = T)
                        Mean_grp2 <- mean(y, na.rm = T)
                        if (is.na(Mean_grp2)){ Mean_grp2 = NA }
                        
                        if (length(x) != 0 && length(y) != 0){
                                wilcTest <- wilcox.test(x = x, y = y, alternative = "two", paired = F, exact = F)
                                pValue <- wilcTest$p.value
                                W <- wilcTest$statistic
                                # calculate standardized rank sum Wilcoxon statistics as in multtest::mt.minP
                                Ranks <- rank(c(x, y))
                                n1 <- length(x)
                                n2 <- length(y)
                                # Wx <- sum(Ranks[1:n1])-(n1*(n1+1)/2) # would be the same as W
                                # how about the other W?
                                # Wy <- sum(Ranks[(n1+1):n2]) - (n2*(n2+1)/2)
                                standStat <- -1*((sum(Ranks[1:n1]) - n1*(n1+n2+1)/2)/sqrt(n1*n2*(n1+n2+1)/12))
                                
                                # # if you want to check that multtest::mt.minP would give the same statistic
                                # mati <- matrix(c(x,y), nrow = 1)
                                # grFac <- c(rep(fac_levels[i], n1), rep(fac_levels[j], n2))
                                # grFac <- factor(grFac, levels = c(fac_levels[i], fac_levels[j]))
                                # standStat2 <- multtest::mt.minP(mati, grFac, test = "wilcoxon")$teststat
                                # # identical(standStat, standStat2) # TRUE
                                # uncomment all with standStat2 to test all the way
                                
                        } else {
                                pValue = NA
                                W <- NA
                                standStat = NA
                                n1 <- length(x)
                                n2 <- length(y)
                                # standStat2 = NA
                        }
                        
                        
                        # -- add fisher exact test of presence differences (should be none in simulation) --
                        fisherMat <- matrix(c(Present_grp1, Zeros_grp1, Present_grp2, Zeros_grp2), 
                                            ncol = 2, dimnames = list(c("Present", "Absent"), c("grp1", "grp2")) )
                        Test <- fisher.test(fisherMat)
                        
                        c(teststat = standStat, p_val = pValue, Median_grp1 = Median_grp1, Median_grp2 = Median_grp2, 
                          Mean_grp1 = Mean_grp1, Mean_grp2 = Mean_grp2, n1 = n1, n2 = n2, Present_grp1 = Present_grp1, 
                          Present_grp2 = Present_grp2, Zeros_grp1 = Zeros_grp1, Zeros_grp2 = Zeros_grp2, 
                          prev_PC_grp1 = prev_PC_grp1, prev_PC_grp2 = prev_PC_grp2, W, 
                          p_val_Fisher = Test$p.value, Test$estimate) #, teststat2 = standStat2 # waste of typing effort:)
                })
                
                res_mat <- t(res_mat)
                colnames(res_mat) <- c("teststat", "p_val", "Median_grp1", "Median_grp2", "Mean_grp1", "Mean_grp2", "n1",
                                       "n2", "Present_grp1", "Present_grp2", "Zeros_grp1", "Zeros_grp2", "prev_PC_grp1", "prev_PC_grp2", "W",
                                       "p_val_Fisher", "oddsRatioFisher") #, , "teststat2"
                
                
                DF <- data.frame(Taxon = rownames(res_mat), res_mat)
                DF$p_val_adj <- p.adjust(DF$p_val, method = p.adjust.method)
                DF$p_val_Fisher_adj <- p.adjust(DF$p_val_Fisher, method = p.adjust.method)
                
                # - add sginificance and direction -
                symnum.args$x <- DF$p_val
                DF$signi. <- do.call(stats::symnum, symnum.args) %>% as.character()
                symnum.args$x <- DF$p_val_adj
                DF$signi_adj <- do.call(stats::symnum, symnum.args) %>% as.character()
                symnum.args$x <- DF$p_val_Fisher
                DF$signif_fisher <- do.call(stats::symnum, symnum.args) %>% as.character()
                symnum.args$x <- DF$p_val_Fisher_adj
                DF$signif_fisher_adj <- do.call(stats::symnum, symnum.args) %>% as.character()
                
                DF$direction <- rep("down", nrow(DF))
                DF$direction[DF$Median_grp1 > DF$Median_grp2] <- "up"
                DF$direction_fisher <- rep("down", nrow(DF))
                DF$direction_fisher[DF$prev_PC_grp1 > DF$prev_PC_grp2] <- "up"
                # --
                
                
                DF <- dplyr::select(DF, 1:3, p_val_adj, signi., signi_adj, direction, p_val_Fisher, p_val_Fisher_adj, signif_fisher,
                                    signif_fisher_adj, direction_fisher, 4:7, Present_grp1, Present_grp2, prev_PC_grp1, prev_PC_grp2)
                # teststat/standStat2 version:
                # DF <- dplyr::select(DF, 1:2, 19, 3, 20:22, 17, 23:25, 4:7, 10:15, 8:9, 16, 18)
                
                DF <- cbind(DF, tax_table(physeq))
                # DF <- dplyr::arrange(DF, desc(abs(teststat)))
                DF <- dplyr::arrange(DF, p_val)
                
                result_list[[k]] <- DF
                names(result_list)[[k]] <- paste(fac_levels[i], "_vs_", fac_levels[j], sep = "")
        }
        
        result_list
}
# --
