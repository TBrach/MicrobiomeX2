

# --
#######################################
#### visualize_filtering
#######################################


visualize_filtering <- function(physeq, prevalence, taxa_sums_quantile, col_vec = NULL){ 
        
        df_ab_prev <- data.frame(Taxon_No = 1:ntaxa(physeq), 
                                 total_counts = taxa_sums(physeq),
                                 prevalence = colSums(as(otu_table(physeq), "matrix") != 0),
                                 sparsity = colSums(as(otu_table(physeq), "matrix") == 0),
                                 mean_count_nonzero = apply(as(otu_table(physeq), "matrix"), 2, function(x){mean(x[x > 0])}),
                                 median_count_nonzero = apply(as(otu_table(physeq), "matrix"), 2, function(x){median(x[x > 0])})
        )
        
        df_ab_prev <- cbind(df_ab_prev, tax_table(physeq))
        
        prev_thresh <- (prevalence/100)*nsamples(physeq)
        abund_thresh <- quantile(taxa_sums(physeq), probs = taxa_sums_quantile/100)
        
        df_ab_prev_filt <- dplyr::filter(df_ab_prev, prevalence > prev_thresh | total_counts > abund_thresh)
        
        no_samples <- nsamples(physeq)
        shade_df <- data.frame(total_counts = 0, prevalence = 0)
        
        
        Tr_prev_vs_log10_ab <- ggplot(df_ab_prev, aes(x = total_counts, y = prevalence))
        Tr_prev_vs_log10_ab <- Tr_prev_vs_log10_ab +
                geom_point(col = cbPalette[2], size = 2, alpha = 0.7) +
                scale_x_log10() +
                geom_rect(data = shade_df, xmin = -Inf, xmax = log10(abund_thresh), ymin = -Inf, ymax = prev_thresh, fill = "#660033", alpha = 0.4) +
                geom_hline(yintercept = (prevalence/100)*nsamples(physeq), col = cbPalette[1], lty = "dashed") +
                geom_vline(xintercept = quantile(taxa_sums(physeq), probs = taxa_sums_quantile/100), col = cbPalette[1], lty = "dashed") +
                xlab("total counts (taxa_sums())") + 
                theme_bw() +
                ggtitle(paste(nrow(df_ab_prev_filt), " of ", nrow(df_ab_prev), " taxa (", round(100*nrow(df_ab_prev_filt)/nrow(df_ab_prev), 1),
                              " %) and ", round(sum(df_ab_prev_filt$total_counts)), " of ", round(sum(df_ab_prev$total_counts)), " counts (",
                              round((sum(df_ab_prev_filt$total_counts)/sum(df_ab_prev$total_counts))*100, 1), " %) remain", sep = ""))
        
        
        df_ab_prev2 <- df_ab_prev
        # - adjust color and order of the phyla in the following plot - 
        col <- "Phylum"
        
        if (!is.null(col_vec)){
                df_ab_prev2[[col]] <- as.character(df_ab_prev2[[col]])
                df_ab_prev2[[col]][is.na(df_ab_prev2[[col]])] <- "NA"
                if (length(col_vec) != length(unique(df_ab_prev2[[col]]))){
                        stop("provided col_vec did not fit in length to the col variable")
                }
                df_ab_prev2[[col]] <- factor(df_ab_prev2[[col]], levels = names(col_vec), ordered = TRUE)
                custom_colors <- col_vec
        } else {
                CountOrder <- dplyr::group_by_(df_ab_prev2, col) %>% dplyr::summarise(total_count_sum = sum(total_counts)) %>% dplyr::arrange(desc(total_count_sum))
                
                CountOrder[[col]] <- as.character(CountOrder[[col]])
                CountOrder[[col]][is.na(CountOrder[[col]])] <- "NA"
                
                if (nrow(CountOrder) <= 15){
                        custom_colors <- make_color_vector(CountOrder[[col]], QuantColors15)
                } else {
                        custom_colors <- make_color_vector(CountOrder[[col]], viridis(nrow(CountOrder)))
                }
                
                df_ab_prev2[[col]] <- as.character(df_ab_prev2[[col]])
                df_ab_prev2[[col]][is.na(df_ab_prev2[[col]])] <- "NA"
                df_ab_prev2[[col]] <- factor(df_ab_prev2[[col]], levels = names(custom_colors), ordered = TRUE)
                
        }
        # --
        
        
        
        Tr_prev_vs_log10_ab_col <- ggplot(df_ab_prev2, aes(x = total_counts, y = prevalence))
        Tr_prev_vs_log10_ab_col <- Tr_prev_vs_log10_ab_col +
                geom_point(aes(col = Phylum), size = 2, alpha = 0.7) +
                scale_color_manual(values = custom_colors) +
                scale_x_log10() +
                geom_rect(data = shade_df, xmin = -Inf, xmax = log10(abund_thresh), ymin = -Inf, ymax = prev_thresh, fill = "#660033", alpha = 0.4) +
                geom_hline(yintercept = (prevalence/100)*nsamples(physeq), col = cbPalette[1], lty = "dashed") +
                geom_vline(xintercept = quantile(taxa_sums(physeq), probs = taxa_sums_quantile/100), col = cbPalette[1], lty = "dashed") +
                xlab("total counts (taxa_sums())") +
                facet_wrap(~Phylum) +
                theme_bw() +
                theme(legend.position = "none")
        
        
        
        out <- list(Tr_prev_vs_log10_ab = Tr_prev_vs_log10_ab,
                    Tr_prev_vs_log10_ab_col = Tr_prev_vs_log10_ab_col)
        
}
# --



# --
#######################################
### plot_sizeFactors
#######################################
# just to check if size factors differ between groups and how they are related with original sample sizes

plot_sizeFactors <- function(physeq, SFs, group_var, color_levels, shape, test = "t.test", p_adjust_method = "fdr",
                             symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                             hide.ns = FALSE){
        
        
        # in case you do not want to see all samples
        if (!all(unique(sample_data(physeq)[[group_var]]) %in% names(color_levels))) {
                keepSamples <- sample_names(physeq)[sample_data(physeq)[[group_var]] %in% names(color_levels)]
                physeq <- prune_samples(samples = keepSamples, physeq)
                sample_data(physeq)[[group_var]] <- factor(sample_data(physeq)[[group_var]], levels = names(color_levels), order = TRUE)
                SFs <- SFs[names(SFs) %in% keepSamples]
        }
        
        DF <- cbind(sample_data(physeq), SF = SFs)
        
        Tr <-  ggplot(DF, aes_string(x = group_var, y = "SF", color = group_var)) 
        Tr <- Tr + 
                geom_boxplot(outlier.color = NA) +
                geom_jitter(aes_string(shape = shape), width = .2, height = 0, alpha = 0.65) +
                xlab("") +
                ylab("size factor") +
                scale_color_manual("", values = color_levels) +
                theme_bw()
        if(is.null(shape)){
                Tr <- Tr + theme(legend.position = "none")
        }
        
        # - since you might have more than two levels in each plot you need to set the comparisons argument in stat_compare_means -
        group_fac <- factor(sample_data(physeq)[[group_var]])
        fac_levels <- levels(group_fac)
        
        comparisonList <- get_unique_facLevel_combinations(fac_levels)
        
        Tr <- Tr + ggpubr::stat_compare_means(comparisons = comparisonList, label = "p.signif", method = test, hide.ns = hide.ns)
        
        formulaa <- as.formula(paste("SF ~", group_var, sep = " "))
        
        pVals <- compare_means(formula = formulaa, data = DF, method = test, p.adjust.method = p_adjust_method, symnum.args = symnum.args)
        
        # - Add total_counts vs SFs plot -
        DF$total_count <- sample_sums(physeq)
        Tr1 <-  ggplot(DF, aes_string(x = "total_count", y = "SF", color = group_var)) 
        Tr1 <- Tr1 + 
                geom_point(aes_string(shape = shape), alpha = 0.65) +
                xlab("sample_sums()") +
                ylab("size factor") +
                scale_color_manual("", values = color_levels) +
                theme_bw()
        if(is.null(shape)){
                Tr1 <- Tr1 + theme(legend.position = "none")
        }
        
        
        
        list(pVals = pVals, Tr = Tr, Tr1 = Tr1)
}
# --









