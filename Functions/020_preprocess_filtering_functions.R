# --
#######################################
### calc_SFs
#######################################
# implementation of DESeq2 library size similar to estimateSizeFactorsForMatrix
# in difference to adjust_LS (obsolete) ignore.zero.ratios is always TRUE (so 0 ratios are always ignored)


calc_SFs <- function(physeq, zeros.count = FALSE, percentile = 50)  {
        
        # ---- Step 1: Calculate Geometric mean for each taxa over all samples ------
        # NB: these GM is basically the reference sample
        if(taxa_are_rows(physeq)){
                GM <- apply(otu_table(physeq), 1, gm_own, zeros.count = zeros.count)   
        } else {
                GM <- apply(otu_table(physeq), 2, gm_own, zeros.count = zeros.count) 
        }
        
        # ---- Step 2: Calculate Size factors --------
        
        # NB: x/y = exp(log(x) - log(y))
        
        if (taxa_are_rows(physeq)) {
                SFs <- apply(as(otu_table(physeq), "matrix"), 2, function(sample_cnts){exp(quantile((log(sample_cnts)-log(GM))[sample_cnts > 0], probs = percentile/100, na.rm = T))})
                SFs <- SFs/exp(mean(log(SFs)))
        } else {
                SFs <- apply(as(otu_table(physeq), "matrix"), 1, function(sample_cnts){exp(quantile((log(sample_cnts)-log(GM))[sample_cnts > 0], probs = percentile/100, na.rm = T))})
                SFs <- SFs/exp(mean(log(SFs)))
        } 
        
        
        if(min(SFs) == 0) { warning("in at least one sample the Size Factor was 0!") }
        
        SFs
}
# --



# --
#######################################
### simply_adjust_LS
#######################################
# implementation of DESeq2 library size similar to estimateSizeFactorsForMatrix
# in difference to adjust_LS (obsolete) ignore.zero.ratios is always TRUE (so 0 ratios are always ignored)

## Input:
# physeq
# zeros.count: zeros.count of gm_own, if TRUE 0 will be considered when calculating the geometric mean
# if FALSE not and thus the geometric means will be bigger (see gm_own)
# percentile: the percentile to select the size factor SF of a sample based on its count ratios to the reference sample. In
# DESeq percentile = 50, i.e. stats::median is used. 
# plots: if TRUE SFs and plots will be given in an output list, otherwise just the adjusted physeq is returned


simply_adjust_LS <- function(physeq, SFs = NULL, zeros.count = FALSE, percentile = 50)  {
        
        # ---- Step 1: Calculate Geometric mean for each taxa over all samples ------
        # NB: these GM is basically the reference sample
        if(taxa_are_rows(physeq)){
                GM <- apply(otu_table(physeq), 1, gm_own, zeros.count = zeros.count)   
        } else {
                GM <- apply(otu_table(physeq), 2, gm_own, zeros.count = zeros.count) 
        }
        
        # ---- Step 2: Calculate Size factors  unless given --------
        
        # NB: x/y = exp(log(x) - log(y))
        
        if (is.null(SFs)){
                
                if (taxa_are_rows(physeq)) {
                        SFs <- apply(as(otu_table(physeq), "matrix"), 2, function(sample_cnts){exp(quantile((log(sample_cnts)-log(GM))[sample_cnts > 0], probs = percentile/100, na.rm = T))})
                        SFs <- SFs/exp(mean(log(SFs)))
                } else {
                        SFs <- apply(as(otu_table(physeq), "matrix"), 1, function(sample_cnts){exp(quantile((log(sample_cnts)-log(GM))[sample_cnts > 0], probs = percentile/100, na.rm = T))})
                        SFs <- SFs/exp(mean(log(SFs)))
                }
                
        }
        
        
        if(min(SFs) == 0) { warning("in at least one sample the Size Factor was 0!") }
        
        
        # --- 3: calculate the new counts and put into a physeq object
        
        if(taxa_are_rows(physeq)){
                if (!identical(colnames(otu_table(physeq)), names(SFs))) {stop("names SFs do not fit to physeq")}
                Mat <- sweep(otu_table(physeq),2,SFs, "/")
                phynew <- physeq
                otu_table(phynew) <- otu_table(Mat, taxa_are_rows = TRUE)
        } else {
                if (!identical(rownames(otu_table(physeq)), names(SFs))) {stop("names SFs do not fit to physeq")}
                Mat <- sweep(otu_table(physeq),1,SFs, "/") 
                phynew <- physeq
                otu_table(phynew) <- otu_table(Mat, taxa_are_rows = FALSE)
        }
        
        
                
        list(physeq = phynew, SFs = SFs)

}
# --



# --
#######################################
#### visualize_filtering
#######################################


visualize_filtering <- function(physeq, prevalence, taxa_sums_quantile){ 
        
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
        CountOrder <- dplyr::group_by_(df_ab_prev2, col) %>% dplyr::summarise(total_count_sum = sum(total_counts)) %>% dplyr::arrange(desc(total_count_sum))
        
        CountOrder[[col]] <- as.character(CountOrder[[col]])
        CountOrder[[col]][is.na(CountOrder[[col]])] <- "NA"
        df_ab_prev2[[col]] <- as.character(df_ab_prev2[[col]])
        df_ab_prev2[[col]][is.na(df_ab_prev2[[col]])] <- "NA"
        df_ab_prev2[[col]] <- factor(df_ab_prev2[[col]], levels = CountOrder[[col]], ordered = TRUE)
        
        if (nrow(CountOrder) <= 15) {
                custom_colors <- QuantColors15[1:nrow(CountOrder)]
                
        } else {
                custom_colors <- viridis(nrow(CountOrder))
                
        }
        
        names(custom_colors) <- levels(df_ab_prev2[[col]])
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



# --
#######################################
#### plot_sample_bars_compare
#######################################
# generates abundance barplots (see plot_bar_own) to compare ps to ps_tca, i.e. to see how SFs adjustment affects the abundances 

plot_sample_bars_compare <- function(physeq, physeq2, x = "Sample", y = "Abundance", group_var, color_levels, fill = NULL,
                                 color_sample_names = TRUE){
        
        # - prepare mdf of ps physeq -
        if(taxa_are_rows(physeq)) { physeq <- t(physeq) }
        
        if (!is.factor(sample_data(physeq)[[group_var]])) {sample_data(physeq)[[group_var]] <- as.factor(sample_data(physeq)[[group_var]])}
        
        if (is.null(fill)) { fill = "Phylum"}
        
        mdf <- phyloseq::psmelt(physeq)
        
        # order samples according to levels
        LookUpDF <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
        LookUpDF <- LookUpDF[order(match(LookUpDF$Group, levels(LookUpDF$Group))), ]
        mdf$Sample <- factor(mdf$Sample, levels = LookUpDF$Sample, ordered = TRUE)
        # mdf$Sample <- factor(mdf$Sample, levels = c("A-15A", "A-5A", "A-2A", "A-1A", "B-15A", "B-5A", "B-2A", "B-1A"), ordered = TRUE)
        
        # order fill levels according to abundance over all samples
        mdf[, fill] <- as.character(mdf[, fill])
        mdf[is.na(mdf[, fill]), fill] <- "NA"
        sums <- group_by_(mdf, fill) %>% summarise(sum_abundance = sum(Abundance)) %>% arrange(sum_abundance)
        mdf[, fill] <- factor(mdf[, fill], levels = as.character(sums[[1]]), ordered = TRUE)
        
        # --
        
        # - prepare mdf of ps_tca physeq -
        if(taxa_are_rows(physeq2)) { physeq2 <- t(physeq2) }
        
        if (!is.factor(sample_data(physeq2)[[group_var]])) {sample_data(physeq2)[[group_var]] <- as.factor(sample_data(physeq2)[[group_var]])}
        
        if (is.null(fill)) { fill = "Phylum"}
        
        mdf2 <- phyloseq::psmelt(physeq2)
        
        # order samples according to levels
        LookUpDF <- data.frame(Sample = sample_names(physeq2), Group = sample_data(physeq2)[[group_var]])
        LookUpDF <- LookUpDF[order(match(LookUpDF$Group, levels(LookUpDF$Group))), ]
        mdf2$Sample <- factor(mdf2$Sample, levels = LookUpDF$Sample, ordered = TRUE)
        # mdf2$Sample <- factor(mdf2$Sample, levels = c("A-15A", "A-5A", "A-2A", "A-1A", "B-15A", "B-5A", "B-2A", "B-1A"), ordered = TRUE)
        
        # order fill levels according to abundance over all samples
        mdf2[, fill] <- as.character(mdf2[, fill])
        mdf2[is.na(mdf2[, fill]), fill] <- "NA"
        sums <- group_by_(mdf2, fill) %>% summarise(sum_abundance = sum(Abundance)) %>% arrange(sum_abundance)
        mdf2[, fill] <- factor(mdf2[, fill], levels = as.character(sums[[1]]), ordered = TRUE)
        # --
        
        mdf$Typer <- "before"
        mdf2$Typer <- "after SF adjustment"
        
        mdf <- rbind(mdf, mdf2)
        mdf$Typer <- factor(mdf$Typer, levels = c("before", "after SF adjustment"), ordered = TRUE)
        
        # - define names of x axis using color_levels (which must be a named character vector) -        
        colxaxis <- color_levels[LookUpDF$Group]
        # --
        
        
        
        if (length(levels(mdf[, fill])) <= 15) {
                fill_colors <- QuantColors15[1:length(levels(mdf[, fill]))]
                names(fill_colors) <- rev(levels(mdf[, fill]))
        } else {
                fill_colors <- rev(viridis(length(levels(mdf[, fill]))))
                names(fill_colors) <- rev(levels(mdf[, fill]))
        }
        
        
        Tr <- ggplot(mdf, aes_string(x = x, y = y, fill = fill))
        Tr <- Tr + 
                geom_bar(stat = "identity", position = "stack") +
                theme_bw() +
                scale_fill_manual(values = fill_colors) +
                xlab("") +
                facet_wrap(~ Typer, ncol = 1) +
                theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0, colour = colxaxis))
        
        Tr
}
# --





