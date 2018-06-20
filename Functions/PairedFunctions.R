# --
#######################################
### calc_pVals_alphdiv
#######################################
# Output:
# a DF with the p-values of all pairwise comparisons between the levels in group_var for the different alpha diversity measures
# NB: uses ggpubr: compare_means

calc_pVals_alphdiv_paired <- function(DF_alpha, measures, group_var, compare, paired_var, test = "t.test", 
                               symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                               p.adjust.method = "BH", paired = TRUE){
        
        if(! group_var %in% colnames(DF_alpha)) {
                stop("The given group_var is not a variable in the DF_alpha.")
        }
        
        if(! paired_var %in% colnames(DF_alpha)) {
                stop("The given paired_var is not a variable in the DF_alpha.")
        }
        
        if (!all(compare %in% unique(DF_alpha[[group_var]]))) {
                stop("Not all given compare levels are actually levels in the group_var column.")
        }
        
        DF_alpha <- DF_alpha[DF_alpha[[group_var]] %in% compare,]
        DF_alpha[[group_var]] <- factor(DF_alpha[[group_var]], levels = compare, ordered = T)
        
        # - order samples based on group_var and paired_var -
        # kick out all samples that are not in exactly twice
        Twice <- names(table(DF_alpha[[paired_var]])[table(DF_alpha[[paired_var]]) == 2])
        DF_alpha <- DF_alpha[DF_alpha[[paired_var]] %in% Twice,]
        DF_alpha <- arrange_(DF_alpha, group_var, paired_var)
        # --
        
        
        y_columns <- which(colnames(DF_alpha) %in% c(measures, paste0(measures, "_resids", sep = "")))
        ttestList <- list()
        for (i in 1:length(y_columns)){
                comparison <- as.formula(paste(colnames(DF_alpha)[y_columns[i]], " ~ ", group_var, sep = ""))
                ttestList[[i]] <- ggpubr::compare_means(formula = comparison, data = DF_alpha, method = test, p.adjust.method = p.adjust.method, symnum.args = symnum.args, paired = paired)
        }
        alpha_div_ttests <- do.call("rbind", ttestList) %>% arrange(.y.)
}
# --






# --
#######################################
### boxplots_alphdiv_paired
#######################################
# Output:
# a list of boxplots for the alpha diversity measures given (including resid plots), the plots indicate the p-values between all levels in the group variable

# NB: makes use of ggpubr::stat_compare_means that allows you to adjust more settings, e.g. plotting symbols instead of numbers (change label to "p.signif")

boxplots_alphdiv_paired <- function(DF_alpha, measures, group_var, paired_var, shape, color_levels, test = "t.test", symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")), hide.ns = FALSE, paired = TRUE){
        
        if(! group_var %in% colnames(DF_alpha)) {
                stop("The given group_var is not a variable in the DF_alpha.")
        }
        
        if(! paired_var %in% colnames(DF_alpha)) {
                stop("The given paired_var is not a variable in the DF_alpha.")
        }
        
        if (!all(names(color_levels) %in% unique(DF_alpha[[group_var]]))) {
                stop("Not all given color levels are actually levels in the group_var column.")
        }
        
        DF_alpha <- DF_alpha[DF_alpha[[group_var]] %in% names(color_levels),]
        DF_alpha[[group_var]] <- factor(DF_alpha[[group_var]], levels = names(color_levels), ordered = T)
        
        # - order samples based on group_var and paired_var -
        # kick out all samples that are not in exactly twice
        Twice <- names(table(DF_alpha[[paired_var]])[table(DF_alpha[[paired_var]]) == 2])
        DF_alpha <- DF_alpha[DF_alpha[[paired_var]] %in% Twice,]
        DF_alpha <- arrange_(DF_alpha, group_var, paired_var)
        # --
        
        
        boxplotList <- list()
        fac_levels <- levels(DF_alpha[[group_var]])
        comparisonList <- get_unique_facLevel_combinations(fac_levels)
        
        y_columns <- which(colnames(DF_alpha) %in% c(measures, paste0(measures, "_resids", sep = "")))
        
        for (i in 1:length(y_columns)){
                aes_map <- aes_string(x = group_var, y = colnames(DF_alpha)[y_columns[i]], color = group_var)
                Tr <-  ggplot(DF_alpha, aes_map) + 
                        geom_boxplot(na.rm = TRUE, outlier.color = NA) +
                        geom_line(aes_string(group = paired_var), col = cbPalette[1]) +
                        geom_point() +
                        xlab("") +
                        scale_color_manual("", values = color_levels) +
                        theme_bw()
                Tr <- Tr + stat_compare_means(comparisons = comparisonList, method = test, label = "p.signif", symnum.args = symnum.args, hide.ns = hide.ns, paired = paired) # "p.signif"
                boxplotList[[i]] <- Tr
                names(boxplotList)[i] <- colnames(DF_alpha)[y_columns[i]]
        }
        boxplotList
}
# --











# --
#######################################
### boxplots_alphdiv_paired
#######################################
# Output:
# a list of boxplots for the alpha diversity measures given (including resid plots), the plots indicate the p-values between all levels in the group variable

# NB: makes use of ggpubr::stat_compare_means that allows you to adjust more settings, e.g. plotting symbols instead of numbers (change label to "p.signif")

boxplots_alphdiv_diff <- function(DF_alpha, measures, group_var, paired_var, color_levels){
        
        if(! group_var %in% colnames(DF_alpha)) {
                stop("The given group_var is not a variable in the DF_alpha.")
        }
        
        if(! paired_var %in% colnames(DF_alpha)) {
                stop("The given paired_var is not a variable in the DF_alpha.")
        }
        
        if (!all(names(color_levels) %in% unique(DF_alpha[[group_var]]))) {
                stop("Not all given color levels are actually levels in the group_var column.")
        }
        
        if (length(names(color_levels)) != 2) {
                stop("For this paired comparison only two groups allowed.")
        }
        
        
        level1 <- names(color_levels)[1]
        level2 <- names(color_levels)[2]
        
        
        DF_alpha <- DF_alpha[DF_alpha[[group_var]] %in% names(color_levels),]
        DF_alpha[[group_var]] <- factor(DF_alpha[[group_var]], levels = names(color_levels), ordered = T)
        
        
        # - order samples based on group_var and paired_var -
        # kick out all samples that are not in exactly twice
        Twice <- names(table(DF_alpha[[paired_var]])[table(DF_alpha[[paired_var]]) == 2])
        DF_alpha <- DF_alpha[DF_alpha[[paired_var]] %in% Twice,]
        DF_alpha <- arrange_(DF_alpha, group_var, paired_var)
        # --
        
        
        
        # - calculate differences - 
        # keep only relevant columns
        DF_diff <- DF_alpha[, colnames(DF_alpha) %in% c(measures, paste0(measures, "_resids", sep = ""), group_var, paired_var, shape)]
        y_columns <- which(colnames(DF_diff) %in% c(measures, paste0(measures, "_resids", sep = "")))
        # since order is equal, calculate Diffes easily
        First <- DF_diff[DF_diff[[group_var]] == level1, y_columns]
        Second <- DF_diff[DF_diff[[group_var]] == level2, y_columns]
        DF_Diff <- Second-First
        
        DF_Diff$Diff <- paste("diff. ", level1, " to ", level2, sep = "") 

        # --
        
        
        boxplotList <- list()
        
        y_columns <- which(colnames(DF_Diff) %in% c(measures, paste0(measures, "_resids", sep = "")))
        
        for (i in 1:length(y_columns)){
                aes_map <- aes_string(x = "Diff", y = colnames(DF_Diff)[y_columns[i]])
                Tr <-  ggplot(DF_Diff, aes_map) + 
                        geom_boxplot(na.rm = TRUE, outlier.color = NA) +
                        geom_jitter(width = .2, height = 0, alpha = 0.65, col = cbPalette[2]) +
                        xlab("") +
                        theme_bw()
                tt <- t.test(DF_Diff[[y_columns[i]]])
                Tr <- Tr +
                        ggtitle(paste("p.value t.test: ", format(tt$p.value, digits = 5), sep = ""))
                boxplotList[[i]] <- Tr
                names(boxplotList)[i] <- colnames(DF_alpha)[y_columns[i]]
        }
        boxplotList
}
# --









# --
#######################################
### plot_hittaxa_boxAndviolin_paired ##
#################


plot_hittaxa_boxAndviolin_paired <- function(physeq, group_var, color_levels, paired_var, taxa_info_df = NULL, taxa_annotation = NULL, 
                                      facet_cols = 5, shape = NULL, excludeZero = FALSE, logTransform = FALSE){
        
        # - make sure group_var and shape are in the sample_data -
        if (!group_var %in% colnames(sample_data(physeq))) {
                stop("group_var must be a variable in sample_data(physeq)")
        }
        
        if(! paired_var %in% colnames(sample_data(physeq))) {
                stop("The given paired_var is not a variable in sample_data(physeq).")
        }
        
        if (!is.null(shape) && !shape %in% colnames(sample_data(physeq))) {
                shape = NULL
        }
        # --
        
        
        # - show only samples "represented" in color levels -
        if (!all(unique(sample_data(physeq)[[group_var]]) %in% names(color_levels))) {
                keepSamples <- sample_names(physeq)[sample_data(physeq)[[group_var]] %in% names(color_levels)]
                physeq <- prune_samples(samples = keepSamples, physeq)
        }
        
        
        # make sure the order of the group levels is as defined by color_levels
        sample_data(physeq)[[group_var]] <- factor(sample_data(physeq)[[group_var]], levels = names(color_levels), order = TRUE)
        # --
        
        # - keep only taxa in taxa_info_df in physeq -
        if (!is.null(taxa_info_df)) {
                
                if (!all(rownames(taxa_info_df) %in% taxa_names(physeq))) {
                        stop("not all taxa in taxa_info_df are taxa in physeq!")
                }
                pruned_ps <- prune_taxa(rownames(taxa_info_df), physeq)
        } else {
                pruned_ps <- physeq
        }
        # --
        
        # - check or set taxa_annotation -
        if (is.null(taxa_annotation)){
                if (!is.null(taxa_info_df)) {
                        taxa_annotation <- paste(taxa_info_df$Taxon, "_", taxa_info_df$p_val_adj, sep = "")
                } else {
                        taxa_annotation <- taxa_names(pruned_ps) # or you could use 
                        # taxa_annotation <- get_taxon_names(as.data.frame(tax_table(pruned_ps)))
                }
        }
        
        if (length(taxa_annotation) != ntaxa(pruned_ps)) {
                warning("taxa_annotation did not fit in length to nrow(taxa_info_df) or ntaxa(physeq), used taxa_names")
                taxa_annotation <- taxa_names(pruned_ps)
        }
        
        taxa_annotation <- as.character(taxa_annotation)
        
        taxa_annotation[is.na(taxa_annotation)] <- "NA"
        
        taxa_annotation <- make.unique(taxa_annotation)
        
        if (!is.null(taxa_info_df)) {
                taxa_lookup <- data.frame(Taxon = taxa_info_df$Taxon, Annotation = taxa_annotation)
        } else {
                taxa_lookup <- data.frame(Taxon = taxa_names(pruned_ps), Annotation = taxa_annotation)
        }
        # --
        
        # - generate data frame for plotting, add annotation, and make sure order of annotation fits -
        DF_CT <- psmelt(pruned_ps)
        
        DF_CT$Annotation <- taxa_lookup$Annotation[match(DF_CT$OTU, taxa_lookup$Taxon)]
        DF_CT$Annotation <- factor(DF_CT$Annotation, levels = taxa_annotation, ordered = TRUE) #
        # --
        
        # - exclude zeros and log transform if asked for -
        if (excludeZero) {
                DF_CT <- filter(DF_CT, Abundance != 0)
        }
        
        if (logTransform) {
                
                minValue <- min(DF_CT$Abundance)
                
                if (minValue < 0){
                        stop("you asked for log_transform but the lowest count in your data is already below 0!")
                }
                if (minValue == 0) {
                        pseudocount <- min(DF_CT$Abundance[DF_CT$Abundance > 0])/2
                        DF_CT$Abundance[DF_CT$Abundance == 0] <- pseudocount
                }
                DF_CT$Abundance <- log10(DF_CT$Abundance)
                
        }
        # --
        
        
        
        Tr <- ggplot(DF_CT, aes_string(x = group_var, y = "Abundance", col = group_var))
        Tr <- Tr +
                geom_boxplot(outlier.color = NA, na.rm = TRUE) +
                geom_line(aes_string(group = paired_var), col = cbPalette[1]) +
                geom_point(aes_string(shape = shape), size = 1, alpha = 0.6) +
                facet_wrap(~ Annotation, ncol = facet_cols, scales = "free_y") +
                scale_color_manual("", values = color_levels) +
                xlab("") +
                ylab("abundance") +
                theme_bw()
        
        if (is.null(shape)) {
                Tr <- Tr + theme(legend.position = "none")
        }
        
        
        Tr1 <- ggplot(DF_CT, aes_string(x = group_var, y = "Abundance", col = group_var))
        Tr1 <- Tr1 +
                geom_violin(fill = NA) +
                geom_line(aes_string(group = paired_var), col = cbPalette[1]) +
                geom_point(aes_string(shape = shape), size = 1, alpha = 0.6) +
                facet_wrap(~ Annotation, ncol = facet_cols, scales = "free_y") +
                scale_color_manual("", values = color_levels) +
                xlab("") +
                ylab("abundance") +
                theme_bw()
        
        if (is.null(shape)) {
                Tr1 <- Tr1 + theme(legend.position = "none")
        }
        
        
        
        list(Boxplot = Tr, ViolinPlot = Tr1)   
}             

# --






# --
####################################
## plot_taxa_ratios_paired 
###################################
# see plot_taxa_ratios_levelPairs: Here you directly calculate the count by count ratio matrix only for the tax_nom, you still facet by taxa_den
# (denominator) but you keep all levels in all plots. Makes extensive use of ggpubr, NB: ggpubr is so smart to adjust p-values when you use
# scale_y_log10

## Output: 
# - list(pVals = pVals, Tr = Tr, Tr1 = Tr1, pValsLog = pValsLog, Tr2 = Tr2, Tr3 = Tr3)


plot_taxa_ratios_paired <- function(physeq, group_var, color_levels, paired_var, tax_names = NULL,
                                       taxa_nom = "Firmicutes", taxa_den = NULL, test = "t.test", p_adjust_method = "fdr",
                                       tax_order = NULL,
                                       symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) {
        
        
        if (taxa_are_rows(physeq)) { 
                physeq <- t(physeq)
        }
        
        if(! group_var %in% colnames(sample_data(physeq))) {
                stop("The given group_var is not a variable in the sample data of the loaded phyloseq object.")
        }
        
        if (!all(names(color_levels) %in% unique(sample_data(physeq)[[group_var]]))) {
                stop("Not all names in names(color_levels)are actually levels in the group_var column.")
        }
        
        
        if(! paired_var %in% colnames(sample_data(physeq))) {
                stop("The given paired_var is not a variable in sample_data(physeq).")
        }
        
        if (!(test %in% c("t.test", "wilcox.test"))) {
                stop("test should be t.test or wilcox.test")
        }
        
        group_fac <- factor(sample_data(physeq)[[group_var]])
        
        compare <- names(color_levels)
        
        if (!is.null(compare)) {
                group_var_levels <- compare
        } else {
                group_var_levels <- levels(group_fac)
        }
        
        if (length(group_var_levels) != 2) {
                stop(paste0("sorry, compare (group_var_levels) must consist of two groups - you asked for ",
                            paste(group_var_levels, collapse = ", ")))
        }
        
        
        # - check that given tax_names fit to physeq and change taxa_names of physeq -
        if (is.null(tax_names)){
                tax_names <- paste("T", 1:ntaxa(physeq), sep = "_")
        } 
        
        if(!identical(ntaxa(physeq), length(tax_names))){stop("tax_names do not fit in length to physeq")}
        
        tax_names <- make.unique(tax_names)
        taxa_names(physeq) <- tax_names
        # --
        
        # - calculate the matrix taxa_nom/(all other taxa) -         
        CT <- t(as(otu_table(physeq), 'matrix')) # now taxa are rows and samples are columns
        # ONLY keep samples defined by group_var_levels (= names(color_levels))
        CT <- CT[, group_fac %in% group_var_levels]
        
        
        i <- which(rownames(CT) == taxa_nom)
        if (length(i) != 1) {stop("taxa_nom not found in tax_names or tax_names not unique!")}
        
        
        TbTmatrix <- apply(CT, 2, function(samp_cnts){samp_cnts[i]/samp_cnts})
        # produces for each taxon (= host taxon) a TbTMatrix
        # NB: there are possibly Inf, and NaN values in the matrix, specifically
        # 0/x = 0, x/0 = Inf; 0/0 = NaN!
        # --
        
        TbT_DF <- as.data.frame(TbTmatrix)
        TbT_DF$Taxon <- rownames(TbT_DF)
        # - use taxa_den (denominator) to restrict the taxa to which taxa_nom is compared to -
        if (is.null(taxa_den)) {taxa_den <- tax_names}
        TbT_DF <- TbT_DF[TbT_DF$Taxon %in% taxa_den, ]
        # --
        
        # - change to long DF -
        TbT_DF_l <- gather(TbT_DF, key = Sample, value = Ratio, -Taxon)
        # --
        
        # - add the sample_data information  -
        SS <- as(sample_data(physeq), "data.frame")
        SS$Sample <- rownames(SS)
        TbT_DF_l <- merge(TbT_DF_l, SS, by = "Sample")
        # --
        
        # - change all ratios where either nominator taxon or denominator taxon had count = 0 to NA -
        # remember: 0/0 = NaN (not finite), 0/x = 0, x/0 = Inf
        TbT_DF_l$Ratio[!is.finite(TbT_DF_l$Ratio) | TbT_DF_l$Ratio == 0] <- NA
        # --
        
        # - kick out all NA cases -
        TbT_DF_l <- TbT_DF_l[!is.na(TbT_DF_l$Ratio), ]
        
        # - calculate paired p.values -
        allTaxa <- unique(TbT_DF_l$Taxon)
        pVals <- vector(mode = "numeric", length = length(allTaxa))
        pValslog10 <- vector(mode = "numeric", length = length(allTaxa))
        for (i in 1:length(allTaxa)){
                currentTaxon <- allTaxa[i]
                currentDF <- filter(TbT_DF_l, Taxon == currentTaxon)
                # keep only samples present in both levels
                currentDF <- currentDF[currentDF[[paired_var]] %in% names(table(currentDF[[paired_var]])[table(currentDF[[paired_var]]) == 2]),]
                # order
                currentDF <- arrange_(currentDF, group_var, paired_var)
                
                if (nrow(currentDF) > 3){
                        if (test == "t.test") {
                                Tester <- t.test(x = currentDF$Ratio[currentDF[[group_var]] == group_var_levels[1]], y = currentDF$Ratio[currentDF[[group_var]] == group_var_levels[2]], paired = TRUE)
                                Testerlog10 <- t.test(x = log10(currentDF$Ratio[currentDF[[group_var]] == group_var_levels[1]]), y = log10(currentDF$Ratio[currentDF[[group_var]] == group_var_levels[2]]), paired = TRUE)
                        } else if (test == "wilcox.test"){
                                Tester <- wilcox.test(x = currentDF$Ratio[currentDF[[group_var]] == group_var_levels[1]], y = currentDF$Ratio[currentDF[[group_var]] == group_var_levels[2]], paired = TRUE)
                                Testerlog10 <- wilcox.test(x = log10(currentDF$Ratio[currentDF[[group_var]] == group_var_levels[1]]), y = log10(currentDF$Ratio[currentDF[[group_var]] == group_var_levels[2]]), paired = TRUE)
                        } else {
                                stop("test must be t.test or wilcox.test")
                        }
                        pVals[i] <- Tester$p.value
                        pValslog10[i] <- Testerlog10$p.value
                } else {
                        pVals[i] <- NA
                        pValslog10[i] <- NA
                }
                
                
        }
        
        pVals <- data.frame(Taxon = allTaxa, p = pVals, p.adj = p.adjust(pVals, method = p_adjust_method)) %>% arrange(p)
        symnum.args$x <- pVals$p.adj
        pVals$signi <- do.call(stats::symnum, symnum.args) %>% as.character()
        pValslog10 <- data.frame(Taxon = allTaxa, p = pValslog10, p.adj = p.adjust(pValslog10, method = p_adjust_method)) %>% arrange(p)
        symnum.args$x <- pValslog10$p.adj
        pValslog10$signi <- do.call(stats::symnum, symnum.args) %>% as.character()
        # --
        
        
        
        # - order taxa based on pVals result or based on tax_order -
        if (is.null(tax_order)){
                TbT_DF_l$Taxon <- factor(TbT_DF_l$Taxon, levels = unique(pVals$Taxon), ordered = TRUE)
        } else {
                if(!all(unique(pVals$Taxon) %in% tax_order)){
                        stop("given tax_order does not fit to tax_names")
                }
                TbT_DF_l$Taxon <- factor(TbT_DF_l$Taxon, levels = tax_order, ordered = TRUE)
        }
        # --
        

        Tr <- ggplot(TbT_DF_l, aes_string(x = group_var, y = "Ratio", col = group_var))
        Tr <- Tr +
                geom_violin() +
                geom_line(aes_string(group = paired_var), col = cbPalette[1]) +
                geom_point(aes_string(shape = shape), size = 1, alpha = 0.6) +
                facet_wrap(~ Taxon,  scales = "free_y") +
                scale_color_manual("", values = color_levels) +
                xlab("") +
                ylab(paste("abundance ratio of", taxa_nom, "to stated taxon")) +
                theme_bw()
        
        if (is.null(shape)) {
                Tr <- Tr + theme(legend.position = "none")
        }

        
        Tr1 <- ggplot(TbT_DF_l, aes_string(x = group_var, y = "Ratio", col = group_var))
        Tr1 <- Tr1 +
                geom_boxplot(outlier.color = NA, na.rm = TRUE) +
                geom_line(aes_string(group = paired_var), col = cbPalette[1]) +
                geom_point(aes_string(shape = shape), size = 1, alpha = 0.6) +
                facet_wrap(~ Taxon,  scales = "free_y") +
                scale_color_manual("", values = color_levels) +
                xlab("") +
                ylab(paste("abundance ratio of", taxa_nom, "to stated taxon")) +
                theme_bw()
        
        if (is.null(shape)) {
                Tr1 <- Tr1 + theme(legend.position = "none")
        }
        
        Tr2 <- Tr + scale_y_log10()
        
        Tr3 <- Tr1 + scale_y_log10()
        
        
        list(pVals = pVals, Tr = Tr, Tr1 = Tr1, pValsLog = pValslog10, Tr2 = Tr2, Tr3 = Tr3)
        
}
#--



