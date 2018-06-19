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




