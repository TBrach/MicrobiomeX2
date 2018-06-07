# --
#######################################
### calc_alphadiv_plusLmResids
#######################################
# Output:
# outlist: 
#    [[1]]: DF of alpha_diversity values, plus residuals of the values to a linear fit against total counts (sample_sums)
#    [[2]]: list of the fit objects derived from lm of alpha diversity values against sample_sums()

# NB: is based on estimate_richness function from phyloseq
# You could easily calculate Shannon, Chao1, Observed self, see alphaDiversityMeasures.Rmd
# estimate_richness itself uses functions from the vegan package

calc_alphadiv_plusLmResids <- function(physeq, measures = c("Observed", "Shannon")) {
        
        DF_alpha <- suppressWarnings(phyloseq::estimate_richness(physeq, measures = measures))
        
        DF_alpha$Total <- sample_sums(physeq)
        
        rownames(DF_alpha) <- sample_names(ps)
        # because linear fits of alpha diversity measures to total counts are often highly significant, I add the residuals of these
        # linear fits. 

        fitlist <- list()
        ncol_df_alpha <- ncol(DF_alpha)
        
        for (i in 1:length(measures)){
                fit <- lm(DF_alpha[,measures[i]] ~ DF_alpha[,"Total"])
                fitlist[[i]] <- fit
                DF_alpha[, ncol_df_alpha + i] <- residuals(fit)
                colnames(DF_alpha)[ncol_df_alpha + i] <- paste(measures[i], "resids", sep = "_")
        }
        
        names(fitlist) <- measures
        
        
        DF_alpha <- data.frame(DF_alpha, sample_data(physeq))
        
        outlist <- list(DF_alpha = DF_alpha, fitlist = fitlist)
}
# --



# --
#######################################
### calc_pVals_alphdiv
#######################################
# Output:
# a DF with the p-values of all pairwise comparisons between the levels in group_var for the different alpha diversity measures
# NB: uses ggpubr: compare_means

calc_pVals_alphdiv <- function(DF_alpha, measures, group, test = "t.test", 
                               symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                               p.adjust.method = "BH"){
        
        y_columns <- which(colnames(DF_alpha) %in% c(measures, paste0(measures, "_resids", sep = "")))
        ttestList <- list()
        for (i in 1:length(y_columns)){
                comparison <- as.formula(paste(colnames(DF_alpha)[y_columns[i]], " ~ ", group, sep = ""))
                ttestList[[i]] <- ggpubr::compare_means(formula = comparison, data = DF_alpha, method = test, p.adjust.method = p.adjust.method, symnum.args = symnum.args)
        }
        alpha_div_ttests <- do.call("rbind", ttestList) %>% arrange(.y.)
}
# --



# --
#######################################
### boxplots_alphdiv
#######################################
# Output:
# a list of boxplots for the alpha diversity measures given (including resid plots), the plots indicate the p-values between all levels in the group variable

# NB: makes use of ggpubr::stat_compare_means that allows you to adjust more settings, e.g. plotting symbols instead of numbers (change label to "p.signif")

boxplots_alphdiv <- function(DF_alpha, measures, group, shape, color_levels, test = "t.test", hide.ns = FALSE){
        
        boxplotList <- list()
        fac_levels <- levels(DF_alpha[[group]])
        
        y_columns <- which(colnames(DF_alpha) %in% c(measures, paste0(measures, "_resids", sep = "")))
        
        # - get the level combis -
        i_sAndj_s <- get_unique_facLevel_combis(fac_levels)
        
        i_s <- i_sAndj_s[["i_s"]]
        j_s <- i_sAndj_s[["j_s"]]
        # --
        
        
        comparisonList <- list()
        for (k in seq_along(i_s)){
                comparisonList[[k]] <- c(fac_levels[i_s[k]], fac_levels[j_s[k]])
        }
        
        for (i in 1:length(y_columns)){
                aes_map <- aes_string(x = group, y = colnames(DF_alpha)[y_columns[i]], color = group)
                Tr <-  ggplot(DF_alpha, aes_map) + 
                        geom_boxplot(na.rm = TRUE, outlier.color = NA) +
                        geom_jitter(aes_string(shape = shape), width = .2, height = 0, alpha = 0.65) +
                        xlab("") +
                        scale_color_manual("", values = color_levels) +
                        theme_bw()
                Tr <- Tr + stat_compare_means(comparisons = comparisonList, method = "t.test", label = "p.signif", hide.ns = hide.ns) # "p.signif"
                boxplotList[[i]] <- Tr
                names(boxplotList)[i] <- colnames(DF_alpha)[y_columns[i]]
        }
        boxplotList
}
# --



# --
#######################################
### lmPlots_alphdiv
#######################################
# Output:
# a list of linear fit plots for the alpha diversity measures given against Total = sample_sums()

# NB: requires the lmp function!

lmPlots_alphdiv <- function(DF_alpha, lm_fitlist, measures, group, shape, color_levels, test = "t.test"){
        
        lm_plotList <- list()
        for (i in 1:length(measures)){
                aes_map <- aes_string(x = "Total", y = measures[i], color = group, shape = shape)
                Tr <-  ggplot(DF_alpha, aes_map) + 
                        geom_point(na.rm = TRUE) +
                        xlab("total counts (sample_sums())") +
                        scale_color_manual("", values = color_levels) +
                        theme_bw()
                fit <- lm_fitlist[[measures[i]]]
                adjR2 <- round(summary(fit)$adj.r.squared,3)
                if (adjR2 != 0){
                        pfit <- lmp(fit)
                } else {
                        pfit <- NA
                }
                Tr <- Tr + geom_smooth(aes_string(x = "Total", y = measures[i]), method = "lm", se = TRUE, inherit.aes = F)
                
                Tr <- Tr + ggtitle(paste("lm fit: p-value: ", format(pfit, digits = 4), ", adj.r2: ", adjR2, sep = ""))
                lm_plotList[[i]] <- Tr
                names(lm_plotList)[i] <- measures[i]
                
        }
        lm_plotList
}
# --



# --
#######################################
### raref_curve_richness
#######################################
# read rarefaction_curve_own, this fast version uses the vegan::rarefy function, so it only generates rarefaction
# for richness. NB: vegan::rarefy used here does averaging (read help)


raref_curve_richness <- function(physeq, group_var = NULL, max_total = NULL, step_size = 200, color_levels, seed = 123) {
        
        
        if (taxa_are_rows(physeq)) {
                seqtab <- t(as(otu_table(physeq), "matrix"))
        } else {
                seqtab <- as(otu_table(physeq), "matrix")
        }
        
        
        if (is.null(max_total)) {
                max_total <- quantile(sample_sums(physeq), probs = .25)
        }
        
        steps <- seq(from = 0, to = max_total, by = step_size)
        NoSamples <- nsamples(physeq)
        NoSteps <- length(steps)
        richness_matrix <- matrix(data = NA, nrow = NoSamples, ncol = NoSteps)
        totalamplicons <- sample_sums(physeq)
        
        set.seed(seed)
        
        # ptm <- Sys.time()
        for (i in 1:length(steps)) {
                richness_matrix[,i] <- vegan::rarefy(seqtab, sample = steps[i], se = FALSE)
        }
        # Sys.time() - ptm
        
        # set richness for samples at which steps > totalamplicons to NA
        for (i in 1:NoSamples) {
                NaIndex <- which(totalamplicons[i] < steps)[1]
                if (!is.na(NaIndex)){
                        richness_matrix[i, NaIndex:ncol(richness_matrix)] <- NA
                }
        }
        
        rownames(richness_matrix) <- rownames(seqtab)
        colnames(richness_matrix) <- paste("step_", steps, sep = "")
        richness_df <- as.data.frame(richness_matrix)
        
        
        plot_div_df3 <- function (div_df, type = "richness") {
                
                div_df$Sample <- rownames(div_df)
                div_df$Total <- totalamplicons
                div_df <- tidyr::gather(div_df, key = step, value = div, -Sample, -Total)
                div_df <- tidyr::separate(div_df, col = step, into = c("term", "step"), sep = "_")
                div_df$step <- as.numeric(div_df$step)
                
                Tr <- ggplot(div_df, aes(x = step, y = div, group = Sample, col = Total))
                Tr <- Tr +
                        geom_line() +
                        scale_color_gradient2("total counts bef.", low = cbPalette[6], mid = cbPalette[1], high = cbPalette[2], midpoint = median(div_df$Total)) +
                        xlab("total counts in sample") +
                        ylab(type) +
                        theme_bw()
        }
        
        Tr_richness_grad <- plot_div_df3(richness_df, type = "richness")
        
        
        
        if (!is.null(group_var)){
                Group <- sample_data(physeq)[[group_var]] 
                if (!is.null(Group) && !is.factor(Group)) {
                        Group <- as.factor(Group)
                }
        } else { Group <- NULL}
        
        
        if (!is.null(Group)) {
                
                richness_df$Group <- Group
                
                # pairwise.tt_richness <- lapply(richness_df[, -ncol(richness_df)], function(step){ 
                #         ptt <- pairwise.t.test(x = step, g = richness_df$Group, alternative = "two", p.adjust.method = "none", var.equal = F, pool.sd = F)
                #         ptt$p.value
                # })
                
                pairwise.tt_richness <- NA # NB: I currently do not even show these p-values and it is currently a source of error as soon as you do not
                # have enough samples within a group at a certain step for a t.test (i.e. only 1 sample in each group left.)
                
                plot_div_df2 <- function (div_df, type = "richness") {
                        
                        div_df$Sample <- rownames(div_df)
                        div_df <- tidyr::gather(div_df, key = step, value = div, -Sample, -Group)
                        div_df <- tidyr::separate(div_df, col = step, into = c("term", "step"), sep = "_")
                        div_df$step <- as.numeric(div_df$step)
                        
                        Tr <- ggplot(div_df, aes(x = step, y = div, group = Sample, col = Group))
                        Tr <- Tr +
                                geom_line() +
                                xlab("total counts in sample") +
                                scale_color_manual("", values = color_levels) +
                                ylab(type) +
                                theme_bw(12)
                }
                
                plot_div_df_group <- function (div_df, type = "alpha diversity") {
                        div_df <- tidyr::gather(div_df, key = step, value = div, - Group)
                        div_df <- tidyr::separate(div_df, col = step, into = c("term", "step"), sep = "_")
                        div_df$step <- as.numeric(div_df$step)
                        div_df <- group_by(div_df, Group, step)
                        div_df <- dplyr::summarise(div_df, Mean_div = mean(div, na.rm = T), SD_div = sd(div, na.rm = T), n = n(), SE_div = SD_div/sqrt(n))
                        
                        Tr <- ggplot(div_df, aes(x = step, y = Mean_div, col = Group))
                        Tr <- Tr + 
                                geom_line() +
                                geom_point(size = 1) +
                                geom_errorbar(aes(ymin = Mean_div-SE_div, ymax = Mean_div+SE_div)) +
                                scale_color_manual("", values = color_levels) +
                                xlab("total counts in sample") +
                                ylab(paste(type, "(group mean +- SE)")) +
                                theme_bw(12)
                }
                Tr_richness_col <- plot_div_df2(richness_df, type = "richness")
                
                Tr_richness_group <- plot_div_df_group(richness_df, type = "richness")
                
        } else {
                Tr_richness_col <- NA
                Tr_richness_group <- NA
                pairwise.tt_richness <- NA
        }
        
        outlist <- list(rarefaction_richness_df = richness_df, Tr_richness_grad = Tr_richness_grad, Tr_richness_col = Tr_richness_col,
                        Tr_richness_group = Tr_richness_group, pairwise.tt_richness = pairwise.tt_richness)
        
}
# --

