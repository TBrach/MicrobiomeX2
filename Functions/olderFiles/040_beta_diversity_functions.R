# --
#######################################
### FUNCTION: calc_beta_div_distances
#######################################

calc_beta_div_distances <- function(physeq, dist_methods = c("bray")) {
        
        dist_list <- vector("list", length(dist_methods))
        names(dist_list) = dist_methods
        
        for (i in dist_methods) {
                iDist <- phyloseq::distance(physeq, method=i)
                dist_list[[i]] = iDist
        }
        
        return(dist_list)
        
}
# --



# --
######################################
### FUNCTION: compare_beta_div_distances_directly
#######################################
# INPUT:
# dist_list: named list of dist objects (so for each distance one object)
# physeq: phyloseq object used to make the dist objects
# group_var: character identifying the grouping variable in the sample_data of the phyloseq
# OUTPUT:
# list of three named lists, one with the plots and one with the data_frames of the p_values from pairwise t.tests, and one (important if more than two levels in group_var) where
# the levels in group_var are compared individually always within_gr1 vs between vs within_gr2

compare_beta_div_distances_directly <- function(dist_list, physeq, group_var, test = "t.test", symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),
                               p.adjust.method = "BH", hide.ns = FALSE) {
        
        TrList <- vector(mode = "list", length = length(dist_list))
        pValListAll <- vector(mode = "list", length = length(dist_list))
        pValListDirect <- vector(mode = "list", length = length(dist_list)) # only within to between group comparisons for each level combination
        # for example: you have 3 levels in group_var, in Direct you do not compare distances 1 to 2 vs distances 3 to 3. In ALL you do. 
        
        for (i in 1:length(dist_list)){
                DistMat <- as(dist_list[[i]], "matrix")
                # - Transform the lower triangel of Distmatrix into a data frame: Sample1 (Row), Sample2 (Col), Distance -
                rowCol <- expand.grid(rownames(DistMat), colnames(DistMat))
                labs <- rowCol[as.vector(lower.tri(DistMat, diag=F)),]
                df <- cbind(labs, DistMat[lower.tri(DistMat, diag=F)])
                colnames(df) <- c("Row","Col","Distance")
                # NB: we exclude here distances of samples with themselves, therefore for x samples: (x-1)*x/2 distances = nrow(df)
                # --
                
                # - add for each distance the group_var levels the samples are from -
                samdf <- data.frame(Sample = sample_names(physeq), Group = sample_data(physeq)[[group_var]])
                df$Row_Group <- samdf$Group[match(df$Row, samdf$Sample)]
                df$Col_Group <- samdf$Group[match(df$Col, samdf$Sample)]
                # -- 
                
                
                # - add GroupToGroup using matching to levels (and ordering! important whne there are more than 2 levels!) to make sure that for example "Old to Young" and "Young to Old" both become "Young vs Old"
                df$GroupToGroup <- apply(dplyr::select(df, Row_Group, Col_Group), 1, function(x){
                        paste(x[order(match(x, levels(samdf$Group)))], collapse = " to ")
                })
                # --

                
                # - change GroupToGroup into an ordered factor so in the plots later you see always distances within first level, distances between first and
                # second level, distances within second level - 
                df$GroupToGroupOrder <- apply(df[, c("Row_Group", "Col_Group")], 1, function(x){
                        x1 <- match(x[1], levels(samdf$Group))
                        x2 <- match(x[2], levels(samdf$Group))
                        if (x1 <= x2) {
                                as.numeric(paste(x1, x2, sep = ""))
                        } else {
                                as.numeric(paste(x2, x1, sep = ""))
                        }
                })
                
                df$Type <- "between"
                df$Type[df$Row_Group == df$Col_Group] <- "within"
                # NB: for within group distances, you have samples*(samples-1)/2 distances, for between group distances
                # you have samplesgrp1 * samplesgrp2 distances.
                
                GToGLeveldf <- unique(df[, c("GroupToGroup", "GroupToGroupOrder", "Type")])
                GroupToGroupLevels <- GToGLeveldf[order(GToGLeveldf$GroupToGroupOrder), ]$GroupToGroup
                GToGLeveldfBetween <- GToGLeveldf[GToGLeveldf$Type == "between", ] # needed for facetting the plots later, see below
                GroupToGroupLevelsBetween <- GToGLeveldfBetween[order(GToGLeveldfBetween$GroupToGroupOrder), ]$GroupToGroup
                df$GroupToGroup <- factor(df$GroupToGroup, levels = GroupToGroupLevels, ordered = TRUE)
                # --
                
                # - use compare_means to calculate p-values for ALL pairwise GroupToGroup distance comparisons --
                df_p <- ggpubr::compare_means(formula = Distance ~ GroupToGroup, data = df, method = test, p.adjust.method = p.adjust.method, symnum.args = symnum.args)
                pValListAll[[i]] <- df_p
                # --
                
                # in case there are more than two levels in group_var I want a faceted plot, i.e. one facet for each GroupToGroupLevelsBetween, 
                # for this I need to duplicate the within data
                
                if (length(GroupToGroupLevelsBetween) > 1) {
                        # - generate new df with required data duplicated -
                        df_list <- lapply(GroupToGroupLevelsBetween, function(level){
                                grps <- unlist(strsplit(level, " to "))
                                df_current <- df[df$Row_Group %in% grps & df$Col_Group %in% grps, ]
                                df_current$Level <- level
                                df_current
                        })
                        df_plot <- do.call("rbind", df_list)
                        # --
                        
                        
                        df_plot$Level <- factor(df_plot$Level, levels = GroupToGroupLevelsBetween, ordered = TRUE)
                        
                        df_pDirect <- ggpubr::compare_means(formula = Distance ~ GroupToGroup, data = df_plot, group.by = "Level", method = test, p.adjust.method = p.adjust.method, symnum.args = symnum.args)
                        pValListDirect[[i]] <- df_pDirect
                        
                        # This worked fine for calculating the p-values but it did not work with Tr + stat_compare_means, therefore other way below
                        # Tr <- ggplot(df_plot, aes(x = GroupToGroup, y = Distance, col = Type))
                        # Tr <- Tr + geom_boxplot(outlier.color = NA) + 
                        #         geom_jitter(position = position_jitter(width = 0.3, height = 0), alpha = 0.65) +
                        #         facet_wrap(~ Level, ncol = 2, scales = "free_x") +
                        #         xlab("") +
                        #         ylab(paste(names(dist_list)[i], "distance", sep = " ")) +
                        #         # scale_color_manual("", values = c("within" = cbPalette[6], "between" = cbPalette[7])) +
                        #         theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                        #               legend.position = "none")
                        # 
                        # Tr + stat_compare_means(label = "p.signif", method = test, group.by = "Level")
                        
                        # - generate a new variable Type2 so you have always the same groups (within1 between within2) for each factet (needed for stat_compare_menas) -
                        df_plot$Type2 <- df_plot$Type
                        df_plot <- separate(df_plot, col = Level, into = c("Level1", "Level2"), sep = " to ", remove = FALSE)
                        df_plot$Type2[df_plot$Type == "within"] <- "within_gr1"
                        df_plot$Type2[df_plot$Type == "within" & (df_plot$Row_Group == df_plot$Level2)] <- "within_gr2"
                        df_plot$Type2 <- factor(df_plot$Type2, levels = c("within_gr1", "between", "within_gr2"), ordered = TRUE)
                        # --
                        
                        Tr <- ggplot(df_plot, aes(x = Type2, y = Distance, col = Type))
                        Tr <- Tr + geom_boxplot(outlier.color = NA)
                        if (nsamples(physeq) < 100) {
                                Tr <- Tr + geom_jitter(position = position_jitter(width = 0.3, height = 0), alpha = 0.65)
                        }
                        Tr <- Tr +
                                facet_wrap(~ Level, ncol = 2, scales = "free_x") +
                                xlab("") +
                                ylab(paste(names(dist_list)[i], "distance", sep = " ")) +
                                scale_color_manual("", values = c("within" = cbPalette[6], "between" = cbPalette[7])) +
                                theme_bw() +
                                theme(legend.position = "none")
                        
                        my_comparisons <- list(c("within_gr1", "between"), c("within_gr1", "within_gr2"), c("within_gr2", "between"))
                        
                        Tr <- Tr + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = test, hide.ns = hide.ns) # alternative: "p.format"
                
                        
                } else {
                        
                        Tr <- ggplot(df, aes(x = GroupToGroup, y = Distance, col = Type))
                        Tr <- Tr + geom_boxplot(outlier.color = NA)
                        if (nsamples(physeq) < 100) {
                                Tr <- Tr + geom_jitter(position = position_jitter(width = 0.3, height = 0), alpha = 0.65)
                        }
                        Tr <- Tr +
                                xlab("") +
                                ylab(paste(names(dist_list)[i], "distance", sep = " ")) + 
                                scale_color_manual("", values = c("within" = cbPalette[6], "between" = cbPalette[7])) +
                                theme_bw() +
                                theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                                      legend.position = "none")
                        
                        # there is probably a less clumsy way to get to the my_comparisons
                        comparisonList <- get_unique_facLevel_combinations(GroupToGroupLevels)
                        
                        Tr <- Tr + stat_compare_means(comparisons = comparisonList, label = "p.signif", method = test, hide.ns = hide.ns)
                        
                        pValListDirect[[i]] <- df_p
                        
                }
                
                TrList[[i]] <- Tr
        }
        
        names(TrList) <- names(dist_list)
        names(pValListAll) <- names(dist_list)
        names(pValListDirect) <- names(dist_list)
        out <- list(DistanceBoxplots = TrList, DistancePValuesAll = pValListAll, DistancePValuesDirect = pValListDirect)
        
}
# --



# --
#######################################
### FUNCTION: loop_vegan_adonis
#######################################
# Function is very much based on pairwise.perm.manova {RVAideMemoire}
# but it also records R2 while looping through vegan::adonis, and generates a
# result data frame in which the results are shown in the order of the group_fac levels
# INPUT:
# dist_obj: dist object
# group_fac: the factor that groups the samples
# nperm: permutations in adonis
# p.adj.method: method to adjust p.values
# symnum.args: symbols linked to p-value cutpoints
# OUTPUT:
# data.frame showing p.values and R2 and adjusted p.values for the different between group comparisons


loop_vegan_adonis <- function(dist_obj, group_fac, nperm = 999, 
                                     p.adj.method = "none", symnum.args = list(cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) {
        
        if (!("dist" %in% class(dist_obj))){
                stop("dist_obj must be of class dist")
        }
        
        group_fac <- factor(group_fac)
        
        # - get all level combinations within group_fac - # changed this so that I think order is more logical, check regularly if correct
        fac_levels <- levels(group_fac)
        fac_levels_num <- setNames(seq_along(fac_levels), fac_levels)
        i_sAndj_s <- get_unique_facLevel_combis(fac_levels)
        i_s <- i_sAndj_s[["i_s"]]
        j_s <- i_sAndj_s[["j_s"]]
        # -- 
        
        p_vals <- vector(mode = "numeric", length = length(i_s))
        r2s <- vector(mode = "numeric", length = length(i_s))
        
        # - loop vegan::adonis for all level combinations and generate df -
        for (k in seq_along(i_s)){
                i <- i_s[k]
                j <- j_s[k]
                group_fac2 <- droplevels(group_fac[as.numeric(group_fac) %in% c(i, j)])
                dist_obj_mat <- as.matrix(dist_obj)
                rows <- which(group_fac %in% levels(group_fac2))
                dist_obj2 <- as.dist(dist_obj_mat[rows, rows])
                fit <- vegan::adonis(dist_obj2 ~ group_fac2, permutations = nperm)
                p_vals[k] <- fit$aov.tab[1, "Pr(>F)"]
                r2s[k] <- fit$aov.tab[1, "R2"]
        }
        
        
        
        result_df <- data.frame(Comparison = paste0(names(fac_levels_num[i_s]), "_to_", names(fac_levels_num[j_s]), sep = ""),
                                adonis_p_value = p_vals, adonis_R2 = r2s, p_val_adj = p.adjust(p_vals, p.adj.method))
        # --
        
        # - add vegan::adonis including all group levels -
        fit <- vegan::adonis(dist_obj ~ group_fac, permutations = nperm)
        df <- data.frame(Comparison = "Overall",
                         adonis_p_value = fit$aov.tab[1, "Pr(>F)"], 
                         adonis_R2 = fit$aov.tab[1, "R2"], 
                         p_val_adj = fit$aov.tab[1, "Pr(>F)"])
        
        result_df <- rbind(df, result_df)
        # -- 
        
        # - add significance levels -
        symnum.args$x <- result_df$adonis_p_value
        result_df$significance <- do.call(stats::symnum, symnum.args) %>% as.character()
        result_df <- mutate(result_df, adonis_R2_PC = round(100*adonis_R2, 2))
        # --
        
}
# --



# --
#######################################
### FUNCTION: calc_ordination_from_distances
#######################################
# fully based on phyloseq::ordinate and phyloseq::plot_ordination


calc_ordination_from_distances <- function(physeq, dist_list, color_levels, ordination_type = "PCoA", group_var = NULL, shape = NULL, coord_cor = FALSE){
        
        ordination_list <- vector("list", length(dist_list))
        DFList <- vector("list", length(dist_list))
        DF_taxa_List <- vector("list", length(dist_list))
        # TrList <- vector("list", length(dist_list))
        TrList_own <- vector("list", length(dist_list))
        TrList_taxa <- vector("list", length(dist_list))
        
        axes <- 1:2 # currently only allowed to plot first and second
        
        for (i in seq_along(dist_list)) {
                
                ordination <- phyloseq::ordinate(physeq, method = ordination_type, distance = dist_list[[i]])
                ordination_list[[i]] <- ordination
                DF <- phyloseq::plot_ordination(physeq, ordination_list[[i]], color = group_var, justDF = TRUE)
                DFList[[i]] <- DF # just the first two axes cbind to sample_data in physeq
                
                x = colnames(DF)[1]
                y = colnames(DF)[2]
                Tr <- ggplot(DF, aes_string(x = x, y = y, col = group_var, shape = shape)) 
                Tr <- Tr + geom_point() + 
                        scale_color_manual("", values = color_levels) +
                        theme_bw(12) +
                        ggtitle(names(dist_list)[i])
                
                # for labelling axes
                if (ordination_type == "PCoA" || ordination_type == "NMDS") {
                        if (ordination_type == "PCoA") {
                                eigvec <- phyloseq:::extract_eigenvalue.pcoa(ordination)
                        } else {
                                eigvec <- phyloseq:::extract_eigenvalue.default(ordination)
                        }
                        
                        if (length(eigvec[axes]) > 0){
                                fracvar = eigvec[axes]/sum(eigvec)
                                percvar = round(100 * fracvar, 1)
                                strivar = as(c(Tr$label$x, Tr$label$y), "character")
                                strivar = paste0(strivar, " (", percvar, " %)")
                                Tr <- Tr + xlab(strivar[1]) + ylab(strivar[2]) 
                        }
                        
                        if (!is.null(eigvec) && coord_cor) {
                                Tr <- Tr + coord_fixed(sqrt(eigvec[2] / eigvec[1]))
                        }
                        
                }
                
                TrList_own[[i]] <- Tr
                rm(Tr)
                
                
                # TrList[[i]] <- phyloseq::plot_ordination(physeq, ordination_list[[i]], color = group_var) + ggtitle(names(dist_list)[i])
                
                DF_taxa <- phyloseq::plot_ordination(physeq, ordination_list[[i]], type = "taxa", color = "Phylum", justDF = TRUE)
                DF_taxa_List[[i]] <- DF_taxa
                
                x = colnames(DF_taxa)[1]
                y = colnames(DF_taxa)[2]
                Tr <- ggplot(DF_taxa, aes_string(x = x, y = y, col = "Phylum")) 
                Tr <- Tr + geom_point() +
                        # scale_color_manual("", values = cbPalette[2:8]) +
                        theme_bw(12) +
                        ggtitle(names(dist_list)[i])
                
                # for labelling axes
                if (ordination_type == "PCoA" || ordination_type == "NMDS") {
                        if (ordination_type == "PCoA") {
                                eigvec <- phyloseq:::extract_eigenvalue.pcoa(ordination)
                        } else {
                                eigvec <- phyloseq:::extract_eigenvalue.default(ordination)
                        }
                        
                        if (length(eigvec[axes]) > 0){
                                fracvar = eigvec[axes]/sum(eigvec)
                                percvar = round(100 * fracvar, 1)
                                strivar = as(c(Tr$label$x, Tr$label$y), "character")
                                strivar = paste0(strivar, " (", percvar, " %)")
                                Tr <- Tr + xlab(strivar[1]) + ylab(strivar[2]) 
                        }
                        
                        if (!is.null(eigvec) && coord_cor) {
                                Tr <- Tr + coord_fixed(sqrt(eigvec[2] / eigvec[1]))
                        }
                        
                }
                
                TrList_taxa[[i]] <- Tr
                rm(Tr)
                
        }
        
        names(ordination_list) <- names(TrList_taxa) <- names(DFList) <- names(DF_taxa_List) <- names(TrList_own) <- names(dist_list)
        out <- list(ordination_list = ordination_list, DFList = DFList, DF_taxa_List = DF_taxa_List, ordination_Tr_samples = TrList_own, ordination_Tr_taxa = TrList_taxa)
}
# --




