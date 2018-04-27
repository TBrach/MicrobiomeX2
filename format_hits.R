function (results, p.adjust.threshold = 0.1, p.adjust.method = NULL) 
{
        df <- results$table
        df <- dplyr::filter(df, !is.na(pvalue) & !is.na(padj))
        if (!is.null(p.adjust.method)) {
                df$padj = p.adjust(df$pvalue, method = p.adjust.method)
        }
        df <- dplyr::filter(df, padj <= p.adjust.threshold)
        df$Annotation <- apply(df, 1, get_pretty_taxon_name)
        specific_fields <- NULL
        if (!is.null(df$OR)) {
                df$OR_lb <- as.numeric(df$OR_lb)
                df <- dplyr::arrange(df, Direction, desc(OR_lb), padj, 
                                     OR, Taxon)
                df$OR_CE <- paste0("[", df$OR_lb, ", ", df$OR_ub, "]")
                specific_fields = c("OR", "OR_CE")
        }
        if (!is.null(df$log2FC)) {
                df <- dplyr::arrange(df, Direction, desc(abs(log2FC)), 
                                     padj, Taxon)
                specific_fields = c("log2FC")
        }
        if (!is.null(df$W)) {
                df <- dplyr::arrange(df, Direction, desc(abs(teststat)), 
                                     padj, Med_1, Taxon)
                specific_fields = c("Med_1", "Med_2")
        }
        if (!is.null(df$H)) {
                df <- dplyr::arrange(df, Direction, desc(abs(H)), padj, 
                                     Taxon)
                specific_fields = c("H")
        }
        df <- dplyr::select(df, Taxon, Annotation, padj, Significance, 
                            Direction, specific_fields)
        rownames(df) <- df$Taxon
        my_results <- results
        my_results$table <- df
        my_results$nhits <- dim(df)[1]
        return(my_results)
}