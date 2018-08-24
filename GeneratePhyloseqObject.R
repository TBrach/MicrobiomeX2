library(lubridate)


# - Load all dada data in -
datapath <- "/Users/jvb740/MarieCurie_Work/NormalizationProjectNastya/Results/16S_Sequencing_Pre/Dada_Analysis"
load(file.path(datapath, "Dada_Data/DenoisedData.RData")) # laod seqtab.nochim
load(file.path(datapath, "Dada_Taxonomy/Silva_v128/Taxonomy.RData")) # load taxa.species
tree_list <- readRDS(file.path(datapath, "Dada_phylogenetic_tree/phylog_tree.rds")) # load tree
# --


# prepare sample_data -
samdf <- read.csv2(file = file.path(datapath, "PhyseqAnalysis/SampleInfoSheetForBGI.csv"), header = TRUE,
                   stringsAsFactors = FALSE)


samdf$Date <- lubridate::parse_date_time(samdf$Date, orders = "dmy", tz = "CET")
rownames(samdf) <- rownames(seqtab.nochim) # check that it fits first!
# --



# - generate phyloseq object --
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), 
               sample_data(samdf), 
               tax_table(taxa.species),
               phy_tree(tree_list[["fitGTR"]]$tree))

# --

# - add numbers of filtered reads -
# the number of Filtered reads are added to later check the relation of alpha diversity measures to the number of filtered reads
FilteredReads <- ReadSummary[,c("Sample","FilteredReads")]
# usually should be in same order but better check:
FilteredReads <- FilteredReads[match(sample_names(ps), FilteredReads$Sample),]
sample_data(ps)$FilteredReads <- FilteredReads$FilteredReads
# --

# - save the phyloseq object -
saveRDS(object = ps, file = file.path(datapath, "PhyseqAnalysis/phyloseq_objects/ps_normalization.rds"))


# --
