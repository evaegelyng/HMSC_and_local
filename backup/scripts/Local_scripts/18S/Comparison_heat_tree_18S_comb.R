library(phyloseq)
library(metacoder)
library(dplyr)
library(ggplot2)
library(ggpubr)

###Import phyloseq object
DADAwang1<-readRDS("../../RDS/18S_no_c2_3reps.rds")

## Get cleaned OTU table
otus<-data.frame(otu_table(DADAwang1),check.names=F)

## Get taxonomy
tax<-data.frame(tax_table(DADAwang1))

# Load table with corrected phylum names and marine/non-marine
tax_cur<-read.table("../../Tekstfiler/18S/18S_classified_phy_class_curated.tsv",sep="\t",header=T)

# Load table with supergroups and uni-/multicellular
sgroups <- read.table("../../Tekstfiler/Across_barcodes/Supergroups_and_cellularity.tsv", sep='\t', header=T, comment="")

## Add curated phylum names
tax$new_phylum<-tax_cur$new_phylum[match(tax$class,tax_cur$class)]

## Add supergroup
tax$supergroup<-sgroups$supergroup[match(tax$new_phylum,sgroups$new_phylum)]

# Extract needed columns
tax <- tax[,c("supergroup","new_phylum","class")]

# Check that there are no "NAs"
unique(tax$supergroup)

## Get sample data
sam<-data.frame(sample_data(DADAwang1))

#transpose the df
otus_t <- as.data.frame(t(otus))

#create columns of rownames
otus_t$otus <- rownames(otus_t)
tax$otus <- rownames(tax)

#Add domain column to "artificially" fix root of plot
otus_t$domain <- "Eukaryota"

#add lineage from tax table
otus_tax <- otus_t %>% 
  left_join(tax, by = c("otus" = "otus")) %>% 
  relocate(otus, .after = class)

otus_full_tax <- otus_tax
otus_full_tax[,1:310] <- as.numeric(unlist(otus_full_tax[,1:310]))

#-------------------------------------------------------------------------------
#              HEAT TREE COMPARISON  relative read abundance                 
#-------------------------------------------------------------------------------

#Create df of sample ids, habitat and season
samples <- data.frame(sample_id = colnames(otus_full_tax[,1:310]))
samples$Habitat <- sam$habitat[match(samples$sample_id,row.names(sam))]
samples$Season <- sam$season[match(samples$sample_id,row.names(sam))]

#vector of order of habitats
target <- c("rocks", "eelgrass", "sand")

#arrange samples so rocks comes first, important for placement of comparison trees
samples <- samples %>% arrange(factor(Habitat, levels = target))


#### Autumn season ####
#filter for season 
id <- samples %>% pull(sample_id)
habitat <- samples %>% pull(Habitat)

#create meatacoder object with classification until class level
obj <- parse_tax_data(otus_full_tax,
                          class_cols = c("domain", "supergroup", "new_phylum","class"))

#calculate taxonomic abundance
obj$data$tax_table <- calc_taxon_abund(obj, data = "tax_data", cols=id)

#calculate differences between habitats
obj$data$diff_table <- compare_groups(obj, data = "tax_table",
                                          cols = id, # What columns of sample data to use
                                          groups = habitat) # What category each sample is assigned to
print(obj$data$diff_table)

set.seed(1)
hab_tree<-heat_tree_matrix(obj,
                          data = "diff_table",
                          label_small_trees = TRUE,
                          node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                          node_label = taxon_names,
                          node_color = log2_median_ratio, # A column from obj$data$diff_table
                          node_color_range = diverging_palette(), # The built-in palette for diverging data
                          node_color_trans = "linear", # The default is scaled by circle area
                          node_color_interval = c(-3, 3), # The range of mean_diff to display
                          edge_color_interval = c(-3, 3), # The range of mean_diff to display
                          node_size_axis_label = "OTU count",
                          node_color_axis_label = "Log2 median ratio",
                          layout = "davidson-harel", # The primary layout algorithm
                          initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations

hab_tree

ggsave(hab_tree,file="../../Plots/Heattrees/Heattrees_18S_comb.pdf")                                                     
