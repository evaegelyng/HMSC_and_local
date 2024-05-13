library(phyloseq)
library(metacoder)
library(dplyr)
library(ggplot2)
library(ggpubr)

## Get phyloseq object
COSQ_rare2_70 <- readRDS("../RDS/COI_no_c2_3reps.rds")

## Get OTU table
otu_70<-data.frame(otu_table(COSQ_rare2_70),check.names=F)

## Get taxonomy
tax_70<-data.frame(tax_table(COSQ_rare2_70))

## Get sample data
sam<-data.frame(sample_data(COSQ_rare2_70))

## Add substrate to OTU table
otu_70$Substrate<-sam$substrate_type[match(row.names(otu_70),row.names(sam))]

#### Sediment ####
#Subset data to sediment samples
sed<-otu_70[otu_70$Substrate=="sediment",]
sed<-within(sed,rm(Substrate))

## Remove empty samples
sed<-sed[rowSums(sed[])>0,] 

# Remove empty colums (otus)
sed<-sed[colSums(sed)>0]

#transpose the sediment df
sed_t <- as.data.frame(t(sed))

#create columns of rownames
sed_t$otus <- rownames(sed_t)
tax_70$otus <- rownames(tax_70)

#Add domain column to "artificially" fix root of plot
sed_t$domain <- "Eukaryota"

#add lineage from tax table
sed_tax <- sed_t %>% 
  left_join(tax_70, by = c("otus" = "otus")) %>% 
  relocate(otus, .after = species)

#remove NAs (note that they are not seen as NAs by R)
sed_full_tax <- sed_tax %>% 
  filter(order != "NA") %>% 
  filter(class != "NA") 

sed_full_tax[,1:155] <- as.numeric(unlist(sed_full_tax[,1:155]))

#-------------------------------------------------------------------------------
#              HEAT TREE COMPARISON  relative read abundance                 
#-------------------------------------------------------------------------------

#Create df of sample ids, habitat and season
sed_samples <- data.frame(sample_id = colnames(sed_full_tax[,1:155]))
sed_samples$Habitat <- sam$habitat[match(sed_samples$sample_id,row.names(sam))]
sed_samples$Season <- sam$season[match(sed_samples$sample_id,row.names(sam))]

#vector of order of habitats
target <- c("rocks", "eelgrass", "sand")

#arrange samples so rocks comes first, important for placement of comparison trees
sed_samples <- sed_samples %>% arrange(factor(Habitat, levels = target))


#### Autumn season ####
#filter for season 
autumn_id <- sed_samples %>% filter(Season == "autumn") %>% pull(sample_id)
autumn_habitat <- sed_samples %>% filter(Season == "autumn") %>% pull(Habitat)
sed_sam_tax <- sed_full_tax %>% select(all_of(autumn_id)) %>% bind_cols(sed_full_tax[,156:160])


#create meatacoder object with classification until order level
obj_sed <- parse_tax_data(sed_sam_tax,
                          class_cols = c("domain", "kingdom", "phylum","class","order"))

#calculate taxonomic abundance
obj_sed$data$tax_table <- calc_taxon_abund(obj_sed, data = "tax_data", cols = autumn_id)

#calculate differences between habitats
obj_sed$data$diff_table <- compare_groups(obj_sed, data = "tax_table",
                                          cols = autumn_id, # What columns of sample data to use
                                          groups = autumn_habitat) # What category each sample is assigned to
print(obj_sed$data$diff_table)

set.seed(1)
sed_aut<-heat_tree_matrix(obj_sed,
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


#### Spring season ####
#filter for season 
spring_id <- sed_samples %>% filter(Season == "spring") %>% pull(sample_id)
spring_habitat <- sed_samples %>% filter(Season == "spring") %>% pull(Habitat)
sed_sam_tax <- sed_full_tax %>% select(all_of(spring_id)) %>% bind_cols(sed_full_tax[,156:160])


#create meatacoder object with classification until order level
obj_sed <- parse_tax_data(sed_sam_tax,
                          class_cols = c("domain", "kingdom", "phylum","class","order"))

#calculate taxonomic abundance
obj_sed$data$tax_table <- calc_taxon_abund(obj_sed, data = "tax_data", cols = spring_id)

#calculate differences between habitats
obj_sed$data$diff_table <- compare_groups(obj_sed, data = "tax_table",
                                          cols = spring_id, # What columns of sample data to use
                                          groups = spring_habitat) # What category each sample is assigned to
print(obj_sed$data$diff_table)

set.seed(1)
sed_spr<-heat_tree_matrix(obj_sed,
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


#### Water ####
#Subset data to water samples
wat<-otu_70[otu_70$Substrate=="water",]
wat<-within(wat,rm(Substrate))

## Remove empty samples
wat<-wat[rowSums(wat[])>0,] 

# Remove empty colums (otus)
wat<-wat[colSums(wat)>0]

#transpose the water df
wat_t <- as.data.frame(t(wat))

#create columns of rownames
wat_t$otus <- rownames(wat_t)
tax_70$otus <- rownames(tax_70)

#Add domain column to "artificially" fix root of plot
wat_t$domain <- "Eukaryota"

#add lineage from tax table
wat_tax <- wat_t %>% 
  left_join(tax_70, by = c("otus" = "otus")) %>% 
  relocate(otus, .after = species)

#remove NAs (note that they are not seen as NAs by R)
wat_full_tax <- wat_tax %>% 
  filter(order != "NA") %>% 
  filter(class != "NA") 

wat_full_tax[,1:151] <- as.numeric(unlist(wat_full_tax[,1:151]))

#-------------------------------------------------------------------------------
#              HEAT TREE COMPARISON  relative read abundance                 
#-------------------------------------------------------------------------------

#Create df of sample ids, habitat and season
wat_samples <- data.frame(sample_id = colnames(wat_full_tax[,1:151]))
wat_samples$Habitat <- sam$habitat[match(wat_samples$sample_id,row.names(sam))]
wat_samples$Season <- sam$season[match(wat_samples$sample_id,row.names(sam))]

#vector of order of habitats
target <- c("rocks", "eelgrass", "sand")

#arrange samples so rocks comes first, important for placement of comparison trees
wat_samples <- wat_samples %>% arrange(factor(Habitat, levels = target))


#### Autumn season ####
#filter for season 
autumn_id <- wat_samples %>% filter(Season == "autumn") %>% pull(sample_id)
autumn_habitat <- wat_samples %>% filter(Season == "autumn") %>% pull(Habitat)
wat_sam_tax <- wat_full_tax %>% select(all_of(autumn_id)) %>% bind_cols(wat_full_tax[,152:156])


#create meatacoder object with classification until order level
obj_wat <- parse_tax_data(wat_sam_tax,
                          class_cols = c("domain", "kingdom", "phylum","class","order"))

#calculate taxonomic abundance
obj_wat$data$tax_table <- calc_taxon_abund(obj_wat, data = "tax_data", cols = autumn_id)

#calculate differences between habitats
obj_wat$data$diff_table <- compare_groups(obj_wat, data = "tax_table",
                                          cols = autumn_id, # What columns of sample data to use
                                          groups = autumn_habitat) # What category each sample is assigned to
print(obj_wat$data$diff_table)

set.seed(1)
wat_aut<-heat_tree_matrix(obj_wat,
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


#### Spring season ####
#filter for season 
spring_id <- wat_samples %>% filter(Season == "spring") %>% pull(sample_id)
spring_habitat <- wat_samples %>% filter(Season == "spring") %>% pull(Habitat)
wat_sam_tax <- wat_full_tax %>% select(all_of(spring_id)) %>% bind_cols(wat_full_tax[,152:156])


#create meatacoder object with classification until order level
obj_wat <- parse_tax_data(wat_sam_tax,
                          class_cols = c("domain", "kingdom", "phylum", "class","order"))

#calculate taxonomic abundance
obj_wat$data$tax_table <- calc_taxon_abund(obj_wat, data = "tax_data", cols = spring_id)

#calculate differences between habitats
obj_wat$data$diff_table <- compare_groups(obj_wat, data = "tax_table",
                                          cols = spring_id, # What columns of sample data to use
                                          groups = spring_habitat) # What category each sample is assigned to
print(obj_wat$data$diff_table)

set.seed(1)
wat_spr<-heat_tree_matrix(obj_wat,
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


# Combine all four plots

plots<-list()

plots[[1]]<-sed_spr
plots[[2]]<-sed_aut
plots[[3]]<-wat_spr
plots[[4]]<-wat_aut

p<-ggarrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],ncol=2,nrow=2,align="h")
ggsave(file="../Plots/Heattrees_COI.pdf",annotate_figure(p,top="Spring                                                                                 Autumn",
                                                         left="Water                                                                Sediment"))
