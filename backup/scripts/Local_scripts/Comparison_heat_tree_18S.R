library(phyloseq)
library(metacoder)
library(dplyr)
library(ggplot2)
library(ggpubr)

###Import phyloseq object
DADAwang1<-readRDS("../RDS/18S_no_c2_3reps.rds")

## Get cleaned OTU table
otus<-data.frame(otu_table(DADAwang1),check.names=F)

## Get taxonomy
tax<-data.frame(tax_table(DADAwang1))

## Get sample data
sam<-data.frame(sample_data(DADAwang1))

otus$Substrate<-sam$substrate_type[match(row.names(otus),row.names(sam))]

#Subset to sediment samples
sed<-otus[otus$Substrate=="sediment",]
sed<-within(sed,rm(Substrate))

## Remove empty samples
sed<-sed[rowSums(sed[])>0,] 

# Remove empty colums (otus)
sed<-sed[colSums(sed)>0]

#transpose the sediment df
sed_t <- as.data.frame(t(sed))

#create columns of rownames
sed_t$otus <- rownames(sed_t)
tax$otus <- rownames(tax)

#Add domain column to "artificially" fix root of plot
sed_t$domain <- "Eukaryota"

#add lineage from tax table
sed_tax <- sed_t %>% 
  left_join(tax, by = c("otus" = "otus")) %>% 
  relocate(otus, .after = Genus)

sed_full_tax <- sed_tax
sed_full_tax[,1:149] <- as.numeric(unlist(sed_full_tax[,1:149]))

#-------------------------------------------------------------------------------
#              HEAT TREE COMPARISON  relative read abundance                 
#-------------------------------------------------------------------------------

#Create df of sample ids, habitat and season
sed_samples <- data.frame(sample_id = colnames(sed_full_tax[,1:149]))
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
sed_sam_tax <- sed_full_tax %>% select(all_of(autumn_id)) %>% bind_cols(sed_full_tax[,150:158])


#create meatacoder object with classification until class level
obj_sed <- parse_tax_data(sed_sam_tax,
                          class_cols = c("domain", "Division", "Phylum","Class"))

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
sed_sam_tax <- sed_full_tax %>% select(all_of(spring_id)) %>% bind_cols(sed_full_tax[,150:158])


#create meatacoder object with classification until class level
obj_sed <- parse_tax_data(sed_sam_tax,
                          class_cols = c("domain", "Division", "Phylum","Class"))

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
wat<-otus[otus$Substrate=="water",]
wat<-within(wat,rm(Substrate))

## Remove empty samples
wat<-wat[rowSums(wat[])>0,] 

# Remove empty colums (otus)
wat<-wat[colSums(wat)>0]

#transpose the water df
wat_t <- as.data.frame(t(wat))

#create columns of rownames
wat_t$otus <- rownames(wat_t)
tax$otus <- rownames(tax)

#Add domain column to "artificially" fix root of plot
wat_t$domain <- "Eukaryota"

#add lineage from tax table
wat_tax <- wat_t %>% 
  left_join(tax, by = c("otus" = "otus")) %>% 
  relocate(otus, .after = Genus)

wat_full_tax <- wat_tax
wat_full_tax[,1:149] <- as.numeric(unlist(wat_full_tax[,1:149]))

#-------------------------------------------------------------------------------
#              HEAT TREE COMPARISON  relative read abundance                 
#-------------------------------------------------------------------------------

#Create df of sample ids, habitat and season
wat_samples <- data.frame(sample_id = colnames(wat_full_tax[,1:149]))
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
wat_sam_tax <- wat_full_tax %>% select(all_of(autumn_id)) %>% bind_cols(wat_full_tax[,150:158])


#create meatacoder object with classification until class level
obj_wat <- parse_tax_data(wat_sam_tax,
                          class_cols = c("domain", "Division", "Phylum","Class"))

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
wat_sam_tax <- wat_full_tax %>% select(all_of(spring_id)) %>% bind_cols(wat_full_tax[,150:158])


#create meatacoder object with classification until order level
obj_wat <- parse_tax_data(wat_sam_tax,
                          class_cols = c("domain", "Division", "Phylum", "Class"))

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
ggsave(file="../Plots/Heattrees_18S.pdf",annotate_figure(p,top="Spring                                                                                 Autumn",
                                                         left="Water                                                                Sediment"))

