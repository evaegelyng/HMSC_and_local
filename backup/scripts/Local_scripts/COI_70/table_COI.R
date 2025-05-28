
#load packages
library(dplyr)
library(ggplot2)
library(tidyverse)
library(phyloseq)

# Load table with supergroups
sgroups <- read.table("../../Tekstfiler/Across_barcodes/Supergroups_and_cellularity.tsv", sep='\t', header=T, comment="")

# Calculate no. of OTUs per phylum after rarefaction

## Get phyloseq object
COSQ_rare2_70 <- readRDS("../../RDS/COI_no_c2_3reps.rds")

## Get OTU table
otu_70<-data.frame(otu_table(COSQ_rare2_70),check.names=F)

## Get taxonomy
taxa_rare<-data.frame(tax_table(COSQ_rare2_70))

## Add supergroup
taxa_rare$supergroup<-sgroups$supergroup[match(taxa_rare$new_phylum,sgroups$new_phylum)]

## Add supergroup to Bacillariophyta_NA
taxa_rare$supergroup[taxa_rare$new_phylum=="Bacillariophyta_NA"]<-"SAR_Stramenopiles"

# Extract needed columns
taxa_rare <- taxa_rare[,c("supergroup","new_phylum","new_class")]

#Count no.of OTUs per phylum
n_otus_rare<-taxa_rare %>% group_by(new_phylum) %>% summarise(n=n())

#check total no. of OTUs
sum(n_otus_rare$n)
#6356

#rename columns
colnames(n_otus_rare)[1]="Phylum"
colnames(n_otus_rare)[2]="No. of MOTUs"

# Calculate no. of clusters per phylum

##Get metadata
metadata<-as.data.frame(sample_data(COSQ_rare2_70))

###Make phyloseq object
tax_mat_b<-as.matrix(taxa_rare)
OTU = otu_table(otu_70, taxa_are_rows = FALSE)
TAX_b = tax_table(tax_mat_b)
samples = sample_data(metadata)
p_NCBI = phyloseq(OTU, TAX_b, samples)

## Aggregate OTUs from the same phylum
p_NCBI_agg<-tax_glom(p_NCBI,taxrank = "new_phylum")

## Merge samples from the same sampling cluster to allow estimating frequency of taxa (in terms of presence in X no. of clusters)

merged = merge_samples(p_NCBI_agg, "cluster")
SD = merge_samples(sample_data(p_NCBI_agg), "cluster")

sample_names(p_NCBI_agg)
sample_names(merged)
identical(SD, sample_data(merged)) #TRUE

## Count for each taxon the number of clusters where it is present
no_clusters<-colSums(merged@otu_table != 0)
no_clusters_df <- as.data.frame(no_clusters)

## Get new tax table
tax_merged<-data.frame(tax_table(merged),check.names=F)

# Merge all three counts in one table

## Add phylum name to list of cluster counts
no_clusters_df$phylum <- tax_merged$new_phylum[match(row.names(no_clusters_df),row.names(tax_merged))]

## Add no. of clusters to table of OTU counts
n_otus_rare$"No. of clusters" <- no_clusters_df$no_clusters[match(n_otus_rare$Phylum,no_clusters_df$phylum)]

#Add supergroup column
n_otus_rare$Supergroup <- taxa_rare$supergroup[match(n_otus_rare$Phylum,taxa_rare$new_phylum)]

#Sort by supergroup, then by no. of OTUs
n_otus_rare <- n_otus_rare[order(n_otus_rare$Supergroup,-(n_otus_rare$"No. of MOTUs")),]

# Reorder columns by name manually
new_order = c("Supergroup","Phylum","No. of MOTUs","No. of clusters")
n_otus_rare <- n_otus_rare[, new_order]

## Write to file
write.table(n_otus_rare,file="../../Tekstfiler/COI/COI_70/Table_COI.tsv",sep="\t")
