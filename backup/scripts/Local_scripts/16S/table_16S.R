
#load packages
library(dplyr)
library(ggplot2)
library(tidyverse)
library(phyloseq)

# Calculate no. of OTUs per phylum after rarefaction

## Get tax table
prok_rare<-readRDS("../../RDS/16S_no_c2_3reps.rds")
taxa_rare<-as.data.frame(tax_table(prok_rare))

#extract kingdom, phylum columns
head(taxa_rare)
taxa1_rare <- taxa_rare[,c(1,2)]

#Count no.of OTUs per phylum
n_otus_rare<-taxa1_rare %>% group_by(Phylum) %>% summarise(n=n())

#check total no. of OTUs
sum(n_otus_rare$n)
#54418

#rename column 2 
colnames(n_otus_rare)[2]="No. of MOTUs"

# Calculate no. of clusters per phylum

## Aggregate OTUs from the same phylum
p_NCBI_agg<-tax_glom(prok_rare,taxrank = "Phylum")

## Merge samples from the same sampling cluster to allow estimating frequency of taxa (in terms of presence in X no. of clusters)

merged = merge_samples(p_NCBI_agg, "cluster")
SD = merge_samples(sample_data(p_NCBI_agg), "cluster")
#print(SD)
#print(merged)
sample_names(p_NCBI_agg)
sample_names(merged)
identical(SD, sample_data(merged)) #TRUE

## Count for each taxon the number of clusters where it is present
#!=0 fjerner alle celler med værdien 0, altså tæller kun de clusters sammen, hvor værdien er over nul og den givne phylum altså er fundet
no_clusters<-colSums(merged@otu_table != 0)
no_clusters_df <- as.data.frame(no_clusters)

## Get new tax table
tax_merged<-data.frame(tax_table(merged),check.names=F)

# Merge all three counts in one table

## Add phylum name to list of cluster counts
no_clusters_df$Phylum <- tax_merged$Phylum[match(row.names(no_clusters_df),row.names(tax_merged))]

## Add no. of clusters to table of OTU counts
n_otus_rare$"No. of stations" <- no_clusters_df$no_clusters[match(n_otus_rare$Phylum,no_clusters_df$Phylum)]

#Add kingdom and phylum columns and remove combine kingdom+phylum column
n_otus_rare$Kingdom <- tax_merged$Kingdom[match(n_otus_rare$Phylum,tax_merged$Phylum)]

#Sort by kingdom, then by no. of OTUs
n_otus_rare <- n_otus_rare[order(n_otus_rare$Kingdom,-(n_otus_rare$"No. of MOTUs")),]

# Reorder columns by name manually
new_order = c("Kingdom","Phylum","No. of MOTUs","No. of stations")
n_otus_rare <- n_otus_rare[, new_order]

## Write to file
write.table(n_otus_rare,file="../../Tekstfiler/16S/Table_16S.tsv",sep="\t")
