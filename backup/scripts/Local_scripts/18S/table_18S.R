
#load packages
library(dplyr)
library(ggplot2)
library(tidyverse)
library(phyloseq)

# Calculate no. of OTUs per phylum before rarefaction

## Import taxonomy table
tax <- read.delim("../../Tekstfiler/18S/18S_classified_all.tsv", sep="\t", row.names=NULL)

# Load table with corrected phylum names and marine/non-marine
tax_cur<-read.table("../../Tekstfiler/18S/18S_classified_phy_class_curated.tsv",sep="\t",header=T)

# Remove non-marine classes from taxonomy table
tax_mar <- tax[!tax$class %in% tax_cur$class[tax_cur$marine=="no"],]

# Load table with supergroups and uni-/multicellular
sgroups <- read.table("../../Tekstfiler/Across_barcodes/Supergroups_and_cellularity.tsv", sep='\t', header=T, comment="")

## Add curated phylum names
tax_mar$new_phylum<-tax_cur$new_phylum[match(tax_mar$class,tax_cur$class)]

## Add supergroup
tax_mar$supergroup<-sgroups$supergroup[match(tax_mar$new_phylum,sgroups$new_phylum)]

# Extract needed columns
tax_mar <- tax_mar[,c("supergroup","new_phylum","class")]

# Check that there are no "NAs"
unique(tax_mar$supergroup)
unique(tax_mar$new_phylum)

#Count no.of OTUs per phylum
n_otus<-tax_mar %>% group_by(new_phylum) %>% summarise(n=n())

#check total no. of OTUs
sum(n_otus$n)
#10673

#Sort phyla by no. of OTUs
n_otus <- n_otus[order(-(n_otus$n)),]

#rename columns
colnames(n_otus)[1]="Phylum"
colnames(n_otus)[2]="No. of MOTUs"

# Calculate no. of OTUs per phylum after rarefaction and removal of cluster 2 and samples w. only 2 reps

## Get tax table
euk_rare<-readRDS("../../RDS/18S_no_c2_3reps.rds")
taxa_rare<-as.data.frame(tax_table(euk_rare))

## Add curated phylum names
taxa_rare$new_phylum<-tax_cur$new_phylum[match(taxa_rare$class,tax_cur$class)]

## Add supergroup
taxa_rare$supergroup<-sgroups$supergroup[match(taxa_rare$new_phylum,sgroups$new_phylum)]

# Extract needed columns
taxa_rare <- taxa_rare[,c("supergroup","new_phylum","class")]

#Count no.of OTUs per phylum
n_otus_rare<-taxa_rare %>% group_by(new_phylum) %>% summarise(n=n())

#check total no. of OTUs
sum(n_otus_rare$n)
#10672

# Add no. of MOTUs after rarefaction to previous table
n_otus$"No. of MOTUs after rarefaction"<-n_otus_rare$n[match(n_otus$Phylum,n_otus_rare$new_phylum)]

# Calculate no. of clusters per phylum

## Get OTU table
otus1<-as.data.frame(otu_table(euk_rare))

##Get metadata
metadata<-as.data.frame(sample_data(euk_rare))

###Make phyloseq object
tax_mat_b<-as.matrix(taxa_rare)
OTU = otu_table(otus1, taxa_are_rows = FALSE)
TAX_b = tax_table(tax_mat_b)
samples = sample_data(metadata)
p_NCBI = phyloseq(OTU, TAX_b, samples)

## Aggregate OTUs from the same phylum
p_NCBI_agg<-tax_glom(p_NCBI,taxrank = "new_phylum")

## Merge samples from the same sampling cluster to allow estimating frequency of taxa (in terms of presence in X no. of clusters)

merged = merge_samples(p_NCBI_agg, "cluster")
SD = merge_samples(sample_data(p_NCBI_agg), "cluster")
#print(SD)
#print(merged)
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
n_otus$"No. of clusters" <- no_clusters_df$no_clusters[match(n_otus$Phylum,no_clusters_df$phylum)]

#Add supergroup column
n_otus$Supergroup <- tax$supergroup[match(n_otus$Phylum,tax$new_phylum)]

#Sort by supergroup, then by no. of OTUs
n_otus <- n_otus[order(n_otus$Supergroup,-(n_otus$"No. of MOTUs")),]

# Reorder columns by name manually
new_order = c("Supergroup","Phylum","No. of MOTUs","No. of MOTUs after rarefaction","No. of clusters")
n_otus <- n_otus[, new_order]

## Write to file
write.table(n_otus,file="../../Tekstfiler/18S/Table_18S.tsv",sep="\t")
