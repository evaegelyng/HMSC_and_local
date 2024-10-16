
#load packages
library(dplyr)
library(ggplot2)
library(tidyverse)
library(phyloseq)

# Calculate no. of OTUs per phylum before rarefaction

## Get tax table
taxa<-read.table("../Tekstfiler/18S/cleaned_tax_pident90_lulu_97.txt", sep="\t", row.names = 1, header=T)

#extract kingdom, phylum columns
head(taxa)
taxa1 <- taxa[,c(1,2)]
#Add combined kingdom and phylum column to account for "NA" in phylum column
taxa1$kin_phy<-paste(taxa1$kingdom,taxa1$phylum,sep="_")

#Count no.of OTUs per phylum
n_otus<-taxa1 %>% group_by(kin_phy) %>% summarise(n=n())

#check total no. of OTUs
sum(n_otus$n)
#10884

#Sort phyla by no. of OTUs
n_otus <- n_otus[order(-(n_otus$n)),]

#rename column 2 
colnames(n_otus)[2]="No. of MOTUs"

# Calculate no. of OTUs per phylum after rarefaction

## Get tax table
taxa_rare<-read.table("../Tekstfiler/18S/tax_rarefy_97.txt", sep="\t", row.names = 1, header=T)

#extract kingdom, phylum columns
head(taxa_rare)
taxa1_rare <- taxa_rare[,c(1,2)]
#Add combined kingdom and phylum column to account for "NA" in phylum column
taxa1_rare$kin_phy<-paste(taxa1_rare$kingdom,taxa1_rare$phylum,sep="_")

#Count no.of OTUs per phylum
n_otus_rare<-taxa1_rare %>% group_by(kin_phy) %>% summarise(n=n())

#check total no. of OTUs
sum(n_otus_rare$n)
#10883

# Add no. of MOTUs after rarefaction to previous table
n_otus$"No. of MOTUs after rarefaction"<-n_otus_rare$n[match(n_otus$kin_phy,n_otus_rare$kin_phy)]

# Calculate no. of clusters per phylum

## Get OTU table
otus1<-read.table("../Tekstfiler/18S/COSQ_Curated_LULU_97.tsv", sep="\t", header=T, row.names=1,check.names=F)
otus1<-within(otus1,rm(total_reads))

##Get metadata
metadata<-read.table("../Tekstfiler/18S/cleaned_metadata_pident90.txt", sep="\t", header=T)

###Make phyloseq object
tax_mat_b<-as.matrix(taxa1)
OTU = otu_table(otus1, taxa_are_rows = TRUE)
TAX_b = tax_table(tax_mat_b)
samples = sample_data(metadata)
p_NCBI = phyloseq(OTU, TAX_b, samples)

## Aggregate OTUs from the same phylum
p_NCBI_agg<-tax_glom(p_NCBI,taxrank = "kin_phy")

## Merge samples from the same sampling cluster to allow estimating frequency of taxa (in terms of presence in X no. of clusters)

merged = merge_samples(p_NCBI_agg, "cluster")
SD = merge_samples(sample_data(p_NCBI_agg), "cluster")
#print(SD)
#print(merged)
sample_names(p_NCBI_agg)
sample_names(merged)
identical(SD, sample_data(merged)) #TRUE

## Count for each taxon the number of clusters where it is present
## Optionally, subset to metazoans first
## merged<-subset_taxa(merged,kingdom=="Metazoa")
#!=0 fjerner alle celler med værdien 0, altså tæller kun de clusters sammen, hvor værdien er over nul og den givne phylum altså er fundet
no_clusters<-colSums(merged@otu_table != 0)
no_clusters_df <- as.data.frame(no_clusters)

## Get new tax table
tax_merged<-data.frame(tax_table(merged),check.names=F)


# Merge all three counts in one table

## Add phylum name to list of cluster counts
no_clusters_df$kin_phy <- tax_merged$kin_phy[match(row.names(no_clusters_df),row.names(tax_merged))]

## Add no. of clusters to table of OTU counts
n_otus$"No. of clusters" <- no_clusters_df$no_clusters[match(n_otus$kin_phy,no_clusters_df$kin_phy)]

#Add kingdom and phylum columns and remove combine kingdom+phylum column
n_otus$Kingdom <- taxa1$kingdom[match(n_otus$kin_phy,taxa1$kin_phy)]
n_otus$Phylum <- taxa1$phylum[match(n_otus$kin_phy,taxa1$kin_phy)]
n_otus<-within(n_otus,rm(kin_phy))

#Sort by kingdom, then by no. of OTUs
n_otus <- n_otus[order(n_otus$Kingdom,-(n_otus$"No. of MOTUs")),]

# Reorder columns by name manually
new_order = c("Kingdom","Phylum","No. of MOTUs","No. of MOTUs after rarefaction","No. of clusters")
n_otus <- n_otus[, new_order]

## Write to file
write.table(n_otus,file="../Tekstfiler/18S/Table_18S.tsv",sep="\t")
