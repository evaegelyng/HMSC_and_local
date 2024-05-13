
#load packages
library(dplyr)
library(ggplot2)
library(tidyverse)
library(phyloseq)

# Calculate no. of OTUs per phylum before rarefaction

## Get tax table
taxa<-read.table("../Tekstfiler/cleaned_ncbi_12_12_22.txt")

# Add kingdoms
taxa<-taxa %>%
  mutate(Kingdom = case_when(Division == "Viridiplantae" ~ "Plantae",
                             Division %in% c("Alveolata", "Rhizaria", "Stramenopiles", "Cryptista", "Haptista") ~ 'Chromista',
                             Division %in% c("Rhodophyta","Choanoflagellata","Amoebozoa_IS","Apusomonadida","Ichthyosporea","Ancyromonadida","Discoba","Telonemida","Filasterea","Rotosphaerida","Breviatea","Mantamonadida","Rhodelphea","Metamonada","Tunicaraptor") ~ 'Protozoa',
                             Division ==  "Metazoa" ~ 'Metazoa',
                             Division ==  "Fungi" ~ 'Fungi',TRUE ~ NA_character_),
         Kingdom = factor(Kingdom, levels = c('Plantae',"Chromista","Protozoa","Metazoa","Fungi")))

#extract id, kingdom, phylum columns
taxa$id<-row.names(taxa)
taxa1 <- taxa[,c(3,8,9)]

#aggregate
phyla_OTU1 <- aggregate(.~Kingdom+Phylum,data=taxa1 , FUN=function(x) length(unique(x)))

#check no. of OTUs
sum(phyla_OTU1$id)
#4703

#Sort phyla by no. of OTUs
phyla_OTU1 <- phyla_OTU1[order(-(phyla_OTU1$id)),]

#rename column 3 and make phylum column into row names
colnames(phyla_OTU1)[3]="No. of MOTUs"
row.names(phyla_OTU1)<-phyla_OTU1$Phylum


# Calculate no. of OTUs per phylum after rarefaction

## Get tax table
taxa<-read.table("../Tekstfiler/f_tax_ncbi_5_01_22.txt")

# Add kingdoms
taxa<-taxa %>%
  mutate(Kingdom = case_when(Division %in% c("Viridiplantae", "Rhodophyta","Rhodelphea") ~ "Plantae",
                             Division %in% c("Alveolata", "Rhizaria", "Stramenopiles", "Cryptista", "Haptista", "Telonemida") ~ 'Chromista',
                             Division %in% c("Choanoflagellata","Amoebozoa_IS","Apusomonadida","Ichthyosporea","Ancyromonadida","Discoba","Filasterea","Rotosphaerida","Breviatea","Mantamonadida","Metamonada","Tunicaraptor") ~ 'Protozoa',
                             Division ==  "Metazoa" ~ 'Metazoa',
                             Division ==  "Fungi" ~ 'Fungi',TRUE ~ NA_character_),
         Kingdom = factor(Kingdom, levels = c('Plantae',"Chromista","Protozoa","Metazoa","Fungi")))

#extract id, kingdom, phylum columns
taxa$id<-row.names(taxa)
taxa1 <- taxa[,c(3,8,9)]

#aggregate
phyla_OTU <- aggregate(.~Kingdom+Phylum,data=taxa1 , FUN=function(x) length(unique(x)))

#check no. of OTUs
sum(phyla_OTU$id)
#4607

# Add no. of MOTUs after rarefaction to previous table
phyla_OTU1$"No. of MOTUs after rarefaction"<-phyla_OTU$id[match(row.names(phyla_OTU1),phyla_OTU$Phylum)]

# Calculate no. of clusters per phylum

## Get OTU table
otus1<-read.table("../Tekstfiler/cleaned_tax_ncbi_12_12_22.txt", sep="\t", header=T, row.names=1,check.names=F)

## Sum read counts per sample
sums <- colSums(otus1)
## Check range and mean/median values
summary(sums)
## Check std
sd(sums)

##Get metadata
metadata<-read.table("../Tekstfiler/cleaned_ncbi_metadata_12_12_22.txt", sep="\t", header=T)

###Make phyloseq object
tax_mat_b<-as.matrix(taxa)
OTU = otu_table(otus1, taxa_are_rows = TRUE)
TAX_b = tax_table(tax_mat_b)
samples = sample_data(metadata)
p_NCBI = phyloseq(OTU, TAX_b, samples)

## Aggregate OTUs from the same phylum
p_NCBI_agg<-tax_glom(p_NCBI,taxrank = "Phylum")

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
no_clusters_df$phylum <- tax_merged$Phylum[match(row.names(no_clusters_df),row.names(tax_merged))]

## Add no. of clusters to table of OTU counts
phyla_OTU1$"No. of clusters" <- no_clusters_df$no_clusters[match(row.names(phyla_OTU1),no_clusters_df$phylum)]

## Write to file
write.table(phyla_OTU1,file="../Tekstfiler/TableS8.tsv",sep="\t")

