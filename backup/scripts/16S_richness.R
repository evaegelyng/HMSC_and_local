#Load libraries
library("phyloseq")
library("dplyr")

###Make phyloseq object from raw data
otu_mat<-as.matrix(read.table("data/f_otu_silva_5_01_22.txt", sep="\t", header=T, row.names=1,check.names=F))

taxonomy_silva<-read.table("data/f_tax_silva_5_01_22.txt", sep='\t', header=T, comment="")
tax_mat_s<-as.matrix(taxonomy_silva)

OTUSILVA = otu_table(otu_mat, taxa_are_rows = FALSE)
TAX_S = tax_table(tax_mat_s)
p_SILVA = phyloseq(OTUSILVA, TAX_S)

#Load metadata
metadata<-read.table("data/f_silva_metadata_5_01_22.txt", sep="\t", header=T)

#Create extra variables
metadata$sshc<-paste(metadata$substrate_type, metadata$season, metadata$habitat, metadata$cluster, sep="_")
##Including the number of field reps per sample
bio_reps<-metadata %>% group_by(sshc) %>% summarise(n=n())
metadata$bio_reps<-bio_reps$n[match(metadata$sshc,bio_reps$sshc)]

sampledata = sample_data(data.frame(metadata, row.names=metadata$root, stringsAsFactors=FALSE))
DADAwang1 = merge_phyloseq(p_SILVA, sampledata)
DADAwang1

##Remove cluster 2, which was only sampled in one season
DADAwang1<-subset_samples(DADAwang1, !cluster==2)
DADAwang1 # 6 field reps removed

#Select only samples with 3 replicates
DADAwang1<-subset_samples(DADAwang1, bio_reps==3)
DADAwang1 # 4 field reps (2 samples) removed

##Merge field replicates
DT1<-merge_samples(DADAwang1, "sshc", fun = mean)

## Extract OTU table
OTU_agg = data.frame(otu_table(DT1),check.names=F)

##Convert to presence/absence
OTU_agg[OTU_agg>0]<-1

## Calculate OTU richness per sample
rich<-rowSums(OTU_agg)
rich<-as.data.frame(rich)

write.table(rich,file="results/bact_rich.tsv",sep="\t",row.names=T)

#Sample data needs to be remade after merging of field replicates
d<-data.frame(sample_data(DT1)[,c("cluster","season","habitat")])
#Splits rowname into habitat from third argument
d$habitat<- sapply(strsplit(as.character(rownames(d)), "_"), function(x) unlist(strsplit(x[3], "_"))[1])
#Splits rowname into season from second argument
d$season<- sapply(strsplit(as.character(rownames(d)), "_"), function(x) unlist(strsplit(x[2], "_"))[1])
d$sshc<-rownames(d)
d$sch<-paste(d$season,d$cluster,d$habitat,sep="_")
d$substrate_type<-metadata$substrate_type[match(rownames(d),metadata$sshc)]

sample_data(DT1)<-d[,c("sshc","sch","cluster","season","habitat","substrate_type")]

#Export phyloseq object as rds file
saveRDS(DT1,file="results/16S_no_c2_3reps.rds")

