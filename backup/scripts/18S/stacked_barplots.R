library("phyloseq")
library("ggplot2")
library("vegan")
library("reshape2")
library("plyr")
library("scales")
library("stringr")
library("RColorBrewer")
sessionInfo()

#Load tables

###Make phyloseq object from raw data
otu_mat<-as.matrix(read.table("data/f_otu_ncbi_5_01_22.txt", sep="\t", header=T, row.names=1,check.names=F))

taxonomy_ncbi<-read.table("data/f_tax_ncbi_5_01_22.txt", sep='\t', header=T, comment="")
tax_mat_b<-as.matrix(taxonomy_ncbi)

OTU = otu_table(otu_mat, taxa_are_rows = FALSE)
TAX_b = tax_table(tax_mat_b)
p_SILVA = phyloseq(OTU, TAX_b)

#Load metadata
metadata<-read.table("data/f_ncbi_metadata_5_01_22.txt", sep="\t", header=T)

#Create extra variables
metadata$sshc<-paste(metadata$substrate_type, metadata$season, metadata$habitat, metadata$cluster, sep="_")
metadata$ssc<-paste(metadata$substrate_type, metadata$season, metadata$cluster, sep="_")
metadata$stc<-paste(metadata$substrate_type, metadata$cluster, sep="_")

sampledata = sample_data(data.frame(metadata, row.names=metadata$sample_root, stringsAsFactors=FALSE))

DADAwang1 = merge_phyloseq(p_SILVA, sampledata)
DADAwang1

# Remove cluster 2, which was only sampled in spring
COSQ_no_c2<-subset_samples(DADAwang1, !cluster==2)
COSQ_no_c2 # Six samples were removed

#####
#####
#Taxonomy plots

#For stacked plots, use "tudao" as starting object:

tudao<-merge_samples(COSQ_no_c2, "sshc")
tudao = filter_taxa(tudao, function(x) sum(x) > 0, TRUE)
tudao = prune_samples(sample_sums(tudao)>0,tudao)

d<-data.frame(sample_data(tudao)[,c("cluster","season","habitat","substrate_type")])
#d$habitat<-cut(d$habitat, 3, labels=c("eelgrass","rocks","sand"))
d$habitat<-sampledatas$habitat[match(row.names(d),sampledatas$sshc)]
#d$substrate_type<-cut(d$substrate_type, 2, labels=c("sediment","water"))
d$substrate_type<-sampledatas$substrate_type[match(row.names(d),sampledatas$sshc)]
#d$season<-cut(d$season, 2, labels=c("autumn","spring"))
d$season<-sampledatas$season[match(row.names(d),sampledatas$sshc)]
d$sshc<-rownames(d)
head(d)
sample_data(tudao)<-d[,c("cluster","season","habitat","substrate_type","sshc")]

#Non-metazoans
no<-subset_taxa(tudao, !Division=="Metazoa")
datarg = transform_sample_counts(no, function(x) x/sum(x))

datag = tax_glom(datarg, "Phylum")
TopNOTUs = names(sort(taxa_sums(datag), TRUE)[1:6])
data9 = prune_taxa(TopNOTUs, datag)

d<-data.frame(sample_data(data9))

tax<-data.frame(tax_table(data9), stringsAsFactors=FALSE)
otu<-otu_table(data9)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Phylum
tab<-tab[,-1:-7]

pttab<-t(data.frame(tab,check.names=F))
ttab<-data.frame(pttab)
ttab$Others<-NA
for (i in 1:nrow(ttab))
{
  ttab[i,"Others"]<-1-sum(ttab[i,1:6])
}

a<-rownames(ttab)
b<-rownames(d)

if(identical(a,b)==TRUE) {combined<-cbind(ttab, d)}

taxa<-colnames(ttab)
design<-colnames(d)

b<-melt(combined, id=design, measure=taxa)

cdata2 <- ddply(b, c("season", "substrate_type", "habitat", "cluster", "variable"), summarise,
               N    = length(value),
               mean = mean(value)
)

#Check the colour scale from COI plots
brewer.pal(colourCount, "Set2")
# "#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494"

ggplot(data=cdata2, aes(x=as.factor(cluster), y=mean, fill=variable)) + 
  geom_bar(stat="identity", size=0.05, width=1, colour="black") + 
  scale_fill_manual(breaks = c("Bacillariophyta", "Chlorophyta", "Ciliophora", "Dinozoa", "Ochrophyta", "Rhodophyta_IS", "Others"), values = c("#66C2A5", "#FC8D62", "forestgreen", "deepskyblue1", "#E78AC3", "#FFD92F", "#E5C494")) + 
  facet_grid(substrate_type+season~habitat, scale="free_x", space="free_x")+ labs(title="", x ="Cluster number", y = "Relative abundance", fill = "") + theme_bw() + 
  scale_y_continuous(limits=c(0, 1.02), expand = c(0, 0)) +
  theme(legend.position = "top", axis.text.y = element_text(size = 8), axis.text.x = element_text(angle = 90, hjust = 1, size=6, vjust=0.5), strip.text = element_text(size=10), legend.text=element_text(size=10), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm")) +
  guides(fill = guide_legend(nrow = 1))
ggsave("results/Both/18S_barplot_stacked_Phylum_nop_ncbi_5_1_22.pdf")

#Phylum - Metazoa

no<-subset_taxa(tudao, Division=="Metazoa")
datarg = transform_sample_counts(no, function(x) x/sum(x))

datag = tax_glom(datarg, "Phylum")
TopNOTUs = names(sort(taxa_sums(datag), TRUE)[1:6])
data9 = prune_taxa(TopNOTUs, datag)

d<-data.frame(sample_data(data9))

tax<-data.frame(tax_table(data9), stringsAsFactors=FALSE)
otu<-otu_table(data9)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Phylum
tab<-tab[,-1:-7]

pttab<-t(data.frame(tab,check.names=F))
ttab<-data.frame(pttab)
ttab$Others<-NA
for (i in 1:nrow(ttab))
{
  ttab[i,"Others"]<-1-sum(ttab[i,1:6])
}

a<-rownames(ttab)
b<-rownames(d)

if(identical(a,b)==TRUE) {combined<-cbind(ttab, d)}

taxa<-colnames(ttab)
design<-colnames(d)

b<-melt(combined, id=design, measure=taxa)

cdata2 <- ddply(b, c("season", "substrate_type", "habitat", "cluster", "variable"), summarise,
               N    = length(value),
               mean = mean(value)
)

# Check colour palette from COI plot
brewer.pal(7,"Paired")

ggplot(data=cdata2, aes(x=as.factor(cluster), y=mean, fill=variable)) + 
  geom_bar(stat="identity", size=0.05, width=1, colour="black") + 
  #scale_fill_manual(breaks = c("Annelida", "Arthropoda", "Chordata", "Cnidaria", "Gastrotricha", "Mollusca", "Nematoda", "Platyhelminthes", "Others"), values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#6A3D9A", "#E31A1C", "#FFFF99", "#FDBF6F", "#CAB2D6")) + 
  scale_fill_manual(breaks = c("Annelida", "Arthropoda", "Cnidaria", "Mollusca", "Nematoda", "Platyhelminthes", "Others"), values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "brown", "#FB9A99", "#FDBF6F")) + 
  facet_grid(substrate_type+season~habitat, scale="free_x", space="free_x")+ 
  labs(title="", x ="Cluster number", y = "", fill = "") + 
  theme_bw() + scale_y_continuous(limits=c(0, 1.02), expand = c(0, 0)) +
  theme(legend.position = "top", axis.text.y = element_text(size = 8), axis.text.x = element_text(angle = 90, hjust = 1, size=6, vjust=0.5), strip.text = element_text(size=10), legend.text=element_text(size=10), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"))+
  guides(fill = guide_legend(nrow = 1))
ggsave("results/Both/18S_barplot_stacked_Phylum_m_ncbi_5_1_22.pdf")
