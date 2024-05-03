library("phyloseq")
library("ggplot2")
library("vegan")
library("reshape2")
library("plyr")
library("scales")
library("stringr")
library("RColorBrewer")
sessionInfo()


# Load rarefied 70%-similarity dataset (MOTU dataset)
COSQ_rare2<-readRDS("data/COSQ_rare2_correct_score_tax.rds")

# Load metadata
metadatas<-sample_data(COSQ_rare2)

#Create extra metadata variables
metadatas$sshc<-paste(metadatas$substrate_type, metadatas$season, metadatas$habitat, metadatas$cluster, sep="_")
metadatas$ssc<-paste(metadatas$substrate_type, metadatas$season, metadatas$cluster, sep="_")
metadatas$stc<-paste(metadatas$substrate_type, metadatas$cluster, sep="_")

sampledatas = sample_data(data.frame(metadatas, row.names=metadatas$root, stringsAsFactors=FALSE))

# Make new phyloseq object with extended metadata
OTU_COI = otu_table(COSQ_rare2, taxa_are_rows = TRUE)
TAX_S = tax_table(COSQ_rare2)
p_COI = phyloseq(OTU_COI, TAX_S)
COSQ_new = merge_phyloseq(p_COI, sampledatas)

# Removing cluster 2, which was only sampled in spring
COSQ_no_c2<-subset_samples(COSQ_new, !cluster==2)

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
no<-subset_taxa(tudao, !kingdom=="Metazoa")
datarg = transform_sample_counts(no, function(x) x/sum(x))

datag = tax_glom(datarg, "phylum")
TopNOTUs = names(sort(taxa_sums(datag), TRUE)[1:6])
data9 = prune_taxa(TopNOTUs, datag)

d<-data.frame(sample_data(data9))

tax<-data.frame(tax_table(data9), stringsAsFactors=FALSE)
otu<-otu_table(data9)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$phylum
tab<-tab[,-1:-9]

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

colourCount = length(unique(cdata2$variable))

ggplot(data=cdata2, aes(x=as.factor(cluster), y=mean, fill=variable)) + 
  geom_bar(stat="identity", size=0.05, width=1, colour="black") + 
  scale_fill_manual(breaks = c("Bacillariophyta", "Chlorophyta", "Discosea", "Ochrophyta", "Oomycota", "Rhodophyta", "Others"), values = brewer.pal(colourCount, "Set2")) +
  facet_grid(substrate_type+season~habitat, scale="free_x", space="free_x")+ 
  labs(title="", x ="", y = "Relative abundance", fill = "") + 
  theme_bw() + scale_y_continuous(limits=c(0, 1.02), expand = c(0, 0)) +
  theme(legend.position = "top", axis.text.y = element_text(size = 8), axis.text.x = element_text(angle = 90, hjust = 1, size=6, vjust=0.5), strip.text = element_text(size=10), legend.text=element_text(size=10), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm")) +
  guides(fill = guide_legend(nrow = 1))
ggsave("results/Both/COI_barplot_stacked_Phylum_nop_ncbi_5_1_22.pdf")

#Phylum - Metazoa

no<-subset_taxa(tudao, kingdom=="Metazoa")
datarg = transform_sample_counts(no, function(x) x/sum(x))

datag = tax_glom(datarg, "phylum")
TopNOTUs = names(sort(taxa_sums(datag), TRUE)[1:6])
data9 = prune_taxa(TopNOTUs, datag)

d<-data.frame(sample_data(data9))

tax<-data.frame(tax_table(data9), stringsAsFactors=FALSE)
otu<-otu_table(data9)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$phylum
tab<-tab[,-1:-9]

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

colourCount = length(unique(cdata2$variable))

ggplot(data=cdata2, aes(x=as.factor(cluster), y=mean, fill=variable)) + 
  geom_bar(stat="identity", size=0.05, width=1, colour="black") + 
  scale_fill_manual(values=brewer.pal(colourCount,"Paired")) +
  facet_grid(substrate_type+season~habitat, scale="free_x", space="free_x")+ 
  labs(title="", x ="", y = "", fill = "") + 
  theme_bw() + scale_y_continuous(limits=c(0, 1.02), expand = c(0, 0)) +
  theme(legend.position = "top", axis.text.y = element_text(size = 8), axis.text.x = element_text(angle = 90, hjust = 1, size=6, vjust=0.5), strip.text = element_text(size=10), legend.text=element_text(size=10), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"))+
  guides(fill = guide_legend(nrow = 1))
ggsave("results/Both/COI_barplot_stacked_Phylum_m_ncbi_5_1_22.pdf")