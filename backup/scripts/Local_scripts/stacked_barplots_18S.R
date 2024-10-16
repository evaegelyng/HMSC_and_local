library("phyloseq")
library("ggplot2")
library("reshape2")
library("plyr")
library("RColorBrewer")

#Load rarefied dataset
COSQ_rare<-readRDS("../RDS/18S_no_c2_3reps_pident90.rds")

#Taxonomy plots

#For stacked plots, use "tudao" as starting object:
tudao = filter_taxa(COSQ_rare, function(x) sum(x) > 0, TRUE)
tudao = prune_samples(sample_sums(tudao)>0,tudao)

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
head(tab)
#Remove taxonomy columns
tab<-tab[,-1:-3]

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

ggplot(data=cdata2, aes(x=as.factor(cluster), y=mean, fill=variable)) + 
  geom_bar(stat="identity", linewidth=0.05, width=1, colour="black") + 
  scale_fill_manual(breaks = c("Bacillariophyta", "Cercozoa", "Chlorophyta", "Ciliophora", "Myzozoa", "PX.clade", "Others"), values = c("#66C2A5", "#FC8D62", "forestgreen", "deepskyblue1", "#E78AC3", "#FFD92F", "#E5C494")) + 
  facet_grid(substrate_type+season~habitat, scale="free_x", space="free_x")+ labs(title="", x ="", y = "Relative abundance", fill = "") + theme_bw() + 
  scale_y_continuous(limits=c(0, 1.02), expand = c(0, 0)) +
  theme(legend.position = "top", axis.title = element_text(size=16), axis.text.y = element_text(size = 9), axis.text.x = element_text(angle = 90, hjust = 1, size=8, vjust=0.5), strip.text = element_text(size=16), legend.text=element_text(size=16), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm")) +
  guides(fill = guide_legend(nrow = 1))
ggsave("../Plots/Stacked_barplots/18S_barplot_stacked_Phylum_nop.pdf")

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
tab<-tab[,-1:-3]

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
  geom_bar(stat="identity", linewidth=0.05, width=1, colour="black") + 
  scale_fill_manual(breaks = c("Annelida", "Arthropoda", "Cnidaria", "Mollusca", "Nematoda", "Platyhelminthes", "Others"), values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "brown", "#FB9A99", "#FDBF6F")) + 
  facet_grid(substrate_type+season~habitat, scale="free_x", space="free_x")+ 
  labs(title="", x ="", y = "", fill = "") + 
  theme_bw() + scale_y_continuous(limits=c(0, 1.02), expand = c(0, 0)) +
  theme(legend.position = "top", axis.title = element_text(size=16),axis.text.y = element_text(size = 9), axis.text.x = element_text(angle = 90, hjust = 1, size=8, vjust=0.5), strip.text = element_text(size=16), legend.text=element_text(size=16), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"))+
  guides(fill = guide_legend(nrow = 1))
ggsave("../Plots/Stacked_barplots/18S_barplot_stacked_Phylum_m.pdf")
