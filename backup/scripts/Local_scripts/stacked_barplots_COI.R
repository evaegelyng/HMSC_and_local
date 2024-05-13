library("phyloseq")
library("ggplot2")
library("reshape2")
library("plyr")
library("RColorBrewer")

# Load rarefied 70%-similarity dataset (MOTU dataset)
COSQ_rare<-readRDS("../RDS/COI_no_c2_3reps.rds")

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
  theme(legend.position = "top", axis.title = element_text(size=16), axis.text.y = element_text(size = 9), axis.text.x = element_text(angle = 90, hjust = 1, size=8, vjust=0.5), strip.text = element_text(size=16), legend.text=element_text(size=16), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm")) +
  guides(fill = guide_legend(nrow = 1))
ggsave("../Plots/COI_barplot_stacked_Phylum_nop.pdf")

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
  theme(legend.position = "top", axis.title = element_text(size=16), axis.text.y = element_text(size = 9), axis.text.x = element_text(angle = 90, hjust = 1, size=8, vjust=0.5), strip.text = element_text(size=16), legend.text=element_text(size=16), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"))+
  guides(fill = guide_legend(nrow = 1))
ggsave("../Plots/COI_barplot_stacked_Phylum_m.pdf")