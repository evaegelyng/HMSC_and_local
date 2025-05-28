library("phyloseq")
library("ggplot2")
library("ggpubr")
library("reshape2")
library("plyr")
library("RColorBrewer")

#Load rarefied dataset
COSQ_rare<-readRDS("../../RDS/16S_no_c2_3reps.rds")

#Taxonomy plots

#For stacked plots, use "tudao" as starting object:
tudao = filter_taxa(COSQ_rare, function(x) sum(x) > 0, TRUE)
tudao = prune_samples(sample_sums(tudao)>0,tudao)

#Bacteria
no<-subset_taxa(tudao, Kingdom=="Bacteria")
datarg = transform_sample_counts(no, function(x) x/sum(x))

datag = tax_glom(datarg, "Family")
TopMOTUs = names(sort(taxa_sums(datag), TRUE)[1:8])
data9 = prune_taxa(TopMOTUs, datag)

d<-data.frame(sample_data(data9))

tax<-data.frame(tax_table(data9), stringsAsFactors=FALSE)
otu<-otu_table(data9)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Family
head(tab)
#Remove taxonomy columns
tab<-tab[,-1:-6]

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
taxa
design<-colnames(d)

b<-melt(combined, id=design, measure=taxa)

cdata2 <- ddply(b, c("season", "substrate_type", "habitat", "cluster", "variable"), summarise,
               N    = length(value),
               mean = mean(value)
)

bac <- ggplot(data=cdata2, aes(x=as.factor(cluster), y=mean, fill=variable)) + 
  geom_bar(stat="identity", linewidth=0.05, width=1, colour="black") + 
  facet_grid(substrate_type+season~habitat, scale="free_x", space="free_x")+ labs(title="", x ="Sampling station", y = "Relative abundance", fill = "") + theme_bw() + 
  scale_y_continuous(limits=c(0, 1.02), expand = c(0, 0)) +
  theme(legend.position = "top", axis.title = element_text(size=16), axis.text.y = element_text(size = 9), axis.text.x = element_text(angle = 90, hjust = 1, size=8, vjust=0.5), strip.text = element_text(size=16), legend.text=element_text(size=14), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm")) +
  guides(fill = guide_legend(nrow = 1))

bac

#Archaea

arch<-subset_taxa(tudao, Kingdom=="Archaea")
tudao = prune_samples(sample_sums(arch)>0,arch)
datarg = transform_sample_counts(tudao, function(x) x/sum(x))

datag = tax_glom(datarg, "Family")
TopMOTUs = names(sort(taxa_sums(datag), TRUE)[1:8])
data9 = prune_taxa(TopMOTUs, datag)

d<-data.frame(sample_data(data9))

tax<-data.frame(tax_table(data9), stringsAsFactors=FALSE)
otu<-otu_table(data9)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$Family
#Remove taxonomy columns
tab<-tab[,-1:-6]

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
taxa
design<-colnames(d)

b<-melt(combined, id=design, measure=taxa)

cdata2 <- ddply(b, c("season", "substrate_type", "habitat", "cluster", "variable"), summarise,
               N    = length(value),
               mean = mean(value)
)

arch_plot <- ggplot(data=cdata2, aes(x=as.factor(cluster), y=mean, fill=variable)) + 
  geom_bar(stat="identity", linewidth=0.05, width=1, colour="black") + 
  facet_grid(substrate_type+season~habitat, scale="free_x", space="free_x")+ 
  labs(title="", x ="Sampling station", y = "", fill = "") + 
  theme_bw() + scale_y_continuous(limits=c(0, 1.02), expand = c(0, 0)) +
  theme(legend.position = "top", axis.title = element_text(size=16),axis.text.y = element_text(size = 9), axis.text.x = element_text(angle = 90, hjust = 1, size=8, vjust=0.5), strip.text = element_text(size=16), legend.text=element_text(size=14), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"))+
  guides(fill = guide_legend(nrow = 1))

arch_plot

ggsave(file="../../Plots/Barplots/16S_stacked_fam.png",ggarrange(bac,arch_plot,nrow=2,ncol=1),width=45,height=45,units="cm")
