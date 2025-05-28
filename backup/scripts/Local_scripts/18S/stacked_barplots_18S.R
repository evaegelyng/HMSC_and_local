library("phyloseq")
library("ggplot2")
library("reshape2")
library("plyr")
library("RColorBrewer")
library("ggpubr")

#Load rarefied dataset
COSQ_rare<-readRDS("../../RDS/18S_no_c2_3reps.rds")

## Get cleaned OTU table
otus<-data.frame(otu_table(COSQ_rare),check.names=F)

## Get sample data
sam<-data.frame(sample_data(COSQ_rare),check.names=F)

## Get taxonomy
tax<-data.frame(tax_table(COSQ_rare))

# Load table with corrected phylum names
tax_cur<-read.table("../../Tekstfiler/18S/18S_classified_phy_class_curated.tsv",sep="\t",header=T)

## Add curated phylum names
tax$new_phylum<-tax_cur$new_phylum[match(tax$class,tax_cur$class)]

## Import overview of uni-/multicellularity for eukaryotes
cell<-read.table("../../Tekstfiler/Across_barcodes/Supergroups_and_cellularity.tsv",sep="\t", header=T)

# Add cellularity
tax$uni_multi <- cell$uni_multi[match(tax$class,cell$class)]

# Extract needed columns
tax <- tax[,c("uni_multi","new_phylum")]

###Make new phyloseq object
tax_mat_b<-as.matrix(tax)
OTU = otu_table(otus, taxa_are_rows = FALSE)
TAX_b = tax_table(tax_mat_b)
SAM = sample_data(sam)
p_NCBI = phyloseq(OTU, TAX_b, SAM)

#Taxonomy plots

#For stacked plots, use "tudao" as starting object:
tudao = filter_taxa(p_NCBI, function(x) sum(x) > 0, TRUE)
tudao = prune_samples(sample_sums(tudao)>0,tudao)

#Unicellular
uni<-subset_taxa(tudao, uni_multi=="U")
datarg = transform_sample_counts(uni, function(x) x/sum(x))

datag = tax_glom(datarg, "new_phylum")
TopMOTUs = names(sort(taxa_sums(datag), TRUE)[1:6])
data9 = prune_taxa(TopMOTUs, datag)

d<-data.frame(sample_data(data9))

tax<-data.frame(tax_table(data9), stringsAsFactors=FALSE)
otu<-otu_table(data9)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$new_phylum
head(tab)
#Remove taxonomy columns
tab<-tab[,-1:-2]

pttab<-t(data.frame(tab,check.names=F))
ttab<-data.frame(pttab,check.names=F)
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

uni_plot<-ggplot(data=cdata2, aes(x=as.factor(cluster), y=mean, fill=variable)) + 
  geom_bar(stat="identity", linewidth=0.05, width=1, colour="black") + 
  scale_fill_manual(breaks = c("Bacillariophyceae*","Cercozoa","Ciliophora","Coscinodiscophyceae*","Dinophyceae*", "Fragilariophyceae*","Others"), values = c("deepskyblue1","#66C2A5", "forestgreen","#FC8D62","#FFD92F","brown","#E5C494")) + 
  facet_grid(substrate_type+season~habitat, scale="free_x", space="free_x")+ labs(title="", x ="", y = "Relative abundance", fill = "") + theme_bw() + 
  scale_y_continuous(limits=c(0, 1.02), expand = c(0, 0)) +
  theme(legend.position = "top", axis.title = element_text(size=16), axis.text.y = element_text(size = 9), axis.text.x = element_text(angle = 90, hjust = 1, size=8, vjust=0.5), strip.text = element_text(size=16), legend.text=element_text(size=14), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm")) +
  guides(fill = guide_legend(nrow = 2, byrow=TRUE))

uni_plot

#Multicellular (incl. classes that include both unicellular and multicellular groups)

multi<-subset_taxa(tudao, uni_multi!="U")
datarg = transform_sample_counts(multi, function(x) x/sum(x))

datag = tax_glom(datarg, "new_phylum")
TopMOTUs = names(sort(taxa_sums(datag), TRUE)[1:6])
data9 = prune_taxa(TopMOTUs, datag)

d<-data.frame(sample_data(data9))

tax<-data.frame(tax_table(data9), stringsAsFactors=FALSE)
otu<-otu_table(data9)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)
rownames(tab)<-tab$new_phylum
tab2<-tab[,-1:-2]

pttab<-t(data.frame(tab2,check.names=F))
ttab<-data.frame(pttab,check.names=F)
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

## Add phylum names
cdata2$phylum<-tab$new_phylum[match(cdata2$variable,row.names(tab))]
cdata2$phylum[cdata2$variable=="Others"]<-"Others"

colourCount = length(unique(cdata2$variable))
brewer.pal(colourCount,"Paired")

m<-ggplot(data=cdata2, aes(x=as.factor(cluster), y=mean, fill=phylum)) + 
  geom_bar(stat="identity", linewidth=0.05, width=1, colour="black") + 
  scale_fill_manual(breaks = c("Annelida", "Arthropoda", "Mollusca", "Nematoda", "Phaeophyceae*","Platyhelminthes","Others"), values = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C","#FF7F00", "brown", "#E5C494")) + 
  facet_grid(substrate_type+season~habitat, scale="free_x", space="free_x")+ 
  labs(title="", x ="", y = "", fill = "") + 
  theme_bw() + scale_y_continuous(limits=c(0, 1.02), expand = c(0, 0)) +
  theme(legend.position = "top", axis.title = element_text(size=16),axis.text.y = element_text(size = 9), axis.text.x = element_text(angle = 90, hjust = 1, size=8, vjust=0.5), strip.text = element_text(size=16), legend.text=element_text(size=14), axis.ticks.length=unit(.04, "cm"), legend.key.size = unit(0.4, "cm"))+
  guides(fill = guide_legend(nrow = 2, byrow=TRUE))

m

ggsave(file="../../Plots/Barplots/18S_stacked.png",ggarrange(uni_plot,m,nrow=1,ncol=2),width=45,height=15,units="cm")
