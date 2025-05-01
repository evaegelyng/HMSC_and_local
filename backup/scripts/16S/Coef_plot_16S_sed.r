# Activate hmsc environment:
# conda activate hmsc

rm(list = ls())

library(tidyverse)
library(patchwork)
library(fields)
library(Hmsc)
library(RColorBrewer)
library(phyloseq)

model <- readRDS("results/sediment/16S/fitTF.rds")

################################################################################
##  Specify taxonomic reference database     ###################################
################################################################################
#Load rarefied dataset
COSQ_rare<-readRDS("data/16S_no_c2_3reps.rds")
taxa <- as.data.frame(tax_table(COSQ_rare))

taxonomy <- as.data.frame(colnames(model$Y))
names(taxonomy) <- "class"
taxonomy$class_name <- sub("_[^_]+$", "", taxonomy$class)
taxonomy$kingdom<-taxa$Kingdom[match(taxonomy$class_name,taxa$Class)]
taxonomy$phylum<-taxa$Phylum[match(taxonomy$class_name,taxa$Class)]

################################################################################
##  The coefficient plot for habitat types   ###################################
################################################################################

#Credibility interval addition
postBeta = getPostEstimate(model,q =c(0.025,0.975), parName="Beta")

supportLevel=0.95
mbeta<- postBeta$mean
betaP=postBeta$support
toPlot = mbeta
# Support level operations. if the logical expression on the right is False it will set to 0. 
#The beta supportlevel paramater is whether the mean beta is significant.
toPlot = toPlot * ((betaP>supportLevel) + (betaP<(1-supportLevel))>0)

toPlot = data.frame(toPlot)

rownames(toPlot)<- colnames(model$XScaled)

#Credibility interval addition
lower_cred <- postBeta$q[1,,]* ((betaP>supportLevel) + (betaP<(1-supportLevel))>0)
upper_cred <- postBeta$q[2,,]* ((betaP>supportLevel) + (betaP<(1-supportLevel))>0)

colnames(lower_cred)<- colnames(toPlot)
colnames(upper_cred)<- colnames(toPlot)
rownames(lower_cred)<- paste(rownames(toPlot),"lower_cred",sep="_")
rownames(upper_cred)<- paste(rownames(toPlot),"upper_cred",sep="_")

toPlot <- rbind(toPlot, lower_cred,upper_cred)

toPlot$variable<- rownames(toPlot)


length(toPlot)

coef_plot<- toPlot %>% pivot_longer(
  cols=1:(length(toPlot)-1),
  values_to="Betapar",
  names_to="class",
  names_repair = "minimal"
)


coef_plot$phylum<-taxonomy$phylum[match(coef_plot$class,taxonomy$class)] 
coef_plot$kingdom<-taxonomy$kingdom[match(coef_plot$class,taxonomy$class)] 
coef_plot<- coef_plot[order(coef_plot$kingdom,coef_plot$phylum),]
coef_plot$class<- factor(coef_plot$class, levels=unique(coef_plot$class))
coef_plot$Betapar<- as.numeric(coef_plot$Betapar) 

require(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category =="qual",]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

head(coef_plot)

# Extract presence/absence data
coef_plot_pa <- coef_plot[grep(".PA", coef_plot$class), ]

# Remove ".PA" from class names
coef_plot_pa$class <- sub("_[^_]+$", "", coef_plot_pa$class)

# Plot for habitat type
hab<- coef_plot_pa %>% group_by(variable) %>% filter(str_starts(variable,"habitat",negate=F))

kingdomOrder = c("Archaea","Bacteria")

hab <- hab[order(factor(hab$kingdom,levels=kingdomOrder),hab$phylum,hab$class),]

phylumOrder <- as.vector(hab$phylum)
phylumOrder <- phylumOrder[!duplicated(phylumOrder)]

hab <- within(hab, phylum <- factor(phylum, levels = phylumOrder))

#Take credibility intervals out of "hab" object
lower_cred<- hab[grepl("lower_cred",hab$variable),]
upper_cred <-hab[grepl("upper_cred",hab$variable),]

#leave only coeficients in "hab" object
hab <-  hab[!grepl("lower_cred",hab$variable),]
hab <-  hab[!grepl("upper_cred",hab$variable),]


p <- hab%>% 
  filter(variable != "(Intercept)")%>%
  ggplot(aes(x=class, y=Betapar, color=variable)) +
  facet_grid(~phylum + kingdom, 
             scales = "free_x", # Let the x axis vary across facets.
             space = "free_x",
             switch = "x")+
  scale_color_manual(values=c("yellowgreen","cornflowerblue"))+
  theme_bw() + theme(panel.spacing = unit(0, "lines"), strip.background = element_blank(), strip.placement = "outside", axis.text.x = element_text(angle = 90, hjust = 1, size=8, vjust=0), strip.text.x = element_text(angle = 90, size = 8), legend.position="top",legend.box = "horizontal", panel.grid.major.y = element_blank(),legend.text = element_text(size=10)) +
  xlab("") +
  # add in a dotted line at zero
  geom_hline(yintercept = 0)+
  geom_point(size=3,alpha=0.9)+
  labs(y ="Estimated effect", title = "")+
  geom_errorbar(aes(ymin =lower_cred$Betapar, 
                    ymax = upper_cred$Betapar, 
                    color = variable), width = 0.5, linewidth = 0.5)

ggsave(p,file="results/sediment/16S/16S_coef_sed_hab.png",height=9,width=15)


# Extract presence only (richness) data
coef_plot_p <- coef_plot[!grepl(".PA", coef_plot$class), ]

# Remove ".P" from class names
coef_plot_p$class <- sub("_[^_]+$", "", coef_plot_p$class)

# Plot for habitat type
hab<- coef_plot_p %>% group_by(variable) %>% filter(str_starts(variable,"habitat",negate=F))

hab <- hab[order(factor(hab$kingdom,levels=kingdomOrder),hab$phylum,hab$class),]

phylumOrder <- as.vector(hab$phylum)
phylumOrder <- phylumOrder[!duplicated(phylumOrder)]

hab <- within(hab, phylum <- factor(phylum, levels = phylumOrder))

#Take credibility intervals out of "hab" object
lower_cred<- hab[grepl("lower_cred",hab$variable),]
upper_cred <-hab[grepl("upper_cred",hab$variable),]

#leave only coeficients in "hab" object
hab <-  hab[!grepl("lower_cred",hab$variable),]
hab <-  hab[!grepl("upper_cred",hab$variable),]


p <- hab%>% 
  filter(variable != "(Intercept)")%>%
  ggplot(aes(x=class, y=Betapar, color=variable)) +
  facet_grid(~phylum + kingdom, 
             scales = "free_x", # Let the x axis vary across facets.
             space = "free_x",
             switch = "x")+
  scale_color_manual(values=c("yellowgreen","cornflowerblue"))+
  theme_bw() + theme(panel.spacing = unit(0, "lines"), strip.background = element_blank(), strip.placement = "outside", axis.text.x = element_text(angle = 90, hjust = 1, size=8, vjust=0), strip.text.x = element_text(angle = 90, size = 8), legend.position="top",legend.box = "horizontal", panel.grid.major.y = element_blank(),legend.text = element_text(size=10)) +
  xlab("") +
  # add in a dotted line at zero
  geom_hline(yintercept = 0)+
  geom_point(size=3,alpha=0.9)+
  labs(y ="Estimated effect", title = "")+
  geom_errorbar(aes(ymin =lower_cred$Betapar, 
                    ymax = upper_cred$Betapar, 
                    color = variable), width = 0.5, linewidth = 0.5)

ggsave(p,file="results/sediment/16S/16S_coef_sed_hab_rich.png",height=9,width=15)

# Plot for salinity
sal <- coef_plot_p %>% group_by(variable) %>% filter(str_starts(variable,"poly",negate=F))

sal <- sal[order(factor(sal$kingdom,levels=kingdomOrder),sal$phylum,sal$class),]

phylumOrder <- as.vector(sal$phylum)
phylumOrder <- phylumOrder[!duplicated(phylumOrder)]

sal <- within(sal, phylum <- factor(phylum, levels = phylumOrder))

#Take credibility intervals out of "hab" object
lower_cred<- sal[grepl("lower_cred",sal$variable),]
upper_cred <-sal[grepl("upper_cred",sal$variable),]

#leave only coeficients in "hab" object
sal <-  sal[!grepl("lower_cred",sal$variable),]
sal <-  sal[!grepl("upper_cred",sal$variable),]


p <- sal%>% 
  filter(variable != "(Intercept)")%>%
  ggplot(aes(x=class, y=Betapar, color=variable)) +
  facet_grid(~phylum + kingdom, 
             scales = "free_x", # Let the x axis vary across facets.
             space = "free_x",
             switch = "x")+
  scale_color_manual(values=c("red","black"))+
  theme_bw() + theme(panel.spacing = unit(0, "lines"), strip.background = element_blank(), strip.placement = "outside", axis.text.x = element_text(angle = 90, hjust = 1, size=8, vjust=0), strip.text.x = element_text(angle = 90, size = 8), legend.position="top",legend.box = "horizontal", panel.grid.major.y = element_blank(),legend.text = element_text(size=10)) +
  xlab("") +
  # add in a dotted line at zero
  geom_hline(yintercept = 0)+
  geom_point(size=3,alpha=0.9)+
  labs(y ="Estimated effect", title = "")+
  geom_errorbar(aes(ymin =lower_cred$Betapar, 
                    ymax = upper_cred$Betapar, 
                    color = variable), width = 0.5, linewidth = 0.5)

ggsave(p,file="results/sediment/16S/16S_coef_sed_sal_rich.png",height=9,width=15)