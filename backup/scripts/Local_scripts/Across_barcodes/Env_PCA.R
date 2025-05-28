#Variable plots 

# Load packages
library("tidyverse")
library("vegan")
library("fastDummies")
library(ggfortify)
library(ggpubr)

# Read variable dataframes:
dare<- read.table("../Tekstfiler/metadata_spat_extended.txt" ,sep='\t',header=TRUE)
fsed<- read.table("../Tekstfiler/sed_metadata.txt", sep="\t", header=T)
fwat<-read.table("../Tekstfiler/wat_metadata.txt", sep="\t", header=T)
#dare$Depth<-fwat$CTD_depth[match(row.names(dare),fwat$snch)]
dare$d18O<- fwat$d18O[match(row.names(dare),fwat$snch)]
dare$d2H<- fwat$d12H[match(row.names(dare),fwat$snch)]
dare$Density<- fsed$Density[match(row.names(dare),fsed$snch)]

# Rename isotope variables to accepted naming
dare$d13C<- dare$d12C_13C
dare$cube_d15N<-dare$cube_d14N_15N

# remove Cluster 2
dare<- dare %>% filter(cluster!=2)

# Remove unwanted variables for sediment PCA
PCA_df_sed<- dare  %>%
  select(-c("Longitude","Latitude","Traffic","cell","log_PO4",
            "log_Si","log_NO2","log_NO3",
            "log_NH3","log_DN","Time","Latitude.1",
            "Longitude.1","sal_group","Prop_gravel_and_coarse_sand" 
            ,"Prop_mud_and_sandy_mud", "Prop_muddy_sand", "Prop_quaternary_clay_and_silt",
            "Prop_sand", "Prop_sedimentary_rock","Prop_till_diamicton",
            "Fisheries","Traffic","Distance_To_Traffic","log_Chlorophyll","log_CN_ratio",
            "log_NP_ratio","NO2","NO3","DN","NH3","PO4","Si","cluster",
            "d18O","habitat","season","d2H","Oxygen.depletion","cube_d14N_15N",
            "Salinity","Temperature", "d12C_13C","N","C","Organic_content","Inorganic_content"
  )) %>%
  na.omit()

# Check that the correct variables remain
names(PCA_df_sed)


# Remove unwanted variables for water PCA
PCA_df_Water<-dare  %>%
  select(-c("Longitude","Latitude","Traffic","cell",
            "log_N","log_C","log_Organic_content",
            "log_Inorganic_content","log_TP",
            "log_DN","Time","Latitude.1",
            "Longitude.1","sal_group","Prop_gravel_and_coarse_sand" 
            ,"Prop_mud_and_sandy_mud", "Prop_muddy_sand", "Prop_quaternary_clay_and_silt",
            "Prop_sand", "Prop_sedimentary_rock","Prop_till_diamicton",
            "Fisheries","Traffic","Distance_To_Traffic","log_CN_ratio",
            "log_NP_ratio","cube_d14N_15N","Water_content","Dry_matter",
            "Grain_size","Inorganic_content","d12C_13C","cube_d15N",
            "cluster","C","N","Organic_content","habitat","season","DN","Oxygen.depletion",
            "Density","d12C_13C","d13C","PO4","Si","NO3","NO2","NH3"
  )) %>%
  na.omit()

# Check that the correct variables remain
names(PCA_df_Water)


##PCA plots.

#Make qualitative variables into binary data with fastdummies

PCA<- prcomp(PCA_df_sed,scale=TRUE)

summary(PCA)

pcaData <- as.data.frame(PCA$x[, 1:2])
rownames(pcaData)<- rownames(PCA_df_sed)
pcaData$season <-dare$season[match(rownames(pcaData),rownames(dare))]
pcaData$habitat <-dare$habitat[match(rownames(pcaData),rownames(dare))]

pcaVars <- PCA$rotation[,c(1,2)]
pcaLabs <- as.data.frame(pcaVars*5.5)
pcaLabs$"var" <- rownames(pcaLabs)
colnames(pcaLabs) <- c("lPC1","lPC2","var")
pcaVars <- pcaVars*5

pcaVars <- cbind(pcaVars,pcaLabs)


Habitat_colour=c("yellowgreen", "cornflowerblue","thistle3")



p<- ggplot(pcaData, aes(x=PC1,y=PC2)) +
  geom_point(size=3,aes(colour=habitat,shape=season)) +
  geom_segment(data = pcaVars, aes(x = 0, xend = PC1, y = 0, yend = PC2), arrow = arrow(length = unit(0.25, "cm")), colour = "blue", alpha=0.4) +
  geom_text(data = pcaVars, aes(x = lPC1, y = lPC2, label = var), size = 3, fontface="bold") +
  xlab("PC1 (48.5%)") +
  ylab("PC2 (13.6%)") +
  theme_classic()+
  theme(
    legend.position = ""
  )+
  scale_colour_manual(values=Habitat_colour)


#################PCA WATER

PCA<- prcomp(PCA_df_Water,scale=TRUE)

summary(PCA)

pcaData <- as.data.frame(PCA$x[, 1:2])
rownames(pcaData)<- rownames(PCA_df_Water)
pcaData$season <-dare$season[match(rownames(pcaData),rownames(dare))]
pcaData$habitat <-dare$habitat[match(rownames(pcaData),rownames(dare))]

pcaVars <- PCA$rotation[,c(1,2)]
pcaLabs <- as.data.frame(pcaVars*5.5)
pcaLabs$"var" <- rownames(pcaLabs)
colnames(pcaLabs) <- c("lPC1","lPC2","var")
pcaVars <- pcaVars*5

pcaVars <- cbind(pcaVars,pcaLabs)

#Move some labels a little bit
pcaVars$lPC2[c(9)]<-pcaVars$lPC2[c(9)]-0.2
pcaVars$lPC1[c(6)]<-pcaVars$lPC1[c(6)]-0.3
pcaVars$lPC2[c(8)]<-pcaVars$lPC2[c(8)]+0.2

q<- ggplot(pcaData, aes(x=PC1,y=PC2)) +
  geom_point(size=3,aes(colour=habitat,shape=season)) +
  geom_segment(data = pcaVars, aes(x = 0, xend = PC1, y = 0, yend = PC2), arrow = arrow(length = unit(0.25, "cm")), colour = "blue", alpha=0.4) +
  geom_text(data = pcaVars, aes(x = lPC1, y = lPC2, label = var), size = 3, fontface="bold") +
  xlab("PC1 (44.3%)") +
  ylab("PC2 (20.0%)") +
  theme_classic()+
  theme(
    legend.position = "top"
  )+
  scale_colour_manual(values=Habitat_colour)

plots<-list()
plots[[1]]<-p
plots[[2]]<-q

ggsave(file="../Plots/Env_plots/PCA.png",ggarrange(plots[[1]],plots[[2]],ncol=2,nrow=1,align="h",common.legend=TRUE,legend="top"),width=30,height=15,units="cm")


#################################
#Heatmaps
#################################

# Load libraries
library(fastDummies)
library(reshape)


# Sediment
cor.cont.data<-cor(PCA_df_sed)

# elongated the data frame for plotting in ggplot
cor.cont.data<- melt(cor.cont.data)

q<-cor.cont.data %>% ggplot(aes(x=X1,y=X2, fill=value))+
  geom_tile()+
  geom_text(size=2,aes(label = round(value, 1)))+
  xlab("")+
  ylab("")+
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") +
  coord_fixed()+
  theme(panel.spacing = unit(0, "lines"), strip.background = element_blank(),
        strip.placement = "outside", 
        axis.text.x = element_text(angle = 90, hjust = 1, size=10, vjust=0),
        strip.text.x = element_text(angle = 90, size = 6),
        panel.grid.major.y = element_blank(),
        legend.text = element_text(size=10)
        )


# Water

cor.cont.data<-cor(PCA_df_Water)

# Elongated the data frame for plotting in ggplot
cor.cont.data<- melt(cor.cont.data)


p<-cor.cont.data %>% ggplot(aes(x=X1,y=X2, fill=value))+
  geom_tile()+
  geom_text(size=2,aes(label = round(value, 1)))+
  xlab("")+
  ylab("")+
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") +
  coord_fixed()+
  theme(panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), 
        strip.placement = "outside", 
        axis.text.x = element_text(angle = 90, hjust = 1, size=10, vjust=0), 
        strip.text.x = element_text(angle = 90, size = 6),
        legend.box = "vertical",
        panel.grid.major.y = element_blank(),
        legend.text = element_text(size=10),
        legend.position="right")

ggarrange(plotlist=list(q,p), common.legend = TRUE, legend = "right")


