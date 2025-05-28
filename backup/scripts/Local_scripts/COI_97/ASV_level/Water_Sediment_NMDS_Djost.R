# libraries

library(vegan)
library(ggplot2)

# Water

four_palette <- c('#1B9E77', '#E7298A', '#A6761D', '#666666')

nmds_water <- read.table("../../../Tekstfiler/COI/COI_ASV/wat_Djost.txt", sep=",", header=TRUE)
rownames(nmds_water) <- nmds_water$X
nmds_water <- nmds_water[,-1] # remove "X" column
nmds_water[is.na(nmds_water)] = 0 # set NA values to 0

set.seed(123)
nmds_water <- metaMDS(nmds_water, k=4, maxit=999, trymax=250)
metadata<-read.table("../../../Tekstfiler/COI/COI_ASV/Mads_metadata_water2_250515.txt", sep="\t", header=TRUE)

site_scrs <- as.data.frame(scores(nmds_water, display = "sites"))
site_scrs <- cbind(site_scrs, Sub_Area = metadata$Sub_Area)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ Sub_Area, data=site_scrs, FUN = mean)
segs <- merge(site_scrs, setNames(cent, c("Sub_Area", "oNMDS1", "oNMDS2")), by = "Sub_Area", sort=FALSE)

site_scrs <- cbind(site_scrs, label = metadata$Community)

ggplot(site_scrs, aes(x=NMDS1, y=NMDS2, colour = Sub_Area)) + 
  geom_segment(data = segs, mapping = aes (xend = oNMDS1, yend= oNMDS2), alpha = 0.4) +
  scale_color_manual(values = four_palette,
                     name = "", breaks = unique(metadata$Sub_Area),
                     labels = unique(metadata$Sub_Area)) +
  geom_point(data = cent, size = 5) +
  geom_point(alpha=0.4) +
  coord_fixed() +
  theme_classic() + 
  theme(legend.position = "right", legend.text = element_text(size = 8), legend.title = element_text(size = 8), axis.text = element_text(size = 6))

ggsave("../../../Plots/NMDS/Water_NMDS_dim1_dim2.pdf", width=14, height=8)

## Now the same NMDS, but colored by salinity

## Get salinities
coord<-read.table("../../../Tekstfiler/Across_barcodes/merged_metadata_230427.txt", sep="\t", row.names=1, header=T, check.names=FALSE)
## Extract only rock and sand habitat (eelgrass sites not included in ASV-level analyses)
coord_sub<-coord[which(coord$habitat!="eelgrass"),c("Salinity","cluster")]
## Aggregate samples across seasons and habitats, keeping the mean salinity
coord_agg<-aggregate(. ~ cluster, coord_sub, mean)

## Add cluster column to site_scrs to allow matching with metadata
site_scrs$cluster <- lapply(site_scrs$label, function(x) gsub("water_","",x))
site_scrs$salinity <- coord_agg$Salinity[match(site_scrs$cluster, coord_agg$cluster)]

# set breaks and colours
breaks <- c(7,11,17,22,27,33)
breaks

cols <- RColorBrewer::brewer.pal(9, "RdYlBu")

nmds_sal<-ggplot(site_scrs, aes(x=NMDS1, y=NMDS2)) + 
  scale_colour_stepsn(colours = rev(cols),
                      breaks = breaks,
                      name = "Salinity") +
  geom_point(alpha=1, size=4, aes(colour = salinity)) +
  coord_fixed() +
  theme_classic() + 
  theme(legend.position = "right", legend.text = element_text(size = 14), legend.title = element_text(size = 14), axis.text = element_text(size = 14),axis.title.x = element_text(size = 14), axis.title.y = element_blank(),plot.margin = unit(c(0,0.4,0,0.4), 'lines'))+
  ggtitle("") +
  #set breaks on x-axis
  scale_x_continuous(breaks = c(-0.5,-0.25,0,0.25,0.5))+
  #set breaks on y-axis
  scale_y_continuous(breaks = c(-0.25,0,0.25))

saveRDS(nmds_sal,"../../../Plots/NMDS/NMDS_ASVs_sal_wat.rds")


# Sediment

nmds_sediment <- read.table("../../../Tekstfiler/COI/COI_ASV/Sed_Djost.txt", sep=",", header=TRUE)
rownames(nmds_sediment) <- nmds_sediment$X
nmds_sediment <- nmds_sediment[,-1] # remove "X" column
nmds_sediment[is.na(nmds_sediment)] = 0 # set NA values to 0

set.seed(123)
nmds_sediment <- metaMDS(nmds_sediment, k=4, maxit=999, trymax=250)
metadata<-read.table("../../../Tekstfiler/COI/COI_ASV/Mads_metadata_sed2_250515.txt", sep="\t", header=TRUE)

site_scrs <- as.data.frame(scores(nmds_sediment, display = "sites"))
site_scrs <- cbind(site_scrs, Sub_Area = metadata$Sub_Area)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ Sub_Area, data=site_scrs, FUN = mean)
segs <- merge(site_scrs, setNames(cent, c("Sub_Area", "oNMDS1", "oNMDS2")), by = "Sub_Area", sort=FALSE)

site_scrs <- cbind(site_scrs, label = metadata$Community)

ggplot(site_scrs, aes(x=NMDS1, y=NMDS2, colour = Sub_Area)) + 
  geom_segment(data = segs, mapping = aes (xend = oNMDS1, yend= oNMDS2), alpha = 0.4) +
  scale_color_manual(values = four_palette,
                     name = "", breaks = unique(metadata$Sub_Area),
                     labels = unique(metadata$Sub_Area)) +
  geom_point(data = cent, size = 5) +
  geom_point(alpha=0.4) +
  coord_fixed() +
  theme_classic() + 
  theme(legend.position = "right", legend.text = element_text(size = 8), legend.title = element_text(size = 8), axis.text = element_text(size = 6))

ggsave("../../../Plots/NMDS/Sediment_NMDS_dim1_dim2.pdf", width=14, height=8)

## Now the same NMDS, but colored by salinity
## Add cluster column to site_scrs to allow matching with metadata
site_scrs$cluster <- lapply(site_scrs$label, function(x) gsub("sediment_","",x))
site_scrs$salinity <- coord_agg$Salinity[match(site_scrs$cluster, coord_agg$cluster)]

# set breaks and colours
breaks <- c(7,11,17,22,27,33)
breaks

cols <- RColorBrewer::brewer.pal(9, "RdYlBu")

nmds_sal<-ggplot(site_scrs, aes(x=NMDS1, y=NMDS2)) + 
  scale_colour_stepsn(colours = rev(cols),
                      breaks = breaks,
                      name = "Salinity") +
  geom_point(alpha=1, size=4, aes(colour = salinity)) +
  coord_fixed() +
  theme_classic() + 
  theme(legend.position = "NA", legend.text = element_text(size = 14), legend.title = element_text(size = 14), axis.text = element_text(size = 14), axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),plot.margin = unit(c(0,0.4,0,0.4), 'lines'))+
  ggtitle("") +
  #set breaks on x-axis
  scale_x_continuous(breaks = c(-0.5,-0.25,0,0.25,0.5))+
  #set breaks on y-axis
  scale_y_continuous(breaks = c(-0.25,0,0.25))

saveRDS(nmds_sal,"../../../Plots/NMDS/NMDS_ASVs_sal_sed.rds")
