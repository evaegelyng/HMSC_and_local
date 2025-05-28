library(rnaturalearth)
library(rnaturalearthhires)
library(sf)
library(sp)
library(ggplot2)
library(ggspatial)
library(cowplot)
library(grid)

denmark <- ne_countries(scale = "large", returnclass = 'sf')
#denmark <- ne_countries(country = 'denmark', scale="large", returnclass = 'sf') - Use this instead if you don't want Sweden/Germany shown

# Eelgrass data points
clusters_E <- data.frame(
  station = c("8", "9", "10", "16", "17", "19","20", "21", "24", "25", "26", "27", "28", "29", "30", "31"),
  latitude = c(57.0103, 56.71281, 56.55926, 55.49418, 55.4722, 55.06489, 54.88502, 54.84749, 55.1359, 55.56932, 55.75147, 56.04451, 55.96898, 55.68357, 55.75885, 55.34755),
  longitude = c(9.6348, 9.209116, 8.615045, 10.618222, 9.724384, 10.21874, 11.013504, 11.479693, 12.107394, 12.288887, 12.593566, 12.609195, 11.846132, 12.073107, 11.349378, 11.101177)
)

# Non-eelgrass data points
clusters_NE <- data.frame(
  station = c("1", "3", "4", "5", "6", "7", "11", "12", "13", "14", "15", "18", "23", "32", "33"),
  latitude = c(55.46984, 55.98989, 57.12644, 57.59213, 57.31915, 57.32298, 56.7972, 56.60705, 56.44816, 56.12285, 55.70038, 55.11826, 54.9502, 55.25714, 55.02055),
  longitude = c(8.31535, 8.120517, 8.615189, 9.947258, 10.532144, 11.132922, 10.28607, 10.297121, 10.955088, 10.223454, 9.710981, 9.488224, 12.475707, 14.819494, 14.91195)
)

# All used sampling stations
clusters_ALL <- rbind(clusters_E, clusters_NE)


# Clusters 2 and 22
clusters_unused <- data.frame(
  station = c("2", "22"),
  latitude = c(55.993522, 54.757024),
  longitude = c(8.310490, 11.862938)
)
 
# Getting the projections right
coordinates(clusters_E)<-~longitude + latitude
proj4string(clusters_E)<-"+proj=longlat +datum=WGS84"
crs_clusters_E <- CRS(SRS_string = "EPSG:3035")
clusters_E<-as.data.frame(spTransform(clusters_E, crs_clusters_E))

coordinates(clusters_NE)<-~longitude + latitude
proj4string(clusters_NE)<-"+proj=longlat +datum=WGS84"
crs_clusters_NE <- CRS(SRS_string = "EPSG:3035")
clusters_NE<-as.data.frame(spTransform(clusters_NE, crs_clusters_NE))

coordinates(clusters_unused)<-~longitude + latitude
proj4string(clusters_unused)<-"+proj=longlat +datum=WGS84"
crs_clusters_unused <- CRS(SRS_string = "EPSG:3035")
clusters_unused<-as.data.frame(spTransform(clusters_unused, crs_clusters_unused))

coordinates(clusters_ALL)<-~longitude + latitude
proj4string(clusters_ALL)<-"+proj=longlat +datum=WGS84"
crs_clusters_ALL <- CRS(SRS_string = "EPSG:3035")
clusters_ALL<-as.data.frame(spTransform(clusters_ALL, crs_clusters_ALL))

# clusters_text - this is just if you want to be able to change where the text is placed (currently inside the circles) 
clusters_text <- data.frame(
  station = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33"),
  latitude = c(55.46984, 55.993522, 55.98989, 57.12644, 57.59213, 57.31915, 57.32298, 57.0103, 56.71281, 56.55926, 56.7972, 56.60705, 56.44816, 56.12285, 55.70038, 55.49418, 55.4722, 55.11826, 55.06489, 54.88502, 54.84749, 54.757024, 54.9502, 55.1359, 55.56932, 55.75147, 56.04451, 55.96898, 55.68357, 55.75885, 55.34755, 55.25714, 55.02055),
  longitude = c(8.31535, 8.310490, 8.120517, 8.615189, 9.947258, 10.532144, 11.132922, 9.6348, 9.209116, 8.615045, 10.28607, 10.297121, 10.955088, 10.223454, 9.710981, 10.618222, 9.724384, 9.488224, 10.21874, 11.013504, 11.479693, 11.862938, 12.475707, 12.107394, 12.288887, 12.593566, 12.609195, 11.846132, 12.073107, 11.349378, 11.101177, 14.819494, 14.91195)
)

# Getting text projection right
coordinates(clusters_text)<-~longitude + latitude
proj4string(clusters_text)<-"+proj=longlat +datum=WGS84"
crs_clusters_text <- CRS(SRS_string = "EPSG:3035")
clusters_text<-as.data.frame(spTransform(clusters_text, crs_clusters_text))

#Add asterisks to stations where sampling was not complete
clusters_unused$station <- c("2*","22**")
clusters_text$station[2] <- " 2*"
clusters_text$station[22] <- "   22**"

# Plot it all with north arrow and annotation scale
denmark_plot<-ggplot(denmark) +
  geom_sf(aes(fill = continent), color = 'black') +
  geom_point(data = clusters_ALL, aes(x = coords.x1, y = coords.x2, color = "Sampled locations"), size = 8,alpha=0.8) + # Change color and size here
  geom_point(data = clusters_unused, aes(x = coords.x1, y = coords.x2, color = "Planned, but not sampled"), size = 8,alpha=0.8) + # Change color and size here
  geom_text(data = clusters_text, aes(x = coords.x1, y = coords.x2, label = station), size = 5, fontface="bold") +
  coord_sf(crs = st_crs(3035),
           xlim = c(4200000, 4630000),
           ylim = c(3500000, 3850000)) +
  scale_fill_manual(values = c(NA, NA, NA, '#dfd294', NA, NA, NA, NA), guide = 'none', na.value = 'white') + # Note that the NA's specify colors of other continents, only color I chose is for Europe (other continents aren't shown anyway)
  theme(panel.background = element_rect(fill = '#dceced'), panel.grid.major = element_line(linewidth = 0.1, color = '#80808080')) +
  xlab("") +
  ylab("") +
  annotation_north_arrow(pad_x = unit(0.2, "cm"), pad_y = unit(0.7, "cm"), height = unit(0.7, "cm"), width = unit(0.7, "cm")) + # North arrow
  annotation_scale() + # Scale bar
  scale_color_manual(name='',
                     breaks=c("Sampled locations","Planned, but not sampled"),
                     values=c("Sampled locations"="white", "Planned, but not sampled"="white"))+ 
  theme(legend.position = "none", text = element_text(size = 18),legend.background=element_blank())

denmark_plot

ggsave("../Images_figure_1/Map_hab_emil.png", denmark_plot, height = 16, width = 22, units = "cm")
