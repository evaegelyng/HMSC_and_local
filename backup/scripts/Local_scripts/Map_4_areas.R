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

# North Sea data points
clusters_NS <- data.frame(
  station = c("1", "3", "4", "5"),
  latitude = c(55.46984, 55.98989, 57.12644, 57.59213),
  longitude = c(8.31535, 8.120517, 8.615189, 9.947258)
)

# Limfjord data points
clusters_LF <- data.frame(
  station = c("8", "9", "10"),
  latitude = c(57.0103, 56.71281, 56.55926),
  longitude = c(9.6348, 9.209116, 8.615045)
)

# Kattegat data points
clusters_KG <- data.frame(
  station = c("6", "7", "11", "12", "13", "14", "15", "16", "17", "18", "19", "26", "27", "28", "29", "30", "31"),
  latitude = c(57.31915, 57.32298, 56.7972, 56.60705, 56.44816, 56.12285, 55.70038, 55.49418, 55.4722, 55.11826, 55.06489, 55.75147, 56.04451, 55.96898, 55.68357, 55.75885, 55.34755),
  longitude = c(10.532144, 11.132922, 10.28607, 10.297121, 10.955088, 10.223454, 9.710981, 10.618222, 9.724384, 9.488224, 10.21874, 12.593566, 12.609195, 11.846132, 12.073107, 11.349378, 11.101177)
)

# Baltic Sea data points
clusters_BS <- data.frame(
  station = c("20", "21", "23", "24", "25", "32", "33"),
  latitude = c(54.88502, 54.84749, 54.9502, 55.1359, 55.56932, 55.25714, 55.02055),
  longitude = c(11.013504, 11.479693, 12.475707, 12.107394, 12.288887, 14.819494, 14.91195)
)

# Cluster 2 and 22
clusters_unused <- data.frame(
  station = c("2", "22"),
  latitude = c(55.993522, 54.757024),
  longitude = c(8.310490, 11.862938)
)
 
# Getting the projections right
coordinates(clusters_NS)<-~longitude + latitude
proj4string(clusters_NS)<-"+proj=longlat +datum=WGS84"
crs_clusters_NS <- CRS(SRS_string = "EPSG:3035")
clusters_NS<-as.data.frame(spTransform(clusters_NS, crs_clusters_NS))

coordinates(clusters_LF)<-~longitude + latitude
proj4string(clusters_LF)<-"+proj=longlat +datum=WGS84"
crs_clusters_LF <- CRS(SRS_string = "EPSG:3035")
clusters_LF<-as.data.frame(spTransform(clusters_LF, crs_clusters_LF))

coordinates(clusters_KG)<-~longitude + latitude
proj4string(clusters_KG)<-"+proj=longlat +datum=WGS84"
crs_clusters_KG <- CRS(SRS_string = "EPSG:3035")
clusters_KG<-as.data.frame(spTransform(clusters_KG, crs_clusters_KG))

coordinates(clusters_BS)<-~longitude + latitude
proj4string(clusters_BS)<-"+proj=longlat +datum=WGS84"
crs_clusters_BS <- CRS(SRS_string = "EPSG:3035")
clusters_BS<-as.data.frame(spTransform(clusters_BS, crs_clusters_BS))

coordinates(clusters_unused)<-~longitude + latitude
proj4string(clusters_unused)<-"+proj=longlat +datum=WGS84"
crs_clusters_unused <- CRS(SRS_string = "EPSG:3035")
clusters_unused<-as.data.frame(spTransform(clusters_unused, crs_clusters_unused))

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

# Plot it all with north arrow and annotation scale
denmark_plot<-ggplot(denmark) + 
  geom_sf(aes(fill = continent), color = 'black') +
  geom_point(data = clusters_NS, aes(x = coords.x1, y = coords.x2, color = "North Sea"), size = 7,alpha=0.5) + # Change color and size here
  geom_point(data = clusters_LF, aes(x = coords.x1, y = coords.x2, color = "Limfjord"), size = 7,alpha=0.5) + # Change color and size here
  geom_point(data = clusters_KG, aes(x = coords.x1, y = coords.x2, color = "Kattegat"), size = 7,alpha=0.5) + # Change color and size here
  geom_point(data = clusters_BS, aes(x = coords.x1, y = coords.x2, color = "Baltic Sea"), size = 7,alpha=0.5) + # Change color and size here
  geom_point(data = clusters_unused, aes(x = coords.x1, y = coords.x2, color = "Unused"), size = 7,stroke=2) + # Change color and size here
  geom_text(data = clusters_text, aes(x = coords.x1, y = coords.x2, label = station), size = 6) + # Change text size here
  coord_sf(crs = st_crs(3035),
           xlim = c(4200000, 4630000),
           ylim = c(3500000, 3850000)) +
  scale_fill_manual(values = c(NA, NA, NA, 'lightyellow1', NA, NA, NA, NA), guide = 'none', na.value = 'white') + # Note that the NA's specify colors of other continents, only color I chose is for Europe (other continents aren't shown anyway)
  theme(panel.background = element_rect(fill = 'mintcream'), panel.grid.major = element_line(linewidth = 0.1, color = '#80808080')) +
  xlab("") +
  ylab("") +
  annotation_north_arrow(pad_x = unit(0.2, "cm"), pad_y = unit(0.7, "cm"), height = unit(0.7, "cm"), width = unit(0.7, "cm")) + # North arrow
  annotation_scale() + # Scale bar
  scale_color_manual(name='',
                     breaks=c("Baltic Sea","Kattegat","Limfjord","North Sea"),
                     values=c('North Sea'='#666666', 'Limfjord'="#A6761D", 'Kattegat'="#E7298A", "Baltic Sea"="#1B9E77", "Unused"="white"))

ggsave("../Kort/Map_240305.png", denmark_plot, height = 16, width = 22, units = "cm")
