#FINAL MAPS OF SPECIES RICHNESS IN DK

# load libraries:
library(crosstalk)
library(leaflet)
library(tidyverse)
library(magrittr)
library(phyloseq)
library(dplyr)
library(tibble)

#### 1. Create a df with mean no. of species pr sample + coordinate data####

# Load rds file with normalized 97-data
COSQ_rare2_97 <- readRDS("../RDS/COSQ_rare2_correct_pident97.rds")
COSQ_rare2_97
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 489 taxa and 931 samples ]
#sample_data() Sample Data:       [ 931 samples by 7 sample variables ]
#tax_table()   Taxonomy Table:    [ 489 taxa by 9 taxonomic ranks ]

## Remove cluster 2 (which was only sampled in 1 season)
COSQ_no_c2<-subset_samples(COSQ_rare2_97,!cluster==2)

# Extract sample data
meta<-data.frame(sample_data(COSQ_no_c2), check.names=F)

# Load further metadata
#metadata <- read.delim("../Tekstfiler/no_control_no_sing_samples_cleaned_metadata_ASV_wise.txt")
coord <- read.delim("../Tekstfiler/merged_metadata_230427.txt",row.names=1)

# merge taxonomy and otu table from the phyloseq object
otu_tab<-merge(COSQ_no_c2@tax_table,COSQ_no_c2@otu_table,by="row.names")

# keep only final.id of the taxonomy columns
otu_tab1 <- otu_tab[,c(-1:-9)]

# Count number of species
length(unique(otu_tab$final.id))
#474

## Aggregate MOTUs from the same species, summing read counts
species_agg<-aggregate(. ~ final.id, otu_tab1, sum)
row.names(species_agg) <- species_agg$final.id
species_agg <- within(species_agg,rm(final.id))

##Convert to presence/absence
species_agg[species_agg>0]<-1

## Calculate total richness pr sample, sum of counts in each sample pr phylum
sample_rich<-species_agg %>%
  colSums
sample_rich<-data.frame(sample_rich)

## Add richness to metadata object
meta$rich<-sample_rich$sample_rich[match(row.names(sample_rich),meta$root)]

## Remove unnecessary columns
meta<-within(meta,rm("field_replicate","over_median","root","season"))

## Aggregate sample replicates, calculating the mean richness per sample replicate across seasons
meta_agg<-aggregate(. ~ cluster + habitat + substrate_type, meta, mean)


##### 2. Create a map #####
shared <- SharedData$new(meta_agg) #shared function is used to be able to fit large datasets into memory - easier to work with

## Subset metadata to spring season to get unique coordinates per cluster+habitat
coord_spring<-coord[which(coord$season=="spring"),]

## Make new copy of dataframe
both_seasons <- meta_agg
## Make column to allow matching with metadata
both_seasons$new_name<-paste("spring","_",both_seasons$cluster,"_",both_seasons$habitat,sep="")
## Add coordinates from metadata to ASV table
both_seasons$Latitude<-coord$Latitude[match(both_seasons$new_name,row.names(coord))]
both_seasons$Longitude<-coord$Longitude[match(both_seasons$new_name,row.names(coord))]
## Make new row names for combined seasons
row.names(both_seasons)<-paste("C",both_seasons$cluster,"_",both_seasons$habitat,"_",both_seasons$substrate_type,sep="")
both_seasons<-within(both_seasons,rm("new_name"))

#round up the richness column
both_seasons$richness<-round(both_seasons$rich)


#### RICHNESS MAPS - All phyla ####

#1. Substrate = sediment, habitat = rocks

#subset data
SR <- both_seasons[both_seasons$substrate_type=="sediment",]
SR <- SR[SR$habitat=="rocks",]

shared_SR <- SharedData$new(SR)

#map
map_richness_SR<- leaflet(shared_SR) %>%
  addTiles() %>%
  addCircles(~Longitude, ~Latitude, radius = ~richness * 200, # set the size of circles based on mean_richness
             popup = ~as.character(richness),
             fillColor = ~colorNumeric("YlOrRd",richness)(richness),
             fillOpacity = 0.8, stroke = TRUE,color = "black", weight=2) %>%
  addLegend("topright", pal = colorNumeric("YlOrRd", SR$richness),
            values = ~richness, title = "Mean number of species (SR)", opacity = 1)

# render the map using leafletOutput
leafletOutput("map_richness_SR", height = "100%")

print(map_richness_SR)


#2. Substrate_type = sediment, habitat=eelgrass

#subset data
SE <- both_seasons[both_seasons$substrate_type=="sediment",]
SE <- SE[SE$habitat=="eelgrass",]

#obs there are only 16 clusters with this combination (not all clusters had eelgrass)

shared_SE <- SharedData$new(SE)

#map 
map_richness_SE <- leaflet(shared_SE) %>%
  addTiles() %>%
  addCircles(~Longitude, ~Latitude, radius = ~richness * 200, # set the size of circles based on richness
             popup = ~as.character(richness),
             fillColor = ~colorNumeric("YlOrRd",richness)(richness),
             fillOpacity = 0.8, stroke = TRUE,color = "black", weight=2) %>%
  addLegend("topright", pal = colorNumeric("YlOrRd", SE$richness),
            values = ~richness, title = "Mean number of species (SE)", opacity = 1)

# render the map using leafletOutput
leafletOutput("map_richness_SE", height = "100%")

print(map_richness_SE)

#3. Substrate_type = sediment, habitat=sand

#subset data
SS <- both_seasons[both_seasons$substrate_type=="sediment",]
SS <- SS[SS$habitat=="sand",]

#obs there are only 16 clusters with this combination (not all clusters had eelgrass)

shared_SS <- SharedData$new(SS)


#map 
map_richness_SS <- leaflet(shared_SS) %>%
  addTiles() %>%
  addCircles(~Longitude, ~Latitude, radius = ~richness * 200, # set the size of circles based on richness
             popup = ~as.character(richness),
             fillColor = ~colorNumeric("YlOrRd",richness)(richness),
             fillOpacity = 0.8, stroke = TRUE,color = "black", weight=2) %>%
  addLegend("topright", pal = colorNumeric("YlOrRd", SS$richness),
            values = ~richness, title = "Mean number of species (SS)", opacity = 1)


# render the map using leafletOutput
leafletOutput("map_richness_SS", height = "100%")

print(map_richness_SS)


#4. Substrate_type = Water, habitat= rocks

#subset data
WR <- both_seasons[both_seasons$substrate_type=="water",]
WR <- WR[WR$habitat=="rocks",]

#obs there are only 16 clusters with this combination (not all clusters had eelgrass)

shared_WR <- SharedData$new(WR)


#map 
map_richness_WR <- leaflet(shared_WR) %>%
  addTiles() %>%
  addCircles(~Longitude, ~Latitude, radius = ~richness * 200, # set the size of circles based on richness
             popup = ~as.character(richness),
             fillColor = ~colorNumeric("YlOrRd",richness)(richness),
             fillOpacity = 0.8, stroke = TRUE,color = "black", weight=2) %>%
  addLegend("topright", pal = colorNumeric("YlOrRd", WR$richness),
            values = ~richness, title = "Mean number of species (WR)", opacity = 1)


# render the map using leafletOutput
leafletOutput("map_richness_WR", height = "100%")

print(map_richness_WR)


#5. Substrate_type = Water, habitat= eelgrass

#subset data
WE <- both_seasons[both_seasons$substrate_type=="water",]
WE <- WE[WE$habitat=="eelgrass",]

#obs there are only 16 clusters with this combination (not all clusters had eelgrass)

shared_WE <- SharedData$new(WE)


#map 
map_richness_WE <- leaflet(shared_WE) %>%
  addTiles() %>%
  addCircles(~Longitude, ~Latitude, radius = ~richness * 200, # set the size of circles based on richness
             popup = ~as.character(richness),
             fillColor = ~colorNumeric("YlOrRd",richness)(richness),
             fillOpacity = 0.8, stroke = TRUE,color = "black", weight=2) %>%
  addLegend("topright", pal = colorNumeric("YlOrRd", WE$richness),
            values = ~richness, title = "Mean number of species (WE)", opacity = 1)


# render the map using leafletOutput
leafletOutput("map_richness_WE", height = "100%")

print(map_richness_WE)


#6. Substrate_type = Water, habitat= sand

#subset data
WS <- both_seasons[both_seasons$substrate_type=="water",]
WS <- WS[WS$habitat=="sand",]

#obs there are only 16 clusters with this combination (not all clusters had eelgrass)

shared_WS <- SharedData$new(WS)

#map water + sand
map_richness_WS <- leaflet(shared_WS) %>%
  addTiles() %>%
  addCircles(~Longitude, ~Latitude, radius = ~richness * 200, # set the size of circles based on richness
             popup = ~as.character(richness),
             fillColor = ~colorNumeric("YlOrRd",richness)(richness),
             fillOpacity = 0.8, stroke = TRUE,color = "black", weight=2) %>%
  addLegend("topright", pal = colorNumeric("YlOrRd", WS$richness),
            values = ~richness, title = "Mean number of species (WS)", opacity = 1)


# render the map using leafletOutput
leafletOutput("map_richness_WS", height = "100%")

print(map_richness_WS)
