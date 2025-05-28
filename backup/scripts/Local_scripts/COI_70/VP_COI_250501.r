# Load libraries
library("Hmsc")
library("dplyr")
library("RColorBrewer")
library("tidyverse")

'
Variance partitioning, but now a cool function grinding out all the noise and trouble.
 What ease of use!

'


#This first function calculates the variance partitioning
#the grouping is for grouping taxa under higher taxonomic levels. 
compute_var<- function(model,tax_group=NULL,variables=NULL){
  
  VP = computeVariancePartitioning(model)
  
  #Transpose the matrix for visualisation
  VPtable<-data.frame(t(VP$vals))
  
  if(!is.null(variables)){
    if(length(variables)!=length(VPtable)){
      stop("Variable vector misspecified: length does not match")
    }
    
    
    colnames(VPtable)<-variables
  }
  
  
  VPtable$clade<-rownames(VPtable)
 
  #catching the sr2 to pop into the data frame. 
  preds = computePredictedValues(model)
  MF= evaluateModelFit(hM =model, predY = preds)
  
  sp_sr2<- as.data.frame(cbind(model$spNames,MF$SR2))
  
  colnames(sp_sr2)<- c("clade", "SR2")
  VPtable$SR2<-sp_sr2$SR2[match(VPtable$clade,sp_sr2$clade)]
  
  VPtable<- VPtable %>% na.omit()
  
  
  ## we need this table in order to make separate our data and make some calc. We also need something to
  ##plot and summarise in gpglot barplot.
  nVPtable<- VPtable
  VPtable<-within(VPtable,rm("SR2","clade"))
  VPtable<- VPtable*as.numeric(nVPtable[,"SR2"])
  VPtable$clade <- nVPtable$clade
  
  
  #the sr2 correspond to the r2 in terms of residual explanation. It is a standardised form of r2 
  #that relates to different pred etc. 
  #The residuals is the amount of variance that the model cant explain. Therefore we do the following:
  VPtable$residual<- 1-as.numeric(nVPtable$SR2)
  
  if(is.data.frame(tax_group)){
    if(any(is.na(tax_group))){
                stop("NA values are not allowed in taxonomic_grouping")}
      
   
  VPtable$group<-tax_group[match(nVPtable$clade,tax_group[,1]),2]}
  
  
  VPtable<- VPtable %>% na.omit()
  
  nVPtable<- VPtable
 
  #we pivot our data. In this format we can plot it. all the variables is popped into one column and all there
  #explanatory power into the value column

  var_len<-length(rownames(VP$vals))
 
  VPtable<- VPtable %>% pivot_longer(
      cols=c(1:var_len,"residual"),
      names_to="variables",
      values_to="value",
      names_repair="minimal"
  )
  VPtable$clade<- factor(VPtable$clade, levels=unique(VPtable$clade))

  VPtable$value<-as.numeric(VPtable$value)
  
  VPtable$clade<- factor(VPtable$clade, levels=unique(VPtable$clade))
  
  if(is.data.frame(tax_group)){
    VPtable$group<- factor(VPtable$group, levels=unique(VPtable$group))
    }
  
  return(VPtable)
    
}



#this is our plotting function

var_Plot<- function(var_table,title=NULL,x_title=NULL,y_title=NULL,grouping=NULL,color_pallette=NULL){
  var_table$variables<-
    factor(var_table$variables,levels=unique(var_table$variables))
  
  
  
  dat<- var_table %>% 
    pivot_wider(names_from = variables, values_from= value)
  
  summing<- data.frame(
    sum= numeric(length(unique(var_table$variables))))
  
  rownames(summing)<- as.character(unique(var_table$variables))
  means<- c()
  for(i in unique(var_table$variables)){
    means<-c(means,mean(dat[[i]]))
  }
  
  summing$sum<- means
  
  category<-c()
  for(i in rownames(summing)){
    if(grepl("Residual",i)){
      category=c(category,1)
    }
    else {
      if(grepl("random",i)){
        category=c(category,2)}
      
      
      else{
        category=c(category,3)}
    }}
  
  summing$category<-category 
  summing<-summing[order(summing$category,summing$sum),]
  
  var_table$variables<-factor(var_table$variables,levels=rev(rownames(summing)))
  variable<- rev(rownames(summing))
  color_pallette<- color_pallette[order(summing$category,summing$sum)]
  
    require(RColorBrewer)
    n <- 5
    qual_col_pals = brewer.pal.info[brewer.pal.info$category =="qual",]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    if(!is.null(color_pallette)){
      if(length(variable)!=length(color_pallette)){
        stop("Missing colors for color pallette")
      }
      col_vector=color_pallette
    }
    
    
    if(!is.null(grouping)){
    for(i in 2:length(colnames(grouping))){
      clade=colnames(grouping)[i]
      print(clade)
      var_table[[clade]]<- grouping[match(var_table$clade,grouping[,1]),clade]
    } 

      clade=colnames(grouping)[length(colnames(grouping))]
      
      supergroupOrder = c("Discobids","Cryptista","Archaeplastids","Haptista","SAR_Stramenopiles","SAR_Alveolates",
"SAR_Rhizarians","Amoebozoans","Breviates","Apusomonads","Opisthokonts")

      var_table<- var_table[order(factor(var_table$supergroup,levels=supergroupOrder),var_table$phylum,var_table$clade,var_table$variables),]
      
      for(i in 2:length(colnames(grouping))){
        clade=colnames(grouping)[i]
        levels <- unique(var_table[[clade]])
        var_table[[clade]]<- factor(var_table[[clade]], levels=levels)
    
      } 
      }
  var_table$clade<- factor(var_table$clade, levels=rev(unique(var_table$clade)))
  
  # Remove "_DNA.P" from class names
  var_table$clade <- sub("_[^_]+$", "", var_table$clade)
  
  var_table$variables<-factor(var_table$variables,levels=rev(levels(var_table$variables)))
  
  # The scale fill manual labels HAVE TO be organised in correspondence to the dataframe.
  plots<- ggplot(data = var_table, aes(x = clade, y = value, fill = variables)) +
    geom_bar(stat = "identity", width = 1.1, linewidth = 0) +
    theme_bw() +
    labs(title=title)+
    theme(panel.spacing = unit(0, "lines"),
          strip.background = element_blank(), 
          strip.placement = "outside",
          strip.text.x = element_text(angle=0),
          axis.text.x = element_text(angle = 0, hjust = 1, size=8, vjust=0),
          strip.text.y.left = element_text(angle = 0, size = 7),
          axis.text.y=element_text(size=7),
          legend.position="top",
          legend.box = "vertical",
          panel.grid.major.x  = element_blank(),
          legend.text = element_text(size=8)) + 
    xlab(x_title) +
    ylab(y_title)+
    labs(fill="Variables")+
    scale_y_continuous(expand = c(0,0)) + 
    scale_x_discrete(expand = c(0,0))+coord_flip()
  
  
  #creating the labels for the fill manual
  labelling<- character()
  
  for(i in variable){
    summing<- mean(dat[[i]]) 
    
    labelling <- c(labelling, paste(i, round(summing, digits = 2)))
  }
  
  
  plots<- plots+
    scale_fill_manual(
      values=col_vector,
      labels = labelling,
      breaks = variable)

    
    if(!is.null(grouping)){
      plots<- plots+ 
        facet_grid(supergroup+phylum~., 
                 scales = "free", # Let the x axis vary across facets.
                 space = "free",  # Let the width of facets vary and force all bars to have the same width.
                 switch = "y")
    }
    
    return(plots)
  }


model<- readRDS("../../RDS/fitTF_COI_sed.rds")
model_w<- readRDS("../../RDS/fitTF_COI_wat.rds")

# Load table with supergroups
sgroups <- read.table("../../Tekstfiler/Across_barcodes/Supergroups_and_cellularity.tsv", sep='\t', header=T, comment="")

# Make taxonomy table for sediment
taxonomy <- as.data.frame(colnames(model$Y))
names(taxonomy) <- "class"
taxonomy$class_name <- sub("_[^_]+$", "", taxonomy$class)

# Add new phylum names
taxonomy$new_phylum<-sgroups$new_phylum[match(taxonomy$class_name,sgroups$class)]

# Add supergroup to taxonomy table
taxonomy$supergroup <- sgroups$supergroup[match(taxonomy$new_phylum,sgroups$new_phylum)]

## Rename phylum column to match functions above
names(taxonomy)[3] <- "phylum"


# Make taxonomy table for water
taxonomy_w <- as.data.frame(colnames(model_w$Y))
names(taxonomy_w) <- "class"
taxonomy_w$class_name <- sub("_[^_]+$", "", taxonomy_w$class)

# Add new phylum names
taxonomy_w$new_phylum<-sgroups$new_phylum[match(taxonomy_w$class_name,sgroups$class)]

# Add supergroup to taxonomy table
taxonomy_w$supergroup <- sgroups$supergroup[match(taxonomy_w$new_phylum,sgroups$new_phylum)]

## Rename phylum column to match functions above
names(taxonomy_w)[3] <- "phylum"

'
This is all we need. remember if you have a higher taxa level you want to group under, insert the table.
The taxonomic grouping should be in the 2nd column and the first should be the response variable taxa group as illustrated just abve.

function(model,tax_group=NULL,variables=NULL)
The variables function specifies the names of the variables, as it can be inconvenient with salinity or random effects
long names in the model. It has to be in the right order. 
If in doubt of the variable order look at model$covNames.
ex: variables=c("Salinity","time","space").

'

variance<- compute_var(model,tax_group=taxonomy)
variance_w<-compute_var(model_w,tax_group=taxonomy_w)

unique(variance$variables)

variance <- variance %>% mutate(variables=recode(variables,"poly.Salinity..degree...2..raw...TRUE."="Salinity","cube_d14N_15N"="d15N^3",
                    "habitat"="Habitat","log_N"="log(Nitrogen)","log_TP"="log(Total phosporus)",
                    "Grain_size"="Grain size","log_Organic_content"="log(Organic content)",
                    "Random..Time_d"="Time (random)",
                   "Random..space"="Space (random)","residual"="Residual"))

unique(variance_w$variables)

variance_w <- variance_w %>% mutate(variables=recode(variables,"poly.Salinity..degree...2..raw...TRUE."="Salinity","log_Chlorophyll"="log(Chlorophyll)",
                                                 "habitat"="Habitat","log_Si"="log(Silicate)","log_PO4"="log(Phosphate)","log_DN"="log(Dissolved nitrogen)",
                                                 "Random..Time_d"="Time (random)",
                                                 "Random..space"="Space (random)","residual"="Residual"))
                                                 

'
The plotting function. 
Hvis grouping=FALSE kommer der ingen grupperede taxonomiske niveauer. 
function(var_table,title=NULL,x_title=NULL, y_title=NULL, grouping=NULL)
'

#The colour pallette for coherent variable variance partitioning colours.  
salinity1="#386CB0"
CubeN="#FDC086"
Organic_content="#BF5B17"
N="#F0027F"
Habitat="darkgreen"
Grain_size="#666666"
TP= "red"
Random_Time_d="#8DA0CB"
Random_space="black"
residual="#F2F2F2"

# The variables in the pallette currently need to be manually ordered!
color_pallette=c(salinity1,
                Habitat,
                N,
                Organic_content,
                Grain_size,
                CubeN,
                TP,
                Random_Time_d,
                Random_space,
                residual)
                 
#Water colours and pallette
Si="#BEAED4"
PO4="purple"
DN="#FFFF99"
Temperature="#B3CDE3"
Chlorophyll="#7FC97F"

# The variables in the pallette currently need to be manually ordered!

color_pallette_water=c(salinity1,
                       Habitat,
                       Si,
                       Temperature,
                       Chlorophyll,
                       PO4,
                       DN,
                       Random_Time_d,
                       Random_space,
                       residual)


sed_plot<- var_Plot(variance,title="",
x_title="",
y_title="Explained variance",
grouping=taxonomy, color_pallette=color_pallette)

ggsave("../../Plots/Var_par/COI_varpar_sed.png", plot = sed_plot, width=12, height=12)

water_plot<- var_Plot(variance_w,title="",
x_title="",
y_title="Explained variance",
grouping=taxonomy_w, color_pallette=color_pallette_water)

ggsave("../../Plots/Var_par/COI_varpar_wat.png", plot = water_plot, width=12, height=12)
