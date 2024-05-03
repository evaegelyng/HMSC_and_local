# Load libraries
library("Hmsc")
library("vegan")
library("reshape2")
library("plyr")
library("scales")
library("stringr")
library("RColorBrewer")
library("corrplot")
library("tidyverse")
library("phyloseq")

'
Variance partitioning, but now a cool function grinding out all the noise and trouble.
 What ease of use!.

'



#THis first function calcultes the variance partitioning, flips and flops the table for plotting entrance
#the grouping is basically for grouping taxa under higher levels(etc. phylum). 
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
  
  
  #colnames(VPtable)<-c("salinity","cube_d14N_15N","log_Organic_content","log_N","habitat",                                                                "Grain_size","log_TP","Fisheries","Distance_To_Traffic","Oxygen.depletion","time","space")

  VPtable$clade<-rownames(VPtable)
 
  #catching the sr2 to pop into the data frame. 
  preds = computePredictedValues(model)
  MF= evaluateModelFit(hM =model, predY = preds)
  
  sp_sr2<- as.data.frame(cbind(model$spNames,MF$SR2))
  
  colnames(sp_sr2)<- c("clade", "SR2")
  VPtable$SR2<-sp_sr2$SR2[match(VPtable$clade,sp_sr2$clade)]
  
  VPtable<- VPtable %>% na.omit()
  
  
  ## we need this table in order to make seperate our data and make some calc. We also need something to
  ##plot and summarise in gpglot barplot.
  nVPtable<- VPtable
  VPtable<-nVPtable %>% select(!c("SR2","clade"))
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
# VPtable<-  VPtable[order(VPtable$class),]
  VPtable$clade<- factor(VPtable$clade, levels=unique(VPtable$clade))

  VPtable$value<-as.numeric(VPtable$value)
  
  VPtable$clade<- factor(VPtable$clade, levels=unique(VPtable$clade))
  
  
  
  if(is.data.frame(tax_group)){
    VPtable$group<- factor(VPtable$group, levels=unique(VPtable$group))
    }
  
  return(VPtable)
    
}



#this is our plotting function, it should not be a hustle to quickly go into the function and add a ylab or something extra if need be :-)

var_Plot<- function(var_table,title=NULL,x_title=NULL,y_title=NULL,grouping=NULL,color_pallette=NULL){
    var_table$variables<- factor(var_table$variables,levels=unique(var_table$variables))
    variable<- unique(var_table$variables)
    
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
      var_table[[clade]]<- tax[match(var_table$clade,grouping[,1]),clade]
    } 

      clade=colnames(grouping)[length(colnames(grouping))]
      
      var_table<- var_table[order(var_table$Division,var_table$Phylum),]
      
      var_table$clade<- factor(var_table$clade, levels=unique(var_table$clade))
  
      for(i in 2:length(colnames(grouping))){
        clade=colnames(grouping)[i]
        levels <- unique(var_table[[clade]])
        var_table[[clade]]<- factor(var_table[[clade]], levels=levels)
    
      } 
      }
  
  # The scale fill manual labels HAVE TO be organised in correspondance to the dataframe.
    plots<- ggplot(data = var_table, aes(x = clade, y = value, fill = variables)) +
    geom_bar(stat = "identity", width = 1.1, linewidth = 0) +
    theme_bw() +
    labs(title=title)+
    theme(panel.spacing = unit(0, "lines"),
    strip.background = element_blank(), 
    strip.placement = "outside",
    axis.text.x = element_text(angle = 90, hjust = 1, size=7, vjust=0),
    strip.text.x = element_text(angle = 90, size = 8),
    legend.position="top",
    legend.box = "horizontal",
    panel.grid.major.y = element_blank(),
    legend.text = element_text(size=10)) + 
      xlab(x_title) +
      ylab(y_title)+
    scale_y_continuous() + 
      scale_x_discrete(expand = c(0,0))


    #creating the labels for the fill manual
    labelling<- character()
    dat<- var_table %>% 
      pivot_wider(names_from = variables, values_from= value)
  
  for(i in variable){
    summing<- mean(dat[[i]]) 

    labelling <- c(labelling, paste(i, round(summing, digits = 2)))
  }
    
    
    plots<- plots+
      scale_fill_manual(
        values=col_vector,
        labels = labelling)
      
    
    if(!is.null(grouping)){
      plots<- plots+ 
        facet_grid(~ Phylum+Division,, 
                 scales = "free_x", # Let the x axis vary across facets.
                 space = "free_x",  # Let the width of facets vary and force all bars to have the same width.
                 switch = "x")
    }
    
    return(plots)
  }


model<- readRDS("results/sediment/18S/18S_rich_p50prevclass_sediment_plus_c2_thin10_s6000_tr02_untrans.rds")
model_w<- readRDS("results/Water/18S/18S_rich_p50prevclass_water_plus_c2_thin10_s6000_tr02_untrans.rds")

taxonomy_ncbi<-read.table("data/f_tax_ncbi_5_01_22.txt", sep='\t', header=T, comment="")

tax<- taxonomy_ncbi[,c("Class","Phylum","Division")]

'
This is all we need. remember if you have a higher taxa level you wanna group under insert the table.
The taxonomic grouping should be in the 2nd column and the first should be the responce variable taxa group as illustrated just abve.

function(model,tax_group=NULL,variables=NULL)
The variables function specifies the names of the variables, as it can be inconvenient with salinity or random effects
long names in the model. It has to be in the right order. 
If in doubt of the variable order look at model$covNames.
ex: variables=c("Salinity","time","space").

'

variance<- compute_var(model,tax_group=tax)
variance_w<-compute_var(model_w,tax_group=tax)


'
The plotting function. 
Hvis grouping=FALSE kommer der ingen grupperede taxonomiske niveauer. 
function(var_table,title=NULL,x_title=NULL, y_title=NULL, grouping=NULL)
'

colnames(variance)
colnames(variance_w)

#The colour pallette for coherent variable variance partitioning colours.  
salinity1="#386CB0"
CubeN="#FDC086"
Organic_content="#BF5B17"
N="#F0027F"
Habitat="#B3CDE3"
Grain_size="#666666"
TP= "red"
Oxygen_Depletion= "#FFFF99"
Random_Time_d="#8DA0CB"
Random_space="#F2F2F2"
residual="black"

color_pallette=c(salinity1,
                 CubeN,
                 N,
                 Organic_content,
                 Habitat,
                 Grain_size,
                 TP,
                 Oxygen_Depletion,
                 Random_Time_d,
                 Random_space,
                 residual)
                 
#Water colours and pallette
Si="#BEAED4"
PO4="#FDC086"
DN="#666666"
Temperature="red"
Chlorophyll="#7FC97F"

color_pallette_water=c(salinity1,
                 Si,
                 PO4,
                 DN,
                 Temperature,
                 Habitat,
                 Chlorophyll,
                 Oxygen_Depletion,
                 Random_Time_d,
                 Random_space,
                 residual)

the_plot<- var_Plot(variance,title="Sediment functional group variance partitioning",
x_title="",
y_title="Explained variance",
grouping=tax, color_pallette=color_pallette)

ggsave("results/sediment/18S/18S_varpar_sed_240318.png", plot = the_plot, width=10, height=9)

water_plot<- var_Plot(variance_w,title="Water functional group variance partitioning",
x_title="",
y_title="Explained variance",
grouping=tax, color_pallette=color_pallette_water)

ggsave("results/Water/18S/18S_varpar_wat_240318.png", plot = water_plot, width=10, height=9)