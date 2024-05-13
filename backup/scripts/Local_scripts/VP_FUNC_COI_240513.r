# Load libraries
library("Hmsc")
library("dplyr")
library("RColorBrewer")
library("tidyverse")
library("phyloseq")

'
Variance partitioning, but now a cool function grinding out all the noise and trouble.
 What ease of use!

'



#THis first function calculates the variance partitioning, flips and flops the table for plotting entrance
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
  VPtable<-within(VPtable,rm("SR2","clade"))
  #VPtable<-nVPtable %>% select(!c("SR2","clade"))
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



#this is our plotting function, it should not be a hustle to quickly go into the function and add a ylab or something extra if need be :-)

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
      var_table[[clade]]<- tax[match(var_table$clade,grouping[,1]),clade]
    } 

      clade=colnames(grouping)[length(colnames(grouping))]
      
      var_table<- var_table[order(var_table$kingdom,var_table$phylum,var_table$class,var_table$clade,var_table$variables),]
      
  
      for(i in 2:length(colnames(grouping))){
        clade=colnames(grouping)[i]
        levels <- unique(var_table[[clade]])
        var_table[[clade]]<- factor(var_table[[clade]], levels=levels)
    
      } 
      }
    var_table$clade<- factor(var_table$clade, levels=rev(unique(var_table$clade)))
  
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
        facet_grid(kingdom+phylum+class~., 
                 scales = "free", # Let the x axis vary across facets.
                 space = "free",  # Let the width of facets vary and force all bars to have the same width.
                 switch = "y")
    }
    
    return(plots)
  }


model<- readRDS("../RDS/COI_240426.rds")
model_w<- readRDS("../RDS/COI_wat_240426.rds")

COSQ_rare2<-readRDS("../RDS/COI_no_c2_3reps.rds")
OTU_COI = otu_table(COSQ_rare2, taxa_are_rows = TRUE)
TAX_S = tax_table(COSQ_rare2)
tax <- data.frame(TAX_S,rownames=F)
tax<- tax[,c("order","class","phylum","kingdom")]

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

variance <- variance %>% mutate(variables=recode(variables,"poly.Salinity..degree...2..raw...TRUE."="Salinity","d14N_15N"="d15N",
                                                 "habitat"="Habitat","N"="Nitrogen","TP"="Total phosporus",
                                                 "Grain_size"="Grain size","Organic_content"="Organic content","Oxygen.depletion"="Oxygen depletion",
                                                 "bac_rich"="Bacterial ASV richness","FishingTrawling"="Fishing","Random..Time_d"="Time (random)",
                                                 "Random..space"="Space (random)","residual"="Residual"))

unique(variance_w$variables)

variance_w <- variance_w %>% mutate(variables=recode(variables,"poly.Salinity..degree...2..raw...TRUE."="Salinity","d14N_15N"="d15N",
                                                     "habitat"="Habitat","Si"="Silicate","PO4"="Phosphate","DN"="Dissolved nitrogen",
                                                     "Oxygen.depletion"="Oxygen depletion","bac_rich"="Bacterial ASV richness","FishingTrawling"="Fishing","Random..Time_d"="Time (random)",
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
Habitat="#B3CDE3"
Grain_size="#666666"
TP= "red"
Oxygen_Depletion= "#FFFF99"
Bac_rich = "darkgreen"
Fishing = "purple"  
Random_Time_d="#8DA0CB"
Random_space="black"
residual="#F2F2F2"

color_pallette=c(salinity1,
                 Habitat,
                 Bac_rich,
                 Organic_content,
                 N,
                 Grain_size,
                 CubeN,
                 TP,
                 Fishing,
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
                 Habitat,
                 Fishing,
                 Temperature,
                 PO4,
                 Chlorophyll,
                 Si,
                 DN,
                 Oxygen_Depletion,
                 Bac_rich,
                 Random_Time_d,
                 Random_space,
                 residual)

the_plot<- var_Plot(variance,title="",
x_title="",
y_title="Explained variance",
grouping=tax, color_pallette=color_pallette)

ggsave("../Plots/COI_varpar_sed_240426.png", plot = the_plot, width=14, height=9)

water_plot<- var_Plot(variance_w,title="",
x_title="",
y_title="Explained variance",
grouping=tax, color_pallette=color_pallette_water)

ggsave("../Plots/COI_varpar_wat_240426.png", plot = water_plot, width=14, height=9)
