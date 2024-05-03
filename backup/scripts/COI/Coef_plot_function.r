library(tidyverse)
library(patchwork)
library(fields)
library(Hmsc)


model<-  readRDS("C:/Users/zoega/OneDrive/Skrivebord/Kandidat_3/projekt/results/COI_rich_sed_FIVAL.rds")
model_water<- readRDS("C:/Users/zoega/OneDrive/Skrivebord/Kandidat_3/projekt/results/COI_rich_p50prevorder_water_FINAL.rds")


coef_beta<- function(model){
  
  postBeta = getPostEstimate(model, parName="Beta")
  
  supportLevel=0.95
  mbeta<- postBeta$mean
  betaP=postBeta$support
  toPlot = mbeta

  toPlot = toPlot * ((betaP>supportLevel) + (betaP<(1-supportLevel))>0)
  toPlot=data.frame(toPlot)
  rownames(toPlot)<- colnames(model$XScaled)
  toPlot$variable<- rownames(toPlot) 
 
  
  len<-length(toPlot)-1
  
  coef_plot<- toPlot %>% pivot_longer(
    cols=1:len,
    values_to="Betapar",
    names_to="Species",
    names_repair = "minimal"
  )
  
  
  coef_plot$Betapar<- as.numeric(coef_plot$Betapar) 
  #coef_plot$Species<- factor(coef_plot$Species, levels=unique(coef_plot$Species))
 
  return(coef_plot)
}


coef_plot<- function(beta,parameter=NULL,grouping=NULL,title=NULL){
  require(RColorBrewer)
  n <- 9
  qual_col_pals = brewer.pal.info[brewer.pal.info$category =="qual",]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  beta<-beta %>% data.frame() %>% filter(variable==parameter)

  if(!is.null(grouping)){
    for(i in 2:length(colnames(grouping))){
      clade=colnames(grouping)[i]
      beta[[clade]]<- tax[match(beta$Species,grouping[,1]),clade]
    } 

      clade=colnames(grouping)[length(colnames(grouping))]
      
      beta<- beta[order(beta$kingdom,beta$phylum,beta$class),]
      
      beta$Species<- factor(beta$Species, levels=unique(beta$Species))
      
      for(i in 2:length(colnames(grouping))){
        clade=colnames(grouping)[i]
        levels <- unique(beta[[clade]])
        beta[[clade]]<- factor(beta[[clade]], levels=levels)
      } 
      }
     beta$Betapar=round(as.numeric(beta$Betapar),2)

  
  plo<-beta %>% 
    filter(variable != "(Intercept)")%>%
    ggplot(aes(x=Species, y=Betapar))+
    scale_color_manual(values=col_vector[11:23])+
    theme_bw()+
    scale_y_continuous(labels = scales::label_number(accuracy = 1))  +
    theme(panel.spacing = unit(0, "lines"), 
          strip.background = element_blank(), 
          strip.placement = "outside",
          axis.text.x = element_text(angle = 90, hjust = 1, size=8, vjust=0),
          strip.text.x = element_text(angle = 90, size = 8),
          legend.position="top",
          legend.box = "horizontal",
          panel.grid.major.y = element_blank(),
          legend.text = element_text(size=10)) +
    xlab("Kingdom") +
    geom_point(size=3,alpha=0.9)+
    # add in a dotted line at zero
    geom_hline(yintercept = 0) +
    labs(
      y ="Estimate of effect of environmental predictor",
      title = paste(parameter,title, sep=" "))
  

  if(!is.null(grouping)){
  plo<- plo+ 
      facet_grid(~ class+phylum+kingdom, 
                 scales = "free_x", # Let the x axis vary across facets.
                 space = "free_x",  # Let the width of facets vary and force all bars to have the same width.
                 switch = "x")}    # }
  
  return(plo)
  
}



'
Specify taxonomic refference database:
'

COSQ_rare2<-readRDS("C:/Users/zoega/OneDrive/Skrivebord/Kandidat_3/projekt/Data/COSQ_rare2_correct_score_tax.rds")
OTU_COI = otu_table(COSQ_rare2, taxa_are_rows = TRUE)
TAX_S = tax_table(COSQ_rare2)
tax <- data.frame(TAX_S,rownames=F)
tax<- tax[,c("order","class","phylum","kingdom")]


#This function estimates beta coefficients and formats data to plot with ggplot
beta<-coef_beta(model)

#EVA REMEMBER TO CHANGE KINGDOM IN THE FUNCTION TO division and to take the 18s dataframe.
####OBS de forskellige parametre har forskellige units!! hvor højt du vil afrunde kan ændres i round(), 
#og hvor mange descimaler i plottet på y aksen kan ændres under scale y continous.  
#ploting function that filters after specified parameter and facets according to refference taxa. 
plot<-coef_plot(beta, parameter=c("Organic_content"),title="my sed coef plot",grouping=tax)

plo<-c()
#Estimating all coefficients:
for(i in unique(beta$variable)){
 plo[[i]]<-coef_plot(beta, parameter=c(i),title="COI Sediment Order coefficient plot",grouping=tax)
  
}

plo[4]
