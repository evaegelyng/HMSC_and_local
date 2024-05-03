library(tidyverse)
library(patchwork)
library(fields)
library(Hmsc)


model<- readRDS("results/Water/18S/18S_rich_p50prevclass_water_plus_c2_thin10_s6000_tr02_untrans.rds")


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
    names_to="Class",
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
      beta[[clade]]<- tax[match(beta$Class,grouping[,1]),clade]
    } 

      clade=colnames(grouping)[length(colnames(grouping))]
      
      beta<- beta[order(beta$Division,beta$Phylum,beta$Class),]
      
      beta$Class<- factor(beta$Class, levels=unique(beta$Class))
      
      for(i in 2:length(colnames(grouping))){
        clade=colnames(grouping)[i]
        levels <- unique(beta[[clade]])
        beta[[clade]]<- factor(beta[[clade]], levels=levels)
      } 
      }
     beta$Betapar=round(as.numeric(beta$Betapar),2)

  
  plo<-beta %>% 
    filter(variable != "(Intercept)")%>%
    ggplot(aes(x=Class, y=Betapar))+
    scale_color_manual(values=col_vector[11:23])+
    theme_bw()+
    scale_y_continuous(labels = scales::label_number(accuracy = 0.01))  +
    #scale_y_continuous(labels = scales::label_number(accuracy = 0.001))  + #For quadratic salinity effect
    theme(panel.spacing = unit(0, "lines"), 
          strip.background = element_blank(), 
          strip.placement = "outside",
          axis.text.x = element_text(angle = 90, hjust = 1, size=8, vjust=0),
          strip.text.x = element_text(angle = 90, size = 8),
          legend.position="top",
          legend.box = "horizontal",
          panel.grid.major.y = element_blank(),
          legend.text = element_text(size=10)) +
    xlab("Division") +
    geom_point(size=3,alpha=0.9)+
    # add in a dotted line at zero
    geom_hline(yintercept = 0) +
    labs(
      y ="Estimated effect",
      title = paste(parameter,title, sep=" "))
  

  if(!is.null(grouping)){
  plo<- plo+ 
      facet_grid(~ Class+Phylum+Division, 
                 scales = "free_x", # Let the x axis vary across facets.
                 space = "free_x",  # Let the width of facets vary and force all bars to have the same width.
                 switch = "x")}    # }
  
  return(plo)
  
}



'
Specify taxonomic reference database:
'
tax<-read.table("data/f_tax_ncbi_5_01_22.txt", sep='\t', header=T, comment="")
tax<- tax[,c("Class","Phylum","Division")]

#This function estimates beta coefficients and formats data to plot with ggplot
beta<-coef_beta(model)

####OBS de forskellige parametre har forskellige units!! hvor højt du vil afrunde kan ændres i round(), 
#og hvor mange descimaler i plottet på y aksen kan ændres under scale y continous.  
#plotting function that filters after specified parameter and facets according to reference taxa. 

#Check variable names:
unique(beta$variable)

sal1<-coef_plot(beta, parameter="poly(Salinity, degree = 2, raw = TRUE)1",title="Water",grouping=tax)
Si<-coef_plot(beta, parameter="Si",title="Water",grouping=tax)
PO4<-coef_plot(beta, parameter="PO4",title="Water",grouping=tax)
habitatrocks<-coef_plot(beta, parameter="habitatrocks",title="Water",grouping=tax)
habitateelgrass<-coef_plot(beta, parameter="habitateelgrass",title="Water",grouping=tax)
DN<-coef_plot(beta, parameter="DN",title="Water",grouping=tax)
Temperature<-coef_plot(beta, parameter="Temperature",title="Water",grouping=tax)
Oxygen.depletion<-coef_plot(beta, parameter="Oxygen.depletion",title="Water",grouping=tax)
Chlorophyll<-coef_plot(beta, parameter="Chlorophyll",title="Water",grouping=tax)

ggsave("results/Water/18S/18S_coeff_wat_sal1.png", sal1, width=10, height=9)
ggsave("results/Water/18S/18S_coeff_wat_Si.png", Si, width=10, height=9)
ggsave("results/Water/18S/18S_coeff_wat_PO4.png", PO4, width=10, height=9)
ggsave("results/Water/18S/18S_coeff_wat_rocks.png", habitatrocks, width=10, height=9)
ggsave("results/Water/18S/18S_coeff_wat_eel.png", habitateelgrass, width=10, height=9)
ggsave("results/Water/18S/18S_coeff_wat_DN.png", DN, width=10, height=9)
ggsave("results/Water/18S/18S_coeff_wat_temp.png", Temperature, width=10, height=9)
ggsave("results/Water/18S/18S_coeff_wat_oxy.png", Oxygen.depletion, width=10, height=9)
ggsave("results/Water/18S/18S_coeff_wat_chl.png", Chlorophyll, width=10, height=9)


#sal2<-coef_plot(beta, parameter="poly(Salinity, degree = 2, raw = TRUE)2",title="Water",grouping=tax)
#ggsave("results/Water/18S/18S_coeff_wat_sal2.png", sal2, width=10, height=9)
