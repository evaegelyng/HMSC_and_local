library(tidyverse)
library(patchwork)
library(fields)
library(Hmsc)
library(phyloseq)


model<-  readRDS("results/sediment/18S/18S_240426.rds")

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
    #scale_y_continuous(labels = scales::label_number(accuracy = 0.00000001))  + #For bacterial richness and N
    theme(panel.spacing = unit(0, "lines"), 
          strip.background = element_blank(), 
          strip.placement = "outside",
          axis.text.x = element_text(angle = 90, hjust = 1, size=8, vjust=0),
          strip.text.x = element_text(angle = 90, size = 8),
          legend.position="top",
          legend.box = "horizontal",
          panel.grid.major.y = element_blank(),
          legend.text = element_text(size=10)) +
    xlab("") +
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
COSQ_rare2<-readRDS("data/18S_no_c2_3reps.rds")
OTU_18S = otu_table(COSQ_rare2, taxa_are_rows = TRUE)
TAX_S = tax_table(COSQ_rare2)
tax <- data.frame(TAX_S,rownames=F)
tax<- tax[,c("Class","Phylum","Division")]

#This function estimates beta coefficients and formats data to plot with ggplot
beta<-coef_beta(model)

####OBS de forskellige parametre har forskellige units!! hvor højt du vil afrunde kan ændres i round(), 
#og hvor mange descimaler i plottet på y aksen kan ændres under scale y continous.  
#plotting function that filters after specified parameter and facets according to reference taxa. 

plo<-list()
#Estimating all coefficients:
for(i in unique(beta$variable)){
 plo[[i]]<-coef_plot(beta, parameter=c(i),title=paste("Sediment,",i,sep=" "),grouping=tax)
  
}

sal1<-coef_plot(beta, parameter="poly(Salinity, degree = 2, raw = TRUE)1",title="Sediment",grouping=tax)
d14N_15N<-coef_plot(beta, parameter="d14N_15N",title="Sediment",grouping=tax)
Organic_content<-coef_plot(beta, parameter="Organic_content",title="Sediment",grouping=tax)
Grain_size<-coef_plot(beta, parameter="Grain_size",title="Sediment",grouping=tax)
TP<-coef_plot(beta, parameter="TP",title="Sediment",grouping=tax)
Oxygen.depletion<-coef_plot(beta, parameter="Oxygen.depletion",title="Sediment",grouping=tax)
No_fishing<-coef_plot(beta, parameter="FishingTrawlingNo fishing/No ban",title="Sediment",grouping=tax)
Yes_fishing<-coef_plot(beta, parameter="FishingTrawlingYes fishing/No ban",title="Sediment",grouping=tax)

ggsave("results/sediment/18S/18S_coeff_sed_sal1.png", sal1, width=10, height=9)
ggsave("results/sediment/18S/18S_coeff_sed_d14N_15N.png", d14N_15N, width=10, height=9)
ggsave("results/sediment/18S/18S_coeff_sed_Organic.png", Organic_content, width=10, height=9)
ggsave("results/sediment/18S/18S_coeff_sed_grain.png", Grain_size, width=10, height=9)
ggsave("results/sediment/18S/18S_coeff_sed_TP.png", TP, width=10, height=9)
ggsave("results/sediment/18S/18S_coeff_sed_oxy.png", Oxygen.depletion, width=10, height=9)
ggsave("results/sediment/18S/18S_coeff_sed_nofish.png", No_fishing, width=10, height=9)
ggsave("results/sediment/18S/18S_coeff_sed_yesfish.png", Yes_fishing, width=10, height=9)

#sal2<-coef_plot(beta, parameter="poly(Salinity, degree = 2, raw = TRUE)2",title="Sediment",grouping=tax)
#N<-coef_plot(beta, parameter="N",title="Sediment",grouping=tax)
#Bac_rich<-coef_plot(beta, parameter="bac_rich",title="Sediment",grouping=tax)
#ggsave("results/sediment/18S/18S_coeff_sed_sal2.png", sal2, width=10, height=9)
#ggsave("results/sediment/18S/18S_coeff_sed_N.png", N, width=10, height=9)
#ggsave("results/sediment/18S/18S_coeff_sed_bac.png", Bac_rich, width=10, height=9)

'
The coefficient plot for habitat types

'
postBeta = getPostEstimate(model, parName="Beta")
supportLevel=0.95
mbeta<- postBeta$mean
betaP=postBeta$support
toPlot = mbeta
# Support level operations. if the logical expression on the right is False it will set to 0. 
#The beta supportlevel paramater is wether the mean beta is significant.
toPlot = toPlot * ((betaP>supportLevel) + (betaP<(1-supportLevel))>0)

toPlot=data.frame(toPlot)

rownames(toPlot)<- colnames(model$XScaled)

toPlot$variable<- rownames(toPlot) 

length(toPlot)

coef_plot<- toPlot %>% pivot_longer(
  cols=1:(length(toPlot)-1),
  values_to="Betapar",
  names_to="Class",
  names_repair = "minimal"
)

coef_plot$Phylum<-tax$Phylum[match(coef_plot$Class,tax$Class)] 
coef_plot$Division<-tax$Division[match(coef_plot$Class,tax$Class)] 
coef_plot<- coef_plot[order(coef_plot$Division,coef_plot$Phylum),]
coef_plot$Class<- factor(coef_plot$Class, levels=unique(coef_plot$Class))
coef_plot$Betapar<- as.numeric(coef_plot$Betapar) 

require(RColorBrewer)
#n <- 9
qual_col_pals = brewer.pal.info[brewer.pal.info$category =="qual",]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


head(coef_plot)

# Plot for habitat type
hab<- coef_plot %>% group_by(variable) %>% filter(str_starts(variable,"habitat",negate=F))
hab<- hab[order(hab$Division,hab$Phylum,hab$Class),]

p<-hab%>% 
  filter(variable != "(Intercept)")%>%
  ggplot(aes(x=Class, y=Betapar, color=variable)) +
  facet_grid(~Class+Phylum+Division, 
             scales = "free_x", # Let the x axis vary across facets.
             space = "free_x",  # Let the width of facets vary and force all bars to have the same width.
             switch = "x")+
  scale_color_manual(values=c("yellowgreen","cornflowerblue"))+
  theme_bw() + theme(panel.spacing = unit(0, "lines"), strip.background = element_blank(), strip.placement = "outside", axis.text.x = element_text(angle = 90, hjust = 1, size=8, vjust=0), strip.text.x = element_text(angle = 90, size = 8), legend.position="top",legend.box = "horizontal", panel.grid.major.y = element_blank(),legend.text = element_text(size=10)) +
  xlab("") +
  geom_point(size=3,alpha=0.9)+
  # add in a dotted line at zero
  geom_hline(yintercept = 0) +
  labs(
    y ="Estimated effect",
    title = "")

ggsave(p,file="results/sediment/18S/coef_sed_hab.png",height=14,width=12)