#Coef plot function

library(tidyverse)
library(patchwork)
library(fields)
library(Hmsc)
library("phyloseq")


model<-  readRDS("C:/Users/zoega/OneDrive/Skrivebord/Kandidat_3/projekt/results/COI_rich_sed_FIVAL.rds") # Just to read in the model

coef_beta<- function(model){
  
  postBeta = getPostEstimate(model, parName="Beta")
  rownames(as.data.frame(postBeta))
  cred_intervals = getPostEstimate(model,q =c(0.025,0.975), parName="Beta")
  
  supportLevel=0.95
  mbeta<- postBeta$mean
  betaP=postBeta$support
  toPlot = mbeta

  toPlot = toPlot * ((betaP>supportLevel) + (betaP<(1-supportLevel))>0)
  toPlot=data.frame(toPlot)
  rownames(toPlot)<- colnames(model$XScaled)
  lower_cred <- cred_intervals$q[1,,]* ((betaP>supportLevel) + (betaP<(1-supportLevel))>0)
  upper_cred <- cred_intervals$q[2,,]* ((betaP>supportLevel) + (betaP<(1-supportLevel))>0)
  
  colnames(lower_cred)<- colnames(toPlot)
  colnames(upper_cred)<- colnames(toPlot)
  rownames(lower_cred)<- paste(rownames(toPlot),"lower_cred",sep="_")
  rownames(upper_cred)<- paste(rownames(toPlot),"upper_cred",sep="_")
  
  
  
  toPlot <- rbind(toPlot, lower_cred,upper_cred)
  toPlot$variable<- rownames(toPlot)
  length(colnames(model$XScaled))

  
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
  
  
  beta<-beta %>% data.frame() %>% filter(grepl(parameter,variable))
  
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
   
  beta$Betapar=round(as.numeric(beta$Betapar),5)
  lower_cred<- beta[grepl("lower_cred",beta$variable),]
  upper_cred <-beta[grepl("upper_cred",beta$variable),]
  
  
  beta <-  beta[!grepl("lower_cred",beta$variable),]
  beta <-  beta[!grepl("upper_cred",beta$variable),]

  
  plo<-beta %>% 
    filter(variable != "(Intercept)")%>%
    ggplot(aes(x=Species, y=Betapar))+
    geom_point(size=3,alpha=0.9,aes(color=variable))+
    scale_color_manual(values=c("#B3DE69","#80B1D3"))+
    theme_bw()+
    #scale_y_continuous(limits=c(min(beta$Betapar),max(beta$Betapar)),labels = scales::label_number(accuracy = 0.001))  + 
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
    # add in a dotted line at zero
    geom_hline(yintercept = 0) +
    labs(
      y ="Estimate of effect",
      title =NULL #paste(parameter,title, sep=" ")
      #title="Salinity COI Sediment Order coefficient plot"
      )+
    geom_errorbar(aes(ymin =lower_cred$Betapar, ymax = upper_cred$Betapar), width = 0.5, color = "blue")
  
  if(!is.null(grouping)){
   
    
  plo<- plo+ 
      facet_grid(~ class+phylum+kingdom, 
                 scales = "free_x", # Let the x axis vary across facets.
                 space = "free_x",  # Let the width of facets vary and force all bars to have the same width.
                 switch = "x")}    # }
  
  return(plo)
  
}

COSQ_rare2<-readRDS("C:/Users/zoega/OneDrive/Skrivebord/Kandidat_3/projekt/Data/COSQ_rare2_correct_score.rds")
OTU_COI = otu_table(COSQ_rare2, taxa_are_rows = TRUE)
TAX_S = tax_table(COSQ_rare2)
tax <- data.frame(TAX_S,rownames=F)
tax<- tax[,c("order","class","phylum","kingdom")]

beta<-coef_beta(model)


plot<-coef_plot(beta, parameter=c("Salinity"),title=NULL,grouping=tax)



#As adviced under Dereks Sargacious mathmatical mind
#credibility intervals of a distribution is made directly from the posterior distribution (Beta posterior in this case. THe beta is made from the link function of the poisson. I know very confusing!).

#This function is the same that we used to estimate the mean(coef value). It can also give the credibility intervals
# by adding the "q" argument.
postBeta = getPostEstimate(model,q =c(0.025,0.975), parName="Beta")
View(postBeta$q)

#We now have the intervals but we need to plot it!!
#Intervals are nested in the postBeta under q in a 3 dimensional matrix.
#if you want to acess the intervals of each variable variable it is done by entering the variable number in the second dimension, as demonstrated below.
postBeta$q[,10,] #Here we get every credible interval for the 10th variable

#So how do we plot this? Unfortunately i do not have access to the server so you will have to do the bulk yourself :)
#i will try an demonstrate how below just add it to the scripts you have on the server.  

unique(beta$variable) # <- this is the number of the variables

#YOU NEED TO BE WARY OF THE FOLLOWING:
#1. The errorbars need to be added in the same order as species plotting order. Otherwise you get the wrong errorbars for each species.
#2 If you set the coefficient to 0 due to low support level you also need to do it for the error bars.

###########################################
#####Implementation intstruction###########
###########################################

#add the credibility intervals to the beta dataframe (you can allways check how i did above in the coef beta):
lower_cred <- cred_intervals$q[1,,]* ((betaP>supportLevel) + (betaP<(1-supportLevel))>0)
upper_cred <- cred_intervals$q[2,,]* ((betaP>supportLevel) + (betaP<(1-supportLevel))>0)

colnames(lower_cred)<- colnames(toPlot)
colnames(upper_cred)<- colnames(toPlot)
rownames(lower_cred)<- paste(rownames(toPlot),"lower_cred",sep="_")
rownames(upper_cred)<- paste(rownames(toPlot),"upper_cred",sep="_")

toPlot <- rbind(toPlot, lower_cred,upper_cred)

############################################
###TO add the credibility intervals#########

#change this: beta<-beta %>% data.frame() %>% filter(variable in parameter)
# to this: beta<-beta %>% data.frame() %>% filter(grepl(parameter,variable))

#seperate the credibility intervals in the original dataframe (check how i did it in coef_plot:
lower_cred<- beta[grepl("lower_cred",beta$variable),]
upper_cred <-beta[grepl("upper_cred",beta$variable),]

beta <-  beta[!grepl("lower_cred",beta$variable),]
beta <-  beta[!grepl("upper_cred",beta$variable),]

#add the errorbar:
geom_errorbar(aes(ymin =lower_cred$Betapar, ymax = upper_cred$Betapar), width = 0.5, color = "blue")





                                     