library(tidyverse)
library(patchwork)
library(fields)
library(Hmsc)
library(phyloseq)

model<- readRDS("results/sediment/COI/COI_sed_241007.rds")

#YOU NEED TO BE WARY OF THE FOLLOWING:
#1. The errorbars need to be added in the same order as species plotting order. Otherwise you get the wrong errorbars for each species.
#2 If you set the coefficient to 0 due to low support level you also need to do it for the error bars.

coef_beta<- function(model){
  
  postBeta = getPostEstimate(model, parName="Beta")
  rownames(as.data.frame(postBeta))

  #credibility intervals are made directly from the posterior distribution (Beta posterior in this case. The beta is made from the link function of the poisson.)
  cred_intervals = getPostEstimate(model,q =c(0.025,0.975), parName="Beta")

  #Intervals are nested in the postBeta under q in a 3 dimensional matrix.
  #if you want to access the intervals of each variable, it is done by entering the variable number in the second dimension, as demonstrated below.
  #cred_intervals$q[,10,] #Here we get every credible interval for the 10th variable

  # Define the minimum accepted support level for a coefficient
  supportLevel=0.95
  mbeta<- postBeta$mean
  betaP=postBeta$support
  toPlot = mbeta

  # Extract only coefficients with support above 0.95 or below 0.05 (why the latter?)
  toPlot = toPlot * ((betaP>supportLevel) + (betaP<(1-supportLevel))>0)

  # Change from matrix to dataframe
  toPlot=data.frame(toPlot)

  # Put env variables as rownames and as a separate column
  rownames(toPlot)<- colnames(model$XScaled)

  #add the credibility intervals to the beta dataframe:
  lower_cred <- cred_intervals$q[1,,]* ((betaP>supportLevel) + (betaP<(1-supportLevel))>0)
  upper_cred <- cred_intervals$q[2,,]* ((betaP>supportLevel) + (betaP<(1-supportLevel))>0)
  
  colnames(lower_cred)<- colnames(toPlot)
  colnames(upper_cred)<- colnames(toPlot)
  rownames(lower_cred)<- paste(rownames(toPlot),"lower_cred",sep="_")
  rownames(upper_cred)<- paste(rownames(toPlot),"upper_cred",sep="_")
  

  toPlot <- rbind(toPlot, lower_cred,upper_cred)
  toPlot$variable<- rownames(toPlot)
  length(colnames(model$XScaled))

  # Count number of orders (response variables)
  len<-length(toPlot)-1
  
  # Change format from list to table
  coef_plot<- toPlot %>% pivot_longer(
    cols=1:len,
    values_to="Betapar",
    names_to="order",
    names_repair = "minimal"
  )
  
  # Change values to numeric format
  coef_plot$Betapar<- as.numeric(coef_plot$Betapar) 
 
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
      beta[[clade]]<- tax[match(beta$order,grouping[,1]),clade]
    } 

      clade=colnames(grouping)[length(colnames(grouping))]
      
      beta<- beta[order(beta$supergroup,beta$phylum,beta$class),]
      
      beta$order<- factor(beta$order, levels=unique(beta$order))
      
      for(i in 2:length(colnames(grouping))){
        clade=colnames(grouping)[i]
        levels <- unique(beta[[clade]])
        beta[[clade]]<- factor(beta[[clade]], levels=levels)
      } 
      }
  
  beta$Betapar=round(as.numeric(beta$Betapar),5)

  #separate the credibility intervals in the original dataframe
  lower_cred<- beta[grepl("lower_cred",beta$variable),]
  upper_cred <-beta[grepl("upper_cred",beta$variable),]
  
  
  beta <-  beta[!grepl("lower_cred",beta$variable),]
  beta <-  beta[!grepl("upper_cred",beta$variable),]
  
  plo<-beta %>% 
    filter(variable != "(Intercept)")%>%
    ggplot(aes(x=order, y=Betapar))+
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
    xlab("") +
    # add in a dotted line at zero
    geom_hline(yintercept = 0) +
    labs(
      y ="Estimated effect",
 title =NULL #paste(parameter,title, sep=" ")
      )+
    geom_errorbar(aes(ymin =lower_cred$Betapar, ymax = upper_cred$Betapar), width = 0.5, color = "blue")

  

  if(!is.null(grouping)){
  plo<- plo+ 
      facet_grid(~ class+phylum+supergroup, 
                 scales = "free_x", # Let the x axis vary across facets.
                 space = "free_x",  # Let the width of facets vary and force all bars to have the same width.
                 switch = "x")}    # }
  
  return(plo)
  
}


'
Specify taxonomic reference database:
'

#COSQ_rare2<-readRDS("data/COI_no_c2_3reps.rds")
#TAX_S = tax_table(COSQ_rare2)
#tax <- data.frame(TAX_S,rownames=F)
#tax<- tax[,c("order","class","phylum","supergroup")]

## Get curated taxonomy
tax<-read.table("data/curated_taxonomy.tsv",sep="\t",header=T)

## Rename phylum column to match functions above
names(tax)[names(tax) == "new_phylum"] <- "phylum"

#This function estimates beta coefficients and formats data to plot with ggplot
beta<-coef_beta(model)

#Check variable names:
unique(beta$variable)

# Plot coefficients including credible intervals
sal1<-coef_plot(beta, parameter="poly(Salinity, degree = 2, raw = TRUE)1",title="Sediment",grouping=tax)
d14N_15N<-coef_plot(beta, parameter="d14N_15N",title="Sediment",grouping=tax)
Organic_content<-coef_plot(beta, parameter="Organic_content",title="Sediment",grouping=tax)
Grain_size<-coef_plot(beta, parameter="Grain_size",title="Sediment",grouping=tax)
TP<-coef_plot(beta, parameter="TP",title="Sediment",grouping=tax)
sal2<-coef_plot(beta, parameter="poly(Salinity, degree = 2, raw = TRUE)2",title="Sediment",grouping=tax)
N<-coef_plot(beta, parameter="N",title="Sediment",grouping=tax)

ggsave("results/sediment/COI/COI_coeff_sed_sal1.png", sal1, width=10, height=9)
ggsave("results/sediment/COI/COI_coeff_sed_d14N_15N.png", d14N_15N, width=10, height=9)
ggsave("results/sediment/COI/COI_coeff_sed_Organic.png", Organic_content, width=10, height=9)
ggsave("results/sediment/COI/COI_coeff_sed_grain.png", Grain_size, width=10, height=9)
ggsave("results/sediment/COI/COI_coeff_sed_TP.png", TP, width=10, height=9)
ggsave("results/sediment/COI/COI_coeff_sed_sal2.png", sal2, width=10, height=9)
ggsave("results/sediment/COI/COI_coeff_sed_N.png", N, width=10, height=9)

'
The coefficient plot for habitat type

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
  names_to="order",
  names_repair = "minimal"
)

coef_plot$phylum<-tax$phylum[match(coef_plot$order,tax$order)] 
coef_plot$supergroup<-tax$supergroup[match(coef_plot$order,tax$order)] 
coef_plot$class<-tax$class[match(coef_plot$order,tax$order)] 
coef_plot<- coef_plot[order(coef_plot$supergroup,coef_plot$phylum,coef_plot$class),]
coef_plot$order<- factor(coef_plot$order, levels=unique(coef_plot$order))
coef_plot$Betapar<- as.numeric(coef_plot$Betapar) 

require(RColorBrewer)
#n <- 9
qual_col_pals = brewer.pal.info[brewer.pal.info$category =="qual",]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


head(coef_plot)

# Plot for habitat type
hab<- coef_plot %>% group_by(variable) %>% filter(str_starts(variable,"habitat",negate=F))
hab<- hab[order(hab$supergroup,hab$phylum,hab$class,hab$order),]

p<-hab%>% 
  filter(variable != "(Intercept)")%>%
  ggplot(aes(x=order, y=Betapar, color=variable)) +
  facet_grid(~class+phylum+supergroup, 
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

ggsave(p,file="results/sediment/COI/coef_sed_hab.png",height=9,width=10)