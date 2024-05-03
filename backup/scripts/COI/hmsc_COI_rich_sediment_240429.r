# First, create conda environment with the needed packages:
# conda create -n hmsc bioconductor-phyloseq r-ggplot2 r-hmsc r-vegan r-reshape2 r-plyr r-scales r-stringr r-RColorBrewer r-corrplot
# Activate environment: 
#conda activate hmsc 

# Load packages
library("phyloseq")
library("Hmsc")
library("dplyr")
library("tidyverse")

# Load rarefied 70%-similarity dataset (MOTU dataset)
COSQ_rare<-readRDS("data/COI_no_c2_3reps.rds")

#Subset to sediment substrate
COSQ_s<-subset_samples(COSQ_rare, substrate_type=="sediment")

#Load env data
pc_bs<-read.table("data/merged_metadata_230427.txt", sep="\t", header=T, row.names=1)
fsed<-read.table("data/sed_metadata.txt", sep="\t", header=T)
spat<-read.table("data/Spatial_values_240430.csv", sep=",", header=T, row.names=1)
pc_bs$TP<-fsed$TP[match(row.names(pc_bs),fsed$snch)]
pc_bs$d14N_15N<-fsed$d14N_15N[match(row.names(pc_bs),fsed$snch)]
pc_bs$Oxygen.depletion<-spat$Oxygen.depletion[match(row.names(pc_bs),row.names(spat))]
pc_bs$FishingTrawling<-spat$FishingTrawling[match(row.names(pc_bs),row.names(spat))]

# Load bacterial richness data
bac<-read.table("results/bact_rich.tsv", sep="\t", header=T, row.names=1)
pc_bs$sshc<-paste("sediment",pc_bs$season,pc_bs$habitat,pc_bs$cluster,sep="_")
pc_bs$bac_rich<-bac$rich[match(pc_bs$sshc,rownames(bac))]

#Load distance matrix
distsea<-read.table("data/dist_by_sea.txt", sep="\t", header=T)

#merge at order level.
DT1.2<-tax_glom(COSQ_s, taxrank="order")

#creating abundance dataframe from phyloseq object
dataYabund<-data.frame(otu_table(DT1.2))
taxa<-data.frame(tax_table(DT1.2))
colnames(dataYabund)<-taxa$order

# Preparing richness calculations
otuo<-data.frame(otu_table(COSQ_s))
taxo<-data.frame(tax_table(COSQ_s), stringsAsFactors=FALSE)
colnames(otuo)<-taxo$order
do<-data.frame(sample_data(COSQ_s))

#Categorisation into taxonomic level. Richness is the sum of unique OTU's within each order in a specific sample.
clades<-levels(factor(taxo$order))

# Prepare for calculating richness of each order in each sample
otuo2<-t(data.frame(otuo, check.names=F))
tabr<-cbind(taxo, otuo2)
# Check tax label columns and remove unnecessary ones
tabr[1:2,1:9]
tabr<-within(tabr,rm("kingdom","class","phylum","family","genus","species","score.id.pident.70"))
ch<-do$sshc
z<-expand.grid("sshc"=ch, stringsAsFactors=FALSE)

#Filtering out all "order scores" (proportion of best BLAST matches that belonged to the assigned order) beneath 85 percent.  
tabr<- tabr %>% filter(as.numeric(order_score)>85)
tabr<-subset(tabr,select=-order_score)


#Sum columns of all samples within the specified taxonomic level acquiring OTU richness.  
for(i in 1:length(clades))
{
gtr<-subset(tabr, order==clades[i])
rich<-colSums(gtr[,-1] != 0)
z<-cbind(z, rich)
t<-1+i
colnames(z)[t]<-clades[i]
}

#This is a new dataframe where we remove the habitat column. 
rich_asv<-z[,-1]
dataYrich<-rich_asv

# Calculating relative read abundance
#Transforms sample counts of a taxa abundance matrix by sum. 
datarg = transform_sample_counts(DT1.2, function(x) x/sum(x))
tax<-data.frame(tax_table(datarg), stringsAsFactors=FALSE)
otu<-otu_table(datarg)
otu<-t(data.frame(otu,check.names=F))
tab<-cbind(tax, otu)

#Checking for NA values in the order column
tab[tab$order == "NA",]
#removing any non classified taxa
tab<- tab[tab$order != "NA",]
                                 
# Check tax label columns and remove unnecessary ones
tab[1:2,1:9]
tab<-within(tab,rm("kingdom","phylum","class","order","family","genus","species","score.id.pident.70"))
ttab<-t(data.frame(tab, check.names=F))
class(ttab) <- "numeric"


###################################################

#Make dist data rectangle
#The distance between each coordinate measured using euclidean distance
ds<-distsea[with(distsea, order(SiteA, SiteB)), ]
ds2 <- with(ds, Mads_calc) 
                                 
#takes all the values of the variables site A and B,  Concatenates them.
#ss is a vector with all the values in site A and B
ss <- with(ds, unique(c(as.character(SiteA), as.character(SiteB))))

#This assigns attributes                                
attributes(ds2) <- with(ds, list(Size = length(ss),
                                  Labels = ss,
                                  Diag = FALSE,
                                  Upper = FALSE,
                                  method = "user"))
class(ds2) <- "dist"

#Format date as days from first sampling day
dsds<-sort(as.Date(pc_bs$Time, format="%d-%m-%y"))
pc_bs$Time_d<-as.integer(abs(as.Date(pc_bs$Time, format="%d-%m-%y")-dsds[1]))

#Prepare environment data
#For HMSC it is required that categorical variables are factors.
pc_bs$season <- factor(pc_bs$season, levels = c("spring","autumn"))
pc_bs$habitat <- factor(pc_bs$habitat, levels = c("sand","rocks","eelgrass"))
pc_bs$sc<-paste(pc_bs$season, pc_bs$cluster, sep="_")
pc_bs$sch<-rownames(pc_bs)
pc_bs$sshc<-paste("sediment",pc_bs$season,pc_bs$habitat, pc_bs$cluster, sep="_")

head(pc_bs)
#What variables should be kept?
pc_bs2<-pc_bs[,c("Oxygen.depletion","Salinity","d14N_15N","Organic_content","habitat","season","sshc","Grain_size","TP","N","bac_rich","FishingTrawling")]

# Identify samples with complete metadata
print("Total number of samples")
nrow(pc_bs2)
good_samples<-rownames(pc_bs2)[rowSums(is.na(pc_bs2)) == 0]
print("Number of samples with complete metadata")
length(good_samples)
pc_bs3<-pc_bs2
rownames(pc_bs3)<-pc_bs3$sshc
good_samples_sp<-rownames(pc_bs3)[rowSums(is.na(pc_bs3)) == 0]

# Subset data to samples with complete metadata
dataYrich <- dataYrich[rownames(dataYrich) %in% good_samples_sp, ]
dataYabund <- dataYabund[rownames(dataYabund) %in% good_samples_sp, ]

pc_bs2<-na.omit(pc_bs2)
pc_bs2$sch<-rownames(pc_bs2)
Time_d<-data.frame(pc_bs[,c("season","Time_d")])
Time_d$season<-as.numeric(Time_d$season)

#Forcing change in time in order to make unique season_date combinations per cluster
Time_d[row.names(Time_d)=="autumn_10_rocks",2]<-131
Time_d[row.names(Time_d)=="autumn_3_rocks",2]<-145
Time_d[row.names(Time_d)=="spring_8_rocks",2]<-5

	
#Include only samples with non missing data
Time_d <- Time_d[good_samples, ]
Time_d2<-as.matrix(unique(Time_d))
btr<-pc_bs
btr<-btr[good_samples, ]

#Making Time data rownames match the rest of the data.
rownames(Time_d2)<-unique(btr$sc)


######  Preparing alternative environment data for analysis with substrate type. Popping everything into one data frame 
dare<-data.frame(sample_data(DT1.2))
dare <- dare[rownames(dare) %in% good_samples_sp, ]
dare$sch<-paste(dare$season,dare$cluster, dare$habitat, sep="_")
dare$sc<-paste(dare$season,dare$cluster, sep="_")
dare$habitat <- factor(dare$habitat, levels = c("sand","rocks","eelgrass"))
dare$cluster<-as.factor(dare$cluster)
dare$season<-as.factor(dare$season)
dare$sch<-as.factor(dare$sch)
dare$sshc<-as.factor(dare$sshc)
dare$sc<-as.factor(dare$sc)
dare$Salinity<-pc_bs2$Salinity[match(dare$sch, pc_bs2$sch)]
dare$d14N_15N<-pc_bs2$d14N_15N[match(dare$sch, pc_bs2$sch)]
dare$Organic_content<-pc_bs2$Organic_content[match(dare$sch, pc_bs2$sch)]
dare$Oxygen.depletion<- pc_bs2$Oxygen.depletion[match(dare$sch, pc_bs2$sch)]
dare$Grain_size<- pc_bs2$Grain_size[match(dare$sch, pc_bs2$sch)]
dare$TP<- pc_bs2$TP[match(dare$sch, pc_bs2$sch)]
dare$N<-pc_bs2$N[match(dare$sch, pc_bs2$sch)]
dare$bac_rich<-pc_bs2$bac_rich[match(dare$sch, pc_bs2$sch)]
dare$FishingTrawling<-as.factor(pc_bs2$FishingTrawling[match(dare$sch, pc_bs2$sch)])

#Variables to include in the final dataframe
dare<-dare[,c("Oxygen.depletion","sch","cluster","season","sc","sshc","habitat","Salinity","d14N_15N","Organic_content","Grain_size","TP","N","bac_rich","FishingTrawling")]

#Filtering by rich, abund and sites occurring
#I.e. picking out the most rich, abundant and frequent classes of species. 
#Too many rare taxa result in too much noise in the modelling                                  
pabund<-as.data.frame(cbind(log10(colMeans(dataYabund)), colnames(dataYabund)))
colnames(pabund)[2]<-"clades_p"
prich<-as.data.frame(cbind(log10(colMeans(dataYrich)), colnames(dataYrich)))
colnames(prich)[2]<-"clades_p"
trich<-t(dataYrich)
prich$sqrt_sites_occur <-sqrt(rowSums(trich != 0))

#filtering orders based on no. of OTUs, reads and clusters                                 
abund_rich_summary <- as.data.frame(merge(prich, pabund, by="clades_p"))
colnames(abund_rich_summary)<-c("clades_p","log10_mean_rich","sqrt_sites_occur","log10_mean_rel_abund")
abund_rich_summary$log10_mean_rich<-as.numeric(abund_rich_summary$log10_mean_rich)
abund_rich_summary$log10_mean_rel_abund<-as.numeric(abund_rich_summary$log10_mean_rel_abund)
abundsorted<-abund_rich_summary[order( abund_rich_summary[,"log10_mean_rel_abund"] ),]
richsorted<-abund_rich_summary[order( abund_rich_summary[,"log10_mean_rich"] ),]
sitessorted<-abund_rich_summary[order( abund_rich_summary[,"sqrt_sites_occur"] ),]

#Get 50% most prevalent phyla
print("Number of classes before prevalence filtering")
nrow(abund_rich_summary)
rich_p50<-abund_rich_summary[abund_rich_summary$sqrt_sites_occur > quantile(abund_rich_summary$sqrt_sites_occur, 0.5), ]
print("Number of classes after prevalence filtering")
nrow(rich_p50)
tax_to_keep2<-as.character(rich_p50$clades_p)

#Update tables                                 
dataYrich <- dataYrich[, tax_to_keep2]
dataYabund <- dataYabund[, tax_to_keep2]

#Removing undefined clades
dataYrich<- dataYrich[, !grepl("NA", names(dataYrich))]

write.table(dare, file = "tmp/Dare_COI_sed.tsv", sep = "\t", row.names = FALSE)
write.table(dataYrich, file="tmp/Species_data_COI_sed.tsv", sep= "\t")

#Check prevalence and abundance
P = colMeans(dataYrich>0)
A = colSums(dataYrich)

#Histogram of scaled richness densities
pdf("results/sediment/COI/hist_COI_rich_p50prevorder_sediment_plus.pdf")
par(mfrow=c(1,2))
hist(P,xlim=c(0,1),breaks = seq(from=0,to=1,by=0.1), col = "grey", xlab = "Prevalence")
hist(log(A,base=10),breaks=10, col = "grey", xlab = "log10 Richness")
dev.off()

P = colMeans(dataYabund>0)
A = colSums(dataYabund)/sum(dataYabund)

#histogram of scaled abundance densities
pdf("results/sediment/COI/hist_COI_abund_p50prevorder_sediment_plus.pdf")
par(mfrow=c(1,2))
hist(P,xlim=c(0,1),breaks = seq(from=0,to=1,by=0.1), col = "grey", xlab = "Prevalence")
hist(log(A,base=10),breaks=10, col = "grey", xlab = "log10 Abundance")
dev.off()

#Removing the zero inflated classes. 
num<- length(names(dataYrich))

zero.inflated<- data.frame(matrix(NA, nrow = nrow(dataYrich), ncol = 0))
rownames(zero.inflated)<-rownames(dataYrich)

for(i in 1:num){
  min_1=FALSE
  for(j in 2:5){
    col_index=6*i+j-6
    if(summary(dataYrich)[col_index]=="1st Qu.:0.0000  "){
      min_1=TRUE}
    

    if(5+6*(i-1)==col_index && min_1 && summary(dataYrich)[col_index]=="3rd Qu.:0.0000  "){
      
      print("zero inflated")
      
      zero.inflated<- cbind(zero.inflated,dataYrich[i])
      
      }
    
    }
  
}

names(zero.inflated)
dataYrich1<- dataYrich[, !names(dataYrich) %in% names(zero.inflated)]

write.table(dataYrich1, file="tmp/Species_data_COI_filt_sed.tsv", sep= "\t")

###############################
#HMSC study design. 
#Study design needs to include comparable columns to the randomeffectlevel and predictor dataframe(dare).
#See it as intermediate dataframe between the randomeffect levels and the rest of the data.
studyDesign = data.frame(
    season = as.factor(dare$season),
    sch = as.factor(dare$sch),
    habitat = as.factor(dare$habitat),
    cluster = as.factor(dare$cluster),
    Time_d = as.factor(dare$sc),
    space = as.factor(dare$cluster),
    sshc = as.factor(dare$sshc)
)

#Create random effects structure
rl1 = HmscRandomLevel(sData = Time_d2)
rl2 = HmscRandomLevel(distMat = ds2)

#Setting priors to a default 5. 
rl1 = setPriors(rl1,nfMax=5)
rl2 = setPriors(rl2,nfMax=5)

rnd_ef<-list("Time_d"=rl1,"space"=rl2)

#Create formula
XFormula = ~ poly(Salinity, degree = 2, raw = TRUE) + d14N_15N + N + Organic_content + habitat + Grain_size + TP + Oxygen.depletion + bac_rich + FishingTrawling


#Set models; abund or rich, evaluate different distributions
#Lognormal poisson allows more flexible poisson modelling that doesnt need standard deviation being equal to the mean.
m.full.r.lp = Hmsc(Y=dataYrich1, XData=dare, XFormula=XFormula,
studyDesign=studyDesign, ranLevels=rnd_ef, distr="lognormal poisson")


#Run
nChains = 2
nParallel = nChains
thin = 40
samples = 6000
transient = 0.2*thin*samples
verbose = samples/10


models = sampleMcmc(m.full.r.lp, thin = thin,
samples = samples, transient = transient,
nChains = nChains, nParallel = nParallel, verbose = verbose)

#Saving model
saveRDS(models, file = "results/sediment/COI/COI_240426.rds")

#Diagnostic
modelsII<-models

#Widely Applicable information Criteria
WAIC2 = computeWAIC(hM=modelsII, byColumn=FALSE)
preds = computePredictedValues(modelsII)

#To get the sr2 the predicted values are used to estimates amount of explained variance.
MF= evaluateModelFit(hM = modelsII, predY = preds)

#Converting model object to r interpretable coda
mpost = convertToCodaObject(modelsII, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
ess.beta = effectiveSize(mpost$Beta)
gd<-gelman.diag(mpost$Beta, multivariate=T)

print("WAIC2")
WAIC2

#Explanatory power
#For Poisson models, a pseudo-R2 is computed as squared spearman correlation between observed and
# predicted values, times the sign of the correlation (SR2).
print("Mean of SR2")
mean(MF$SR2)

print("Std dev of SR2")
sd(MF$SR2)

#Predictive power
print("Mean of ESS")
mean(ess.beta)

print("Std dev of ESS")
sd(ess.beta)

#Potential scale reduction factor gives an estimate of model parameter convergence. This number should be <1.1
print("Mean of PSRF")
mean(gd$psrf)

print("Std dev of PSRF")
sd(gd$psrf)

#Histogram of effective sample size
pdf("results/sediment/COI/diag_240426.pdf")
par(mfrow=c(1,2))
hist(effectiveSize(mpost$Beta), main="ess(beta)")
hist(gelman.diag(mpost$Beta, multivariate=F)$psrf, main="psrf(beta)")
dev.off()

#Posterior beta estimation over different iterations gives an idea of parameter convergence for each parameter.
pdf("results/sediment/COI/diag_beta_240426.pdf")
plot(mpost$Beta)
dev.off()

head(modelsII$X)



