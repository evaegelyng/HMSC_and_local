# First, create conda environment with the needed packages:
# Note you need python<3.12 for tensorflow compatibility
#   conda create -n hmsc-hpc bioconductor-phyloseq r-dplyr r-devtools python=3.12
# Activate environment: 
#   conda activate hmsc-hpc 
# Then install hmsc from pip (python package)
#   pip install git+https://github.com/hmsc-r/hmsc-hpc.git@main
# Install tensorflow and cuda
#   pip install tensorflow[and-cuda]
# open R and install Hmsc package from GitHub
#   library(devtools)
#   install_github("hmsc-r/HMSC")

# Load packages
library("phyloseq")
library("Hmsc")
library("dplyr")
library("ape")
library("jsonify")

#Load rarefied dataset
COSQ_rare<-readRDS("data/18S_no_c2_3reps.rds")

# Subset to sediment substrate
COSQ_s<-subset_samples(COSQ_rare, substrate_type=="sediment")

# Load env data
pc_bs<-read.table("data/merged_metadata_230427.txt", sep="\t", header=T, row.names=1)

# Load distance matrix
distsea<-read.table("data/dist_by_sea.txt", sep="\t", header=T)

# Merge at class level 
DT1.2<-tax_glom(COSQ_s, taxrank="class")

# Extracting taxonomy table from phyloseq object
taxa<-data.frame(tax_table(DT1.2))

# Preparing richness calculations
otuo<-data.frame(otu_table(COSQ_s))
taxo<-data.frame(tax_table(COSQ_s), stringsAsFactors=FALSE)
colnames(otuo)<-taxo$class
do<-data.frame(sample_data(COSQ_s))

#Categorisation into taxonomic level. Richness is the sum of unique OTU's within each class in a specific sample.
clades<-levels(factor(taxo$class))

# Prepare for calculating richness of each class in each sample
otuo2<-t(data.frame(otuo, check.names=F))
tabr<-cbind(taxo, otuo2)
# Check tax label columns and remove unnecessary ones
tabr[1:2,1:7]
tabr<-within(tabr,rm("phylum"))
ch<-do$sshc
z<-expand.grid("sshc"=ch, stringsAsFactors=FALSE)

#Sum columns of all samples within the specified taxonomic level acquiring OTU richness.  
for(i in 1:length(clades))
{
gtr<-subset(tabr, class==clades[i])
rich<-colSums(gtr[,-1] != 0)
z<-cbind(z, rich)
t<-1+i
colnames(z)[t]<-clades[i]
}

#This is a new dataframe where we remove the sample column (sample names saved in row names). 
rich_asv<-z[,-1]
dataYrich<-rich_asv

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
pc_bs2<-pc_bs[,c("Salinity","cube_d14N_15N","log_Organic_content","habitat","season","sshc","Grain_size","log_TP","log_N")]

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

pc_bs2<-na.omit(pc_bs2)

pc_bs2$sch<-rownames(pc_bs2)
Time_d<-data.frame(pc_bs[,c("season","Time_d")])
Time_d$season<-as.numeric(Time_d$season)

##Forcing change in time in order to make unique season_date combinations per cluster
Time_d[row.names(Time_d)=="autumn_10_rocks",2]<-131
Time_d[row.names(Time_d)=="autumn_3_rocks",2]<-145
Time_d[row.names(Time_d)=="spring_8_rocks",2]<-5

#Include only samples with non-missing data
Time_d <- Time_d[good_samples, ]
Time_d2<-as.matrix(unique(Time_d))
btr<-pc_bs
btr<-btr[good_samples, ]

#Making Time data rownames match the rest of the data.
rownames(Time_d2)<-unique(btr$sc)

######    Preparing alternative environment data for analysis with substrate type. Popping everything into one data frame
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
dare$cube_d14N_15N<-pc_bs2$cube_d14N_15N[match(dare$sch, pc_bs2$sch)]
dare$log_Organic_content<-pc_bs2$log_Organic_content[match(dare$sch, pc_bs2$sch)]
dare$Grain_size<- pc_bs2$Grain_size[match(dare$sch, pc_bs2$sch)]
dare$log_TP<- pc_bs2$log_TP[match(dare$sch, pc_bs2$sch)]
dare$log_N<-pc_bs2$log_N[match(dare$sch, pc_bs2$sch)]

#Variables to include in the final dataframe
dare<-dare[,c("sch","cluster","season","sc","sshc","habitat","Salinity","cube_d14N_15N","log_Organic_content","Grain_size","log_TP","log_N")]

#Check for "NA"s in class names
dataYrich[, grepl("NA", names(dataYrich))]
# No NAs in class names

write.table(dare, file = "tmp/18S_sed_env.tsv", sep = "\t", row.names = FALSE)
write.table(dataYrich, file="tmp/18S_sed_rich.tsv", sep= "\t")

#Check richness values
A = colSums(dataYrich)
A

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

#Setting priors to a default of 5. 
rl1 = setPriors(rl1,nfMax=5)
rl2 = setPriors(rl2,nfMax=5)

## Create a presence/absence matrix of classes
dataYrich = as.matrix(dataYrich)
dataYrich.pa = 1*(dataYrich>0)
head(dataYrich.pa)

## Create a richness matrix dependent on presence
dataYrich.p = dataYrich
dataYrich.p[dataYrich.p==0] = NA

# Check for skewness of data
jpeg("results/sediment/18S/hist_18S_sed.jpg", width = 350, height = "350")
hist(dataYrich.p)
dev.off()

# Log-transform richness data
dataYrich.p = log(dataYrich.p)

# Check for skewness of data
jpeg("results/sediment/18S/hist_18S_sed_log.jpg", width = 350, height = "350")
hist(dataYrich.p)
dev.off()

Y = cbind(dataYrich.pa,dataYrich.p)
sps = colnames(dataYrich.pa)
ns = length(sps)
sps.subset = colnames(dataYrich.p)
ns.subset = length(sps.subset)
colnames(Y) = c(paste0(sps,"_DNA.PA"),paste0(sps.subset,"_DNA.P"))
colnames(Y) <- noquote(colnames(Y))
head(Y)

# Add taxonomic tree
## Load table with corrected phylum names and marine/non-marine
tax_cur<-read.table("data/18S_classified_phy_class_curated.tsv",sep="\t",header=T)
## Add curated phylum names
taxonomy <- as.data.frame(colnames(Y))
names(taxonomy) <- "Y_names"
taxonomy$class <- sub("_[^_]+$", "", taxonomy$Y_names)
taxonomy$new_phylum<-tax_cur$new_phylum[match(taxonomy$class,tax_cur$class)]

# Load table with supergroups and uni-/multicellular
sgroups <- read.table("data/Supergroups_and_cellularity.tsv", sep='\t', header=T, comment="")
## Add supergroup
taxonomy$supergroup<-sgroups$supergroup[match(taxonomy$new_phylum,sgroups$new_phylum)]

## Define taxonomic hierarchy
frm <- ~supergroup/new_phylum/Y_names
taxonomy$supergroup <- as.factor(taxonomy$supergroup)
taxonomy$new_phylum <- as.factor(taxonomy$new_phylum)
taxonomy$Y_names <- as.factor(taxonomy$Y_names)
tr <- as.phylo(frm, data = taxonomy, collapse=FALSE)
tr$edge.length <- rep(1, nrow(tr$edge))
#plot(tr, show.node.label=TRUE)
Nnode(tr)

# Add a table of data type (p/a and presence-dependent richness)
TrData = data.frame(datatype=as.factor(c(rep("DNA_pa",ns),rep("DNA_p",ns))))
rownames(TrData) = colnames(Y)
head(TrData)

#Create formulas
XFormula = ~ poly(Salinity, degree = 2, raw = TRUE) + cube_d14N_15N + log_N + log_Organic_content + habitat + Grain_size + log_TP

TrFormula = ~ datatype

#Create vector of distributions, setting probit for the p/a data and lognormal poisson for the richness data
my.distr = c(rep("probit",ns),rep("lognormal poisson",ns.subset))

#Set models
#Lognormal poisson allows more flexible poisson modelling that doesnt need standard deviation being equal to the mean.
m = Hmsc(Y=Y, phyloTree = tr, XData=dare, XFormula=XFormula, 
TrData = TrData, TrFormula = TrFormula,
#studyDesign=studyDesign, ranLevels=list("Time_d"=rl1,"space"=rl2), distr=my.distr)
studyDesign=studyDesign, ranLevels=list("Time_d"=rl1), distr=my.distr)

# SAVING MODEL
models = list(m)
names(models) = c("Sediment 18S")
save(models, file = "results/sediment/18S/unfitted_models.RData")

# TESTING THAT MODELS FIT WITHOUT ERRORS
for(i in 1:length(models)){
  print(i)
  sampleMcmc(models[[i]],samples=2)
}

#Run
nSamples = 250
thin = 100
nChains = 4
verbose = 100
transient = ceiling(0.5*nSamples*thin)

# export the model object, using the setting ´engine="HPC"´, which denotes that we are not interested in sampling the model, but only to initialize the sampling.
init_obj = sampleMcmc(m, samples=nSamples, thin=thin,
                                            transient=transient, nChains=nChains,
                                            verbose=verbose, engine="HPC")

init_file_path = "results/sediment/18S/init_file.rds"
saveRDS(to_json(init_obj), file=init_file_path)

# Now run from terminal:
# conda activate hmsc-hpc
# sbatch scripts/hmsc_18S_sed.sh