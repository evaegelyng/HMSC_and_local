
# http://citeseerx.ist.psu.edu/viewdoc/download;jsessionid=2E55EEDF17FC2211F307B3B703D6AF8E?doi=10.1.1.694.9409&rep=rep1&type=pdf
#Haploid Data
#The explanation given above using heterozygosity for determining FST for diploid loci is all well and
#good, but how do you define FST for a haploid locus, where heterozygosity is relatively meaningless? This
#has been approached in a number of ways (as usual), but the simplest way of reconceptualising FST for
#haploid loci is to think in terms of haplotype diversity instead of heterozygosity. Haplotype diversity
#(conveniently, also referred to as H) is a measure of the degree of variation in haplotypes found within a
#population, and is calculated as:

# H = 1- sum(p²) haplotipic diversity
#Fst = (Ht - mean(Hs)) / Ht
#Ht -> Htotal
#Hs -> Hsubpob

# according to Carles Carreres with the same definition what we do is calculate the Djost which has the same def as the Fst

library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(tidyverse)
library(adegenet)
library(mmod)
library(future.apply)
library(phyloseq)

# error standard
errstd <- function(x) sd(x, na.rm = T)/sqrt(length(x[!is.na(x)]))
sq <- function(x,a) {
  if (x==0) x<- 'A' else
    if (x<=a[1]) x<-'B' else
      if (x<=a[2]) x<-'C' else
        x<-'D'
}
semi_quant <- function(x){
  a <- unname(quantile(mitjanes_Djost_dins_fora$mean_Djost, probs = c(.1,.9)))
  x<- sapply(x, sq, a = a)
}
semi_quantb <- function(x){
  a <- unname(quantile(mitjanes_Djost_dins_fora$total_reads, probs = c(.25,.75)))
  x<- sapply(x, sq, a = a)
}

correr_sig=FALSE


# import data
initial_data<-read.table("../../../Tekstfiler/COI/COI_ASV/sed_ASVs_250512.txt", sep="\t", header=T, check.names = F)
motu<-read.table("../../../Tekstfiler/COI/COI_ASV/COSQ_final_ASV.tsv", sep="\t", header=T, check.names=F)
row.names(motu)<-motu$id

initial_data$motu<-motu$motu[match(row.names(initial_data),row.names(motu))]
initial_data$seq <- "ACGT" # random sequence
initial_data$count <- "5" # random number

initial_data <- data.frame(motu=initial_data$final.id, id=initial_data$id,seq=initial_data$seq, count=initial_data$count, initial_data[,3:(ncol(initial_data)-3)], check.names = F)
initial_data_temp <- initial_data[,5:ncol(initial_data)]

initial_data_temp[initial_data_temp>=0.75] <- 4
initial_data_temp[initial_data_temp>=0.50 & initial_data_temp<0.75] <- 3
initial_data_temp[initial_data_temp>=0.25 & initial_data_temp<0.50] <- 2
initial_data_temp[initial_data_temp>0 & initial_data_temp<0.25] <- 1
initial_data_temp[initial_data_temp==0 ] <- 0

initial_data <- data.frame(initial_data[,1:4], initial_data_temp, check.names = F)

all_data_reads <- c(rep(5, length(unique(initial_data$motu))))
names(all_data_reads) <- unique(initial_data$motu)

# I reworked the metadatafile in excel
metadata<-read.table("../../../Tekstfiler/COI/COI_ASV/Mads_metadata_sed2_250515.txt", sep="\t", header=TRUE) # check if metadata file is correct

n_sites<-ncol(initial_data_temp)

Djost_df_3D <- matrix(NA,nrow = n_sites, ncol = n_sites)
Djost_df_3D_fronts <- Djost_df_3D
Djost_sig_3D <- Djost_df_3D

for (row in 1:n_sites) {
  for (column in 1:n_sites) {
    if (row == column) {
      Djost_df_3D_fronts[row, column] <- 0
    } else if (abs(metadata$Geographic_Order[row] - metadata$Geographic_Order[column]) == 1 & metadata$Sub_Area[row] == metadata$Sub_Area[column]) {
      # } else if (abs(sample_metadata_sorted_sumLS$Geographic_Order[row] - sample_metadata_sorted_sumLS$Geographic_Order[column]) == 1 & sample_metadata_sorted_sumLS$Sub_Area[row] == sample_metadata_sorted_sumLS$Sub_Area[column]) {
      Djost_df_3D_fronts[row, column] <- 1
    } else if (abs(metadata$Geographic_Order[row] - metadata$Geographic_Order[column]) == 1 & metadata$Sub_Area[row] != metadata$Sub_Area[column]) {
      # } else if (abs(sample_metadata_sorted_sumLS$Geographic_Order[row] - sample_metadata_sorted_sumLS$Geographic_Order[column]) == 1 & sample_metadata_sorted_sumLS$Sub_Area[row] != sample_metadata_sorted_sumLS$Sub_Area[column]) {
      Djost_df_3D_fronts[row, column] <- 2
    } else  {
      Djost_df_3D_fronts[row, column] <- 0
    }
  }
}


run <- 0

motus <- unique(initial_data$motu)

sample_cols <- c(5:ncol(initial_data))

for(id in motus){
  # id <- str_extract(id, 'PHY1_.........')
  # I import the MOTUs that will be analyzed
  # to the files the haplotypes and to the columns you show them
  
  data_inicial <- initial_data[initial_data$motu==id,]
  # each MOTU has been connected to the DnoisE program https://github.com/adriantich/DnoisE
  # first hi has a minimal abundance filter that is fet in the script filtering_form_DnoisE.R
  # The remaining haplotypes are corrected
  
  
  # we pass the data from the reads to relative values to be able to do the transformation
  # after values from 0 to 4
  rel_data <- data_inicial[,sample_cols]
  positive_samples <- colSums(rel_data)>0
  
  df4genind <- c()
  df4genind_pop <- c()
  for (i in 1:dim(rel_data)[1]) {
    for (j in which(positive_samples, TRUE)) {
      individuals <- data.frame('loc1' = rep(as.character(data_inicial$id[i]),rel_data[i,j]))
      df4genind <- rbind(df4genind,individuals)
      df4genind_pop <- c(df4genind_pop,rep(colnames(rel_data)[j],rel_data[i,j]))
    }
  }
  genindobj <- df2genind(df4genind,pop = df4genind_pop, ploidy = 1)
  Djost_df <- as.matrix(mmod::pairwise_D(genindobj))
  Djost_df[Djost_df<0] <- 0
  
  # add absent samples
  Djost_df <- rbind(Djost_df,
                    matrix(NA,
                           nrow = (ncol(rel_data) - dim(Djost_df)[1]),
                           ncol = (dim(Djost_df)[2])))
  Djost_df <- cbind(Djost_df,
                    matrix(NA,
                           nrow = (dim(Djost_df)[1]),
                           ncol = (ncol(rel_data) - dim(Djost_df)[2])))
  colnames(Djost_df)[colnames(Djost_df)==''] <- names(which(!positive_samples))
  rownames(Djost_df)[rownames(Djost_df)==''] <- names(which(!positive_samples))
  Djost_df <- Djost_df[,as.character(metadata$Community)]
  Djost_df <- Djost_df[as.character(metadata$Community),]
  
  
  # Djost_df_sum[positive_samples, positive_samples] <- Djost_df_sum[positive_samples, positive_samples] + Djost_df[positive_samples, positive_samples]
  Djost_df_3D <- abind::abind(Djost_df_3D,Djost_df, along = 3)
  
  if (correr_sig) {
    # Djost_df_recount[positive_samples, positive_samples] <- Djost_df_recount[positive_samples, positive_samples] + 1
    # compute significance from randomization
    genindobj <- df2genind(df4genind,pop = df4genind_pop, ploidy = 1)
    Djost_df <- as.matrix(mmod::pairwise_D(genindobj))
    plan(multisession)
    rndmlist <- future_lapply(1:1000, FUN=function(x){
      as.matrix(
        mmod::pairwise_D(
          df2genind(
            as.data.frame(df4genind[sample(1:length(df4genind$loc1)),]),
            pop = df4genind_pop, ploidy = 1)))
    })
    d3array <- abind::abind(rndmlist,along = 3)
    CIarray <- apply(d3array, 1:2, function(x){
      x[x<0] <- 0
      quantile(x,.975)}
    )
    Djost_sig <- CIarray<Djost_df
    Djost_sig <- rbind(Djost_sig,
                       matrix(NA,
                              nrow = (n_sites - dim(Djost_sig)[1]),
                              ncol = (dim(Djost_sig)[2])))
    Djost_sig <- cbind(Djost_sig,
                       matrix(NA,
                              nrow = (dim(Djost_sig)[1]),
                              ncol = (n_sites - dim(Djost_sig)[2])))
    colnames(Djost_sig)[colnames(Djost_sig)==''] <- names(which(!positive_samples))
    rownames(Djost_sig)[rownames(Djost_sig)==''] <- names(which(!positive_samples))
    Djost_sig <- Djost_sig[,as.character(metadata$Community)]
    Djost_sig <- Djost_sig[as.character(metadata$Community),]
    
    Djost_sig_3D <- abind::abind(Djost_sig_3D,Djost_sig, along = 3)
  } 
  
  # here I am creating a data frame with the averages and errors for each motu of the Djost values within front and between fronts
  if (run == 0) {
    mitjanes_Djost_dins_fora <- data.frame("id" = id, "mean_Djost" = mean(Djost_df[Djost_df_3D_fronts > 0],na.rm = TRUE), 
                                           "intra_front" = mean(Djost_df[Djost_df_3D_fronts == 1],na.rm = TRUE), 
                                           "inter_front" = mean(Djost_df[Djost_df_3D_fronts == 2],na.rm = TRUE),
                                           "err_Djost" = errstd(Djost_df[Djost_df_3D_fronts > 0]),"err_intra_front" = errstd(Djost_df[Djost_df_3D_fronts == 1]), 
                                           "err_inter_front" = errstd(Djost_df[Djost_df_3D_fronts == 2]),
                                           'total_reads'=all_data_reads[grep(id, names(all_data_reads))])
  } else {
    mitjanes_Djost_dins_fora <- rbind(mitjanes_Djost_dins_fora,
                                      data.frame("id" = id, "mean_Djost" = mean(Djost_df[Djost_df_3D_fronts > 0],na.rm = TRUE), 
                                                 "intra_front" = mean(Djost_df[Djost_df_3D_fronts == 1],na.rm = TRUE), 
                                                 "inter_front" = mean(Djost_df[Djost_df_3D_fronts == 2],na.rm = TRUE),
                                                 "err_Djost" = errstd(Djost_df[Djost_df_3D_fronts > 0]),"err_intra_front" = errstd(Djost_df[Djost_df_3D_fronts == 1]), 
                                                 "err_inter_front" = errstd(Djost_df[Djost_df_3D_fronts == 2]),
                                                 'total_reads'=all_data_reads[grep(id, names(all_data_reads))]))
    
  }
  
  run <- run +1 
  print(paste((run*100/length(motus)),'% done'))
  save.image('../../../RDS/data_during_run_presabs_loc.RData')
}
save.image('../../../RDS/data_after_run_presabs_loc.RData')

#load('data_after_run_presabs_loc.RData')

#UNTIL HERE

mitjanes_Djost_dins_fora <- cbind(mitjanes_Djost_dins_fora,data.frame("sep"=NA))
mitjanes_Djost_dins_fora$sep[c(mitjanes_Djost_dins_fora$intra_front + mitjanes_Djost_dins_fora$err_intra_front)<c(mitjanes_Djost_dins_fora$inter_front - mitjanes_Djost_dins_fora$err_inter_front)] <- "a"
mitjanes_Djost_dins_fora$sep[!c(c(mitjanes_Djost_dins_fora$intra_front + mitjanes_Djost_dins_fora$err_intra_front)<c(mitjanes_Djost_dins_fora$inter_front - mitjanes_Djost_dins_fora$err_inter_front))] <- "b"
# mitjanes_Djost_dins_fora$sep[c(mitjanes_Djost_dins_fora$intra_front)<c(mitjanes_Djost_dins_fora$inter_front)] <- "a"
# mitjanes_Djost_dins_fora$sep[!c(c(mitjanes_Djost_dins_fora$intra_front )<c(mitjanes_Djost_dins_fora$inter_front ))] <- "b"
# mitjanes_Djost_dins_fora$sep[c(mitjanes_Djost_dins_fora$intra_front<0.25)&c(mitjanes_Djost_dins_fora$inter_front>0.5)] <- "a"
# mitjanes_Djost_dins_fora$sep[!c(c(mitjanes_Djost_dins_fora$intra_front<0.25)&c(mitjanes_Djost_dins_fora$inter_front>0.5))] <- "b"

mitjanes_Djost_dins_fora <- mitjanes_Djost_dins_fora[!is.nan(mitjanes_Djost_dins_fora$inter_front) & !is.nan(mitjanes_Djost_dins_fora$intra_front),]

mitjanes_Djost_dins_fora$semiquant_mean <- sapply(mitjanes_Djost_dins_fora$mean_Djost, FUN = semi_quant)
mitjanes_Djost_dins_fora$semiquant_read <- sapply(mitjanes_Djost_dins_fora$total_reads, FUN = semi_quantb)


# Djost_df_mean <- Djost_df_sum / Djost_df_recount
Djost_df_mean <- apply(Djost_df_3D,1:2,FUN = function(x) mean(x,na.rm = T))
colnames(Djost_df_mean) <- metadata$Community
rownames(Djost_df_mean) <- metadata$Community
if (correr_sig) {
  Djost_sig_sum <- apply(Djost_sig_3D,1:2,FUN = function(x) sum(x,na.rm = T))
  colnames(Djost_sig_sum) <- metadata$Community
  rownames(Djost_sig_sum) <- metadata$Community
  Djost_sig_mean <- apply(Djost_sig_3D,1:2,FUN = function(x) mean(x,na.rm = T))
  colnames(Djost_sig_mean) <- metadata$Community
  rownames(Djost_sig_mean) <- metadata$Community
  Djost_sig_count <- apply(Djost_sig_3D,1:2,FUN = function(x) sum(!is.na(x),na.rm = T))
  colnames(Djost_sig_count) <- metadata$Community
  rownames(Djost_sig_count) <- metadata$Community
}
Djost_df_median <- apply(Djost_df_3D,1:2,FUN = function(x) median(x,na.rm = T))
colnames(Djost_df_median) <- metadata$Community
rownames(Djost_df_median) <- metadata$Community

Djost_df_3D_sig <- Djost_df_3D
Djost_df_3D_sig[Djost_sig_3D==FALSE] <- NA
Djost_df_mean_sig <- apply(Djost_df_3D_sig,1:2,FUN = function(x) mean(x,na.rm = T))
colnames(Djost_df_mean_sig) <- metadata$Community
rownames(Djost_df_mean_sig) <- metadata$Community
Djost_df_median_sig <- apply(Djost_df_3D_sig,1:2,FUN = function(x) median(x,na.rm = T))
colnames(Djost_df_median_sig) <- metadata$Community
rownames(Djost_df_median_sig) <- metadata$Community



# mean Djost values
Djost_df_no0 <- Djost_df_mean
sample_x <- rownames(Djost_df_no0)
Djost_df_no0 <- cbind(sample_x,as.data.frame(Djost_df_no0))
Djost_df_no0$sample_x <- factor(Djost_df_no0$sample_x, levels = as.character(Djost_df_no0$sample_x))

if(correr_sig){
  # sum of number of significative Djost
  Djost_sum_no0 <- Djost_sig_sum
  sample_x <- rownames(Djost_sum_no0)
  Djost_sum_no0 <- cbind(sample_x,as.data.frame(Djost_sum_no0))
  Djost_sum_no0$sample_x <- factor(Djost_sum_no0$sample_x, levels = as.character(Djost_sum_no0$sample_x))
  
  # mean of number of significative Djost among all counts of Djost
  Djost_mean_no0 <- Djost_sig_mean
  sample_x <- rownames(Djost_mean_no0)
  Djost_mean_no0 <- cbind(sample_x,as.data.frame(Djost_mean_no0))
  Djost_mean_no0$sample_x <- factor(Djost_mean_no0$sample_x, levels = as.character(Djost_mean_no0$sample_x))
  
  # number of Djost
  Djost_count_no0 <- Djost_sig_count
  sample_x <- rownames(Djost_count_no0)
  Djost_count_no0 <- cbind(sample_x,as.data.frame(Djost_count_no0))
  Djost_count_no0$sample_x <- factor(Djost_count_no0$sample_x, levels = as.character(Djost_count_no0$sample_x))
  
}

# median Djost values
Djost_median_no0 <- Djost_df_median
sample_x <- rownames(Djost_median_no0)
Djost_median_no0 <- cbind(sample_x,as.data.frame(Djost_median_no0))
Djost_median_no0$sample_x <- factor(Djost_median_no0$sample_x, levels = as.character(Djost_median_no0$sample_x))

# mean of significative Djost values
Djost_mean_sig_no0 <- Djost_df_mean_sig
sample_x <- rownames(Djost_mean_sig_no0)
Djost_mean_sig_no0 <- cbind(sample_x,as.data.frame(Djost_mean_sig_no0))
Djost_mean_sig_no0$sample_x <- factor(Djost_mean_sig_no0$sample_x, levels = as.character(Djost_mean_sig_no0$sample_x))

# median of significative Djost values
Djost_median_sig_no0 <- Djost_df_median_sig
sample_x <- rownames(Djost_median_sig_no0)
Djost_median_sig_no0 <- cbind(sample_x,as.data.frame(Djost_median_sig_no0))
Djost_median_sig_no0$sample_x <- factor(Djost_median_sig_no0$sample_x, levels = as.character(Djost_median_sig_no0$sample_x))

Djost_gathered <- gather(Djost_df_no0, "sample_y", "Djost", -sample_x)
if(correr_sig){
  Djost_gathered_sig <- gather(Djost_sum_no0, "sample_y", "Djost_sig", -sample_x)
  Djost_gathered_sig_mean <- gather(Djost_mean_no0, "sample_y", "Djost_sig_mean", -sample_x)
  Djost_gathered_sig_count <- gather(Djost_count_no0, "sample_y", "Djost_sig_count", -sample_x)
}
Djost_gathered_median <- gather(Djost_median_no0, "sample_y", "Djost_median", -sample_x)
Djost_gathered_mean_sig <- gather(Djost_mean_sig_no0, "sample_y", "Djost_mean_sig", -sample_x)
Djost_gathered_median_sig <- gather(Djost_median_sig_no0, "sample_y", "Djost_median_sig", -sample_x)

Djost_gathered$sample_x <- factor(Djost_gathered$sample_x, levels = as.character(Djost_df_no0$sample_x))
Djost_gathered$sample_y <- factor(Djost_gathered$sample_y, levels = as.character(Djost_df_no0$sample_x))
if(correr_sig){
  Djost_gathered$Djost_sig <- Djost_gathered_sig$Djost_sig
  Djost_gathered$Djost_sig_mean <- Djost_gathered_sig_mean$Djost_sig_mean
  Djost_gathered$Djost_sig_count <- Djost_gathered_sig_count$Djost_sig_count
}
Djost_gathered$Djost_median <- Djost_gathered_median$Djost_median
Djost_gathered$Djost_mean_sig <- Djost_gathered_mean_sig$Djost_mean_sig
Djost_gathered$Djost_median_sig <- Djost_gathered_median_sig$Djost_median_sig

# Djost_gathered <- Djost_gathered[as.numeric(Djost_gathered$sample_x)>=as.numeric(Djost_gathered$sample_y),]

# Djost_gathered <- cbind(Djost_gathered, "Sd" = Djost_gathered_sd$Sd)
Djost_gathered$com_x <- NA
Djost_gathered$com_y <- NA
Djost_gathered$reg_x <- NA
Djost_gathered$reg_y <- NA
Djost_gathered$dist_x <- NA
Djost_gathered$dist_y <- NA
Djost_gathered$fronts <- as.character(0)


# for (com in levels(metadata$Community)) {
for (com in metadata$Community) {
  Djost_gathered$com_x[grep(com,Djost_gathered$sample_x)] <- com
  Djost_gathered$com_y[grep(com,Djost_gathered$sample_y)] <- com
  # reg <- as.character(metadata$Sub_Area[grep(com,metadata$Community)][1])
  # distance <- metadata$dist_from_C33.km.[grep(com,metadata$Community)][1]
  reg <- as.character(metadata$Sub_Area[grep(com,metadata$Community)][1])
  distance <- metadata$dist_from_C33.km[grep(com,metadata$Community)][1]
  Djost_gathered$reg_x[grep(com,Djost_gathered$sample_x)] <- reg
  Djost_gathered$reg_y[grep(com,Djost_gathered$sample_y)] <- reg
  Djost_gathered$dist_x[grep(com,Djost_gathered$sample_x)] <- distance
  Djost_gathered$dist_y[grep(com,Djost_gathered$sample_y)] <- distance
  
}

Djost_gathered$same_com <- Djost_gathered$com_x == Djost_gathered$com_y
Djost_gathered$same_reg <- Djost_gathered$reg_x == Djost_gathered$reg_y
Djost_gathered$one_front <- Djost_gathered$reg_x != Djost_gathered$reg_y & c(Djost_gathered$reg_x == 'Canal_Ibiza_S' | Djost_gathered$reg_y == 'Canal_Ibiza_S')
Djost_gathered$two_front <- Djost_gathered$reg_x != Djost_gathered$reg_y & c(Djost_gathered$reg_x != 'Canal_Ibiza_S' & Djost_gathered$reg_y != 'Canal_Ibiza_S')
Djost_gathered$fronts[Djost_gathered$one_front] <- 1
Djost_gathered$fronts[Djost_gathered$two_front] <- 2
Djost_gathered$dist_xy <- abs(c(Djost_gathered$dist_x - Djost_gathered$dist_y))

########
## sd ##
########

Djost_sd <- apply(Djost_df_3D, c(1,2), errstd )
Djost_df <- as.data.frame(Djost_df)
Djost_df_no0_sd <- Djost_sd
sample_x <- rownames(Djost_df_no0_sd)
Djost_df_no0_sd <- cbind(sample_x,as.data.frame(Djost_df_no0_sd))
Djost_df_no0_sd$sample_x <- factor(Djost_df_no0_sd$sample_x, levels = as.character(Djost_df_no0_sd$sample_x))

Djost_gathered_sd <- gather(Djost_df_no0_sd, "sample_y", "Sd", -sample_x)
Djost_gathered_sd$sample_x <- factor(Djost_gathered_sd$sample_x, levels = as.character(Djost_df_no0_sd$sample_x))
Djost_gathered_sd$sample_y <- factor(Djost_gathered_sd$sample_y, levels = as.character(Djost_df_no0_sd$sample_x))
Djost_gathered_sd$com_x <- NA
Djost_gathered_sd$com_y <- NA
Djost_gathered_sd$reg_x <- NA
Djost_gathered_sd$reg_y <- NA
Djost_gathered_sd$dist_x <- NA
Djost_gathered_sd$dist_y <- NA
Djost_gathered_sd$fronts <- as.character(0)


# for (com in levels(metadata$Community)) {
for (com in metadata$Community) {
  Djost_gathered_sd$com_x[grep(com,Djost_gathered_sd$sample_x)] <- com
  Djost_gathered_sd$com_y[grep(com,Djost_gathered_sd$sample_y)] <- com
  # reg <- as.character(metadata$Sub_Area[grep(com,metadata$Community)][1])
  # distance <- metadata$dist_from_C33.km.[grep(com,metadata$Community)][1]
  reg <- as.character(metadata$Sub_Area[grep(com,metadata$Community)][1])
  distance <- metadata$dist_from_C33.km[grep(com,metadata$Community)][1]
  Djost_gathered_sd$reg_x[grep(com,Djost_gathered_sd$sample_x)] <- reg
  Djost_gathered_sd$reg_y[grep(com,Djost_gathered_sd$sample_y)] <- reg
  Djost_gathered_sd$dist_x[grep(com,Djost_gathered_sd$sample_x)] <- distance
  Djost_gathered_sd$dist_y[grep(com,Djost_gathered_sd$sample_y)] <- distance
  
}

Djost_gathered_sd$same_com <- Djost_gathered_sd$com_x == Djost_gathered_sd$com_y
Djost_gathered_sd$same_reg <- Djost_gathered_sd$reg_x == Djost_gathered_sd$reg_y
Djost_gathered_sd$one_front <- Djost_gathered_sd$reg_x != Djost_gathered_sd$reg_y & c(Djost_gathered_sd$reg_x == 'Canal_Ibiza_S' | Djost_gathered_sd$reg_y == 'Canal_Ibiza_S')
Djost_gathered_sd$two_front <- Djost_gathered_sd$reg_x != Djost_gathered_sd$reg_y & c(Djost_gathered_sd$reg_x != 'Canal_Ibiza_S' & Djost_gathered_sd$reg_y != 'Canal_Ibiza_S')
Djost_gathered_sd$fronts[Djost_gathered_sd$one_front] <- 1
Djost_gathered_sd$fronts[Djost_gathered_sd$two_front] <- 2
Djost_gathered_sd$dist_xy <- abs(c(Djost_gathered_sd$dist_x - Djost_gathered_sd$dist_y))

# Dropped the below, as there are no negative values
#Djost_gathered$Djost[Djost_gathered$Djost<0] <- Djost_gathered$Djost[Djost_gathered$Djost<0] * (-1)



# save(sq,semi_quant,semi_quantb,mitjanes_Djost_dins_fora,Djost_df_3D, Djost_df_recount, Djost_df_sum, Djost_gathered, Djost_df_mean, file="Djost_meanABC.RData")
if(correr_sig){
    save(sq,semi_quant,semi_quantb,mitjanes_Djost_dins_fora,Djost_df_3D,
       Djost_gathered, Djost_df_mean, Djost_sig_3D,Djost_sig_sum,
       Djost_df_median,Djost_df_mean_sig, Djost_df_median_sig,Djost_sig_mean,
       file="../../../RDS/Djost_meanABC_presabs_loc.RData")
} else {
  save(sq,semi_quant,semi_quantb,mitjanes_Djost_dins_fora,Djost_df_3D,
       Djost_gathered, Djost_df_mean,
       Djost_df_median,Djost_df_mean_sig, Djost_df_median_sig,
       file="../../../RDS/Djost_meanABC_presabs_loc.RData")
}

write.csv(Djost_df_mean,'../../../Tekstfiler/COI/COI_ASV/Sed_Djost.txt')

# Calculate number of comparisons (species) per site pair
test <- Djost_df_3D
test[!is.na(test)] <- 1 # Everything not NA will get the value 1
out=matrix(NA,n_sites,n_sites) 
for(n in 1:n_sites) { 
  for(i in 1:n_sites){out[i,n]=sum(test[i,n,],na.rm=T)}} 
diag(out)=NA
row.names(out)<-row.names(Djost_df)
colnames(out)<-colnames(Djost_df)
hist(out,50)
mean(out,na.rm=T)
write.table(out, "../../../Tekstfiler/COI/COI_ASV/Sediment_no_of_comp.txt",sep="\t", quote=F) 
