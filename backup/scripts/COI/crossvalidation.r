'
crossvalidation
'
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


#water_order<- readRDS("results/COI_rich_p50prevorder_water_plus_c2_thin10_s2000_tr02.rds")
#sed_order_20<- readRDS("results/COI_rich_p50prevorder_sediment_plus_c2_thin10_s2000_tr02.rds")

sed_order_thin40_s60k <- readRDS("results/sediment/COI_rich_p50prevorder_sediment_plus_c2_thin40_s6000_tr03.rds")
water_order_thin40_s60k<- readRDS("results/Water/order/COI_rich_p50prevorder_water_plus_c2_thin40_s5000_tr03_notnorm")


model<-sed_order_thin40_s60k

#models to cv.
models<- list(model)

cross_val<- sapply(models,
    function(model, frame){
        #the number of folds:
        k=2
        #we start by creating the data partition( note that the computational time will increase k by amount of k folds.)
        folds<- createPartition(model,nfolds=k)
        #We make our predictions using the different data partitions. 
        preds = computePredictedValues(model,partition=folds)
        MF= evaluateModelFit(hM =model, predY = preds)

        saveRDS(MF, file =paste0("results/","COI_rich_p50prevorder_water_plus_c2_thin10_s2000_tr02_notnorm",model,".rds"))
        print("Mean of SR2")
        mean(MF$SR2)

        print("Std dev of SR2")
        sd(MF$SR2)

        print("Mean of PSRF")
        mean(gd$psrf)

        frame<- data.frame(
            "pred_RMSE"=mean(MF$SR2),
            "Std dev of SR2"=sd(MF$SR2),
            "mean psrg"=mean(gd$psrf)
            )
        frame
        }
        )
        


cross_val