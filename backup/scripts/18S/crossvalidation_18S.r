'
crossvalidation
'
library("Hmsc")

#model <- readRDS("results/sediment/18S/18S_rich_p50prevclass_sediment_plus_c2_thin10_s6000_tr02_untrans.rds")
model<- readRDS("results/Water/18S/18S_rich_p50prevclass_water_plus_c2_thin10_s6000_tr02_untrans.rds")

#we start by creating the data partition( note that the computational time will increase k by amount of k folds.)
folds<- createPartition(model,nfolds=2,column="Time_d")

#We make our predictions using the different data partitions. 
preds = computePredictedValues(model,partition=folds)
MF= evaluateModelFit(hM = model, predY = preds)

#saveRDS(MF, file = "results/sediment/18S/18S_rich_p50prevclass_sediment_untrans_cv_time_d.rds")
saveRDS(MF, file = "results/Water/18S/18S_rich_p50prevclass_water_untrans_cv_time_d.rds")

#Explanatory power
#For Poisson models, a pseudo-R2 is computed as squared spearman correlation between observed and
# predicted values, times the sign of the correlation (SR2).
print("Mean of SR2")
mean(MF$SR2)
print("Std dev of SR2")
sd(MF$SR2)