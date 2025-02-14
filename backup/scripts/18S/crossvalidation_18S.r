'
crossvalidation
'
library("Hmsc")

model <- readRDS("results/sediment/18S/18S_sed_241007.rds")
#model<- readRDS("results/Water/18S/18S_wat_241007.rds")

#we start by creating the data partition( note that the computational time will increase k by amount of k folds.)
folds<- createPartition(model,nfolds=10,column="Time_d")

#We make our predictions using the different data partitions. 
preds = computePredictedValues(model,partition=folds)
MF= evaluateModelFit(hM = model, predY = preds)

saveRDS(MF, file = "results/sediment/18S/18S_sed_241007_cv.rds")
#saveRDS(MF, file = "results/Water/18S/18S_wat_241007_cv.rds")

#Explanatory power
#For Poisson models, a pseudo-R2 is computed as squared spearman correlation between observed and
# predicted values, times the sign of the correlation (SR2).
print("Mean of SR2")
mean(MF$SR2)
print("Std dev of SR2")
sd(MF$SR2)