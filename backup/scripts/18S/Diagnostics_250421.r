library(jsonlite)
library("Hmsc")

# Diagnostics
#unfit_m <- readRDS(file = "results/Water/18S/unfitted_models.rds") # doesnt work!
#unfit_m <- load(file = "results/Water/18S/unfitted_models.RData") # doesnt work!

chain1 <- from_json(readRDS(file="results/Water/18S/post_file00.rds")[[1]])
chain2 <- from_json(readRDS(file="results/Water/18S/post_file01.rds")[[1]])
chain3 <- from_json(readRDS(file="results/Water/18S/post_file02.rds")[[1]])
chain4 <- from_json(readRDS(file="results/Water/18S/post_file03.rds")[[1]])

cat(sprintf("fitting time %.1f sec\n", chain1[[2]]))
cat(sprintf("fitting time %.1f sec\n", chain2[[2]]))
cat(sprintf("fitting time %.1f sec\n", chain3[[2]]))
cat(sprintf("fitting time %.1f sec\n", chain4[[2]]))

postList = list(chain1[[1]],chain2[[1]],chain3[[1]],chain4[[1]])

fitTF = importPosteriorFromHPC(m, postList, nSamples, thin, transient)

jpeg("results/Water/18S/var_par_18S_wat.jpg")
plotVariancePartitioning(fitTF, computeVariancePartitioning(fitTF), args.legend=list(x="bottomright"))
dev.off()

preds = computePredictedValues(fitTF)

#To get the sr2, the predicted values are used to estimates amount of explained variance.
MF= evaluateModelFit(hM = fitTF, predY = preds)

#Converting model object to r interpretable coda
mpost = convertToCodaObject(fitTF, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
ess.beta = effectiveSize(mpost$Beta)

#Explanatory power
#For Poisson models, a pseudo-R2 is computed as squared spearman correlation between observed and
# predicted values, times the sign of the correlation (SR2).
print("Mean of SR2")
mean(MF$SR2,na.rm = TRUE)
print("Std dev of SR2")
sd(MF$SR2, na.rm = TRUE)

#Predictive power
print("Mean of ESS")
mean(ess.beta)

print("Std dev of ESS")
sd(ess.beta)

#Histogram of effective sample size
pdf("results/Water/18S/diag_250421.pdf")
par(mfrow=c(1,2))
hist(effectiveSize(mpost$Beta), main="ess(beta)")
hist(gelman.diag(mpost$Beta, multivariate=F)$psrf, main="psrf(beta)")
dev.off()

#Posterior beta estimation over different iterations gives an idea of parameter convergence for each parameter.
pdf("results/Water/18S/diag_beta_250421.pdf")
plot(mpost$Beta)
dev.off()

head(fitTF$X)