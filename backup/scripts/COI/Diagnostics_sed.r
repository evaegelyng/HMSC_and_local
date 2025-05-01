# Activate conda environment before starting R:
# conda activate hmsc-hpc
# Run from results folder

library("jsonify")
library("Hmsc")
library(colorspace)
library(vioplot)

load(file = "sediment/COI/unfitted_models.RData")

chain1 <- from_json(readRDS(file="sediment/COI/post_file00.rds")[[1]])
chain2 <- from_json(readRDS(file="sediment/COI/post_file01.rds")[[1]])
chain3 <- from_json(readRDS(file="sediment/COI/post_file02.rds")[[1]])
chain4 <- from_json(readRDS(file="sediment/COI/post_file03.rds")[[1]])

cat(sprintf("fitting time %.1f sec\n", chain1[[2]]))
cat(sprintf("fitting time %.1f sec\n", chain2[[2]]))
cat(sprintf("fitting time %.1f sec\n", chain3[[2]]))
cat(sprintf("fitting time %.1f sec\n", chain4[[2]]))

postList = list(chain1[[1]],chain2[[1]],chain3[[1]],chain4[[1]])

nSamples = 250
thin = 1000
transient = nSamples*thin

fitTF = importPosteriorFromHPC(models[[1]], postList, nSamples, thin, transient)
saveRDS(fitTF, file = "sediment/COI/fitTF.rds")

# MAKE THE SCRIPT REPRODUCIBLE 
set.seed(1)

#### Evaluate convergence #### 

# SETTING COMMONLY ADJUSTED PARAMETERS 
showBeta = TRUE
showGamma = TRUE
showOmega = TRUE
maxOmega = 100
showRho = TRUE
showAlpha = TRUE

text.file = file.path("sediment/COI/MCMC_convergence.txt")
cat("MCMC Convergence statistics\n\n",file=text.file,sep="")

ma.beta = NULL
na.beta = NULL
ma.gamma = NULL
na.gamma = NULL
ma.omega= NULL
na.omega = NULL
ma.alpha = NULL
na.alpha = NULL  
ma.rho = NULL
na.rho = NULL

mpost = convertToCodaObject(fitTF, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
nr = fitTF$nr
if(showBeta){
  psrf = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
  tmp = summary(psrf)
  cat("\nbeta\n\n",file=text.file,sep="",append=TRUE)
  cat(tmp[,1],file=text.file,sep="\n",append=TRUE)
  if(is.null(ma.beta)){
    ma.beta = psrf[,1]
    na.beta = paste0(as.character(thin),",",as.character(nSamples))
  } else {
    ma.beta = cbind(ma.beta,psrf[,1])
    na.beta = c(na.beta,paste0(as.character(thin),",",as.character(nSamples)))
  }
}
if(showGamma){
  psrf = gelman.diag(mpost$Gamma,multivariate=FALSE)$psrf
  tmp = summary(psrf)
  cat("\ngamma\n\n",file=text.file,sep="",append=TRUE)
  cat(tmp[,1],file=text.file,sep="\n",append=TRUE)
  if(is.null(ma.gamma)){
    ma.gamma = psrf[,1]
    na.gamma = paste0(as.character(thin),",",as.character(nSamples))
  } else {
    ma.gamma = cbind(ma.gamma,psrf[,1])
    if(j==1){
      na.gamma = c(na.gamma,paste0(as.character(thin),",",as.character(nSamples)))
    } else {
      na.gamma = c(na.gamma,"")
    }
  }
}
if(showRho & !is.null(mpost$Rho)){
  psrf = gelman.diag(mpost$Rho,multivariate=FALSE)$psrf
  cat("\nrho\n\n",file=text.file,sep="",append=TRUE)
  cat(psrf[1],file=text.file,sep="\n",append=TRUE)
}
if(showOmega & nr>0){
  cat("\nomega\n\n",file=text.file,sep="",append=TRUE)
  for(k in 1:nr){
    cat(c("\n",names(fitTF$ranLevels)[k],"\n\n"),file=text.file,sep="",append=TRUE)
    tmp = mpost$Omega[[k]]
    z = dim(tmp[[1]])[2]
    if(z > maxOmega){
      sel = sample(1:z, size = maxOmega)
      for(i in 1:length(tmp)){
        tmp[[i]] = tmp[[i]][,sel]
      }
    }
    psrf = gelman.diag(tmp, multivariate = FALSE)$psrf
    tmp = summary(psrf)
    cat(tmp[,1],file=text.file,sep="\n",append=TRUE)
    if(is.null(ma.omega)){
      ma.omega = psrf[,1]
      na.omega = paste0(as.character(thin),",",as.character(nSamples))
    } else {
      ma.omega = cbind(ma.omega,psrf[,1])
      na.omega = c(na.omega,paste0(as.character(thin),",",as.character(nSamples)))
    }
  }
}
if(showAlpha & nr>0){
  for(k in 1:nr){
    if(fitTF$ranLevels[[k]]$sDim>0){
      cat("\nalpha\n\n",file=text.file,sep="\n",append=TRUE)
      cat(c("\n",names(fitTF$ranLevels)[k],"\n\n"),file=text.file,sep="",append=TRUE)
      psrf = gelman.diag(mpost$Alpha[[k]],multivariate = FALSE)$psrf
      cat(psrf[,1],file=text.file,sep="\n",append=TRUE)            
    }
  }
}

nm = length(models)

pdf(file= file.path("sediment/COI/MCMC_convergence.pdf"))
if(showBeta){
  par(mfrow=c(2,1))
  vioplot(ma.beta,col=rainbow_hcl(nm),names=na.beta,ylim=c(0,max(ma.beta)),main="psrf(beta)")
  legend("topright",legend = names(models), fill=rainbow_hcl(nm))
  vioplot(ma.beta,col=rainbow_hcl(nm),names=na.beta,ylim=c(0.9,1.1),main="psrf(beta)")
}
if(showGamma){
  par(mfrow=c(2,1))
  vioplot(ma.gamma,col=rainbow_hcl(nm),names=na.gamma,ylim=c(0,max(ma.gamma)),main="psrf(gamma)")
  legend("topright",legend = names(models), fill=rainbow_hcl(nm))
  vioplot(ma.gamma,col=rainbow_hcl(nm),names=na.gamma,ylim=c(0.9,1.1),main="psrf(gamma)")
}
if(showOmega & !is.null(ma.omega)){
  par(mfrow=c(2,1))
  vioplot(ma.omega,col=rainbow_hcl(nm),names=na.omega,ylim=c(0,max(ma.omega)),main="psrf(omega)")
  legend("topright",legend = names(models), fill=rainbow_hcl(nm))
  vioplot(ma.omega,col=rainbow_hcl(nm),names=na.omega,ylim=c(0.9,1.1),main="psrf(omega)")
}
dev.off()

#### Compute model fit ####

nfolds = 2 #change the number of CV-folds
nParallel = 4 #Set to 1 to disable parallel computing

preds = computePredictedValues(fitTF)
MF = evaluateModelFit(hM=fitTF, predY=preds)
partition = createPartition(fitTF, nfolds = nfolds, column="Time_d") 
preds = computePredictedValues(fitTF,partition=partition, nParallel = nParallel)
MFCV = evaluateModelFit(hM=fitTF, predY=preds)
WAIC = computeWAIC(fitTF)
      
save(MF,MFCV,WAIC,file = "sediment/COI/model_fit.Rdata")