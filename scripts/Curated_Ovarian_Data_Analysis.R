# Analysis of PVR/CD274 survival associations using curatedOvarianData package
# Mostly followed package vignette
# Updated to reflect reviewer comments
# removed datasets from PVR analysis so same datasets are used for both PVR and CD274

library(data.table)
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(cowplot)
library(curatedOvarianData)
library(logging)
library(metafor)

# source data
source(system.file("extdata", "patientselection.config", package="curatedOvarianData"))

sapply(ls(), function(x) if(!x %in% c("remove.samples", "duplicates")) print(get(x)))
source(system.file("extdata", "createEsetList.R", package = "curatedOvarianData"))

# modify forestplot function
# require PD1/CD274 presence 
forestplot <- function(esets, y="y", probeset, formula=y~probeset,
                       mlab="Overall", rma.method="FE", at=NULL,xlab="Hazard Ratio",...) {
  require(metafor)
  esets <- esets[sapply(esets, function(x) c('CD274') %in% featureNames(x))]
  coefs <- sapply(1:length(esets), function(i) {
    tmp   <- as(phenoData(esets[[i]]), "data.frame")
    tmp$y <- esets[[i]][[y]]
    tmp$probeset <- exprs(esets[[i]])[probeset,]
    
    summary(coxph(formula, data=tmp))$coefficients[1, c(1,3)]
  })  
  
  res.rma <- metafor::rma(yi = coefs[1,], sei = coefs[2,], 
                          method=rma.method)
  
  if (is.null(at)) at <- c(-1,0,1)
  forest.rma(res.rma, xlab=xlab, slab=gsub("_eset$", "", names(esets)), at=at, mlab=mlab,...)
  return(res.rma)
}


# look at associations in all datasets
res <- forestplot(esets=esets, probeset="PVR")
res$pval

# screen to keep datasets: debulking of most interest
# as it often strongly predicts survival
idx.debulking <- sapply(esets, function(X) 
  sum(X$debulking=="suboptimal", na.rm=TRUE)) > 0


# include debulking in forestplot survival model
# CD155 / PVR
res <- forestplot(esets=esets[idx.debulking],
                  probeset="PVR",
                  formula=y~probeset+debulking, 
                  xlab = 'log hazard ratio')
res$pval
dev.copy2pdf(file='figures/CD155_Ovarian_Survival_Revised.pdf', height=4.75, width=6)

# PD-L1 / CD274
res <- forestplot(esets=esets[idx.debulking],
                  probeset="CD274",formula=y~probeset+debulking,
                  xlab = 'log hazard ratio')
res$pval # pval
dev.copy2pdf(file='figures/PDL1_Ovarian_Survival_Revised.pdf', height=4, width=6)
