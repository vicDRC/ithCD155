## Evaluation of CD155 expression in in-house IROC cohort
# includes data preprocessing: processed data provided in data folder

# required packages
library(pheatmap)
library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)
library(limma)

# set wd
wd <- '~/Documents/R/CD155_final_draft/'
setwd(wd)

# #--------------------------------------------#
# #----------- tabulate patient data ----------#
 pdat <- fread('data/patient_characteristics.csv')
# 
# ageMean <- mean(pdat$Age)
# ageSD <- sd(pdat$Age)
# table(pdat$Stage) # stages
# dim(pdat) # 51 patients
# 
# #-----------------------------------------------#
# #- Nanostring expression - normalize with Voom -#
# #-----------------------------------------------#
# 
# expr_dat <- fread('~/Documents/R/CD155_final_draft/data/nano_expression.csv') %>% 
#   data.frame(row.names=1) %>% na.omit
# vexpr <- limma::voom(expr_dat)
# vexpr <- data.matrix(vexpr)
# 
# vnames <- gsub("\\D", "", gsub('[0-9]$', '', colnames(vexpr))) %>% 
#   as.numeric %>% 
#   sprintf("%02d", .)
# 
# # take average by patient
# vexpr <- setkey(data.table(t(vexpr), pid=vnames)[, lapply(.SD, mean), by=pid], pid)
# pvr <- vexpr[, c('pid', 'PVR', 'CD274')]
# 
# #-----------------------------------------------#
# #----------------- IHC expression --------------#
# #-----------------------------------------------#
# 
# iroc <- fread('~/Documents/R/CD155_final_draft/data/iroc_ihc_database.csv')
# iroc$IROC_ID <- sprintf("%02d", iroc$IROC_ID)
# iroc <- setkey(iroc, IROC_ID)
# 
# iroc_expr <- iroc[pvr, ][ , c('IROC_ID', 'PVR', 'CD274', 'Total_CD8_epithelial', 'Total_CD8_stroma')]
# 
# # write out 
# write.csv(iroc_expr, file='data/IROC_expression_and_TIL_counts.csv', row.names=FALSE)
 
# --- note: manually removed two cases for consistency with CD155 IHC case set --- #

#-----------------------------------------------#
#------------ begin analysis here --------------#
#-----------------------------------------------#

# load data
 iroc_dat <- fread('data/IROC_expression_and_TIL_counts.csv')

# make plots - use total T cells on log10 scale
p1 <- ggplot(iroc_dat, aes(x=PVR, y=log10(Total_CD8_epithelial + Total_CD8_stroma + 1))) + 
  geom_point(color=rgb(0,.5,1), size=3) + geom_smooth(method='lm', fill=NA) +
  ylab(expression(log[10]~CD8+~(cells/mm^2))) +
  xlab(expression(log[2]~italic(PVR)~expression))

p2 <- ggplot(iroc_dat, aes(x=CD274, y=log10(Total_CD8_epithelial + Total_CD8_stroma + 1))) + 
  geom_point(color=rgb(0,.5,1), size=3) + geom_smooth(method='lm', fill=NA) +
  ylab(expression(log[10]~CD8+~(cells/mm^2))) +
  xlab(expression(log[2]~italic(CD*"274")~expression)) 

plot_grid(p2, p1)
dev.copy2pdf(file='~/Documents/R/CD155_final_draft/figures/IROC_CD155_PDL1_NanoString_July10.pdf', height=3, width=6)

# Get Pvals from correlation
with( iroc_dat, cor.test(PVR, log10(Total_CD8_epithelial + Total_CD8_stroma + 1), method='sp') )
with( iroc_dat, cor.test(CD274, log10(Total_CD8_epithelial + Total_CD8_stroma + 1), method='sp') )

#-----------------------------------------------#
#------------ load and process H scores --------#
#-----------------------------------------------#
# # read and merge H score data
# hscores <- fread('data/Hscores CD155 CD112 PDL1 compiled 10Jul19.csv')
# 
# hscores <- hscores[ ,c('Unique ID', 'Average CD155', 'Average cd112', 'Average pdl1') ]
# hscores$'Unique ID' <- substr(hscores$'Unique ID', 1, 3) %>%
#   as.numeric %>% 
#   sprintf("%02d", .)
# hscores <- setkey(hscores, 'Unique ID')
# 
# iroc_dat$IROC_ID <- as.character(iroc_dat$IROC_ID) 
# iroc_dat <- setkey(iroc_dat, 'IROC_ID')
# 
# hmerge <- hscores[iroc_dat, ]
# write.csv(hscores, file='data/CD155_and_CD112_Hscores.csv')

#------ plots ------#
hmerge <- read.csv(file='data/CD155_and_CD112_Hscores.csv')

p1 <- ggplot(hmerge, aes(x=`Average CD155`, y=log10(Total_CD8_epithelial + Total_CD8_stroma + 1))) + 
  geom_point(color=rgb(0,.5,1), size=3) + geom_smooth(method='lm', fill=NA) +
  ylab(expression(log[10]~CD8+~(cells/mm^2))) +
  xlab('CD155 H-score')

p2 <- ggplot(hmerge, aes(x=`Average pdl1`, y=log10(Total_CD8_epithelial + Total_CD8_stroma + 1))) + 
  geom_point(color=rgb(0,.5,1), size=3) + geom_smooth(method='lm', fill=NA) +
  ylab(expression(log[10]~CD8+~(cells/mm^2))) +
  xlab('PD-L1 H-score')

plot_grid(p2, p1)
dev.copy2pdf(file='figures/IROC_CD155_PDL1_IHC_July10.pdf', height=3, width=6)

#------ pvals ------#
with( hmerge, cor.test(`Average CD155`,log10(Total_CD8_epithelial + Total_CD8_stroma + 1), method='sp') )
with( hmerge, cor.test(`Average pdl1`,log10(Total_CD8_epithelial + Total_CD8_stroma + 1), method='sp') )
