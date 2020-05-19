# Analysis of outcome associations with PVR and PD-L1
# requires data from PanCanAtlas -- includes clinical, RNAseq, and Immunological data
# requires Ovarian TCGA subtypes -- provided in data folder

# updated for reviewer comments
library(data.table)
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(cowplot)
library(survival)
library(purrr)

# setwd
wd <- '~/Documents/R/CD155_final_draft/'
setwd(wd)

#----- read PanCan data - accessed Spring 2019 -----#
clins <- fread('data/pancan_clins.csv') # clinical data renamed
ovsubs <- fread('data/TCGA_OV_subtypes.csv') # ov subtypes
expr <- fread('data/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv') # RNAseq
pcimm <- fread('data/pancan_imm.csv') # immune

#----- convert and clean expression -----#
# drop duplicates - second occurence of geneid
expr <- data.frame(expr, row.names=1)
expr <- expr[!grepl('\\?', rownames(expr)), ]
gids <- gsub('\\|.*$', '', rownames(expr))
expr <- expr[!duplicated(gids), ]
rownames(expr) <- gsub('\\|.*$', '', rownames(expr)) # rename

# remove patients with duplicate samples
# crop to 12 digit codes
pats <- substr(colnames(expr), 1, 12)
dups <- pats[duplicated(pats)]
expr <- expr[ , !pats %in% dups]
colnames(expr) <- gsub('\\.', '-', substr(colnames(expr), 1, 12))

# align clinical and expression data
pat <- intersect(colnames(expr), clins$bcr_patient_barcode)
clins <- data.frame(clins, row.names = 2)
clins <- clins[pat, ]
expr <- expr[ , pat]

# subset genes of interest
# log transform
eset <- log2(t(expr[c('PVR', 'GZMA', 'PRF1', 'TIGIT', 'CD274', 'PDCD1', 'CD8A', 'PVRL2', 'CD276'), ]) + 1) %>% 
  data.frame

# add cytolytic score and cancer type
eset <- data.frame(eset, cyt=(eset$GZMA+eset$PRF1)/2, type=clins$type)

# loop-works for now to build df
# convert to lapply later
cors <- data.frame(row.names=levels(eset$type))
count <- 1
for(s in levels(eset$type)) {

  cors$pval[count] <- signif(cor.test(eset[eset$type==s, 'PVR'], 
                                  as.numeric(eset[eset$type==s, 'cyt']), 
                                  method='sp')$p.value, 2
  )
  cors$coef[count]<- signif(cor.test(eset[eset$type==s, 'PVR'], 
                                 as.numeric(eset[eset$type==s, 'cyt']), 
                                 method='sp')$estimate, 2
  )
  count <- count+1
}
cors$padj <- p.adjust(cors$pval, method='BH')
cors <- cors[order(cors$coef), ]
write.csv(cors,'figures/PanCancer_Cyt_PVR_Spearman_Correlations.csv')


#-------- Survival --------$
# overall survival across cancers
os <- Surv(clins$OS.time, clins$OS)
cmod <- coxph(os ~ eset$PVR + strata(clins$type))

# plot assumption - seems ok
cox.zph(cmod)
plot(cox.zph(cmod))

# function to get hazard ratio and confints
sfunc <- function(site, gene) {
  sos <- os[clins$type==site]
  spvr <- eset[clins$type==site, gene]
  cmod <- summary(coxph(sos ~ spvr))
  t(cmod$conf.int)
}

# loop for cancers
ss <- unique(clins$type)

# run for PVR
smods <- t(purrr::map_df(setNames(ss, ss), sfunc, gene='PVR'))
smods <- log(smods[rev(order(smods[,1])), -2])
colnames(smods) <- c('coef', 'lci', 'uci')
smods <- data.frame(smods, rank=rev(rank(smods[,'coef'])), 
                    type=factor(rownames(smods), levels=rownames(smods)))

# then PDL1
smodsPDL1 <- t(purrr::map_df(setNames(ss, ss), sfunc, gene='CD274'))
smodsPDL1 <- log(smodsPDL1[rev(order(smodsPDL1[,1])), -2])
colnames(smodsPDL1) <- c('coef', 'lci', 'uci')
smodsPDL1 <- data.frame(smodsPDL1, rank=rev(rank(smodsPDL1[,'coef'])), 
                    type=factor(rownames(smodsPDL1), levels=rownames(smodsPDL1)))

# merge PVR and CD274 hazards
smerge <- merge(smods, smodsPDL1, by.x='type', by.y='type')

# test association
cor.test(smerge$coef.x, smerge$coef.y, method='sp')

# plot 2-D hazards
ggplot(smerge, aes(x=coef.x, y=coef.y)) + 
  geom_hline(yintercept = 0, linetype = 'dashed', color=rgb(.5,.5,.5,0.5)) +
  geom_vline(xintercept = 0, linetype = 'dashed', color=rgb(.5,.5,.5,0.5)) +
  geom_errorbar(aes(ymin=lci.y, ymax=uci.y), col=rgb(0,.5,1,0.2)) +
  geom_errorbarh(aes(xmin=lci.x, xmax=uci.x), col=rgb(0,.5,1,0.2)) +
  geom_point() +
  xlab(expression(italic(PVR)~log~hazard~ratio)) +
  ylab(expression(italic(CD274)~log~hazard~ratio)) +
  annotate('text', label=expression(rho==0.04), x=2.5, y=1.25)
  

dev.copy2pdf(file='figures/TCGA_PDL1_CD155_hazards.pdf')

#----------------------------------------------------------#
#-------- Analysis of Residual Expression in OV -----------#
#----------------------------------------------------------#

# define new set with relevant markers from Danaher et al.
eset <- log2(t(expr[c('PVR', 'GZMA', 'PRF1', 'TIGIT', 
                      'CD274', 'PDCD1', 'CD8A','CD8B', 
                      'PVRL2', 'CD276', 'XCL1', 'XCL2', 'NCR1',
                      'KIR2DL3', 'KIR3DL1', 'KIR3DL2', 'IL21R'),  
                    ]) + 1) %>% data.frame

# restrict to subtyped OV cases
ovsubs <- data.frame(ovsubs, row.names=1)
ovs <- intersect(rownames(ovsubs), rownames(eset))

ovexp <- eset[ovs, ]
ovsubs <- ovsubs[ovs, ]

# update with Danaher CD8 score
# use geo mean of CD8A and CD8B
tigit_res <- resid(lm(TIGIT ~ I((CD8A+CD8B)/2), data=ovexp))
pd1_res <- resid(lm(PDCD1 ~ I((CD8A+CD8B)/2), data=ovexp))

# compute CD8 and 2 NK signatures
cd8 <- with(ovexp, (CD8A+CD8B)/2)
nk1 <- with(ovexp, (XCL1 + XCL2 + NCR1)/3)
nk2 <- with(ovexp, (KIR2DL3 + KIR3DL1 + KIR3DL2 + IL21R)/4)

# test combined signal of co-expression controlling for NKs and T cells
# association between PD1 and TIGIT remains strong
mod <- lm(ovexp$TIGIT ~ ovexp$PDCD1 + nk2 + nk1 + cd8)

cor.test(tigit_res, pd1_res, method='sp')

# Make Figure
sub_pal <- rainbow(4, s=0.7, v=0.6)
ggdat <- data.frame(tigit_res, pd1_res, subtype=ovsubs$SUBTYPE)
ggplot(ggdat, aes(x=pd1_res, y=tigit_res, color=subtype)) + 
  geom_point(size=3) +
  scale_color_manual('', values=sub_pal) + 
  xlab(expression(italic(PDCD*"1")~residual~expression)) +
  ylab(expression(italic(TIGIT)~residual~expression)) +
  annotate('text', 3, 3, label=expression(P<10^-12))

dev.copy2pdf(file='figures/TIGIT_PD1_TCGA_residuals.pdf', width=6, height=4)


