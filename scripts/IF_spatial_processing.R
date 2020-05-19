# This script will run hotspot analysis for a set QuPath detection files
# These are quite large: summary file provided in data

# dependencies
library(data.table)
library(spdep)
library(spatstat)
library(dplyr)
library(ggplot2)
library(cowplot)

theme_set(theme_cowplot())

## run hotspot analysis -- first source functions
source('~/Documents/R/ithCD155/scripts/Hotspot_Functions.R')

# # get file list of detections
# qp_files <- list.files('data/IF_data_restain', full.names = TRUE, pattern = 'detection')

# quick wrapper
# hwrap <- function(x) run_hotspot(generate_ppp(qp_files[x]))

# get Gi* and Fcs
# params as set in function scripts - require 5 of each cell type

# runs slowly--save and source
#test_gstats <- lapply(seq_along(qp_files), plyr::failwith(NULL, hwrap))
#saveRDS(test_gstats, file='~/Documents/R/ithCD155/data/Gstats_for_IROC_April30.rds')

#----------------------------------------------------------------#
#----------------------------------------------------------------#
# NOTE: Update April 2020: updated hotspot generating function for
# compatibility/robustness with new versions of dependencies. 
# Masking strategy appeared to have broken 
# (due to an update in spatstat::Smooth causing underflow?). 
# Set sigma for bandwidth fn to bw.ppl with padding (1.5).
# Results in small differences from previous analysis: see below.
#----------------------------------------------------------------#
#----------------------------------------------------------------#

# load g statistics
trds <- readRDS(file='data/Gstats_for_IROC_April30.rds')

# get Fc stats
fset <- lapply(trds, getFcs)

# relevant vals for PD1 Fc
pdl1_fcs <- lapply(seq_along(fset), function(x) fset[[x]]$f_pdl1_pd1$c)
cd155_fcs <- lapply(seq_along(fset), function(x) fset[[x]]$f_cd155_pd1$c)

# relevant vals for CD8 Fc
pdl1_cd8_fcs <- lapply(seq_along(fset), function(x) fset[[x]]$f_pdl1_cd8$c)
cd155_cd8_fcs <- lapply(seq_along(fset), function(x) fset[[x]]$f_cd155_cd8$c)

# map core IDs to patient IDs
# get TMA map and clean for matching
tmap <- read.csv('data/ID position List-Table 1.csv', 
                  stringsAsFactors = FALSE)
tmap$TMA.Location <- gsub('\\[|\\]', '', tmap$TMA.Location)
qp_map <- sapply(strsplit(qp_files, '\\[|\\]'), '[', 2)

pid <- tmap[match(qp_map, tmap$TMA.Location), 'Patient.ID']
pid <- substr(pid, 1, 3)

# compile and extract relevant stats
fcdf <- data.table(pdl1=unlist(pdl1_fcs),
                   cd155=unlist(cd155_fcs),
                   pdl1_cd8=unlist(pdl1_cd8_fcs),
                   cd155_cd8=unlist(cd155_cd8_fcs),
                   pid)
fcdf <- na.omit(fcdf[, lapply(.SD, mean, na.rm=TRUE), by=pid])
fcdf <- fcdf[pid!='Ton']

# write mapped stats
write.csv(fcdf, 'data/Mapped_Fc_Stats.csv', row.names = FALSE)

fcdf <- fread('data/Mapped_Fc_Stats.csv')
fcdf$pid <- as.factor(fcdf$pid)

# compile separately for PD-1 and CD8 overlaps
fcdfm <- melt(fcdf[, c('pid', 'pdl1', 'cd155')])
fcdfm_cd8 <- melt(fcdf[, c('pid', 'pdl1_cd8', 'cd155_cd8')])

# Generate plots
g1 <- ggplot(fcdfm, aes(x=variable, y=value, group=pid)) + 
  geom_line(color='gray') + 
  geom_point() + xlab('') +
  ylab(expression(F[c]))  + ggtitle('PD1') +
  scale_x_discrete(labels= c('PD-L1', 'CD155')) + 
  theme(axis.text.x = element_text(angle = 90)) +
  theme_cowplot()

g2 <- ggplot(fcdfm_cd8, aes(x=variable, y=value, group=pid)) + 
  geom_line(color='gray') + 
  geom_point() + xlab('') +
  ylab(expression(F[c]))  + ggtitle('CD8') +
  scale_x_discrete(labels= c('PD-L1', 'CD155')) + 
  theme(axis.text.x = element_text(angle = 90)) +
  theme_cowplot()

plot_grid(g1, g2)

# export
dev.copy2pdf(file='figures/Hotspot_Figure_April29.pdf', height=4, width=4)

# calculate P values -- paired T-tests
t.test(fcdf$pdl1, fcdf$cd155, paired = TRUE)
t.test(fcdf$pdl1_cd8, fcdf$cd155_cd8, paired = TRUE)



