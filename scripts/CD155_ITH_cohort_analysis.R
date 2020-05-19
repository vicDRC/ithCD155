# CD155 Analysis: using Zhang et al NanoString Data

# dependencies
library(data.table)
library(dplyr)
library(lme4)
library(ggplot2)
library(cowplot)
library(ggtext)
library(ggpubr)
library(emmeans)

# set cowplot theme
theme_set(theme_cowplot())

# set dir
wd <- '~/Documents/R/CD155_final_draft/'
setwd(wd)

# define simple function for plots
# parse emmeans output for contrast plotting
getSigs <- function(emm, 
                    ypos, 
                    stagger, 
                    hide.nsig=FALSE) {
  
  pdf <- data.frame(emm[[2]])
  rtab <- read.table(text=as.character(pdf$contrast))
  outdf <- data.frame(group1=as.character(rtab[, 1]),
                      group2=as.character(rtab[, 3]),
                      p.value = signif(pdf$p.value, 2),
                      y.position = ypos + seq(0, 
                                              (nrow(pdf)-1) * stagger, 
                                              by = stagger),
                      'p.stars' = cut(pdf$p.value, 
                                      breaks = c(0, 0.001, 0.01, 0.05, 1), 
                                      labels = c('***', '**', '*', 'n.s') )
                     )
  
  if(hide.nsig) outdf <- outdf[outdf$p.value < 0.05, ]
  
  return(outdf)
  
}

#---------------------- datasets ----------------------#
# these are taken from Supp.of Zhang et al. 2018, Cell # 
#------------------------------------------------------#

# IHC cell counts
cells <- fread('data/ITH_cell_densities.csv') %>% data.frame(row.names=1)

# nanostring data
nano <- fread('data/nanostring_data_ids_mapped.tsv') %>% data.frame(row.names=2) %>% .[, -c(1,2)]
colnames(nano) <- gsub('^X', '', colnames(nano)) 

# addition for reviewer comment - site mapping
maps_red <- fread('data/site_map.csv') %>% data.frame(., row.names='sample_id')

# clean up
cells[cells=='#N/A'] <- NA
cells_na <- na.omit(cells[, 2:13])

#get intersecting patients
nano_com <- Reduce(intersect, list(colnames(nano), rownames(cells)))
nano_int <- intersect(nano_com, rownames(maps_red))

#extract patient number
pat_a <- as.factor(gsub('_.*$', '', nano_com))

# mixed models
# use CD8 as response var
mod <- lmer(as.numeric(log2(cells[nano_com, 'Overall.CD8..density']+1)) ~ t(nano['PVR', nano_com]) + (1|pat_a))
car::Anova(mod)

# add in sample site
# borderline effect of site on CD8--reacpitulates what we expected
moddf <- data.frame(logCD8=log2(cells[nano_com, 'Overall.CD8..density']+1),
                    PVR=t(nano['PVR', nano_com]),
                    tissue=factor(maps_red[nano_com, 'tissue'], levels=c('Omentum', 'Ovary', 'Other')),
                    pat=pat_a
                  )
mod1 <- lmer(logCD8 ~ PVR * tissue + (1|pat), data=moddf)
car::Anova(mod1)

# look at PVR as a function of site alone: appears significant
mod2 <- lmer(PVR ~ tissue + (1|pat), data=moddf)
car::Anova(mod2)
em <- emmeans(mod2, list(pairwise ~ tissue), adjust = "tukey")

# build manual sig df for plotting
plotsigs <- data.frame(group1=c('Omentum', 'Omentum', 'Ovary'),
                       group2=c('Ovary', 'Other', "Other"),
                       p=c('***', '**', 'ns'),
                       y.position=c(11.5, 11.75, 12))


# plot PVR ~ site
ggplot(moddf, aes(x=tissue, y=PVR)) + geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width=0.1, alpha=0.7, size=2.5, aes(color=as.factor(pat))) +
  scale_color_discrete('Patient', rainbow(20)) +
  xlab('') + ylab("*PVR* (log2)") + theme(axis.title.y = element_markdown()) +
  stat_pvalue_manual(plotsigs, tip.length = 0)
dev.copy2pdf(file='figures/PVR_across_disease_sites_March2020.pdf', height=4, width=4)


# overall plotting df for other data
ggdat <- na.omit(data.frame(cd8=log10(as.numeric(cells[nano_com, 'Overall.CD8..density'])+1), 
                    pvr=t(nano['PVR', nano_com]), patient=as.factor(pat_a), 
                    til_cluster=cells[nano_com, 'til_cluster'], 
                    pdl1=t(nano['CD274', nano_com])))
ggdat$til_cluster <- factor(ggdat$til_cluster, levels=c('N-TIL', 'S-TIL', 'ES-TIL'))


# pvalues
mod <- lmer(PVR ~ til_cluster + (1|patient), data=ggdat)
emPVR <- emmeans(mod, list(pairwise ~ til_cluster), adjust = "tukey")
mod <- lmer(CD274 ~ til_cluster + (1|patient), data=ggdat)
em274 <- emmeans(mod, list(pairwise ~ til_cluster), adjust = "tukey")

# plots #
p1 <- ggplot(ggdat, aes(x=PVR, y=cd8, color=patient)) + geom_point(size=3) + 
  geom_smooth(method='lm',  fill=NA, size=.5, alpha=0.3) + 
  ylim(0,3.5) + xlab(expression(paste('log'[2],' ',italic(PVR),' expression'))) + 
  ylab(expression(paste('log'[10],' CD8+ (cells/mm'^2,')'))) +
  theme(legend.position = "none")  
p2 <- ggplot(ggdat, aes(x=CD274, y=cd8, color=patient)) + geom_point(size=3) + 
  geom_smooth(method='lm',  fill=NA, size=.5, alpha=0.3) + 
  ylim(0,3.5) + xlab(expression(paste('log'[2],' ',italic(CD*"274"),' expression'))) + 
  ylab(expression(paste('log'[10],' CD8+ (cells/mm'^2,')'))) + theme(legend.position = "none") 


p3 <- ggplot(ggdat, aes(y=PVR, x=til_cluster)) + geom_boxplot() + 
  ylab(expression(paste('log'[2],' ',italic(PVR),' expression'))) +
  xlab('') + 
  stat_pvalue_manual(getSigs(emPVR, ypos=11.5, stagger=0.25, hide.nsig = TRUE), label = 'p.stars',  tip.length= 0)
p4 <- ggplot(ggdat, aes(y=CD274, x=til_cluster)) + geom_boxplot() + 
  ylab(expression(paste('log'[2],' ',italic(CD*"274"),' expression'))) +
  xlab('') +
  stat_pvalue_manual(getSigs(em274, ypos=9, stagger=0.35, hide.nsig = TRUE), label = 'p.stars',  tip.length= 0)

#------------ write plots ------------#
dev.off()
plot_grid(p2, p1, align='hv', axis='tblr')
dev.copy2pdf(file='~/Documents/R/CD155_final_draft/figures/CD8_vs_PVR_ITH_cohort_with_lines.pdf',
             width=6, height=3)

plot_grid(p4, p3, align='hv', axis='tblr')
dev.copy2pdf(file='~/Documents/R/CD155_final_draft/figures/PVR_and_CD274_boxplot_ITH_cohort_March2020.pdf', 
             width=6, height=3)




