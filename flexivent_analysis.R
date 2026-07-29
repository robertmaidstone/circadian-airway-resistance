library(tidyverse)
library(openxlsx)
library(grid)
library(patchwork)
library(ggtext)
library(drc)
library(lme4)
library(broom)
library(rstatix)
library(purrr)

# load data ---------------------------------------------------------------
setwd("~/Amlan/AirwayResistance/")
source("functions_davidcomments.R")
read.xlsx("data/AC_FVData_Complete.xlsx",sheet=5,fillMergedCells = TRUE) -> LF_data

LF_data %>% filter(Mch_conc!=0) -> LF_data

# mann whitney u on pairwise zts ----------------------------------

mw_results<-dr_MWU_pairwise(LF_data)

# dose response model -----------------------------------------------------

# param_formodel <- dr_fit(LF_data) #uncomment these lines to rerun dr curve fits - warning takes some time
# save(param_formodel,file = "data/drcmodelparams_flex.RData")
load("data/drcmodelparams_flex.RData")

anova_pvals_upper <- dr_anova(param_formodel,"Upper","log10(params) ~ ZT * Treatment") # not used

anova_pvals_slope <- dr_anova(param_formodel,"Slope","log10(-params) ~ ZT * Treatment")

anova_pvals_slope_geno <- dr_anova_gen(param_formodel,"Slope","log10(-params) ~ ZT*Genotype")

anova_pvals_slope %>% filter(WT<0.05|KO<0.05)
anova_pvals_slope_geno %>% filter(PBS<0.05|HDM<0.05)

# plotting dose response curve --------------------------------------------

p1<-dr_plot(LF_data,anova_pvals_slope,mw_results,c(0.7,5.7),
            y_lab="Mean Airway Resistance R<sub>rs</sub>(cm.H<sub>2</sub>O.s.ml<sup>-1</sup>)",
            x_lab=expression("Methacholine Concentration (mg.mL"^"-1"*")"))
p1

png("plots_v2/flex_meth_dose_response.png", width = 3000, height = 1500, res = 300)  # adjust size/res as needed
grid.draw(
  p1
)
#grid.text( expression("Methacholine Concentration (mg.mL"^"-1"*")"), y = unit(0.03, "npc"), gp = gpar(fontsize = 12))
dev.off()


# AUC sinusoidal analysis -----------------------------------------------------

rhy_plot_bar(LF_data,"AUC",y_lim=c(-6,50),y_lab="AUC of Airway Resistance R<sub>rs</sub>(cm.H<sub>2</sub>O.s.ml<sup>-1</sup>)") -> analysis_out

p_auc <- analysis_out$combined
#ggsave(p_auc,filename="plots_v2/flex_AUC_bar_v2.png",width=10,height=5)

png("plots_v2/flex_AUC_bar_v3.png", width = 3000, height = 1500, res = 300)  # adjust size/res as needed
grid.draw(
  p_auc
)
grid.text( "WT", y = unit(0.03, "npc"),x = unit(0.175, "npc"), gp = gpar(fontsize = 12))
grid.text( "CCSP-Reverba KO", y = unit(0.03, "npc"),x = unit(0.4, "npc"), gp = gpar(fontsize = 12))
grid.text( "WT", y = unit(0.03, "npc"),x = unit(0.66, "npc"), gp = gpar(fontsize = 12))
grid.text( "CCSP-Reverba KO", y = unit(0.03, "npc"),x = unit(0.875, "npc"), gp = gpar(fontsize = 12))
dev.off()


# Max sinusoidal analysis -----------------------------------------------------

rhy_plot(LF_data,"Max",y_lim=c(-.1,1.3),y_lab="Max Airway Resistance R<sub>rs</sub>(cm.H<sub>2</sub>O.s.ml<sup>-1</sup>)") -> analysis_out

p_max <- analysis_out$combined
#ggsave(p_max,filename="plots_v2/flex_max.png",width=10,height=5)

png("plots_v2/flex_max.png", width = 3000, height = 1500, res = 300)  # adjust size/res as needed
grid.draw(
  p_max
)
grid.text( "WT", y = unit(0.03, "npc"),x = unit(0.175, "npc"), gp = gpar(fontsize = 12))
grid.text( "CCSP-Reverba KO", y = unit(0.03, "npc"),x = unit(0.4, "npc"), gp = gpar(fontsize = 12))
grid.text( "WT", y = unit(0.03, "npc"),x = unit(0.66, "npc"), gp = gpar(fontsize = 12))
grid.text( "CCSP-Reverba KO", y = unit(0.03, "npc"),x = unit(0.875, "npc"), gp = gpar(fontsize = 12))
dev.off()