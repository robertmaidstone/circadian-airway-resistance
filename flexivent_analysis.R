library(tidyverse)
library(openxlsx)
library(grid)
library(patchwork)
library(ggtext)
library(drc)
library(lme4)

library(rstatix)
library(purrr)

# load data ---------------------------------------------------------------
setwd("~/Amlan/AirwayResistance/")
source("functions.R")
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

# plotting dose response curve --------------------------------------------

p1<-dr_plot(LF_data,anova_pvals_slope,mw_results,c(0.7,5.7),y_lab="Maximum Airway Resistance R<sub>rs</sub>(cm.H<sub>2</sub>O.s.ml<sup>-1</sup>)")
p1

png("plots/flex_meth_dose_response.png", width = 2400, height = 1500, res = 300)  # adjust size/res as needed
grid.draw(
  p1
)
grid.text( expression("Methacholine Concentration (mg.mL"^"-1"*")"), y = unit(0.03, "npc"), gp = gpar(fontsize = 10))
dev.off()


# AUC sinusoidal analysis -----------------------------------------------------

rhy_plot(LF_data,"AUC",y_lim=c(-6,50),y_lab="AUC Airway Resistance R<sub>rs</sub>(cm.H<sub>2</sub>O.s.ml<sup>-1</sup>)") -> analysis_out

p_auc <- analysis_out$combined
ggsave(p_auc,filename="plots/flex_AUC.png",width=8,height=5)

# Max sinusoidal analysis -----------------------------------------------------

rhy_plot(LF_data,"Max",y_lim=c(-.1,1.3),y_lab="Max Airway Resistance R<sub>rs</sub>(cm.H<sub>2</sub>O.s.ml<sup>-1</sup>)") -> analysis_out

p_max <- analysis_out$combined
ggsave(p_max,filename="plots/flex_max.png",width=8,height=5)
