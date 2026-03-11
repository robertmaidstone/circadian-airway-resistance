library(tidyverse)
library(openxlsx)
library(grid)
library(patchwork)
library(ggtext)
library(drc)
library(lme4)

# load data ---------------------------------------------------------------
setwd("~/Amlan/AirwayResistance/")
source("functions.R")
read.xlsx("data/AC_FVData_Complete.xlsx",sheet=5,fillMergedCells = TRUE) -> LF_data

LF_data %>% filter(Mch_conc!=0) -> LF_data

# anovas replace with tests as discussed ----------------------------------

model <- lmer(Value ~ ZT + Mch_conc + (1 | Sample), data = LF_data%>% filter(Treatment=="HDM",Genotype=="WT"))
anova(model)


lm(Value~ZT+Mch_conc,data=LF_data%>% filter(Treatment=="HDM",Genotype=="WT")) %>% anova
lm(Value~Treatment*ZT,data=LF_data%>% filter(Genotype=="KO")) %>% anova

lm(Value~ZT+Mch_conc,data=LF_data%>% filter(Genotype=="WT",Treatment=="HDM")) %>% anova
lm(Value~ZT,data=LF_data%>% filter(Genotype=="KO",Treatment=="HDM")) %>% anova

# dose response model -----------------------------------------------------

# param_formodel <- dr_fit(LF_data) #uncomment these lines to rerun dr curve fits - warning takes some time
# save(param_formodel,file = "data/drcmodelparams_flex.RData")
load("data/drcmodelparams_flex.RData")

anova_pvals_upper <- dr_anova(param_formodel,"Upper","log10(params) ~ ZT * Treatment")

anova_pvals_slope <- dr_anova(param_formodel,"Slope","log10(-params) ~ ZT * Treatment")

# plotting dose response curve --------------------------------------------

dr_plot(LF_data,anova_pvals_slope,anova_pvals_upper,c(0.5,6),y_lab="Maximum Airway Resistance R<sub>rs</sub>(cm.H<sub>2</sub>O.s.ml<sup>-1</sup>)")

png("plots/flex_meth_dose_response.png", width = 2400, height = 1500, res = 300)  # adjust size/res as needed
grid.draw(
  dr_plot(LF_data,anova_pvals_slope,anova_pvals_upper,c(0.5,6),y_lab="Maximum Airway Resistance R<sub>rs</sub>(cm.H<sub>2</sub>O.s.ml<sup>-1</sup>)")
)
grid.text( expression("Methacholine Concentration (mg.mL"^"-1"*")"), y = unit(0.03, "npc"), gp = gpar(fontsize = 10))
dev.off()


# AUC sinusoidal analysis -----------------------------------------------------

rhy_plot(LF_data,"AUC",y_lim=c(-6,50)) -> analysis_out

p_auc <- analysis_out$combined
ggsave(p_auc,filename="plots/flex_AUC.png",width=8,height=5)

# Max sinusoidal analysis -----------------------------------------------------

rhy_plot(LF_data,"Max",y_lim=c(-.1,1.3)) -> analysis_out

p_max <- analysis_out$combined
ggsave(p_max,filename="plots/flex_max.png",width=8,height=5)
