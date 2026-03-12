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

read.xlsx("data/LungFunctionMasterDataSet.xlsx",sheet=2) -> animal_data
read.xlsx("data/LungFunctionMasterDataSet.xlsx",sheet=1,fillMergedCells = TRUE) -> LF_data
###
LF_data %>% filter(Parameter=="DCP_sGAW") %>% pivot_longer(cols=4:39,names_to = "Animal.ID") %>%
  group_by(`Animal.ID`,ZT) %>%
  filter(Mch_Conc!=0) %>%
  merge(animal_data,by="Animal.ID")-> LF_data

LF_data %>% rename(Sample=Animal.ID,Mch_conc=Mch_Conc,Value=value) %>%
  dplyr::select(-Parameter,-Cull_time) %>%
  dplyr::mutate(Genotype=ifelse(Genotype=="HET","KO","WT"))-> LF_data

# anovas replace with tests as discussed ----------------------------------

model <- lmer(Value ~ ZT + Mch_conc + (1 | Sample), data = LF_data%>% filter(Treatment=="HDM",Genotype=="WT"))
anova(model)


lm(Value~ZT+Mch_conc,data=LF_data%>% filter(Treatment=="HDM",Genotype=="WT")) %>% anova
lm(Value~Treatment*ZT,data=LF_data%>% filter(Genotype=="KO")) %>% anova

lm(Value~ZT+Mch_conc,data=LF_data%>% filter(Genotype=="WT",Treatment=="HDM")) %>% anova
lm(Value~ZT,data=LF_data%>% filter(Genotype=="KO",Treatment=="HDM")) %>% anova

# dose response model -----------------------------------------------------
param_formodel <- dr_fit_sep(LF_data) #uncomment these lines to rerun dr curve fits - warning takes some time
save(param_formodel,file = "data/drcmodelparams_sGAW.RData")
load("data/drcmodelparams_sGAW.RData")

anova_pvals_upper <- dr_anova(param_formodel,"Upper","log10(params+1) ~ ZT * Treatment")

anova_pvals_slope <- dr_anova(param_formodel,"Slope","params ~ ZT * Treatment")

# plotting dose response curve --------------------------------------------

dr_plot(LF_data,anova_pvals_slope,anova_pvals_upper,c(-.1,.5),y_lab="Median EF50 (ml.sec<sup>-1</sup>)")

png("plots/ef50_meth_dose_response.png", width = 2400, height = 1500, res = 300)  # adjust size/res as needed
grid.draw(
  dr_plot(LF_data,anova_pvals_slope,anova_pvals_upper,c(1,2.75),y_lab="Median EF50 (ml.sec<sup>-1</sup>)")
)
grid.text( expression("Methacholine Concentration (mg.mL"^"-1"*")"), y = unit(0.03, "npc"), gp = gpar(fontsize = 10))
dev.off()


# AUC sinusoidal analysis -----------------------------------------------------

rhy_plot(LF_data,"AUC",y_lim=c(-17.5,0)) -> analysis_out

analysis_out$combined
p_auc <- analysis_out$combined
ggsave(p_auc,filename="plots/sGAW_AUC.png",width=8,height=5)

# Max sinusoidal analysis -----------------------------------------------------

rhy_plot(LF_data,"Min",y_lim=c(-1,-.1)) -> analysis_out

analysis_out$combined
p_max <- analysis_out[["plot_pbs"]] + analysis_out[["plot_hdm"]]
ggsave(p_max,filename="plots/sGAW_min.png",width=8,height=5)
