library(tidyverse)
library(openxlsx)
library(grid)
library(patchwork)
library(ggtext)
library(drc)
library(lme4)
library(broom)

# load data ---------------------------------------------------------------
setwd("~/Amlan/AirwayResistance/")
source("functions.R")

read.xlsx("data/LungFunctionMasterDataSet.xlsx",sheet=2) -> animal_data
read.xlsx("data/LungFunctionMasterDataSet.xlsx",sheet=1,fillMergedCells = TRUE) -> LF_data
###
LF_data %>% filter(Parameter=="DCP_sRAW") %>% pivot_longer(cols=4:39,names_to = "Animal.ID") %>%
  group_by(`Animal.ID`,ZT) %>%
  filter(Mch_Conc!=0) %>%
  merge(animal_data,by="Animal.ID")-> LF_data

LF_data %>% rename(Sample=Animal.ID,Mch_conc=Mch_Conc,Value=value) %>%
  dplyr::select(-Parameter,-Cull_time) %>%
  dplyr::mutate(Genotype=ifelse(Genotype=="HET","KO","WT"))-> LF_data

# mann whitney u on pairwise zts ----------------------------------

mw_results<-dr_MWU_pairwise(LF_data)

# dose response model -----------------------------------------------------

# param_formodel <- dr_fit(LF_data) #uncomment these lines to rerun dr curve fits - warning takes some time
# save(param_formodel,file = "data/drcmodelparams_sRAW.RData")
load("data/drcmodelparams_sRAW.RData")

anova_pvals_upper <- dr_anova(param_formodel,"Upper","log10(params+1) ~ ZT * Treatment")

anova_pvals_slope <- dr_anova(param_formodel,"Slope","params ~ ZT * Treatment")

anova_pvals_slope_geno <- dr_anova_gen(param_formodel,"Slope","log10(-params) ~ ZT*Genotype")

anova_pvals_slope %>% filter(WT<0.05|KO<0.05)
anova_pvals_slope_geno %>% filter(PBS<0.05|HDM<0.05)

# plotting dose response curve --------------------------------------------

p1 <- dr_plot(LF_data,anova_pvals_slope,mw_results,c(2,7.5),y_lab="Mean sRAW (cm.H<sub>2</sub>O.sec)")
p1

png("plots/sRAW_meth_dose_response.png", width = 3000, height = 1500, res = 300)  # adjust size/res as needed
grid.draw(
    patchworkGrob(p1)
)
grid.text(
  "HDM slope\nchange by\ngenotype *",
  x = unit(0.89, "npc"),   # horizontal position
  y = unit(0.23, "npc"),   # vertical position
  just = "left",
  gp = gpar(fontsize = 12)
)
grid.text( expression("Methacholine Concentration (mg.mL"^"-1"*")"), y = unit(0.015, "npc"), gp = gpar(fontsize = 12))
dev.off()


# AUC sinusoidal analysis -----------------------------------------------------

rhy_plot(LF_data,"AUC",y_lim=c(5,22),y_lab="AUC of sRAW (cm.H<sub>2</sub>O.sec)") -> analysis_out

analysis_out$combined
p_auc <- analysis_out$combined
ggsave(p_auc,filename="plots/sRAW_AUC.png",width=10,height=5)

# Max sinusoidal analysis -----------------------------------------------------

rhy_plot(LF_data,"Max",y_lim=c(0.4,1.2),y_lab="Max sRAW (cm.H<sub>2</sub>O.sec)") -> analysis_out

analysis_out$combined
p_max <- analysis_out$combined
ggsave(p_max,filename="plots/sRAW_max.png",width=10,height=5)
