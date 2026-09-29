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
source("functions_davidcomments.R")

read.xlsx("data/LungFunctionMasterDataSet.xlsx",sheet=2) -> animal_data
read.xlsx("data/LungFunctionMasterDataSet.xlsx",sheet=1,fillMergedCells = TRUE) -> LF_data
###
LF_data %>% filter(Parameter=="DCP_EF50") %>% pivot_longer(cols=4:39,names_to = "Animal.ID") %>%
  group_by(`Animal.ID`,ZT) %>%
  #filter(Mch_Conc!=0) %>%
  #mutate(value=log10(value)) %>%
  merge(animal_data,by="Animal.ID")-> LF_data

LF_data %>% rename(Sample=Animal.ID,Mch_conc=Mch_Conc,Value=value) %>%
  filter(Parameter=="DCP_EF50") %>%
  dplyr::select(-Parameter,-Cull_time) %>%
  dplyr::mutate(Genotype=ifelse(Genotype=="HET","KO","WT"))-> LF_data

# mann whitney u on pairwise zts ----------------------------------

mw_results<-dr_MWU_pairwise(LF_data)

# dose response model -----------------------------------------------------

# param_formodel <- dr_fit(LF_data) #uncomment these lines to rerun dr curve fits - warning takes some time
# save(param_formodel,file = "data/drcmodelparams_ef50.RData")
load("data/drcmodelparams_ef50.RData")

anova_pvals_upper <- dr_anova(param_formodel,"Upper","log10(params+1) ~ ZT * Treatment")

anova_pvals_slope <- dr_anova(param_formodel,"Slope","params ~ ZT * Treatment")

anova_pvals_slope_geno <- dr_anova_gen(param_formodel,"Slope","log10(-params) ~ ZT*Genotype")

anova_pvals_slope %>% filter(WT<0.05|KO<0.05)
anova_pvals_slope_geno %>% filter(PBS<0.05|HDM<0.05)

# plotting dose response curve --------------------------------------------


p1<-dr_plot(LF_data,anova_pvals_slope,mw_results%>% mutate(p.adj=1),c(1,2.5),
            y_lab="Mean EF50 (ml.sec<sup>-1</sup>)",
            x_lab=expression("Methacholine Concentration (mg.mL"^"-1"*")"))
p1
png("plots_v4/ef50_meth_dose_response.png", width = 3000, height = 1500, res = 300)  # adjust size/res as needed
grid.draw(
  p1
)
#grid.text( expression("Methacholine Concentration (mg.mL"^"-1"*")"), y = unit(0.015, "npc"), gp = gpar(fontsize = 12))
dev.off()

p1<-dr_plot(LF_data,anova_pvals_slope,mw_results%>% mutate(p.adj=1),c(0.5,3),
            y_lab="Mean EF50 (ml.sec<sup>-1</sup>)",
            x_lab=expression("Methacholine Concentration (mg.mL"^"-1"*")"),errorbar = T)
p1
png("plots_v4/ef50_meth_dose_response_eb.png", width = 3000, height = 1500, res = 300)  # adjust size/res as needed
grid.draw(
  p1
)
#grid.text( expression("Methacholine Concentration (mg.mL"^"-1"*")"), y = unit(0.015, "npc"), gp = gpar(fontsize = 12))
dev.off()
# AUC sinusoidal analysis -----------------------------------------------------

anova_box(LF_data,"AUC") -> anova_auc

plot_box(LF_data,"AUC",y_lim=c(-11,20),y_lab="AUC of EF50 (ml.sec<sup>-1</sup>)",anova_auc) -> analysis_out

p_auc <- analysis_out
#ggsave(p_auc,filename="plots_v2/flex_AUC_bar_v2.png",width=10,height=5)

png("plots_v4/EF50_AUC_bar_v3.png", width = 2500, height = 1250, res = 300)  # adjust size/res as needed
grid.draw(
  p_auc
)
dev.off()


# Max sinusoidal analysis -----------------------------------------------------
anova_box(LF_data,"Max") -> anova_auc

plot_box(LF_data,"Max",y_lim=c(-.1,.9),y_lab="Max EF50 (ml.sec<sup>-1</sup>)",anova_auc) -> analysis_out

p_auc <- analysis_out
#ggsave(p_auc,filename="plots_v2/flex_AUC_bar_v2.png",width=10,height=5)

png("plots_v4/EF50_max_bar_v3.png", width = 2500, height = 1250, res = 300)  # adjust size/res as needed
grid.draw(
  p_auc
)
dev.off()

