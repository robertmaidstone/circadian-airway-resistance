library(tidyverse)
library(openxlsx)
library(grid)
library(patchwork)
library(ggtext)
library(drc)
library(lme4)

library(rstatix)
library(purrr)

#par_val<-"Cull_FV_Rrs"
par_val<-"Cull_FV_RN"
#par_val<-"Cull_FV_Rrs"

# load data ---------------------------------------------------------------
read.xlsx("data/LungFunctionMasterDataSet.xlsx",sheet=2) -> animal_data # need this as some animals not in LFmaster
read.xlsx("data/LFmaster.xlsx",sheet=2) -> animal_data_2
read.xlsx("data/LFmaster.xlsx",sheet=1,fillMergedCells = TRUE) -> LF_data2

rbind(animal_data %>% dplyr::select(Animal.ID,Treatment,Genotype) %>% unique,animal_data_2 %>% dplyr::select(Animal.ID=Sample,Treatment,Genotype) %>% unique) %>% unique() %>%
  dplyr::mutate(Genotype=ifelse(Genotype=="HET","KO",Genotype)) %>% 
  unique -> animal_data_m
###
LF_data2 %>% filter(Parameter==par_val) %>% pivot_longer(cols=4:99,names_to = "Animal.ID") %>%
  group_by(`Animal.ID`,ZT) %>%
  filter(Mch_Conc!=0) %>%
  merge(animal_data_m ,by="Animal.ID")-> LF_data2

LF_data2 %>% rename(Sample=Animal.ID,Mch_conc=Mch_Conc,Value=value) %>%
  dplyr::select(-Parameter) -> LF_data2

LF_data2 %>% filter(!is.na(Value)) -> LF_data

# mann whitney u on pairwise zts ----------------------------------

mw_results<-dr_MWU_pairwise(LF_data)

# dose response model -----------------------------------------------------

param_formodel <- dr_fit(LF_data) #uncomment these lines to rerun dr curve fits - warning takes some time
save(param_formodel,file = "data/drcmodelparams_flex_LFmaster_RN.RData")
load("data/drcmodelparams_flex_LFmaster_RN.RData")

anova_pvals_upper <- dr_anova(param_formodel,"Upper","log10(params) ~ ZT * Treatment") # not used

anova_pvals_slope <- dr_anova(param_formodel,"Slope","log10(-params) ~ ZT * Treatment")

# plotting dose response curve --------------------------------------------

p1<-dr_plot(LF_data,anova_pvals_slope,mw_results,c(0,2.5),y_lab="Maximum Airway Resistance R<sub>rs</sub>(cm.H<sub>2</sub>O.s.ml<sup>-1</sup>)")
p1

png("plots/flex_meth_dose_response_LFmaster_RN.png", width = 2400, height = 1500, res = 300)  # adjust size/res as needed
grid.draw(
  p1
)
grid.text( expression("Methacholine Concentration (mg.mL"^"-1"*")"), y = unit(0.03, "npc"), gp = gpar(fontsize = 10))
dev.off()


# AUC sinusoidal analysis -----------------------------------------------------

rhy_plot(LF_data,"AUC",y_lim=c(-25,25),y_lab="AUC Airway Resistance R<sub>rs</sub>(cm.H<sub>2</sub>O.s.ml<sup>-1</sup>)") -> analysis_out

p_auc <- analysis_out$combined
p_auc
ggsave(p_auc,filename="plots/flex_AUC_LFmaster_RN.png",width=8,height=5)

# Max sinusoidal analysis -----------------------------------------------------

rhy_plot(LF_data,"Max",y_lim=c(-.75,.75),y_lab="Max Airway Resistance R<sub>rs</sub>(cm.H<sub>2</sub>O.s.ml<sup>-1</sup>)") -> analysis_out

p_max <- analysis_out$combined
p_max
ggsave(p_max,filename="plots/flex_max_LFmaster_RN.png",width=8,height=5)
