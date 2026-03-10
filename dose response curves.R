library(tidyverse)
library(openxlsx)
library(grid)
library(patchwork)
library(ggtext)
library(drc)


# load data ---------------------------------------------------------------
setwd("~/Amlan/AirwayResistance/")
source("functions.R")
read.xlsx("data/AC_FVData_Complete.xlsx",sheet=5,fillMergedCells = TRUE) -> LF_data

LF_data %>% filter(Mch_conc!=0) -> LF_data

# anovas replace with tests as discussed ----------------------------------

library(lme4)
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

# ANOVA for Emax
anova_pvals_upper <- dr_anova(param_formodel,"Upper","log10(params) ~ ZT * Treatment")

# ANOVA for Slope
anova_pvals_slope <- dr_anova(param_formodel,"Slope","log10(-params) ~ ZT * Treatment")

# plotting dose response curve --------------------------------------------

dr_plot(LF_data,anova_pvals_slope,anova_pvals_upper)

gen_labels <- c("CCSP-Reverbα WT","CCSP-Reverbα KO")
sig_text <- "Max Dose Treatment **"
#sig_text <- ""

t_data %>%
  filter(Genotype=="WT") %>%
  mutate(ZT=as.character(ZT)) %>%
  ggplot(aes(x=Mch_conc,y=Med_Value)) + geom_line(aes(linetype = Treatment, color = ZT))+
  geom_point(aes(color = ZT))+
  scale_x_continuous(breaks=c(0,3.12,6.25,12.5,25,50))+
  theme_bw()+
  ylab("Maximum Airway Resistance R<sub>rs</sub>(cm.H<sub>2</sub>O.s.ml<sup>-1</sup>)") +
  ylim(c(0.5,6))+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(),
        plot.title = element_markdown(),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_manual(values = c("#0072B2","#E69F00","#D55E00", "#009E73")) +
  annotate("text", x = 1, y = 6, label = sig_text, hjust = 0, vjust = 1, size = 4)+
  ggtitle(gen_labels[1])-> p1

sig_text <- "Max Dose Treatment ***"

t_data %>%
  filter(Genotype=="KO") %>%
  mutate(ZT=as.character(ZT)) %>%
  ggplot(aes(x=Mch_conc,y=Med_Value)) + geom_line(aes(linetype = Treatment, color = ZT))+
  geom_point(aes(color = ZT))+
  scale_x_continuous(breaks=c(0,3.12,6.25,12.5,25,50))+
  theme_bw() +
  xlab("")+ 
  ylim(c(0.5,6))+
  theme(axis.title.y = element_blank(),
        plot.title = element_markdown(),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_manual(values = c("#0072B2","#E69F00","#D55E00", "#009E73")) +
  annotate("text", x = 1, y = 6, label = sig_text, hjust = 0, vjust = 1, size = 4)+
  ggtitle(gen_labels[2])->p2

png("~/Amlan/Flexivent_Jan25/meth_dose_response_v3.png", width = 2400, height = 1500, res = 300)  # adjust size/res as needed
#grid.newpage()
grid.draw(
  patchworkGrob(p1 | p2)
)
grid.text( expression("Methacholine Concentration (mg.mL"^"-1"*")"), y = unit(0.03, "npc"), gp = gpar(fontsize = 10))
dev.off()
