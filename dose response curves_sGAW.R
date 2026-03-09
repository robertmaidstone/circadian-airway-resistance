library(tidyverse)
library(openxlsx)
library(grid)
library(patchwork)
library(ggtext)


read.xlsx("~/Amlan/DCP/LungFunctionMasterDataSet.xlsx",sheet=2) -> animal_data
read.xlsx("~/Amlan/DCP/LungFunctionMasterDataSet.xlsx",sheet=1,fillMergedCells = TRUE) -> LF_data
###
LF_data %>% filter(Parameter=="DCP_sGAW") %>% pivot_longer(cols=4:39,names_to = "Animal.ID") %>%
  group_by(`Animal.ID`,ZT) %>%
  filter(Mch_Conc!=0) %>%
  mutate(value=log10(value)) %>%
  merge(animal_data,by="Animal.ID")-> LF_data


LF_data %>% 
  group_by(ZT,Genotype,Treatment,Mch_Conc) %>%
  mutate(Med_Value=mean(value,na.rm=T)) %>%
  dplyr::select(-Animal.ID,-value) %>%
  unique %>% ungroup -> t_data

t_data %>%
  #filter(Genotype=="WT") %>%
  mutate(ZT=as.character(ZT)) %>%
  ggplot(aes(x=Mch_Conc,y=Med_Value)) + geom_line(aes(linetype = Treatment, color = ZT))+
  geom_point(aes(color = ZT))+
  scale_x_continuous(breaks=c(0,3.12,6.25,12.5,25,50))+
  theme_bw() +
  facet_grid(~Genotype)


library(drc)

# Fit dose response model for each geotype, treatment, time triplicate
LF_data$Group <- interaction(LF_data$Genotype, LF_data$Treatment,LF_data$ZT)


model <- drm(value ~ Mch_Conc, data = LF_data,
             curveid = Group,
             fct = LL.4(names = c("Slope", "Lower", "Upper", "EC50")))

summary(model)


# Extract coefficients
params <- coef(model)

# Convert to data frame
param_df <- as.data.frame(params)

# Add group names as a column
param_df$Group <- rownames(param_df)

# Separate Group into Time and Condition
param_df <- param_df %>%
  separate(Group, into = c("Genotype", "Treatment","Time"), sep = "\\.") %>%
  separate(Genotype, into = c("Parameter","Genotype"), sep = ":")


param_df$Time <- factor(param_df$Time)
param_df$Treatment <- factor(param_df$Treatment)
param_df$Genotype <- factor(param_df$Genotype)

# ANOVA for Emax

aov_emax <- aov(params ~ Time * Treatment, data = param_df %>% filter(Genotype=="WT",Parameter=="Upper"))
summary(aov_emax)

# ANOVA for Slope
aov_slope <-  aov(params ~ Time * Treatment, data = param_df %>% filter(Genotype=="WT",Parameter=="Slope"))
summary(aov_slope)

#######

# Fit dose response model
LF_data$Group <- interaction(LF_data$Animal.ID,LF_data$ZT)


model <- drm(value ~ Mch_Conc, data = LF_data,
             curveid = Group,
             fct = LL.4(names = c("Slope", "Lower", "Upper", "EC50")))

summary(model)


# Extract coefficients
params <- coef(model)

# Convert to data frame
param_df <- as.data.frame(params)

# Add group names as a column
param_df$Group <- rownames(param_df)

# Separate Group into Time and Condition
param_df <- param_df %>%
  separate(Group, into = c("Parameter","Sample"), sep = ":")

param_df %>% merge(LF_data %>% dplyr::select(Animal.ID,Genotype,Treatment,ZT) %>% unique %>% mutate(Sample=paste(Animal.ID,ZT,sep = ".")),by="Sample")-> param_formodel

param_formodel$ZT <- factor(param_formodel$ZT)
param_formodel$Treatment <- factor(param_formodel$Treatment)
param_formodel$Genotype <- factor(param_formodel$Genotype)

#save(param_formodel,file = "drcmodelparams_EF50.RData")

# ANOVA for Emax
param_formodel %>% filter(Genotype=="WT",Parameter=="Upper") %>% dplyr::select(params) %>% unlist %>% log10 %>% hist
aov_emax <- aov(log10(params+1) ~ ZT * Treatment, data = param_formodel %>% filter(Genotype=="WT",Parameter=="Upper"))
summary(aov_emax)
plot(aov_emax)

log10((param_formodel %>% filter(Genotype=="HET",Parameter=="Upper") %>% dplyr::select(params) %>% unlist )+1) %>% hist
aov_emax <- aov(log10(params+1) ~ ZT * Treatment, data = param_formodel %>% filter(Genotype=="HET",Parameter=="Upper"))
summary(aov_emax)
plot(aov_emax)

ggplot(param_formodel %>% filter(Parameter == "Upper", Genotype == "HET"),
       aes(x = ZT, y = log10(params), color = Treatment)) +
  geom_point(position = position_jitter(width = 0.2)) +
  stat_summary(fun = mean, geom = "line", aes(group = Treatment)) +
  labs(title = "Emax by ZT and Treatment (KO)")

# ANOVA for Slope
param_formodel %>% filter(Genotype=="WT",Parameter=="Slope") %>% dplyr::select(params) %>% unlist %>% (function(x){log10((x))}) %>%hist
aov_slope <-  aov((params) ~ ZT * Treatment, data = param_formodel %>% filter(Genotype=="WT",Parameter=="Slope"))
summary(aov_slope)
plot(aov_slope)

param_formodel %>% filter(Genotype=="HET",Parameter=="Slope") %>% dplyr::select(params) %>% unlist %>% (function(x){log10((x*-1))}) %>%hist
aov_slope <-  aov(params ~ ZT * Treatment, data = param_formodel %>% filter(Genotype=="HET",Parameter=="Slope"))
summary(aov_slope)
plot(aov_slope)

ggplot(param_formodel %>% filter(Parameter == "Slope", Genotype == "KO"),
       aes(x = ZT, y = log10(-params), color = Treatment)) +
  geom_point(position = position_jitter(width = 0.2)) +
  stat_summary(fun = mean, geom = "line", aes(group = Treatment)) +
  labs(title = "Slope by ZT and Treatment (KO)")


gen_labels <- c("CCSP-Reverbα WT","CCSP-Reverbα KO")
sig_text <- "Max Dose Treatment *"

t_data %>%
  filter(Genotype=="WT") %>%
  mutate(ZT=as.character(ZT)) %>%
  ggplot(aes(x=Mch_Conc,y=Med_Value)) + geom_line(aes(linetype = Treatment, color = ZT))+
  geom_point(aes(color = ZT))+
  scale_x_continuous(breaks=c(0,3.12,6.25,12.5,25,50))+
  theme_bw()+
  ylab("Median sGAW (cm.H<sub>2</sub>O.sec<sup>-1</sup>)") +
  ylim(c(-0.75,-0.25))+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(),
        plot.title = element_markdown(),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_manual(values = c("#0072B2","#E69F00","#D55E00", "#009E73")) +
  annotate("text", x = 1, y = -0.7, label = sig_text, hjust = 0, vjust = 1, size = 4)+
  ggtitle(gen_labels[1])-> p1

sig_text <- ""

t_data %>%
  filter(Genotype=="HET") %>%
  mutate(ZT=as.character(ZT)) %>%
  ggplot(aes(x=Mch_Conc,y=Med_Value)) + geom_line(aes(linetype = Treatment, color = ZT))+
  geom_point(aes(color = ZT))+
  scale_x_continuous(breaks=c(0,3.12,6.25,12.5,25,50))+
  theme_bw() +
  xlab("")+ 
  ylim(c(-0.75,-0.25))+
  theme(axis.title.y = element_blank(),
        plot.title = element_markdown(),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_manual(values = c("#0072B2","#E69F00","#D55E00", "#009E73")) +
  annotate("text", x = 1, y = .85, label = sig_text, hjust = 0, vjust = 1, size = 4)+
  ggtitle(gen_labels[2])->p2

png("~/Amlan/plots_200126/meth_dose_response_sGAW.png", width = 2400, height = 1500, res = 300)  # adjust size/res as needed
#grid.newpage()
grid.draw(
  patchworkGrob(p1 | p2)
)
grid.text( expression("Methacholine Concentration (mg.mL"^"-1"*")"), y = unit(0.03, "npc"), gp = gpar(fontsize = 10))
dev.off()


### anova on mch_conc 50


lm(value~Treatment*ZT,data=LF_data%>% filter(Mch_Conc==25,Genotype=="WT")) %>% anova
lm(value~Treatment*ZT,data=LF_data%>% filter(Mch_Conc==25,Genotype=="HET")) %>% anova

lm(value~ZT,data=LF_data%>% filter(Mch_Conc==25,Genotype=="WT",Treatment=="HDM")) %>% anova
lm(value~ZT,data=LF_data%>% filter(Mch_Conc==25,Genotype=="HET",Treatment=="HDM")) %>% anova
