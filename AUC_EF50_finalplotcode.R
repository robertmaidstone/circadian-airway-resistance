library(tidyverse)
library(openxlsx)
library(lmerTest)
library(factoextra)

read.xlsx("~/Amlan/DCP/LungFunctionMasterDataSet.xlsx",sheet=2) -> animal_data
read.xlsx("~/Amlan/DCP/LungFunctionMasterDataSet.xlsx",sheet=1,fillMergedCells = TRUE) -> LF_data

###

LF_data %>% filter(Parameter=="DCP_EF50") %>% pivot_longer(cols=4:39,names_to = "Animal.ID") %>%
  group_by(`Animal.ID`,ZT) %>%
  filter(Mch_Conc!=0) %>%
  mutate(value=log10(value)) %>% #is this still right for DCP??
  mutate(Max_Value=max(value,na.rm=T)) %>%
  mutate(AUC=sum(diff(Mch_Conc)*(value[-1]+value[-4])/2,na.rm=T)) %>%
  dplyr::select(`Animal.ID`,ZT,Max_Value,AUC) %>%
  unique %>% ungroup %>%
  #mutate(AUC=Max_Value) %>%
  merge(animal_data,by="Animal.ID")-> sum_data_sGAW

lmer(AUC~1+Genotype*Treatment*sin(2*pi/24*ZT) + Genotype*Treatment*cos(2*pi/24*ZT)+Genotype*Treatment+(1|Animal.ID),data=sum_data_sGAW) -> lm_sGAW
summary(lm_sGAW)
anova(lm_sGAW)

lm(AUC~1+Genotype*Treatment*sin(2*pi/24*ZT) + Genotype*Treatment*cos(2*pi/24*ZT)+Genotype*Treatment,data=sum_data_sGAW) -> lm_sGAW
summary(lm_sGAW)

newdata <- expand.grid(
  ZT = seq(from=0,to=24,length.out=12),
  Treatment = unique(sum_data_sGAW$Treatment),
  Genotype = unique(sum_data_sGAW$Genotype)
  #Animal.ID = unique(sum_data_sGAW$Animal.ID)
)

pred_resp <- predict.lm(lm_sGAW, newdata = newdata, se.fit=TRUE, interval="confidence", level=0.95)
newdata$Predicted_Response <- pred_resp$fit[,1]
newdata$Predicted_Response_LCI <- pred_resp$fit[,2]
newdata$Predicted_Response_UCI <- pred_resp$fit[,3]

newdata -> predict_values_sGAW

### nls

nls_model <- nls(AUC ~ I + A * sin(2*pi/24*ZT + phi) + (Genotype=="HET")*(I1+ A1 * sin(2*pi/24*ZT + phi1)) + (Treatment=="HDM")*(I2+ A2 * sin(2*pi/24*ZT + phi2)), data = sum_data_sGAW,
                 start = list(A = 1, phi = 0,A1 = 0.1, phi1 = 0.1,A2 = 0.1, phi2 = 0.1, I=0,I1=0,I2=0))
summary(nls_model)

### plotting
sum_data_sGAW$AUC %>% summary

col_vec <- c("#0072B2","#E69F00")
y_limit <- c(-11,12)

sig_text <- "Mesor Treatment **\nPhase Genotype ***\nPhase Treatment ***"
gen_labels <- c("CCSP-Reverbα WT","CCSP-Reverbα KO")

sum_data_sGAW %>%
  mutate(Genotype=factor(Genotype,levels=c("WT","HET"),labels=c("WT","KO"))) %>%
  filter(Treatment=="PBS") %>%
  ggplot(aes(x = ZT, y = AUC, colour = Genotype)) +
  geom_point() +
  geom_ribbon(aes(ymax=Predicted_Response_UCI,ymin=Predicted_Response_LCI,y=NULL,fill=Genotype,colour=NULL),
              data=predict_values_sGAW %>%filter(Treatment=="PBS") %>% mutate(Genotype=factor(Genotype,levels=c("WT","HET"),
                                                                  labels=c("WT","KO"))),size=1,alpha=0.2) +
  geom_line(aes(y=Predicted_Response),
            data=predict_values_sGAW %>%filter(Treatment=="PBS") %>% mutate(Genotype=factor(Genotype,levels=c("WT","HET"),
                                                                labels=c("WT","KO"))),size=1,linetype="dashed") + 
  ylab("AUC EF50 (ml.sec<sup>-1</sup>)") + scale_colour_manual(values=col_vec) + scale_fill_manual(values=col_vec) +
  scale_x_continuous(breaks=c(0,6,12,18,24))+
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_markdown()) +
  ggtitle("PBS")+
  annotate("text", x = 1, y = -7, label = sig_text, hjust = 0, vjust = 1, size = 4)+
  ylim(y_limit)-> p_pbs

sum_data_sGAW %>%
  mutate(Genotype=factor(Genotype,levels=c("WT","HET"),labels=gen_labels)) %>%
  filter(Treatment=="HDM") %>%
  ggplot(aes(x = ZT, y = AUC, colour = Genotype)) +
  geom_point() +
  geom_ribbon(aes(ymax=Predicted_Response_UCI,ymin=Predicted_Response_LCI,y=NULL,fill=Genotype,colour=NULL),
              data=predict_values_sGAW %>%filter(Treatment=="HDM") %>% mutate(Genotype=factor(Genotype,levels=c("WT","HET"),
                                                                                              labels=gen_labels)),size=1,alpha=0.2) +
  geom_line(aes(y=Predicted_Response),
            data=predict_values_sGAW %>%filter(Treatment=="HDM") %>% mutate(Genotype=factor(Genotype,levels=c("WT","HET"),
                                                                                            labels=gen_labels)),size=1,linetype="dashed") + 
  ylab("AUC") + scale_colour_manual(values=col_vec) + scale_fill_manual(values=col_vec) +
  scale_x_continuous(breaks=c(0,6,12,18,24))+
  theme_bw() +
  theme(legend.position=c(0.75, 0.9),
        legend.text = element_markdown(),
        legend.title = element_blank(),
        legend.background = element_blank(),
        axis.title.y = element_blank())+
  guides(fill = guide_legend(override.aes = list(colour = NA,alpha=1)))+
  ggtitle("HDM")+
  ylim(y_limit)-> p_hdm


library(patchwork)
p_pbs + p_hdm

ggsave(p_pbs + p_hdm,filename="~/Amlan/DCP/p_AUC_EF50_final.png",width=8,height=5)

#################################################################

######KO. HDM
test_data <- sum_data_sGAW %>% filter(Genotype=="HET",Treatment=="HDM")
nls_model <- nls(AUC ~ I + A * sin(2*pi/24*ZT + phi), data = test_data,
                 start = list(A = 1, phi = 0, I=0))
#summary(nls_model)
lm_model <- lm(AUC ~ 1 , data = test_data)
#summary(lm_model)
logLik(lm_model)
lmtest::lrtest(nls_model,lm_model)

newdata <- expand.grid(
  ZT = seq(from=0,to=24,length.out=12),
  Treatment = "HDM",
  Genotype = "HET"
  #Animal.ID = unique(sum_data_sGAW$Animal.ID)
)

pred_resp <- predict.lm(lm_model, newdata = newdata, se.fit=TRUE, interval="confidence", level=0.95)
newdata$Predicted_Response <- pred_resp$fit[,1]
newdata$Predicted_Response_LCI <- pred_resp$fit[,2]
newdata$Predicted_Response_UCI <- pred_resp$fit[,3]
newdata_mast <- newdata

######KO. PBS
test_data <- sum_data_sGAW %>% filter(Genotype=="HET",Treatment=="PBS")
nls_model <- nls(AUC ~ I + A * sin(2*pi/24*ZT + phi), data = test_data,
                 start = list(A = 1, phi = 0, I=0))
#summary(nls_model)
lm_model <- lm(AUC ~ 1 , data = test_data)
#summary(lm_model)
logLik(lm_model)
lmtest::lrtest(nls_model,lm_model)

newdata <- expand.grid(
  ZT = seq(from=0,to=24,length.out=12),
  Treatment = "PBS",
  Genotype = "HET"
  #Animal.ID = unique(sum_data_sGAW$Animal.ID)
)

pred_resp <- predict.lm(lm_model, newdata = newdata, se.fit=TRUE, interval="confidence", level=0.95)
newdata$Predicted_Response <- pred_resp$fit[,1]
newdata$Predicted_Response_LCI <- pred_resp$fit[,2]
newdata$Predicted_Response_UCI <- pred_resp$fit[,3]
newdata_mast <- rbind(newdata,newdata_mast)

######WT PBS
test_data <- sum_data_sGAW %>% filter(Genotype=="WT",Treatment=="PBS")
nls_model <- nls(AUC ~ I + A * sin(2*pi/24*ZT + phi), data = test_data,
                 start = list(A = 1, phi = 0, I=0))
#summary(nls_model)
lm_model <- lm(AUC ~ 1 , data = test_data)
#summary(lm_model)
logLik(lm_model)
lmtest::lrtest(nls_model,lm_model)

newdata <- expand.grid(
  ZT = seq(from=0,to=24,length.out=12),
  Treatment = "PBS",
  Genotype = "WT"
  #Animal.ID = unique(sum_data_sGAW$Animal.ID)
)

pred_resp <- predict.lm(lm_model, newdata = newdata, se.fit=TRUE, interval="confidence", level=0.95)
newdata$Predicted_Response <- pred_resp$fit[,1]
newdata$Predicted_Response_LCI <- pred_resp$fit[,2]
newdata$Predicted_Response_UCI <- pred_resp$fit[,3]
newdata_mast <- rbind(newdata,newdata_mast)

######WT HDM
test_data <- sum_data_sGAW %>% filter(Genotype=="WT",Treatment=="HDM")
nls_model <- nls(AUC ~ I + A * sin(2*pi/24*ZT + phi), data = test_data,
                 start = list(A = 1, phi = 0, I=0))
#summary(nls_model)
lm_model <- lm(AUC ~ 1 , data = test_data)
#summary(lm_model)
logLik(lm_model)
lmtest::lrtest(nls_model,lm_model)

newdata <- expand.grid(
  ZT = seq(from=0,to=24,length.out=12),
  Treatment = "HDM",
  Genotype = "WT"
  #Animal.ID = unique(sum_data_sGAW$Animal.ID)
)

pred_resp <- predict.lm(lm_model, newdata = newdata, se.fit=TRUE, interval="confidence", level=0.95)
newdata$Predicted_Response <- pred_resp$fit[,1]
newdata$Predicted_Response_LCI <- pred_resp$fit[,2]
newdata$Predicted_Response_UCI <- pred_resp$fit[,3]
newdata_mast <- rbind(newdata,newdata_mast)


sum_data_sGAW %>%
  mutate(Genotype=factor(Genotype,levels=c("WT","HET"),labels=c("WT","KO"))) %>%
  filter(Treatment=="PBS") %>%
  ggplot(aes(x = ZT, y = AUC, colour = Genotype)) +
  geom_point() +
  geom_ribbon(aes(ymax=Predicted_Response_UCI,ymin=Predicted_Response_LCI,y=NULL,fill=Genotype,colour=NULL),
              data=newdata_mast %>%filter(Treatment=="PBS") %>% mutate(Genotype=factor(Genotype,levels=c("WT","HET"),
                                                                                       labels=c("WT","KO"))),size=1,alpha=0.2) +
  geom_line(aes(y=Predicted_Response),
            data=newdata_mast %>%filter(Treatment=="PBS") %>% mutate(Genotype=factor(Genotype,levels=c("WT","HET"),
                                                                                     labels=c("WT","KO"))),size=1) + 
  ylab("AUC EF50 (ml.sec<sup>-1</sup>)")  + scale_colour_manual(values=col_vec) + scale_fill_manual(values=col_vec) +
  scale_x_continuous(breaks=c(0,6,12,18,24))+
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_markdown()) +
  ggtitle("PBS")+
  annotate("text", x = 1, y = -6, label = sig_text, hjust = 0, vjust = 1, size = 4)+
  ylim(y_limit)-> p_pbs

sum_data_sGAW %>%
  mutate(Genotype=factor(Genotype,levels=c("WT","HET"),labels=gen_labels)) %>%
  filter(Treatment=="HDM") %>%
  ggplot(aes(x = ZT, y = AUC, colour = Genotype)) +
  geom_point() +
  geom_ribbon(aes(ymax=Predicted_Response_UCI,ymin=Predicted_Response_LCI,y=NULL,fill=Genotype,colour=NULL),
              data=newdata_mast %>%filter(Treatment=="HDM") %>% mutate(Genotype=factor(Genotype,levels=c("WT","HET"),
                                                                                       labels=gen_labels)),size=1,alpha=0.2) +
  geom_line(aes(y=Predicted_Response),
            data=newdata_mast %>%filter(Treatment=="HDM") %>% mutate(Genotype=factor(Genotype,levels=c("WT","HET"),
                                                                                     labels=gen_labels)),size=1) + 
  ylab("Min Value") + scale_colour_manual(values=col_vec) + scale_fill_manual(values=col_vec) +
  scale_x_continuous(breaks=c(0,6,12,18,24))+
  theme_bw() +
  theme(legend.position=c(0.9, 0.9),
        legend.title = element_blank(),
        legend.background = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_markdown())+
  guides(fill = guide_legend(override.aes = list(colour = NA,alpha=1)))+
  ggtitle("HDM")+
  ylim(y_limit)-> p_hdm


library(patchwork)
p_pbs + p_hdm
ggsave(p_pbs + p_hdm,filename="~/Amlan/plots_200126/p_AUC_EF50_final.png",width=8,height=5)

