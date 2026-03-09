library(tidyverse)
library(openxlsx)
library(lmerTest)
library(factoextra)
library(ggtext)
gen_labels <- c("CCSP-Reverbα WT","CCSP-Reverbα KO")
read.xlsx("~/Amlan/Flexivent_Jan25/AC_FVData_Complete.xlsx",sheet=5,fillMergedCells = TRUE) -> LF_data

###
LF_data %>% filter(Mch_conc!=0) %>%
  mutate(Value=log10(Value)) -> LF_data

###

LF_data %>% 
  group_by(Sample,ZT,Genotype,Treatment) %>%
  mutate(Max_Value=max(Value,na.rm=T)) %>%
  mutate(AUC=sum(diff(Mch_conc)*(Value[-1]+Value[-5])/2)) %>%
  dplyr::select(Sample,Max_Value,AUC) %>%
  mutate(AUC=Max_Value) %>%
  unique %>% ungroup -> sum_data

lm(AUC~1+Genotype*Treatment*sin(2*pi/24*ZT) + Genotype*Treatment*cos(2*pi/24*ZT)+Genotype*Treatment,data=sum_data) -> lm_sGAW
summary(lm_sGAW)

newdata <- expand.grid(
  ZT = seq(from=0,to=24,length.out=12),
  Treatment = unique(sum_data$Treatment),
  Genotype = unique(sum_data$Genotype)
  #Animal.ID = unique(sum_data_sGAW$Animal.ID)
)

pred_resp <- predict.lm(lm_sGAW, newdata = newdata, se.fit=TRUE, interval="confidence", level=0.95)
newdata$Predicted_Response <- pred_resp$fit[,1]
newdata$Predicted_Response_LCI <- pred_resp$fit[,2]
newdata$Predicted_Response_UCI <- pred_resp$fit[,3]

newdata %>%   mutate(Genotype=factor(Genotype,levels=c("WT","KO"),labels=gen_labels)) -> predict_values_sGAW

### nls

nls_model <- nls(AUC ~ I + A * sin(2*pi/24*ZT + phi) + (Genotype=="KO")*(I1+ A1 * sin(2*pi/24*ZT + phi1)) + (Treatment=="HDM")*(I2+ A2 * sin(2*pi/24*ZT + phi2)), data = sum_data,
                 start = list(A = 1, phi = 0,A1 = 0.1, phi1 = 0.1,A2 = 0.1, phi2 = 0.1, I=0,I1=0,I2=0))
summary(nls_model)
summary(nls_model)[[10]][,4] %>% round(4)
nls_model_o <- nls(AUC ~ I , data = sum_data,
                 start = list(I=0))
summary(nls_model_o)
lmtest::lrtest(nls_model,nls_model_o)

##PBS only
nls_model <- nls(AUC ~ I + A * sin(2*pi/24*ZT + phi) + (Genotype=="KO")*(I1+ A1 * sin(2*pi/24*ZT + phi1)), data = sum_data %>% filter(Treatment=="PBS"),
                 start = list(A = 1, phi = 0,A1 = 0.1, phi1 = 0.1, I=0,I1=0))
summary(nls_model)
summary(nls_model)[[10]][,4] %>% round(4)
##HDM only
nls_model <- nls(AUC ~ I + A * sin(2*pi/24*ZT + phi) + (Genotype=="KO")*(I1+ A1 * sin(2*pi/24*ZT + phi1)), data = sum_data %>% filter(Treatment=="HDM"),
                 start = list(A = 1, phi = 0,A1 = 0.1, phi1 = 0.1, I=0,I1=0))
summary(nls_model)
summary(nls_model)[[10]][,4] %>% round(4)
### plotting


col_vec <- c("#0072B2","#E69F00")
y_limit <- c(-.1,1.3)

sig_text <- ""
#gen_labels <- c("Cre<sup>-/-</sup>","Cre<sup>+/-</sup>")
gen_labels <- c("CCSP-Reverbα WT","CCSP-Reverbα KO")

sum_data %>%
  mutate(Genotype=factor(Genotype,levels=c("WT","KO"),labels=gen_labels)) %>%
  filter(Treatment=="PBS") %>%
  ggplot(aes(x = ZT, y = AUC, colour = Genotype)) +
  geom_point() +
  geom_ribbon(aes(ymax=Predicted_Response_UCI,ymin=Predicted_Response_LCI,y=NULL,fill=Genotype,colour=NULL),
              data=predict_values_sGAW %>%filter(Treatment=="PBS"),size=1,alpha=0.2) +
  geom_line(aes(y=Predicted_Response),
            data=predict_values_sGAW %>%filter(Treatment=="PBS"),size=1,linetype="dashed") +
  ylab("Max Airway Resistance R<sub>rs</sub>(cm.H<sub>2</sub>O.s.ml<sup>-1</sup>)") + scale_colour_manual(values=col_vec) + scale_fill_manual(values=col_vec) +
  scale_x_continuous(breaks=c(0,6,12,18,24),limits = c(0,24))+
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_markdown()) +
  ggtitle("PBS")+
  annotate("text", x = 1, y = max(y_limit), label = sig_text, hjust = 0, vjust = 1, size = 4)+
  ylim(y_limit)-> p_pbs

sig_text <- "Amplitude **\nPhase****"

sum_data %>%
  mutate(Genotype=factor(Genotype,levels=c("WT","KO"),labels=gen_labels)) %>%
  filter(Treatment=="HDM") %>%
  ggplot(aes(x = ZT, y = AUC, colour = Genotype)) +
  geom_point() +
  geom_ribbon(aes(ymax=Predicted_Response_UCI,ymin=Predicted_Response_LCI,y=NULL,fill=Genotype,colour=NULL),
              data=predict_values_sGAW %>%filter(Treatment=="HDM"),size=1,alpha=0.2) +
  geom_line(aes(y=Predicted_Response),
            data=predict_values_sGAW %>%filter(Treatment=="HDM",Genotype=="CCSP-Reverbα WT"),size=1) + 
  geom_line(aes(y=Predicted_Response),
            data=predict_values_sGAW %>%filter(Treatment=="HDM",Genotype=="CCSP-Reverbα KO"),size=1,linetype="dashed") +
  ylab("Max Value") + scale_colour_manual(values=col_vec,drop=FALSE) + scale_fill_manual(values=col_vec,drop=FALSE) +
  scale_x_continuous(breaks=c(0,6,12,18,24),limits = c(0,24))+
  theme_bw() +
  theme(legend.position=c(0.75, 0.9),
        legend.title = element_blank(),
        legend.text = element_markdown(),
        legend.background = element_blank(),
        axis.title.y = element_blank())+
  guides(fill = guide_legend(override.aes = list(colour = col_vec,alpha=1),drop=FALSE))+
  ggtitle("HDM")+
  annotate("text", x = 1, y = max(y_limit), label = sig_text, hjust = 0, vjust = 1, size = 4)+
  ylim(y_limit)-> p_hdm


library(patchwork)
p_pbs + p_hdm

ggsave(p_pbs + p_hdm,filename="~/Amlan/Flexivent_Jan25/p_Max_final_v4b.png",width=8,height=5)


######KO. HDM
test_data <- sum_data %>% filter(Genotype=="KO",Treatment=="HDM")
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
  Genotype = "KO"
  #Animal.ID = unique(sum_data_sGAW$Animal.ID)
)

pred_resp <- predict.lm(lm_model, newdata = newdata, se.fit=TRUE, interval="confidence", level=0.95)
newdata$Predicted_Response <- pred_resp$fit[,1]
newdata$Predicted_Response_LCI <- pred_resp$fit[,2]
newdata$Predicted_Response_UCI <- pred_resp$fit[,3]
newdata_mast <- newdata

######KO. PBS
test_data <- sum_data %>% filter(Genotype=="KO",Treatment=="PBS")
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
  Genotype = "KO"
  #Animal.ID = unique(sum_data_sGAW$Animal.ID)
)

pred_resp <- predict.lm(lm_model, newdata = newdata, se.fit=TRUE, interval="confidence", level=0.95)
newdata$Predicted_Response <- pred_resp$fit[,1]
newdata$Predicted_Response_LCI <- pred_resp$fit[,2]
newdata$Predicted_Response_UCI <- pred_resp$fit[,3]
newdata_mast <- rbind(newdata,newdata_mast)


######WT. HDM
test_data <- sum_data %>% filter(Genotype=="WT",Treatment=="HDM")
nls_model <- nls(AUC ~ I + A * sin(2*pi/24*ZT + phi), data = test_data,
                 start = list(A = 1, phi = 0, I=0))
#summary(nls_model)
lm_model <- lm(AUC ~ 1 , data = test_data)
#summary(lm_model)
logLik(lm_model)
lmtest::lrtest(nls_model,lm_model)

lm(AUC~1+sin(2*pi/24*ZT) + cos(2*pi/24*ZT),data=test_data) -> lm_sGAW
summary(lm_sGAW)
lmtest::lrtest(lm_sGAW,lm_model)

newdata <- expand.grid(
  ZT = seq(from=0,to=24,length.out=12),
  Treatment = "HDM",
  Genotype = "WT"
  #Animal.ID = unique(sum_data_sGAW$Animal.ID)
)

pred_resp <- predict.lm(lm_sGAW, newdata = newdata, se.fit=TRUE, interval="confidence", level=0.95)
newdata$Predicted_Response <- pred_resp$fit[,1]
newdata$Predicted_Response_LCI <- pred_resp$fit[,2]
newdata$Predicted_Response_UCI <- pred_resp$fit[,3]
newdata_mast <- rbind(newdata,newdata_mast)

#WT. PBS
test_data <- sum_data %>% filter(Genotype=="WT",Treatment=="PBS")
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

newdata_mast <- newdata_mast %>%  mutate(Genotype=factor(Genotype,levels=c("WT","KO"),labels=gen_labels)) 

sum_data %>%
  mutate(Genotype=factor(Genotype,levels=c("WT","KO"),labels=gen_labels)) %>%
  filter(Treatment=="PBS") %>%
  ggplot(aes(x = ZT, y = AUC, colour = Genotype)) +
  geom_point() +
  geom_ribbon(aes(ymax=Predicted_Response_UCI,ymin=Predicted_Response_LCI,y=NULL,fill=Genotype,colour=NULL),
              data=newdata_mast %>%filter(Treatment=="PBS") ,size=1,alpha=0.2) +
  geom_line(aes(y=Predicted_Response),
            data=newdata_mast %>%filter(Treatment=="PBS") ,size=1) +
  ylab("Max Airway Resistance R<sub>rs</sub>(cm.H<sub>2</sub>O.s.ml<sup>-1</sup>)") + scale_colour_manual(values=col_vec) + scale_fill_manual(values=col_vec) +
  scale_x_continuous(breaks=c(0,6,12,18,24),limits = c(0,24))+
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y = element_markdown()) +
  ggtitle("PBS")+
  annotate("text", x = 1, y = max(y_limit), label = sig_text, hjust = 0, vjust = 1, size = 4)+
  ylim(y_limit)-> p_pbs

sum_data %>%
  mutate(Genotype=factor(Genotype,levels=c("WT","KO"),labels=gen_labels)) %>%
  filter(Treatment=="HDM") %>%
  ggplot(aes(x = ZT, y = AUC, colour = Genotype)) +
  geom_point() +
  geom_ribbon(aes(ymax=Predicted_Response_UCI,ymin=Predicted_Response_LCI,y=NULL,fill=Genotype,colour=NULL),
              data=newdata_mast %>%filter(Treatment=="HDM") ,size=1,alpha=0.2) +
  geom_line(aes(y=Predicted_Response),
            data=newdata_mast %>%filter(Treatment=="HDM") ,size=1)+ 
  ylab("Max Value") + scale_colour_manual(values=col_vec,drop=FALSE) + scale_fill_manual(values=col_vec,drop=FALSE) +
  scale_x_continuous(breaks=c(0,6,12,18,24),limits = c(0,24))+
  theme_bw() +
  theme(legend.position=c(0.9, 0.9),
        legend.title = element_blank(),
        legend.text = element_markdown(),
        legend.background = element_blank(),
        axis.title.y = element_blank())+
  guides(fill = guide_legend(override.aes = list(colour = col_vec,alpha=1),drop=FALSE))+
  ggtitle("HDM")+
  ylim(y_limit)-> p_hdm


library(patchwork)
p_pbs + p_hdm

ggsave(p_pbs + p_hdm,filename="~/Amlan/Flexivent_Jan25/p_Max_final_v3.png",width=8,height=5)
