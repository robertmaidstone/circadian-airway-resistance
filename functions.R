

# function to fit dose response curve and manipulate parameters -------------------------------------

dr_fit <- function(LFdata) {
  LF_data$Group <- LF_data$Sample
  
  
  model <- drm(
    Value ~ Mch_conc,
    data = LF_data,
    curveid = Group,
    fct = LL.4(names = c("Slope", "Lower", "Upper", "EC50"))
  )
  
  params <- coef(model)
  param_df <- as.data.frame(params)
  param_df$Group <- rownames(param_df)
  param_df <- param_df %>%
    separate(Group,
             into = c("Parameter", "Sample"),
             sep = ":")
  
  param_df %>% merge(LF_data %>% dplyr::select(Sample, Genotype, Treatment, ZT) %>% unique,
                     by = "Sample") -> param_formodel
  
  param_formodel$ZT <- factor(param_formodel$ZT)
  param_formodel$Treatment <- factor(param_formodel$Treatment)
  param_formodel$Genotype <- factor(param_formodel$Genotype)
  return(param_formodel)
}


# function to run anova on dr_fit output ----------------------------------

dr_anova <- function(param_formodel, parameter, formula_string){
  # Check parameter exists
  if(!parameter %in% param_formodel$Parameter){
    stop("Error: unknown parameter")
  }
  # Convert the character string into a formula
  model_formula <- as.formula(formula_string)
  # WT model
  aov_out_WT <- summary(
    aov(model_formula,
        data = param_formodel %>% 
          dplyr::filter(Genotype == "WT", Parameter == parameter))
  )
  # KO model
  aov_out_KO <- summary(
    aov(model_formula,
        data = param_formodel %>% 
          dplyr::filter(Genotype == "KO", Parameter == parameter))
  )
  # Return p-values in a tidy data frame
  out <- data.frame(
    Variable = row.names(aov_out_WT[[1]]),
    WT = aov_out_WT[[1]]$`Pr(>F)`,
    KO = aov_out_KO[[1]]$`Pr(>F)`
  )
  # Remove residual row
  out <- out[-nrow(out), ]
  return(out)
}


# p to stars --------------------------------------------------------------

p_to_stars <- function(p){
  dplyr::case_when(
    p < 0.0001 ~ "****",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
}


# extract sig text --------------------------------------------------------

extract_sig_text <- function(df, Gen, prefix = "Upper") {
  
  sigp <- df %>%
    dplyr::select(Variable, all_of(Gen)) %>%
    dplyr::filter(.data[[Gen]] < 0.05)
  
  if (nrow(sigp) == 0) return("")
  
  sigp %>%
    dplyr::mutate(
      stars = p_to_stars(.data[[Gen]]),
      text  = paste0(prefix, Variable, " ", stars)
    ) %>%
    dplyr::pull(text) %>%
    paste(collapse = "\n")
}

# dose response curve plot ------------------------------------------------

dr_plot <- function(LF_data,anova_pvals_slope,anova_pvals_upper){

anova_pvals_slope$Variable <- c("ZT","Treatment","Interaction")
anova_pvals_upper$Variable <- c("ZT","Treatment","Interaction")
  
gen_labels <- c("WT","CCSP-Reverbα KO")

LF_data %>% 
  group_by(ZT,Genotype,Treatment,Mch_conc) %>%
  mutate(Med_Value=mean(Value,na.rm=T)) %>%
  dplyr::select(-Sample,-Value) %>%
  unique %>% ungroup -> t_data

for(Gen in c("WT","KO")){
  upper_text <- extract_sig_text(anova_pvals_upper, Gen, prefix = "Upper ")
  slope_text <- extract_sig_text(anova_pvals_slope, Gen, prefix = "Slope ")
  
  sig_text <- paste(upper_text, slope_text, sep = "\n")
  
t_data %>%
  filter(Genotype==Gen) %>%
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
  ggtitle(gen_labels[which(Gen==c("WT","KO"))])-> p1
 assign(paste0("p_",Gen),value = p1)
}

  return(p_WT|p_KO)
}


fit_nls_model <- function(df) {
  nls(
    AUC ~ I + A * sin(2*pi/24*ZT + phi) +
      (Genotype=="KO")*(I1 + A1 * sin(2*pi/24*ZT + phi1)),
    data = df,
    start = list(A = 1, phi = 0, A1 = 0.1, phi1 = 0.1, I = 0, I1 = 0)
  )
}

fit_lrtest <- function(df) {
  nls_model <- nls(
    AUC ~ I + A * sin(2*pi/24*ZT + phi),
    data = df,
    start = list(A = 1, phi = 0, I = 0)
  )
  
  lm_model <- lm(AUC ~ 1, data = df)
  
  # return the LR test p-value
  lmtest::lrtest(nls_model, lm_model)$`Pr(>Chisq)`[2]
}

plot_rhy_funcs <- function(df,predict_values,annot_pvals,sig_line_pvals,Tr){
  col_vec <- c("#0072B2","#E69F00")
  y_limit <- c(-6,50)
  gen_labels <- c("CCSP-Reverbα WT","CCSP-Reverbα KO")
  
  sig_text <- extract_sig_text(annot_pvals, Tr, prefix = "")
  
  pred_df <- merge(
    predict_values,
    sig_line_pvals %>%
      pivot_longer(
        cols = -1,
        names_to = "Treatment",
        values_to = "Linetype"
      ),
    by = c("Treatment", "Genotype"),
    all.x = TRUE
  ) %>%
    mutate(Linetype = Linetype < 0.05)
  
  sum_data %>%
    mutate(Genotype=factor(Genotype,levels=c("WT","KO"),labels=gen_labels)) %>%
    filter(Treatment==Tr) %>%
    ggplot(aes(x = ZT, y = AUC, colour = Genotype)) +
    geom_point() +
    geom_ribbon(aes(ymax=Predicted_Response_UCI,ymin=Predicted_Response_LCI,y=NULL,fill=Genotype,colour=NULL),
                data=pred_df %>%filter(Treatment==Tr) %>% mutate(Genotype=factor(Genotype,levels=c("WT","KO"),
                                                                                 labels=gen_labels)),size=1,alpha=0.2) +
    geom_line(aes(y=Predicted_Response,linetype=Linetype),
              data=pred_df %>%filter(Treatment==Tr) %>% mutate(Genotype=factor(Genotype,levels=c("WT","KO"),
                                                                               labels=gen_labels)),size=1) +
    ylab("AUC Airway Resistance R<sub>rs</sub>(cm.H<sub>2</sub>O.s.ml<sup>-1</sup>)") + scale_colour_manual(values=col_vec) + scale_fill_manual(values=col_vec) +
    scale_x_continuous(breaks=c(0,6,12,18,24))+
    scale_linetype_manual(values=c("dashed","solid"))+
    theme_bw() +
    theme(legend.position = "none",
          axis.title.y = element_markdown()) +
    ggtitle(Tr)+
    annotate("text", x = 1, y = max(y_limit), label = sig_text, hjust = 0, vjust = 1, size = 4)+
    ylim(y_limit)-> p1
  return(p1)
}


predict_values<- merge(predict_values,list_out$loglik_pvals%>% 
                         pivot_longer(cols = -1,names_to = "Treatment",values_to = "Linetype"),
                       by=c("Treatment","Genotype"),all.x=TRUE) %>%
  mutate(Linetype=Linetype<0.05)


rhy_plot<-function(LF_data){
  treatments <- c("PBS", "HDM")
  genotypes  <- c("WT", "KO")
  list_out <- list()
  
  ## data manip
  
  LF_data %>% 
    mutate(Value=log10(Value)) %>%
    group_by(Sample,ZT,Genotype,Treatment) %>%
    mutate(Max_Value=max(Value,na.rm=T)) %>%
    mutate(AUC=sum(diff(Mch_conc)*(Value[-1]+Value[-5])/2)) %>%
    dplyr::select(Sample,Max_Value,AUC) %>%
    #mutate(AUC=Max_Value) %>%
    unique %>% ungroup -> sum_data
  
  ## nls model fitting
  pval_mat <- lapply(treatments, function(trt) {
    df_sub <- sum_data %>% dplyr::filter(Treatment == trt)
    model  <- fit_nls_model(df_sub)
    summary(model)$parameters[c("A1","phi1","I1"), "Pr(>|t|)"]
  })
  list_out[["nls_pvals"]] <- do.call(cbind, pval_mat)
  colnames(list_out[["nls_pvals"]]) <- treatments
  list_out[["nls_pvals"]]<-as.data.frame(list_out[["nls_pvals"]]) %>% mutate(Variable=rownames(list_out[["nls_pvals"]]))
  
  ## log likelihood test comparing rhythmic to constant
  results <- expand.grid(Genotype = genotypes,
                         Treatment = treatments,
                         stringsAsFactors = FALSE)
  
  results$p_value <- mapply(function(g, t) {
    df_sub <- sum_data %>% 
      dplyr::filter(Genotype == g, Treatment == t)
    
    fit_lrtest(df_sub)
  }, results$Genotype, results$Treatment)
  results %>% pivot_wider(names_from = Treatment,values_from = p_value) %>%
    as_data_frame()-> results 
  list_out[["loglik_pvals"]] <- results
  
  ## plotting
  
  lm(AUC~1+Genotype*Treatment*sin(2*pi/24*ZT) + Genotype*Treatment*cos(2*pi/24*ZT)+Genotype*Treatment,data=sum_data) -> lm_sGAW
  
  predict_values <- expand.grid(
    ZT = seq(from=0,to=24,length.out=12),
    Treatment = unique(sum_data$Treatment),
    Genotype = unique(sum_data$Genotype)
    #Animal.ID = unique(sum_data_sGAW$Animal.ID)
  )
  
  pred_resp <- predict.lm(lm_sGAW, newdata = predict_values, se.fit=TRUE, interval="confidence", level=0.95)
  predict_values$Predicted_Response <- pred_resp$fit[,1]
  predict_values$Predicted_Response_LCI <- pred_resp$fit[,2]
  predict_values$Predicted_Response_UCI <- pred_resp$fit[,3]
  
  ##
  list_out[["plot_pbs"]] <- plot_rhy_funcs(sum_data,predict_values,annot_pvals = list_out[["nls_pvals"]],sig_line_pvals = list_out[["loglik_pvals"]],"PBS")
  list_out[["plot_hdm"]] <- plot_rhy_funcs(sum_data,predict_values,annot_pvals = list_out[["nls_pvals"]],sig_line_pvals = list_out[["loglik_pvals"]],"HDM")
  
  return(list_out)
}