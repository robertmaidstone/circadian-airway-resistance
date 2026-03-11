

# function to fit dose response curve and manipulate parameters -------------------------------------

dr_fit <- function(LFdata) {
  #LF_data$Group <- LF_data$Sample
  LF_data$Group <- interaction(LF_data$Sample,LF_data$ZT)
  
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
             into = c("Parameter", "Group"),
             sep = ":")
  
  param_df %>% merge(LF_data %>% dplyr::select(Group, Genotype, Treatment, ZT) %>% unique,
                     by = "Group") -> param_formodel
  
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

dr_plot <- function(LF_data,anova_pvals_slope,anova_pvals_upper,y_lim,y_lab){

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
  ylab(y_lab) +
  ylim(y_lim)+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(),
        plot.title = element_markdown(),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_manual(values = c("#0072B2","#E69F00","#D55E00", "#009E73")) +
  annotate("text", x = 1, y = max(y_lim), label = sig_text, hjust = 0, vjust = 1, size = 4)+
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

plot_rhy_funcs <- function(df, predict_values, annot_pvals, sig_line_pvals,
                           Tr, y_axis = TRUE, legend = TRUE, y_lab, y_lim) {
  
  col_vec    <- c("#0072B2", "#E69F00")
  gen_labels <- c("CCSP-Reverbα WT", "CCSP-Reverbα KO")
  y_limit    <- y_lim
  # --- 1. Prepare significance text -----------------------------------------
  sig_text <- extract_sig_text(annot_pvals, Tr, prefix = "")
  # --- 2. Prepare prediction dataframe --------------------------------------
  pred_df <- predict_values %>%
    left_join(
      sig_line_pvals %>%
        pivot_longer(cols = -1, names_to = "Treatment", values_to = "Linetype"),
      by = c("Treatment", "Genotype")
    ) %>%
    mutate(
      Linetype = Linetype < 0.05,
      Genotype = factor(Genotype, levels = c("WT", "KO"), labels = gen_labels)
    )
  # --- 3. Prepare observed data ---------------------------------------------
  df_plot <- df %>%
    filter(Treatment == Tr) %>%
    mutate(Genotype = factor(Genotype, levels = c("WT", "KO"), labels = gen_labels))
  # --- 4. Build the base plot ------------------------------------------------
  p <- ggplot(df_plot, aes(x = ZT, y = AUC, colour = Genotype)) +
    geom_point() +
    geom_ribbon(
      data = pred_df %>% filter(Treatment == Tr),
      aes(ymax = Predicted_Response_UCI,
          ymin = Predicted_Response_LCI,
          y=NULL,
          fill = Genotype,
          colour = NULL),
      alpha = 0.2
    ) +
    geom_line(
      data = pred_df %>% filter(Treatment == Tr),
      aes(y = Predicted_Response, linetype = Linetype),
      linewidth = 1
    ) +
    scale_colour_manual(values = col_vec) +
    scale_fill_manual(values = col_vec) +
    scale_linetype_manual(values = c("dashed", "solid"), guide = "none") +
    scale_x_continuous(breaks = c(0, 6, 12, 18, 24)) +
    guides(linetype = "none") +
    ylab(y_lab) +
    ggtitle(Tr) +
    annotate("text", x = 1, y = max(y_limit), label = sig_text,
             hjust = 0, vjust = 1, size = 4) +
    ylim(y_limit) +
    theme_bw() +
    theme(
      axis.title.y = element_markdown(),
      legend.position = "none"
    )
  # --- 5. Legend logic -------------------------------------------------------
  if (legend) {
    p <- p + theme(
      legend.position = c(0.75, 0.9),
      legend.title = element_blank(),
      legend.text = element_markdown(),
      legend.background = element_blank()
    )
  }
  # --- 6. Y-axis logic -------------------------------------------------------
  if (!y_axis) {
    p <- p + theme(axis.title.y = element_blank())
  }
  return(p)
}


rhy_plot<-function(LF_data,Type,y_lim){
  treatments <- c("PBS", "HDM")
  genotypes  <- c("WT", "KO")
  list_out <- list()
  
  ## data manip
  
  LF_data %>% 
    mutate(Value=log10(Value)) %>%
    arrange(across(all_of(c("Sample", "ZT", "Genotype", "Treatment"))), Mch_conc) %>%
    group_by(Sample,ZT,Genotype,Treatment) %>%
    mutate(Max_Value=max(Value,na.rm=T)) %>%
    mutate(Min_Value=min(Value,na.rm=T)) %>%
    mutate(AUC = sum(diff(Mch_conc) * (Value[-1] + Value[-length(Value)]) / 2)) %>%
    #mutate(AUC = sum(diff(Mch_conc) * (Value[-1] + Value[-4]) / 2,na.rm=T)) %>%
    dplyr::select(Sample,Max_Value,Min_Value,AUC) %>%
    distinct %>% ungroup -> sum_data
  if(Type=="Max"){
    sum_data <- sum_data %>% mutate(AUC=Max_Value)
    y_lab <- "Max Airway Resistance R<sub>rs</sub>(cm.H<sub>2</sub>O.s.ml<sup>-1</sup>)"
  }else if(Type=="Min"){
    sum_data <- sum_data %>% mutate(AUC=Min_Value)
    y_lab<-"Min Airway Resistance R<sub>rs</sub>(cm.H<sub>2</sub>O.s.ml<sup>-1</sup>)"
  }else if(Type=="AUC"){
    y_lab<-"AUC Airway Resistance R<sub>rs</sub>(cm.H<sub>2</sub>O.s.ml<sup>-1</sup>)"
  }else{
    stop("Error type not recognised")
  }
  ## nls model fitting
  pval_mat <- lapply(treatments, function(trt) {
    df_sub <- sum_data %>% dplyr::filter(Treatment == trt)
    model  <- fit_nls_model(df_sub)
    summary(model)$parameters[c("A1","phi1","I1"), "Pr(>|t|)"]
  })
  list_out[["nls_pvals"]] <- do.call(cbind, pval_mat)
  colnames(list_out[["nls_pvals"]]) <- treatments
  list_out[["nls_pvals"]]<-as.data.frame(list_out[["nls_pvals"]]) %>% mutate(Variable=rownames(list_out[["nls_pvals"]])) %>%
    mutate(Variable= dplyr::case_when(
      Variable =="A1" ~ "Amplitude",
      Variable =="phi1" ~ "Phase",
      Variable =="I1" ~ "Mesor"
    ))
  
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
    as.data.frame()-> results 
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
  list_out[["plot_pbs"]] <- plot_rhy_funcs(sum_data,predict_values,
                                           annot_pvals = list_out[["nls_pvals"]],sig_line_pvals = list_out[["loglik_pvals"]],
                                           "PBS",legend=FALSE,y_lab=y_lab,y_lim=y_lim)
  list_out[["plot_hdm"]] <- plot_rhy_funcs(sum_data,predict_values,
                                           annot_pvals = list_out[["nls_pvals"]],sig_line_pvals = list_out[["loglik_pvals"]],
                                           "HDM",y_axis=FALSE,y_lab=y_lab,y_lim=y_lim)
  list_out$combined <- list_out$plot_pbs + list_out$plot_hdm
  return(list_out)
}
