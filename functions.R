

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
      text  = paste0(prefix, " ", Variable, " ", stars)
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
  upper_text <- extract_sig_text(anova_pvals_upper, Gen, prefix = "Upper")
  slope_text <- extract_sig_text(anova_pvals_slope, Gen, prefix = "Slope")
  
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
