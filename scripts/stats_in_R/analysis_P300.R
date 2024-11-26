#task EEG analysis

# preparattions ----
#load libraries
{
  library(readxl)
  library(lmerTest) # for gmler()
  library(car)      # for Anova()
  library(sjPlot)   # for plot_model()
  library(ggplot2)
}

#load and format data
setwd("/data/p_02191/Analysis/Nadine/final_scripts_and_data/") #set working directory
{
per_trial_data <-  read_excel("data/for_analysis/all_data_per_trial_P300.xlsx") #load trial based data
all_data = read_excel("data/for_analysis/all_data_behavioral.xlsx") #load rest of data to get at covariates and DFA and 1/f slope

drop = c("condition","accuracy","RT") #exclude condition, accuracy and RT variable to not confuse the merge
all_data = all_data[ , !(names(all_data) %in% drop)]
all_data = unique(all_data)

data = merge(all_data,per_trial_data, by = "ID") #merge datasets to have trial based data

#format data
{
  data$zWM_tired = scale(data$task_tiredness)
  data$zWM_conc = scale(data$task_concentration)
  data$zIQ = scale(as.numeric(data$IQ))
  data$zAge = scale(as.numeric(data$Age))
  data$zDFS = scale(as.numeric(data$DFS))
  data$zRelAlphaPow_rest = scale(as.numeric(data$relAlphaPow_rest))
  data$zBIS = scale(as.numeric(data$BIS))
  data$zBAS = scale(as.numeric(data$BAS))
  data$zBMI = scale(as.numeric(data$BMI))
  data$zDSBack = scale(as.numeric(data$DSBack))
  data$Gender = as.factor(data$Gender)
  
  data <- data[data$correct != "unknown", ] #remove misses
  data = data[data$reaction_time > 0.2,] #remove false alarm
  data$P3_scale = data$P3b*1000000 #scale to microvolt
}

  #get averaged DFA & 1/f from behavioral analysis
  {  
    result_matrix_DFA <- read_excel("data/for_analysis/result_matrix_DFA.xlsx")
    filtered_electrodes <- result_matrix_DFA[result_matrix_DFA$sig == "***" & result_matrix_DFA$cluster == 2,"elec"]
    eleclist_sig <- as.list(filtered_electrodes)
    eleclist_sig <- unlist(eleclist_sig)
    eleclist_sig_DFA = eleclist_sig[!is.na(eleclist_sig)]
    # Calculate the row means for the selected columns
    data$DFA <- rowMeans(data[, eleclist_sig_DFA])
    
    result_matrix_oneF_AA <- read_excel("data/for_analysis/result_matrix_oneF.xlsx")
    filtered_electrodes <- result_matrix_oneF_AA[result_matrix_oneF_AA$sig == "***","elec"]
    eleclist_sig <- as.list(filtered_electrodes)
    eleclist_sig <- unlist(eleclist_sig)
    eleclist_sig_oneF = eleclist_sig[!is.na(eleclist_sig)] #because both clusters show same behavioral effect, we pooled them together 
    # Calculate the row means for the selected columns
    data$oneF <- rowMeans(data[, eleclist_sig_oneF])
  }
  subset = data[data$condition == 'ignore' | data$condition == 'update',]   #subset data to look at conditions of interest only
  subset = as.data.frame(subset)
  subset$correct = as.factor(subset$correct)
  subset$condition = as.factor(subset$condition)
}  
# Set contrasts
#contrasts(subset$condition) <- matrix(c(-1, 1), ncol = 1)
# Check the new contrasts to confirm
contrasts(subset$condition) #needs to be 0 and 1 for this analysis, because otherwise there will be convergence problems

# stats P3 on behavior ----
{
  model_P3b_behav <- glmer(correct ~ condition*P3_scale + (1|ID),
                           family = "binomial",
                           nAGQ = 1,
                           data = subset,
                           control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e6)))
  Anova(model_P3b_behav, type = 3) #condition:P3b_scale   2.9538  1    0.08567 .
  summary(model_P3b_behav)
  #plot_model(model_P3b_behav, type = "pred", terms = c("P3b_scale","condition"), show.data = F) + theme_classic()
  

  #with control variables
  {
    model_P3b_behav_ctrl <- glmer(correct ~ condition*P3_scale + AAratio + zWM_tired + zWM_conc + zIQ + zAge + zDFS + zRelAlphaPow_rest + zBIS + zBAS + zBMI + zDSBack + (1|ID),
                                  family = "binomial",
                                  nAGQ = 1,
                                  data = subset,
                                  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e7)))
    Anova(model_P3b_behav_ctrl, type = 3) #condition:DFA:AAratio 58.0683  1  2.532e-14 ***
    plot_model(model_P3b_behav_ctrl, type = "pred", terms = c('P3_scale','condition'))
    #with trial numb it does not converge
    
    #write table to word
    {
      library(officer)
      library(flextable)
      library(broom)
      library(dplyr)
      
      anova_results = Anova(model_P3b_behav_ctrl, type = 3)
      anova_df <- broom::tidy(anova_results)
      anova_df <- anova_df %>%
        mutate(across(where(is.numeric), ~ round(., 3)))
      ft <- flextable(anova_df)
      doc <- read_docx()
      doc <- body_add_flextable(doc, value = ft)
      print(doc, target = "~/Work/PhD/WORMCRI/Results/results_P300_cond_ctrl.docx")
      
    }
  }
 
  posthoc <- emtrends(model_P3b_behav_ctrl, ~ condition, var = "P3b_scale")
  posthoc
  #contrast(posthoc, method = "pairwise")
  
  #post hoc
  {
  subset_ign = subset[subset$condition=='ignore',]
  subset_upt = subset[subset$condition=='update',]
  
  model_P3b_behav_ign <- glmer(correct ~ P3b_scale + AAratio + zWM_tired + zWM_conc + zIQ + zAge + zDFS + zRelAlphaPow_rest + zBIS + zBAS + zBMI + zDSBack + (1|ID),
                               family = "binomial",
                               nAGQ = 1,
                               data = subset_ign,
                               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e7)))
  summary(model_P3b_behav_ign)
  
  model_P3b_behav_upt <- glmer(correct ~ P3b_scale + AAratio + zWM_tired + zWM_conc + zIQ + zAge + zDFS + zRelAlphaPow_rest + zBIS + zBAS + zBMI + zDSBack + (1|ID),
                               family = "binomial",
                               nAGQ = 1,
                               data = subset_upt,
                               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e7)))
  summary(model_P3b_behav_upt)
  }
  
  #plot 
  { 
    #keep = c("ID","correct","condition","P3b_scale")
    keep = c("ID","correct","condition","P3_scale","AAratio","zWM_tired","zWM_conc","zIQ","zAge","zDFS","zRelAlphaPow_rest","zBIS","zBAS","zBMI","zDSBack")
    plotdata = subset[ , (names(subset) %in% keep)]
    plotdata = na.omit(plotdata)
    plotdata$fitted = fitted(model_P3b_behav_ctrl)
    plotdata = unique(plotdata)
   
    ggplot(data = plotdata, aes(P3_scale, fitted, colour = condition, shape = condition, fill = condition)) +
      geom_smooth(method = "lm") +
      geom_point(alpha = .4, size = 0.5) +
      scale_color_manual(values = c("#ffc107ff","#8aaedcff")) +
      scale_fill_manual(values = c("#ffc107ff","#8aaedcff")) +
      scale_shape_manual(values = c(24, 22)) +  # Set filled triangle (24) and filled square (22) shapes for groups
      theme_classic() +
      theme(axis.text.x = element_text(size = 14, color = "black")) +
      theme(axis.text.y = element_text(size = 14, color = "black")) +
      labs(y = "probability correct") +
      labs(x = "P300 amplitude (µV)") +
      theme(axis.title.y = element_text(margin = margin(r = 15), size = 14)) +
      theme(axis.title.x = element_text(margin = margin(r = 15), size = 14)) +
      theme(strip.background = element_blank(),  # Remove the background of facet labels
            strip.text = element_text(size = 12)) +
      guides(colour = guide_legend(title = "condition", 
                                   title.position = "left", 
                                   title.hjust = 0,
                                   label.position = "right",
                                   nrow = 1), 
             shape = guide_legend(title = "condition", 
                                  title.position = "left", 
                                  title.hjust = 0,
                                  label.position = "right",
                                  nrow = 1))+
      theme(legend.position = c(1, 1),  # Position the legend in the top-right corner
            legend.justification = c(1, 1),
            legend.text = element_text(size = 10.5))  # Anchor the legend to the top-right corner
  }
}

#is there a difference in amplitude per condition?
{
model_P3_diff = lmer(P3_scale ~ condition + (1+condition|ID), data = subset)
summary(model_P3_diff)

#mean and sd per condition
mean_sd_p3b_per_condition <- subset %>%
  group_by(condition) %>%  # Group data by the 'condition' column
  summarize(
    mean_p3b = mean(P3_scale, na.rm = TRUE),  # Calculate mean, ignoring NA values
    sd_p3b = sd(P3_scale, na.rm = TRUE)        # Calculate standard deviation, ignoring NA values
  )


#plot
{
  subset_w0_out = subset[subset$P3b > min(subset$P3b),]
  ggplot(subset_w0_out, aes(x = condition, y = P3b*1000000, fill = condition, color = condition, shape = condition)) +
  geom_violin(trim = TRUE, alpha = 0.8) +  # Draw violin plots 
  geom_point(alpha = 0.3, size = 2, position = position_jitter(width = 0.2)) +  # Add jittered points
  labs(x = "Condition", y = "P3b (µV)") +
  theme_minimal() +
  scale_color_manual(values = c("#ffc107ff", "#8aaedcff")) +  # Manual color scale for points and outlines
  scale_fill_manual(values = c("#ffc107ff", "#8aaedcff")) +  # Manual fill color scale for violins and boxplots
  scale_shape_manual(values = c(15, 18)) +  # Assign shapes to conditions
  geom_boxplot(width = 0.2, fill = "white", color = "gray") +   
  theme(
    axis.text.x = element_text(size = 18, colour = "black"),  # Adjust x-axis tick label size
    axis.text.y = element_text(size = 18, colour = "black"),  # Adjust y-axis tick label size
    axis.title.x = element_text(size = 18),  # Adjust x-axis label size
    axis.title.y = element_text(size = 18)   # Adjust y-axis label size 
  )
}
}


#### do DFA or oneF*AA influence P300? ----
model_P3_DFA <- lmer(P3_scale ~ DFA*condition + (1|ID), data = subset)
summary(model_P3_DFA)

model_P3_oneF_AA = lmer(P3_scale ~ oneF*AAratio*condition + (1|ID), data = subset)
summary(model_P3_oneF_AA)
plot_model(model_P3_oneF_AA, type = "pred", terms = c("AAratio","condition"))

#write table to word
{
  anova_results = Anova(model_P3_DFA, type = 3)
  anova_df <- broom::tidy(anova_results)
  anova_df <- anova_df %>%
    mutate(across(where(is.numeric), ~ round(., 3)))
  ft <- flextable(anova_df)
  doc <- read_docx()
  doc <- body_add_flextable(doc, value = ft)
  print(doc, target = "~/Work/PhD/WORMCRI/Results/results_P300_DFA.docx")
  
}

#only for updating? since only here P300 amplitude seemed to be behaviorally relevant?
{
  subset_upt = subset[subset$condition == 'update',]
  model_P3_DFA_upt <- lmer(P3_scale ~ DFA + (1|ID), data = subset_upt)
  summary(model_P3_DFA_upt)
  
  model_P3_oneF_upt <- lmer(P3_scale ~ oneF+AAratio + (1|ID), data = subset_upt)
  summary(model_P3_oneF_upt)
  plot_model(model_P3_oneF_upt, type = "pred", terms = c("AAratio"))
  
  model_P3_AA_upt <- lmer(P3_scale ~ AAratio + (1|ID), data = subset_upt)
  summary(model_P3_AA_upt)
  plot_model(model_P3_AA_upt, type = "pred", terms = c("AAratio"))
}
