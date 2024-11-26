#load libraries
{
  library(readxl)
  library(lme4)
  #library(lmerTest)
  library(ggplot2)
  library(sjPlot)   #for plot_model()
  library(car)      #for Anova()
  library(lsmeans)  #for posthoc
  library(dplyr)
  library(openxlsx)
  library(igraph)   #for graph function to find clusters
  
}
data <- read_excel("~/Work/PhD/WORMCRI/Data/all_data.xlsx")
DFA_Slope_allSubs <- read_excel("~/Work/PhD/WORMCRI/Data/source/DFA_Slope_allSubs.xlsx")

#format data
{
  data$zWM_tired = scale(data$task_tiredness)
  data$zWM_conc = scale(data$task_concentration)
  data$zBMI = scale(data$BMI)
  data$zIQ = scale(as.numeric(data$IQ))
  data$zAge = scale(as.numeric(data$Age))
  data$zDFS = scale(as.numeric(data$DFS))
  data$zRelAlphaPow_rest = scale(as.numeric(data$relAlphaPow_rest))
  data$zBIS = scale(as.numeric(data$BIS))
  data$zBAS = scale(as.numeric(data$BAS))
  data$zDSBack = scale(as.numeric(data$DSBack))
  data$Gender = as.factor(data$Gender)
  data$condition <- as.factor(data$condition)
  
  median_AA <- median(data$AAratio, na.rm = TRUE)
  data$AA_group <- ifelse(data$AAratio <= median_AA, "low", "high")
  
  data <- data %>% #rename m2 and m1
    mutate(condition = recode(condition,
                              'm1' = 'control_long',
                              'm2' = 'control_short'))
}

#get averaged DFA & 1/f from behavioral analysis
{  
  result_matrix_DFA <- read_excel("~/Work/PhD/WORMCRI/Data/source/results_matrix_DFA_source.xlsx")
  filtered_electrodes <- result_matrix_DFA[result_matrix_DFA$V3 == "***","V1"]
  eleclist_sig <- as.list(filtered_electrodes)
  eleclist_sig <- unlist(eleclist_sig)
  eleclist_sig_DFA = eleclist_sig[!is.na(eleclist_sig)]
  # Calculate the row means for the selected columns
  # Calculating row means and converting to dataframe
  DFA <- as.data.frame(rowMeans(DFA_Slope_allSubs[, eleclist_sig_DFA]))
  DFA[,2] <- DFA_Slope_allSubs[, 1]
  names(DFA) <- c("source_DFA", "ID")  
  test = merge(data, DFA, by = "ID")

  
  result_matrix_oneF_AA <- read_excel("~/Work/PhD/WORMCRI/Data/source/results_matrix_oneF_source.xlsx")
  filtered_electrodes <- result_matrix_oneF_AA[result_matrix_oneF_AA$V5 == "***","V1"]
  eleclist_sig <- as.list(filtered_electrodes)
  eleclist_sig <- unlist(eleclist_sig)
  eleclist_sig_oneF = eleclist_sig[!is.na(eleclist_sig)] #because both clusters show same behavioral effect, we pooled them together 
  # Calculate the row means for the selected columns
  oneF <- as.data.frame(rowMeans(DFA_Slope_allSubs[, eleclist_sig_oneF]))
  oneF[,2] <- DFA_Slope_allSubs[, 1]
  names(oneF) <- c("source_oneF", "ID")  
  data = merge(test, oneF, by = "ID")
}

model_DFA_source <- lmer(accuracy ~ condition * source_DFA + (1|ID), data = data)
Anova(model_DFA_source, type = 3)
summary(model_DFA_source)

#plot 
{ 
  keep = c("ID","correct","condition","source_DFA")
  plotdata = data[ , (names(data) %in% keep)]
  plotdata = na.omit(plotdata)
  plotdata$fitted = fitted(model_DFA_source)
  plotdata = unique(plotdata)
  
  plotdata$condition <- factor(data$condition, levels = c("ignore", "control_long", "update", "control_short"))
  
  
  ggplot(data = plotdata, aes(source_DFA, fitted, colour = condition, shape = condition, fill = condition)) +
    geom_smooth(method = "lm") +
    geom_point(alpha = .6, size = 3) +
    scale_color_manual(values = c("#ffc107ff","grey80","#8aaedcff","grey40")) +
    scale_fill_manual(values = c("#ffc107ff","grey80","#8aaedcff","grey40")) +
    scale_shape_manual(values = c(15,16,17,18)) +  # Set filled triangle (24) and filled square (22) shapes for groups
    theme_classic() +
    facet_grid(.~ condition)+
    theme(axis.text.x = element_text(size = 14, color = "black", family = "arial")) +
    theme(axis.text.y = element_text(size = 14, color = "black", family = "arial")) +
    labs(y = "predicted accuracy") +
    labs(x = "average DFA (source level)") +
    theme(axis.title.y = element_text(margin = margin(r = 15), size = 14, family = "Arial")) +
    theme(axis.title.x = element_text(margin = margin(r = 15), size = 14, family = "Arial")) +
    theme(strip.background = element_blank(),  # Remove the background of facet labels
          strip.text = element_text(size = 12)) +
    scale_x_continuous(breaks = c(0.65, 0.95))
}

#posthoc
{
  subset_ign = data[data$condition == "ignore",]
  model <- lm(accuracy ~ source_DFA, data = subset_ign)
  summary(model)
  #             Estimate   Std. Error  t value  Pr(>|t|)    
  # (Intercept)   1.05944    0.07111   14.899   <2e-16 ***
  # source_DFA   -0.24754    0.10046   -2.464   0.0162 * 
  
  
  subset_upt = data[data$condition == "update",]
  model <- lm(accuracy ~ source_DFA, data = subset_upt)
  summary(model)
  #             Estimate   Std. Error  t value  Pr(>|t|)    
  # (Intercept)   0.93003    0.07291   12.756   <2e-16 ***
  # source_DFA    0.07742    0.10305   0.751    0.455
  
  
  subset_m1 = data[data$condition == "control_long",]
  model <- lm(accuracy ~ source_DFA, data = subset_m1)
  summary(model)
  #              Estimate   Std. Error  t value  Pr(>|t|)    
  # (Intercept)   1.10408    0.05780    19.101   < 2e-16 ***
  # source_DFA   -0.27913    0.08044    -3.47    0.000889 *** 
  
  
  subset_m2 = data[data$condition == "control_short",]
  model <- lm(accuracy ~ source_DFA, data = subset_m2)
  summary(model)
  #               Estimate   Std. Error  t value  Pr(>|t|)    
  # (Intercept)    1.04081    0.03908    26.631   <2e-16 ***
  # source_DFA    -0.08663    0.05639    -1.536    0.129 
}

#### oneF ----

model_oneF_source <- lmer(accuracy ~ condition * source_oneF * AAratio + (1|ID), data = data)
Anova(model_oneF_source, type = 3)
summary(model_oneF_source)

#plot 
{ 
  keep = c("ID","accuracy","condition","source_oneF", "AA_group")
  plotdata = data[ , (names(data) %in% keep)]
  plotdata = na.omit(plotdata)
  plotdata$fitted = fitted(model_oneF_source)
  plotdata = unique(plotdata)
  
  plotdata$condition <- factor(plotdata$condition, levels = c("ignore", "control_long", "update", "control_short"))
  
  
  ggplot(data = plotdata, aes(source_oneF, fitted, colour = AA_group, shape = AA_group, fill = AA_group)) +
    geom_smooth(method = "lm") +
    geom_point(alpha = .4, size = 2) +
    scale_color_manual(values = c("#fa6640ff","#3ea68fff")) +
    scale_fill_manual(values = c("#fa6640ff","#3ea68fff")) +
    scale_shape_manual(values = c(25, 15)) +  # Set filled triangle (24) and filled square (22) shapes for groups
    theme_classic() +
    facet_grid(.~ condition) +
    theme(axis.text.x = element_text(size = 14, color = "black", family = "arial")) +
    theme(axis.text.y = element_text(size = 14, color = "black", family = "arial")) +
    labs(y = "predicted accuracy") +
    labs(x = "average 1/f slope (source level)") +
    theme(axis.title.y = element_text(margin = margin(r = 15), size = 14, family = "Arial")) +
    theme(axis.title.x = element_text(margin = margin(r = 15), size = 14, family = "Arial")) +
    theme(strip.background = element_blank(),  # Remove the background of facet labels
          strip.text = element_text(size = 12)) 
}
