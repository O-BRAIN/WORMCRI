#load libraries
{
  library(readxl)
  library(lme4)
  library(lmerTest)
  library(ggplot2)
  library(sjPlot)   #for plot_model()
  library(car)      #for Anova()
  library(lsmeans)  #for posthoc
  library(dplyr)
  library(openxlsx)
  library(igraph)   #for graph function to find clusters
  
}
setwd("/data/p_02191/Analysis/Nadine/final_scripts_and_data/data/for_analysis/")
data <- read_excel("all_data_behavioral.xlsx")

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

#read in DFA and oneF (for when you want to run the script without having to run the electrode loop)
{
  result_matrix_oneF <- read_excel("result_matrix_oneF.xlsx")
  filtered_electrodes <- result_matrix_oneF[result_matrix_oneF$sig == "***","elec"]
  eleclist_sig <- as.list(filtered_electrodes)
  eleclist_sig <- unlist(eleclist_sig)
  eleclist_sig_oneF = eleclist_sig[!is.na(eleclist_sig)]
  data$oneF = rowMeans(data[, eleclist_sig_oneF])
  
  result_matrix_DFA <- read_excel("result_matrix_DFA.xlsx")
  filtered_electrodes <- result_matrix_DFA[result_matrix_DFA$sig == "***" & result_matrix_DFA$cluster == 2,"elec"]
  eleclist_sig <- as.list(filtered_electrodes)
  eleclist_sig <- unlist(eleclist_sig)
  eleclist_sig_DFA = eleclist_sig[!is.na(eleclist_sig)]
  data$DFA <- rowMeans(data[, eleclist_sig_DFA])
}

#sample characteristics
{length(unique(data$ID)) #sample size
sum(data$Gender == 1)/4 #numer of females (diveded by 4, because 4 conditions)

mean(as.numeric(data$Age))
sd(as.numeric(data$Age))
min(as.numeric(data$Age))
max(as.numeric(data$Age))

mean(as.numeric(data$IQ))
sd(as.numeric(data$IQ))
min(as.numeric(data$IQ))
max(as.numeric(data$IQ))

mean(as.numeric(data$BMI))
sd(as.numeric(data$BMI))
min(as.numeric(data$BMI))
max(as.numeric(data$BMI))


}

#set contrasts
#options(contrasts = c("contr.treatment","contr.poly"))
contrasts(data$condition) #needs to be treatment, otherwise convergence problems

# create function to check if two electrodes are adjacent
is_adjacent <- function(electrode1, electrode2, adj_matrix) {
  # Extract the row corresponding to the first electrode
  row <- adj_matrix[electrode1, ]
  
  # Extract the value for the column of the second electrode
  adjacency_value <- row[electrode2]
  
  # Return TRUE if adjacent (1), FALSE otherwise
  return(adjacency_value == 1)
}

#stats behavior per condition ----

model <- lmer(accuracy ~ condition + (1|ID), data = data)
Anova(model, type = 3)
summary(model)

emm <- emmeans(model, pairwise ~ condition)
pairs <- emm$contrasts
summary(pairs, adjust = "tukey")  # Adjust using Tukey method

#summary stats
{summary_stats <- data %>%
  group_by(condition) %>%  # Group data by the 'condition' column
  summarize(
    mean_accuracy = mean(accuracy, na.rm = TRUE),  # Calculate mean, ignoring NA values
    sd_accuracy = sd(accuracy, na.rm = TRUE)       # Calculate standard deviation, ignoring NA values
  )
}


#plot
{
  data$condition <- factor(data$condition, levels = c("ignore", "control_long", "update", "control_short"))
  
  
    ggplot(data, aes(x = condition, y = accuracy, fill = condition, color = condition, shape = condition)) +
    geom_violin(trim = TRUE, alpha = 0.3) +  # Draw violin plots
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  # Add boxplots on top
    geom_point(alpha = 0.5, size = 2, position = position_jitter(width = 0.2)) +  # Add jittered points
    stat_summary(fun = median, geom = "point", shape = 23, size = 3, color = "darkgrey", fill = "black") +  # Add black dots for means
    theme_classic() +
    scale_color_manual(values = c("#ffc107ff", "grey80", "#8aaedcff", "grey40")) +  # Manual color scale for points and outlines
    scale_fill_manual(values = c("#ffc107ff", "grey80", "#8aaedcff", "grey40")) +  # Manual fill color scale for violins and boxplots
    scale_shape_manual(values = c(15, 16, 17, 18)) +  # Assign shapes to conditions
    scale_y_continuous(labels = scales::percent) + 
    theme(
      axis.text.x = element_text(size = 18, colour = "black"),  # Adjust x-axis tick label size
      axis.text.y = element_text(size = 18, colour = "black"),
      axis.title.x = element_text(size = 18),  # Adjust x-axis label size
      axis.title.y = element_text(size = 18)   # Adjust y-axis label size 
    )

}


#stats DFA-behavior ----
{

#loop over all electrodes to find in which electrodes DFA relates to condition-dependent performance
{
# Find the index of 'Fp1' in the column names
start_index <- which(colnames(data) == 'DFA_Fp1')
# Find the index of 'P08' in the column names
end_index <- which(colnames(data) == 'DFA_PO8')

eleclist <- colnames(data)[start_index:end_index]
#initiate result_matrix
result_matrix_DFA <- matrix(nrow = length(eleclist), ncol = 8)
colnames(result_matrix_DFA) <- c("elec", "p_DFA_cond", "Chi_DFA_cond", "sig","p_DFA_cond_after_AA", "p_DFA_cond_AA", "Chi_DFA_cond_AA", "sig_DFA_cond_AA")
count = 1

for (elec in eleclist) {
  # Build the formula with the current ELEC variable
  formula1 <- as.formula(paste("accuracy ~ condition *", elec, "+ (1|ID)"))
  formula2 <- as.formula(paste("accuracy ~ condition *", elec, "* AAratio + (1|ID)"))
  print(paste("Processing elec:", elec))
  # Fit the glmer model
  model1 <- lmer(formula1, data = data)
  model2 <- lmer(formula2, data = data)
  
  model_summary1 <- Anova(model1, type = 3)
  model_summary2 <- Anova(model2, type = 3)
  tvals1 = model_summary1[4,1]
  tvals2 = model_summary2[8,1]
  p_DFA <-  model_summary1$Pr[4]
  p_DFA_AA <- model_summary2$Pr[8]
  p_DFA_after_AA <- model_summary2$Pr[5]

  result_matrix_DFA[count,1] = elec
  result_matrix_DFA[count,2] = p_DFA
  result_matrix_DFA[count,3] = tvals1
  if (p_DFA < 0.05){
    result_matrix_DFA[count,4] = '***'}
  if (p_DFA_after_AA < 0.05){
    result_matrix_DFA[count,5] = '***'}
  result_matrix_DFA[count,6] = p_DFA_AA
  result_matrix_DFA[count,7] = tvals2 #store tvalue for the output of model1 condition*DFA effect (for permutation testing later)
  if (p_DFA_AA < 0.05){
    result_matrix_DFA[count,8] = '***'}
  
  count = count + 1
  
}
write.xlsx(as.data.frame(result_matrix_DFA), file = "~/Work/PhD/WORMCRI/Results/result_matrix_DFA.xlsx")

}
  
#do cluster based permutation testing on DFA results
{
#load and format adjacency matrix
{
adj_matrix <- read_excel("~/Work/PhD/WORMCRI/scripts/adjacency_matrix.xlsx")
adj_matrix_df <- as.data.frame(adj_matrix)

# Change column and rownames names
colnames(adj_matrix_df) <- paste("DFA", colnames(adj_matrix), sep = "_")
rownames(adj_matrix_df) <- adj_matrix_df[,1] # Adjust this if your row names are in a different column
adj_matrix_df <- adj_matrix_df[,-1] # Remove the column now serving as row names, if necessary
rownames(adj_matrix_df) <- paste("DFA", rownames(adj_matrix_df), sep = "_")
}

#perm loop
{
N_perm = 1000
sum_perm_tvalue <- matrix(nrow = N_perm, ncol = 1)
  
for (perm in 1:N_perm){
  print(paste("perm # ", perm))
  #shuffle data
  permdata <- data %>%
  mutate(accuracy = sample(accuracy, size = n(), replace = FALSE))
  
  result_matrix_perm <- matrix(nrow = length(eleclist), ncol = 5)
  count = 1
  
  for (elec in eleclist) {
    # Build the formula with the current ELEC variable
    formula <- as.formula(paste("accuracy ~ condition *", elec, "+ (1|ID)"))
    # Fit the glmer model
    model <- lmer(formula, data = permdata)
    model_summary <- Anova(model, type = 3)
    p_values <- model_summary$Pr
    p_DFA <- p_values[4]
    tvals = summary(model)$coefficients[, "t value"]
    tval_upt = tvals[8]
    
    result_matrix_perm[count,1] = elec
    result_matrix_perm[count,2] = p_DFA
    if (p_DFA < 0.05){
      result_matrix_perm[count,3] = '***'}
    result_matrix_perm[count,4] = tval_upt
    count = count + 1
  }
  
  #collect sum T value for sig. cluster of adjacent electrodes
  result_matrix_perm = as.data.frame(result_matrix_perm)
  sig_elecs <- result_matrix_perm[result_matrix_perm$V3 == "***","V1"] #get significant electrodes
  sig_elecs <- na.omit(sig_elecs)
  
  if (length(sig_elecs) == 0) {
    sum_perm_tvalue[perm,] = 0
  } else {
    
    # Check adjacency for all pairs in the list
    length_sig <- length(sig_elecs)
    adjacency_results <- matrix(NA, nrow = length_sig, ncol = length_sig)
    
    for (i in 1:length_sig) {
      for (j in 1:length_sig) {
        if (i != j) {
          adjacency_results[i, j] <- is_adjacent(sig_elecs[i], sig_elecs[j], adj_matrix_df)
        }
      }
    }
    
    rownames(adjacency_results) <- sig_elecs
    colnames(adjacency_results) <- sig_elecs
    
    #find clusters
    {
      # convert the adjacency_results matrix to a graph
      g <- graph_from_adjacency_matrix(adjacency_results, mode = "undirected", diag = FALSE)
      # find connected components of the graph
      components <- components(g)
      # extract the clusters
      clusters <- split(names(V(g)), components$membership)
    }
    
    #calculate sum of T value of biggest cluster
    {
      cluster_sizes <- sapply(clusters, length)
      biggest_cluster_index <- which.max(cluster_sizes) # find the index of the biggest cluster
      biggest_cluster <- clusters[[biggest_cluster_index]] # select the biggest cluster
      
      x <- result_matrix_perm[result_matrix_perm$V1 %in% biggest_cluster, ]
      sum_current_Tval = sum(as.numeric(x$V4))
    }
    sum_perm_tvalue[perm,] = sum_current_Tval
  }
  
    
}  


# Save the current workspace
save.image(file = "~/Work/PhD/WORMCRI/Results/WS_perm_1000_DFA.RData")}

#compare to summed Tvalue original cluster
{
result_matrix_DFA = as.data.frame(result_matrix_DFA)
orig_sig_elecs <- result_matrix_DFA[result_matrix_DFA$V3 == "***","V1"] #get significant electrodes
orig_sig_elecs <- na.omit(orig_sig_elecs)

# Check adjacency for all pairs in the list
length_sig <- length(orig_sig_elecs)
adjacency_results <- matrix(NA, nrow = length_sig, ncol = length_sig)

for (i in 1:length_sig) {
  for (j in 1:length_sig) {
    if (i != j) {
      adjacency_results[i, j] <- is_adjacent(orig_sig_elecs[i], orig_sig_elecs[j], adj_matrix_df)
    }
  }
}

rownames(adjacency_results) <- orig_sig_elecs
colnames(adjacency_results) <- orig_sig_elecs

#find clusters
{
  # convert the adjacency_results matrix to a graph
  library(igraph)
  g <- graph_from_adjacency_matrix(adjacency_results, mode = "undirected", diag = FALSE)
  # find connected components of the graph
  components <- components(g)
  # extract the clusters
  clusters <- split(names(V(g)), components$membership)
}

#calculate sum of T value of biggest cluster
{cluster_sizes <- sapply(clusters, length)
biggest_cluster_index <- which.max(cluster_sizes) # find the index of the biggest cluster
biggest_cluster_DFA <- clusters[[biggest_cluster_index]] # select the biggest cluster
smallest_cluster_index <- which.min(cluster_sizes) # find the index of the smallest cluster
smallest_cluster_DFA <- clusters[[smallest_cluster_index]] # select the smallest cluster

x <- result_matrix_DFA[result_matrix_DFA$V1 %in% biggest_cluster_DFA, ]
sum_orig_Tval_bigClust = sum(as.numeric(x$V6))
x <- result_matrix_DFA[result_matrix_DFA$V1 %in% smallest_cluster_DFA, ]
sum_orig_Tval_smallClust = sum(as.numeric(x$V6))
}

#plot
{
  par(cex = 1.5)  # Adjust the size multiplier as needed
  
  hist(sum_perm_tvalue, xlab = "sum T-value", breaks = 150)
  abline(v=sum_orig_Tval_bigClust, col="red", lwd=2)
  text(x = sum_orig_Tval_bigClust, y = 782, labels = "big cluster", pos = 4, col = "red", cex = 1.2)
  abline(v=sum_orig_Tval_smallClust, col="green4", lwd=2)
  text(x = sum_orig_Tval_smallClust, y = 782, labels = "small cluster", pos = 2, col = "green4", cex = 1.2)
  
  percentile <- quantile(sum_perm_tvalue, 0.95)
  abline(v=percentile, col="blue", lwd=2)
  #text(x = percentile, y = 720, labels = "95th percentile", pos = 2, col = "blue", cex = 1.2)
  
}

{
hist(sum_perm_tvalue, xlab = "sum T-value", breaks = 150)
abline(v=sum_orig_Tval_bigClust, col="red", lwd=2)
text(x = sum_orig_Tval_bigClust, y = 760, labels = "big cluster", pos = 4, col = "red")
abline(v=sum_orig_Tval_smallClust, col="green4", lwd=2)
text(x = sum_orig_Tval_smallClust, y = 760, labels = "small cluster", pos = 2, col = "green4")

percentile <- quantile(sum_perm_tvalue, 0.95)
abline(v=percentile, col="blue", lwd=2)
text(x = percentile, y = 710, labels = "95th percentile", pos = 4, col = "blue")

}
}
# DFA big cluster is significant at 5%
# DFA small cluster not significant
}
  
#model output for model with average DFA
{data$DFA_bigClust <- rowMeans(data[, biggest_cluster_DFA])
data$DFA_smallClust = rowMeans(data[,smallest_cluster_DFA])

model_DFA_bigClust <- lmer(accuracy ~ condition * DFA + (1|ID), data = data)
Anova(model_DFA_bigClust, type = 3)
summary(model_DFA_bigClust)

model_DFA_bigClust_ctrl <- lmer(accuracy ~ condition * DFA + Gender + zWM_tired + zWM_conc + zIQ + zAge + zDFS + zRelAlphaPow_rest + zBIS + zBAS + zBMI + zDSBack + (1|ID), data = data)
Anova(model_DFA_bigClust_ctrl, type = 3)#results hold
summary(model_DFA_bigClust_ctrl)
} 
  
#write table to word
{
  library(officer)
  library(flextable)
  library(broom)
  
  anova_results = Anova(model_DFA_bigClust_ctrl, type = 3)
  anova_df <- broom::tidy(anova_results)
  anova_df <- anova_df %>%
    mutate(across(where(is.numeric), ~ round(., 3)))
  ft <- flextable(anova_df)
  doc <- read_docx()
  doc <- body_add_flextable(doc, value = ft)
  print(doc, target = "~/Work/PhD/WORMCRI/Results/anova_results_DFA_ctrl.docx")
  
}

#plot 
{ 
  keep = c("ID","correct","condition","DFA")
  plotdata = data[ , (names(data) %in% keep)]
  plotdata = na.omit(plotdata)
  plotdata$fitted = fitted(model_DFA_bigClust)
  plotdata = unique(plotdata)
  
  ggplot(data = plotdata, aes(DFA, fitted, colour = condition, shape = condition, fill = condition)) +
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
    labs(x = "average DFA (big cluster)") +
    theme(axis.title.y = element_text(margin = margin(r = 15), size = 14, family = "Arial")) +
    theme(axis.title.x = element_text(margin = margin(r = 15), size = 14, family = "Arial")) +
    theme(strip.background = element_blank(),  # Remove the background of facet labels
          strip.text = element_text(size = 12)) +
    scale_x_continuous(breaks = c(0.65, 0.95))
}

#posthoc
{
subset_ign = data[data$condition == "ignore",]
model_DFA_bigClust_ign <- lm(accuracy ~ DFA, data = subset_ign)
summary(model_DFA_bigClust_ign)
#             Estimate   Std. Error  t value  Pr(>|t|)    
# (Intercept)   1.05944    0.07111   14.899   <2e-16 ***
# DFA_bigClust -0.21258    0.09118   -2.331   0.0226 * 


subset_upt = data[data$condition == "update",]
model_DFA_bigClust_upt <- lm(accuracy ~ DFA, data = subset_upt)
summary(model_DFA_bigClust_upt)
#             Estimate   Std. Error  t value  Pr(>|t|)    
# (Intercept)   0.93003    0.07291   12.756   <2e-16 ***
# DFA_bigClust  0.01894    0.09349   0.203     0.84 


subset_m1 = data[data$condition == "control_long",]
model_DFA_bigClust_m1 <- lm(accuracy ~ DFA, data = subset_m1)
summary(model_DFA_bigClust_m1)
#              Estimate   Std. Error  t value  Pr(>|t|)    
# (Intercept)   1.10408    0.05780    19.101   < 2e-16 ***
# DFA_bigClust -0.22140    0.07412    -2.987   0.00386 ** 


subset_m2 = data[data$condition == "control_short",]
model_DFA_bigClust_m2 <- lm(accuracy ~ DFA, data = subset_m2)
summary(model_DFA_bigClust_m2)
#               Estimate   Std. Error  t value  Pr(>|t|)    
# (Intercept)    1.04081    0.03908    26.631   <2e-16 ***
# DFA_bigClust  -0.11078    0.05011    -2.211   0.0303 * 
}
}#end DFA section


#stats oneF-behavior ----
{
  
  #loop over electrodes
  {# Find the index of 'Fp1' in the column names
  start_index <- which(colnames(data) == 'oneF_Fp1')
  # Find the index of 'P08' in the column names
  end_index <- which(colnames(data) == 'oneF_PO8')
  
  eleclist <- colnames(data)[start_index:end_index]
  result_matrix_oneF <- matrix(nrow = length(eleclist), ncol = 7)
  count = 1
  
  for (elec in eleclist) {
    # Build the formula with the current ELEC variable
    formula1 <- as.formula(paste("accuracy ~ condition *", elec, "+ (1|ID)"))
    formula2 <- as.formula(paste("accuracy ~ condition *", elec, "* AAratio + (1|ID)"))
    print(paste("Processing elec:", elec))
    # Fit the glmer model
    model1 <- lmer(formula1, data = data)
    model2 <- lmer(formula2, data = data)
    
    model_summary1 <- Anova(model1, type = 3)
    model_summary2 <- Anova(model2, type = 3)
    tvals1 = model_summary1[4,1]
    tvals2 = model_summary2[8,1]
    p_values1 <- model_summary1$Pr
    p_DFA <- p_values1[4]
    p_values2 <- model_summary2$Pr
    p_DFA_AA <- p_values2[8]
    
    result_matrix_oneF[count,1] = elec
    result_matrix_oneF[count,2] = p_DFA
    result_matrix_oneF[count,3] = tvals1
    if (p_DFA < 0.05){
      result_matrix_oneF[count,4] = '***'}
    result_matrix_oneF[count,5] = p_DFA_AA
    result_matrix_oneF[count,6] = tvals2 #store tvalue for the output of model1 condition*DFA effect (for permutation testing later)
    if (p_DFA_AA < 0.05){
      result_matrix_oneF[count,7] = '***'}
    
    count = count + 1
  }
  write.xlsx(as.data.frame(result_matrix_oneF), file = "~/Work/PhD/WORMCRI/Results/result_matrix_oneF.xlsx")
}
  
  #perform permutation testing on effect
  {#load and format adjacency matrix
  {
    adj_matrix <- read_excel("~/Work/PhD/WORMCRI/scripts/adjacency_matrix.xlsx")
    adj_matrix_df <- as.data.frame(adj_matrix)
    
    # Change column and rownames names
    colnames(adj_matrix_df) <- paste("oneF", colnames(adj_matrix), sep = "_")
    rownames(adj_matrix_df) <- adj_matrix_df[,1] # Adjust this if your row names are in a different column
    adj_matrix_df <- adj_matrix_df[,-1] # Remove the column now serving as row names, if necessary
    rownames(adj_matrix_df) <- paste("oneF", rownames(adj_matrix_df), sep = "_")
  }
  
  #perm loop
  {
    N_perm = 1000
    sum_perm_tvalue_oneF <- matrix(nrow = N_perm, ncol = 1)
    
    for (perm in 1:N_perm){
      print(paste("perm # ", perm))
      #shuffle data
      permdata <- data %>%
        mutate(accuracy = sample(accuracy, size = n(), replace = FALSE))
      
      result_matrix_perm <- matrix(nrow = length(eleclist), ncol = 5)
      count = 1
      
      for (elec in eleclist) {
        # Build the formula with the current ELEC variable
        formula <- as.formula(paste("accuracy ~ condition * AAratio *", elec, "+ (1|ID)"))
        # Fit the glmer model
        model <- lmer(formula, data = permdata)
        model_summary <- Anova(model, type = 3)
        p_values <- model_summary$Pr
        p_DFA <- p_values[8]
        tvals = summary(model)$coefficients[,"t value"]
        tval_upt = tvals[16]
        
        result_matrix_perm[count,1] = elec
        result_matrix_perm[count,2] = p_DFA
        if (p_DFA < 0.05){
          result_matrix_perm[count,3] = '***'}
        result_matrix_perm[count,4] = tval_upt
        count = count + 1
      }
      #collect sum T value for sig. cluster of adjacent electrodes
      result_matrix_perm = as.data.frame(result_matrix_perm)
      sig_elecs <- result_matrix_perm[result_matrix_perm$V3 == "***","V1"] #get significant electrodes
      sig_elecs <- na.omit(sig_elecs)
      
      if (length(sig_elecs) == 0) {
        sum_perm_tvalue[perm,] = 0
      } else {
        
        # Check adjacency for all pairs in the list
        length_sig <- length(sig_elecs)
        adjacency_results <- matrix(NA, nrow = length_sig, ncol = length_sig)
        
        for (i in 1:length_sig) {
          for (j in 1:length_sig) {
            if (i != j) {
              adjacency_results[i, j] <- is_adjacent(sig_elecs[i], sig_elecs[j], adj_matrix_df)
            }
          }
        }
        
        rownames(adjacency_results) <- sig_elecs
        colnames(adjacency_results) <- sig_elecs
        
        #find clusters
        {
          # convert the adjacency_results matrix to a graph
          library(igraph)
          g <- graph_from_adjacency_matrix(adjacency_results, mode = "undirected", diag = FALSE)
          # find connected components of the graph
          components <- components(g)
          # extract the clusters
          clusters <- split(names(V(g)), components$membership)
        }
        
        #calculate sum of T value of biggest cluster
        {
          cluster_sizes <- sapply(clusters, length)
          biggest_cluster_index <- which.min(cluster_sizes) # find the index of the biggest cluster
          biggest_cluster <- clusters[[biggest_cluster_index]] # select the biggest cluster
          
          x <- result_matrix_perm[result_matrix_perm$V1 %in% biggest_cluster, ]
          sum_current_Tval = sum(as.numeric(x$V4))
        }
        sum_perm_tvalue_oneF[perm,] = sum_current_Tval
      }
      
      
    }  
  
  # Save the current workspace
  save.image(file = "WS_perm_1000_oneF.RData")}
  
  
  #compare to summed Tvalue original cluster
  {
    result_matrix_oneF = as.data.frame(result_matrix_oneF_AA)
    orig_sig_elecs <- result_matrix_oneF[result_matrix_oneF$V5 == "***","V1"] #get significant electrodes
    orig_sig_elecs <- na.omit(orig_sig_elecs)
    
    # Check adjacency for all pairs in the list
    length_sig <- length(orig_sig_elecs)
    adjacency_results <- matrix(NA, nrow = length_sig, ncol = length_sig)
    
    for (i in 1:length_sig) {
      for (j in 1:length_sig) {
        if (i != j) {
          adjacency_results[i, j] <- is_adjacent(orig_sig_elecs[i], orig_sig_elecs[j], adj_matrix_df)
        }
      }
    }
    
    rownames(adjacency_results) <- orig_sig_elecs
    colnames(adjacency_results) <- orig_sig_elecs
    
    #find clusters
    {
      # convert the adjacency_results matrix to a graph
      library(igraph)
      g <- graph_from_adjacency_matrix(adjacency_results, mode = "undirected", diag = FALSE)
      # find connected components of the graph
      components <- components(g)
      # extract the clusters
      clusters <- split(names(V(g)), components$membership)
    }
    
    #calculate sum of T value of biggest cluster
    {cluster_sizes <- sapply(clusters, length)
      biggest_cluster_index <- which.max(cluster_sizes) # find the index of the biggest cluster
      biggest_cluster <- clusters[[biggest_cluster_index]] # select the biggest cluster
      smallest_cluster_index <- which.min(cluster_sizes) # find the index of the biggest cluster
      smallest_cluster <- clusters[[smallest_cluster_index]] # select the biggest cluster
      
      x <- result_matrix_oneF[result_matrix_oneF$V1 %in% biggest_cluster, ]
      sum_orig_Tval_bigClust = sum(as.numeric(x$V6))
      x <- result_matrix_oneF[result_matrix_oneF$V1 %in% smallest_cluster, ]
      sum_orig_Tval_smallClust = sum(as.numeric(x$V6))
    }
    #library(openxlsx)
    #write.xlsx(result_matrix_oneF, file = "result_matrix_oneF.xlsx")
    
    #plot
    {
      par(cex = 1.5)  # Adjust the size multiplier as needed
      
      hist(sum_perm_tvalue_oneF, xlab = "sum T-value", breaks = 150)
      abline(v=sum_orig_Tval_bigClust, col="red", lwd=2)
      text(x = sum_orig_Tval_bigClust, y = 118, labels = "big cluster", pos = 4, col = "red", cex = 1.2)
      abline(v=sum_orig_Tval_smallClust, col="green4", lwd=2)
      text(x = sum_orig_Tval_smallClust, y = 118, labels = "small cluster", pos = 2, col = "green4", cex = 1.2)
      
      percentile <- quantile(sum_perm_tvalue, 0.95)
      abline(v=percentile, col="blue", lwd=2)
      text(x = percentile, y = 108, labels = "95th percentile", pos = 2, col = "blue", cex = 1.2)
      
    }
    #oneF cluster is significant at 99%
  }
}
  #model output for average oneF 
  {
  data$oneF_smallClust <- rowMeans(data[,smallest_cluster])
  data$oneF_bigClust <- rowMeans(data[,biggest_cluster])
  
  model_oneF_smallClust <- lmer(accuracy ~ condition * oneF_smallClust * AAratio + (1|ID), data = data)
  Anova(model_oneF_smallClust, type = 3)
  #plot_model(model_oneF_smallClust, type = "pred", terms = c("oneF_smallClust","AAratio","condition"))+ theme_classic()
  
  model_oneF_bigClust <- lmer(accuracy ~ condition * oneF_bigClust * AAratio + (1|ID), data = data)
  Anova(model_oneF_bigClust, type = 3)
  plot_model(model_oneF_bigClust, type = "pred", terms = c("oneF","AAratio","condition"))+ theme_classic()
  
  model_oneF_AA <- lmer(accuracy ~ condition * oneF * AAratio + (1|ID), data = data)
  Anova(model_oneF_AA, type = 3)
  plot_model(model_oneF_AA, type = "pred", terms = c("oneF","condition"))+ theme_classic()

  
  #plot 
  { 
    keep = c("ID","accuracy","condition","oneF", "AA_group")
    plotdata = data[ , (names(data) %in% keep)]
    plotdata = na.omit(plotdata)
    plotdata$fitted = fitted(model_oneF_bigClust)
    plotdata = unique(plotdata)
    
    ggplot(data = plotdata, aes(oneF, fitted, colour = AA_group, shape = AA_group, fill = AA_group)) +
      geom_smooth(method = "lm") +
      geom_point(alpha = .4, size = 2) +
      scale_color_manual(values = c("#fa6640ff","#3ea68fff")) +
      scale_fill_manual(values = c("#fa6640ff","#3ea68fff")) +
      scale_shape_manual(values = c(25, 15)) +  # Set filled triangle (24) and filled square (22) shapes for groups
      theme_classic() +
      facet_grid(.~ condition) +
      theme(axis.text.x = element_text(size = 14, color = "black")) +
      theme(axis.text.y = element_text(size = 14, color = "black")) +
      labs(y = "predicted accuracy") +
      labs(x = "average 1/f slope") +
      theme(axis.title.y = element_text(margin = margin(r = 15), size = 14)) +
      theme(axis.title.x = element_text(margin = margin(r = 15), size = 14)) +
      theme(strip.background = element_blank(),  # Remove the background of facet labels
            strip.text = element_text(size = 12)) 
  }
  
  #check mean AAratio per group
  {
    mean_AAratio_by_group <- data %>%
      group_by(AA_group) %>%
      summarize(
        mean = if(all(is.na(AAratio))) NA else mean(AAratio, na.rm = TRUE),
        sd = if(all(is.na(AAratio))) NA else sd(AAratio, na.rm = TRUE),
        min = if(all(is.na(AAratio))) NA else min(AAratio, na.rm = TRUE),
        max = if(all(is.na(AAratio))) NA else max(AAratio, na.rm = TRUE)
      )
  }
  
  #posthoc
  {subset_ign = data[data$condition == "ignore",]
    model_oneF_bigClust_ign <- lm(accuracy ~ oneF*AAratio, data = subset_ign)
    summary(model_oneF_bigClust_ign)
    #                         Estimate   Std. Error  t value  Pr(>|t|)    
    # oneF_bigClust:AAratio   -3.5156     1.4934     -2.354   0.0217 *
    
    
    subset_upt = data[data$condition == "update",]
    model_oneF_bigClust_upt <- lm(accuracy ~ oneF*AAratio, data = subset_upt)
    summary(model_oneF_bigClust_upt)
    #                         Estimate   Std. Error  t value  Pr(>|t|)    
    # oneF_bigClust:AAratio   1.3016     1.5573      0.836    0.406 
    
    
    subset_m1 = data[data$condition == "control_long",]
    model_oneF_bigClust_m1 <- lm(accuracy ~ oneF*AAratio, data = subset_m1)
    summary(model_oneF_bigClust_m1)
    #                        Estimate   Std. Error  t value  Pr(>|t|)    
    # oneF_bigClust:AAratio  -1.7344    1.2763      -1.359    0.179
    
    
    subset_m2 = data[data$condition == "control_short",]
    model_oneF_bigClust_m2 <- lm(accuracy ~ oneF*AAratio, data = subset_m2)
    summary(model_oneF_bigClust_m2)
    #                        Estimate   Std. Error  t value  Pr(>|t|)    
    # oneF_bigClust:AAratio  -1.9424    0.8149      -2.384   0.0201 *
  }
  
  
  #check summed over both clusters with control variables
  {
  model_oneF_AA_ctrl <- lmer(accuracy ~ condition * oneF * AAratio  + Gender + zWM_tired + zWM_conc + zIQ + zAge + zDFS + zRelAlphaPow_rest + zBIS + zBAS + zBMI + zDSBack + (1|ID), data = data)
  Anova(model_oneF_AA_ctrl, type = 3)
  #plot_model(model_oneF_AA_ctrl, type = "pred", terms = c("oneF","AAratio", "condition"))
  
  #check co-linerarity of onef and rel alpha power
  {model_oneF_AA_relAlpha <- lmer(accuracy ~ condition + oneF + AAratio + relAlphaPow_rest + (1|ID), data = data)
    vif(model_oneF_AA_relAlpha)}
  }#results hold
  
  #write table to word
  {
    library(officer)
    library(flextable)
    library(broom)
    
    anova_results = Anova(model_oneF_AA_ctrl, type = 3)
    anova_df <- broom::tidy(anova_results)
    anova_df <- anova_df %>%
      mutate(across(where(is.numeric), ~ round(., 3)))
    ft <- flextable(anova_df)
    doc <- read_docx()
    doc <- body_add_flextable(doc, value = ft)
    print(doc, target = "~/Work/PhD/WORMCRI/Results/anova_results_oneF_ctrl.docx")
    
    anova_results = Anova(model_oneF_smallClust_ctrl, type = 3)
    anova_df <- broom::tidy(anova_results)
    ft <- flextable(anova_df)
    doc <- read_docx()
    doc <- body_add_flextable(doc, value = ft)
    print(doc, target = "~/Work/PhD/WORMCRI/Results/anova_results_oneF_ctrl_small.docx")
    
  }
  }
  }#end oneF


#check quadratic effects of AAratio ----

data$AAratio2 = data$AAratio*data$AAratio

#loop over oneF
{# Find the index of 'Fp1' in the column names
  start_index <- which(colnames(data) == 'oneF_Fp1')
  # Find the index of 'P08' in the column names
  end_index <- which(colnames(data) == 'oneF_PO8')
  
  eleclist <- colnames(data)[start_index:end_index]
  result_matrix_oneF <- matrix(nrow = length(eleclist), ncol = 7)
  count = 1
  
  for (elec in eleclist) {
    # Build the formula with the current ELEC variable
    formula1 <- as.formula(paste("accuracy ~ condition *", elec, "+ (1|ID)"))
    formula2 <- as.formula(paste("accuracy ~ condition *", elec, "* AAratio2 + (1|ID)"))
    print(paste("Processing elec:", elec))
    # Fit the glmer model
    model1 <- lmer(formula1, data = data)
    model2 <- lmer(formula2, data = data)
    
    model_summary1 <- Anova(model1, type = 3)
    model_summary2 <- Anova(model2, type = 3)
    tvals1 = summary(model1)$coefficients[8,4]
    tvals2 = summary(model2)$coefficients[16,4]
    p_DFA <- model_summary1$Pr[4]
    p_DFA_AA <-  model_summary2$Pr[8]
    
    result_matrix_oneF[count,1] = elec
    result_matrix_oneF[count,2] = p_DFA
    result_matrix_oneF[count,3] = tvals1
    if (p_DFA < 0.05){
      result_matrix_oneF[count,4] = '***'}
    result_matrix_oneF[count,5] = p_DFA_AA
    result_matrix_oneF[count,6] = tvals2 #store tvalue for the output of model1 condition*DFA effect (for permutation testing later)
    if (p_DFA_AA < 0.05){
      result_matrix_oneF[count,7] = '***'}
    
    count = count + 1
  }
  write.xlsx(as.data.frame(result_matrix_oneF), file = "~/Work/PhD/WORMCRI/Results/result_matrix_oneF_AAquad.xlsx")
  
}
#some more electrodes become significant

#loop DFA
{
  # Find the index of 'Fp1' in the column names
  start_index <- which(colnames(data) == 'DFA_Fp1')
  # Find the index of 'P08' in the column names
  end_index <- which(colnames(data) == 'DFA_PO8')
  
  eleclist <- colnames(data)[start_index:end_index]
  result_matrix_DFA <- matrix(nrow = length(eleclist), ncol = 7)
  count = 1
  
  for (elec in eleclist) {
    formula2 <- as.formula(paste("accuracy ~ condition *", elec, "* AAratio2 + (1|ID)"))
    print(paste("Processing elec:", elec))
    model2 <- lmer(formula2, data = data)
    
    model_summary2 <- Anova(model2, type = 3)
    tvals2 = summary(model2)$coefficients[16,4]
    p_DFA_AA <- model_summary2$Pr[8]
 
    result_matrix_DFA[count,1] = elec
    result_matrix_DFA[count,2] = p_DFA_AA
    result_matrix_DFA[count,3] = tvals2
    if (p_DFA_AA < 0.05){
      result_matrix_DFA[count,4] = '***'}
    count = count + 1
  }
  #AAratio as quadratic does not change anything
}
#for DFA it stays the same...

#AAratio as quadratic does not change anything

#check which model fits better
{model_oneF_AA2 = lmer(accuracy ~ condition * oneF * AAratio2 + (1|ID), data = data)
model_oneF_AA = lmer(accuracy ~ condition * oneF * AAratio + (1|ID), data = data)
BIC(model_oneF_AA, model_oneF_AA2)
anova(model_oneF_AA, model_oneF_AA2)}
#model with quadratic term fits better, but negligibly






#oneF x DFA ----

  #read in results matrix to extract big cluster DFA electrodes
  {
    result_matrix_DFA <- read_excel("~/Work/PhD/WORMCRI/Analyses/Results/result_matrix_DFA.xlsx")
    filtered_electrodes <- result_matrix_DFA[result_matrix_DFA$sig == "***" & result_matrix_DFA$cluster == 2,"elec"]
    eleclist_sig <- as.list(filtered_electrodes)
    eleclist_sig <- unlist(eleclist_sig)
    eleclist_sig_DFA = eleclist_sig[!is.na(eleclist_sig)]
    #data$DFA = rowMeans(data[, eleclist_sig_DFA])
  }
  
  #create average oneF for sig. DFA electrodes
  {
  eleclist_oneF_DFA <- sub("DFA","oneF",eleclist_sig_DFA)
  data$oneF_DFA <- rowMeans(data[, eleclist_oneF_DFA])
  median_oneF = median(data$oneF_DFA, na.rm = T)
  data$oneF_group <- ifelse(data$oneF_DFA <= median_oneF, "flatter", "steeper")
  subset = data[data$condition == "ignore" | data$condition == "update",]
}
  
  model <- lmer(accuracy ~ DFA*oneF_DFA*condition + (1|ID), data = subset)
  summary(model)
  Anova(model, type = 3)
  plot_model(model, type = "pred", terms = c("DFA","oneF_DFA"), show.data = T)+ theme_classic()
  
  # DFA_bigClust:oneF_DFA           5.7406  1    0.01658 *
  # condition:oneF_DFA              0.5203  1    0.47072  
  # DFA_bigClust:condition:oneF_DFA 0.6905  1    0.40599 
  
  #check VIF because DFA and oneF might be co-linear
  {model_w0_inter <- lmer(accuracy ~ DFA+oneF_DFA+condition + (1|ID), data = subset)
  vif(model_w0_inter)
  #      DFA  oneF_DFA condition 
  # 1.108649  1.108649  1.000000
  #none above 5 ... no troublesome co-linearity
  }
  
  #plot 
  { 
    keep = c("ID","accuracy","condition","DFA", "oneF_group")
    plotdata = subset[ , (names(subset) %in% keep)]
    plotdata = na.omit(plotdata)
    plotdata$fitted = fitted(model)
    plotdata = unique(plotdata)
    
    ggplot(data = plotdata, aes(DFA, fitted, colour = condition, shape = condition, fill = condition)) +
      geom_smooth(method = "lm") +
      geom_point(alpha = .4, size = 3) +
      scale_color_manual(values = c("#ffc107ff","#8aaedcff")) +
      scale_fill_manual(values = c("#ffc107ff","#8aaedcff")) +
      scale_shape_manual(values = c(15, 17)) +  # Set filled triangle (24) and filled square (22) shapes for groups
      theme_classic() +
      facet_grid(.~ oneF_group)+
      theme(axis.text.x = element_text(size = 18, color = "black", family = "arial")) +
      theme(axis.text.y = element_text(size = 18, color = "black", family = "arial")) +
      labs(y = "predicted accuracy") +
      labs(x = "average DFA exponent (big cluster)") +
      theme(axis.title.y = element_text(margin = margin(r = 15), size = 18, family = "Arial")) +
      theme(axis.title.x = element_text(margin = margin(r = 15), size = 18, family = "Arial")) +
      theme(strip.background = element_blank(),  # Remove the background of facet labels
            strip.text = element_text(size = 18)) 
  }
  
  #posthoc
  {
  subset_excit = subset[subset$oneF_group == "flatter",]
  model_excit <- lmer(accuracy ~ DFA*condition + (1|ID), data = subset_excit)
  summary(model_excit)
  Anova(model_excit, type = 3)
  #                             Estimate  Std. Error  t value   p
  #(Intercept)                    0.9771     0.1217   8.030
  #DFA_bigClust                  -0.1057     0.1641  -0.644     0.519
  #conditionupdate               -0.1710     0.1176  -1.454
  #DFA_bigClust:conditionupdate   0.2868     0.1586   1.809     0.07052 .
  
  #plot_model(model_excit, type = "pred", terms = c("DFA_bigClust"))+ theme_classic() + ggtitle("excit")
  
  subset_excit_upt = subset_excit[subset_excit$condition == "update",]
  model_excit_upt <- lm(accuracy ~ DFA, data = subset_excit_upt)
  summary(model_excit_upt)
  
  subset_excit_ign = subset_excit[subset_excit$condition == "ignore",]
  model_excit_ign <- lm(accuracy ~ DFA, data = subset_excit_ign)
  summary(model_excit_ign)
  
  subset_inhib = subset[subset$oneF_group == "steeper",]
  model_inhib <- lmer(accuracy ~ DFA*condition + (1|ID), data = subset_inhib)
  summary(model_inhib)
  Anova(model_inhib, type = 3)
  #                             Estimate  Std. Error  t value   p
  #(Intercept)                   1.15720    0.09861   11.735
  #DFA_bigClust                 -0.32847    0.12062   -2.723    0.006***
  #conditionupdate              -0.08709    0.11829   -0.736
  #DFA_bigClust:conditionupdate  0.18045    0.14470    1.247
  
  plot_model(model_inhib, type = "pred", terms = c("DFA_bigClust"))+ theme_classic() + ggtitle("inhib")
  }
  
  #how do oneF and DFA relate to each other?
  {subset_wide = subset[subset$condition == "ignore",] #get dataframe into wide format (one row per participant only)
  model_oneF_DFA = lm(oneF_DFA ~ DFA, data = subset_wide)
  summary(model_oneF_DFA)
  #Coefficients:
  #               Estimate Std. Error t value   Pr(>|t|)    
  # (Intercept)   1.1342     0.2384   4.759     9.96e-06 ***
  # DFA           0.8489     0.3056   2.777     0.007 **
  plot_model(model_oneF_DFA, type = "pred", terms = c("DFA"), show.data = T) + 
    theme_classic() +
    labs(y = "PSD slope") +
    annotate("text", x = 0.92, y = 1.25, label = "p = 0.007", hjust = 0, vjust = 0.2, size = 4, color = "black",fontface = "italic")
  
  
  #effect of mean DFA between the oneF groups?
  t.test(subset_wide$DFA ~ subset_wide$oneF_group, subset_wide)
  #t = -3.5784, df = 69.307, p-value = 0.0006362
  #plot boxplot per group
  {ggplot(subset_wide, aes(x = oneF_group, y = DFA, fill = oneF_group)) +
    geom_violin(trim = TRUE, alpha = 0.8) +  # Draw violin plots 
    geom_boxplot(width = 0.2, fill = "white", color = "black") + 
    geom_point(alpha = 0.8, size = 2, position = position_jitter(width = 0.1)) +  # Add jittered points
    labs(x = "PSD slope group", y = "DFA exponent") +
    theme_minimal() +
    scale_color_manual(values = c("gray", "gray25")) +  # Manual color scale for points and outlines
    scale_fill_manual(values = c("gray", "gray25")) +  # Manual fill color scale for violins and boxplots
    scale_shape_manual(values = c(15, 18)) +  # Assign shapes to conditions
    theme(
      axis.text.x = element_text(size = 18, colour = "black"),  # Adjust x-axis tick label size
      axis.text.y = element_text(size = 18, colour = "black"),  # Adjust y-axis tick label size
      axis.title.x = element_text(size = 18),  # Adjust x-axis label size
      axis.title.y = element_text(size = 18)   # Adjust y-axis label size 
    )}
  
}

#CONCLUSIONS: 
# stats on average accuracy show condition*DFA effect in P6 P4, FC4, C4, P2, P6 etc.
# no effect of condition*DFA*AAratio in any electrode
# when adding AAratio to condition*DFA, the two-interaction of condition*DFA is gone
# so AAratio does explain some variance there?

# effect for 1/f slope*condition*AAratio in many electrodes 
# no 1/f * condition effect

# DFA and oneF interact on accuracy, not different per condition

#additional checks: ----

#does AAratio directly relate to oneF?
{model = lm(AAratio ~ oneF, data = subset_wide)
summary(model)
plot_model(model, type = "pred", show.data = T) + theme_classic()
#               Estimate  Std. Error  t value   Pr(>|t|)    
# (Intercept)  0.234802     0.020800  11.289   <2e-16 ***
# oneF         0.007875     0.011662  0.675     0.502 
}
#no it does not

#does AAratio directly relate to DFA?
{model = lm(AAratio ~ DFA, data = subset_wide)
summary(model)
#              Estimate   Std. Error t value  Pr(>|t|)    
# (Intercept)  0.23499    0.02255    10.422   1.4e-15 ***
# DFA          0.01769    0.02882    0.614    0.541  
plot_model(model, type = "pred", show.data = T) + theme_classic()
}#no it does not

#what are our AAratio ranges?
{AA_low = data[data$AA_group == "low",]
min(AA_low$AAratio, na.rm = T)#0.2048945
mean(AA_low$AAratio, na.rm = T)#0.22955
max(AA_low$AAratio, na.rm = T)#0.2462801

AA_high = data[data$AA_group == "high",]
min(AA_high$AAratio, na.rm = T)#0.247
mean(AA_high$AAratio, na.rm = T)#0.268
max(AA_high$AAratio, na.rm = T)#0.299
}
  
#does the difference in variance in conditions explain any effects?
{
    # Calculate variance in accuracy per condition
    variance_per_condition <- data %>%
      group_by(condition) %>%
      summarize(variance = sd(accuracy))
    
    # condition        variance
    # <fct>            <dbl>
    # 1 ignore         0.00596
    # 2 control_long   0.00411
    # 3 control_short  0.00179
    # 4 update         0.00582
    
    # Perform Levene's Test
    levene_test <- leveneTest(accuracy ~ condition, data = data)
    #there is an overall difference. Posthoc, where exactly:
    {
    subset_ign_upt = data[data$condition == "ignore" | data$condition == "update",]
    leveneTest(accuracy ~ condition, data = subset_ign_upt)
    #       Df F value Pr(>F)
    #group   1  2.4401 0.1205
    
    subset_ign_m1 = data[data$condition == "ignore" | data$condition == "control_long",]
    leveneTest(accuracy ~ condition, data = subset_ign_m1)
    #       Df F value Pr(>F)
    #group   1   2.579 0.1105
    
    subset_ign_m2 = data[data$condition == "ignore" | data$condition == "control_short",]
    leveneTest(accuracy ~ condition, data = subset_ign_m2)
    #       Df F value Pr(>F)
    #group   1  17.721 4.482e-05 ***
    }
    #the variance in control short is significantly different from the other conditions
    
}
  
#does alpha power alone explain condition-specific task performance
{
model_alpha_pow <- lmer(accuracy ~ condition * relAlphaPow_rest *AAratio + (1|ID), data = data)
Anova(model_alpha_pow, type = 3)
summary(model_alpha_pow)

#plot 
{ 
  keep = c("ID","correct","condition","relAlphaPow_rest")
  plotdata = data[ , (names(data) %in% keep)]
  plotdata = na.omit(plotdata)
  plotdata$fitted = fitted(model_alpha_pow)
  plotdata = unique(plotdata)
  plotdata$condition <- factor(plotdata$condition, levels = c("ignore", "control_long", "update", "control_short"))

  
  ggplot(data = plotdata, aes(relAlphaPow_rest, fitted, colour = condition, shape = condition, fill = condition)) +
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
    labs(x = "relAlphaPow_rest") +
    theme(axis.title.y = element_text(margin = margin(r = 15), size = 14, family = "Arial")) +
    theme(axis.title.x = element_text(margin = margin(r = 15), size = 14, family = "Arial")) +
    theme(strip.background = element_blank(),  # Remove the background of facet labels
          strip.text = element_text(size = 12)) 
}

#correlation alpha power and oneF
data_wide = subset[subset$condition == "ignore",]
model_alpha_pow_oneF <- lm(oneF ~ relAlphaPow_rest, data = data_wide)
summary(model_alpha_pow_oneF)
plot_model(model_alpha_pow_oneF, type = "pred", terms = c("relAlphaPow_rest"), show.data = TRUE) + theme_classic()
}


