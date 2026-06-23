rm(list=ls())																							#Remove all variables in the current session
library(readr); library(caret); library("pROC"); library(dplyr); library(gridExtra); library(ggplot2); source("http://peterhaschke.com/Code/multiplot.R"); library(ggpubr)
addTaskCallback(function(...) {set.seed(500);TRUE})
ctrl <- trainControl(method = "repeatedcv", sampling='down', number = 10, repeats = 3, savePredictions = TRUE, 
                     classProbs = TRUE, selectionFunction = "best", summaryFunction = twoClassSummary) 
#Prepare train dataset
ppi_A10ccall <- read.csv("./polyA_model.csv");                                                          #Load dataset     
ppi_A10ccall$ppi <- factor(ppi_A10ccall$ppi,levels = c("0", "1"),labels = c("no","yes"))                #Factor PPI variable
ppi_A10ccall <- ppi_A10ccall[, -c(1:2)]                                                                 #Delete the AC and polyX positions
rf_ppi_A10cc <- train(ppi ~ ., data = ppi_A10ccall, method = "rf",	metric = "ROC", trControl = ctrl)   #Train the model
new_files <- c("A-AA.csv", "A-AS.csv")                                                                  #Files to apply the model

#Function to read, clean, predict and tag
predict_on_file <- function(file_path) {
  tag <- tools::file_path_sans_ext(basename(file_path))
  df  <- read.csv(file_path, stringsAsFactors = FALSE)
  df$ppi <- factor(df$ppi, levels = c(0,1), labels = c("no","yes"))
  df_feat <- df[, !names(df) %in% c("ac", names(df)[2], "ppi")]
  probs   <- predict(rf_ppi_A10cc, newdata = df_feat, type = "prob")
  data.frame(file   = tag, length = df$length, score  = probs[ , "yes"], row.names = NULL)
}
#Apply function to both files and bind into one
results <- do.call(rbind, lapply(new_files, predict_on_file))                                          

#Plot results
cmp <- list(c("A-AA","A-AS"))
a<-ggplot(results, aes(x = file, y = score)) +  geom_boxplot() + theme_classic() +
  stat_compare_means( comparisons = cmp, method = "wilcox.test", label = "p.format", label.y = 0.81 ) +
  labs(x = "PolyX", y = "Predicted PPI Score") + ylim(0,0.85)

b<-ggplot(results, aes(x = length, y = score, colour = file)) +  theme_classic() +
  geom_jitter(width = 0.1, height = 0, size = 2, alpha = 0.7) + ylim(0,0.85) +
  labs(x = "Length", y = "Predicted PPI Score", colour = "PolyX") + theme(legend.position = c(0.8, 0.8))

grid.arrange(a,b, layout_matrix = matrix(c(1,2), ncol=2))
pdf('score.pdf', height=3, width=6); grid.arrange(a,b, layout_matrix = matrix(c(1,2), ncol=2));  dev.off()

