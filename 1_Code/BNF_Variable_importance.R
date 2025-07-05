library(tidyverse)
library(caret)
library(DALEX)
library(GGally)
library(RColorBrewer)
library(plyr)
library(raster)
library(patchwork)
library(tidyverse)
library(randomForest)
library(ppcor)
min_max_norm <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}
library(iml)
library(gridExtra)
library(viridis)
library(ggplot2)
library(pdp)
library(ggpubr)
library(cowplot)
library(DataExplorer)
library(fastshap) 
library(tidyr)
library(dplyr)
library(ggbeeswarm)
library(readxl)
windowsFonts(Font = windowsFont("Times New Roman"))
par(family = 'Font')
setwd('../0_Data')

##1.SNF################################################################
# for SNF
shui_data <- read_excel("BNFdata_use.xlsx", sheet = "Symbiotic")
shui_data = shui_data[,c(-1)]
select=dplyr::select

# use field with all
Variables<-c('BNF','MAT','MAT_season','MAP','MAP_season','AI','AET','VPD','srad','tmax','tmin','CEC',
             'BD','pH','Sand','Silt','Clay','SWC','SOC','TN','TP','MBC','MBN','FB_ratio','NPP','BNPP',
             'NDVI','Ndep','LDMC','SLA','LNC','LPC','Vcmax',
             'EM_tree','AM_tree','Nfixer')

SNF_data = shui_data[,3:38]
head(SNF_data)
nrow(SNF_data)

# Important factors affecting nitrogen fixation in previous studies
var_imp<-c('MAT','MAP','AET','NPP','Nfixer','TN','SOC')
## VIF test
alias(lm(BNF~.,data = SNF_data))
corvalue=car::vif(lm(BNF~.,data =SNF_data));corvalue
cor_variable=names(corvalue[corvalue <= 4])
length(corvalue)
SNF_data2 = SNF_data%>%select(BNF,!!cor_variable,var_imp)
corvalue=car::vif(lm(BNF~.,data =SNF_data2));corvalue
length(corvalue)
SNF_data2 = SNF_data2[,c("BNF","MAP_season","pH","SWC","TP","NPP","BNPP","Ndep","LNC",
                             "MAT","MAP","AET","Nfixer","TN","SOC")]
corvalue=car::vif(lm(BNF~.,data =SNF_data2));corvalue

vif_values=as.data.frame(corvalue)
vif_values$Variable=rownames(vif_values)
colnames(vif_values)[1]=c("VIF")
vif_values=vif_values%>%arrange(VIF) 

VIF2=vif_values%>%mutate(Variable=factor(Variable,levels =vif_values$Variable ))%>%
  ggplot(aes(x=Variable,y=VIF))+
  geom_bar(stat="identity", width=0.7, fill="#0082C6") + 
  scale_y_continuous(limits = c(0,10))+
  coord_flip()+
  geom_hline(aes(yintercept = 10),linetype=2,linewidth=1)+
  labs(x="",y="VIF") 
VIF2

## factors important analysis
rsq<-vector()
impToPlot<-list()
for (i in 1:100) {
  set.seed(i)
  # use random forests
  mean_ran_for<-randomForest(BNF ~., 
                             data = SNF_data2, ntree=300, importance=TRUE,na.action=na.roughfix, mtry=3)
  
  rsq[i]<-mean_ran_for$rsq[length(mean_ran_for$rsq)]
  impToPlot[[i]]<-unname(randomForest::importance(mean_ran_for)[,1])
}
rf_variables5_full_rsq<-mean(rsq)
rf_variables5_full_std<-sd(rsq)

set.seed(999)
mod.formula4 <- as.formula(paste0(colnames(SNF_data2)[1],"~",paste0(colnames(SNF_data2)[2:ncol(SNF_data2)],collapse = "+")))
rf_variables5 <- caret::train(mod.formula4,data = SNF_data2,
                             trControl=trainControl(method="cv", number=10),
                             method = "ranger",
                             metric = "RMSE")
SNF_rf5 <- DALEX::explain(rf_variables5,label = "SNF",
                        data = SNF_data2,
                        y = SNF_data2$BNF)

per_rf5 <- model_performance(SNF_rf5)
per_rf5
plot(per_rf5)
plot(per_rf5,geom = "boxplot")

importance_rf5<-variable_importance(
  SNF_rf5,
  loss_function = loss_root_mean_square
)
plot(importance_rf5)

# Save factor importance results for Fig.3a
importance_rf5_pivoted <- importance_rf5 %>%
  select(variable, permutation, dropout_loss) %>%
  pivot_wider(
    names_from = variable,        
    values_from = dropout_loss,    
    id_cols = permutation          
  )
# write.xlsx(importance_rf5_pivoted, file = "RelativeImportance.xlsx", sheetName = "SNF", rowNames = TRUE)

#################################
# shapley value analysis
str(SNF_data2)
SNF_data2 <- lapply(SNF_data2, as.numeric)
SNF_data2 <- as.data.frame(SNF_data2)

SNF_rf5 <- DALEX::explain(
  model = rf_variables5,
  label = "SNF",
  data = SNF_data2[, -1], # Exclude the response variable (BNF)
  y = SNF_data2$BNF
)

# Compute Shapley values using fastshap
set.seed(321) 
shap_values <- fastshap::explain(
  object = SNF_rf5$model$finalModel, 
  X = as.matrix(SNF_data2[, -1]), 
  nsim = 100, 
  pred_wrapper = function(model, newdata) {
    predict(model, data = newdata)$predictions
  },
  feature_names = colnames(SNF_data2)[-1]
)

# Convert Shapley values to a data frame
shap_df <- as.data.frame(shap_values)
colnames(shap_df) <- colnames(SNF_data2)[-1] 

# Step 4: Create scatter plots for each predictor vs. its Shapley value
# Combine predictor values with Shapley values
predictor_names <- colnames(SNF_data2)[-1]
shap_var_name <- paste0(predictor_names, "_shap")
colnames(shap_df) <- shap_var_name
plot_data <- cbind(SNF_data2[, -1], shap_df)

# write.xlsx(plot_data, file = "SHAP value.xlsx", sheetName = "SNF", rowNames = TRUE)

# Function to create a scatter plot for a single variable
create_shap_scatter <- function(var_name) {
  shap_var_name <- paste0(var_name, "_shap") # Use Shapley value column
  ggplot(plot_data, aes(x = .data[[var_name]], y = .data[[shap_var_name]])) +
    geom_point(color = "#0082C6", alpha = 0.6) +
    geom_smooth(method = "loess", color = "red", se = T, linewidth = 1) + # Add trend line
    labs(x = var_name, y = paste("Shapley Value for", var_name)) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    )
}

# Generate scatter plots for all predictors
shap_scatter_plots <- lapply(predictor_names, create_shap_scatter)

# Combine plots using patchwork
combined_scatter_plots <- wrap_plots(shap_scatter_plots, ncol = 4) +
  plot_annotation(title = "Shapley Value vs Predictor Relationships", 
                  theme = theme(plot.title = element_text(size = 16, hjust = 0.5)))
print(combined_scatter_plots)
# ggsave("shapley_scatter_SNF.png", combined_scatter_plots, width = 12, height = 10)

## NPP-SNF Fig 3c
var_name <- "NPP"
shap_var_name <- paste0(var_name, "_shap")
ggnpp <- ggplot(plot_data, aes(x = .data[[var_name]]/10, y = .data[[shap_var_name]])) +
  geom_point(color = "white", fill = "grey33", shape = 21, alpha = 1, size = 3) +
  geom_smooth(method = "loess", color = "red", fill = "#FF6666", se = TRUE, size = 2) +
  geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 1) +
  labs(x = paste(var_name, "(gC/m2/yr)"), y = "SHAP value") +
  theme_classic() + ylim(-8,20) +
  theme(
    text = element_text(family = "Times New Roman"),
    axis.text = element_text(size = 32, color = "black"),
    axis.title = element_text(size = 48, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 1),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks = element_line(color = "black")
  )
ggnpp
# ggsave("SNF-NPP.png", ggnpp, width = 6, height = 4)

##2.FNF#################################################################
shui_data <- read_excel("BNFdata_use.xlsx", sheet = "Free-living")
shui_data <- shui_data[,c(-1)]
# use field with all
Variables<-c('BNF','MAT','MAT_season','MAP','MAP_season','AI','AET','VPD','srad','tmax','tmin','CEC',
             'BD','pH','Sand','Silt','Clay','SWC','SOC','TN','TP','MBC','MBN','FB_ratio','NPP','BNPP',
             'NDVI','Ndep','LDMC','SLA','LNC','LPC','Vcmax',
             'EM_tree','AM_tree','Nfixer')
FNF_data <- shui_data[,3:38]
nrow(FNF_data)

var_imp<-c('MAT','MAP','AET','NPP','Nfixer','TN','SOC')
alias(lm(BNF~.,data = FNF_data))
corvalue=car::vif(lm(BNF~.,data =FNF_data));corvalue
cor_variable=names(corvalue[corvalue <= 5])
FNF_data2 = FNF_data %>% select(BNF,!!cor_variable,var_imp)
corvalue=car::vif(lm(BNF~.,data =FNF_data2));corvalue
length(corvalue)
FNF_data2 = FNF_data2[,c("BNF","MAP_season","pH","SWC","TP","NPP","BNPP","Ndep","LNC",
                             "MAT","MAP","AET","Nfixer","TN","SOC")]
corvalue=car::vif(lm(BNF~.,data =FNF_data2));corvalue
vif_values=as.data.frame(corvalue)
vif_values$Variable=rownames(vif_values)
colnames(vif_values)[1]=c("VIF")
vif_values=vif_values%>%arrange(VIF)
VIF2=vif_values%>%mutate(Variable=factor(Variable,levels =vif_values$Variable ))%>%
  ggplot(aes(x=Variable,y=VIF))+
  geom_bar(stat="identity", width=0.7, fill="#0082C6") + 
  scale_y_continuous(limits = c(0,10))+
  coord_flip()+
  geom_hline(aes(yintercept=10),linetype=2,size=1)+
  labs(x="",y="VIF") 
VIF2

## factors important analysis
rsq<-vector()
impToPlot<-list()
for (i in 1:100) {
  set.seed(i)
  # use random forests
  mean_ran_for<-randomForest(BNF ~., 
                             data = FNF_data2, ntree=300, importance=TRUE,na.action=na.roughfix, mtry=3)
  
  rsq[i]<-mean_ran_for$rsq[length(mean_ran_for$rsq)]
  impToPlot[[i]]<-unname(randomForest::importance(mean_ran_for)[,1])
}
rf_variables5_full_rsq<-mean(rsq)
rf_variables5_full_std<-sd(rsq)

set.seed(999)
mod.formula4 <- as.formula(paste0(colnames(FNF_data2)[1],"~",paste0(colnames(FNF_data2)[2:ncol(FNF_data2)],collapse = "+")))
rf_variables5 <- caret::train(mod.formula4,data = FNF_data2,
                             trControl=trainControl(method="cv", number=10),
                             method = "ranger",
                             metric = "RMSE")
FNF_rf5 <- DALEX::explain(rf_variables5,label = "FNF",
                        data = FNF_data2,
                        y = FNF_data2$BNF)
per_rf5 <- model_performance(FNF_rf5)
per_rf5
plot(per_rf5)
plot(per_rf5,geom = "boxplot")

importance_rf5<-variable_importance(
  FNF_rf5,
  loss_function = loss_root_mean_square
)
plot(importance_rf5)

# Save factor importance results for Fig.3b
importance_rf5_pivoted <- importance_rf5 %>%
  select(variable, permutation, dropout_loss) %>%
  pivot_wider(
    names_from = variable,        
    values_from = dropout_loss,   
    id_cols = permutation          
  )
# wb <- loadWorkbook("RelativeImportance.xlsx")  
# addWorksheet(wb, "FNF")      # 添加新 sheet
# writeData(wb, "FNF", importance_rf5_pivoted)
# saveWorkbook(wb, "RelativeImportance.xlsx", overwrite = TRUE)

################################
# shapley value analysis
str(FNF_data2)
FNF_data2 <- lapply(FNF_data2, as.numeric)
FNF_data2 <- as.data.frame(FNF_data2)

FNF_rf5 <- DALEX::explain(
  model = rf_variables5,
  label = "FNF",
  data = FNF_data2[, -1], # Exclude the response variable (BNF)
  y = FNF_data2$BNF
)

# Step 1: Compute Shapley values using fastshap
set.seed(123) # For reproducibility
shap_values <- fastshap::explain(
  object = FNF_rf5$model$finalModel, # Ranger model
  X = as.matrix(FNF_data2[, -1]), # Numeric predictor matrix
  nsim = 100, # Number of Monte Carlo simulations
  pred_wrapper = function(model, newdata) {
    predict(model, data = newdata)$predictions
  },
  feature_names = colnames(FNF_data2)[-1]
)

# Convert Shapley values to a data frame
shap_df <- as.data.frame(shap_values)
colnames(shap_df) <- colnames(FNF_data2)[-1] # Assign predictor names

# Combine predictor values with Shapley values
predictor_names <- colnames(FNF_data2)[-1]
shap_var_name <- paste0(predictor_names, "_shap")
colnames(shap_df) <- shap_var_name
plot_data <- cbind(FNF_data2[, -1], shap_df)

# wb <- loadWorkbook("SHAP value.xlsx")  
# addWorksheet(wb, "FNF")      
# writeData(wb, "FNF", plot_data)
# saveWorkbook(wb, "SHAP value.xlsx", overwrite = TRUE)

# Function to create a scatter plot for a single variable
create_shap_scatter <- function(var_name) {
  shap_var_name <- paste0(var_name, "_shap") # Use Shapley value column
  ggplot(plot_data, aes(x = .data[[var_name]], y = .data[[shap_var_name]])) +
    geom_point(color = "#0082C6", alpha = 0.6) +
    geom_smooth(method = "loess", color = "red", se = T, size = 1) + # Add trend line
    labs(x = var_name, y = paste("Shapley Value for", var_name)) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    )
}
shap_scatter_plots <- lapply(predictor_names, create_shap_scatter)

# Combine plots using patchwork
combined_scatter_plots <- wrap_plots(shap_scatter_plots, ncol = 4) +
  plot_annotation(title = "Shapley Value vs Predictor Relationships", 
                  theme = theme(plot.title = element_text(size = 16, hjust = 0.5)))
print(combined_scatter_plots)
# ggsave("shapley_scatter_FNF.png", combined_scatter_plots, width = 12, height = 10)

## Fig 3d SOC-FNF
var_name <- "SOC"
shap_var_name <- paste0(var_name, "_shap")
ggsoc <- ggplot(plot_data, aes(x = .data[[var_name]]/100, y = .data[[shap_var_name]])) +
  geom_point(color = "white", fill = "grey33", shape = 21, alpha = 1, size = 3) +
  geom_smooth(method = "loess", color = "red", fill = "#FF6666", se = TRUE, size = 2) +
  geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 1) +
  labs(x = paste(var_name, "(% of weight)"), y = "SHAP value") +
  theme_classic() + 
  theme(
    text = element_text(family = "Times New Roman"),
    axis.text = element_text(size = 32, color = "black"),
    axis.title = element_text(size = 48, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 1),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks = element_line(color = "black")
  )
ggsoc
# ggsave("FNF-SOC.png", ggsoc, width = 6, height = 4)