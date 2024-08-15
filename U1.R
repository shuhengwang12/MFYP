install.packages("glmnet")
library(glmnet)
library(pROC)
library(sampling)
library(rsample)
library(qiime2R)
library(tidyverse)
library(dplyr)
library(vegan)
library(ggrepel)
library(lme4)
library(tidyr)
library(emmeans)
library (randomForest)
library(rfUtilities)
library (caret)
library(stringr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(randomForestSRC)
getwd()
#grouped by community
set.seed(123)
removezerodata_meta_NA$community.number <- as.factor(removezerodata_meta_NA$community.number)
total_communities <- length(unique(removezerodata_meta_NA$community.number))
test_size <- ceiling(total_communities * 0.7)
test_communities <- sample(levels(removezerodata_meta_NA$community.number), size = test_size)
Community_Group_TV <- removezerodata_meta_NA
Community_Group_TV$set <- "validation"
Community_Group_TV[Community_Group_TV$community.number %in% test_communities, ]$set <- "test"
Community_Group_Test <- Community_Group_TV[Community_Group_TV$set == "test", ]
Community_Group_Validations <- Community_Group_TV[Community_Group_TV$set == "validation", ]


Community_Group_Test_clean_Chem_Com <- Community_Group_Test %>%
  select(-SampleID,  -community.name, 
         -chem.code, -complexity, -well, -plate,
         -set)
Community_Group_Test_clean_Chem_Com$community.number <- as.factor(Community_Group_Test_clean_Chem_Com$community.number)
Community_Group_Test_clean_Chem_Com$Amoxicillin <- as.factor(Community_Group_Test_clean_Chem_Com$Amoxicillin)
Community_Group_Test_clean_Chem_Com$Chlorothalonil <- as.factor(Community_Group_Test_clean_Chem_Com$Chlorothalonil)
Community_Group_Test_clean_Chem_Com$Diflufenican <- as.factor(Community_Group_Test_clean_Chem_Com$Diflufenican)
Community_Group_Test_clean_Chem_Com$Glyphosate <- as.factor(Community_Group_Test_clean_Chem_Com$Glyphosate)
Community_Group_Test_clean_Chem_Com$Imidacloprid <- as.factor(Community_Group_Test_clean_Chem_Com$Imidacloprid)
Community_Group_Test_clean_Chem_Com$Metaldehyde <- as.factor(Community_Group_Test_clean_Chem_Com$Metaldehyde)
Community_Group_Test_clean_Chem_Com$Oxytetracycline <- as.factor(Community_Group_Test_clean_Chem_Com$Oxytetracycline)
Community_Group_Test_clean_Chem_Com$Tebuconazole <- as.factor(Community_Group_Test_clean_Chem_Com$Tebuconazole)
str(Community_Group_Test_clean_Chem_Com)
#Community_group_RFSRC
RFSRC_Community_Chem_Com <- rfsrc(Multivar(community.number, Amoxicillin, Chlorothalonil, Diflufenican, 
                                           Glyphosate, Imidacloprid, Metaldehyde, Oxytetracycline, Tebuconazole) ~ .,
                         data = Community_Group_Test_clean_Chem_Com, ntree = 501, proximity = TRUE)

summary(RFSRC_Community_Chem_Com)
get.mv.error(RFSRC_Community_Chem_Com, standardize = TRUE)


Community_Group_Validations_clean_Chem_Com <- Community_Group_Validations %>%
  select(-SampleID, -community.name, -chem.code, 
         -complexity, -well, -plate,
         -set)
Community_Group_Validations_clean_Chem_Com$community.number <- as.factor(Community_Group_Validations_clean_Chem_Com$community.number)
Community_Group_Validations_clean_Chem_Com$Amoxicillin <- as.factor(Community_Group_Validations_clean_Chem_Com$Amoxicillin)
Community_Group_Validations_clean_Chem_Com$Chlorothalonil <- as.factor(Community_Group_Validations_clean_Chem_Com$Chlorothalonil)
Community_Group_Validations_clean_Chem_Com$Diflufenican <- as.factor(Community_Group_Validations_clean_Chem_Com$Diflufenican)
Community_Group_Validations_clean_Chem_Com$Glyphosate <- as.factor(Community_Group_Validations_clean_Chem_Com$Glyphosate)
Community_Group_Validations_clean_Chem_Com$Imidacloprid <- as.factor(Community_Group_Validations_clean_Chem_Com$Imidacloprid)
Community_Group_Validations_clean_Chem_Com$Metaldehyde <- as.factor(Community_Group_Validations_clean_Chem_Com$Metaldehyde)
Community_Group_Validations_clean_Chem_Com$Oxytetracycline <- as.factor(Community_Group_Validations_clean_Chem_Com$Oxytetracycline)
Community_Group_Validations_clean_Chem_Com$Tebuconazole <- as.factor(Community_Group_Validations_clean_Chem_Com$Tebuconazole)

Predictions_Community_Group_Validations_clean_Chem_Com <- predict(RFSRC_Community_Chem_Com, newdata = Community_Group_Validations_clean_Chem_Com, type = "response")
Predictions_Community_Group_Validations_clean_Chem_Com 
str(Predictions_Community_Group_Validations_clean_Chem_Com )
#Confusion matrix
Response_Community_Group_Chem_Com <- c("community.number", "Amoxicillin", "Chlorothalonil", "Diflufenican", "Glyphosate", "Imidacloprid",
                                       "Metaldehyde", "Oxytetracycline", "Tebuconazole")
for (response_var in Response_Community_Group_Chem_Com) {
  predictions_com_chem_com <- Predictions_Community_Group_Validations_clean_Chem_Com$classOutput[[response_var]]$class
  actual_cm_comg_chem_com <- Community_Group_Validations_clean_Chem_Com[[response_var]]

  confusion_matrix <- confusionMatrix(predictions_com_chem_com, actual_cm_comg_chem_com)
  
  print(paste("Confusion Matrix for", response_var))
  print(confusion_matrix)
} 


#ROC and AUC
print(str(Predictions_Community_Group_Validations_clean_Chem_Com$classOutput))
auc_values <- list()


for (var_name in Response_Community_Group_Chem_Com) {

  actual_value <- as.factor(Community_Group_Validations_clean_Chem_Com[[var_name]])

  prediction_prob <- Predictions_Community_Group_Validations_clean_Chem_Com$classOutput[[var_name]]$predicted

  if (length(actual_value) == nrow(prediction_prob)) {
 
    if (ncol(prediction_prob) == 2) {
      roc_curve <- roc(actual_value, prediction_prob[, 2], levels = rev(levels(actual_value)))
      auc_value <- auc(roc_curve)
    } else {
     
      roc_curve <- multiclass.roc(actual_value, as.numeric(Predictions_Community_Group_Validations_clean_Chem_Com$classOutput[[var_name]]$class))
      auc_value <- auc(roc_curve)
    }
    
  
    auc_values[[var_name]] <- auc_value
    

    print(paste("AUC for", var_name, "is", auc_value))
  } else {
    print(paste("Lengths of actual values and predictions do not match for", var_name))
  }
}
#Remove some cate
auc_values <- list()

for (var_name in Response_Community_Group_Chem_Com) {

  actual_value <- as.factor(Community_Group_Validations_clean_Chem_Com[[var_name]])
  

  prediction_prob <- Predictions_Community_Group_Validations_clean_Chem_Com$classOutput[[var_name]]$predicted
  
  if (length(actual_value) == nrow(prediction_prob)) {
 
    if (ncol(prediction_prob) == 2) {
      roc_curve <- roc(actual_value, prediction_prob[, 2], levels = rev(levels(actual_value)))
      auc_value <- auc(roc_curve)
    } else {

      present_levels <- levels(actual_value)
      
    
      actual_value <- factor(actual_value, levels = intersect(present_levels, colnames(prediction_prob)))
      

      roc_curve <- multiclass.roc(actual_value, as.numeric(Predictions_Community_Group_Validations_clean_Chem_Com$classOutput[[var_name]]$class))
      auc_value <- auc(roc_curve)
    }
    

    auc_values[[var_name]] <- auc_value

    print(paste("AUC for", var_name, "is", auc_value))
  } else {
    print(paste("Lengths of actual values and predictions do not match for", var_name))
  }
}


#No_com
Community_Group_Test_clean_Chem <- Community_Group_Test %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code, -complexity, -well, -plate,
         -set)
view(Community_Group_Test_clean_Chem)
Community_Group_Test_clean_Chem$Amoxicillin <- as.factor(Community_Group_Test_clean_Chem$Amoxicillin)
Community_Group_Test_clean_Chem$Chlorothalonil <- as.factor(Community_Group_Test_clean_Chem$Chlorothalonil)
Community_Group_Test_clean_Chem$Diflufenican <- as.factor(Community_Group_Test_clean_Chem$Diflufenican)
Community_Group_Test_clean_Chem$Glyphosate <- as.factor(Community_Group_Test_clean_Chem$Glyphosate)
Community_Group_Test_clean_Chem$Imidacloprid <- as.factor(Community_Group_Test_clean_Chem$Imidacloprid)
Community_Group_Test_clean_Chem$Metaldehyde <- as.factor(Community_Group_Test_clean_Chem$Metaldehyde)
Community_Group_Test_clean_Chem$Oxytetracycline <- as.factor(Community_Group_Test_clean_Chem$Oxytetracycline)
Community_Group_Test_clean_Chem$Tebuconazole <- as.factor(Community_Group_Test_clean_Chem$Tebuconazole)
str(Community_Group_Test_clean_Chem)
#Community_group_RFSRC
RFSRC_Community_Chem<- rfsrc(Multivar(Amoxicillin, Chlorothalonil, Diflufenican, 
                                           Glyphosate, Imidacloprid, Metaldehyde, Oxytetracycline, Tebuconazole) ~ .,
                                  data = Community_Group_Test_clean_Chem, ntree = 501, proximity = TRUE, importance = TRUE)

summary(RFSRC_Community_Chem)
get.mv.error(RFSRC_Community_Chem, standardize = TRUE)

vimp_RFSRC_Community_Chem <- vimp(RFSRC_Community_Chem)
importance_RFSRC_Community_Chem <- as.data.frame(vimp_RFSRC_Community_Chem)
str(RFSRC_Community_Chem)

print(importance_RFSRC_Community_Chem)
str(vimp_RFSRC_Community_Chem)
Community_Group_Validations_clean_Chem <- Community_Group_Validations %>%
  select(-SampleID, -community.name, -community.number, -chem.code, 
         -complexity, -well, -plate,
         -set)

Community_Group_Validations_clean_Chem$Amoxicillin <- as.factor(Community_Group_Validations_clean_Chem$Amoxicillin)
Community_Group_Validations_clean_Chem$Chlorothalonil <- as.factor(Community_Group_Validations_clean_Chem$Chlorothalonil)
Community_Group_Validations_clean_Chem$Diflufenican <- as.factor(Community_Group_Validations_clean_Chem$Diflufenican)
Community_Group_Validations_clean_Chem$Glyphosate <- as.factor(Community_Group_Validations_clean_Chem$Glyphosate)
Community_Group_Validations_clean_Chem$Imidacloprid <- as.factor(Community_Group_Validations_clean_Chem$Imidacloprid)
Community_Group_Validations_clean_Chem$Metaldehyde <- as.factor(Community_Group_Validations_clean_Chem$Metaldehyde)
Community_Group_Validations_clean_Chem$Oxytetracycline <- as.factor(Community_Group_Validations_clean_Chem$Oxytetracycline)
Community_Group_Validations_clean_Chem$Tebuconazole <- as.factor(Community_Group_Validations_clean_Chem$Tebuconazole)

Predictions_Community_Group_Validations_clean_Chem <- predict(RFSRC_Community_Chem, newdata = Community_Group_Validations_clean_Chem, type = "response")
Predictions_Community_Group_Validations_clean_Chem 


Response_Community_Group_Chem <- c("Amoxicillin", "Chlorothalonil", "Diflufenican", "Glyphosate", "Imidacloprid",
                                   "Metaldehyde", "Oxytetracycline", "Tebuconazole")
auc_values <- list()


for (var_name in Response_Community_Group_Chem) {

  actual_value <- as.factor(Community_Group_Validations_clean_Chem[[var_name]])
  
  prediction_prob <- Predictions_Community_Group_Validations_clean_Chem$classOutput[[var_name]]$predicted

  if (length(actual_value) == nrow(prediction_prob)) {
    
    roc_curve <- roc(actual_value, prediction_prob[, 2], levels = rev(levels(actual_value)), direction = "<")
    auc_value <- auc(roc_curve)

    auc_values[[var_name]] <- auc_value

    print(paste("AUC for", var_name, "is", auc_value))
  }
}
Predictions_Community_Group_Validations_clean_Chem$classOutput$Oxytetracycline$predicted

#CF
for (var in Response_Community_Group_Chem) {
  actual_values <- as.factor(Community_Group_Validations_clean_Chem[[var]])
  predictions_raw <- Predictions_Community_Group_Validations_clean_Chem$classOutput[[var]]$predicted
  num_classes <- length(levels(actual_values))
  if (num_classes > 2) {

    predictions <- factor(predictions_raw, levels = levels(actual_values))
  } else {
    predictions <- ifelse(predictions_raw[, 1] > predictions_raw[, 2], 0, 1)
    predictions <- factor(predictions, levels = levels(actual_values))
  }

  if (length(actual_values) != length(predictions)) {
    print(paste("Length mismatch for", var))
    print(paste("Length of actual_values:", length(actual_values)))
    print(paste("Length of predictions:", length(predictions)))
    next
  }

  confusion_matrix <- confusionMatrix(predictions, actual_values)

  print(paste("Confusion Matrix for", var))
  print(confusion_matrix)
}





#Not delete community.number

RFSRC_Community_Chem_NC<- rfsrc(Multivar(Amoxicillin, Chlorothalonil, Diflufenican, 
                                      Glyphosate, Imidacloprid, Metaldehyde, Oxytetracycline, Tebuconazole) ~ .,
                             data = Community_Group_Test_clean_Chem_Com, ntree = 501, proximity = TRUE)

summary(RFSRC_Community_Chem_NC)
get.mv.error(RFSRC_Community_Chem_NC, standardize = TRUE)

Predictions_Community_Group_Validations_clean_Chem_NC <- predict(RFSRC_Community_Chem_NC, newdata = Community_Group_Validations_clean_Chem_Com, type = "response")
Predictions_Community_Group_Validations_clean_Chem 
auc_values <- list()


for (var_name in Response_Community_Group_Chem) {
 
  actual_value <- as.factor(Community_Group_Validations_clean_Chem_Com[[var_name]])
  

  prediction_prob <- Predictions_Community_Group_Validations_clean_Chem_NC$classOutput[[var_name]]$predicted
  
 
  if (length(actual_value) == nrow(prediction_prob)) {
    
    roc_curve <- roc(actual_value, prediction_prob[, 2], levels = rev(levels(actual_value)))
    auc_value <- auc(roc_curve)

    auc_values[[var_name]] <- auc_value

    print(paste("AUC for", var_name, "is", auc_value))
  }
}



param_grid <- expand.grid(
  ntree = seq(100, 1000, by = 100),
  maxdepth = seq(5, 20, by = 5),
  nodesize = c(1, 5, 10)
)

best_performance <- double.infinity
best_params <- NULL


for (i in 1:nrow(param_grid)) {
  set.seed(123)
  model <- rfsrc(Multivar(community.number, Amoxicillin, Chlorothalonil, Diflufenican, 
                          Glyphosate, Imidacloprid, Metaldehyde, Oxytetracycline, Tebuconazole) ~ .,
                 data = Community_Group_Test_clean_Chem_Com, 
                 ntree = param_grid$ntree[i],
                 maxdepth = param_grid$maxdepth[i],
                 nodesize = param_grid$nodesize[i],
                 proximity = TRUE)
  
  performance <- 1 - model$err.rate[1]  
  

  if (performance < best_performance) {
    best_performance <- performance
    best_params <- param_grid[i, ]
  }
}


print("Best parameters found:")
print(best_params)




