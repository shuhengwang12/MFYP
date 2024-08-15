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
library(reshape2)
library(extrafont)
font_import(pattern = "Arial")
#1
set.seed(123)
removezerodata_meta_NA$community.number <- as.factor(removezerodata_meta_NA$community.number)
removezerodata_meta_NA$chem.code <- as.factor(removezerodata_meta_NA$chem.code)
Chem_Com_group <- removezerodata_meta_NA %>%
  group_by(community.number, chem.code) %>%
  mutate(set = cur_group_id()) %>%
  ungroup() 
G_num <- max(unique(Chem_Com_group$set))
Train_data <- ceiling(G_num * 0.7)
test_community_g <- sample(unique(Chem_Com_group$set), size = Train_data, replace = FALSE)
view(removezerodata_meta)
Community_Chem_G_Set <- Chem_Com_group
Community_Chem_G_Set$group <- "validation"
Community_Chem_G_Set[Community_Chem_G_Set$set %in% test_community_g, ]$group <- "test"
Community_Chem_G_Test <- Community_Chem_G_Set[Community_Chem_G_Set$group == "test", ]
Community_Chem_G_Val <- Community_Chem_G_Set[Community_Chem_G_Set$group == "validation", ]
vrow_indices <- match(Community_Chem_G_Val$SampleID, removezerodata_meta_NA$SampleID)
row_indices <- match(Community_Chem_G_Test$SampleID, removezerodata_meta_NA$SampleID)
t <- removezerodata_meta_NA[row_indices, ]
v <- removezerodata_meta_NA[vrow_indices, ]

Train_set_Chem_Com <- t %>%
  select(-SampleID,  -community.name, 
         -chem.code, -complexity, -well, -plate)
Train_set_Chem_Com$community.number <- as.factor(Train_set_Chem_Com$community.number)
Train_set_Chem_Com$Amoxicillin <- as.factor(Train_set_Chem_Com$Amoxicillin)
Train_set_Chem_Com$Chlorothalonil <- as.factor(Train_set_Chem_Com$Chlorothalonil)
Train_set_Chem_Com$Diflufenican <- as.factor(Train_set_Chem_Com$Diflufenican)
Train_set_Chem_Com$Glyphosate <- as.factor(Train_set_Chem_Com$Glyphosate)
Train_set_Chem_Com$Imidacloprid <- as.factor(Train_set_Chem_Com$Imidacloprid)
Train_set_Chem_Com$Metaldehyde <- as.factor(Train_set_Chem_Com$Metaldehyde)
Train_set_Chem_Com$Oxytetracycline <- as.factor(Train_set_Chem_Com$Oxytetracycline)
Train_set_Chem_Com$Tebuconazole <- as.factor(Train_set_Chem_Com$Tebuconazole)
str(Train_set_Chem_Com)
#Community_group_RFSRC
RFSRC_One_Chem_Com <- rfsrc(Multivar(community.number, Amoxicillin, Chlorothalonil, Diflufenican, 
                                           Glyphosate, Imidacloprid, Metaldehyde, Oxytetracycline, Tebuconazole) ~ .,
                                  data = Train_set_Chem_Com, ntree = 501, proximity = TRUE)

summary(RFSRC_One_Chem_Com)
get.mv.error(RFSRC_One_Chem_Com, standardize = TRUE)
Val_set_Chem_Com <- v %>%
  select(-SampleID, -community.name, -chem.code, 
         -complexity, -well,-plate)
Val_set_Chem_Com$community.number <- as.factor(Val_set_Chem_Com$community.number)
Val_set_Chem_Com$Amoxicillin <- as.factor(Val_set_Chem_Com$Amoxicillin)
Val_set_Chem_Com$Chlorothalonil <- as.factor(Val_set_Chem_Com$Chlorothalonil)
Val_set_Chem_Com$Diflufenican <- as.factor(Val_set_Chem_Com$Diflufenican)
Val_set_Chem_Com$Glyphosate <- as.factor(Val_set_Chem_Com$Glyphosate)
Val_set_Chem_Com$Imidacloprid <- as.factor(Val_set_Chem_Com$Imidacloprid)
Val_set_Chem_Com$Metaldehyde <- as.factor(Val_set_Chem_Com$Metaldehyde)
Val_set_Chem_Com$Oxytetracycline <- as.factor(Val_set_Chem_Com$Oxytetracycline)
Val_set_Chem_Com$Tebuconazole <- as.factor(Val_set_Chem_Com$Tebuconazole)

Predictions_Val_set_Chem_Com <- predict(RFSRC_One_Chem_Com, newdata = Val_set_Chem_Com, type = "response")
Predictions_Val_set_Chem_Com 


#ROC and AUC-factor
auc_values_one <- list()

response_var_one <- c("community.number", "Amoxicillin", "Chlorothalonil", "Diflufenican", "Glyphosate", "Imidacloprid", "Metaldehyde", "Oxytetracycline", "Tebuconazole")

for (var in response_var_one) {
  actual_value_one <- as.factor(Val_set_Chem_Com[[var]])
  if (length(levels(actual_value_one)) == 2) {
    prediction_one <- Predictions_Val_set_Chem_Com$classOutput[[var]]$predicted[, 2]
    roc_curve <- roc(actual_value_one, prediction_one, direction = "<")
    auc_value <- auc(roc_curve)
    auc_values_one[[var]] <- auc_value
    print(paste("AUC for", var, "is", auc_value))
  
    
  }  else {
    prediction_one <- Predictions_Val_set_Chem_Com$classOutput[[var]]$predicted
    roc_list <- lapply(1:ncol(prediction_one), function(i) {
      multiclass.roc(actual_value_one, prediction_one[, i])
    })
    auc_values <- sapply(roc_list, function(roc_obj) auc(roc_obj))
    auc_values_one[[var]] <- auc_values
    levels <- levels(actual_value_one)
    for (i in seq_along(auc_values)) {
      print(paste("AUC for", var, "class", levels[i], "is", auc_values[i]))
    }
  }
}
  
  
#one AUC

auc_values_one <- list()

response_var_one <- c("community.number", "Amoxicillin", "Chlorothalonil", "Diflufenican", "Glyphosate", "Imidacloprid", "Metaldehyde", "Oxytetracycline", "Tebuconazole")

for (var in response_var_one) {
  actual_value_one <- as.factor(Val_set_Chem_Com[[var]])
  
  if (length(levels(actual_value_one)) == 2) {
    prediction_one <- Predictions_Val_set_Chem_Com$classOutput[[var]]$predicted[, 2]
      roc_curve <- roc(actual_value_one, prediction_one,direction = "<")
      auc_value <- auc(roc_curve)
      auc_values_one[[var]] <- auc_value
      print(paste("AUC for", var, "is", auc_value))
    
  } else {
  
    prediction_one <- Predictions_Val_set_Chem_Com$classOutput[[var]]$predicted
      roc_curve <- multiclass.roc(actual_value_one, prediction_one)
      auc_value <- auc(roc_curve)
      auc_values_one[[var]] <- auc_value
      print(paste("AUC for", var, "is", auc_value))
    } 
  }

Predictions_Val_set_Chem_Com$classOutput$Oxytetracycline$predicted
#CF
for (var in response_var_one) {
  actual_values <- as.factor(Val_set_Chem_Com[[var]])
  predictions_raw <- Predictions_Val_set_Chem_Com$classOutput[[var]]$predicted
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

#ACC-Com
actual_classes_one_com <- as.character(Val_set_Chem_Com$community.number)

predicted_classes_com_one <- apply(Predictions_Val_set_Chem_Com$classOutput$community.number$predicted, 1, which.max)
predicted_classes_com_one  <- as.character(predicted_classes_com_one )
adjusted_predicted_classes_com_one <- as.character(ifelse(predicted_classes_com_one  != "1", as.numeric(predicted_classes_com_one) + 2, as.character(predicted_classes_com_one)))
accuracy_one_com <- sum(adjusted_predicted_classes_com_one  == actual_classes_one_com) / length(adjusted_predicted_classes_com_one)




#Cross Val
Cross_Val_One <- removezerodata_meta_NA
Cross_Val_One <- Cross_Val_One %>%
  select(-SampleID,  -community.name, 
         -chem.code, -complexity, -well, -plate)
Cross_Val_One$community.number <- as.factor(Cross_Val_One$community.number)
Cross_Val_One$Amoxicillin <- as.factor(Cross_Val_One$Amoxicillin)
Cross_Val_One$Chlorothalonil <- as.factor(Cross_Val_One$Chlorothalonil)
Cross_Val_One$Diflufenican <- as.factor(Cross_Val_One$Diflufenican)
Cross_Val_One$Glyphosate <- as.factor(Cross_Val_One$Glyphosate)
Cross_Val_One$Imidacloprid <- as.factor(Cross_Val_One$Imidacloprid)
Cross_Val_One$Metaldehyde <- as.factor(Cross_Val_One$Metaldehyde)
Cross_Val_One$Oxytetracycline <- as.factor(Cross_Val_One$Oxytetracycline)
Cross_Val_One$Tebuconazole <- as.factor(Cross_Val_One$Tebuconazole)

train_control_one <- trainControl(method = "cv", number = 10)  
Ten_fold_CV_data_one <- removezerodata_meta_NA  %>%
  select(-SampleID, -community.name, -chem.code, 
         -well, -plate)

set.seed(123)  
RFSRC_One_Chem_Com_cv <- train(Multivar(community.number, Amoxicillin, Chlorothalonil, Diflufenican,
                                              Glyphosate, Imidacloprid, Metaldehyde, Oxytetracycline, Tebuconazole) ~ .,
                                     data = Ten_fold_CV_data_one,
                                     method = "rfsrc", 
                                     ntree = 501,
                                     proximity = TRUE,
                                     trControl = train_control_one)




response_factors <- lapply(Ten_fold_CV_data_one[, c("community.number", "Amoxicillin", "Chlorothalonil", "Diflufenican",
                                                    "Glyphosate", "Imidacloprid", "Metaldehyde", "Oxytetracycline", "Tebuconazole")],
                           as.factor)

rfsrc_custom <- list(
  type = "Classification", 
  library = "randomForestSRC",
  loop = NULL,
  

  fit = function(x, y, wts, param, lev, last, classProbs, ...) {
    rfsrc(y ~ ., data = cbind(y, x), ntree = param$ntree, ...)
  },
  

  predict = function(modelFit, newdata, submodels = NULL) {
    predict(modelFit, newdata)$predicted
  },
  
  
  grid = function(x, y, len = NULL, search = "grid") {
    data.frame(ntree = 501)  
  },
  

  parameters = data.frame(parameter = "ntree", class = "numeric", label = "Number of Trees"),
  

  prob = NULL
)


train_control_one <- trainControl(method = "cv", number = 10)
set.seed(123)


models <- lapply(names(response_factors), function(response_var) {
  predictors <- Ten_fold_CV_data_one[, !(names(Ten_fold_CV_data_one) %in% c("community.number", "Amoxicillin", "Chlorothalonil", "Diflufenican",
                                                                            "Glyphosate", "Imidacloprid", "Metaldehyde", "Oxytetracycline", "Tebuconazole"))]
  model <- train(x = predictors, 
                 y = response_factors[[response_var]],  
                 method = rfsrc_custom, 
                 trControl = train_control_one)
  model
})


print(models[[1]])




#2 Done
#3
Train_data_three <- removezerodata_meta_NA %>%
  filter(complexity %in% c(0, 1, 2))

Test_data_three <- removezerodata_meta_NA %>%
  filter(complexity %in% c(4,8))

Train_set_Chem_Com_three <- Train_data_three %>%
  select(-SampleID,  -community.name, -complexity, -community.number,
         -chem.code,  -well, -plate)
view(Train_set_Chem_Com_three)
Train_set_Chem_Com_three$Amoxicillin <- as.factor(Train_set_Chem_Com_three$Amoxicillin)
Train_set_Chem_Com_three$Chlorothalonil <- as.factor(Train_set_Chem_Com_three$Chlorothalonil)
Train_set_Chem_Com_three$Diflufenican <- as.factor(Train_set_Chem_Com_three$Diflufenican)
Train_set_Chem_Com_three$Glyphosate <- as.factor(Train_set_Chem_Com_three$Glyphosate)
Train_set_Chem_Com_three$Imidacloprid <- as.factor(Train_set_Chem_Com_three$Imidacloprid)
Train_set_Chem_Com_three$Metaldehyde <- as.factor(Train_set_Chem_Com_three$Metaldehyde)
Train_set_Chem_Com_three$Oxytetracycline <- as.factor(Train_set_Chem_Com_three$Oxytetracycline)
Train_set_Chem_Com_three$Tebuconazole <- as.factor(Train_set_Chem_Com_three$Tebuconazole)

str(Train_set_Chem_Com_three)
#Community_group_RFSRC
RFSRC_Three_Chem_Com <- rfsrc(Multivar(Amoxicillin, Chlorothalonil, Diflufenican, 
                                     Glyphosate, Imidacloprid, Metaldehyde, Oxytetracycline, Tebuconazole) ~ .,
                            data = Train_set_Chem_Com_three, ntree = 501, proximity = TRUE)

summary(RFSRC_Three_Chem_Com)
get.mv.error(RFSRC_Three_Chem_Com, standardize = TRUE)
  
  
Val_set_Chem_Com_three <- Test_data_three  %>%
  select(-SampleID, -community.name, -chem.code,  -complexity, -community.number,
         -well, -plate)
view(Val_set_Chem_Com_three)
count_oxy_zeros_three <- sum(Train_set_Chem_Com_three$Oxytetracycline == 0)

Predictions_Val_set_Chem_Com_three <- predict(RFSRC_Three_Chem_Com, newdata = Val_set_Chem_Com_three, type = "response")
Predictions_Val_set_Chem_Com_three 
#ROC and AUC-factor
auc_values_three <- list()

response_var_three <- c("Amoxicillin", "Chlorothalonil", "Diflufenican", "Glyphosate", "Imidacloprid", "Metaldehyde", "Oxytetracycline", "Tebuconazole")

for (var in response_var_three) {
  actual_value_three <- as.factor(Val_set_Chem_Com_three[[var]])
    prediction_three <- Predictions_Val_set_Chem_Com_three$classOutput[[var]]$predicted[,2]
    roc_curve <- roc(actual_value_three, prediction_three, direction = "<")
    auc_value <- auc(roc_curve)
    auc_values_three[[var]] <- auc_value
    print(paste("AUC for", var, "is", auc_value))
    } 

#CF
for (var in response_var_three) {
  actual_values <- as.factor(Val_set_Chem_Com_three[[var]])
  predictions_raw <- Predictions_Val_set_Chem_Com_three$classOutput[[var]]$predicted
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







#3-one and zero
Train_data_three_one <- removezerodata_meta_NA %>%
  filter(complexity %in% c(0, 1))

Test_data_three_one <- removezerodata_meta_NA %>%
  filter(complexity %in% c(2,4,8))

Train_set_Chem_Com_three_one <- Train_data_three_one %>%
  select(-SampleID,  -community.name, -complexity,
         -chem.code,  -well, -plate)
Train_set_Chem_Com_three_one$community.number <- as.factor(Train_set_Chem_Com_three_one$community.number)
Train_set_Chem_Com_three_one$Amoxicillin <- as.factor(Train_set_Chem_Com_three_one$Amoxicillin)
Train_set_Chem_Com_three_one$Chlorothalonil <- as.factor(Train_set_Chem_Com_three_one$Chlorothalonil)
Train_set_Chem_Com_three_one$Diflufenican <- as.factor(Train_set_Chem_Com_three_one$Diflufenican)
Train_set_Chem_Com_three_one$Glyphosate <- as.factor(Train_set_Chem_Com_three_one$Glyphosate)
Train_set_Chem_Com_three_one$Imidacloprid <- as.factor(Train_set_Chem_Com_three_one$Imidacloprid)
Train_set_Chem_Com_three_one$Metaldehyde <- as.factor(Train_set_Chem_Com_three_one$Metaldehyde)
Train_set_Chem_Com_three_one$Oxytetracycline <- as.factor(Train_set_Chem_Com_three_one$Oxytetracycline)
Train_set_Chem_Com_three_one$Tebuconazole <- as.factor(Train_set_Chem_Com_three_one$Tebuconazole)

str(Train_set_Chem_Com_three_one)
#Community_group_RFSRC
RFSRC_Three_Chem_Com_one <- rfsrc(Multivar(Amoxicillin, Chlorothalonil, Diflufenican, 
                                       Glyphosate, Imidacloprid, Metaldehyde, Oxytetracycline, Tebuconazole) ~ .,
                              data = Train_set_Chem_Com_three_one, ntree = 501, proximity = TRUE)

summary(RFSRC_Three_Chem_Com_one)
get.mv.error(RFSRC_Three_Chem_Com_one, standardize = TRUE)


Val_set_Chem_Com_three_one <- Test_data_three_one  %>%
  select(-SampleID, -community.name, -chem.code, -complexity,
         -well, -plate)
Val_set_Chem_Com_three_one$community.number <- as.factor(Val_set_Chem_Com_three_one$community.number)
Val_set_Chem_Com_three_one$Amoxicillin <- as.factor(Val_set_Chem_Com_three_one$Amoxicillin)
Val_set_Chem_Com_three_one$Chlorothalonil <- as.factor(Val_set_Chem_Com_three_one$Chlorothalonil)
Val_set_Chem_Com_three_one$Diflufenican <- as.factor(Val_set_Chem_Com_three_one$Diflufenican)
Val_set_Chem_Com_three_one$Glyphosate <- as.factor(Val_set_Chem_Com_three_one$Glyphosate)
Val_set_Chem_Com_three_one$Imidacloprid <- as.factor(Val_set_Chem_Com_three_one$Imidacloprid)
Val_set_Chem_Com_three_one$Metaldehyde <- as.factor(Val_set_Chem_Com_three_one$Metaldehyde)
Val_set_Chem_Com_three_one$Oxytetracycline <- as.factor(Val_set_Chem_Com_three_one$Oxytetracycline)
Val_set_Chem_Com_three_one$Tebuconazole <- as.factor(Val_set_Chem_Com_three_one$Tebuconazole)

Predictions_Val_set_Chem_Com_three_one <- predict(RFSRC_Three_Chem_Com_one, newdata = Val_set_Chem_Com_three_one, type = "response")
Predictions_Val_set_Chem_Com_three_one
#ROC and AUC-factor
auc_values_three_one <- list()

response_var_three_one <- c("Amoxicillin", "Chlorothalonil", "Diflufenican", "Glyphosate", "Imidacloprid", "Metaldehyde", "Oxytetracycline", "Tebuconazole")

for (var in response_var_three_one) {
  actual_value_three_one <- as.factor(Val_set_Chem_Com_three_one[[var]])
  prediction_three_one <- Predictions_Val_set_Chem_Com_three_one$classOutput[[var]]$predicted[,2]
  roc_curve <- roc(actual_value_three_one, prediction_three_one)
  auc_value <- auc(roc_curve)
  auc_values_three_one[[var]] <- auc_value
  print(paste("AUC for", var, "is", auc_value))
} 


#Importance
var_importance_three_one <- vimp(RFSRC_Three_Chem_Com_one)


importance_three_one_amoxicillin <- var_importance_three_one$classOutput$Amoxicillin$importance

print(importance_three_one_amoxicillin)

importance_three_one_amo_df <- as.data.frame(importance_three_one_amoxicillin)
importance_three_one_amo_df$Variables <- rownames(importance_three_one_amo_df)
print(importance_three_one_amo_df)
importance_three_one_amo_df <- importance_three_one_amo_df[, c("Variables", colnames(importance_three_one_amo_df)[1:(ncol(importance_three_one_amo_df)-1)])]
importance_three_one_amo_df <- importance_three_one_amo_df[order(importance_three_one_amo_df$all, decreasing = TRUE), ]

ggplot(importance_three_one_amo_df, aes(x = reorder(Variables, all), y = all)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(x = "Importance", y = "Variables", title = "Variable Importance for Amoxicillin") +
  theme_minimal()

ggplot(importance_three_one_amo_df, aes(x = reorder(Variables, all, FUN = rev), y = all)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Variable Importance for Amoxicillin",
       x = "Variables",
       y = "Overall Importance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))

str(taxonomy)

asvs_three_one <- colnames(Train_set_Chem_Com_three_one)[1:92]
new_column_names_three_one <- ifelse(substr(asvs_three_one, 1, 1) == "X", substring(asvs_three_one, 2), asvs_three_one)
colnames(Train_set_Chem_Com_three_one)[1:92] <- new_column_names_three_one
asvs_three_one <- colnames(Train_set_Chem_Com_three_one)[1:92]

asv_indices_three_one <- match(asvs_three_one, rownames(taxonomy))

matched_taxonomy_three_one <- taxonomy[asv_indices_three_one, ]
matched_taxonomy_three_one$ASV <- rownames(matched_taxonomy_three_one)
rownames(matched_taxonomy_three_one) <- NULL


asvs_three_one_amo_im <- rownames(importance_three_one_amo_df)
new_column_names_three_one_amo_im <- ifelse(substr(asvs_three_one_amo_im, 1, 1) == "X", substring(asvs_three_one_amo_im, 2), asvs_three_one_amo_im)
rownames(importance_three_one_amo_df) <- new_column_names_three_one_amo_im
importance_three_one_amo_df$Variables <- rownames(importance_three_one_amo_df)




merged_data_three_one_amo <- merge(importance_three_one_amo_df, matched_taxonomy_three_one, by.x = "Variables", by.y = "ASV", all.x = TRUE)

merged_data_three_one_amo$Family[is.na(merged_data_three_one_amo$Family)] <- merged_data_three_one_amo$Variables[is.na(merged_data_three_one_amo$Family)]


amo_three_one_family_range <- merged_data_three_one_amo %>%
  group_by(Family) %>%
  summarise(Min = min(all, na.rm = TRUE), Max = max(all, na.rm = TRUE))


ggplot(amo_three_one_family_range, aes(x = Family, y = 0, color = Family)) +
  geom_point(aes(x = Family, y = Max), size = 2) +
  geom_errorbar(aes(ymin = Min, ymax = Max), width = 0.2, na.rm = TRUE) +
  coord_flip() +
  labs(x = "Family", y = "Importance Scores", title = "Range of Importance Scores by Family") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))











#K-Fold   
set.seed(123)  
train_control_three <- trainControl(method = "cv", number = 5) 
data_three <- removezerodata_meta_NA %>%
  select(-SampleID, -community.name, -chem.code, -complexity,
         -well, -plate)
data_three$community.number <- as.factor(data_three$community.number)
data_three$Amoxicillin <- as.factor(data_three$Amoxicillin)
data_three$Chlorothalonil <- as.factor(data_three$Chlorothalonil)
data_three$Diflufenican <- as.factor(data_three$Diflufenican)
data_three$Glyphosate <- as.factor(data_three$Glyphosate)
data_three$Imidacloprid <- as.factor(data_three$Imidacloprid)
data_three$Metaldehyde <- as.factor(data_three$Metaldehyde)
data_three$Oxytetracycline <- as.factor(data_three$Oxytetracycline)
data_three$Tebuconazole <- as.factor(data_three$Tebuconazole)
response_vars <- c("Amoxicillin", "Chlorothalonil", "Diflufenican", "Glyphosate", 
                   "Imidacloprid", "Metaldehyde", "Oxytetracycline", "Tebuconazole")

rfsrc_wrapper <- function(data_three, indices) {
  # Create training set
  train_data <- data_three[indices, ]
  # Create testing set
  test_data <- data_three[-indices, ]
  
  rfsrc_model <- rfsrc(Multivar(Amoxicillin, Chlorothalonil, Diflufenican, 
                                 Glyphosate, Imidacloprid, Metaldehyde, 
                                 Oxytetracycline, Tebuconazole) ~ .,
                        data = train_data, ntree = 501, proximity = TRUE)
  

  predictions <- predict(rfsrc_model, newdata = test_data)
 
  print(str(predictions))

  error_rate <- sapply(response_vars, function(var) {
    mean(predictions$class[[var]] != test_data[[var]])
  })
  
  return(list(model = rfsrc_model, error_rate = error_rate))
}


train_function <- function(data_three, method = "rfsrc", trControl, response_vars) {
  if (method == "rfsrc") {
    
    folds <- createFolds(data_three[[response_vars[1]]], k = trControl$number)
    resamples <- lapply(folds, function(indices) {
    
      rfsrc_result <- rfsrc_wrapper(data_three, indices)
     
      return(rfsrc_result)
    })
    return(resamples)
  } else {
    stop("Unsupported method")
  }
}



models <- train_function(data = data_three, method = "rfsrc", trControl = train_control_three, response_vars = response_vars)

get_model_performance <- function(models, metric = "error") {
  errors <- sapply(models, function(model_result) {
    if (!is.null(model_result$error_rate)) {
      return(mean(model_result$error_rate))
    } else {
      return(NA)
    }
  })
  print(errors)  
  return(mean(errors, na.rm = TRUE))
}


mean_error <- get_model_performance(models)
print(mean_error)

print(str(predictions))










#4-Oxy res
removezerodata_meta_NA$community.number <- as.factor(removezerodata_meta_NA$community.number)
removezerodata_meta_NA$chem.code <- as.factor(removezerodata_meta_NA$chem.code)
Chem_Com_group <- removezerodata_meta_NA %>%
  group_by(community.number, chem.code) %>%
  mutate(set = cur_group_id()) %>%
  ungroup() 
G_num <- max(unique(Chem_Com_group$set))
Train_data <- ceiling(G_num * 0.7)
test_community_g <- sample(unique(Chem_Com_group$set), size = Train_data, replace = FALSE)

Community_Chem_G_Set <- Chem_Com_group
Community_Chem_G_Set$group <- "validation"
Community_Chem_G_Set[Community_Chem_G_Set$set %in% test_community_g, ]$group <- "test"
Community_Chem_G_Test <- Community_Chem_G_Set[Community_Chem_G_Set$group == "test", ]
Community_Chem_G_Val <- Community_Chem_G_Set[Community_Chem_G_Set$group == "validation", ]
vrow_indices <- match(Community_Chem_G_Val$SampleID, removezerodata_meta_NA$SampleID)
row_indices <- match(Community_Chem_G_Test$SampleID, removezerodata_meta_NA$SampleID)
t_f <- removezerodata_meta_NA[row_indices, ]
v_f <- removezerodata_meta_NA[vrow_indices, ]

Train_set_Oxy <- t_f %>%
  select(-SampleID,  -community.number, -community.name, 
         -chem.code, -Amoxicillin, -Chlorothalonil, -Diflufenican, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, 
         -Tebuconazole, -complexity, -well, -plate)
view(t_f)
Train_set_Oxy$Oxytetracycline <- as.factor(Train_set_Oxy$Oxytetracycline)

str(Train_set_Oxy)
RFSRC_Four_Oxy_Res <- rfsrc(Oxytetracycline ~ ., data = Train_set_Oxy, ntree = 501, importance = TRUE, proximities = TRUE)
get.mv.error(RFSRC_Four_Oxy_Res, standardize = TRUE)

Val_set_Oxy <- v_f %>%
  select(-SampleID,  -community.number, -community.name, 
         -chem.code, -Amoxicillin, -Chlorothalonil, -Diflufenican, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, 
         -Tebuconazole, -complexity, -well, -plate)
Val_set_Oxy$Oxytetracycline <- as.factor(Val_set_Oxy$Oxytetracycline)
str(Val_set_Oxy)
Predictions_Four_Oxy_Res <- predict(RFSRC_Four_Oxy_Res, newdata = Val_set_Oxy, type = "response")
Predictions_Four_Oxy_Res
roc_curve_Oxy <- roc(Val_set_Oxy$Oxytetracycline, Predictions_Four_Oxy_Res$predicted[,2])

# 计算AUC值
auc_value_Oxy <- auc(roc_curve_Oxy)
auc_value_Oxy

#CF
response_var_four_Oxy <- c("Oxytetracycline")
for (var in response_var_four_Oxy) {
  actual_values <- as.factor(Val_set_Oxy[[var]])
  predictions_raw <- Predictions_Four_Oxy_Res$predicted

    predictions <- ifelse(predictions_raw[, 1] > predictions_raw[, 2], 0, 1)
    predictions <- factor(predictions, levels = levels(actual_values))
  }

  confusion_matrix <- confusionMatrix(predictions, actual_values)
  
  print(paste("Confusion Matrix for", var))
  print(confusion_matrix)


#Oxy 
err_rate_four_Oxy <- RFSRC_Four_Oxy_Res$err.rate

err_rate_four_Oxy_df <- as.data.frame(err_rate_four_Oxy)

model_ids_four_Oxy <- dimnames(err_rate_four_Oxy)[[2]]
colnames(err_rate_four_Oxy_df) <- model_ids_four_Oxy

err_rate_four_Oxy_long <- err_rate_four_Oxy_df %>%
  mutate(ntrees = 1:501) %>%
  pivot_longer(-ntrees, names_to = "model_id", values_to = "err_rate")

ggplot(err_rate_four_Oxy_long, aes(x = ntrees, y = err_rate, color = model_id, group = model_id)) +
  geom_line() +
  geom_point()+
  geom_smooth(method = "loess", se = FALSE) +
  labs(
    title = "Community's Error Rate by Number of Trees and Model",
    x = "Number of Trees",
    y = "Error Rate",
    color = "Model"
  ) +
  theme_minimal()


ggplot(err_rate_four_Oxy_long, aes(x = ntrees, y = err_rate, color = model_id, group = model_id)) +
  geom_line() +
  geom_point()+
  geom_smooth(method = "loess", se = FALSE) +
  labs(
    x = "Number of Trees",
    y = "Error Rate",
    color = "Model"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.ticks = element_line(color = "black", size = 1),
    panel.border = element_rect(colour = "black", fill = NA, size = 2),
    axis.ticks.length = unit(0.1, "cm"),
    axis.text = element_text(family = "sans", face = "bold", size = 10, color = "black"),
    axis.title = element_text(family = "sans", face = "bold", size = 11, color = "black"),
    legend.text = element_text(family = "sans", face = "bold", size = 11, color = "black"),
    legend.title = element_text(family = "sans", face = "bold", size = 11, color = "black"),
    legend.position = c(0.95, 0.95),  
    legend.justification = c(1, 1)  
  )









#Com_Oxy
Train_set_Oxy_Com <- t_f %>%
  select(-SampleID,  -community.name, 
         -chem.code, -Amoxicillin, -Chlorothalonil, -Diflufenican, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, 
         -Tebuconazole, -complexity, -well, -plate)

Train_set_Oxy_Com$Oxytetracycline <- as.factor(Train_set_Oxy_Com$Oxytetracycline)
Train_set_Oxy_Com$community.number <- as.factor(Train_set_Oxy_Com$community.number)

RFSRC_Four_Oxy_Res_Com <- rfsrc(Multivar(community.number, Oxytetracycline) ~ ., data = Train_set_Oxy_Com, ntree = 501, importance = TRUE, proximities = TRUE)
get.mv.error(RFSRC_Four_Oxy_Res_Com, standardize = TRUE)

Val_set_Oxy_Com <- v_f %>%
  select(-SampleID,  -community.name, 
         -chem.code, -Amoxicillin, -Chlorothalonil, -Diflufenican, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, 
         -Tebuconazole, -complexity, -well, -plate)
Val_set_Oxy_Com$Oxytetracycline <- as.factor(Val_set_Oxy_Com$Oxytetracycline)
Val_set_Oxy_Com$community.number <- as.factor(Val_set_Oxy_Com$community.number)

Predictions_Four_Oxy_Res_Com <- predict(RFSRC_Four_Oxy_Res_Com, newdata = Val_set_Oxy_Com, type = "response")
Predictions_Four_Oxy_Res_Com
#ROC and AUC-factor
auc_values_Four_Oxy_Com <- list()
response_var_four_Oxy_Com <- c("community.number",  "Oxytetracycline")


for (var in response_var_four_Oxy_Com) {
  actual_value_four_Oxy_Com <- as.factor(Val_set_Oxy_Com[[var]])
  
  if (length(levels(actual_value_four_Oxy_Com )) == 2) {
    prediction_four_Oxy_Com <- Predictions_Four_Oxy_Res_Com$classOutput[[var]]$predicted[, 2]
    roc_curve <- roc(actual_value_four_Oxy_Com, prediction_four_Oxy_Com)
    auc_value <- auc(roc_curve)
    auc_values_Four_Oxy_Com[[var]] <- auc_value
    print(paste("AUC for", var, "is", auc_value))
    
  } else {
    
    prediction_four_Oxy_Com <- Predictions_Four_Oxy_Res_Com$classOutput[[var]]$predicted
    roc_curve <- multiclass.roc(actual_value_four_Oxy_Com, prediction_four_Oxy_Com)
    auc_value <- auc(roc_curve)
    auc_values_Four_Oxy_Com[[var]] <- auc_value
    print(paste("AUC for", var, "is", auc_value))
  } 
}
#CF
for (var in response_var_four_Oxy_Com) {
  actual_values <- as.factor(Val_set_Oxy_Com[[var]])
  predictions_raw <- Predictions_Four_Oxy_Res_Com$classOutput[[var]]$predicted
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
#ACC-Com
actual_classes_oxy_com_com <- as.character(Val_set_Oxy_Com$community.number)

predicted_classes_com_oxy_com <- apply(Predictions_Four_Oxy_Res_Com$classOutput$community.number$predicted, 1, which.max)
predicted_classes_com_oxy_com<- as.character(predicted_classes_com_oxy_com)
adjusted_predicted_classes_com_oxy_com <- as.character(ifelse(predicted_classes_com_oxy_com != "1", as.numeric(predicted_classes_com_oxy_com) + 2, as.character(predicted_classes_com_oxy_com)))
accuracy_oxy_com_com <- sum(adjusted_predicted_classes_com_oxy_com   == actual_classes_oxy_com_com) / length(adjusted_predicted_classes_com_oxy_com )



#Com
err_rate_Oxy_Com_four_com <- RFSRC_Four_Oxy_Res_Com$classOutput$community.number$err.rate

err_rate_Oxy_Com_four_com_df <- as.data.frame(err_rate_Oxy_Com_four_com)

model_ids_Oxy_Com_four_com <- dimnames(err_rate_Oxy_Com_four_com)[[2]]
colnames(err_rate_Oxy_Com_four_com_df) <- model_ids_Oxy_Com_four_com

err_rate_Oxy_Com_four_com_long <- err_rate_Oxy_Com_four_com_df %>%
  mutate(ntrees = 1:501) %>%
  pivot_longer(-ntrees, names_to = "model_id", values_to = "err_rate")

ggplot(err_rate_Oxy_Com_four_com_long, aes(x = ntrees, y = err_rate, color = model_id, group = model_id)) +
  geom_line() +
  geom_point()+
  geom_smooth(method = "loess", se = FALSE) +
  labs(
    x = "Number of Trees",
    y = "Error Rate",
    color = "Model"
  ) +
  main_theme



#Oxy
err_rate_Oxy_Com_four_Oxy <- RFSRC_Four_Oxy_Res_Com$classOutput$Oxytetracycline$err.rate

err_rate_Oxy_Com_four_Oxy_df <- as.data.frame(err_rate_Oxy_Com_four_Oxy)

model_ids_Oxy_Com_four_Oxy <- dimnames(err_rate_Oxy_Com_four_Oxy)[[2]]
colnames(err_rate_Oxy_Com_four_Oxy_df) <- model_ids_Oxy_Com_four_Oxy

err_rate_Oxy_Com_four_Oxy_long <- err_rate_Oxy_Com_four_Oxy_df %>%
  mutate(ntrees = 1:501) %>%
  pivot_longer(-ntrees, names_to = "model_id", values_to = "err_rate")

ggplot(err_rate_Oxy_Com_four_Oxy_long, aes(x = ntrees, y = err_rate, color = model_id, group = model_id)) +
  geom_line() +
  geom_point()+
  geom_smooth(method = "loess", se = FALSE) +
  labs(
    x = "Number of Trees",
    y = "Error Rate",
    color = "Model"
  ) +
  main_theme

ggplot(err_rate_Oxy_Com_four_Oxy_long, aes(x = ntrees, y = err_rate, color = model_id, group = model_id)) +
  geom_line() +
  geom_point()+
  geom_smooth(method = "loess", se = FALSE) +
  labs(
    x = "Number of Trees",
    y = "Error Rate",
    color = "Model"
  ) +
  theme_minimal()


#Single no com 70% 30%

removezerodata_meta_NA$community.number <- as.factor(removezerodata_meta_NA$community.number)
total_communities_sin <- length(unique(removezerodata_meta_NA$community.number))
test_size_sin <- ceiling(total_communities_sin * 0.7)
train_communities_sin <- sample(levels(removezerodata_meta_NA$community.number), size = test_size_sin, replace = FALSE)
Community_Group_TV_sin <- removezerodata_meta_NA
Community_Group_TV_sin$set <- "validation"
Community_Group_TV_sin[Community_Group_TV_sin$community.number %in% train_communities_sin, ]$set <- "train"
Community_Group_Train_sin <- Community_Group_TV_sin[Community_Group_TV_sin$set == "train", ]
Community_Group_Validations_sin <- Community_Group_TV_sin[Community_Group_TV_sin$set == "validation", ]

Community_Group_Train_sin$community.number <- as.factor(Community_Group_Train_sin$community.number)
Community_Group_Train_sin$Amoxicillin <- as.factor(Community_Group_Train_sin$Amoxicillin)
Community_Group_Train_sin$Chlorothalonil <- as.factor(Community_Group_Train_sin$Chlorothalonil)
Community_Group_Train_sin$Diflufenican <- as.factor(Community_Group_Train_sin$Diflufenican)
Community_Group_Train_sin$Glyphosate <- as.factor(Community_Group_Train_sin$Glyphosate)
Community_Group_Train_sin$Imidacloprid <- as.factor(Community_Group_Train_sin$Imidacloprid)
Community_Group_Train_sin$Metaldehyde <- as.factor(Community_Group_Train_sin$Metaldehyde)
Community_Group_Train_sin$Oxytetracycline <- as.factor(Community_Group_Train_sin$Oxytetracycline)
Community_Group_Train_sin$Tebuconazole <- as.factor(Community_Group_Train_sin$Tebuconazole)

Community_Group_Validations_sin$community.number <- as.factor(Community_Group_Validations_sin$community.number)
Community_Group_Validations_sin$Amoxicillin <- as.factor(Community_Group_Validations_sin$Amoxicillin)
Community_Group_Validations_sin$Chlorothalonil <- as.factor(Community_Group_Validations_sin$Chlorothalonil)
Community_Group_Validations_sin$Diflufenican <- as.factor(Community_Group_Validations_sin$Diflufenican)
Community_Group_Validations_sin$Glyphosate <- as.factor(Community_Group_Validations_sin$Glyphosate)
Community_Group_Validations_sin$Imidacloprid <- as.factor(Community_Group_Validations_sin$Imidacloprid)
Community_Group_Validations_sin$Metaldehyde <- as.factor(Community_Group_Validations_sin$Metaldehyde)
Community_Group_Validations_sin$Oxytetracycline <- as.factor(Community_Group_Validations_sin$Oxytetracycline)
Community_Group_Validations_sin$Tebuconazole <- as.factor(Community_Group_Validations_sin$Tebuconazole)

#Amo
Community_Group_Train_sin_amo  <- Community_Group_Train_sin %>%
  select(-SampleID,  -community.name,  -community.number, -Chlorothalonil, -Diflufenican,
         -Glyphosate, -Imidacloprid, -Metaldehyde, -Oxytetracycline, -Tebuconazole,
         -chem.code, -complexity, -well, -plate,
         -set)
str(Community_Group_Train_sin_amo)
RFSRC_nocomsin_amo <- rfsrc(Amoxicillin ~ ., Community_Group_Train_sin_amo, ntree = 501, importance = TRUE, proximities = TRUE)
get.mv.error(RFSRC_nocomsin_amo)
Community_Group_Validations_sin_amo <- Community_Group_Validations_sin %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code,  -Chlorothalonil, -Diflufenican, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate,-set)

str(Community_Group_Validations_sin_amo)
Predictions_nocomsin_amo <- predict(RFSRC_nocomsin_amo, newdata = Community_Group_Validations_sin_amo, type = "response")

roc_curve_nocomsin_amo <- roc(Community_Group_Validations_sin_amo$Amoxicillin, Predictions_nocomsin_amo$predicted[,2])


auc_value_nocomsin_amo <- auc(roc_curve_nocomsin_amo)
auc_value_nocomsin_amo

#CF
0.672

#Chl
Community_Group_Train_sin_chl  <- Community_Group_Train_sin %>%
  select(-SampleID,  -community.name,  -community.number, -Amoxicillin, -Diflufenican,
         -Glyphosate, -Imidacloprid, -Metaldehyde, -Oxytetracycline, -Tebuconazole,
         -chem.code, -complexity, -well, -plate,
         -set)
str(Community_Group_Train_sin_chl)
RFSRC_nocomsin_chl <- rfsrc(Chlorothalonil ~ ., Community_Group_Train_sin_chl, ntree = 501, importance = TRUE, proximities = TRUE)
get.mv.error(RFSRC_nocomsin_chl)
Community_Group_Validations_sin_chl <- Community_Group_Validations_sin %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code,  -Amoxicillin, -Diflufenican, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate, -set)
str(Community_Group_Validations_sin_chl)
Predictions_nocomsin_chl<- predict(RFSRC_nocomsin_chl, newdata = Community_Group_Validations_sin_chl, type = "response")

roc_curve_nocomsin_chl <- roc(Community_Group_Validations_sin_chl$Chlorothalonil, Predictions_nocomsin_chl$predicted[,2])


auc_value_nocomsin_chl <- auc(roc_curve_nocomsin_chl)
auc_value_nocomsin_chl

#CF
0.716

#Dif
Community_Group_Train_sin_dif  <- Community_Group_Train_sin %>%
  select(-SampleID,  -community.name,  -community.number, -Amoxicillin, -Chlorothalonil,
         -Glyphosate, -Imidacloprid, -Metaldehyde, -Oxytetracycline, -Tebuconazole,
         -chem.code, -complexity, -well, -plate,
         -set)
str(Community_Group_Train_sin_dif)
RFSRC_nocomsin_dif <- rfsrc(Diflufenican ~ ., Community_Group_Train_sin_dif, ntree = 501, importance = TRUE, proximities = TRUE)
get.mv.error(RFSRC_nocomsin_dif)
Community_Group_Validations_sin_dif <- Community_Group_Validations_sin %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code,  -Amoxicillin, -Chlorothalonil, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate, -set)
str(Community_Group_Validations_sin_dif)
Predictions_nocomsin_dif <- predict(RFSRC_nocomsin_dif, newdata = Community_Group_Validations_sin_dif, type = "response")

roc_curve_nocomsin_dif <- roc(Community_Group_Validations_sin_dif$Diflufenican, Predictions_nocomsin_dif$predicted[,2])


auc_value_nocomsin_dif <- auc(roc_curve_nocomsin_dif)
auc_value_nocomsin_dif

#CF
0.657


#Gly
Community_Group_Train_sin_gly  <- Community_Group_Train_sin %>%
  select(-SampleID,  -community.name,  -community.number, -Amoxicillin, -Chlorothalonil,
         -Diflufenican, -Imidacloprid, -Metaldehyde, -Oxytetracycline, -Tebuconazole,
         -chem.code, -complexity, -well, -plate,
         -set)
str(Community_Group_Train_sin_gly)
RFSRC_nocomsin_gly <- rfsrc(Glyphosate ~ ., Community_Group_Train_sin_gly, ntree = 501, importance = TRUE, proximities = TRUE)
get.mv.error(RFSRC_nocomsin_gly)
Community_Group_Validations_sin_gly <- Community_Group_Validations_sin %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code,  -Amoxicillin, -Chlorothalonil, 
         -Diflufenican, -Imidacloprid, -Metaldehyde, -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate, -set)
str(Community_Group_Validations_sin_gly)
Predictions_nocomsin_gly <- predict(RFSRC_nocomsin_gly, newdata = Community_Group_Validations_sin_gly, type = "response")

roc_curve_nocomsin_gly <- roc(Community_Group_Validations_sin_gly$Glyphosate, Predictions_nocomsin_gly$predicted[,2])

auc_value_nocomsin_gly <- auc(roc_curve_nocomsin_gly)
auc_value_nocomsin_gly

#CF
0.701


#Imi
Community_Group_Train_sin_Imi  <- Community_Group_Train_sin %>%
  select(-SampleID,  -community.name,  -community.number, -Amoxicillin, -Chlorothalonil,
         -Diflufenican, -Glyphosate, -Metaldehyde, -Oxytetracycline, -Tebuconazole,
         -chem.code, -complexity, -well, -plate,
         -set)
str(Community_Group_Train_sin_Imi)
RFSRC_nocomsin_imi <- rfsrc(Imidacloprid ~ ., Community_Group_Train_sin_Imi, ntree = 501, importance = TRUE, proximities = TRUE)
get.mv.error(RFSRC_nocomsin_imi)
Community_Group_Validations_sin_imi <- Community_Group_Validations_sin %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code,  -Amoxicillin, -Chlorothalonil, 
         -Diflufenican, -Glyphosate, -Metaldehyde, -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate, -set)
str(Community_Group_Validations_sin_imi)
Predictions_nocomsin_imi <- predict(RFSRC_nocomsin_imi, newdata = Community_Group_Validations_sin_imi, type = "response")

roc_curve_nocomsin_imi <- roc(Community_Group_Validations_sin_imi$Imidacloprid, Predictions_nocomsin_imi$predicted[,2])


auc_value_nocomsin_imi <- auc(roc_curve_nocomsin_imi)
auc_value_nocomsin_imi

#CF
0.761

#Met
Community_Group_Train_sin_met  <- Community_Group_Train_sin %>%
  select(-SampleID,  -community.name,  -community.number, -Amoxicillin, -Chlorothalonil,
         -Diflufenican, -Glyphosate, -Imidacloprid, -Oxytetracycline, -Tebuconazole,
         -chem.code, -complexity, -well, -plate,
         -set)
str(Community_Group_Train_sin_met)
RFSRC_nocomsin_met <- rfsrc(Metaldehyde ~ ., Community_Group_Train_sin_met, ntree = 501, importance = TRUE, proximities = TRUE)
get.mv.error(RFSRC_nocomsin_met)
Community_Group_Validations_sin_met <- Community_Group_Validations_sin %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code,  -Amoxicillin, -Chlorothalonil, 
         -Diflufenican, -Glyphosate, -Imidacloprid, -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate, -set)
str(Community_Group_Validations_sin_met)
Predictions_nocomsin_met <- predict(RFSRC_nocomsin_met, newdata = Community_Group_Validations_sin_met, type = "response")

roc_curve_nocomsin_met <- roc(Community_Group_Validations_sin_met$Metaldehyde, Predictions_nocomsin_met$predicted[,2])


auc_value_nocomsin_met <- auc(roc_curve_nocomsin_met)
auc_value_nocomsin_met

#CF
0.716

#Oxy
Community_Group_Train_sin_oxy  <- Community_Group_Train_sin %>%
  select(-SampleID,  -community.name,  -community.number, -Amoxicillin, -Chlorothalonil,
         -Diflufenican, -Glyphosate, -Imidacloprid, -Metaldehyde, -Tebuconazole,
         -chem.code, -complexity, -well, -plate,
         -set)
str(Community_Group_Train_sin_oxy)
RFSRC_nocomsin_oxy <- rfsrc(Oxytetracycline ~ ., Community_Group_Train_sin_oxy, ntree = 501, importance = TRUE, proximities = TRUE)
get.mv.error(RFSRC_nocomsin_oxy)
Community_Group_Validations_sin_oxy <- Community_Group_Validations_sin %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code,  -Amoxicillin, -Chlorothalonil, 
         -Diflufenican, -Glyphosate, -Imidacloprid, -Metaldehyde,
         -Tebuconazole, -complexity, -well, -plate, -set)
str(Community_Group_Validations_sin_oxy)
Predictions_nocomsin_oxy <- predict(RFSRC_nocomsin_oxy, newdata = Community_Group_Validations_sin_oxy, type = "response")

roc_curve_nocomsin_oxy <- roc(Community_Group_Validations_sin_oxy$Oxytetracycline, Predictions_nocomsin_oxy$predicted[,2])


auc_value_nocomsin_oxy <- auc(roc_curve_nocomsin_oxy)
auc_value_nocomsin_oxy

#CF
0.940

#Teb
Community_Group_Train_sin_teb  <- Community_Group_Train_sin %>%
  select(-SampleID,  -community.name,  -community.number, -Amoxicillin, -Chlorothalonil,
         -Diflufenican, -Glyphosate, -Imidacloprid, -Metaldehyde, -Oxytetracycline,
         -chem.code, -complexity, -well, -plate,
         -set)
str(Community_Group_Train_sin_teb)
RFSRC_nocomsin_teb <- rfsrc(Tebuconazole ~ ., Community_Group_Train_sin_teb, ntree = 501, importance = TRUE, proximities = TRUE)
get.mv.error(RFSRC_nocomsin_teb)
Community_Group_Validations_sin_teb <- Community_Group_Validations_sin %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code,  -Amoxicillin, -Chlorothalonil, 
         -Diflufenican, -Glyphosate, -Imidacloprid, -Metaldehyde,
         -Oxytetracycline, -complexity, -well, -plate, -set)
str(Community_Group_Validations_sin_teb)
Predictions_nocomsin_teb <- predict(RFSRC_nocomsin_teb, newdata = Community_Group_Validations_sin_teb, type = "response")

roc_curve_nocomsin_teb <- roc(Community_Group_Validations_sin_teb$Tebuconazole, Predictions_nocomsin_teb$predicted[,2])


auc_value_nocomsin_teb <- auc(roc_curve_nocomsin_teb)
auc_value_nocomsin_teb

#CF
0.701





# Singe_4
Train_set_Com <- t_f %>%
  select(-SampleID,  -community.name, 
         -chem.code, -Amoxicillin, -Chlorothalonil, -Diflufenican, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate)


Train_set_Com$community.number <- as.factor(Train_set_Com$community.number)

RFSRC_Four_Com <- rfsrc(community.number ~ ., data = Train_set_Com, ntree = 501, importance = TRUE, proximities = TRUE)
Val_set_Com <- v_f %>%
  select(-SampleID,  -community.name, 
         -chem.code, -Amoxicillin, -Chlorothalonil, -Diflufenican, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate,-set)

Val_set_Com$community.number <- as.factor(Val_set_Com$community.number)

Predictions_Four_Com <- predict(RFSRC_Four_Com, newdata = Val_set_Com, type = "response")
Predictions_Four_Com
roc_curve_Com <- multiclass.roc(Val_set_Com$community.number, Predictions_Four_Com$predicted)


auc_value_Com <- auc(roc_curve_Com)
auc_value_Com
#ACC-Com
actual_classes_sin_com <- as.character(Val_set_Com$community.number)

predicted_classes_com_sin <- apply(Predictions_Four_Com$predicted, 1, which.max)
predicted_classes_com_sin  <- as.character(predicted_classes_com_sin )
adjusted_predicted_classes_com_sin <- as.character(ifelse(predicted_classes_com_sin  != "1", as.numeric(predicted_classes_com_sin ) + 2, as.character(predicted_classes_com_sin )))
accuracy_sin_com <- sum(adjusted_predicted_classes_com_sin  == actual_classes_sin_com) / length(adjusted_predicted_classes_com_sin)


#Com
err_rate_four_com <- RFSRC_Four_Com$err.rate

err_rate_four_com_df <- as.data.frame(err_rate_four_com)

model_ids_four_com <- dimnames(err_rate_four_com)[[2]]
colnames(err_rate_four_com_df) <- model_ids_four_com

err_rate_four_com_long <- err_rate_four_com_df %>%
  mutate(ntrees = 1:501) %>%
  pivot_longer(-ntrees, names_to = "model_id", values_to = "err_rate")

ggplot(err_rate_four_com_long, aes(x = ntrees, y = err_rate, color = model_id, group = model_id)) +
  geom_line() +
  geom_point()+
  geom_smooth(method = "loess", se = FALSE) +
  labs(
    title = "Community's Error Rate by Number of Trees and Model",
    x = "Number of Trees",
    y = "Error Rate",
    color = "Model"
  ) +
  theme_minimal()




Train_set_Amo <- t_f %>%
  select(-SampleID,  -community.name, -community.number, 
         -chem.code,  -Chlorothalonil, -Diflufenican, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate)


Train_set_Amo$Amoxicillin <- as.factor(Train_set_Amo$Amoxicillin)

RFSRC_Four_Amo <- rfsrc(Amoxicillin ~ ., data = Train_set_Amo, ntree = 501, importance = TRUE, proximities = TRUE)
Val_set_Amo <- v_f %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code,  -Chlorothalonil, -Diflufenican, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate)

Val_set_Amo$Amoxicillin <- as.factor(Val_set_Amo$Amoxicillin)

Predictions_Four_Amo <- predict(RFSRC_Four_Amo, newdata = Val_set_Amo, type = "response")
Predictions_Four_Amo
roc_curve_Amo <- roc(Val_set_Amo$Amoxicillin, Predictions_Four_Amo$predicted[,2])


auc_value_Amo <- auc(roc_curve_Amo)
auc_value_Amo
#CF
response_var_four_Amo <- c("Amoxicillin")
for (var in response_var_four_Amo) {
  actual_values <- as.factor(Val_set_Amo[[var]])
  predictions_raw <- Predictions_Four_Amo$predicted
  
  predictions <- ifelse(predictions_raw[, 1] > predictions_raw[, 2], 0, 1)
  predictions <- factor(predictions, levels = levels(actual_values))
}

confusion_matrix <- confusionMatrix(predictions, actual_values)


print(paste("Confusion Matrix for", var))
print(confusion_matrix)



#Chl
Train_set_Chl <- t_f %>%
  select(-SampleID,  -community.name, -community.number, 
         -chem.code, -Amoxicillin, -Diflufenican, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate)


Train_set_Chl$Chlorothalonil <- as.factor(Train_set_Chl$Chlorothalonil)

RFSRC_Four_Chl <- rfsrc(Chlorothalonil ~ ., data = Train_set_Chl, ntree = 501, importance = TRUE, proximities = TRUE)
Val_set_Chl <- v_f %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code, -Amoxicillin, -Diflufenican, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate)

Val_set_Chl$Chlorothalonil <- as.factor(Val_set_Chl$Chlorothalonil)

Predictions_Four_Chl <- predict(RFSRC_Four_Chl, newdata = Val_set_Chl, type = "response")
Predictions_Four_Chl
roc_curve_Chl <- roc(Val_set_Chl$Chlorothalonil, Predictions_Four_Chl$predicted[,2])


auc_value_Chl <- auc(roc_curve_Chl)
auc_value_Chl

#CF
response_var_four_Chl <- c("Chlorothalonil")
for (var in response_var_four_Chl) {
  actual_values <- as.factor(Val_set_Chl[[var]])
  predictions_raw <- Predictions_Four_Chl$predicted
  
  predictions <- ifelse(predictions_raw[, 1] > predictions_raw[, 2], 0, 1)
  predictions <- factor(predictions, levels = levels(actual_values))
}

confusion_matrix <- confusionMatrix(predictions, actual_values)

print(paste("Confusion Matrix for", var))
print(confusion_matrix)


#Dif
Train_set_Dif <- t_f %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code, -Amoxicillin, -Chlorothalonil, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate)


Train_set_Dif$Diflufenican <- as.factor(Train_set_Dif$Diflufenican)

RFSRC_Four_Dif <- rfsrc(Diflufenican ~ ., data = Train_set_Dif, ntree = 501, importance = TRUE, proximities = TRUE)
Val_set_Dif <- v_f %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code, -Amoxicillin, -Chlorothalonil, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate)

Val_set_Dif$Diflufenican <- as.factor(Val_set_Dif$Diflufenican)

Predictions_Four_Dif <- predict(RFSRC_Four_Dif, newdata = Val_set_Dif, type = "response")
Predictions_Four_Dif
roc_curve_Dif <- roc(Val_set_Dif$Diflufenican, Predictions_Four_Dif$predicted[,2])


auc_value_Dif <- auc(roc_curve_Dif)
auc_value_Dif


#CF
response_var_four_Dif <- c("Diflufenican")
for (var in response_var_four_Dif) {
  actual_values <- as.factor(Val_set_Dif[[var]])
  predictions_raw <- Predictions_Four_Dif$predicted
  
  predictions <- ifelse(predictions_raw[, 1] > predictions_raw[, 2], 0, 1)
  predictions <- factor(predictions, levels = levels(actual_values))
}

confusion_matrix <- confusionMatrix(predictions, actual_values)


print(paste("Confusion Matrix for", var))
print(confusion_matrix)

#Gly
Train_set_Gly <- t_f %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code, -Amoxicillin, -Chlorothalonil, -Diflufenican, 
         -Imidacloprid, -Metaldehyde, -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate)


Train_set_Gly$Glyphosate <- as.factor(Train_set_Gly$Glyphosate)

RFSRC_Four_Gly <- rfsrc(Glyphosate ~ ., data = Train_set_Gly, ntree = 501, importance = TRUE, proximities = TRUE)
Val_set_Gly <- v_f %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code, -Amoxicillin, -Chlorothalonil, -Diflufenican, 
         -Imidacloprid, -Metaldehyde, -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate)

Val_set_Gly$Glyphosate <- as.factor(Val_set_Gly$Glyphosate)

Predictions_Four_Gly <- predict(RFSRC_Four_Gly, newdata = Val_set_Gly, type = "response")
Predictions_Four_Gly
roc_curve_Gly <- roc(Val_set_Gly$Glyphosate, Predictions_Four_Gly$predicted[,2])


auc_value_Gly <- auc(roc_curve_Gly)
auc_value_Gly


#CF
response_var_four_Gly <- c("Glyphosate")
for (var in response_var_four_Gly) {
  actual_values <- as.factor(Val_set_Gly[[var]])
  predictions_raw <- Predictions_Four_Gly$predicted
  
  predictions <- ifelse(predictions_raw[, 1] > predictions_raw[, 2], 0, 1)
  predictions <- factor(predictions, levels = levels(actual_values))
}

confusion_matrix <- confusionMatrix(predictions, actual_values)


print(paste("Confusion Matrix for", var))
print(confusion_matrix)

#Imi
Train_set_Imi <- t_f %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code, -Amoxicillin, -Chlorothalonil, -Diflufenican, 
         -Glyphosate, -Metaldehyde, -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate)


Train_set_Imi$Imidacloprid <- as.factor(Train_set_Imi$Imidacloprid)

RFSRC_Four_Imi <- rfsrc(Imidacloprid ~ ., data = Train_set_Imi, ntree = 501, importance = TRUE, proximities = TRUE)
Val_set_Imi <- v_f %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code, -Amoxicillin, -Chlorothalonil, -Diflufenican, 
         -Glyphosate, -Metaldehyde, -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate)

Val_set_Imi$Imidacloprid <- as.factor(Val_set_Imi$Imidacloprid)

Predictions_Four_Imi <- predict(RFSRC_Four_Imi, newdata = Val_set_Imi, type = "response")
Predictions_Four_Imi
roc_curve_Imi <- roc(Val_set_Imi$Imidacloprid, Predictions_Four_Imi$predicted[,2])


auc_value_Imi <- auc(roc_curve_Imi)

#CF
response_var_four_Imi <- c("Imidacloprid")
for (var in response_var_four_Imi) {
  actual_values <- as.factor(Val_set_Imi[[var]])
  predictions_raw <- Predictions_Four_Imi$predicted
  
  predictions <- ifelse(predictions_raw[, 1] > predictions_raw[, 2], 0, 1)
  predictions <- factor(predictions, levels = levels(actual_values))
}

confusion_matrix <- confusionMatrix(predictions, actual_values)


print(paste("Confusion Matrix for", var))
print(confusion_matrix)



#Met
Train_set_Met <- t_f %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code, -Amoxicillin, -Chlorothalonil, -Diflufenican, 
         -Glyphosate, -Imidacloprid,  -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate)


Train_set_Met$Metaldehyde <- as.factor(Train_set_Met$Metaldehyde)

RFSRC_Four_Met <- rfsrc(Metaldehyde ~ ., data = Train_set_Met, ntree = 501, importance = TRUE, proximities = TRUE)
Val_set_Met <- v_f %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code, -Amoxicillin, -Chlorothalonil, -Diflufenican, 
         -Glyphosate, -Imidacloprid, -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate)

Val_set_Met$Metaldehyde <- as.factor(Val_set_Met$Metaldehyde)

Predictions_Four_Met <- predict(RFSRC_Four_Met, newdata = Val_set_Met, type = "response")
Predictions_Four_Met
roc_curve_Met <- roc(Val_set_Met$Metaldehyde, Predictions_Four_Met$predicted[,2])


auc_value_Met <- auc(roc_curve_Met)
auc_value_Met

#CF
response_var_four_Met <- c("Metaldehyde")
for (var in response_var_four_Met) {
  actual_values <- as.factor(Val_set_Met[[var]])
  predictions_raw <- Predictions_Four_Met$predicted
  
  predictions <- ifelse(predictions_raw[, 1] > predictions_raw[, 2], 0, 1)
  predictions <- factor(predictions, levels = levels(actual_values))
}

confusion_matrix <- confusionMatrix(predictions, actual_values)


print(paste("Confusion Matrix for", var))
print(confusion_matrix)


#Teb
Train_set_Teb <- t_f %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code, -Amoxicillin, -Chlorothalonil, -Diflufenican, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, -Oxytetracycline,
         -complexity, -well, -plate)


Train_set_Teb$Tebuconazole <- as.factor(Train_set_Teb$Tebuconazole)

RFSRC_Four_Teb <- rfsrc(Tebuconazole ~ ., data = Train_set_Teb, ntree = 501, importance = TRUE, proximities = TRUE)
Val_set_Teb <- v_f %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code, -Amoxicillin, -Chlorothalonil, -Diflufenican, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, -Oxytetracycline,
         -complexity, -well, -plate)

Val_set_Teb$Tebuconazole <- as.factor(Val_set_Teb$Tebuconazole)

Predictions_Four_Teb <- predict(RFSRC_Four_Teb, newdata = Val_set_Teb, type = "response")
Predictions_Four_Teb

roc_curve_Teb <- roc(Val_set_Teb$Tebuconazole, Predictions_Four_Teb$predicted[,2])

auc_value_Teb <- auc(roc_curve_Teb)

auc_value_Teb

#CF
response_var_four_Teb <- c("Tebuconazole")
for (var in response_var_four_Teb) {
  actual_values <- as.factor(Val_set_Teb[[var]])
  predictions_raw <- Predictions_Four_Teb$predicted
  
  predictions <- ifelse(predictions_raw[, 1] > predictions_raw[, 2], 0, 1)
  predictions <- factor(predictions, levels = levels(actual_values))
}
confusion_matrix <- confusionMatrix(predictions, actual_values)


print(paste("Confusion Matrix for", var))
print(confusion_matrix)


#Com and chem
#Com_Amo
Train_set_Amo_Com <- t_f %>%
  select(-SampleID,  -community.name, 
         -chem.code, -Chlorothalonil, -Diflufenican, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate)

Train_set_Amo_Com$Amoxicillin <- as.factor(Train_set_Amo_Com$Amoxicillin)
Train_set_Amo_Com$community.number <- as.factor(Train_set_Amo_Com$community.number)

RFSRC_Four_Amo_Res_Com <- rfsrc(Multivar(community.number, Amoxicillin) ~ ., data = Train_set_Amo_Com, ntree = 501, importance = TRUE, proximities = TRUE)
get.mv.error(RFSRC_Four_Amo_Res_Com, standardize = TRUE)

Val_set_Amo_Com <- v_f %>%
  select(-SampleID,  -community.name, 
         -chem.code, -Chlorothalonil, -Diflufenican, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate)
Val_set_Amo_Com$Amoxicillin <- as.factor(Val_set_Amo_Com$Amoxicillin)
Val_set_Amo_Com$community.number <- as.factor(Val_set_Amo_Com$community.number)

Predictions_Four_Amo_Res_Com <- predict(RFSRC_Four_Amo_Res_Com, newdata = Val_set_Amo_Com, type = "response")
Predictions_Four_Amo_Res_Com
#ROC and AUC-factor
auc_values_Four_Amo_Com <- list()
response_var_four_Amo_Com <- c("community.number",  "Amoxicillin")


for (var in response_var_four_Amo_Com) {
  actual_value_four_Amo_Com <- as.factor(Val_set_Amo_Com[[var]])
  
  if (length(levels(actual_value_four_Amo_Com )) == 2) {
    prediction_four_Amo_Com <- Predictions_Four_Amo_Res_Com$classOutput[[var]]$predicted[, 2]
    roc_curve <- roc(actual_value_four_Amo_Com, prediction_four_Amo_Com)
    auc_value <- auc(roc_curve)
    auc_values_Four_Amo_Com[[var]] <- auc_value
    print(paste("AUC for", var, "is", auc_value))
    
  } else {
    
    prediction_four_Amo_Com <- Predictions_Four_Amo_Res_Com$classOutput[[var]]$predicted
    roc_curve <- multiclass.roc(actual_value_four_Amo_Com, prediction_four_Amo_Com)
    auc_value <- auc(roc_curve)
    auc_values_Four_Amo_Com[[var]] <- auc_value
    print(paste("AUC for", var, "is", auc_value))
  } 
}
#CF
for (var in response_var_four_Amo_Com) {
  actual_values <- as.factor(Val_set_Amo_Com[[var]])
  predictions_raw <- Predictions_Four_Amo_Res_Com$classOutput[[var]]$predicted
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
#ACC-Com
actual_classes_amo_com_com <- as.character(Val_set_Amo_Com$community.number)

predicted_classes_com_amo_com <- apply(Predictions_Four_Amo_Res_Com$classOutput$community.number$predicted, 1, which.max)
predicted_classes_com_amo_com<- as.character(predicted_classes_com_amo_com)
adjusted_predicted_classes_com_amo_com <- as.character(ifelse(predicted_classes_com_amo_com != "1", as.numeric(predicted_classes_com_amo_com) + 2, as.character(predicted_classes_com_amo_com)))
accuracy_amo_com_com <- sum(adjusted_predicted_classes_com_amo_com  == actual_classes_amo_com_com) / length(adjusted_predicted_classes_com_amo_com)




#Com_Chl
Train_set_Chl_Com <- t_f %>%
  select(-SampleID,  -community.name, 
         -chem.code, -Amoxicillin,  -Diflufenican, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate)

Train_set_Chl_Com$Chlorothalonil <- as.factor(Train_set_Chl_Com$Chlorothalonil)
Train_set_Chl_Com$community.number <- as.factor(Train_set_Chl_Com$community.number)

RFSRC_Four_Chl_Res_Com <- rfsrc(Multivar(community.number, Chlorothalonil) ~ ., data = Train_set_Chl_Com, ntree = 501, importance = TRUE, proximities = TRUE)
get.mv.error(RFSRC_Four_Chl_Res_Com, standardize = TRUE)

Val_set_Chl_Com <- v_f %>%
  select(-SampleID,  -community.name, 
         -chem.code, -Amoxicillin, -Diflufenican, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate)
Val_set_Chl_Com$Chlorothalonil <- as.factor(Val_set_Chl_Com$Chlorothalonil)
Val_set_Chl_Com$community.number <- as.factor(Val_set_Chl_Com$community.number)

Predictions_Four_Chl_Res_Com <- predict(RFSRC_Four_Chl_Res_Com, newdata = Val_set_Chl_Com, type = "response")
Predictions_Four_Chl_Res_Com
#ROC and AUC-factor
auc_values_Four_Chl_Com <- list()
response_var_four_Chl_Com <- c("community.number",  "Chlorothalonil")


for (var in response_var_four_Chl_Com) {
  actual_value_four_Chl_Com <- as.factor(Val_set_Chl_Com[[var]])
  
  if (length(levels(actual_value_four_Chl_Com )) == 2) {
    prediction_four_Chl_Com <- Predictions_Four_Chl_Res_Com$classOutput[[var]]$predicted[, 2]
    roc_curve <- roc(actual_value_four_Chl_Com, prediction_four_Chl_Com)
    auc_value <- auc(roc_curve)
    auc_values_Four_Chl_Com[[var]] <- auc_value
    print(paste("AUC for", var, "is", auc_value))
    
  } else {
    
    prediction_four_Chl_Com <- Predictions_Four_Chl_Res_Com$classOutput[[var]]$predicted
    roc_curve <- multiclass.roc(actual_value_four_Chl_Com, prediction_four_Chl_Com)
    auc_value <- auc(roc_curve)
    auc_values_Four_Chl_Com[[var]] <- auc_value
    print(paste("AUC for", var, "is", auc_value))
  } 
}
#CF
for (var in response_var_four_Chl_Com) {
  actual_values <- as.factor(Val_set_Chl_Com[[var]])
  predictions_raw <- Predictions_Four_Chl_Res_Com$classOutput[[var]]$predicted
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
#ACC-Com
actual_classes_chl_com_com <- as.character(Val_set_Chl_Com$community.number)

predicted_classes_com_chl_com <- apply(Predictions_Four_Chl_Res_Com$classOutput$community.number$predicted, 1, which.max)
predicted_classes_com_chl_com <- as.character(predicted_classes_com_chl_com )
adjusted_predicted_classes_com_chl_com <- as.character(ifelse(predicted_classes_com_chl_com  != "1", as.numeric(predicted_classes_com_chl_com ) + 2, as.character(predicted_classes_com_chl_com )))
accuracy_chl_com_com <- sum(adjusted_predicted_classes_com_chl_com  == actual_classes_chl_com_com) / length(adjusted_predicted_classes_com_chl_com)



#Com_Dif
Train_set_Dif_Com <- t_f %>%
  select(-SampleID,  -community.name, 
         -chem.code, -Amoxicillin, -Chlorothalonil, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate)

Train_set_Dif_Com$Diflufenican <- as.factor(Train_set_Dif_Com$Diflufenican)
Train_set_Dif_Com$community.number <- as.factor(Train_set_Dif_Com$community.number)

RFSRC_Four_Dif_Res_Com <- rfsrc(Multivar(community.number, Diflufenican) ~ ., data = Train_set_Dif_Com, ntree = 501, importance = TRUE, proximities = TRUE)
get.mv.error(RFSRC_Four_Dif_Res_Com, standardize = TRUE)

Val_set_Dif_Com <- v_f %>%
  select(-SampleID,  -community.name, 
         -chem.code, -Amoxicillin, -Chlorothalonil, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate)
Val_set_Dif_Com$Diflufenican <- as.factor(Val_set_Dif_Com$Diflufenican)
Val_set_Dif_Com$community.number <- as.factor(Val_set_Dif_Com$community.number)

Predictions_Four_Dif_Res_Com <- predict(RFSRC_Four_Dif_Res_Com, newdata = Val_set_Dif_Com, type = "response")
Predictions_Four_Dif_Res_Com
#ROC and AUC-factor
auc_values_Four_Dif_Com <- list()
response_var_four_Dif_Com <- c("community.number",  "Diflufenican")


for (var in response_var_four_Dif_Com) {
  actual_value_four_Dif_Com <- as.factor(Val_set_Dif_Com[[var]])
  
  if (length(levels(actual_value_four_Dif_Com )) == 2) {
    prediction_four_Dif_Com <- Predictions_Four_Dif_Res_Com$classOutput[[var]]$predicted[, 2]
    roc_curve <- roc(actual_value_four_Dif_Com, prediction_four_Dif_Com)
    auc_value <- auc(roc_curve)
    auc_values_Four_Dif_Com[[var]] <- auc_value
    print(paste("AUC for", var, "is", auc_value))
    
  } else {
    
    prediction_four_Dif_Com <- Predictions_Four_Dif_Res_Com$classOutput[[var]]$predicted
    roc_curve <- multiclass.roc(actual_value_four_Dif_Com, prediction_four_Dif_Com)
    auc_value <- auc(roc_curve)
    auc_values_Four_Dif_Com[[var]] <- auc_value
    print(paste("AUC for", var, "is", auc_value))
  } 
}
#CF
for (var in response_var_four_Dif_Com) {
  actual_values <- as.factor(Val_set_Dif_Com[[var]])
  predictions_raw <- Predictions_Four_Dif_Res_Com$classOutput[[var]]$predicted
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
#ACC-Com
actual_classes_dif_com_com <- as.character(Val_set_Dif_Com$community.number)

predicted_classes_com_dif_com <- apply(Predictions_Four_Dif_Res_Com$classOutput$community.number$predicted, 1, which.max)
predicted_classes_com_dif_com<- as.character(predicted_classes_com_dif_com)
adjusted_predicted_classes_com_dif_com <- as.character(ifelse(predicted_classes_com_dif_com != "1", as.numeric(predicted_classes_com_dif_com) + 2, as.character(predicted_classes_com_dif_com)))
accuracy_dif_com_com <- sum(adjusted_predicted_classes_com_dif_com  == actual_classes_dif_com_com) / length(adjusted_predicted_classes_com_dif_com )



#Com_Gly
Train_set_Gly_Com <- t_f %>%
  select(-SampleID,  -community.name, 
         -chem.code, -Amoxicillin, -Chlorothalonil, -Diflufenican,
         -Imidacloprid, -Metaldehyde, -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate)

Train_set_Gly_Com$Glyphosate <- as.factor(Train_set_Gly_Com$Glyphosate)
Train_set_Gly_Com$community.number <- as.factor(Train_set_Gly_Com$community.number)

RFSRC_Four_Gly_Res_Com <- rfsrc(Multivar(community.number, Glyphosate) ~ ., data = Train_set_Gly_Com, ntree = 501, importance = TRUE, proximities = TRUE)
get.mv.error(RFSRC_Four_Gly_Res_Com, standardize = TRUE)

Val_set_Gly_Com <- v_f %>%
  select(-SampleID,  -community.name, 
         -chem.code, -Amoxicillin, -Chlorothalonil, -Diflufenican,
         -Imidacloprid, -Metaldehyde, -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate)
Val_set_Gly_Com$Glyphosate <- as.factor(Val_set_Gly_Com$Glyphosate)
Val_set_Gly_Com$community.number <- as.factor(Val_set_Gly_Com$community.number)

Predictions_Four_Gly_Res_Com <- predict(RFSRC_Four_Gly_Res_Com, newdata = Val_set_Gly_Com, type = "response")
Predictions_Four_Gly_Res_Com
#ROC and AUC-factor
auc_values_Four_Gly_Com <- list()
response_var_four_Gly_Com <- c("community.number",  "Glyphosate")


for (var in response_var_four_Gly_Com) {
  actual_value_four_Gly_Com <- as.factor(Val_set_Gly_Com[[var]])
  
  if (length(levels(actual_value_four_Gly_Com )) == 2) {
    prediction_four_Gly_Com <- Predictions_Four_Gly_Res_Com$classOutput[[var]]$predicted[, 2]
    roc_curve <- roc(actual_value_four_Gly_Com, prediction_four_Gly_Com)
    auc_value <- auc(roc_curve)
    auc_values_Four_Gly_Com[[var]] <- auc_value
    print(paste("AUC for", var, "is", auc_value))
    
  } else {
    
    prediction_four_Gly_Com <- Predictions_Four_Gly_Res_Com$classOutput[[var]]$predicted
    roc_curve <- multiclass.roc(actual_value_four_Gly_Com, prediction_four_Gly_Com)
    auc_value <- auc(roc_curve)
    auc_values_Four_Gly_Com[[var]] <- auc_value
    print(paste("AUC for", var, "is", auc_value))
  } 
}
#CF
for (var in response_var_four_Gly_Com) {
  actual_values <- as.factor(Val_set_Gly_Com[[var]])
  predictions_raw <- Predictions_Four_Gly_Res_Com$classOutput[[var]]$predicted
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
#ACC-Com
actual_classes_gly_com_com <- as.character(Val_set_Gly_Com$community.number)

predicted_classes_com_gly_com <- apply(Predictions_Four_Gly_Res_Com$classOutput$community.number$predicted, 1, which.max)
predicted_classes_com_gly_com<- as.character(predicted_classes_com_gly_com)
adjusted_predicted_classes_com_gly_com <- as.character(ifelse(predicted_classes_com_gly_com != "1", as.numeric(predicted_classes_com_gly_com) + 2, as.character(predicted_classes_com_gly_com)))
accuracy_gly_com_com <- sum(adjusted_predicted_classes_com_gly_com  == actual_classes_gly_com_com) / length(adjusted_predicted_classes_com_gly_com)



#Com_Imidacloprid
Train_set_Imi_Com <- t_f %>%
  select(-SampleID,  -community.name, 
         -chem.code, -Amoxicillin, -Chlorothalonil, -Diflufenican,
         -Glyphosate, -Metaldehyde, -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate)

Train_set_Imi_Com$Imidacloprid <- as.factor(Train_set_Imi_Com$Imidacloprid)
Train_set_Imi_Com$community.number <- as.factor(Train_set_Imi_Com$community.number)

RFSRC_Four_Imi_Res_Com <- rfsrc(Multivar(community.number, Imidacloprid) ~ ., data = Train_set_Imi_Com, ntree = 501, importance = TRUE, proximities = TRUE)
get.mv.error(RFSRC_Four_Imi_Res_Com, standardize = TRUE)

Val_set_Imi_Com <- v_f %>%
  select(-SampleID,  -community.name, 
         -chem.code, -Amoxicillin, -Chlorothalonil, -Diflufenican,
         -Glyphosate, -Metaldehyde, -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate)
Val_set_Imi_Com$Imidacloprid <- as.factor(Val_set_Imi_Com$Imidacloprid)
Val_set_Imi_Com$community.number <- as.factor(Val_set_Imi_Com$community.number)

Predictions_Four_Imi_Res_Com <- predict(RFSRC_Four_Imi_Res_Com, newdata = Val_set_Imi_Com, type = "response")
Predictions_Four_Imi_Res_Com
#ROC and AUC-factor
auc_values_Four_Imi_Com <- list()
response_var_four_Imi_Com <- c("community.number",  "Imidacloprid")


for (var in response_var_four_Imi_Com) {
  actual_value_four_Imi_Com <- as.factor(Val_set_Imi_Com[[var]])
  
  if (length(levels(actual_value_four_Imi_Com )) == 2) {
    prediction_four_Imi_Com <- Predictions_Four_Imi_Res_Com$classOutput[[var]]$predicted[, 2]
    roc_curve <- roc(actual_value_four_Imi_Com, prediction_four_Imi_Com)
    auc_value <- auc(roc_curve)
    auc_values_Four_Imi_Com[[var]] <- auc_value
    print(paste("AUC for", var, "is", auc_value))
    
  } else {
    
    prediction_four_Imi_Com <- Predictions_Four_Imi_Res_Com$classOutput[[var]]$predicted
    roc_curve <- multiclass.roc(actual_value_four_Imi_Com, prediction_four_Imi_Com)
    auc_value <- auc(roc_curve)
    auc_values_Four_Imi_Com[[var]] <- auc_value
    print(paste("AUC for", var, "is", auc_value))
  } 
}
#CF
for (var in response_var_four_Imi_Com) {
  actual_values <- as.factor(Val_set_Imi_Com[[var]])
  predictions_raw <- Predictions_Four_Imi_Res_Com$classOutput[[var]]$predicted
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
#ACC-Com
actual_classes_imi_com_com <- as.character(Val_set_Imi_Com$community.number)

predicted_classes_com_imi_com <- apply(Predictions_Four_Imi_Res_Com$classOutput$community.number$predicted, 1, which.max)
predicted_classes_com_imi_com<- as.character(predicted_classes_com_imi_com)
adjusted_predicted_classes_com_imi_com <- as.character(ifelse(predicted_classes_com_imi_com != "1", as.numeric(predicted_classes_com_imi_com) + 2, as.character(predicted_classes_com_imi_com)))
accuracy_imi_com_com <- sum(adjusted_predicted_classes_com_imi_com  == actual_classes_imi_com_com) / length(adjusted_predicted_classes_com_imi_com)

view(Val_set_Imi_Com)


#Com_Met
Train_set_Met_Com <- t_f %>%
  select(-SampleID,  -community.name, 
         -chem.code, -Amoxicillin, -Chlorothalonil, -Diflufenican,
         -Glyphosate, -Oxytetracycline, -Imidacloprid,
         -Tebuconazole, -complexity, -well, -plate)

Train_set_Met_Com$Metaldehyde <- as.factor(Train_set_Met_Com$Metaldehyde)
Train_set_Met_Com$community.number <- as.factor(Train_set_Met_Com$community.number)

RFSRC_Four_Met_Res_Com <- rfsrc(Multivar(community.number, Metaldehyde) ~ ., data = Train_set_Met_Com, ntree = 501, importance = TRUE, proximities = TRUE)
get.mv.error(RFSRC_Four_Met_Res_Com, standardize = TRUE)

Val_set_Met_Com <- v_f %>%
  select(-SampleID,  -community.name, 
         -chem.code, -Amoxicillin, -Chlorothalonil, -Diflufenican,
         -Glyphosate, -Imidacloprid, -Oxytetracycline,
         -Tebuconazole, -complexity, -well, -plate)
Val_set_Met_Com$Metaldehyde <- as.factor(Val_set_Met_Com$Metaldehyde)
Val_set_Met_Com$community.number <- as.factor(Val_set_Met_Com$community.number)

Predictions_Four_Met_Res_Com <- predict(RFSRC_Four_Met_Res_Com, newdata = Val_set_Met_Com, type = "response")
Predictions_Four_Met_Res_Com
#ROC and AUC-factor
auc_values_Four_Met_Com <- list()
response_var_four_Met_Com <- c("community.number",  "Metaldehyde")


for (var in response_var_four_Met_Com) {
  actual_value_four_Met_Com <- as.factor(Val_set_Met_Com[[var]])
  
  if (length(levels(actual_value_four_Met_Com )) == 2) {
    prediction_four_Met_Com <- Predictions_Four_Met_Res_Com$classOutput[[var]]$predicted[, 2]
    roc_curve <- roc(actual_value_four_Met_Com, prediction_four_Met_Com)
    auc_value <- auc(roc_curve)
    auc_values_Four_Met_Com[[var]] <- auc_value
    print(paste("AUC for", var, "is", auc_value))
    
  } else {
    
    prediction_four_Met_Com <- Predictions_Four_Met_Res_Com$classOutput[[var]]$predicted
    roc_curve <- multiclass.roc(actual_value_four_Met_Com, prediction_four_Met_Com)
    auc_value <- auc(roc_curve)
    auc_values_Four_Met_Com[[var]] <- auc_value
    print(paste("AUC for", var, "is", auc_value))
  } 
}
#CF
for (var in response_var_four_Met_Com) {
  actual_values <- as.factor(Val_set_Met_Com[[var]])
  predictions_raw <- Predictions_Four_Met_Res_Com$classOutput[[var]]$predicted
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
#ACC-Com
actual_classes_met_com_com <- Val_set_Met_Com$community.number

predicted_classes_com_met_com <- apply(Predictions_Four_Met_Res_Com$classOutput$community.number$predicted, 1, which.max)
actual_classes_met_com_com <- as.character(actual_classes_met_com_com)
predicted_classes_com_met_com <- as.character(predicted_classes_com_met_com)
adjusted_predicted_classes_com_met_com <- as.character(ifelse(predicted_classes_com_met_com != "1", as.numeric(predicted_classes_com_met_com) + 2, as.character(predicted_classes_com_met_com)))
accuracy_met_com_com <- sum(adjusted_predicted_classes_com_met_com  == actual_classes_met_com_com) / length(adjusted_predicted_classes_com_met_com)


#Com_Teb
Train_set_Teb_Com <- t_f %>%
  select(-SampleID,  -community.name, 
         -chem.code, -Amoxicillin, -Chlorothalonil, -Diflufenican,
         -Glyphosate, -Oxytetracycline, -Imidacloprid,
         -Metaldehyde, -complexity, -well, -plate)

Train_set_Teb_Com$Tebuconazole <- as.factor(Train_set_Teb_Com$Tebuconazole)
Train_set_Teb_Com$community.number <- as.factor(Train_set_Teb_Com$community.number)
str(Train_set_Teb_Com)
RFSRC_Four_Teb_Res_Com <- rfsrc(Multivar(community.number, Tebuconazole) ~ ., data = Train_set_Teb_Com, ntree = 501, importance = TRUE, proximities = TRUE)
get.mv.error(RFSRC_Four_Teb_Res_Com, standardize = TRUE)

Val_set_Teb_Com <- v_f %>%
  select(-SampleID,  -community.name, 
         -chem.code, -Amoxicillin, -Chlorothalonil, -Diflufenican,
         -Glyphosate, -Imidacloprid, -Oxytetracycline,
         -Metaldehyde, -complexity, -well, -plate)
Val_set_Teb_Com$Tebuconazole <- as.factor(Val_set_Teb_Com$Tebuconazole)
Val_set_Teb_Com$community.number <- as.factor(Val_set_Teb_Com$community.number)

Predictions_Four_Teb_Res_Com <- predict(RFSRC_Four_Teb_Res_Com, newdata = Val_set_Teb_Com, type = "response")
Predictions_Four_Teb_Res_Com
#ROC and AUC-factor
auc_values_Four_Teb_Com <- list()
response_var_four_Teb_Com <- c("community.number",  "Tebuconazole")


for (var in response_var_four_Teb_Com) {
  actual_value_four_Teb_Com <- as.factor(Val_set_Teb_Com[[var]])
  
  if (length(levels(actual_value_four_Teb_Com )) == 2) {
    prediction_four_Teb_Com <- Predictions_Four_Teb_Res_Com$classOutput[[var]]$predicted[, 2]
    roc_curve <- roc(actual_value_four_Teb_Com, prediction_four_Teb_Com)
    auc_value <- auc(roc_curve)
    auc_values_Four_Teb_Com[[var]] <- auc_value
    print(paste("AUC for", var, "is", auc_value))
    
  } else {
    
    prediction_four_Teb_Com <- Predictions_Four_Teb_Res_Com$classOutput[[var]]$predicted
    roc_curve <- multiclass.roc(actual_value_four_Teb_Com, prediction_four_Teb_Com)
    auc_value <- auc(roc_curve)
    auc_values_Four_Teb_Com[[var]] <- auc_value
    print(paste("AUC for", var, "is", auc_value))
  } 
}
#CF

for (var in response_var_four_Teb_Com) {
  actual_values <- as.factor(Val_set_Teb_Com[[var]])
  predictions_raw <- Predictions_Four_Teb_Res_Com$classOutput[[var]]$predicted
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
#ACC-Com
actual_classes_teb_com_com <- Val_set_Teb_Com$community.number

predicted_classes_com_teb_com <- apply(Predictions_Four_Teb_Res_Com$classOutput$community.number$predicted, 1, which.max)
actual_classes_teb_com_com<- as.character(actual_classes_teb_com_com)
predicted_classes_com_teb_com <- as.character(predicted_classes_com_teb_com)
adjusted_predicted_classes_com_teb_com<- as.character(ifelse(predicted_classes_com_teb_com != "1", as.numeric(predicted_classes_com_teb_com) + 2, as.character(predicted_classes_com_teb_com)))
accuracy_teb_com_com <- sum(adjusted_predicted_classes_com_teb_com  == actual_classes_teb_com_com) / length(adjusted_predicted_classes_com_met_com)







err_rate_Teb_Com_four_com <- RFSRC_Four_Teb_Res_Com$classOutput$community.number$err.rate

err_rate_Teb_Com_four_com_df <- as.data.frame(err_rate_Teb_Com_four_com)

model_ids_Teb_Com_four_com <- dimnames(err_rate_Teb_Com_four_com)[[2]]
colnames(err_rate_Teb_Com_four_com_df) <- model_ids_Teb_Com_four_com
nrow(err_rate_Teb_Com_four_com_df)

err_rate_Teb_Com_four_com_long <- err_rate_Teb_Com_four_com_df %>%
  mutate(ntrees = 1:501) %>%
  pivot_longer(-ntrees, names_to = "model_id", values_to = "err_rate")

ggplot(err_rate_Teb_Com_four_com_long, aes(x = ntrees, y = err_rate, color = model_id, group = model_id)) +
  geom_line() +
  geom_point()+
  geom_smooth(method = "loess", se = FALSE) +
  labs(
    title = "Community's Error Rate by Number of Trees and Model",
    x = "Number of Trees",
    y = "Error Rate",
    color = "Model"
  ) +
  theme_minimal()





#5 eachcommunity


Qfive_ASV_meta_no_NA <- removezerodata_meta_NA
Qfive_ASV_meta_no_NA$Amoxicillin <- as.factor(Qfive_ASV_meta_no_NA$Amoxicillin)
Qfive_ASV_meta_no_NA$Chlorothalonil <- as.factor(Qfive_ASV_meta_no_NA$Chlorothalonil)
Qfive_ASV_meta_no_NA$Diflufenican <- as.factor(Qfive_ASV_meta_no_NA$Diflufenican)
Qfive_ASV_meta_no_NA$Glyphosate <- as.factor(Qfive_ASV_meta_no_NA$Glyphosate)
Qfive_ASV_meta_no_NA$Imidacloprid <- as.factor(Qfive_ASV_meta_no_NA$Imidacloprid)
Qfive_ASV_meta_no_NA$Metaldehyde <- as.factor(Qfive_ASV_meta_no_NA$Metaldehyde)
Qfive_ASV_meta_no_NA$Oxytetracycline <- as.factor(Qfive_ASV_meta_no_NA$Oxytetracycline)
Qfive_ASV_meta_no_NA$Tebuconazole <- as.factor(Qfive_ASV_meta_no_NA$Tebuconazole)


Com_one <- Qfive_ASV_meta_no_NA [Qfive_ASV_meta_no_NA $community.number == 1, ]
Train_set_chem_one <- subset(Com_one, chem.code %in% c('CDGO', 'CM', 'A', 'IT', 'G'))

Train_set_chem_one <- Train_set_chem_one %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code, -complexity, -well, -plate)
str(Train_set_chem_one )
RFSRC_chem_one <- rfsrc(Multivar(Amoxicillin, Chlorothalonil, Diflufenican, 
                                           Glyphosate, Imidacloprid, Metaldehyde, Oxytetracycline, Tebuconazole) ~ .,
                                  data = Train_set_chem_one, ntree = 501, proximity = TRUE)


get.mv.error(RFSRC_chem_one, standardize = TRUE)

view(Val_set_Chem_one)
v_f_one <- subset(Com_one, chem.code %in% c('M', 'CGMT', 'AD'))


Val_set_Chem_one <- v_f_one  %>%
  select(-SampleID, -community.name, -community.number, -chem.code, -complexity,
         -well, -plate)

Predictions_Val_Chem_one <- predict(RFSRC_chem_one, newdata = Val_set_Chem_one, type = "response")
Predictions_Val_Chem_one


#Amo
roc_amo_com_one <- roc(Val_set_Chem_one$Amoxicillin, Predictions_Val_Chem_one$classOutput$Amoxicillin$predicted[,2])
AUC_amo_com_one <- auc(roc_amo_com_one)
AUC_amo_com_one
#Chl
roc_chl_com_one <- roc(Val_set_Chem_one$Chlorothalonil, Predictions_Val_Chem_one$classOutput$Chlorothalonil$predicted[,2])
AUC_chl_com_one <- auc(roc_chl_com_one)
AUC_chl_com_one
#Dif
roc_dif_com_one <- roc(Val_set_Chem_one$Diflufenican, Predictions_Val_Chem_one$classOutput$Diflufenican$predicted[,2])
AUC_dif_com_one <- auc(roc_dif_com_one)
AUC_dif_com_one
#Gly
roc_gly_com_one <- roc(Val_set_Chem_one$Glyphosate, Predictions_Val_Chem_one$classOutput$Glyphosate$predicted[,2])
AUC_gly_com_one <- auc(roc_gly_com_one)
AUC_gly_com_one
#Met
roc_met_com_one <- roc(Val_set_Chem_one$Metaldehyde, Predictions_Val_Chem_one$classOutput$Metaldehyde$predicted[,2])
AUC_met_com_one <- auc(roc_met_com_one)
AUC_met_com_one
#Teb
roc_teb_com_one <- roc(Val_set_Chem_one$Tebuconazole, Predictions_Val_Chem_one$classOutput$Tebuconazole$predicted[,2])
AUC_teb_com_one <- auc(roc_teb_com_one)
AUC_teb_com_one

#CF
response_com_one <- c ("Amoxicillin", "Chlorothalonil", "Diflufenican", "Glyphosate", "Metaldehyde", "Tebuconazole")
for (var in response_com_one) {
  actual_values <- as.factor(Val_set_Chem_one[[var]])
  predictions_raw <- Predictions_Val_Chem_one$classOutput[[var]]$predicted
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

#Com4-RF

Com_four <- Qfive_ASV_meta_no_NA [Qfive_ASV_meta_no_NA $community.number == 4, ]
Chem_Com_four <- Com_four %>%
  group_by(community.number, chem.code) %>%
  mutate(set = cur_group_id()) %>%
  ungroup() 
G_num_four <- max(unique(Chem_Com_four$set))
Train_data_Com_four <- ceiling(G_num_four * 0.7)
train_community_four <- sample(unique(Chem_Com_four$set), size = Train_data_Com_four, replace = FALSE)

Community_Chem_G_Set_four <- Chem_Com_four
Community_Chem_G_Set_four$group <- "validation"
Community_Chem_G_Set_four[Community_Chem_G_Set_four$set %in% train_community_four, ]$group <- "train"

Community_Chem_G_Train_four <- Community_Chem_G_Set_four[Community_Chem_G_Set_four$group == "train", ]
Community_Chem_G_Val_four <- Community_Chem_G_Set_four[Community_Chem_G_Set_four$group == "validation", ]
vrow_indices_four <- match(Community_Chem_G_Val_four$SampleID, Qfive_ASV_meta_no_NA$SampleID)
row_indices_four <- match(Community_Chem_G_Train_four$SampleID, Qfive_ASV_meta_no_NA$SampleID)
t_f_four <- Qfive_ASV_meta_no_NA[row_indices_four, ]
v_f_four <- Qfive_ASV_meta_no_NA[vrow_indices_four, ]

Train_set_chem_four <- t_f_four %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code, -complexity, -well, -plate)

view(Train_set_chem_four)

RFSRC_chem_four<- rfsrc(Multivar(Amoxicillin, Chlorothalonil, Diflufenican, 
                                 Glyphosate, Imidacloprid, Metaldehyde, Oxytetracycline, Tebuconazole) ~ .,
                        data = Train_set_chem_four, ntree = 501, proximity = TRUE)


get.mv.error(RFSRC_chem_four, standardize = TRUE)


Val_set_Chem_four <- v_f_four  %>%
  select(-SampleID, -community.name, -community.number, -chem.code, -complexity,
         -well, -plate)

Predictions_Val_Chem_four <- predict(RFSRC_chem_four, newdata = Val_set_Chem_four, type = "response")
Predictions_Val_Chem_four

auc_values_five_four <- list()

response_var_five_four <- c("Amoxicillin", "Chlorothalonil", "Diflufenican", "Glyphosate", "Imidacloprid", "Metaldehyde", "Oxytetracycline", "Tebuconazole")

for (var in response_var_five_four) {
  actual_value_five_four <- as.factor(Val_set_Chem_four[[var]])
  prediction_five_four <- Predictions_Val_Chem_four$classOutput[[var]]$predicted[,2]
  roc_curve <- roc(actual_value_five_four, prediction_five_four)
  auc_value <- auc(roc_curve)
  auc_values_five_four[[var]] <- auc_value
  print(paste("AUC for", var, "is", auc_value))
} 

#CF
for (var in response_var_five_four) {
  actual_values <- as.factor(Val_set_Chem_four[[var]])
  predictions_raw <- Predictions_Val_Chem_four$classOutput[[var]]$predicted
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


#Com5-RF

Com_five <- Qfive_ASV_meta_no_NA [Qfive_ASV_meta_no_NA $community.number == 5, ]
Chem_Com_five <- Com_five  %>%
  group_by(community.number, chem.code) %>%
  mutate(set = cur_group_id()) %>%
  ungroup() 
G_num_five <- max(unique(Chem_Com_five$set))
Train_data_Com_five <- ceiling(G_num_five * 0.7)

Community_Chem_G_Train_five <- Chem_Com_five[!Chem_Com_five$set %in% c(2, 5, 6), ]
Community_Chem_G_Val_five <- Chem_Com_five[Chem_Com_five$set %in% c(2, 5, 6), ]
vrow_indices_five <- match(Community_Chem_G_Val_five$SampleID, Qfive_ASV_meta_no_NA$SampleID)
row_indices_five <- match(Community_Chem_G_Train_five$SampleID, Qfive_ASV_meta_no_NA$SampleID)
t_f_five <- Qfive_ASV_meta_no_NA[row_indices_five, ]
v_f_five <- Qfive_ASV_meta_no_NA[vrow_indices_five, ]

Train_set_chem_five <- t_f_five %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code, -complexity, -well, -plate)

str(Train_set_chem_five)

RFSRC_chem_five <- rfsrc(Multivar(Amoxicillin, Chlorothalonil, Diflufenican, 
                                 Glyphosate, Imidacloprid, Metaldehyde, Oxytetracycline, Tebuconazole) ~ .,
                        data = Train_set_chem_five, ntree = 501, proximity = TRUE)


get.mv.error(RFSRC_chem_five, standardize = TRUE)


Val_set_Chem_five <- v_f_five  %>%
  select(-SampleID, -community.name, -community.number, -chem.code, -complexity,
         -well, -plate)


Predictions_Val_Chem_five <- predict(RFSRC_chem_five, newdata = Val_set_Chem_five, type = "response")
Predictions_Val_Chem_five

auc_values_five_five <- list()

response_var_five_five <- c("Amoxicillin", "Chlorothalonil", "Diflufenican", "Glyphosate", "Imidacloprid", "Metaldehyde", "Oxytetracycline", "Tebuconazole")

for (var in response_var_five_five) {
  actual_value_five_five <- as.factor(Val_set_Chem_five[[var]])
  prediction_five_five <- Predictions_Val_Chem_five$classOutput[[var]]$predicted[,2]
  roc_curve <- roc(actual_value_five_five, prediction_five_five, direction = "<")
  auc_value <- auc(roc_curve)
  auc_values_five_five[[var]] <- auc_value
  print(paste("AUC for", var, "is", auc_value))
} 
Predictions_Val_Chem_five$classOutput$Tebuconazole$predicted
Val_set_Chem_five$Tebuconazole
#CF
Val_set_Chem_five$Amoxicillin
Predictions_Val_Chem_five$classOutput$Amoxicillin$predicted
for (var in response_var_five_five) {
  actual_values <- as.factor(Val_set_Chem_five[[var]])
  predictions_raw <- Predictions_Val_Chem_five$classOutput[[var]]$predicted
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



#Com6-RF
Com_six <- Qfive_ASV_meta_no_NA [Qfive_ASV_meta_no_NA $community.number == 6, ]
Chem_Com_six <- Com_six %>%
  group_by(community.number, chem.code) %>%
  mutate(set = cur_group_id()) %>%
  ungroup() 
G_num_six <- max(unique(Chem_Com_six$set))
Train_data_Com_six <- ceiling(G_num_six * 0.7)
train_community_six <- sample(unique(Chem_Com_six$set), size = Train_data_Com_six, replace = FALSE)

Community_Chem_G_Set_six <- Chem_Com_six
Community_Chem_G_Set_six$group <- "validation"
Community_Chem_G_Set_six[Community_Chem_G_Set_six$set %in% train_community_six, ]$group <- "train"

Community_Chem_G_Train_six <- Community_Chem_G_Set_six[Community_Chem_G_Set_six$group == "train", ]
Community_Chem_G_Val_six<- Community_Chem_G_Set_six[Community_Chem_G_Set_six$group == "validation", ]
vrow_indices_six <- match(Community_Chem_G_Val_six$SampleID, Qfive_ASV_meta_no_NA$SampleID)
row_indices_six <- match(Community_Chem_G_Train_six$SampleID, Qfive_ASV_meta_no_NA$SampleID)
t_f_six <- Qfive_ASV_meta_no_NA[row_indices_six, ]
v_f_six <- Qfive_ASV_meta_no_NA[vrow_indices_six, ]


Train_set_chem_six <- t_f_six %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code, -complexity, -well, -plate)



str(Train_set_chem_six)

RFSRC_chem_six <- rfsrc(Multivar(Amoxicillin, Chlorothalonil, Diflufenican, 
                                 Glyphosate, Imidacloprid, Metaldehyde, Oxytetracycline, Tebuconazole) ~ .,
                        data = Train_set_chem_six, ntree = 501, proximity = TRUE)


get.mv.error(RFSRC_chem_six, standardize = TRUE)


Val_set_Chem_six <- v_f_six  %>%
  select(-SampleID, -community.name, -community.number, -chem.code, -complexity,
         -well, -plate)

Predictions_Val_Chem_six <- predict(RFSRC_chem_six, newdata = Val_set_Chem_six, type = "response")
Predictions_Val_Chem_six

auc_values_five_six <- list()

response_var_five_six <- c("Amoxicillin", "Chlorothalonil", "Diflufenican", "Glyphosate", "Imidacloprid", "Metaldehyde", "Oxytetracycline", "Tebuconazole")

for (var in response_var_five_six) {
  actual_value_five_six <- as.factor(Val_set_Chem_six[[var]])
  prediction_five_six <- Predictions_Val_Chem_six$classOutput[[var]]$predicted[,2]
  roc_curve <- roc(actual_value_five_six, prediction_five_six)
  auc_value <- auc(roc_curve)
  auc_values_five_six[[var]] <- auc_value
  print(paste("AUC for", var, "is", auc_value))
} 

#CF
for (var in response_var_five_six) {
  actual_values <- as.factor(Val_set_Chem_six[[var]])
  predictions_raw <- Predictions_Val_Chem_six$classOutput[[var]]$predicted
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

#Com7-RF
Com_seven <- Qfive_ASV_meta_no_NA [Qfive_ASV_meta_no_NA $community.number == 7, ]
Chem_Com_seven <- Com_seven %>%
  group_by(community.number, chem.code) %>%
  mutate(set = cur_group_id()) %>%
  ungroup() 
G_num_seven <- max(unique(Chem_Com_seven$set))
Train_data_Com_seven <- ceiling(G_num_seven * 0.7)
train_community_seven <- sample(unique(Chem_Com_seven$set), size = Train_data_Com_seven, replace = FALSE)

Community_Chem_G_Set_seven <- Chem_Com_seven
Community_Chem_G_Set_seven$group <- "validation"
Community_Chem_G_Set_seven[Community_Chem_G_Set_seven$set %in% train_community_seven, ]$group <- "train"

Community_Chem_G_Train_seven <- Community_Chem_G_Set_seven[Community_Chem_G_Set_seven$group == "train", ]
Community_Chem_G_Val_seven<- Community_Chem_G_Set_seven[Community_Chem_G_Set_seven$group == "validation", ]
vrow_indices_seven <- match(Community_Chem_G_Val_seven$SampleID, Qfive_ASV_meta_no_NA$SampleID)
row_indices_seven <- match(Community_Chem_G_Train_seven$SampleID, Qfive_ASV_meta_no_NA$SampleID)
t_f_seven <- Qfive_ASV_meta_no_NA[row_indices_seven, ]
v_f_seven <- Qfive_ASV_meta_no_NA[vrow_indices_seven, ]


Train_set_chem_seven <- t_f_seven %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code, -complexity, -well, -plate)


str(Train_set_chem_seven)

RFSRC_chem_seven<- rfsrc(Multivar(Amoxicillin, Chlorothalonil, Diflufenican, 
                                 Glyphosate, Imidacloprid, Metaldehyde, Oxytetracycline, Tebuconazole) ~ .,
                        data = Train_set_chem_seven, ntree = 501, proximity = TRUE)


get.mv.error(RFSRC_chem_seven, standardize = TRUE)



Val_set_Chem_seven <- v_f_seven  %>%
  select(-SampleID, -community.name, -community.number, -chem.code, -complexity,
         -well, -plate)

Predictions_Val_Chem_seven <- predict(RFSRC_chem_seven, newdata = Val_set_Chem_seven, type = "response")
Predictions_Val_Chem_seven

auc_values_five_seven <- list()

response_var_five_seven <- c("Amoxicillin", "Chlorothalonil", "Diflufenican", "Glyphosate", "Imidacloprid", "Metaldehyde", "Oxytetracycline", "Tebuconazole")

for (var in response_var_five_seven) {
  actual_value_five_seven <- as.factor(Val_set_Chem_seven[[var]])
  prediction_five_seven <- Predictions_Val_Chem_seven$classOutput[[var]]$predicted[,2]
  roc_curve <- roc(actual_value_five_seven, prediction_five_seven)
  auc_value <- auc(roc_curve)
  auc_values_five_seven[[var]] <- auc_value
  print(paste("AUC for", var, "is", auc_value))
} 
view(Val_set_Chem_seven)

#CF
for (var in response_var_five_seven) {
  actual_values <- as.factor(Val_set_Chem_seven[[var]])
  predictions_raw <- Predictions_Val_Chem_seven$classOutput[[var]]$predicted
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


#Com8-RF
Com_eight <- Qfive_ASV_meta_no_NA [Qfive_ASV_meta_no_NA $community.number == 8, ]
Chem_Com_eight <- Com_eight %>%
  group_by(community.number, chem.code) %>%
  mutate(set = cur_group_id()) %>%
  ungroup() 
G_num_eight <- max(unique(Chem_Com_eight$set))
Train_data_Com_eight <- ceiling(G_num_eight * 0.7)
train_community_eight <- sample(unique(Chem_Com_eight$set), size = Train_data_Com_eight, replace = FALSE)

Community_Chem_G_Set_eight <- Chem_Com_eight
Community_Chem_G_Set_eight$group <- "validation"
Community_Chem_G_Set_eight[Community_Chem_G_Set_eight$set %in% train_community_eight, ]$group <- "train"

Community_Chem_G_Train_eight <- Community_Chem_G_Set_eight[Community_Chem_G_Set_eight$group == "train", ]
Community_Chem_G_Val_eight<- Community_Chem_G_Set_eight[Community_Chem_G_Set_eight$group == "validation", ]
vrow_indices_eight <- match(Community_Chem_G_Val_eight$SampleID, Qfive_ASV_meta_no_NA$SampleID)
row_indices_eight <- match(Community_Chem_G_Train_eight$SampleID, Qfive_ASV_meta_no_NA$SampleID)
t_f_eight<- Qfive_ASV_meta_no_NA[row_indices_eight, ]
v_f_eight <- Qfive_ASV_meta_no_NA[vrow_indices_eight, ]


Train_set_chem_eight <- t_f_eight %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code, -complexity, -well, -plate)


str(Train_set_chem_eight)

RFSRC_chem_eight <- rfsrc(Multivar(Amoxicillin, Chlorothalonil, Diflufenican, 
                                 Glyphosate, Imidacloprid, Metaldehyde, Oxytetracycline, Tebuconazole) ~ .,
                        data = Train_set_chem_eight, ntree = 501, proximity = TRUE)


get.mv.error(RFSRC_chem_eight, standardize = TRUE)


Val_set_Chem_eight <- v_f_eight %>%
  select(-SampleID, -community.name, -community.number, -chem.code, -complexity,
         -well, -plate)

Predictions_Val_Chem_eight <- predict(RFSRC_chem_eight, newdata = Val_set_Chem_eight, type = "response")
Predictions_Val_Chem_eight

auc_values_five_eight <- list()

response_var_five_eight <- c("Amoxicillin", "Chlorothalonil", "Diflufenican", "Glyphosate", "Imidacloprid", "Metaldehyde", "Oxytetracycline", "Tebuconazole")

for (var in response_var_five_eight) {
  actual_value_five_eight <- as.factor(Val_set_Chem_eight[[var]])
  prediction_five_eight <- Predictions_Val_Chem_eight$classOutput[[var]]$predicted[,2]
  roc_curve <- roc(actual_value_five_eight, prediction_five_eight, direction = "<")
  auc_value <- auc(roc_curve)
  auc_values_five_eight[[var]] <- auc_value
  print(paste("AUC for", var, "is", auc_value))
} 

#CF
for (var in response_var_five_eight) {
  actual_values <- as.factor(Val_set_Chem_eight[[var]])
  predictions_raw <- Predictions_Val_Chem_eight$classOutput[[var]]$predicted
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


#Com9-RF
Com_nine <- Qfive_ASV_meta_no_NA [Qfive_ASV_meta_no_NA $community.number == 9, ]
Chem_Com_nine <- Com_nine %>%
  group_by(community.number, chem.code) %>%
  mutate(set = cur_group_id()) %>%
  ungroup() 
G_num_nine <- max(unique(Chem_Com_nine$set))
Train_data_Com_nine <- ceiling(G_num_nine * 0.7)
train_community_nine <- sample(unique(Chem_Com_nine$set), size = Train_data_Com_nine, replace = FALSE)

Community_Chem_G_Set_nine <- Chem_Com_nine
Community_Chem_G_Set_nine$group <- "validation"
Community_Chem_G_Set_nine[Community_Chem_G_Set_nine$set %in% train_community_nine, ]$group <- "train"

Community_Chem_G_Train_nine <- Community_Chem_G_Set_nine[Community_Chem_G_Set_nine$group == "train", ]
Community_Chem_G_Val_nine<- Community_Chem_G_Set_nine[Community_Chem_G_Set_nine$group == "validation", ]
vrow_indices_nine <- match(Community_Chem_G_Val_nine$SampleID, Qfive_ASV_meta_no_NA$SampleID)
row_indices_nine <- match(Community_Chem_G_Train_nine$SampleID, Qfive_ASV_meta_no_NA$SampleID)
t_f_nine<- Qfive_ASV_meta_no_NA[row_indices_nine, ]
v_f_nine <- Qfive_ASV_meta_no_NA[vrow_indices_nine, ]


Train_set_chem_nine <- t_f_nine %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code, -complexity, -well, -plate)


str(Train_set_chem_nine)

RFSRC_chem_nine <- rfsrc(Multivar(Amoxicillin, Chlorothalonil, Diflufenican, 
                                 Glyphosate, Imidacloprid, Metaldehyde, Oxytetracycline, Tebuconazole) ~ .,
                        data = Train_set_chem_nine, ntree = 501, proximity = TRUE)


get.mv.error(RFSRC_chem_nine, standardize = TRUE)

Val_set_Chem_nine <- v_f_nine  %>%
  select(-SampleID, -community.name, -community.number, -chem.code, -complexity,
         -well, -plate)

Predictions_Val_Chem_nine <- predict(RFSRC_chem_nine, newdata = Val_set_Chem_nine, type = "response")
Predictions_Val_Chem_nine

auc_values_five_nine <- list()

response_var_five_nine <- c("Amoxicillin", "Chlorothalonil", "Diflufenican", "Glyphosate", "Imidacloprid", "Metaldehyde", "Oxytetracycline", "Tebuconazole")

for (var in response_var_five_nine) {
  actual_value_five_nine <- as.factor(Val_set_Chem_nine[[var]])
  prediction_five_nine <- Predictions_Val_Chem_nine$classOutput[[var]]$predicted[,2]
  roc_curve <- roc(actual_value_five_nine, prediction_five_nine)
  auc_value <- auc(roc_curve)
  auc_values_five_nine[[var]] <- auc_value
  print(paste("AUC for", var, "is", auc_value))
} 

#CF
for (var in response_var_five_nine) {
  actual_values <- as.factor(Val_set_Chem_nine[[var]])
  predictions_raw <- Predictions_Val_Chem_nine$classOutput[[var]]$predicted
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

#Com10-RF
Com_ten <- Qfive_ASV_meta_no_NA [Qfive_ASV_meta_no_NA $community.number == 10, ]
Chem_Com_ten <- Com_ten %>%
  group_by(community.number, chem.code) %>%
  mutate(set = cur_group_id()) %>%
  ungroup() 
G_num_ten <- max(unique(Chem_Com_ten$set))
Train_data_Com_ten <- ceiling(G_num_ten * 0.7)
train_community_ten <- sample(unique(Chem_Com_ten$set), size = Train_data_Com_ten, replace = FALSE)

Community_Chem_G_Set_ten <- Chem_Com_ten
Community_Chem_G_Set_ten$group <- "validation"
Community_Chem_G_Set_ten[Community_Chem_G_Set_ten$set %in% train_community_ten, ]$group <- "train"

Community_Chem_G_Train_ten <- Community_Chem_G_Set_ten[Community_Chem_G_Set_ten$group == "train", ]
Community_Chem_G_Val_ten<- Community_Chem_G_Set_ten[Community_Chem_G_Set_ten$group == "validation", ]
vrow_indices_ten <- match(Community_Chem_G_Val_ten$SampleID, Qfive_ASV_meta_no_NA$SampleID)
row_indices_ten <- match(Community_Chem_G_Train_ten$SampleID, Qfive_ASV_meta_no_NA$SampleID)
t_f_ten<- Qfive_ASV_meta_no_NA[row_indices_ten, ]
v_f_ten <- Qfive_ASV_meta_no_NA[vrow_indices_ten, ]


Train_set_chem_ten <- t_f_ten %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code, -complexity, -well, -plate)


str(Train_set_chem_ten)

RFSRC_chem_ten <- rfsrc(Multivar(Amoxicillin, Chlorothalonil, Diflufenican, 
                                 Glyphosate, Imidacloprid, Metaldehyde, Oxytetracycline, Tebuconazole) ~ .,
                        data = Train_set_chem_ten, ntree = 501, proximity = TRUE)


get.mv.error(RFSRC_chem_ten, standardize = TRUE)


Val_set_Chem_ten <- v_f_ten  %>%
  select(-SampleID, -community.name, -community.number, -chem.code, -complexity,
         -well, -plate)

Predictions_Val_Chem_ten <- predict(RFSRC_chem_ten, newdata = Val_set_Chem_ten, type = "response")
Predictions_Val_Chem_ten

auc_values_five_ten <- list()

response_var_five_ten <- c("Amoxicillin", "Chlorothalonil", "Diflufenican", "Glyphosate", "Imidacloprid", "Metaldehyde", "Oxytetracycline", "Tebuconazole")

for (var in response_var_five_ten) {
  actual_value_five_ten <- as.factor(Val_set_Chem_ten[[var]])
  prediction_five_ten <- Predictions_Val_Chem_ten$classOutput[[var]]$predicted[,2]
  roc_curve <- roc(actual_value_five_ten, prediction_five_ten)
  auc_value <- auc(roc_curve)
  auc_values_five_ten[[var]] <- auc_value
  print(paste("AUC for", var, "is", auc_value))
} 

#CF
for (var in response_var_five_ten) {
  actual_values <- as.factor(Val_set_Chem_ten[[var]])
  predictions_raw <- Predictions_Val_Chem_ten$classOutput[[var]]$predicted
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


#Com11-RF
Com_eleven <- Qfive_ASV_meta_no_NA [Qfive_ASV_meta_no_NA $community.number == 11, ]
Chem_Com_eleven <- Com_eleven %>%
  group_by(community.number, chem.code) %>%
  mutate(set = cur_group_id()) %>%
  ungroup() 
G_num_eleven <- max(unique(Chem_Com_eleven$set))
Train_data_Com_eleven <- ceiling(G_num_eleven * 0.7)
train_community_eleven <- sample(unique(Chem_Com_eleven$set), size = Train_data_Com_eleven, replace = FALSE)

Community_Chem_G_Set_eleven <- Chem_Com_eleven
Community_Chem_G_Set_eleven$group <- "validation"
Community_Chem_G_Set_eleven[Community_Chem_G_Set_eleven$set %in% train_community_eleven, ]$group <- "train"

Community_Chem_G_Train_eleven <- Community_Chem_G_Set_eleven[Community_Chem_G_Set_eleven$group == "train", ]
Community_Chem_G_Val_eleven<- Community_Chem_G_Set_eleven[Community_Chem_G_Set_eleven$group == "validation", ]
vrow_indices_eleven <- match(Community_Chem_G_Val_eleven$SampleID, Qfive_ASV_meta_no_NA$SampleID)
row_indices_eleven <- match(Community_Chem_G_Train_eleven$SampleID, Qfive_ASV_meta_no_NA$SampleID)
t_f_eleven<- Qfive_ASV_meta_no_NA[row_indices_eleven, ]
v_f_eleven <- Qfive_ASV_meta_no_NA[vrow_indices_eleven, ]


Train_set_chem_eleven <- t_f_eleven %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code, -complexity, -well, -plate)


str(Train_set_chem_eleven)

RFSRC_chem_eleven <- rfsrc(Multivar(Amoxicillin, Chlorothalonil, Diflufenican, 
                                 Glyphosate, Imidacloprid, Metaldehyde, Oxytetracycline, Tebuconazole) ~ .,
                        data = Train_set_chem_eleven, ntree = 501, proximity = TRUE)


get.mv.error(RFSRC_chem_eleven, standardize = TRUE)


Val_set_Chem_eleven <- v_f_eleven  %>%
  select(-SampleID, -community.name, -community.number, -chem.code, -complexity,
         -well, -plate)

Predictions_Val_Chem_eleven <- predict(RFSRC_chem_eleven, newdata = Val_set_Chem_eleven, type = "response")
Predictions_Val_Chem_eleven

auc_values_five_eleven <- list()

response_var_five_eleven <- c("Amoxicillin", "Chlorothalonil", "Diflufenican", "Glyphosate", "Imidacloprid", "Metaldehyde", "Oxytetracycline", "Tebuconazole")

for (var in response_var_five_eleven) {
  actual_value_five_eleven <- as.factor(Val_set_Chem_eleven[[var]])
  prediction_five_eleven <- Predictions_Val_Chem_eleven$classOutput[[var]]$predicted[,2]
  roc_curve <- roc(actual_value_five_eleven, prediction_five_eleven, diection = "<")
  auc_value <- auc(roc_curve)
  auc_values_five_eleven[[var]] <- auc_value
  print(paste("AUC for", var, "is", auc_value))
} 

#CF
for (var in response_var_five_eleven) {
  actual_values <- as.factor(Val_set_Chem_eleven[[var]])
  predictions_raw <- Predictions_Val_Chem_eleven$classOutput[[var]]$predicted
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


#Com12-RF
Com_twelve <- Qfive_ASV_meta_no_NA [Qfive_ASV_meta_no_NA $community.number == 12, ]
Chem_Com_twelve <- Com_twelve %>%
  group_by(community.number, chem.code) %>%
  mutate(set = cur_group_id()) %>%
  ungroup() 
G_num_twelve <- max(unique(Chem_Com_twelve$set))
Train_data_Com_twelve <- ceiling(G_num_twelve * 0.7)
train_community_twelve <- sample(unique(Chem_Com_twelve$set), size = Train_data_Com_twelve, replace = FALSE)

Community_Chem_G_Set_twelve <- Chem_Com_twelve
Community_Chem_G_Set_twelve$group <- "validation"
Community_Chem_G_Set_twelve[Community_Chem_G_Set_twelve$set %in% train_community_twelve, ]$group <- "train"

Community_Chem_G_Train_twelve <- Community_Chem_G_Set_twelve[Community_Chem_G_Set_twelve$group == "train", ]
Community_Chem_G_Val_twelve<- Community_Chem_G_Set_twelve[Community_Chem_G_Set_twelve$group == "validation", ]
vrow_indices_twelve <- match(Community_Chem_G_Val_twelve$SampleID, Qfive_ASV_meta_no_NA$SampleID)
row_indices_twelve <- match(Community_Chem_G_Train_twelve$SampleID, Qfive_ASV_meta_no_NA$SampleID)
t_f_twelve<- Qfive_ASV_meta_no_NA[row_indices_twelve, ]
v_f_twelve <- Qfive_ASV_meta_no_NA[vrow_indices_twelve, ]

Train_set_chem_twelve  <- t_f_twelve  %>%
  select(-SampleID,  -community.name, -community.number,
         -chem.code, -complexity, -well, -plate)


str(Train_set_chem_twelve)

RFSRC_chem_twelve <- rfsrc(Multivar(Amoxicillin, Chlorothalonil, Diflufenican, 
                                 Glyphosate, Imidacloprid, Metaldehyde, Oxytetracycline, Tebuconazole) ~ .,
                        data = Train_set_chem_twelve, ntree = 501, proximity = TRUE)


get.mv.error(RFSRC_chem_twelve, standardize = TRUE)


Val_set_Chem_twelve <- v_f_twelve  %>%
  select(-SampleID, -community.name, -community.number, -chem.code, -complexity,
         -well, -plate)

Predictions_Val_Chem_twelve <- predict(RFSRC_chem_twelve, newdata = Val_set_Chem_twelve, type = "response")
Predictions_Val_Chem_twelve

auc_values_five_twelve <- list()

response_var_five_twelve <- c("Amoxicillin", "Chlorothalonil", "Diflufenican", "Glyphosate", "Imidacloprid", "Metaldehyde", "Oxytetracycline", "Tebuconazole")

for (var in response_var_five_twelve) {
  actual_value_five_twelve <- as.factor(Val_set_Chem_twelve[[var]])
  prediction_five_twelve <- Predictions_Val_Chem_twelve$classOutput[[var]]$predicted[,2]
  roc_curve <- roc(actual_value_five_twelve, prediction_five_twelve, diection = "<")
  auc_value <- auc(roc_curve)
  auc_values_five_twelve[[var]] <- auc_value
  print(paste("AUC for", var, "is", auc_value))
} 
Predictions_Val_Chem_twelve$classOutput$Tebuconazole$predicted
Val_set_Chem_twelve$Tebuconazole
?roc
#CF
Val_set_Chem_twelve$Chlorothalonil
Predictions_Val_Chem_twelve$classOutput$Chlorothalonil$predicted
for (var in response_var_five_twelve) {
  actual_values <- as.factor(Val_set_Chem_twelve[[var]])
  predictions_raw <- Predictions_Val_Chem_twelve$classOutput[[var]]$predicted
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






#No-Oxy
removezerodata_meta_NA$community.number <- as.factor(removezerodata_meta_NA$community.number)
removezerodata_meta_NA$chem.code <- as.factor(removezerodata_meta_NA$chem.code)
removezerodata_meta_NA_no_oxy <- removezerodata_meta_NA[removezerodata_meta_NA$Oxytetracycline == 0, ]
Chem_Com_group_no_oxy <- removezerodata_meta_NA_no_oxy %>%
  group_by(community.number, chem.code) %>%
  mutate(set = cur_group_id()) %>%
  ungroup() 
G_num_no_oxy <- max(unique(Chem_Com_group_no_oxy$set))
Train_data_no_oxy <- ceiling(G_num_no_oxy * 0.7)
test_community_g_no_oxy <- sample(unique(Chem_Com_group_no_oxy$set), size = Train_data_no_oxy, replace = FALSE)

Community_Chem_G_Set_no_oxy <- Chem_Com_group_no_oxy
Community_Chem_G_Set_no_oxy$group <- "validation"

Community_Chem_G_Set_no_oxy[Community_Chem_G_Set_no_oxy$set %in% test_community_g_no_oxy, ]$group <- "test"
Community_Chem_G_Test_no_oxy <- Community_Chem_G_Set_no_oxy[Community_Chem_G_Set_no_oxy$group == "test", ]
Community_Chem_G_Val_no_oxy <- Community_Chem_G_Set_no_oxy[Community_Chem_G_Set_no_oxy$group == "validation", ]


val_sample_ids_no_oxy <- Community_Chem_G_Val_no_oxy$SampleID
meta_sample_ids_no_oxy <- removezerodata_meta_NA_no_oxy$SampleID


common_sample_ids_no_oxy_val <- intersect(val_sample_ids_no_oxy, meta_sample_ids_no_oxy)


v_f_no_oxy <- removezerodata_meta_NA_no_oxy[removezerodata_meta_NA_no_oxy$SampleID %in% common_sample_ids_no_oxy, ]

test_sample_ids_no_oxy <- Community_Chem_G_Test_no_oxy$SampleID


common_sample_ids_no_oxy_test <- intersect(test_sample_ids_no_oxy, meta_sample_ids_no_oxy)


t_f_no_oxy <- removezerodata_meta_NA_no_oxy[removezerodata_meta_NA_no_oxy$SampleID %in% common_sample_ids_no_oxy_test, ]


Train_set_chem_no_oxy_com  <- t_f_no_oxy  %>%
  select(-SampleID,  -community.name, 
         -chem.code, -Oxytetracycline, -complexity, -well, -plate)

Train_set_chem_no_oxy_com$community.number <- as.factor(Train_set_chem_no_oxy_com$community.number)
Train_set_chem_no_oxy_com$Amoxicillin <- as.factor(Train_set_chem_no_oxy_com$Amoxicillin)
Train_set_chem_no_oxy_com$Chlorothalonil <- as.factor(Train_set_chem_no_oxy_com$Chlorothalonil)
Train_set_chem_no_oxy_com$Diflufenican <- as.factor(Train_set_chem_no_oxy_com$Diflufenican)
Train_set_chem_no_oxy_com$Glyphosate <- as.factor(Train_set_chem_no_oxy_com$Glyphosate)
Train_set_chem_no_oxy_com$Imidacloprid <- as.factor(Train_set_chem_no_oxy_com$Imidacloprid)
Train_set_chem_no_oxy_com$Metaldehyde <- as.factor(Train_set_chem_no_oxy_com$Metaldehyde)
Train_set_chem_no_oxy_com$Tebuconazole <- as.factor(Train_set_chem_no_oxy_com$Tebuconazole)

str(Train_set_chem_no_oxy_com)
str(Val_set_chem_no_oxy_com)
RFSRC_chem_no_oxy_com <- rfsrc(Multivar(community.number, Amoxicillin, Chlorothalonil, Diflufenican, 
                                    Glyphosate, Imidacloprid, Metaldehyde, Tebuconazole) ~ .,
                           data = Train_set_chem_no_oxy_com, ntree = 501, proximity = TRUE)

RFSRC_chem_no_oxy_com_new <- rfsrc(Multivar(community.number, Amoxicillin, Chlorothalonil, Diflufenican, 
                                        Glyphosate, Imidacloprid, Metaldehyde, Tebuconazole) ~ .,
                               data = Train_set_chem_no_oxy_com, ntree = 200, proximity = TRUE)

get.mv.error(RFSRC_chem_no_oxy_com_new, standardize = TRUE)


Val_set_chem_no_oxy_com <- v_f_no_oxy  %>%
  select(-SampleID, -community.name,  -chem.code, -complexity, -Oxytetracycline,
         -well, -plate)
Val_set_chem_no_oxy_com$community.number <- as.factor(Val_set_chem_no_oxy_com$community.number)
Val_set_chem_no_oxy_com$Amoxicillin <- as.factor(Val_set_chem_no_oxy_com$Amoxicillin)
Val_set_chem_no_oxy_com$Chlorothalonil <- as.factor(Val_set_chem_no_oxy_com$Chlorothalonil)
Val_set_chem_no_oxy_com$Diflufenican <- as.factor(Val_set_chem_no_oxy_com$Diflufenican)
Val_set_chem_no_oxy_com$Glyphosate <- as.factor(Val_set_chem_no_oxy_com$Glyphosate)
Val_set_chem_no_oxy_com$Imidacloprid <- as.factor(Val_set_chem_no_oxy_com$Imidacloprid)
Val_set_chem_no_oxy_com$Metaldehyde <- as.factor(Val_set_chem_no_oxy_com$Metaldehyde)
Val_set_chem_no_oxy_com$Tebuconazole <- as.factor(Val_set_chem_no_oxy_com$Tebuconazole)

Predictions_Val_Chem_no_oxy_com <- predict(RFSRC_chem_no_oxy_com, newdata = Val_set_chem_no_oxy_com, type = "response")
Predictions_Val_Chem_no_oxy_com
Predictions_Val_Chem_no_oxy_com_new <- predict(RFSRC_chem_no_oxy_com_new, newdata = Val_set_chem_no_oxy_com, type = "response")

auc_values_five_Chem_no_oxy_com <- list()

response_var_five_Chem_no_oxy_com <- c("community.number", "Amoxicillin", "Chlorothalonil", "Diflufenican", "Glyphosate", "Imidacloprid", "Metaldehyde", "Tebuconazole")

for (var in response_var_five_Chem_no_oxy_com) {
  actual_value_five_Chem_no_oxy_com <- as.factor(Val_set_chem_no_oxy_com[[var]])
  if (length(levels(actual_value_five_Chem_no_oxy_com)) == 2) {
  prediction_five_Chem_no_oxy_com <- Predictions_Val_Chem_no_oxy_com$classOutput[[var]]$predicted[,2]
  roc_curve <- roc(actual_value_five_Chem_no_oxy_com, prediction_five_Chem_no_oxy_com, direction = "<")
  auc_value <- auc(roc_curve)
  auc_values_five_Chem_no_oxy_com[[var]] <- auc_value
  print(paste("AUC for", var, "is", auc_value))
} else {
    prediction_five_Chem_no_oxy_com <- Predictions_Val_Chem_no_oxy_com$classOutput[[var]]$predicted
    roc_curve <- multiclass.roc(actual_value_five_Chem_no_oxy_com, prediction_five_Chem_no_oxy_com)
    auc_value <- auc(roc_curve)
    auc_values_five_Chem_no_oxy_com[[var]] <- auc_value
    print(paste("AUC for", var, "is", auc_value))
  } 
}
#CF


for (var in response_var_five_Chem_no_oxy_com) {
  actual_values <- as.factor(Val_set_chem_no_oxy_com[[var]])

  predictions_raw <- Predictions_Val_Chem_no_oxy_com$classOutput[[var]]$predicted
  

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
#ACC-Com
Predictions_Val_Chem_no_oxy_com$classOutput$community.number$predicted

actual_classes_no_oxy_com <- Val_set_chem_no_oxy_com$community.number

predicted_classes_no_oxy_com <- apply(Predictions_Val_Chem_no_oxy_com$classOutput$community.number$predicted, 1, which.max)
actual_classes_no_oxy_com <- as.character(actual_classes_no_oxy_com)
predicted_classes_no_oxy_com <- as.character(predicted_classes_no_oxy_com)





#No-OxyCom
removezerodata_meta_NA$community.number <- as.factor(removezerodata_meta_NA$community.number)
removezerodata_meta_NA_no_oxy <- removezerodata_meta_NA[removezerodata_meta_NA$Oxytetracycline == 0, ]
total_communities <- length(unique(removezerodata_meta_NA_no_oxy$community.number))
test_size <- ceiling(total_communities * 0.7)
test_communities <- sample(levels(removezerodata_meta_NA_no_oxy$community.number), size = test_size, replace = FALSE)
Community_Group_TV <- removezerodata_meta_NA_no_oxy
Community_Group_TV$set <- "validation"
Community_Group_TV[Community_Group_TV$community.number %in% test_communities, ]$set <- "train"
Community_Group_Train <- Community_Group_TV[Community_Group_TV$set == "train", ]
Community_Group_Validations <- Community_Group_TV[Community_Group_TV$set == "validation", ]


Train_set_chem_no_oxycom  <- Community_Group_Train %>%
  select(-SampleID,  -community.name,  -community.number, -Oxytetracycline,
         -chem.code, -complexity, -well, -plate,
         -set)


Train_set_chem_no_oxycom$Amoxicillin <- as.factor(Train_set_chem_no_oxycom$Amoxicillin)
Train_set_chem_no_oxycom$Chlorothalonil <- as.factor(Train_set_chem_no_oxycom$Chlorothalonil)
Train_set_chem_no_oxycom$Diflufenican <- as.factor(Train_set_chem_no_oxycom$Diflufenican)
Train_set_chem_no_oxycom$Glyphosate <- as.factor(Train_set_chem_no_oxycom$Glyphosate)
Train_set_chem_no_oxycom$Imidacloprid <- as.factor(Train_set_chem_no_oxycom$Imidacloprid)
Train_set_chem_no_oxycom$Metaldehyde <- as.factor(Train_set_chem_no_oxycom$Metaldehyde)
Train_set_chem_no_oxycom$Tebuconazole <- as.factor(Train_set_chem_no_oxycom$Tebuconazole)



RFSRC_chem_no_oxycom <- rfsrc(Multivar(Amoxicillin, Chlorothalonil, Diflufenican, 
                                        Glyphosate, Imidacloprid, Metaldehyde, Tebuconazole) ~ .,
                               data = Train_set_chem_no_oxycom, ntree = 501, proximity = TRUE)


get.mv.error(RFSRC_chem_no_oxycom, standardize = TRUE)


Val_set_chem_no_oxycom <- Community_Group_Validations  %>%
  select(-SampleID, -community.name,  -chem.code, -complexity, -Oxytetracycline, -community.number,
         -well, -plate)

Val_set_chem_no_oxycom$Amoxicillin <- as.factor(Val_set_chem_no_oxycom$Amoxicillin)
Val_set_chem_no_oxycom$Chlorothalonil <- as.factor(Val_set_chem_no_oxycom$Chlorothalonil)
Val_set_chem_no_oxycom$Diflufenican <- as.factor(Val_set_chem_no_oxycom$Diflufenican)
Val_set_chem_no_oxycom$Glyphosate <- as.factor(Val_set_chem_no_oxycom$Glyphosate)
Val_set_chem_no_oxycom$Imidacloprid <- as.factor(Val_set_chem_no_oxycom$Imidacloprid)
Val_set_chem_no_oxycom$Metaldehyde <- as.factor(Val_set_chem_no_oxycom$Metaldehyde)
Val_set_chem_no_oxycom$Tebuconazole <- as.factor(Val_set_chem_no_oxycom$Tebuconazole)

Predictions_Val_Chem_no_oxycom <- predict(RFSRC_chem_no_oxycom, newdata = Val_set_chem_no_oxycom, type = "response")
Predictions_Val_Chem_no_oxycom

auc_values_five_Chem_no_oxycom <- list()

response_var_five_Chem_no_oxycom <- c("Amoxicillin", "Chlorothalonil", "Diflufenican", "Glyphosate", "Imidacloprid", "Metaldehyde", "Tebuconazole")

for (var in response_var_five_Chem_no_oxycom) {
  actual_value_five_Chem_no_oxycom <- as.factor(Val_set_chem_no_oxycom[[var]])
  if (length(levels(actual_value_five_Chem_no_oxycom)) == 2) {
    prediction_five_Chem_no_oxycom <- Predictions_Val_Chem_no_oxycom$classOutput[[var]]$predicted[,2]
    roc_curve <- roc(actual_value_five_Chem_no_oxycom, prediction_five_Chem_no_oxycom, direction = "<")
    auc_value <- auc(roc_curve)
    auc_values_five_Chem_no_oxycom[[var]] <- auc_value
    print(paste("AUC for", var, "is", auc_value))
  } else {
    prediction_five_Chem_no_oxycom <- Predictions_Val_Chem_no_oxycom$classOutput[[var]]$predicted
    roc_curve <- multiclass.roc(actual_value_five_Chem_no_oxycom, prediction_five_Chem_no_oxycom)
    auc_value <- auc(roc_curve)
    auc_values_five_Chem_no_oxycom[[var]] <- auc_value
    print(paste("AUC for", var, "is", auc_value))
  } 
}
#CF


for (var in response_var_five_Chem_no_oxycom) {
  actual_values <- as.factor(Val_set_chem_no_oxycom[[var]])
  

  predictions_raw <- Predictions_Val_Chem_no_oxycom$classOutput[[var]]$predicted
  

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
