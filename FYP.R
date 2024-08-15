install.packages("BiocManager")
install.packages("remotes")
remotes::install_github("jbisanz/qiime2R", force = TRUE)

install.packages("Rhdf5lib")
install.packages("rhdf5filters")
install.packages("rhdf5")
install.packages("biomformat")
install.packages("phyloseq")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TreeSummarizedExperiment")
BiocManager::install()
BiocManager::install("qiime2R")
BiocManager::install("biomformat")
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("jbisanz/qiime2R")
find.package("qiime2R")
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
install.packages("rsample")
install.packages("sampling")
install.packages("patchwork")
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
library(openxlsx)
library(gridExtra)
library(patchwork)


rarefied_table <- read_qza("D:/Documents/FYP/Data/rarefied_table.qza")
str(rarefied_table)
head(rarefied_table$data)
rarefied_table_data <- rarefied_table$data
rarefied_td_t <- t(rarefied_table_data)

metadata <- read.csv("D:/Documents/FYP/Data/metadata.tsv", sep = "\t") %>%
  rename(SampleID = sample.id)
metadata[metadata$chem.code == "",]$chem.code <- "Control"
taxonomy <- read_qza("D:/Documents/FYP/Data/taxonomy.qza")$data %>% parse_taxonomy()
tree <- read_qza("D:/Documents/FYP/Data/rooted_tree.qza")$data
main_theme <-   theme_minimal() +
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
    legend.position = "right"  
  )
write.xlsx(taxonomy, file = "D:/Documents/FYP/Data/taxonomy.xlsx", colNames = TRUE, rowNames = TRUE)

write.xlsx(merged_data_three_one_amo, file = "D:/Documents/FYP/Data/merged_data_three_one_amo.xlsx", colNames = TRUE, rowNames = TRUE)
write.xlsx(removezerodata_meta_NA, file = "D:/Documents/FYP/Data/removezerodata_meta_NA.xlsx", colNames = TRUE, rowNames = TRUE)
#No_frozen
metadata_noFrozen <-  metadata %>%
  filter(!str_detect(chem.code, "Frozen Starting Community"))




#NMDS

Dis <- vegdist(rarefied_td_t)
nmds_result <- metaMDS(Dis)
plot(nmds_result)
nmds_coordinates <- nmds_result$points
nmds_df <- as.data.frame(nmds_coordinates)
nmds_df$stress <- nmds_result$stress
nmds_df$SampleID<- rownames(nmds_df)
nmds.df <- left_join(nmds_df, metadata)
options(repr.plot.width = 10, repr.plot.height = 6)
nmds_plot <- ggplot(nmds.df %>% filter(chem.code != "Frozen Starting Community"), aes(x = MDS1, y = MDS2, fill = as.character(community.number))) +
  geom_point(aes(shape = as.character(Oxytetracycline)), size = 3) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_viridis_d() +
  labs(x = "nMDS 1", y = "nMDS 2", fill = "Community", shape = "Oxytetracycline") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  main_theme
nmds_plot

nmds_plot <- nmds_plot + 
  stat_ellipse(aes(x = MDS1, y = MDS2, group = as.character(community.number)), 
               type = "norm", linetype = 2, alpha = 0.5)





nmds_plot_com <- ggplot(nmds.df %>% filter(chem.code != "Frozen Starting Community"), aes(x = MDS1, y = MDS2, fill = as.character(complexity))) +
  geom_point(aes(shape = as.character(Oxytetracycline)), size = 3) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_viridis_d() +
  labs(x = "nMDS 1", y = "nMDS 2", fill = "Complexity", shape = "Oxytetracycline") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  main_theme
nmds_plot_com


ggplot(nmds.df %>% filter(community.number == 4 & chem.code != "Frozen Starting Community"), 
       aes(x = MDS1, y = MDS2, fill = as.character(community.number))) +
  geom_point(aes(shape = as.character(Oxytetracycline)), size = 3) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_viridis_d() +
  labs(x = "nMDS 1", y = "nMDS 2", fill = "Community", shape = "Oxytetracycline") +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  main_theme




ggplot(nmds.df %>% filter(community.number == 4 & chem.code != "Frozen Starting Community"), 
       aes(x = MDS1, y = MDS2, fill = as.character(complexity))) + 
  geom_point(shape = 21, size = 3) +
  geom_label_repel(aes(label = chem.code),
                   segment.color = 'grey50',
                   max.overlaps = 20) +
  scale_fill_viridis_d() +
  labs(x = "nMDS 1", y = "nMDS 2", fill = "# Chemicals", title = "Community 4") +
  main_theme

nmds_plot_f <- ggplot(nmds.df %>% filter(community.number == 4 & chem.code != "Frozen Starting Community"), 
       aes(x = MDS1, y = MDS2, fill = as.character(complexity), shape = as.character(Oxytetracycline))) + 
  geom_point(size = 3) +  # 移除了shape参数，因为我们将在scale_shape_manual中定义它
  geom_label_repel(aes(label = chem.code),
                   segment.color = 'grey50',
                   max.overlaps = 20) +  
  scale_fill_viridis_d() +
  scale_shape_manual(values = c(21, 24)) +  
  labs(x = "nMDS 1", y = "nMDS 2", fill = "Chemicals", shape = "Oxytetracycline", title = "Community 4") +
  main_theme

nmds_plot_n <- ggplot(nmds.df %>% filter(community.number == 9 & chem.code != "Frozen Starting Community"), 
       aes(x = MDS1, y = MDS2, fill = as.character(complexity), shape = as.character(Oxytetracycline))) + 
  geom_point(size = 3) +  # 移除了shape参数，因为我们将在scale_shape_manual中定义它
  geom_label_repel(aes(label = chem.code),
                   segment.color = 'grey50',
                   max.overlaps = 20) +  
  scale_fill_viridis_d() +
  scale_shape_manual(values = c(21, 24)) +  
  labs(x = "nMDS 1", y = "nMDS 2", fill = "Chemicals", shape = "Oxytetracycline", title = "Community 9") +
  main_theme




nmds_plot_t <- ggplot(nmds.df %>% filter(community.number == 10 & chem.code != "Frozen Starting Community"), 
aes(x = MDS1, y = MDS2, fill = as.character(complexity), shape = as.character(Oxytetracycline))) + 
  geom_point(size = 3) +
  geom_label_repel(aes(label = chem.code),
                   segment.color = 'grey50',
                   max.overlaps = 20) +
  scale_fill_viridis_d() +
  scale_shape_manual(values = c(21, 24)) +
  labs(x = "nMDS 1", y = "nMDS 2", fill = "Chemicals", shape = "Oxytetracycline", title = "Community 10") +
  main_theme

nmds_plot_tw <- ggplot(nmds.df %>% filter(community.number == 12 & chem.code != "Frozen Starting Community"), 
       aes(x = MDS1, y = MDS2, fill = as.character(complexity), shape = as.character(Oxytetracycline))) + 
  geom_point(size = 3) +
  geom_label_repel(aes(label = chem.code),
                   segment.color = 'grey50',
                   max.overlaps = 20) +
  scale_fill_viridis_d() +
  scale_shape_manual(values = c(21, 24)) + 
  labs(x = "nMDS 1", y = "nMDS 2", fill = "Chemicals", shape = "Oxytetracycline", title = "Community 12") +
  main_theme

nmds_plot_titled <- nmds_plot + ggtitle("A")
nmds_plot_com_titled <- nmds_plot_com + ggtitle("B")
nmds_plot_f <- nmds_plot_f + ggtitle("C")
nmds_plot_n <- nmds_plot_n + ggtitle("D")
nmds_plot_t <-nmds_plot_t + ggtitle("E")
nmds_plot_tw <- nmds_plot_tw + ggtitle("F")


combined_plots <- nmds_plot_titled + 
  nmds_plot_com_titled + 
  nmds_plot_f + 
  nmds_plot_n + 
  nmds_plot_t + 
  nmds_plot_tw +
  plot_layout(ncol = 2) 


print(combined_plots)


# ANOSIM
valid_samples <- metadata$chem.code != "Frozen Starting Community"
anosim_result <- anosim(dist_matrix, groups)



rarefied_td_t_df <- as.data.frame(rarefied_td_t)

rownames(rarefied_td_t_df) <- rownames(rarefied_td_t)
colnames(rarefied_td_t_df) <- colnames(rarefied_td_t)


rarefied_td_t_df <- data.frame(SampleID = rownames(rarefied_td_t_df), rarefied_td_t_df)
combined_rarefiedmeta <- left_join(rarefied_td_t_df, metadata, by = "SampleID")

filtered_combined <- filter(combined_rarefiedmeta, chem.code != "Frozen Starting Community")
valid_samples <- combined_rarefiedmeta$chem.code != "Frozen Starting Community"

length(valid_samples)

Dis <- as.matrix(Dis)
dim(Dis)

valid_indices <- which(valid_samples)

Dis_filtered <- Dis[valid_indices, valid_indices]
chem.code_filtered <- filtered_combined$chem.code
anosim_result <- anosim(Dis_filtered, chem.code_filtered)
print(anosim_result)
community.name_filtered <- filtered_combined$community.name
anosim_result2 <- anosim(Dis_filtered, community.name_filtered)
print(anosim_result2)
#anosim
rows_with_community_number_one <- filtered_combined[filtered_combined$community.number == 1, ]
rarefied_one <- rows_with_community_number_one [,2:810]
Dis_one <- vegdist(rarefied_one, method = "bray")
chem_code_one <- as.factor(rows_with_community_number_one$chem.code)
anosim_result_one <- anosim(Dis_one, chem_code_one)
anosim_result_one

rarefied_all <- filtered_combined[,2:810]
Dis_all <- vegdist(rarefied_all, method = "bray")
chem_code_all <- as.factor(filtered_combined$chem.code)
anosim_result_all <- anosim(Dis_all, chem_code_all, strata= filtered_combined$community.number)
anosim_result_all

community_all <- as.factor(filtered_combined$community.number)
anosim_result_all_com <- anosim(Dis_all, community_all)
anosim_result_all_com

#com4
rows_with_community_number_four <- filtered_combined[filtered_combined$community.number == 4, ]
rarefied_four <- rows_with_community_number_four [,2:810]
Dis_four <- vegdist(rarefied_four, method = "bray")
chem_code_four <- as.factor(rows_with_community_number_four$chem.code)
anosim_result_four <- anosim(Dis_four, chem_code_four)
anosim_result_four


#5
rows_with_community_number_five <- filtered_combined[filtered_combined$community.number == 5, ]
rarefied_five <- rows_with_community_number_five [,2:810]
Dis_five <- vegdist(rarefied_five, method = "bray")
chem_code_five <- as.factor(rows_with_community_number_five$chem.code)
anosim_result_five <- anosim(Dis_five, chem_code_five)
anosim_result_five



#6
rows_with_community_number_six <- filtered_combined[filtered_combined$community.number == 6, ]
rarefied_six <- rows_with_community_number_six [,2:810]
Dis_six <- vegdist(rarefied_six, method = "bray")
chem_code_six <- as.factor(rows_with_community_number_six$chem.code)
anosim_result_six <- anosim(Dis_six, chem_code_six)
anosim_result_six


#7
rows_with_community_number_seven <- filtered_combined[filtered_combined$community.number == 7, ]
rarefied_seven <- rows_with_community_number_seven [,2:810]
Dis_seven <- vegdist(rarefied_seven, method = "bray")
chem_code_seven <- as.factor(rows_with_community_number_seven$chem.code)
anosim_result_seven <- anosim(Dis_seven, chem_code_seven)
anosim_result_seven


#8
rows_with_community_number_eight <- filtered_combined[filtered_combined$community.number == 8, ]
rarefied_eight <- rows_with_community_number_eight [,2:810]
Dis_eight <- vegdist(rarefied_eight, method = "bray")
chem_code_eight <- as.factor(rows_with_community_number_eight$chem.code)
anosim_result_eight <- anosim(Dis_eight, chem_code_eight)
anosim_result_eight




#9
rows_with_community_number_nine <- filtered_combined[filtered_combined$community.number == 9, ]
rarefied_nine <- rows_with_community_number_nine [,2:810]
Dis_nine <- vegdist(rarefied_nine, method = "bray")
chem_code_nine <- as.factor(rows_with_community_number_nine$chem.code)
anosim_result_nine <- anosim(Dis_nine, chem_code_nine)
anosim_result_nine



#10
rows_with_community_number_ten <- filtered_combined[filtered_combined$community.number == 10, ]
rarefied_ten <- rows_with_community_number_ten [,2:810]
Dis_ten <- vegdist(rarefied_ten, method = "bray")
chem_code_ten <- as.factor(rows_with_community_number_ten$chem.code)
anosim_result_ten <- anosim(Dis_ten, chem_code_ten)
anosim_result_ten



#11
rows_with_community_number_eleven <- filtered_combined[filtered_combined$community.number == 11, ]
rarefied_eleven <- rows_with_community_number_eleven [,2:810]
Dis_eleven <- vegdist(rarefied_eleven, method = "bray")
chem_code_eleven <- as.factor(rows_with_community_number_eleven$chem.code)
anosim_result_eleven <- anosim(Dis_eleven, chem_code_eleven)
anosim_result_eleven



#12
rows_with_community_number_tweleve <- filtered_combined[filtered_combined$community.number == 12, ]
rarefied_tweleve <- rows_with_community_number_tweleve [,2:810]
Dis_tweleve <- vegdist(rarefied_tweleve, method = "bray")
chem_code_tweleve <- as.factor(rows_with_community_number_tweleve$chem.code)
anosim_result_tweleve <- anosim(Dis_tweleve, chem_code_tweleve)
anosim_result_tweleve

#PERMANOVA

rarefied_td_t_f <- rarefied_td_t_df[filtered_combined$SampleID,]
row.names(rarefied_td_t_f) == filtered_combined$SampleID
rare_meta <- filtered_combined[, c(1, 811:824)]
row.names(rarefied_td_t_f) == rare_meta$SampleID
view(rarefied_td_t_f)
rarefied_td_t_f <- subset(rarefied_td_t_f, select = -SampleID)
adonis2(rarefied_td_t_f ~ community.number + (Amoxicillin + Chlorothalonil + Diflufenican +
                                              Glyphosate + Imidacloprid + Metaldehyde + Tebuconazole) * Oxytetracycline ,
        data = rare_meta, method = "bray")
adonis2(rarefied_td_t_f ~ community.number * (Amoxicillin + Chlorothalonil + Diflufenican +
                                              Glyphosate + Imidacloprid + Metaldehyde + Oxytetracycline + Tebuconazole),
        data = rare_meta, method = "bray")
adonis2(rarefied_td_t_f ~ community.number + Amoxicillin + Chlorothalonil + Diflufenican +
          Glyphosate + Imidacloprid + Metaldehyde + Oxytetracycline + Tebuconazole,
        data = rare_meta, method = "bray")
adonis2(rarefied_td_t_f ~ community.number * chem.code,
         data = rare_meta, method = "bray" )
view(filtered_combined)
filtered_combined_no_oxy <- filtered_combined[filtered_combined$Oxytetracycline == 0, ]
rare_meta_no_oxy <- filtered_combined_no_oxy[, c(1, 811:824)]
rare_no_oxy <- filtered_combined_no_oxy[, 2:810]
adonis2(rare_no_oxy ~ community.number *( Amoxicillin + Chlorothalonil + Diflufenican +
                                                Glyphosate + Imidacloprid + Metaldehyde + Tebuconazole),
        data = rare_meta_no_oxy, method = "bray")

filtered_combined$community.number <- as.factor(filtered_combined$community.number)
adonis_result_CC <- adonis2(Dis_filtered ~ community.number + chem.code, data = filtered_combined, permutations = 999)

Frozen_rarefiedmeta <- combined_rarefiedmeta[combined_rarefiedmeta$chem.code == "Frozen Starting Community", ]
rarefied_td_t_frozen <- Frozen_rarefiedmeta[, 2:810]
metadata_frozen <- Frozen_rarefiedmeta[, c(1, 811:824)]
metadata_frozen$community.number <- as.factor(metadata_frozen$community.number)
adonis2(rarefied_td_t_frozen ~ community.number,
        data = metadata_frozen, method = "bray")

metadata_frozen$community.number <- as.factor(metadata_frozen$community.number)


filtered_rarefied_td_t_df <- rarefied_td_t_df[valid_samples, ]
long_data <- gather(filtered_rarefied_td_t_df, key = "ASV", value = "Abundance", -SampleID)
long_meta <- left_join(long_data, rarefied_metadata, by = "SampleID")
lmer_model <- lmer(Abundance ~ Amoxicillin + Chlorothalonil + Diflufenican + Glyphosate + Imidacloprid + Metaldehyde + Tebuconazole + Oxytetracycline + (1|community.name),
                   data = long_meta)
emmeans(lmer_model, ~ Amoxicillin * Chlorothalonil * Diflufenican * Glyphosate * Imidacloprid * Metaldehyde * Tebuconazole)
plot(emmeans(lmer_model, pairwise ~ .), lwd = 2)
emm_results <- emmeans(lmer_model, ~ fixed_effect_variable)

##RF
dim(rarefied_table_data)
ASV_nonzero_counts <- apply(rarefied_table_data, 1, function(y) sum(length(which(y > 0))))
hist(ASV_nonzero_counts, breaks=100, col="grey", main="", ylab="Number of ASVs", xlab="Number of Non-Zero Values")
remove_rare <- function( table , cutoff_pro ) {
  row2keep <- c()
  cutoff <- ceiling( cutoff_pro * ncol(table) )  
  for ( i in 1:nrow(table) ) {
    row_nonzero <- length( which( table[ i , ]  > 0 ) ) 
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  return( table [ row2keep , , drop=F ])
}

ASV_table_rare_removed <- remove_rare(table=rarefied_table_data, cutoff_pro=0.1) # remove OTUs not present in 10% of samples
dim(ASV_table_rare_removed) # we've excluded almost 90% of our OTUs.. that seems a bit concerning
# prepare data for model
# just try for oxytet first
asv_table_oxytet <- data.frame(t(ASV_table_rare_removed)) 
asv_table_oxytet$SampleID <- rownames(asv_table_oxytet)
asv_table_oxytet <- left_join(asv_table_oxytet,
                              metadata %>%
                                dplyr::select(SampleID, Oxytetracycline)) %>%
  dplyr::select(-SampleID) %>%
  filter(!is.na(Oxytetracycline))

# needs to be a factor for classification (rather than regression)
asv_table_oxytet$Oxytetracycline <- as.factor(asv_table_oxytet$Oxytetracycline)

# try to run the model
RF_state_classify <- randomForest(Oxytetracycline ~., data = asv_table_oxytet, ntree=501, importance=TRUE, proximities=TRUE )
RF_state_classify

#Selection
removezerodata <- data.frame(t(ASV_table_rare_removed))
removezerodata$SampleID <- rownames(removezerodata)
removezerodata_meta <- left_join(removezerodata, metadata, by = "SampleID")


# 移除指定列
removezerodata_meta_NA <- removezerodata_meta %>%
  filter(!is.na(Oxytetracycline))
removezerodata_meta_clean <- removezerodata_meta_NA %>%
  select(-c(SampleID, community.number, community.name, 
            Amoxicillin, Chlorothalonil, Diflufenican, 
            Glyphosate, Imidacloprid, Metaldehyde, 
            Oxytetracycline, Tebuconazole, complexity, well, plate))
#CommunityRF-Rep
set.seed(123) 
removezerodata_meta_clean$chem.code <- as.factor(removezerodata_meta_clean$chem.code)
trainIndex <- createDataPartition(removezerodata_meta_clean$chem.code, p = 0.7, list = FALSE)
train_data <- removezerodata_meta_clean[trainIndex, ]
test_data <- removezerodata_meta_clean[-trainIndex, ]
RF_model <- randomForest(chem.code ~ ., data = train_data, ntree=501, importance=TRUE)
RF_model


##No-Rep
set.seed(123) 
#pergroup-percode
removezerodata_meta_NA_g <- removezerodata_meta_NA %>%
  group_by(chem.code) %>%
  mutate(group_id = cur_group_id()) %>%  
  ungroup()
g_max_value <- max(removezerodata_meta_NA_g$group_id)
num_test_groups <- ceiling(0.3 * g_max_value)
group_id_unique <- unique(removezerodata_meta_NA_g$group_id)
group_id_shuffled <- sample(group_id_unique, replace = FALSE)
test_groups <- head(group_id_shuffled, num_test_groups)
validation_groups <- setdiff(group_id_shuffled, test_groups)

group_assignment <- tibble(
  group_id = group_id_shuffled,  
  set = ifelse(group_id_shuffled %in% test_groups, "test", "validation")  
)
removezerodata_meta_NA_g_tv <- left_join(removezerodata_meta_NA_g, group_assignment, by = "group_id")
NoRep_test_data <- removezerodata_meta_NA_g_tv %>% filter(set == "test")
NoRep_validation_data <- removezerodata_meta_NA_g_tv %>% filter(set == "validation")

# Remove some col
NoRep_test_data_clean <- NoRep_test_data %>%
  select(-SampleID, -community.number, -community.name, 
         -Amoxicillin, -Chlorothalonil, -Diflufenican, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, 
         -Oxytetracycline, -Tebuconazole, -complexity, -well, -plate,
         -group_id, -set)
NoRep_test_data_clean$chem.code <- as.factor(NoRep_test_data_clean$chem.code)
RF_model_NoRep_c <- randomForest(chem.code ~ ., data = NoRep_test_data_clean, ntree=501, importance=TRUE)
RF_model_NoRep_c

importance(RF_model_NoRep_c)
varImpPlot(RF_model_NoRep_c)
importance_matrix <- importance(RF_model_NoRep_c)
print(importance_matrix)
NoRep_test_data_clean_Oxy <- NoRep_test_data %>%
  select(-SampleID, -community.number, -community.name, 
         -chem.code, -Amoxicillin, -Chlorothalonil, -Diflufenican, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, 
         -Tebuconazole, -complexity, -well, -plate,
         -group_id, -set)
NoRep_test_data_clean_Oxy$Oxytetracycline <- as.factor(NoRep_test_data_clean_Oxy$Oxytetracycline)
RF_model_NoRep_Oxy <- randomForest(Oxytetracycline ~., data = NoRep_test_data_clean_Oxy, ntree=501, importance=TRUE, proximities=TRUE )
RF_model_NoRep_Oxy

IRF_NR_O_df_MDA_Sorted<- arrange( IRF_NR_O_df  , desc(MeanDecreaseAccuracy)  )
IRF_NR_O_df_MDA_top_features <- IRF_NR_O_df_MDA_Sorted[1:10,"Variable"]
taxonomy %>% filter(row.names(taxonomy) %in% c(IRF_NR_O_df_MDA_top_features,
                                               "6377e612d1fdfdadf18255fc486ce608",
                                               "5203c525a0d0a6c831885ac86a569081",
                                               "4c294f9a88e58b0e40a2b19f4658c209"))


IRF_NR_O <- importance(RF_model_NoRep_Oxy )
IRF_NR_C <- importance(RF_model_NoRep_c)

IRF_NR_O_RNA <- setdiff(names(NoRep_test_data_clean_Oxy), "Oxytetracycline")
IRF_NR_C_RNA <- setdiff(names(NoRep_test_data_clean), "chem.code")  
str(IRF_NR_O)
IRF_NR_O_df <- as.data.frame(IRF_NR_O)
IRF_NR_C_df <- as.data.frame(IRF_NR_C)
IRF_NR_O_df$Variable <- IRF_NR_O_RNA
IRF_NR_C_df$Variable <- IRF_NR_C_RNA
IRF_NR_C_sorted <- IRF_NR_C_df[order(-IRF_NR_C_df$MeanDecreaseGini), ]
IRF_NR_O_sorted <- IRF_NR_O_df[order(-IRF_NR_O_df$MeanDecreaseGini), ]
print(IRF_NR_O_sorted)
print(IRF_NR_C_sorted)

ggplot(IRF_NR_C_df, aes(x = Variable, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity") +
  labs(title = "Variable Importance - IRF_NR_O",
       x = "Variable",
       y = "MeanDecreaseGini")
ggplot(IRF_NR_O_df, aes(x = Variable, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity") +
  labs(title = "Variable Importance - IRF_NR_O",
       x = "Variable",
       y = "MeanDecreaseGini")

IRF_NR_O_df_model <- IRF_NR_O_df
IRF_NR_O_df_model$Model <- "ModelO"
IRF_NR_C_df_model <- IRF_NR_C_df
IRF_NR_C_df_model$Model <- "ModelC"
IRF_NR_O_df_model <- IRF_NR_O_df_model[, c("MeanDecreaseGini", "MeanDecreaseAccuracy", "Variable", "Model")]
IRF_NR_C_df_model <- IRF_NR_C_df_model[, c("MeanDecreaseGini", "MeanDecreaseAccuracy", "Variable", "Model")]
combined_importance_IRF_NR_CO <- rbind(IRF_NR_C_df_model, IRF_NR_O_df_model)
ggplot(combined_importance_IRF_NR_CO, aes(x=Variable, y=MeanDecreaseGini, fill=Model)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(title="Comparison of Variable Importance Across Two Models",
       x="Variable",
       y="Mean Decrease in Gini",
       fill="Model")

ranking_IRF_NR_C <- rank(-IRF_NR_C_df$MeanDecreaseGini, ties.method = "first")
ranking_IRF_NR_O<- rank(-IRF_NR_O_df$MeanDecreaseGini, ties.method = "first") 
comparison <- data.frame(Variable = IRF_NR_C_sorted$Variable,
                         Ranking_C = ranking_IRF_NR_C,
                         Ranking_O = ranking_IRF_NR_O)

comparison$Rank_Difference <- comparison$Ranking_C - comparison$Ranking_O                                             
comparison$Variable <- as.factor(comparison$Variable)
ggplot(comparison, aes(x = Ranking_C, y = Ranking_O, color = Variable)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  theme_minimal() +
  labs(title = "Variable Importance Ranking Comparison",
       x = "Ranking in ranking_IRF_NR_C",
       y = "Ranking in ranking_IRF_NR_O") +
  theme(legend.position = "top") 

NoRep_test_data_c_c <- NoRep_test_data %>%
  select(-SampleID, -community.name, 
         -Amoxicillin, -Chlorothalonil, -Diflufenican, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, 
         -Oxytetracycline, -Tebuconazole, -complexity, -well, -plate,
         -group_id, -set)
NoRep_test_data_c_c$chem.code <- as.factor(NoRep_test_data_c_c$chem.code)
NoRep_test_data_c_c$community.number <- as.factor(NoRep_test_data_c_c$community.number)
RF_model_NoRep_c_c <- randomForest(chem.code ~ ., data = NoRep_test_data_c_c, ntree=501, importance=TRUE)
RF_model_NoRep_c_c


cor_matrix <- cor(removezerodata_meta_NA[, -which(names(removezerodata_meta_NA) == "chem.code")], use = "pairwise.complete.obs")

# 找出与每个RNA群落特征相关性最高的化学物质
cor_chem <- apply(cor_matrix, 2, function(x) sort(x, decreasing = TRUE)[1])



#Oxytetracycline+community.name
#chem.code+community

#RF-DvaraiblesImportance
#NoRep_test_data_clean_DVI <- NoRep_test_data %>%
  select(-SampleID, -community.name, 
         -Amoxicillin, -Chlorothalonil, -Diflufenican, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, 
         -Oxytetracycline, -Tebuconazole, -complexity, -well, -plate,
         -group_id, -set)
NoRep_test_data_clean_DVI$chem.code <- as.factor(NoRep_test_data_clean$chem.code)
NoRep_test_data_clean_DVI$community.number <- as.factor(NoRep_test_data_clean_DVI$community.number)


#RF-NoRep_chemicals
#Oxy
set.seed(123)
NoRep_test_data_clean_Oxy <- NoRep_test_data %>%
  select(-SampleID, -community.number, -community.name, 
         -chem.code, -Amoxicillin, -Chlorothalonil, -Diflufenican, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, 
         -Tebuconazole, -complexity, -well, -plate,
         -group_id, -set)
NoRep_test_data_clean_Oxy$Oxytetracycline <- as.factor(NoRep_test_data_clean_Oxy$Oxytetracycline)
RF_model_NoRep_Oxy <- randomForest(Oxytetracycline ~., data = NoRep_test_data_clean_Oxy, ntree=501, importance=TRUE, proximities=TRUE )
RF_model_NoRep_Oxy
#Amo
NoRep_test_data_clean_Amo <- NoRep_test_data %>%
  select(-SampleID, -community.number, -community.name, 
         -chem.code, -Chlorothalonil, -Diflufenican, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, 
         -Oxytetracycline, -Tebuconazole, -complexity, -well, -plate,
         -group_id, -set)
NoRep_test_data_clean_Amo$Amoxicillin <- as.factor(NoRep_test_data_clean_Amo$Amoxicillin)
RF_model_NoRep_Amo <- randomForest(Amoxicillin ~., data = NoRep_test_data_clean_Amo, ntree=501, importance=TRUE, proximities=TRUE )
RF_model_NoRep_Amo
#Chl
NoRep_test_data_clean_Chl <- NoRep_test_data %>%
  select(-SampleID, -community.number, -community.name, 
         -chem.code, -Amoxicillin, -Diflufenican, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, 
         -Oxytetracycline, -Tebuconazole, -complexity, -well, -plate,
         -group_id, -set)
NoRep_test_data_clean_Chl$Chlorothalonil <- as.factor(NoRep_test_data_clean_Chl$Chlorothalonil)
RF_model_NoRep_Chl <- randomForest(Chlorothalonil ~., data = NoRep_test_data_clean_Chl, ntree=501, importance=TRUE, proximities=TRUE )
RF_model_NoRep_Chl
#Dif
NoRep_test_data_clean_Dif <- NoRep_test_data %>%
  select(-SampleID, -community.number, -community.name, 
         -chem.code, -Amoxicillin, -Chlorothalonil,  
         -Glyphosate, -Imidacloprid, -Metaldehyde, 
         -Oxytetracycline, -Tebuconazole, -complexity, -well, -plate,
         -group_id, -set)
NoRep_test_data_clean_Dif$Diflufenican <- as.factor(NoRep_test_data_clean_Dif$Diflufenican)
RF_model_NoRep_Dif <- randomForest(Diflufenican ~., data = NoRep_test_data_clean_Dif, ntree=501, importance=TRUE, proximities=TRUE )
RF_model_NoRep_Dif
#Gly
NoRep_test_data_clean_Gly <- NoRep_test_data %>%
  select(-SampleID, -community.number, -community.name, 
         -chem.code, -Amoxicillin, -Chlorothalonil, -Diflufenican, 
         -Imidacloprid, -Metaldehyde, 
         -Oxytetracycline, -Tebuconazole, -complexity, -well, -plate,
         -group_id, -set)
NoRep_test_data_clean_Gly$Glyphosate <- as.factor(NoRep_test_data_clean_Gly$Glyphosate)
RF_model_NoRep_Gly <- randomForest(Glyphosate ~., data = NoRep_test_data_clean_Gly, ntree=501, importance=TRUE, proximities=TRUE )
RF_model_NoRep_Gly
#Imi
NoRep_test_data_clean_Imi <- NoRep_test_data %>%
  select(-SampleID, -community.number, -community.name, 
         -chem.code, -Amoxicillin, -Chlorothalonil, -Diflufenican, 
         -Glyphosate, -Metaldehyde, 
         -Oxytetracycline, -Tebuconazole, -complexity, -well, -plate,
         -group_id, -set)
NoRep_test_data_clean_Imi$Imidacloprid <- as.factor(NoRep_test_data_clean_Imi$Imidacloprid)
RF_model_NoRep_Imi <- randomForest(Imidacloprid ~., data = NoRep_test_data_clean_Imi, ntree=501, importance=TRUE, proximities=TRUE )
RF_model_NoRep_Imi
#Met
NoRep_test_data_clean_Met <- NoRep_test_data %>%
  select(-SampleID, -community.number, -community.name, 
         -chem.code, -Amoxicillin, -Chlorothalonil, -Diflufenican, 
         -Glyphosate, -Imidacloprid, 
         -Oxytetracycline, -Tebuconazole, -complexity, -well, -plate,
         -group_id, -set)
NoRep_test_data_clean_Met$Metaldehyde <- as.factor(NoRep_test_data_clean_Met$Metaldehyde)
RF_model_NoRep_Met <- randomForest(Metaldehyde ~., data = NoRep_test_data_clean_Met, ntree=501, importance=TRUE, proximities=TRUE )
RF_model_NoRep_Met
#Teb
NoRep_test_data_clean_Teb <- NoRep_test_data %>%
  select(-SampleID, -community.number, -community.name, 
         -chem.code, -Amoxicillin, -Chlorothalonil, -Diflufenican, 
         -Glyphosate, -Imidacloprid, 
         -Oxytetracycline, -Metaldehyde, -complexity, -well, -plate,
         -group_id, -set)
NoRep_test_data_clean_Teb$Tebuconazole <- as.factor(NoRep_test_data_clean_Teb$Tebuconazole)
RF_model_NoRep_Teb <- randomForest(Tebuconazole ~., data = NoRep_test_data_clean_Teb, ntree=501, importance=TRUE, proximities=TRUE )
RF_model_NoRep_Teb

#Permutation analysis
RF_permute_Norep_Oxy <- rf.significance( x=RF_model_NoRep_Oxy ,  xdata=NoRep_test_data_clean_Oxy[,1:(ncol(NoRep_test_data_clean_Oxy)-1)] , nperm=1000 , ntree=501 )  
RF_permute_Norep_Oxy

NoRep_validation_data_clean_Oxy <- NoRep_validation_data %>%
  select(-SampleID, -community.number, -community.name, 
         -chem.code, -Amoxicillin, -Chlorothalonil, -Diflufenican, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, 
         -Tebuconazole, -complexity, -well, -plate,
         -group_id, -set)
NoRep_validation_data_clean_Oxy$Oxytetracycline <- as.factor(NoRep_validation_data_clean_Oxy$Oxytetracycline)
Predictions_NoRep_validation_data_clean_Oxy <- predict(RF_model_NoRep_Oxy, NoRep_validation_data_clean_Oxy)


str(Predictions_NoRep_validation_data_clean_Oxy)
Accuracy_Predictions_NoRep_validation_data_clean_Oxy <- sum(Predictions_NoRep_validation_data_clean_Oxy == NoRep_validation_data_clean_Oxy$Oxytetracycline) / nrow(NoRep_validation_data_clean_Oxy)
levels(Predictions_NoRep_validation_data_clean_Oxy)

##Leave-one-out cross-validation

fit_control <- trainControl(method = "LOOCV")
NoRep_test_data_clean_Oxy_Com <- NoRep_test_data %>%
  select(-SampleID, -community.name, 
         -chem.code, -Amoxicillin, -Chlorothalonil, -Diflufenican, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, 
         -Tebuconazole, -complexity, -well, -plate,
         -group_id, -set)
NoRep_test_data_clean_Oxy_Com$community.number <- as.factor(NoRep_test_data_clean_Oxy_Com$community.number)
unique_communities <- unique(NoRep_test_data_clean_Oxy_Com$community.number)

NoRep_test_data_clean_Oxy_Com$Oxytetracycline <- as.factor(NoRep_test_data_clean_Oxy_Com$Oxytetracycline)
NoRep_test_data_clean_Oxy_Com_Loocv_Results <- list()

for (i in seq_along(unique_communities)) {
 
  test_data <- NoRep_test_data_clean_Oxy_Com[NoRep_test_data_clean_Oxy_Com$community.number == unique_communities[i], ]
  train_data <- NoRep_test_data_clean_Oxy_Com[!NoRep_test_data_clean_Oxy_Com$community.number %in% unique_communities[i], ]

NoRep_test_data_clean_Oxy_Com_Loocv <- train(
    x = train_data[, -which(names(train_data) %in% c("Oxytetracycline", "community.number"))],
    y = train_data$Oxytetracycline,
    method = "rf",
    ntree = 501,
    trControl = fit_control
  )


NoRep_test_data_clean_Oxy_Com_Loocv_Results[[unique_communities[i]]] <- NoRep_test_data_clean_Oxy_Com_Loocv
}
NoRep_test_data_clean_Oxy_Com_Loocv_Results
#RFSRC
NoRep_test_data_clean_Chem_Com <- NoRep_test_data %>%
  select(-SampleID, -community.name, 
         -chem.code, -complexity, -well, -plate,
         -group_id, -set)
NoRep_test_data_clean_Chem <- NoRep_test_data %>%
  select(-SampleID, -community.name, -community.number,
         -chem.code, -complexity, -well, -plate,
         -group_id, -set)

NoRep_test_data_clean_Chem$Amoxicillin <- as.factor(NoRep_test_data_clean_Chem$Amoxicillin)
NoRep_test_data_clean_Chem$Chlorothalonil <- as.factor(NoRep_test_data_clean_Chem$Chlorothalonil)
NoRep_test_data_clean_Chem$Diflufenican <- as.factor(NoRep_test_data_clean_Chem$Diflufenican)
NoRep_test_data_clean_Chem$Glyphosate <- as.factor(NoRep_test_data_clean_Chem$Glyphosate)
NoRep_test_data_clean_Chem$Imidacloprid <- as.factor(NoRep_test_data_clean_Chem$Imidacloprid)
NoRep_test_data_clean_Chem$Metaldehyde <- as.factor(NoRep_test_data_clean_Chem$Metaldehyde)
NoRep_test_data_clean_Chem$Oxytetracycline <- as.factor(NoRep_test_data_clean_Chem$Oxytetracycline)
NoRep_test_data_clean_Chem$Tebuconazole <- as.factor(NoRep_test_data_clean_Chem$Tebuconazole)

NoRep_Test_data_Chem_RFSRC <- rfsrc(Multivar( Amoxicillin, Chlorothalonil, Diflufenican, Glyphosate, Imidacloprid,
                                                 Metaldehyde, Oxytetracycline, Tebuconazole) ~ ., 
                                        data = NoRep_test_data_clean_Chem, improtance = TRUE, proximity = TRUE, ntree = 501)





Response_NoRep_Test_data_Chem_Com <- c("community.number", "Amoxicillin", "Chlorothalonil", "Diflufenican", "Glyphosate", "Imidacloprid",
                                       "Metaldehyde", "Oxytetracycline", "Tebuconazole")

NoRep_test_data_clean_Chem_Com$community.number <- as.factor(NoRep_test_data_clean_Chem_Com$community.number)
NoRep_test_data_clean_Chem_Com$Amoxicillin <- as.factor(NoRep_test_data_clean_Chem_Com$Amoxicillin)
NoRep_test_data_clean_Chem_Com$Chlorothalonil <- as.factor(NoRep_test_data_clean_Chem_Com$Chlorothalonil)
NoRep_test_data_clean_Chem_Com$Diflufenican <- as.factor(NoRep_test_data_clean_Chem_Com$Diflufenican)
NoRep_test_data_clean_Chem_Com$Glyphosate <- as.factor(NoRep_test_data_clean_Chem_Com$Glyphosate)
NoRep_test_data_clean_Chem_Com$Imidacloprid <- as.factor(NoRep_test_data_clean_Chem_Com$Imidacloprid)
NoRep_test_data_clean_Chem_Com$Metaldehyde <- as.factor(NoRep_test_data_clean_Chem_Com$Metaldehyde)
NoRep_test_data_clean_Chem_Com$Oxytetracycline <- as.factor(NoRep_test_data_clean_Chem_Com$Oxytetracycline)
NoRep_test_data_clean_Chem_Com$Tebuconazole <- as.factor(NoRep_test_data_clean_Chem_Com$Tebuconazole)
str(NoRep_test_data_clean_Chem_Com)
NoRep_Test_data_Chem_Com_RFSRC <- rfsrc(Multivar(community.number, Amoxicillin, Chlorothalonil, Diflufenican, Glyphosate, Imidacloprid,
                                                 Metaldehyde, Oxytetracycline, Tebuconazole) ~ ., 
                                        data = NoRep_test_data_clean_Chem_Com, ntree = 501)
NoRep_test_data_clean_Oxy <- NoRep_test_data_clean_Chem_Com%>%
  select(-community.number,-Amoxicillin, -Chlorothalonil, -Diflufenican, 
         -Glyphosate, -Imidacloprid, -Metaldehyde, 
         -Tebuconazole)
NoRep_test_data_clean_Oxy$Oxytetracycline <- as.factor(NoRep_test_data_clean_Oxy$Oxytetracycline)
NoRep_Test_data_Oxy_RFSRC <- rfsrc(Oxytetracycline ~ ., 
                                     data = NoRep_test_data_clean_Oxy, ntree = 501)
str(NoRep_test_data_clean_Oxy) 

summary(NoRep_Test_data_Chem_Com_RFSRC)
removezerodata_meta_NA_clean4 <- removezerodata_meta_NA %>%
  select(-SampleID, -community.name, 
         -chem.code, -complexity, -well, -plate)
removezerodata_meta_NA_clean4$community.number <- as.factor(removezerodata_meta_NA_clean4$community.number)
removezerodata_meta_NA_clean4$Amoxicillin <- as.factor(removezerodata_meta_NA_clean4$Amoxicillin)
removezerodata_meta_NA_clean4$Chlorothalonil <- as.factor(removezerodata_meta_NA_clean4$Chlorothalonil)
removezerodata_meta_NA_clean4$Diflufenican <- as.factor(removezerodata_meta_NA_clean4$Diflufenican)
removezerodata_meta_NA_clean4$Glyphosate <- as.factor(removezerodata_meta_NA_clean4$Glyphosate)
removezerodata_meta_NA_clean4$Imidacloprid <- as.factor(removezerodata_meta_NA_clean4$Imidacloprid)
removezerodata_meta_NA_clean4$Metaldehyde <- as.factor(removezerodata_meta_NA_clean4$Metaldehyde)
removezerodata_meta_NA_clean4$Oxytetracycline <- as.factor(removezerodata_meta_NA_clean4$Oxytetracycline)
removezerodata_meta_NA_clean4$Tebuconazole <- as.factor(removezerodata_meta_NA_clean4$Tebuconazole)
NoRep_Test_data_Oxy_RFSRC <- rfsrc(Oxytetracycline ~ ., 
                                   data = removezerodata_meta_NA_clean4, ntree = 501)
str(removezerodata_meta_NA_clean4)
#Predicition-RFSRC
NoRep_validation_data_clean_Chem_Com <- NoRep_validation_data %>%
  select(-SampleID, -community.name, -chem.code, 
          -complexity, -well, -plate,
         -group_id, -set)
Predictions_NoRep_validation_data_clean_Chem_Com <- predict(NoRep_Test_data_Chem_Com_RFSRC, NoRep_validation_data_clean_Chem_Com)

NoRep_validation_data_clean_Chem_Com$community.number <- as.factor(NoRep_validation_data_clean_Chem_Com$community.number)
NoRep_validation_data_clean_Chem_Com$Amoxicillin <- as.factor(NoRep_validation_data_clean_Chem_Com$Amoxicillin)
NoRep_validation_data_clean_Chem_Com$Chlorothalonil <- as.factor(NoRep_validation_data_clean_Chem_Com$Chlorothalonil)
NoRep_validation_data_clean_Chem_Com$Diflufenican <- as.factor(NoRep_validation_data_clean_Chem_Com$Diflufenican)
NoRep_validation_data_clean_Chem_Com$Glyphosate <- as.factor(NoRep_validation_data_clean_Chem_Com$Glyphosate)
NoRep_validation_data_clean_Chem_Com$Imidacloprid <- as.factor(NoRep_validation_data_clean_Chem_Com$Imidacloprid)
NoRep_validation_data_clean_Chem_Com$Metaldehyde <- as.factor(NoRep_validation_data_clean_Chem_Com$Metaldehyde)
NoRep_validation_data_clean_Chem_Com$Oxytetracycline <- as.factor(NoRep_validation_data_clean_Chem_Com$Oxytetracycline)
NoRep_validation_data_clean_Chem_Com$Tebuconazole <- as.factor(NoRep_validation_data_clean_Chem_Com$Tebuconazole)
str(NoRep_validation_data_clean_Chem_Com)
Predictions_NoRep_validation_data_clean_Chem_Com <- predict(NoRep_Test_data_Chem_Com_RFSRC, NoRep_validation_data_clean_Chem_Com)
Predictions_NoRep_validation_data_clean_Chem_Com

NoRep_validation_data_clean_Chem <- NoRep_validation_data %>%
  select(-SampleID, -community.name, -community.number, -chem.code, 
         -complexity, -well, -plate,
         -group_id, -set)
Predictions_NoRep_validation_data_clean_Chem <- predict(NoRep_Test_data_Chem_RFSRC, NoRep_validation_data_clean_Chem)

#Confusion Matrix
for (response_var in Response_NoRep_Test_data_Chem_Com) {
     response_values <- as.factor(NoRep_validation_data_clean_Chem_Com[[response_var]])
     predictions_m <- round(Predictions_NoRep_validation_data_clean_Chem_Com$regrOutput[[response_var]]$predicted)
  
     predictions_m <- factor(predictions_m, levels = levels(response_values))
   
    confusion_matrix <- confusionMatrix(predictions_m, response_values)
    
    print(paste("Confusion Matrix for", response_var))
    print(confusion_matrix)
  } 
#Confusion Matrix-Accuracy
accuracies <- c()

for (var_name in Response_NoRep_Test_data_Chem_Com) {

  actual_values <- as.factor(NoRep_validation_data_clean_Chem_Com[[var_name]])

  predictions <- round(Predictions_NoRep_validation_data_clean_Chem_Com$regrOutput[[var_name]]$predicted)
  

  predictions <- factor(predictions, levels = levels(actual_values))

  conf_matrix <- confusionMatrix(predictions, actual_values)

  accuracy <- conf_matrix$overall['Accuracy']
  cat(sprintf("Variable: %s, Accuracy: %f\n", var_name, accuracy))
}


##Leave-one-out cross-validation
fit_control <- trainControl(method = "LOOCV")
NoRep_test_data_clean_Chem_Com <- NoRep_test_data %>%
  select(-SampleID, -community.name, 
         -chem.code, -complexity, -well, -plate,
         -group_id, -set)
NoRep_test_data_clean_Oxy_Com$community.number <- as.factor(NoRep_test_data_clean_Oxy_Com$community.number)
unique_communities <- unique(NoRep_test_data_clean_Oxy_Com$community.number)

NoRep_test_data_clean_Oxy_Com$Oxytetracycline <- as.factor(NoRep_test_data_clean_Oxy_Com$Oxytetracycline)
NoRep_test_data_clean_Oxy_Com_Loocv_Results <- list()

for (i in seq_along(unique_communities)) {

  test_data <- NoRep_test_data_clean_Oxy_Com[NoRep_test_data_clean_Oxy_Com$community.number == unique_communities[i], ]
  train_data <- NoRep_test_data_clean_Oxy_Com[!NoRep_test_data_clean_Oxy_Com$community.number %in% unique_communities[i], ]

  NoRep_test_data_clean_Oxy_Com_Loocv <- train(
    x = train_data[, -which(names(train_data) %in% c("Oxytetracycline", "community.number"))],
    y = train_data$Oxytetracycline,
    method = "rf",
    ntree = 501,
    trControl = fit_control
  )
  

  NoRep_test_data_clean_Oxy_Com_Loocv_Results[[unique_communities[i]]] <- NoRep_test_data_clean_Oxy_Com_Loocv
}
NoRep_test_data_clean_Oxy_Com_Loocv_Results

#ROC and AUC-factor
auc_values <- list()

response_var <- c("community.number", "Amoxicillin", "Chlorothalonil", "Diflufenican", "Glyphosate", "Imidacloprid", "Metaldehyde", "Oxytetracycline", "Tebuconazole")
for (var in response_var) {
    
    actual_value <- as.factor(NoRep_validation_data_clean_Chem_Com[[var_name]])
       
        prediction <- Predictions_NoRep_validation_data_clean_Chem_Com$regrOutput[[var_name]]$predicted
        roc_curve <- roc(actual_value, prediction)
        auc_value <- auc(roc_curve)
        
    
          auc_values[[var]] <- auc_value
          
       
            print(paste("AUC for", var, "is", auc_value))
}














