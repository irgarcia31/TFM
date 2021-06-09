# Path where the results of the project are stored
RESULTS_DIR <- "/home/user14/Nextcloud/share2/share/20201229_TFM/results/Validation_set_2/"

# Load final predictions
final_predictions_xgb <- read.csv(file = file.path(RESULTS_DIR, "final_prediction_xgb_set_2.csv"), row.names = 1, stringsAsFactors = FALSE)
final_predictions_rf <- read.csv(file = file.path(RESULTS_DIR, "final_prediction_rf_set_2.csv"), row.names = 1, stringsAsFactors = FALSE)
final_predictions_lda <- read.csv(file = file.path(RESULTS_DIR, "final_prediction_lda_set_2.csv"), row.names = 1, stringsAsFactors = FALSE)

library(dplyr)

# Add 
final_results_lda_set_2 <- final_predictions_lda %>%
  mutate(Sample_Name = gsub("_[0-9]+_R[0-9]+C[0-9]+$", "", rownames(final_predictions_lda))) %>%
  mutate(max_val_lda = apply(final_predictions_lda, 1, max)) %>%
  mutate(class_lda = colnames(final_predictions_lda)[apply(final_predictions_lda, 1, which.max)])

final_results_rf_set_2 <- final_predictions_rf %>%
  mutate(Sample_Name = gsub("_[0-9]+_R[0-9]+C[0-9]+$", "", rownames(final_predictions_rf))) %>%
  mutate(max_val_rf = apply(final_predictions_rf, 1, max)) %>%
  mutate(class_rf = colnames(final_predictions_rf)[apply(final_predictions_rf, 1, which.max)])

final_results_xgb_set_2 <- final_predictions_xgb %>%
  mutate(Sample_Name = gsub("_[0-9]+_R[0-9]+C[0-9]+$", "", rownames(final_predictions_xgb))) %>%
  mutate(max_val_xgb = apply(final_predictions_xgb, 1, max)) %>%
  mutate(class_xgb = colnames(final_predictions_xgb)[apply(final_predictions_xgb, 1, which.max)])

# Generate a final data frame with max_val and class columns
ss_output <- Reduce(function(x, y) merge(x, y, all=TRUE, by = "Sample_Name"), 
                  list(
                    ss, 
                    final_results_lda_set_2[ , c("Sample_Name", "max_val_lda", "class_lda")],
                    final_results_rf_set_2[ , c("Sample_Name", "max_val_rf", "class_rf")],
                    final_results_xgb_set_2[ , c("Sample_Name", "max_val_xgb", "class_xgb")]))

# Save th output
write.csv(ss_output, file = file.path(RESULTS_DIR, "Validation_set_2_predictions.csv"))
