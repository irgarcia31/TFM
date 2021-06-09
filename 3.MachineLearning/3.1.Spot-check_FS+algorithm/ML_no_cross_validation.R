#===============================================================================
# Machine Learning Workflow
# Author: Inés Rodríguez.
# Date: 2021/01/21.
# Description: Standardized methodology to perform supervised machine learning 
# on a dataset
# Usage: 
#===============================================================================

# PACKAGES
library(data.table) # load dataset
library(RColorBrewer) # color palette
library(Rtsne) # data visualization
library(caret) # machine learning
library(limma) # Differential Methylation Analysis
library(randomForest) # Recursive Feature Elimination
library(doParallel) # paralellize

# GLOBAL VARIABLES
# Path where the data of the project is stored
DATA_DIR <- "/home/user14/Nextcloud/share2/share/20201229_TFM/data/"
# Path where the results of the project are stored
RESULTS_DIR <- "/home/user14/Nextcloud/share2/share/20201229_TFM/results/"
# Path where the code of the project is stored
CODE_DIR <- "/home/user14/Nextcloud/share2/share/20201229_TFM/code/"
# CSV file with samples information, it MUST have a column with tissue site 
# information for supervised classification
SAMPLESHEET <- "samplesheet_ines.csv"
# Name of the column with tissue site information
SAMPLE_GROUP <- "Sample_Group"
# CSV file with b-values
BVALUES <- "all_samples.csv"
seed <- 1234
# A numeric value of frecuency of missing values to filter out probes
ProbeCutoff <- 0.1
# ML model evaluation metrics
metric <- "Accuracy"
# accuracy and kappa: default metrics on binary and multi-class datasets (kappa 
# is more suitable for unbalanced data).
# RMSE and R^2: default metrics on regression datasets.
# AUROC: only for binary datasets
# LogLoss: useful in multi-class datasets
# For Cross-Validation
# Select number of repeats for cross-validation
r <- 3
# Select number of folds for cross-validation
k <- 10
# Select number of features for limma selection
n <- 500

# FUNCTION
source(file = file.path(CODE_DIR, "limma_fs.R"))

################################################################################### 
##STEP 1: LOAD DATA
###################################################################################
# print("1. DATA LOADING")
# # for testing
dataset <- fread(file.path(DATA_DIR, "subset_30_5.csv"))
dataset <- as.data.frame(dataset)
dataset[, SAMPLE_GROUP] <- as.factor(dataset[, SAMPLE_GROUP])

################################################################################### 
##STEP 2: DATA EXPLORATION
###################################################################################
print("2. DATA EXPLORATION")

# Colour palette and shape for visualiazing 32 tumoral types
pal <- rep(brewer.pal(8, "Dark2"), 4)
pch <- c(rep(15,8), rep(16,8), rep(17, 8), rep(18,8))

# t-SNE = t-Distributed Stochastic Neighbor Embedding (of 50000 random CpGs)
# ---------------------------------------------------
# is a non-linear technique for dimensionality reduction that is particularly well 
# suited for the visualization of high-dimensional datasets.
print("Performing t-Distributed Stochastic Neighbor Embedding (t-SNE)...")
set.seed(1234)
dataset_na <- dataset[ , apply(dataset, 2, function(x) !any(is.na(x)))]
dataset_tsne_subset <- sample(dataset_na[-length(dataset_na)], 50000)
tsne <- Rtsne(as.matrix(dataset_tsne_subset), 
              dims = 2, 
              perplexity=6, 
              verbose=TRUE, 
              max_iter = 5000,
              check_duplicates = FALSE)
png(file = file.path(RESULTS_DIR, paste0("tsne_labels.png")),
    res = 300,
    height = 2700,
    width = 2700)
plot(tsne$Y, 
     type = "n", 
     cex = 1 , 
     main = paste0("t-sne"), 
     col = pal[factor(dataset_na[, SAMPLE_GROUP])],
     xlab = "First dimension", 
     ylab = "Second dimension")
text(tsne$Y, 
     labels = rownames(dataset),
     cex = 0.5, 
     col = pal[factor(dataset_na[, SAMPLE_GROUP])])
legend("topright", legend=levels(factor(dataset_na[, SAMPLE_GROUP])), text.col=pal, col = pal,
       bg="white", pch = pch, cex=0.5)
dev.off()

################################################################################
## STEP 3: PREPARE DATA
################################################################################
# Split-out validation dataset
set.seed(seed)
trainIndex <- createDataPartition(dataset[, SAMPLE_GROUP], p=0.70, list=FALSE)
dataTrain <- dataset[ trainIndex,]
dataTest <- dataset[-trainIndex,]

################################################################################
## STEP 4: FEATURE SELECTION
################################################################################
# Remove probes with more than 10% of missing values
keep <- mclapply(dataTrain, function(x) {
    percentage_NA <-sum(is.na(x)) / dim(dataTrain)[1]
    return(ifelse(percentage_NA > ProbeCutoff, FALSE, TRUE))
})
dataTrain <- dataTrain[, unlist(keep)]
dataTest <- dataTest[, unlist(keep)]

# Impute NAs
preProcessImpute <- preProcess(dataTrain, method = "medianImpute")
dataTrain <- predict(preProcessImpute, dataTrain)
dataTest <- predict(preProcessImpute, dataTest)

# DMPs
set.seed(seed)
system.time(keep <- limma_fs(dataTrain[, -length(dataTrain)], factor(dataTrain[, length(dataTrain)]), n))
# create a new dataset with DMPs
keep <- c(keep, SAMPLE_GROUP)
print(paste("Features:", length(keep)))
dataTrain <- dataTrain[, keep]
dataTest <- dataTest[, keep]
dataTrain <- dataTrain[, - grep("[.][1-9]$", colnames(dataTrain))]
dataTest <- dataTest[, - grep("[.][1-9]$", colnames(dataTest))]

# # Rank Features by Importance
# # ensure results are repeatable
# set.seed(seed)
# # prepare training scheme
# control <- trainControl(method="repeatedcv", number=10, repeats=3)
# # train the model
# system.time(model <- train(dataTrain[,-length(dataTrain)], 
#                as.factor(dataTrain[,length(dataTrain)]), 
#                method="lvq",
#                trControl=control,
#                tuneGrid = data.frame(size = 3, k = 1:2)))
# # estimate variable importance
# importance <- varImp(model, scale=FALSE)
# # summarize importance
# print(importance, top = 50)
# # plot importance
# png(file = file.path(RESULTS_DIR, paste0("TRIAL_2_RI_50.png")),
#     res = 300,
#     height = 3500,
#     width = 3000)
# plot(importance, top = 50)
# dev.off()
# # list the chosen features
# keep <- rownames(importance$importance)[1:50]
# # create a new dataset with Most Important Features
# keep <- c(keep, SAMPLE_GROUP)
# print(paste("Features:", length(keep)))
# dataTrain <- dataTrain[, keep]
# dataTest <- dataTest[, keep]
# 
# Wrapped methods: Recursive Feature Selection
set.seed(seed)
# define the control using a random forest selection function
control <- rfeControl(functions=rfFuncs, method="cv", number=10, repeats = 5)
# define a vector with the number of features that should be trained with
subset <- seq(1, ncol(dataTrain), by = 1)
# run the RFE algorithm
system.time(results <- rfe(dataTrain[,-length(dataTrain)],
               as.factor(dataTrain[,length(dataTrain)]),
               size = subset,
               metric = "Kappa",
               rfeControl = control))
# summarize the results
print(results)
# list the chosen features
predictors <- predictors(results)
# # plot the results
# png(file = file.path(RESULTS_DIR, paste0("TRIAL_6_rfe.png")),
#     res = 300,
#     height = 1500,
#     width = 3000)
# plot(results, type=c("g", "o"))
# dev.off()
keep <- colnames(dataTrain) %in% predictors
keep[length(keep)] <- TRUE
print(paste("Features:", sum(keep)))
dataTrain <- dataTrain[, keep]
dataTest <- dataTest[, keep]

# Handle class imbalance
model_weights <- NULL
for (m in dataTrain[, SAMPLE_GROUP]){
    model_weights <- c(model_weights, (1/table(dataTrain[, SAMPLE_GROUP])[m]) * 0.5)
}

################################################################################
## STEP 5: EVALUATE ALGORITHMS
################################################################################
# 10-fold cross validation with 3 repeats
trainControl <- trainControl(method="repeatedcv", number=k, repeats=r)
set.seed(seed)
# Train the model with the training data
# LDA
fit.lda <- train(as.data.frame(dataTrain[, -length(dataTrain)]), 
                 as.factor(dataTrain[, length(dataTrain)]),
                 method = "lda",
                 metric = metric,
                 weights = model_weights,
                 trControl = trainControl)
print(fit.lda)
# estimate skill of LDA on the validation dataset
predictions <- predict(fit.lda, dataTest)
confusionMatrix(predictions, dataTest[, SAMPLE_GROUP])
# Naive Bayes
fit.nb <- train(as.data.frame(dataTrain[, -length(dataTrain)]), 
                as.factor(dataTrain[, length(dataTrain)]),
                method = "nb",
                metric = metric,
                weights = model_weights,
                trControl = trainControl)
print(fit.nb)
# estimate skill of NB on the validation dataset
predictions <- predict(fit.nb, dataTest)
confusionMatrix(predictions, dataTest[, SAMPLE_GROUP])
# SVM
fit.svm <- train(as.data.frame(dataTrain[, -length(dataTrain)]), 
                 as.factor(dataTrain[, length(dataTrain)]),
                 method = "svmRadial",
                 metric = metric,
                 weights = model_weights,
                 trControl = trainControl)
print(fit.svm)
# estimate skill of SVM on the validation dataset
predictions <- predict(fit.svm, dataTest)
confusionMatrix(predictions, dataTest[, SAMPLE_GROUP])
# CART (built-in feature selection)
fit.cart <- train(as.data.frame(dataTrain[, -length(dataTrain)]), 
                  as.factor(dataTrain[, length(dataTrain)]),
                  method = "rpart",
                  metric = metric,
                  weights = model_weights,
                  trControl = trainControl)
print(fit.cart)
# estimate skill of CART on the validation dataset
predictions <- predict(fit.cart, dataTest)
confusionMatrix(predictions, dataTest[, SAMPLE_GROUP])
# KNN
fit.knn <- train(as.data.frame(dataTrain[, -length(dataTrain)]), 
                 as.factor(dataTrain[, length(dataTrain)]),
                 method = "knn",
                 metric = metric,
                 weights = model_weights,
                 trControl = trainControl)
print(fit.knn)
# estimate skill of KNN on the validation dataset
predictions <- predict(fit.knn, dataTest)
confusionMatrix(predictions, dataTest[, SAMPLE_GROUP])
# RF (built-in feature selection)
fit.rf <- train(as.data.frame(dataTrain[, -length(dataTrain)]), 
                as.factor(dataTrain[, length(dataTrain)]),
                method = "rf",
                metric = metric,
                weights = model_weights,
                trControl = trainControl)
print(fit.rf)
# estimate skill of RF on the validation dataset
predictions <- predict(fit.rf, dataTest)
confusionMatrix(predictions, dataTest[, SAMPLE_GROUP])
# C5.0 (built-in feature selection)
fit.c5.0 <- train(as.data.frame(dataTrain[, -length(dataTrain)]), 
                  as.factor(dataTrain[, length(dataTrain)]),
                  method = "C5.0",
                  metric = metric,
                  weights = model_weights,
                  trControl = trainControl)
print(fit.c5.0)
# estimate skill of C5.0 on the validation dataset
predictions <- predict(fit.c5.0, dataTest)
confusionMatrix(predictions, dataTest[, SAMPLE_GROUP])
# TREEBAG
fit.treebag <- train(as.data.frame(dataTrain[, -length(dataTrain)]), 
                     as.factor(dataTrain[, length(dataTrain)]),
                     method = "treebag",
                     metric = metric,
                     weights = model_weights,
                     trControl = trainControl)
print(fit.treebag)
# estimate skill of TREEBAG on the validation dataset
predictions <- predict(fit.treebag, dataTest)
confusionMatrix(predictions, dataTest[, SAMPLE_GROUP])
# GBM (built-in feature selection)
fit.gbm <- train(as.data.frame(dataTrain[, -length(dataTrain)]), 
                 as.factor(dataTrain[, length(dataTrain)]),
                 method = "gbm",
                 metric = metric,
                 weights = model_weights,
                 trControl = trainControl)
print(fit.gbm)
# estimate skill of GBM on the validation dataset
predictions <- predict(fit.gbm, dataTest)
confusionMatrix(predictions, dataTest[, SAMPLE_GROUP])
# Compare algorithms
results <- resamples(list(LDA=fit.lda, NB=fit.nb, SVM=fit.svm,
                          CART=fit.cart, KNN=fit.knn, RF=fit.rf, 
                          C5.0=fit.c5.0, TREEBAG=fit.treebag, GBM=fit.gbm))
summary(results)
png(file = file.path(RESULTS_DIR, paste0("TRIAL_6_dotplot_fs_limma_rf.png")),
    res = 300,
    height = 2700,
    width = 2700)
dotplot(results)
dev.off()

################################################################################
## STEP 6: IMPROVE ACCURACY
################################################################################
# Tunning LDA
# no tunning parametres

# Tunning RF
# Random Search
trainControl <- trainControl(method="repeatedcv", number=10, repeats=3, search="random")
set.seed(seed)
mtry <- sqrt(ncol(as.data.frame(dataTrain[, -length(dataTrain)])))
rfRandom <- train(as.data.frame(dataTrain[, -length(dataTrain)]), 
                  as.factor(dataTrain[, length(dataTrain)]), 
                  method="rf", 
                  metric=metric, 
                  tuneLength=15,
                  trControl=trainControl)
print(rfRandom)
png(file = file.path(RESULTS_DIR, paste0("TUNE_6_limma_rfe_rf_random.png")),
    res = 300,
    height = 2700,
    width = 2700)
plot(rfRandom)
dev.off()

# Grid Search
trainControl <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
set.seed(seed)
tunegrid <- expand.grid(.mtry=c(1:15))
rfGrid <- train(as.data.frame(dataTrain[, -length(dataTrain)]), 
                as.factor(dataTrain[, length(dataTrain)]),
                method="rf", 
                metric=metric, 
                tuneGrid=tunegrid,
                trControl=trainControl)
print(rfGrid)
png(file = file.path(RESULTS_DIR, paste0("TUNE_6_limma_rfe_rf_grid.png")),
    res = 300,
    height = 2700,
    width = 2700)
plot(rfGrid)
dev.off()

# Tunning GBM
fitControl <- trainControl(method = "repeatedcv", number = 10, repeats = 10)
gbmGrid <-  expand.grid(interaction.depth = (1:5) * 2, 
                        n.trees = (1:10)*25, 
                        shrinkage = 0.1,
                        n.minobsinnode = 20)
set.seed(seed)
gbmFit <- train(as.data.frame(dataTrain[, -length(dataTrain)]), 
                 as.factor(dataTrain[, length(dataTrain)]), 
                 method = "gbm", 
                 trControl = fitControl, 
                 verbose = FALSE, 
                 tuneGrid = gbmGrid)
gbmFit
png(file = file.path(RESULTS_DIR, paste0("TUNE_6_limma_rfe_gbm_1_accuracy.png")),
    res = 300,
    height = 2700,
    width = 2700)
trellis.par.set(caretTheme())
plot(gbmFit)
dev.off()

png(file = file.path(RESULTS_DIR, paste0("TUNE_6_limma_rfe_gbm_1_kappa.png")),
    res = 300,
    height = 2700,
    width = 2700)
trellis.par.set(caretTheme())
plot(gbmFit, metric = "Kappa")
dev.off()

set.seed(seed)
gbmGrid2 <-  expand.grid(interaction.depth = c(1, 5, 9), 
                        n.trees = (1:30)*50, 
                        shrinkage = 0.1,
                        n.minobsinnode = 20)
gbmFit2 <- train(as.data.frame(dataTrain[, -length(dataTrain)]), 
                 as.factor(dataTrain[, length(dataTrain)]), 
                 method = "gbm", 
                 trControl = fitControl, 
                 verbose = FALSE, 
                 tuneGrid = gbmGrid2)
gbmFit2
png(file = file.path(RESULTS_DIR, paste0("TUNE_6_limma_rfe_gbm_2_accuracy.png")),
    res = 300,
    height = 2700,
    width = 2700)
trellis.par.set(caretTheme())
plot(gbmFit2)
dev.off()

png(file = file.path(RESULTS_DIR, paste0("TUNE_6_limma_rfe_gbm_2_kappa.png")),
    res = 300,
    height = 2700,
    width = 2700)
trellis.par.set(caretTheme())
plot(gbmFit2, metric = "Kappa")
dev.off()

################################################################################
## STEP 7: FINALIZE MODEL
################################################################################
finalModel <- train(as.data.frame(dataTrain[, -length(dataTrain)]), 
                    as.factor(dataTrain[, length(dataTrain)]),
                    method = "lda",
                    metric = metric,
                    weights = model_weights)
finalPredictions <- predict(finalModel, dataTest)
confusionMatrix(finalPredictions, dataTest[, SAMPLE_GROUP])
saveRDS(finalModel, "./finalModel.rds")

# Save selected features
selected_features <- colnames(dataTrain[-length(dataTrain)])
saveRDS(selected_features, "./selected_features.RData")

superModel <- readRDS("./finalModel.rds")
print(superModel)
finalPredictions <- predict(superModel, dataTest[, -length(dataTest)])
print(finalPredictions)
confusionMatrix(finalPredictions, dataTest[, SAMPLE_GROUP])
