#===============================================================================
# Machine Learning Workflow
# Author: Inés Rodríguez.
# Date: 2021/01/21.
# Description: Standardized methodology to perform supervised machine learning 
# on a dataset with limma as Feature Selection, and LDA, RF and GBM as ML
# evaluated algorithms
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
DATASET <- "/illumina/runs3/dpineyro/students/2020_Ines_Rodriguez/Work/Epicup_remake/20201217_all_dataset_bvals/result/20201214_betas/all_samples.csv"
# Path where the results of the project are stored
RESULTS_DIR <- "/illumina/runs3/dpineyro/students/2020_Ines_Rodriguez/Work/Epicup_remake/david_tests/20210212_test_ML_no_cv_limma_rfe"
# Path where the code of the project is stored
CODE_DIR <- "/illumina/runs3/dpineyro/students/2020_Ines_Rodriguez/Work/Epicup_remake/scripts"
# Sample sheet
SAMPLE_SHEET <- "/illumina/runs3/dpineyro/students/2020_Ines_Rodriguez/Work/Epicup_remake/20201217_all_dataset_bvals/data/samplesheet_ines.csv"
# Name of the column with tissue site information
SAMPLE_GROUP <- "Sample_Group"
# seed 
seed <- 1234
# A numeric value of frecuency of missing values to filter out probes
ProbeCutoff <- 0.1
# Pearson's correlation threshold above which probes are considered to be higly 
# correlated
CorrelatedCutoff <- 0.9
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
n <- 9600
# Multiprocessing cores
CORES <- 48

# FUNCTION
# Imported the function caused a segmentation fault in Dori (I don't know why...).
# I copied the function in the same script.
#source(file = file.path(CODE_DIR, "limma_fs.R"))

#' Impute NAs using the column median
#'
#' Giving two data.frames, this function will compute the column median values
#' and then it will impute both data.frame NA values with the corresponding
#' column median of the first data.frame.
#'
#' @param data_train : data.frame. A data.frame from which the column medians
#'   will be calculated. It is normally a data.frame with samples as rows and
#'   features as columns.
#' @param data_test : data.frame. A second data.frame to impute NA values.
#'
#' @return A list with two elements:
#'   data_train_imp : data.frame. This is the same `data_train` data.frame
#'     but with all the NA values imputed using its column median values.
#'   data_test_imp : data.frame. This is the same `data_test` data.frame
#'     but with all the NA values imputed using the column median values of
#'     `data_train`.
#'
median.imputer <- function(data_train, data_test) {
    stopifnot(dim(data_train)[2] == dim(data_test)[2])
    medians_train <- apply(data_train, 2, median, na.rm = TRUE)
    d_train_imp <- lapply(1:dim(data_train)[2], function(x) {
        imp_col <- data_train[, x]
        imp_col[is.na(imp_col)] <- medians_train[x]
        return(imp_col)
    })
    d_train_imp <- as.data.frame(d_train_imp,
                                 col.names = colnames(data_train))
    d_test_imp <- lapply(1:dim(data_test)[2], function(x) {
        imp_col <- data_test[, x]
        imp_col[is.na(imp_col)] <- medians_train[x]
        return(imp_col)
    })
    d_test_imp <- as.data.frame(d_test_imp,
                                col.names = colnames(data_test))
    return(list(data_train_imp = d_train_imp,
                data_test_imp = d_test_imp))
}

#' Select CpGs by differential methylation using limma
#'
#' This funtion uses limma to select top differentially methylated probes
#' (features) between groups (sample labels factor). The number of class labels
#' can be 2 or more, but no paired designs are implemented. It is equivalent to
#' anova_fs, but much faster. Only applicable to methylation array data.
#'
#' @param input_data [data.frame] a \code{data.frame} rows represent samples and
#'   columns represent features (beta values or methylation ratios).
#' @param sample_labels [factor] a \code{factor} of length == nrow(input_data),
#'   containing the sample labels.
#' @param n [numeric] the number of features (top differentially methylated CpGs
#'   to select).
#' @param remove_dup [bool] Whether to return a unique set of features or not.
#'   This is important when there are more than 2 classes and any feature is
#'   repeatedly selected for more than one class. Default: TRUE.
#'
#' @return [character] A character vector with the \code{n} top differentially
#'   methylated (sorted by p-value) feature names, for all contrasts tested.
#'   For more than two groups in \code{sample_labels}, each contrast is designed
#'   as one Vs the others.
#'
#' @examples
#' set.seed(1234)
#' input_data <- data.frame(a = runif(10, min = 0, max = 1),
#'                          b = runif(10, min = 0, max = 1),
#'                          c = runif(10, min = 0, max = 1),
#'                          d = runif(10, min = 0, max = 1))
#' sample_labels <- as.factor(c(rep("classA", 5), rep("classB", 5)))
#' features_selected <- limma_fs(input_data, sample_labels, 2)
#' print(features_selected)
#'
#' @export
limma_fs <- function(input_data, sample_labels, n, remove_dupl = TRUE) {
  print("[MSG limma_fs]")
  # First, make sure sample_labels has valid names
  sample_labels <- factor(make.names(sample_labels))
  # First, check whether n is < total number of features.
  #if (dim(input_data)[2] <= n) {
  #  warning("The number of features to select is not smaller than the total
  #          number of features. No features were selected.")
  #  return(colnames(input_data))
  #} else if (sum(input_data == 0) | sum(input_data == 1)) {
  #  # Check whether absolute 0 or 1 values exists. So it is usually a sign of
  #  # sequencing data.
  #  warning("0 and/or 1 are found in your data. Is this data generated by
  #          bisulfite sequencing?. Limma filtering is not applied.")
  #  return(colnames(input_data))
  #} else {
    # Calculate M-values using M=log2(Beta/(1-Beta)). All statistics will be
    # performed on M-values. Transpose to have rows as features and columns as
    # samples.
    print("[MSG limma_fs] Starting limma FS")
    sample_names <- rownames(input_data)
    input_data <- t(input_data)
    m_vals <- log2(input_data / (1-input_data))
    # Create a targets dataframe.
    pheno_data <- data.frame(sample_names, sample_labels)
    rownames(pheno_data) <- sample_names
    targets <- stats::model.frame(sample_names ~ sample_labels, pheno_data)
    # Design matrix (only unpaired test supported).
    design <- stats::model.matrix(~0+sample_labels, data=targets)
    colnames(design) <- levels(sample_labels)
    # Contrast matrix (one vs the others).
    print("[MSG limma_fs] Creating contrast matrix")
    if (nlevels(sample_labels) == 2) {
      contr <- paste0(levels(sample_labels)[1], "-", levels(sample_labels)[2])
      contMatrix <- limma::makeContrasts(contrasts = contr, levels = design)
    } else {
      # More than 2 groups. Using One Vs the Others contrasts.
      i <- 1
      contr <- character()
      while (i <= nlevels(sample_labels)) {
        one <- levels(sample_labels)[i]
        the_others <- levels(sample_labels)[-i]
        # We will make contrasts such as "A-(B+C)/2".
        the_others_mean <- paste0("(", paste(the_others, collapse = "+"),
                                  ")/", length(the_others))
        contr <- c(contr, paste0(one, "-", the_others_mean))
        i <- i + 1
      }
      contMatrix <- limma::makeContrasts(contrasts = contr, levels = design)
    }
    print("[MSG limma_fs] Fitting limma")
    # fit the linear model
    fit <- limma::lmFit(m_vals, design)
    # fit the contrasts
    fit2 <- limma::contrasts.fit(fit, contMatrix)
    print("[MSG limma_fs] eBayes")
    fit2 <- limma::eBayes(fit2)
    # Toptable of all contrasts, selected by p-value.
    if (nlevels(sample_labels) < 3) {
      list_of_results <- lapply(1:(nlevels(sample_labels) - 1), function(x) {
        limma::topTable(fit2, coef = x, number=Inf, sort.by = "P")
      })
    } else {
      list_of_results <- lapply(1:(nlevels(sample_labels)), function(x) {
        limma::topTable(fit2, coef = x, number=Inf, sort.by = "P")
      })
    }
    # Selecting the "best" CpGs of each contrast. If
    # n/nlevels(sample_labels) has no integer result, round
    # approximation is taken.
    if (nlevels(sample_labels) < 3) {
      each_contrast_n <- round(n/(nlevels(sample_labels) - 1), 0)
    } else {
      each_contrast_n <- round(n/(nlevels(sample_labels)), 0)
    }
    feature_selection <- character()
    for (r in list_of_results) {
      feature_selection <- c(feature_selection, rownames(r)[1:each_contrast_n])
    }
    if (remove_dupl){
      return(unique(feature_selection))
    } else {
      return(feature_selection)
    }
  #}
}

################################################################################### 
##STEP 1: LOAD DATA
###################################################################################
# print("1. DATA LOADING")
# # for testing
dataset <- fread(DATASET, data.table = FALSE)
rownames(dataset) <- dataset[, 1]
dataset <- dataset[, -1]
colnames(dataset) <-  gsub("_[\\dRC]+", "",  colnames(dataset), perl = TRUE)

ss <- read.table(SAMPLE_SHEET, header = TRUE,
                      sep = ",", skip = 7)
# Order ss as our dataset.
ss <- ss[match(ss$Sample_Name, colnames(dataset)), ]
stopifnot(all(ss$Sample_Name == colnames(dataset)))
targets <- factor(ss[, SAMPLE_GROUP])

# Transpose dataset
dataset <- t(dataset)  # Samples as rows, features as columns.

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
#print("Performing t-Distributed Stochastic Neighbor Embedding (t-SNE)...")
#set.seed(seed)
#dataset_na <- dataset[, apply(dataset, 2, function(x) !any(is.na(x)))]
#tsne_sel <- sample(colnames(dataset_na), 50000)
#dataset_tsne_subset <- dataset_na[, tsne_sel]
#tsne <- Rtsne(dataset_tsne_subset,
#              dims = 2,
#              perplexity=6,
#              verbose=TRUE,
#              max_iter = 5000,
#              num_threads = 8,
#              check_duplicates = FALSE)
#png(file = file.path(RESULTS_DIR, paste0("tsne.png")),
#    res = 300,
#    height = 2700,
#    width = 2700)
#plot(tsne$Y,
#     pch = pch[targets],
#     cex = 1 ,
#     main = paste0("t-sne"),
#     col = pal[targets],
#     xlab = "First dimension",
#     ylab = "Second dimension")
#legend("topright", legend=levels(targets), text.col=pal, col = pal,
#       bg="white", pch = pch, cex=0.5)
#dev.off()

################################################################################
## STEP 3: PREPARE DATA
################################################################################
# Split-out validation dataset
print("Splitting data set into training / test")
set.seed(seed)
trainIndex <- createDataPartition(targets, p=0.70, list=FALSE)
dataTrain <- dataset[ trainIndex,]
targetsTrain <- targets[trainIndex]
dataTest <- dataset[-trainIndex,]
targetsTest <- targets[-trainIndex]

################################################################################
## STEP 4: FEATURE SELECTION
################################################################################
cl <- makePSOCKcluster(CORES)
registerDoParallel(cl)

print("Performing FS:")
print("...removing features by missing value percentage")

# Remove probes with more than 10% of missing values
nas <- apply(dataTrain, 2, function(x) sum(is.na(x)))
percent_NA <- nas / dim(dataTrain)[1]
keep <- ifelse(percent_NA > ProbeCutoff, FALSE, TRUE)
dataTrain <- dataTrain[, keep]
dataTest <- dataTest[, keep]
print("...impute NA values")
imputation <- median.imputer(dataTrain, dataTest)
dataTrain <- imputation[['data_train_imp']]
dataTest <- imputation[['data_test_imp']]
#preProcessImpute <- preProcess(dataTrain, method = "medianImpute")
#dataTrain <- predict(preProcessImpute, dataTrain)
#dataTest <- predict(preProcessImpute, dataTest)

print("...selecting features using limma")
# Filter method: DMPs
keep <- limma_fs(dataTrain, targetsTrain, n)
# create a new dataset with DMPs
#keep <- c(keep, SAMPLE_GROUP)
print(paste("Features:", length(keep)))
dataTrain <- dataTrain[, keep]
dataTest <- dataTest[, keep]

print("...removing highly correlated features")
# find attributes that are highly correlated and remove them
set.seed(seed)
correlations <- cor(dataTrain)
highlyCorrelated <- findCorrelation(correlations, cutoff = CorrelatedCutoff)
print(paste("Features:", dim(dataTrain)[2] - length(highlyCorrelated)))
# create a new dataset without highly correlated features
dataTrain <- dataTrain[,-highlyCorrelated]
dataTest <- dataTest[,-highlyCorrelated]

print("...selecting features using RFE")
# Wrapped methods: Recursive Feature Selection
set.seed(seed)
# define the control using a random forest selection function
control <- rfeControl(functions = rfFuncs, method = "cv", number = 10, repeats = 5)
# define a vector with the number of features that should be trained with
subset <- 1:10
for (n in 1:ncol(dataTrain)) {
  m <- subset[length(subset)] + n
  subset <- c(subset, m)
  if (subset[length(subset)] > ncol(dataTrain)){
    break
  }
}
# run the RFE algorithm
results <- rfe(as.data.frame(dataTrain), 
               targetsTrain,
               size = subset,
               metric = "Kappa",
               rfeControl = control)
# summarize the results
print(results)
# list the chosen features
predictors <- predictors(results)
# plot the results
png(file = file.path(RESULTS_DIR, paste0("TRIAL_7_dotplot_fs_limma_rfe.png")),
    res = 300,
    height = 1500,
    width = 3000)
plot(results, type=c("g", "o"))
dev.off()
keep <- colnames(dataTrain) %in% predictors
keep[length(keep)] <- TRUE
print(paste("Features:", sum(keep)))
dataTrain <- dataTrain[, keep]
dataTest <- dataTest[, keep]

# Handle class imbalance
model_weights <- NULL
for (m in targetsTrain){
    model_weights <- c(model_weights, (1/table(targetsTrain)[m]))
}
print("FS completed")
################################################################################
## STEP 5: EVALUATE ALGORITHMS
################################################################################
# 10-fold cross validation with 3 repeats
trainControl <- trainControl(method="repeatedcv", number=k, repeats=r)
set.seed(seed)
# Train the model with the training data
# LDA
print("Training LDA")
fit.lda <- train(as.data.frame(dataTrain), 
                 targetsTrain,
                 method = "lda",
                 metric = metric,
                 weights = model_weights,
                 trControl = trainControl)
print(fit.lda)
# estimate skill of LDA on the validation dataset
predictions <- predict(fit.lda, dataTest)
confusionMatrix(predictions, targetsTest)

# RF (built-in feature selection)
print("Training RF")
fit.rf <- train(as.data.frame(dataTrain), 
                targetsTrain,
                method = "rf",
                metric = metric,
                weights = model_weights,
                trControl = trainControl)
print(fit.rf)
# estimate skill of RF on the validation dataset
predictions <- predict(fit.rf, dataTest)
confusionMatrix(predictions, targetsTest)

# GBM (built-in feature selection)
print("Training GBM")
fit.gbm <- train(as.data.frame(dataTrain), 
                 targetsTrain,
                 method = "gbm",
                 metric = metric,
                 weights = model_weights,
                 trControl = trainControl)
print(fit.gbm)
# estimate skill of GBM on the validation dataset
predictions <- predict(fit.gbm, dataTest)
confusionMatrix(predictions, targetsTest)

# Compare algorithms
results <- resamples(list(LDA=fit.lda, RF=fit.rf, GBM=fit.gbm))
summary(results)
png(file = file.path(RESULTS_DIR, paste0("TRIAL_7_dotplot_fs_limma_rfe.png")),
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
rfRandom <- train(as.data.frame(dataTrain), 
                  targetsTrain, 
                  method="rf", 
                  metric=metric, 
                  tuneLength=15,
                  weights = model_weights,
                  trControl=trainControl)
print(rfRandom)

png(file = file.path(RESULTS_DIR, paste0("TRIAL_7_limma_rfe_rf_random.png")),
    res = 300,
    height = 2700,
    width = 2700)
plot(rfRandom)
dev.off()

# Tunning GBM
fitControl <- trainControl(method = "repeatedcv", number = 10, repeats = 10)
gbmGrid <-  expand.grid(interaction.depth = (1:5) * 2, 
                        n.trees = (1:10)*25, 
                        shrinkage = 0.1,
                        n.minobsinnode = 20)
set.seed(seed)
gbmFit <- train(as.data.frame(dataTrain), 
                targetsTrain, 
                method = "gbm", 
                trControl = fitControl, 
                weights = model_weights,
                verbose = FALSE, 
                tuneGrid = gbmGrid)
gbmFit

png(file = file.path(RESULTS_DIR, paste0("TRIAL_7_limma_rfe_gbm_accuracy.png")),
    res = 300,
    height = 2700,
    width = 2700)
trellis.par.set(caretTheme())
plot(gbmFit)
dev.off()

png(file = file.path(RESULTS_DIR, paste0("TRIAL_7_limma_rfe_gbm_kappa.png")),
    res = 300,
    height = 2700,
    width = 2700)
trellis.par.set(caretTheme())
plot(gbmFit, metric = "Kappa")
dev.off()

stopCluster(cl)
