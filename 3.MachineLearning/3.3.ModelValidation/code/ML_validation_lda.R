#===============================================================================
# Machine Learning Validation Workflow
# Author: Inés Rodríguez
# Date: 2021/02/01
# Description: Classification of new samples using a Machine Learning model
# previously generated.
# Usage:
#===============================================================================

# PACKAGES
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(data.table)
library(caret)

# GLOBAL VARIABLES
# Raw data folder, with the samplesheet and the idat files
IDAT_DIR <- "/home/user14/Nextcloud/share2/share/20201229_TFM/data/Validation_set_2/Validation_set_2_arrays/"
# Path where the data of the project is stored
DATA_DIR <- "/home/user14/Nextcloud/share2/share/20201229_TFM/data/Validation_set_2/"
# Path where the temporary files and final b-values should be stored
BVALS_DIR <- "/home/user14/Nextcloud/share2/share/20201229_TFM/data/Validation_set_2/Validation_set_2_b_vals/"
# Path where the model is stored
MODEL_DIR <- "/home/user14/Nextcloud/share2/share/20201229_TFM/david_runs/20210308_test_ML_no_cv_limma_rfe_joint_classes_XGB"
# Path where the results of the project are stored
RESULTS_DIR <- "/home/user14/Nextcloud/share2/share/20201229_TFM/results/"
# Path where the code of the project is stored
CODE_DIR <- "/home/user14/Nextcloud/share2/share/20201229_TFM/code/"
# Sample sheet file name.
# It requires sample names in a column named "Sample_Name" and the following
# columns: "Sample_Plate", "Sentrix_ID", "Sentrix_Position".
SAMPLESHEET <- "Validation_set_2.csv"
# csv files containing probes to exclude as they are cross-reactive.
POLYMORPHIC450k <- "/home/user14/Nextcloud/share2/share/documentation/illumina450k_filtering-master/48639-non-specific-probes-Illumina450k.csv"
POLYMORPHICEPIC <- "/home/user14/Nextcloud/share2/share/documentation/illumina450k_filtering-master/EPIC/13059_2016_1066_MOESM1_ESM.csv"
FAILED2020 <- "/home/user14/Nextcloud/share2/share/20201229_TFM/code/EPIC_MFG_Change_List_of_Affected_Probes_091419.csv"
# R object containing the final selected features used by the ML model
SELECTED_FEATURES <- "predictors.csv"
# R object containing the final Machine Learning Model
ML_MODEL <- "trained_algorithms_CV.RData"
# Quality filters:
# Remove poor quality samples based on proportion of failed probes.
max_failed_probes <- 0.1
# Threshold to consider failed probe.
probe_detP_max <- 0.01
# Whether to make failed probes NA.
mask_failed_probes <- TRUE
# If your data includes males and females, remove probes on the sex chromosomes
remove_sex_chr <- TRUE
# Remove probes with SNPs at CpG site or at a single nucleotide extension.
remove_SNP_probes <- TRUE
# Exclude cross reactive probes
remove_cross_reactive <- TRUE
# Load preloaded model arrays to be use in 450k-EPIC array combinations.
load("/home/user14/Nextcloud/share2/share/20201229_TFM/code/model_450k_EPIC.RData")
# For the few parallel computations.
CORES <- 1

#===============================================================================
### FUNCTIONS.

#' Drop CpGs by SNP annotation
#'
#' This is a custom version of minfi::dropLociWithSnps intented to be used
#' with a beta-values or m-values matrix as input, instead of a
#' GenomicRatioSet object, despite of it also accepts GenomicRatioSet objects,
#' in which case returns behaves as minfi::dropLociWithSnps and returns a
#' filtered GenomicRatioSet
#'
#' @param a_matrix a beta-values or m-values matrix. A GenomicRatioSet object
#'                 is also possible, in which case this function behaves as
#'                 minfi::dropLociWithSnps
#' @param snpAnno a DataFrame object generated with minfi::getSnpInfo function
#'                containing the SNP annotation of the `a_matrix` probes
#' @param snps a character vector with the filtering to be used.
#'             Default: c("CpG", "SBE")
#' @param maf a numeric value of minor allele frequency to filter out probes.
#'            Default: 0
#'
#' @details Three filter criteria can be used: "Probe", "CpG" and "SBE"
#'          correspond the SNPs present inside the probe body, at the CpG
#'          interrogation and at the single nucleotide extension respectively
#'
#' @return The same matrix as `a_matrix` but with the SNP probes filtered out.
#'         If a_matrix is a GenomicRatioSet object, this same object class is
#'         returned
dropLociWithSnps_custom <- function(a_matrix,
                                    snpAnno,
                                    snps = c("CpG", "SBE"),
                                    maf = 0) {
  res <- a_matrix
  to_remove <- character(0)
  # Drop "Probe" SNPs
  if ("Probe" %in% snps) {
    to_remove <- c(to_remove,
                   rownames(snpAnno)[!(is.na(snpAnno$Probe_rs)) &
                                       snpAnno$Probe_maf > maf])
  }
  
  # Drop "CpG" SNPs
  if ("CpG" %in% snps) {
    to_remove <- c(to_remove,
                   rownames(snpAnno)[!(is.na(snpAnno$CpG_rs)) &
                                       snpAnno$CpG_maf > maf])
  }
  
  # Drop "SBE" SNPs
  if ("SBE" %in% snps) {
    to_remove <- c(to_remove,
                   rownames(snpAnno)[!(is.na(snpAnno$SBE_rs)) &
                                       snpAnno$SBE_maf > maf])
  }
  to_keep <- !(rownames(res) %in% to_remove)
  res <- res[to_keep, , drop = FALSE]
  return(res)
}

#' Load samples by samplesheet into a 450k-like virtual array
#'
#' This function loads 450k and/or EPIC arrays from a sample sheet using minfi
#' internal functions with the added improvement that it will combine samples
#' from 450k and EPIC array into a 450k-like virtual array. It's also
#' implemented to perform operations in parallel, to improve performance.
#'
#' @param ss A data.frame containing the sample_sheet. Mandatory columns are:
#'   'Sample_Name' & 'Basename'. Very importantly, 'Basename'
#'   should contain the full path to each idat array file, but neither the final
#'   file extension nor the channel extension ('_Red.idat' or '_Grn.idat').
#' @param annot An array annotation object. E.g.
#'   annot <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
#' @param virtual_450k_model an \code{\link[minfi]{RGChannelSet}} of a virtual
#'   450k array already obtained by \code{\link[minfi]{combineArrays}} using
#'   a 450k model and an the same EPIC model expected in \code{model_EPIC}.
#'   It's used to produce always the same \code{\link[minfi]{combineArrays}}
#'   results.
#' @param out_dir path to directory where each sample b-values will be
#'   stored. Important: be sure that the 'out_dir' is empty and will only
#'   contain the output of this function, in other to merge them correctly.
#'   Default: dir.create(file.path(getwd(), 'combined_arrays')).
#' @param cores number of cores to be used in parallel computations. Default:
#'   1L.
#' @param max_failed_probes max proportion of failed probes to keep a sample.
#'   Default: 1L (no sample discarded).
#' @param mask_failed_probes whether to mask ('NA') failed probes or not.
#'   Default: TRUE.
#' @param remove_cpgs a character vector with the CpG ids to be removed (e.g.
#'   polymorphic probes). Default: NULL.
#' @param remove_sex whether to remove probes in sex chromosomes.
#'   Default: FALSE.
#' @param remove_snp whether to remove proves based on a known SNP in the
#'   interrogated CpG or at 1 bp. Default = TRUE.
#' @param probe_detP_max the max detection P value to consider a probe not
#'   failed. Default: 0.01.
#'
#' @details This function uses the "Pool_ID" column from the sample sheet which
#'   should specify either "450k" or "Epic", to select samples. Using this info,
#'   it will load 450k and EPIC samples separatedly and then combine them with a
#'   model array from the other type into a 450k-like arrays. Then, both
#'   450k-like arrays will be merged into a single
#'   \code{\link[minfi]{RGChannelSet}} object, discarding the model array. In
#'   order to mantain probes constant, in the case of only 450k (or only EPIC)
#'   samples, they will also be merged with the model array of the other type
#'   to generate the 450k-like array.
#'
#' @return This function writes all the samples in a separate file, consisting
#'  of a single column with all the beta-values and using the sample_name and
#'  the array basename as colname. In addition, also returns a named bool
#'  vector with TRUE / FALSE depending on the sample passed QC or not.
#'
load.combined.450k.EPIC <- function(
  ss,
  annot,
  virtual_450k_model,
  out_dir = dir.create(file.path(getwd(), 'combined_arrays')),
  cores = 1L,
  max_failed_probes = 1L,
  mask_failed_probes = TRUE,
  remove_cpgs = NULL,
  remove_sex = FALSE,
  remove_snp = TRUE,
  probe_detP_max = 0.01) {
  sample_processed <- mclapply(1:nrow(ss), function(i) {
    sn <- as.character(ss[i, "Sample_Name"])
    bc <- as.character(basename(as.character(ss[i, "Basename"])))
    print(paste0("[MSG] Processing sample ", sn, " with barcode: ", bc))
    rgSet <- minfi::read.metharray.exp(targets = ss[i, ],
                                       force = TRUE)
    # Combine with model array.
    combined <- minfi::combineArrays(rgSet, virtual_450k_model)
    combined <- combined[, 1]
    mSet <- minfi::preprocessNoob(combined, dyeMethod = "single")
    b_vals <- minfi::getBeta(mSet)
    detP <- minfi::detectionP(combined)
    stopifnot(all(rownames(b_vals) == rownames(detP)))
    if (discard.sample.qc(detP, max_failed_probes, probe_detP_max)) {
      print(paste0("[MSG] Sample ", sn, " with barcode: ", bc, " QC rejected!."))
      return(FALSE)
    }
    if (mask_failed_probes) {
      b_vals <- NA.failed.probes(b_vals, detP, probe_detP_max)
    }
    if (!is.null(remove_cpgs)) {
      stopifnot(is.character(remove_cpgs))
      keep <- !(rownames(b_vals) %in% remove_cpgs)
      b_vals <- b_vals[keep, , drop = FALSE]
    }
    if (remove_sex_chr) {
      keep <- !(rownames(b_vals) %in% annot$Name[annot$chr %in% c("chrX","chrY")])
      b_vals <- b_vals[keep, , drop = FALSE]
    }
    if (remove_snp) {
      snp_info <- minfi::getSnpInfo(mSet)
      snp_info <- snp_info[rownames(b_vals), ]
      b_vals <- dropLociWithSnps_custom(b_vals, snpAnno = snp_info,
                                        snps = c("CpG", "SBE"), maf = 0)
    }
    
    full_sample_name <- paste(as.character(ss[i, "Sample_Name"]),
                              basename(as.character(ss[i, "Basename"])),
                              sep = "_")
    total_n_samples <- dim(ss)[1]
    target_id_filename <- paste(
      formatC(0, width = nchar(as.character(total_n_samples)),
              format = "d", flag = "0"),
      "aaa_target_id.csv", sep = "_")
    sample_filename <- paste0(
      formatC(i, width = nchar(as.character(total_n_samples)),
              format = "d", flag = "0"),
      "_", full_sample_name, ".csv")
    
    colnames(b_vals) <- full_sample_name
    target_id <- as.data.frame(rownames(b_vals))
    colnames(target_id) <- "Target_ID"
    
    if (!file.exists(file.path(out_dir, target_id_filename))){
      write.table(target_id,
                  file = file.path(out_dir, target_id_filename),
                  quote = FALSE,
                  sep = ",",
                  row.names = FALSE,
                  col.names = TRUE
      )
    }
    write.table(b_vals,
                file = file.path(out_dir, sample_filename),
                quote = FALSE,
                sep = ",",
                row.names = FALSE,
                col.names = TRUE
    )
    print(paste0("[MSG] Sample ", sn, " with barcode: ", bc, " COMPLETED!"))
    return(TRUE)
  }, mc.cores = cores)
  names(sample_processed) <- ss$Sample_Name
  return(unlist(sample_processed))
}


#' Discard samples by max failed probe proportion
#'
#' @description Determines whether a sample should be discarded by proportion
#' of failed probes or not.
#'
#' @param detP a matrix with the detection P values.
#' @param max_failed_probes the max proportion (between 0 and 1) of failed
#'   probes allowed.
#' @param probe_detP_max the max detection P value to consider a probe not
#'   failed. Default: 0.01.
#' @return a logical vector to use as logical mask to filter out samples.
#'
discard.sample.qc <- function(detP, max_failed_probes, probe_detP_max = 0.01) {
  keep <- apply(detP, 2, function(x) {
    # All probes detP == 0 means totally failed array.
    if (sum(x) == 0) {
      return(FALSE)
    } else {
      failed_probes <- sum(x >= probe_detP_max)
      failed_sample_prop <- failed_probes / length(x)
      return(ifelse(failed_sample_prop < max_failed_probes, TRUE, FALSE))
    }
  })
  return(!keep)
}

#' Mask failed probes
#'
#' @description Mask (change to `NA`) beta values based on detection P value.
#'
#' @param b_vals a beta values matrix.
#' @param detP a detection P values matrix.
#' @param probe_detP_max the max detection P value to consider a probe not
#'   failed. Default: 0.01.
#'
#' @return a b_vals matrix with al the failed probes masked by NAs.
#'
NA.failed.probes <- function(b_vals, detP, probe_detP_max = 0.01){
  # Mark all probes that have failed (detection pval > threshold) as NAs.
  failed_cpgs <- which(detP > probe_detP_max, arr.ind = TRUE)
  b_vals[failed_cpgs] <- NA
  return(b_vals)
}
#===============================================================================
### MAIN.
################################################################################
## Step 1: data preprocessing
################################################################################
# Load Illumina annotation
annot <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# Determine probes to exclude
if (remove_cross_reactive) {
  poly_450k <- read.csv(POLYMORPHIC450k, stringsAsFactors = FALSE)
  poly_EPIC <- read.csv(POLYMORPHICEPIC, stringsAsFactors = FALSE)
  failed2020 <- read.csv(FAILED2020, stringsAsFactors = FALSE)
  probes_to_exclude <- unique(c(poly_450k[,1], poly_EPIC[,1], failed2020[, 1]))
} else {
  probes_to_exclude <- NULL
}
# Load sample sheet
ss <- read.csv(file.path(DATA_DIR, SAMPLESHEET),
               skip = 7,
               stringsAsFactors = FALSE)
# Generate/correct Basename column, as it's the one required by minfi.
ss$Basename <- file.path(IDAT_DIR,
                         #ss$Sample_Plate,
                         #ss$Sentrix_ID,
                         paste(ss$Sentrix_ID,
                               ss$Sentrix_Position,
                               sep = "_")
)
samples_processed <- load.combined.450k.EPIC(
  ss = ss,
  annot = annot,
  virtual_450k_model = virtual_450k_model,
  out_dir = BVALS_DIR,
  cores = CORES,
  max_failed_probes = max_failed_probes,
  mask_failed_probes = mask_failed_probes,
  remove_cpgs = probes_to_exclude,
  remove_sex = remove_sex_chr,
  remove_snp = remove_SNP_probes,
  probe_detP_max = probe_detP_max)

print("[MSG] Processing completed!")
print(paste0("[MSG] Original samples:", dim(ss)[1]))
print(paste0("[MSG] Total of samples passing QC: ", sum(samples_processed)))

print("[MSG] Joining all the files together using linux paste")
system2("paste", args = c("-d", ",",
                          file.path(BVALS_DIR, "*.csv"),
                          ">",
                          file.path(BVALS_DIR, "all_samples.csv"))
)
print("[MSG] Joined completed!")

################################################################################
## Step 2: data loading
################################################################################
# load preprocessed validation dataset
validation <- fread(file.path(BVALS_DIR, "all_samples.csv"))
validation <- as.data.frame(validation)
rownames(validation) <- validation$Target_ID
validation$Target_ID <- NULL

################################################################################
## Step 3     : predictions
################################################################################
# load the model
superModel <- load(file.path(MODEL_DIR, ML_MODEL))
print(fit.lda)

# select features of the model
predictors <- as.character(predictors(fit.lda))
validation <- as.data.frame(validation[predictors, , FALSE])

# Impute NAs
preProcessImpute <- preProcess(validation, method = "medianImpute")
validation <- predict(preProcessImpute, validation)

# transpose the dataset so samples are in rows and features in columns
validation <- as.data.frame(t(validation))

# make predictions
finalPredictions <- predict(fit.lda, newdata = validation, type = "prob")
print(finalPredictions)

# export CSV file
write.csv(finalPredictions, file = file.path(RESULTS_DIR, "final_prediction_lda_set_2.csv"))
