# ==================================================================================

# Metlylation array analysis for DEA (Data Exploration Analysis)
# Author: Inés Rodríguez García
# Date: 2020/09/21
# Description: Standardized methodology to analyze methylation arrays (450k and
#   EPIC) and perform DEA (Data Exploration Analysis)
# Usage: Change variables at GLOBAL VARIABLES section and execute the script.
#   [WARNING] Be careful with variable_of_interest, they'll depend on your 
#   samplesheet columns.

# ==================================================================================
##PACKAGES
library(tidyr)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylationEPICmanifest)
library(RColorBrewer)
library(limma)
library(Rtsne)
library(umap)
library(factoextra)

library(stringr)
arg <- commandArgs(trailingOnly = TRUE)
print(paste(arg))
abbreviation <- substring(str_extract(arg, "\\-.+\\."), 2, nchar(str_extract(arg, "\\-.+\\.")) -1)
print(paste(abbreviation))
dir.create(file.path("/mnt/ElRaid/irodriguez/Nextcloud/share/20201005_IRodriguez_all_TCGA/results", abbreviation))

##GLOBAL VARIABLES
# Path where the data of the project is stored
DATA_DIR <- file.path("/mnt/ElRaid/irodriguez/Nextcloud/share/20201005_IRodriguez_all_TCGA")
# Path where the results are supossed to be saved
RESULTS_DIR <- file.path(paste0("/mnt/ElRaid/irodriguez/Nextcloud/share/20201005_IRodriguez_all_TCGA/results/", abbreviation)) # ir creando carpetas
print(RESULTS_DIR)
# CSV file containing one line per sample , with a number of columns describing each sample
SAMPLESHEET <- arg
# Path to the cross-reactive probes CSV
CROSSREACTIVE450K <- file.path("/mnt/ElRaid/irodriguez/Nextcloud/share/20200930_IRodriguez_TCGA_CESC_PAAN/data/48639-non-specific-probes-Illumina450k.csv")
CROSSREACTIVEEPIC <- file.path("/mnt/ElRaid/irodriguez/Nextcloud/share/documentation/illumina450k_filtering-master/EPIC/13059_2016_1066_MOESM1_ESM.csv")
# Sample sheet column name to use as in differential methylation tests.
class_label <- "tumor_tissue_site"
# Sample sheet column that identify individuals (or any other paired info), for paired tests. 
# If not paired test, let "".
# individual_label <- "Cell_line"
individual_label <- ""
# Contrast to perform (to use in makeContrasts).
contr <- c(
  ""
)
# Pearson's correlation threshold above which a sample is considered a possible
# duplicated sample, taking into account array genotyping snps.
duplicated_snp_thr <- 0.9
# Threshold to consider failed probe.
probe_detP_max <- 0.01
# Remove poor quality samples based on mean detection pval.
sample_mean_detP_max <- 0.01
# Remove poor quality samples based on proportion of failed probes.
max_failed_probes <- 0.05
# Decide whether to remove probes failed in any sample or to mark them as NAs
#remove_or_NA_probes <- "remove"
remove_or_NA_probes <- "NA"
# If your data includes normal tissue, remove normal samples
remove_normal <- TRUE
normal_tissue <- "Solid Tissue Normal"
# If your data includes recurrent tumor, remove recurrent tumor samples
remove_recurrent <- TRUE
recurrent_tumor <- "Recurrent Tumor"
# If your data includes additional - new primary, remove those samples
remove_additional <- TRUE
additional_newprimary <- "Additional - New Primary"
# If your data includes metastasis, remove metastatic samples
# If your intention is to use this data for ML, it would be interesting to train the 
# classifier with and without metastatic samples
remove_metastatic = TRUE # If FALSE be careful with duplicates
metastatic <- "Metastatic"
# If your data includes males and females, remove probes on the sex chromosomes
remove_sex_chr <- TRUE
# Remove probes with SNPs at CpG site or at a single nucleotide extension.
remove_SNP_probes <- TRUE
# Exclude cross reactive probes 
remove_cross_reactive <- TRUE
# Normalization method
# Recommendations:
#   Quantile: use it when you do not expect global methylation differences 
#     between samples, for instance: when all samples are from the same tissue.
#   Funnorm: use it when you expect global methylation differences between  
#     samples, for instance: tumor/normal or different tissues.
#   Illumina: use it when you want to reproduce GenomeStudio results.
#   Noob: use it to allow further sample adding. It is also a good default
#     normalization method.
norm_method <- "Noob"
# Quality filters:
# Perform sex prediction (TRUE or FALSE).
# CAUTION!! It requires a column named "Gender" in the sample sheet. The
# gender code should be "Male/Female" or "M/F". Missing values should be noted
# as NA.
check_sex <- TRUE
# Specify for which variables you want to perform the Data Exploration c("project", "tumor_tissue_site","histological_type", "gender", "age_dichotomic", "race_list")
variables_of_interest <- c("tumor_tissue_site","histological_type", "gender", "age_dichotomic", "race_list")

##FUNCTIONS

# DropLociWithSnps customize so it can use bValue/mValue matrix as input, instead of
# GenomicRatioSet object (although it can be used as input too)

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
  res <- res[to_keep,]
  return(res)
}

#==================================================================================
###################################################################################
## STEP 1: DATA LOADING
###################################################################################
print("1. DATA LOADING")
print("Reading in the samplesheet of the experiment...")
# Read in the sample sheet for the experiment
targets <- read.metharray.sheet(DATA_DIR, pattern=SAMPLESHEET)
# This generates a data.frame with one row for each sample and several columns
# The function uses the specified path and other information from the sample sheet 
# to create a column called Basename which specifies the location of each individual 
# IDAT file in the experiment.
dim(targets)[1]
print(paste0("Number of samples:", dim(targets)[1]))

# If necessary delete Normal Samples
print("Removing Normal Samples if present...")
n <- dim(targets)[1]
if (remove_normal == TRUE) {
  targets <- targets[targets$sample_type != normal_tissue, ]
}
print(paste0(n -  dim(targets)[1], " Normal Samples removed. Total number of samples: ", dim(targets)[1]))

# If necessary delete Recurrent Tumor Samples
print("Removing Recurrent Tumor Samples if present...")
n <- dim(targets)[1]
if (remove_recurrent == TRUE) {
  targets <- targets[targets$sample_type != recurrent_tumor, ]
}
print(paste0(n -  dim(targets)[1], " Recurrent Tumor samples removed. Total number of samples: ", dim(targets)[1]))

# If necessary delete Additional - New Primary
n <- dim(targets)[1]
print("Removing Aditional - New Primary Samples if present...")
if (remove_additional == TRUE) {
  targets <- targets[targets$sample_type != additional_newprimary, ]
}
print(paste0(n -  dim(targets)[1], " Aditional - New Primary Samples removed. Total number of samples: ", dim(targets)[1]))

# If necessary delete Metastatic Samples
n <- dim(targets)[1]
print("Removing Metastatic Samples if present...")
if (remove_metastatic == TRUE) {
  targets <- targets[targets$sample_type != metastatic, ]
}
print(paste0(n -  dim(targets)[1], " Metastatic Samples removed. Total number of samples: ", dim(targets)[1]))

# Read in the raw data from de IDAT files
print("\nReading in the raw data from the IDAT files...")
rgSet <- read.metharray.exp(targets=targets)
# This generates a S4 Object of class RGChannelSet that contains all the raw intensity 
# data, from both the red and green colour channels, for each of the samples

sampleNames(rgSet) <- targets$sample_name

# Detection of Illumina platform
print("Detecting Illumina platform...")
platform_array <- substring(rgSet@annotation[1], 25, nchar(rgSet@annotation[1]))

# Get annotation data and probes to exclude
if (platform_array == "450k") {
  ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  print("The array platform is 450k")
} else if (platform_array == "EPIC") {
  ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  print("The array platform is EPIC")
} else {
  print("ERROR: Problem detecting array platform")
}

# If necessary, modify some variables for the Data Exploration
targets$age_dichotomic[targets$age_at_initial_pathologic_diagnosis < 65] <- "Young"
targets$age_dichotomic[targets$age_at_initial_pathologic_diagnosis >= 65] <- "Old"

###################################################################################
##STEP 2: QUALITY CONTROL
###################################################################################
print("2. QUALITY CONTROL")
# Genotyping SNPs QC to check sample origin. Useful for experiments with more than 
# one array per patient.
print("Genotyping SNPs QC to check sample origin...")
snps <- getSnpBeta(rgSet) 
# getSnpBeta: Retrives the measurements of the 65 SNP probes located on the array. 
# These SNP probes are intended to be used for sample tracking and sample mixups. 
# The return value is a matrix of beta values. Each SNP probe ought to have values 
# clustered around 3 distinct values corresponding to homo-, and hetero-zygotes.
snps_cor <- cor(snps)
possible_duplicates <- which(snps_cor > duplicated_snp_thr & snps_cor < 1, arr.ind = TRUE)
# which() Give the TRUE indices of a logical object, allowing for array indices.
# Check whether any correlation is > than duplicated_snp_thr(0.9), sugesting possible duplicated
if (dim(possible_duplicates)[1] == 0) {
  print(paste0("No duplicated samples detected by SNPs analysis (Pearson's correlation threshold at ", duplicated_snp_thr, ")"))
} else {
  print(paste0("[WARNING] ", (dim(possible_duplicates)[1]) , " Possible duplicated samples detected (Pearson's correlation threshold at " , duplicated_snp_thr, ". Remove them manually!!)"))
  print(paste0(rownames(possible_duplicates)))
  write.csv(snps_cor[rownames(possible_duplicates), 
                     rownames(possible_duplicates)], 
            file = file.path(RESULTS_DIR, paste0("possible_duplicated_samples_", abbreviation ,".csv")))
}
# paste0() = paste(..., sep="", collapse) Concatenate vectors after converting to character.

# Calculate the detection p-values (for every CpG in every sample, which is 
# indicative of the quality of the signal)
print("Calculating detection p-values...")
detP <- detectionP(rgSet)
# Compares the total signal (M+U) for each probe to the background signal level, 
# which is estimated from the negative control probes. 
# Very small p-values are indicative of a reliable signal whilst large p-values, 
# for example >0.01, generally indicate a poor quality signal.

# Examine mean detection p-values across all samples to identify any failed samples
pal <- brewer.pal(8,"Dark2")
png(file = file.path(RESULTS_DIR, paste0("detection_pvals_", abbreviation, ".png")),
    height = 1000,
    width = 1000)
barplot(colMeans(detP), col=pal[factor(targets[ ,class_label])], 
        las=2, 
        cex.names=0.8, 
        ylim = c(0, 0.08),
        ylab="Mean detection p-values")
abline(h=0.05,col="red")
abline(h= 0.01, col = "blue")
legend("topright", legend=levels(factor(targets[ ,class_label])), fill=pal,
       bg="white")
dev.off()

# Examine the proportion of failed probes (detP >= probe_detP_max)
print(paste0("Examining proportion of failed probes (detP >=", probe_detP_max, ")..."))
failed_proportion <- apply(detP, 2, function(x) {
  if (sum(x) == 0) {
    print("Failed array") # If summatory of all probes = 0 -> failed array
    return(NA)
  } else {
    failed_probes <- sum(x >= probe_detP_max)
    return(failed_probes/length(x))
  }
})
# apply(matrix, margin(1->rows, 2->columns), function) -> function to manipulate 
# slices of data from matrices, arrays, lists and dataframes in a repetitive way
# In this case no specific function has been defined
# If the summatory of the column is 0 -> NA
# If the summatory of the column is != 0 then calculate how many P-value >= probe_detP_max
# and / between the nº of total p-values in that column
# Now we know the proportion of failed probes in each site
png(file = file.path(RESULTS_DIR, paste0("proportion_failed_probes_", abbreviation, ".png")),
    height = 1000,
    width = 1000)
barplot(failed_proportion, col=pal[factor(targets[, class_label])], 
        ylim = c(0, 0.12),
        las=2, 
        cex.names=0.8, ylab="Proportion of failed probes")
abline(h=0.1,col="red")
abline(h=0.05,col="blue")
abline(h= 0.01, col = "green")
dev.off()

# Generate QC Report
qcReport(rgSet, 
         sampNames=targets$ID, 
         sampGroups=targets[, class_label], 
         pdf=file.path(RESULTS_DIR, paste0("qcReport_", abbreviation, ".pdf")))

# # Remove possible duplicates based on probes quality (from rgSet, targets and detP)
# print("Removing possible duplicates based on probes quality...")
# # Recommendation: do this manually so you can check them before deleting them 
# dup_failed_probes <- failed_proportion[rownames(possible_duplicates)]
# erase <- NULL
# for (i in seq(1, length(dup_failed_probes), by=2)){
#   if (dup_failed_probes[i] >= dup_failed_probes[i+1]) {
#     erase <- append(erase, dup_failed_probes[i])
#   }
# }
# print(paste0(erase))
# for (i in names(erase)) {
#   targets <- targets[targets$sample_name != i, ]
# }
# for (i in names(erase)) {
#   rgSet <- rgSet[ ,rgSet@colData@listData$sample_name != i]
# }
# for (i in names(erase)) {
#   detP <- detP[ , colnames(detP) != i]
# }
# print(paste0("Total number of samples: ", dim(targets)[1]))

print("Removing samples based on probes quality...")
# Remove poor quality samples based on mean detection p-value from rgSet, targets 
# and detection p-value table 
keep <- colMeans(detP) < sample_mean_detP_max
# ColMeans() is equivalent to use of apply with FUN = mean, but are a lot faster.
print(paste0("Samples to remove based on mean detP threshold of ",
             sample_mean_detP_max,
             ": ",
             paste0(targets$sample_name[!keep],";", colnames(detP)[!keep])))
rgSet <- rgSet[,keep]
targets <- targets[keep,]
detP <- detP[,keep]
print(paste0("Total number of samples: ", dim(targets)[1]))

print("Removing samples based on proportion of failed probes...")
# Remove poor quality samples based on proportion of failed probes.
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
# ifelse(..., TRUE, FALSE)
print(paste0("Samples to remove based on proportion of failed probes over ",
             max_failed_probes,
             ": ",
             paste0(targets$Sample_Name[!keep],";", colnames(detP)[!keep])))
rgSet <- rgSet[,keep]
targets <- targets[keep,]
detP <- detP[,keep]
print(paste0("Total number of samples: ", dim(targets)[1]))

# Covert array 450k or EPIC to a pseudo450-EPIC


###################################################################################
##STEP 3: NORMALIZATION
###################################################################################
# To minimise the unwanted variation within and between samples, various data 
# normalisations can be applied. Although there is no single normalisation method 
# that is universally considered best, a recent study by Fortin et al. (2014) has 
# suggested that a good rule of thumb within the minfi framework is that the 
# preprocessFunnorm (Fortin et al. 2014) function is most appropriate for datasets 
# with global methylation differences such as cancer/normal or vastly different 
# tissue types, whilst the preprocessQuantile function (Touleimat and Tost 2012) is 
# more suited for datasets where you do not expect global differences between your 
# samples, for example a single tissue.
# After normalisation, the data is housed in a GenomicRatioSet object. This is a 
# much more compact representation of the data as the colour channel information has
# been discarded and the M and U intensity information has been converted to M-values
# and beta values, together with associated genomic coordinates.
print("3. NORMALIZATION")
print(paste0("Performing normalization by ", norm_method, "..."))
if (norm_method == "Quantile"){
  mSetSq <- preprocessQuantile(rgSet)
} else if (norm_method == "Funnorm") {
  mSetSq <- preprocessFunnorm(rgSet)
} else if (norm_method == "Noob") {
  mSetSq <- preprocessNoob(rgSet, dyeMethod = "single")
}

# visualise what the data looks like before and after normalisation
png(file = file.path(RESULTS_DIR, paste0("betas_before_after_norm_", abbreviation, ".png")),
    height = 1000,
    width = 2000,
    res = 200)
par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=targets[, class_label], main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets[, class_label])), 
       text.col=brewer.pal(8,"Dark2"), cex = 0.3)
densityPlot(getBeta(mSetSq), sampGroups=targets[, class_label],
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets[, class_label])), 
       text.col=brewer.pal(8,"Dark2"), cex = 0.3)
dev.off()

# SEX PREDICTION
print("Predicting sex...")
if (check_sex) {
  ### Sex prediction. Detect sex inconsistencies between predicted and reported.
  sex_idx <- !is.na(targets$gender) # Remove NAs from sex analysis
  gmSetSq <- mapToGenome(mSetSq, mergeManifest = FALSE) # mapToGenome() Map Ilumina
  # methylation array data to the genome using an annotation package
  mSetSq_sex <- gmSetSq[, sex_idx]
  samples_sex <- targets$sample_name[sex_idx]
  sex_prediction <- getSex(mSetSq_sex)
  predicted_sex <- sex_prediction$predictedSex
  reported_sex <- targets$gender[sex_idx]
  # recode to M/F.
  reported_sex <- ifelse(reported_sex == "MALE" | reported_sex == "M", "M", "F")
  inconsistent <- predicted_sex != reported_sex
  if (sum(inconsistent) == 0) {
    print("No gender inconsistencies were found.")
  } else {
    print("Gender inconsistency for samples:")
    print(paste0(samples_sex[inconsistent]))
  }
  # Save sex prediction results as a data.frame.
  df_inconsistencies <- data.frame(list(
    Sample_Name = samples_sex,
    xMed = sex_prediction$xMed,
    yMed = sex_prediction$yMed,
    Reported_gender = reported_sex,
    Predicted_gender = predicted_sex,
    Inconsistency_found = inconsistent
  ))
  write.csv(df_inconsistencies, file = file.path(RESULTS_DIR,
                                                 paste0("sex_predictions_", abbreviation,".csv")))
  # Plot to check.
  plotSex(addSex(mSetSq_sex))
  # Save the plot to check.
  png(file = file.path(RESULTS_DIR, paste0("sex_prediction_", abbreviation, ".png")),
      height = 1200,
      width = 2000,
      res = 160)
  par(mfrow=c(1,2))
  # Plot Reported.
  plot(yMed ~ xMed,
       data = df_inconsistencies,
       pch = "",
       main = "Gender Reported",
       ylab = "Y chr; median total intensity (log2)",
       xlab = "X chr; median total intensity (log2)")
  text(yMed ~ xMed,
       data = df_inconsistencies,
       labels = df_inconsistencies$Sample_Name,
       col = c("Black", "Red")[df_inconsistencies$Reported_gender],
       cex = 0.5)
  legend("topleft",
         legend = levels(factor(df_inconsistencies$Reported_gender)),
         text.col = labels(factor(df_inconsistencies$Reported_gender)))
  # Plot Predicted.
  plot(yMed ~ xMed,
       data = df_inconsistencies,
       pch = "",
       main = "Gender Predicted",
       ylab = "Y chr; median total intensity (log2)",
       xlab = "X chr; median total intensity (log2)")
  text(yMed ~ xMed,
       data = df_inconsistencies,
       labels = df_inconsistencies$Sample_Name,
       col = c("Black", "Red")[df_inconsistencies$Predicted_gender],
       cex = 0.5)
  legend("topleft",
         legend = levels(factor(df_inconsistencies$Predicted_gender)),
         text.col = labels(factor(df_inconsistencies$Predicted_gender)))
  dev.off()
}

###################################################################################
##STEP 4: FILTERING
###################################################################################
print("4. FILTERING")
# Poor performing probes are generally filtered out prior to differential methylation 
# analysis. As the signal from these probes is unreliable, by removing them we perform 
# fewer statistical tests and thus incur a reduced multiple testing penalty. We filter
# out probes that have failed in one or more samples based on detection p-value.

# calculate M-values for statistical analysis
print("Calculating B and M values...")
bVals_original <- getBeta(mSetSq)
mVals_original <- getM(mSetSq)
# For each CpG, there are two measurements: a methylated intensity (denoted by M) and 
# an unmethylated intensity (denoted by U). These intensity values can be used to 
# determine the proportion of methylation at each CpG locus. Methylation levels are 
# commonly reported as either beta values (β=M/(M+U)) or M-values (Mvalue=log2(M/U))
# Beta values are generally preferable for describing the level of methylation at a 
# locus or for graphical presentation because percentage methylation is easily 
# interpretable. However, due to their distributional properties, M-values are more 
# appropriate for statistical testing (Du et al. 2010).

# Ensure probes are in the same order in the bVals and detP objects
detP <- detP[match(rownames(bVals_original),rownames(detP)),]

# Probes that have failed in one or more samples should be treated to get better results.
# There are 2 options:
#   a) Remove failed probes. Advantages: reduced CpGs set therefore fewer statistical 
#       analysis; no "NA".
#   b) Mark as "NA" failed probes. Advantages: keep more CpGs. THIS IS THE CORRECT 
#       OPTION IF YOUR GOING TO USE THE DATA FOR ML.
if (remove_or_NA_probes == "NA") {
  failed_cpgs <- which(detP > probe_detP_max, arr.ind = TRUE)
  bVals <- bVals_original
  bVals[failed_cpgs] <- NA
  mVals <- mVals_original
  mVals[failed_cpgs] <- NA
  print(paste("NAs before mark failed probes:",
              sum(is.na(bVals_original))))
  print(paste("NAs after marking failed probes:",
              sum(is.na(bVals))))
} else if (remove_or_NA_probes == "remove") {
  keep <- rowSums(detP < probe_detP_max) == ncol(mSetSq)
  # if the probe has < probe_detP_max in some sample then it would be excluded and the 
  # summatory will be < ncol(mSetSq)
  bVals <- bVals_original[keep, ]
  mVals <- mVals_original[keep, ]
  print("Probes before removing failed CpGs.")
  print(paste(dim(bVals_original)[1]))
  print("Probes after removing failed CpGs.")
  print(paste(dim(bVals)[1]))
} else {
  stop("ERROR: remove_or_NA_probes must be set to either NA or remove")
}

# if your data includes males and females, remove probes on the sex chromosomes
if (remove_sex_chr) {
  keep <- !(rownames(bVals) %in% ann$Name[ann$chr %in% c("chrX","chrY")])
  # %in% This operator is used to identify if an element belongs to a vector. 
  # If ann$Name$chr is in chrX or chrY then TRUE
  table(keep)
  bVals <- bVals[keep, ]
  mVals <- mVals[keep, ]
  paste0("Probes after removing sex chr CpGs: ", dim(bVals)[1])
}

# remove probes with SNPs at CpG site
if (remove_SNP_probes) {
  # Get snp info for all our probes.
  snp_info <- getSnpInfo(mSetSq)
  # Subset for the ones still present.
  snp_info <- snp_info[rownames(bVals), ]
  bVals <- dropLociWithSnps_custom(bVals, snpAnno = snp_info, 
                                   snps = c("CpG", "SBE"), maf = 0)
  mVals <- dropLociWithSnps_custom(mVals, snpAnno = snp_info, 
                                   snps = c("CpG", "SBE"), maf = 0)
  print("Probes after removing probes with SNPs at CpG site or at a single nucleotide extension.")
  print(paste(dim(bVals)[1]))
}

# exclude cross reactive probes of 450k
if (remove_cross_reactive) {
  xReactiveProbes <- read.csv(CROSSREACTIVE450K, stringsAsFactors=FALSE)
  keep <- !(rownames(bVals) %in% xReactiveProbes[, 1])
  table(keep)
  bVals <- bVals[keep, ]
  mVals <- mVals[keep, ]
  paste0("Probes after removing cross-reactive probes of 450K: ", dim(bVals)[1])
}
# exclude cross reactive probes of EPIC
if (remove_cross_reactive) {
  xReactiveProbes <- read.csv(CROSSREACTIVEEPIC, stringsAsFactors=FALSE)
  keep <- !(rownames(bVals) %in% xReactiveProbes[, 1])
  table(keep)
  bVals <- bVals[keep, ]
  mVals <- mVals[keep, ]
  paste0("Probes after removing cross-reactive probes of EPIC: ", dim(bVals)[1])
}

# AFTER BVALS
# If there is inconsistency between predicted and reported sex, plot bVals of the samples to check quality
print(paste("Samples with sex inconsisency", samples_sex[inconsistent]))
sex_inconsistency <- samples_sex[inconsistent] # add manually the samples
for (i in sex_inconsistency){
  png(file = file.path(RESULTS_DIR, paste0("samples_with_sex_inconsistency_", i , ".png")),
      height = 1000,
      width = 1000,
      res = 200)
  densityPlot(bVals[, i, drop = FALSE],
              main = paste0(i),
              sampGroups = factor(targets[, class_label]), 
              legend = FALSE, 
              xlab="Beta values")
  legend("top", 
         legend = levels(factor(targets[, class_label])), 
         text.col = brewer.pal(8,"Dark2"))
  dev.off()
}

# Write b-vals and m-vals.
# write.csv(bVals, file = file.path(RESULTS_DIR, paste0("b_values_", abbreviation, ".csv")))
# write.csv(mVals, file = file.path(RESULTS_DIR, paste0("m_values_", abbreviation, ".csv")))

# If you are interested in performing supervised ML you will have to add an 
# attribute to the b/mVals table specifying classes to the classifier.
tbVals_df <- as.data.frame(t(bVals))
tbVals_tumor <- cbind(tbVals_df, as.factor(targets[, class_label]))
colnames(tbVals_tumor)[dim(tbVals_tumor)[2]] <- c("tumor_tissue_site")
saveRDS(tbVals_tumor, file = file.path(RESULTS_DIR, paste0("tb_values_tumor_", abbreviation, ".rds")))
saveRDS(bVals, file = file.path(RESULTS_DIR, paste0("b_values_", abbreviation, ".rds")))
saveRDS(mVals, file = file.path(RESULTS_DIR, paste0("m_values_", abbreviation, ".rds")))
saveRDS(targets, file = file.path(RESULTS_DIR, paste0("targets_", abbreviation, ".rds")))
################################################################################### 
##STEP 5: DATA EXPLORATION
###################################################################################
print("5. DATA EXPLORATION")
## Dimensions of the data
## ----------------------
# If you have a lot of instances, you may need to work with a smaller sample of the 
# data so that model training and evaluation is computationally tractable. 
# If you have a vast number of attributes, you may need to select those that are 
# most relevant (Feature Selection). 
# If you have more attributes than instances you may need to select specific 
# modeling techniques.
print(paste0("The dataset has ", dim(bVals)[2], " instances and ", dim(bVals)[1], " atributes."))

## Multi-dimensional scaling (MDS) plots are excellent for visualising data, and 
## -------------------------------
# are usually some of the first plots that should be made when exploring the data. 
# MDS plots are based on principal components analysis and are an unsupervised method 
# for looking at the similarities and differences between the various samples. Samples 
# that are more similar to each other should cluster together, and samples that are 
# very different should be further apart on the plot.
print("Performing multi-dimensional scaling (MDS)...")
# MDS plots to look at largest sources of variation (uncomment for before and after
# filtering)
pal <- brewer.pal(8, "Dark2")
for (i in variables_of_interest) {
  png(file = file.path(RESULTS_DIR, paste0("MDS_", i , ".png")),
      height = 2000,
      width = 2000,
      res = 200)
  plotMDS(mVals, top=1000, gene.selection="common", cex = 0.5,
          main = paste0("MDS_", i),
          col=pal[factor(targets[, i])])
  legend("bottomleft", legend=levels(factor(targets[, i])), text.col=pal,
          bg="white", cex=0.7)
  dev.off()
}

## PCA = Principal Component Analysis is a linear technique for dimensionality reduction.
## ----------------------------------
print("Performing PCA...")
bVals_na <- bVals[complete.cases(bVals), ]
paste0("Probes before removing probes with NA: ", dim(bVals)[1])
paste0("Probes after removing probes with NA: ", dim(bVals_na)[1])
tbVals_na <- t(bVals_na)
pca.out <- prcomp(tbVals_na, center = TRUE, scale=TRUE)
# it is necessary to do the transpose so rows=samples; columns=CpGs
# prcomp() : Performs a principal components analysis on the given data matrix and 
# returns the results as an object of class prcomp.

# Visualize eigenvalues/variances
png(filename = file.path(RESULTS_DIR, paste0("pca_eigenvalues.png")))
fviz_screeplot(pca.out, addlabels = TRUE, ylim = c(0, 50))
dev.off()
pca.out.names <- cbind(targets$sample_name, pca.out$x)
pca.out.names <- pca.out.names[,-1]
cp1 <- pca.out$x [ ,1]
cp2 <- pca.out$x [ ,2]
# cor(cp1, cp2)
cp3 <- pca.out$x [ ,3]
# cor(cp2, cp3)
cp4 <- pca.out$x [ ,4]
# cor(cp3, cp4)

for (i in variables_of_interest) {
  png(filename = file.path(RESULTS_DIR, paste0("PCA_cp1_cp2_", i , ".png")),
      res = 300,
      height = 2000,
      width = 2000)
  plot( cp1, cp2, 
        xlab=paste0("CP1 (", summary(pca.out)$importance[2, 1] * 100 , "% of total variation)"), 
        ylab=paste0("CP2 (", summary(pca.out)$importance[2, 2] * 100, "% of total variation)"), 
        col=pal[factor(targets[, i])], 
        cex = 0.3,
        type = "n")
  text( cp1, cp2, 
        labels = targets$sample_name,
        cex = 0.3, 
        col=pal[factor(targets[, i])])
  legend("topleft", legend=levels(factor(targets[, i])), text.col=pal,
         bg="white", cex=0.3)
  dev.off()
  png(filename = file.path(RESULTS_DIR, paste0("PCA_cp2_cp3_", i , ".png")),
      res = 300,
      height = 2000,
      width = 2000)
  plot( cp2, cp3, 
        xlab=paste0("CP2 (", summary(pca.out)$importance[2, 2] * 100, "% of total variation)"), 
        ylab=paste0("CP3 (", summary(pca.out)$importance[2, 3] * 100, "% of total variation)"), 
        col=pal[factor(targets[, i])], 
        cex = 0.3,
        type = "n" )
  text( cp2, cp3, 
        labels = targets$sample_name, 
        cex = 0.3, 
        col=pal[factor(targets[, i])])
  legend("topleft", legend=levels(factor(targets[, i])), text.col=pal,
         bg="white", cex=0.3)
  dev.off()
}

## t-SNE = t-Distributed Stochastic Neighbor Embedding (of 50000 random CpGs)
## ---------------------------------------------------
# is a non-linear technique for dimensionality reduction that is particularly well 
# suited for the visualization of high-dimensional datasets.
print("Performing t-Distributed Stochastic Neighbor Embedding (t-SNE)...")
set.seed(1234)
bVals_tsne_subset <- sample(which(complete.cases(bVals)), 50000)
tsne <- Rtsne(t(bVals[bVals_tsne_subset, ]), 
              dims = 2, 
              perplexity=6, 
              verbose=TRUE, 
              max_iter = 5000)

for (i in variables_of_interest) {
  png(file = file.path(RESULTS_DIR, paste0("tsne_", i , ".png")),
    res = 300,
    height = 2700,
    width = 2700)
plot(tsne$Y, t='p', type= "n", cex = 0.5 , main = paste0("tsne_", i), col=pal[factor(targets[, i])],
     xlab = "First dimension", ylab = "Second dimension")
text(tsne$Y, 
     labels = targets$sample_name,
     cex = 0.5, 
     col=pal[factor(targets[, i])])
legend("bottomright", legend=levels(factor(targets[, i])), text.col=pal,
       bg="white", cex=0.5)
dev.off()
}

## UMAP = Uniform Manifold Approximation and Projection is an algorithm for dimensional reduction
## ----------------------------------------------------
print("Performing Uniform Manifold Approximation and Projection (UMAP)...")
umap.out <- umap(tbVals_na, random_state=123)
# The UMAP algorithm uses random numbers and thus results may vary from run to run. 
# It is possible to obtain a minimal level of reproducibility by setting seeds (random_state). 
# Add sample names
umap.out.names <- cbind(targets$sample_name, umap.out$layout)
for (i in variables_of_interest) {
png(filename = file.path(RESULTS_DIR, paste0("UMAP_", i , ".png")),
    res = 300,
    height = 2700,
    width = 2700)
plot(umap.out$layout,
     main=paste0("UMAP_", i),
     col=pal[factor(targets[, i])],
     type = "n")
text(umap.out$layout, labels = rownames(umap.out$layout), cex = 0.5, col=pal[factor(targets[, i])])
legend("topleft", legend=levels(factor(targets[, i])), text.col=pal,
       bg="white", cex=0.5)
dev.off()
}

print("THE END!! :)")
