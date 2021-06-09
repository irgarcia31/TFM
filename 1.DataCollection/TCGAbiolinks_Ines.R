# ==================================================================================

# TCGA Biolinks is able to accesss de NCI to SEARCH, DOWNLOAD and PREPARE relevant 
# data for analysis in R
# Author: Inés Rodríguez García
# Date: 2020/09/30
# WARNING: it is necessary to have version 4.0 of R, otherwise dependecy "purrrogress" 
# won't be available.
# Usage: For downloading specific projects use the script step-by-step.
# 	For downloading all TCGA IDAT files, run the EXAMPLE at the end of this script.

# ==================================================================================
# If necessary install packages -> uncomment
# If not -> comment
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")

## PACKAGES
library(TCGAbiolinks)
library(tidyr)
library(dplyr)

## GLOBAL VARIABLES

## FUNCTIONS

# =================================================================================
###################################################################################
## SEARCHING GDC DATABASE
###################################################################################
query <- GDCquery(project = c(""), #
                  data.category = "Raw microarray data",
                  data.type = "Raw intensities",
                  legacy = TRUE,
                  platform = c("Illumina Human Methylation 450"),
                  file.type="idat",
                  experimental.strategy = "Methylation array"
)

# http://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/query.html
# extended fields information
# 2 available ources to download GDC data:
# - GDC legacy archive: unmodifies copy of the data (hg19, hg18)
# - GDC Harmonize Database: data harmonized against hg38

###################################################################################
## DOWNLOADING AND PREPARING FILES FOR ANALYSIS
###################################################################################
tryCatch(GDCdownload(query, method = "api", files.per.chunk = 10), error = function(e) GDCdownload(query, method = "client"))
data <- GDCprepare(query)

# 2 methods to download GDC data:
# - client: more reliable but slower -> small datasets
# - api: the data is compressed into a tar.gz file
#   If the size and the number of the files are too big this tar.gz will be too big 
#   which might have a high probability of download failure.
#   To solve that -> files.per.chunk argument which will split the files into small chunks.
#   If no files.per.chunk is specified but the files are too big, it will select
#   the number of chunks by default.

# With tryCatch() if an error occurs downloading by api's method, it changes to 
# client's method.

match.file.cases <- getResults(query,cols=c("cases","file_name"))
# This will create a map between idat file name, cases (barcode) and project
readr::write_tsv(match.file.cases, path =  "idat_filename_case.txt")
# code to move all files to local folder
for(file in dir(".",pattern = ".idat", recursive = T)){
  TCGAbiolinks:::move(file,basename(file))
}

###################################################################################
## GET CLINICAL DATA AND GENERATE SAMPLESHEET
###################################################################################
projects <- TCGAbiolinks:::getGDCprojects()$project_id
projects <- projects[grepl('^TCGA',projects,perl=T)]
for (proj in projects){
    query <- GDCquery(project = proj,
                    data.category = "Raw microarray data",
                    data.type = "Raw intensities", 
                    experimental.strategy = "Methylation array", 
                    legacy = TRUE,
                    file.type = ".idat",
                    platform = "Illumina Human Methylation 450")
  query.clinical <- GDCquery(project = proj, 
                  data.category = "Clinical", 
                  file.type = "xml")
  GDCdownload(query.clinical)
  clinical <- GDCprepare_clinic(query.clinical, clinical.info = "patient")
  if ("tumor_tissue_site" %in% colnames(clinical)) {
    clinical_sub <- clinical[, c("bcr_patient_barcode", "tumor_tissue_site", "histological_type", "gender", "age_at_initial_pathologic_diagnosis", "race_list")]
  } else if ("primary_pathology_tumor_tissue_site" %in% colnames(clinical)) {
    clinical_sub <- clinical[, c("bcr_patient_barcode", "primary_pathology_tumor_tissue_site", "primary_pathology_histological_type", "gender", "primary_pathology_age_at_initial_pathologic_diagnosis", "race_list")]
  } else {
    print("ERROR: undefined columns selected")
  }
  samplesheet <- getResults(query,cols=c("cases.submitter_id","file_name","platform","project","center_short_name","sample_type"))
  colnames(samplesheet)[1] <- c("sample_name")
  colnames(clinical_sub)[1] <- c("sample_name")
  final.samplesheet <- merge(samplesheet, clinical_sub, by="sample_name")
  final.samplesheet$file_name <- gsub("*_....idat", "", final.samplesheet$file_name)
  final.samplesheet <- final.samplesheet[!(duplicated(final.samplesheet)), ]
  final.samplesheet <- separate(data = final.samplesheet, col = file_name, into = c("Sentrix_ID", "Sentrix_Position"), sep = "_")
  names(final.samplesheet) <- c("sample_name", "Sentrix_ID", "Sentrix_Position", "platform", "project", "center_short_name", "sample_type", "tumor_tissue_site", "histological_type", "gender", "age_at_initial_pathologic_diagnosis", "race_list")
  readr::write_csv(final.samplesheet,path = paste0("samplesheet_", proj, ".csv"))
}

# all_clin$tissue_source_site ???

# all_clin$days_to_birth -> Number of days between the date used for index and the date 
# from a person's date of birth represented as a calculated negative number of days.
# all_clin$gender -> Text designations that identify gender. Gender is described as 
# the assemblage of properties that distinguish people on the basis of their societal 
# roles. [Explanatory Comment 1: Identification of gender is based upon self-report 
# and may come from a form, questionnaire, interview, etc.]
# (female, male, unknown, unspecified, not reported)

# all_clin$age_at_initial_pathologic_diagnosis

# all_clin$ethnicity -> An individual's self-described social and cultural grouping, 
# specifically whether an individual describes themselves as Hispanic or Latino. 
# The provided values are based on the categories defined by the U.S. Office of 
# Management and Business and used by the U.S. Census Bureau.
# (hispanic or latino, not hispanic or latino, Unknown, not reported, not allowed to collect)

# all_clin$race_list -> An arbitrary classification of a taxonomic group that is 
# a division of a species. It usually arises as a consequence of geographical 
# isolation within a species and is characterized by shared heredity, physical 
# attributes and behavior, and in the case of humans, by common history, nationality, 
# or geographic distribution. The provided values are based on the categories defined 
# by the U.S. Office of Management and Business and used by the U.S. Census Bureau.
# (white, american indian or alaska native, black or african american, asian, 
# native hawaiian or other pacific islander, other, Unknown, not reported, 
# not allowed to collect)

# all_clin$tobacco_smoking_history -> Category describing current smoking status and
# smoking history as self-reported by a patient.
# (1, 2, 3, 4, 5, 6, 7, Unknown, Not Reported)

# all_clin$number_pack_years_smoked -> Numeric computed value to represent lifetime 
# tobacco exposure defined as number of cigarettes smoked per day x number of years 
# smoked divided by 20.

# all_clin$maximum_tumor_dimension

# For further information check the following link:
# https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-entity-list&anchor=clinical

###################################################################################
## EXAMPLE: DNA methylation: Get all TCGA IDAT files
###################################################################################
# Recomendation: execute it in your data directory
projects <- TCGAbiolinks:::getGDCprojects()$project_id
projects <- projects[grepl('^TCGA',projects,perl=T)]
match.file.cases.all <- NULL
for(proj in projects){
  print(proj)
  query <- GDCquery(project = proj,
                    data.category = "Raw microarray data",
                    data.type = "Raw intensities", 
                    experimental.strategy = "Methylation array", 
                    legacy = TRUE,
                    file.type = ".idat",
                    platform = "Illumina Human Methylation 450")
  match.file.cases <- getResults(query,cols=c("cases","file_name"))
  match.file.cases$project <- proj
  match.file.cases.all <- rbind(match.file.cases.all,match.file.cases)
  tryCatch(GDCdownload(query, method = "api", files.per.chunk = 20),
           error = function(e) GDCdownload(query, method = "client"))
}
# This will create a map between idat file name, cases (barcode) and project
readr::write_tsv(match.file.cases.all, path =  "idat_filename_case.txt")
# code to move all files to local folder
for(file in dir(".",pattern = ".idat", recursive = T)){
  TCGAbiolinks:::move(file,basename(file))
}

# =================================================================================
