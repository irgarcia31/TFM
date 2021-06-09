#===============================================================================
# Generation of a methylation array CSV samplesheet
# Author: Inés Rodríguez
# Date: 2021/03/15
# Description: 
# Usage: Use this script in order to generate a samplesheet CSV file. The only
# requirement is to dispose all the IDAT files in the same directory.
#===============================================================================

# GLOBAL VARIABLES
IDAT_FILES <- "/home/user14/Nextcloud/share2/share/20201229_TFM/data/Validation_set_2/Validation_set_2_arrays"
RESULTS_DIR <- "/home/user14/Nextcloud/share2/share/20201229_TFM/data/Validation_set_2/Validation_set_2_arrays"
SAMPLESHEET_NAME <- "samplesheet_validation_set_2.csv"

#===============================================================================

# MAIN
Sentrix_ID <- character(0)
Sentrix_Position <- character(0)
for (i in list.files(IDAT_FILES, pattern = "*_Red.idat")) {
  id <- strsplit(i, "_")[[1]]
  Sentrix_ID <- c(Sentrix_ID, id[1])
  Sentrix_Position <- c(Sentrix_Position, id[2])
}
Sample_Name <- as.character(seq(1, length(Sentrix_ID)))
Sample_Plate <- NA

samplesheet_df <- data.frame(Sample_Name, Sample_Plate, Sentrix_ID, Sentrix_Position)
write.csv(samplesheet_df, file = file.path(RESULTS_DIR, SAMPLESHEET_NAME), row.names = FALSE)
