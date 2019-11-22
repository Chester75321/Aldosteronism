##############################
#Created on Feb 19 2019

#@author: Chester
##############################

##############################
# install & import libraries
##############################
library(devtools)
library(rhdf5)
library(sleuth)

##############################
# get argument
##############################
str_fileName_metadata <- "../Week_16/metadata.txt"
str_filePath_kallisto <- "../Week_16/Kallisto/"
str_filePath_output <- "../Week_16/"

##############################
# main function
##############################
# Get file paths
list_filePath_kallisto <- list.dirs(str_filePath_kallisto, recursive=FALSE)
# Get sample names
list_sampleName <- gsub("/", "", gsub(str_filePath_kallisto, "", list_filePath_kallisto))
dt_filePath <- data.frame(sample=list_sampleName, path=list_filePath_kallisto)

# Read metadata
dt_meta <- read.table(file.path(str_fileName_metadata), header=TRUE, stringsAsFactors=FALSE)
dt_meta <- merge.data.frame(dt_meta, dt_filePath)
dt_meta$path <- as.vector(dt_meta$path)

# Initial sleuth object
so <- sleuth_prep(dt_meta, ~condition, extra_bootstrap_summary=TRUE, read_bootstrap_tpm=TRUE)
# Fit the full model
so <- sleuth_fit(so)
# Fit the reduce model
so <- sleuth_fit(so, ~1, 'reduced')
# Run a likelihood ratio test (LRT) between the two models
so <- sleuth_lrt(so, 'reduced', 'full')

# Run Wald test to get beta
model_so <- models(so)
str_intercept <- colnames(model_so$full$design_matrix)[2]
so <- sleuth_wt(so, str_intercept, 'full')

# Get result
dt_result <- sleuth_results(so, 'reduced:full', 'lrt', show_all=TRUE)
dt_result_beta <- dplyr::select(sleuth_results(so, str_intercept, test_type='wt', show_all=TRUE), target_id, b)
dt_result <- dplyr::inner_join(dt_result, dt_result_beta, by='target_id')
dt_result <- dplyr::arrange(dt_result, qval)
# Filter result
dt_result_significant <- dplyr::filter(dt_result, qval <= 0.05)
dt_result_significant <- dplyr::arrange(dt_result_significant, qval)

# Output result table
write.table(dt_result, file.path(str_filePath_output, "Sleuth_Result.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
# Plot significant transcripts
if (dim(dt_result_significant)[1] > 0){
    for (idx in 1:dim(dt_result_significant)[1]){
        png(file.path(str_filePath_output, paste0("Sleuth_Result_", dt_result$target_id[idx], ".png")))
        plot(plot_bootstrap(so, dt_result$target_id[idx], units="tpm", color_by="condition"))
        dev.off()
    }
}
