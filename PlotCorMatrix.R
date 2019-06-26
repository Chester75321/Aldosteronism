##############################
#Created on Nov 2018

#@author: Chester
##############################

##############################
# import libraries
##############################
library(corrplot)
library(RColorBrewer)

##############################
# set parameters
##############################
str_inputFilePath <- "../28Mouse/mouse_28_kallisto_20170507/"
list_num_week <- c("00", "01", "04", "10", "16")
list_name_0_week0 <- c("C0_1", "C0_2", "A0_1", "A0_2")
list_name_1_week1 <- c("C1_1", "C1_2", "C1_3", "A1_1", "A1_2", "A1_3")
list_name_2_week4 <- c("C4_1", "C4_2", "C4_3", "A4_1", "A4_2", "A4_3")
list_name_3_week10 <- c("C10_1", "C10_2", "C10_3", "A10_1", "A10_2", "A10_3")
list_name_4_week16 <- c("C16_1", "C16_2", "C16_3", "A16_1", "A16_2", "A16_3")
list_name_all_week <- c(list_name_0_week0, list_name_1_week1, list_name_2_week4, list_name_3_week10, list_name_4_week16)

##############################
# load data
##############################
df_raw_tpm <- data.frame()
for(idx_n in 1:length(list_name_all_week)){
    df_temp <- read.table(paste0(str_inputFilePath, list_name_all_week[idx_n], ".tsv"), header=TRUE, sep='\t')
    df_temp <- data.frame(df_temp[, c("target_id", "tpm")])
    colnames(df_temp) <- c("target_id", list_name_all_week[idx_n])
    if(nrow(df_raw_tpm) == 0){
        df_raw_tpm <- df_temp
    }else{
        df_raw_tpm <- merge(x=df_raw_tpm, y=df_temp, by=c("target_id"), all.x=T, sort=F)
    }
}

##############################
# main function
##############################
# plot corMatrix
list_func_cor <- c("pearson", "spearman")
for(idx_t in 1:length(list_func_cor)){
    cov_matrix <- cov(df_raw_tpm[,-1], method=list_func_cor[idx_t])
    cor_matrix <- cov2cor(cov_matrix)
    cor_matrix <- floor(cor_matrix*100)/100
    col_palette_sea <- colorRampPalette(brewer.pal(9, "GnBu"))
    png(paste0(str_inputFilePath, "28Mouse_", list_func_cor[idx_t], ".png"), width=1560, height=1560, pointsize=15)
    corrplot(cor_matrix, method="shade", title=paste0(list_func_cor[idx_t], ": ", nrow(df_raw_tpm)), tl.offset=.2, tl.cex=1, tl.srt=45, tl.col="#000000", type="upper", addCoef.col="white", addCoefasPercent=FALSE, col=col_palette_sea(10), is.corr=FALSE, cl.lim = c(0.5,1), mar = c(0,1,2.5,2))
    dev.off()
}
