
library(openxlsx)


# Create a directory to save figures and tables

mainDir <- "../Supp_Tab"
dir.create(mainDir, showWarnings = TRUE,recursive = TRUE)
tab_folder <- paste(mainDir,sep="/")



############Supplementary Fig 2#################
# Import all the files for Fig 2

# LME Analysis tables

lme_files = list.files(path = "../Tables/Tabs_2",
                       pattern="*.csv",
                       full.names = T)
lme_list = lapply(lme_files, read.csv,header = T)

names(lme_list) <- c("LME_InvSimpson","LME_Shannon","LME_TTP")

lme_list  <-  lapply(names(lme_list), function(i){
                     x <- lme_list[[ i ]]
                     # set 1st column to a new name
                      names(x)[1] <- "Variable"
                     # return
                     x
               })

names(lme_list) <- c("LME_InvSimpson","LME_Shannon","LME_TTP")

# Merge into one big dataframe using a "cool" Reduce function 
write.xlsx(lme_list, file = paste0(tab_folder,"/Supp_Tab_2.xlsx"))



############Supplementary Tables Differential Analysis #################################

# Import Differential analysis from Clinical Trial

diff_clinical_mic_files <-  list.files(path = "../Tables/Limma_mic_drug",
                                       pattern="*.csv",
                                       full.names = T)
fnames <- list.files(path = "../Tables/Limma_mic_drug",
                     pattern="*.csv",
                     full.names = F)
fnames <-  gsub(".csv","",fnames)

diff_mic_clinical_list = lapply(diff_clinical_mic_files, read.csv,header = T)

names(diff_mic_clinical_list) <- fnames

diff_mic_clinical_list  <-  lapply(names(diff_mic_clinical_list), function(i){
  x <- diff_mic_clinical_list[[ i ]]
  # set 1st column to a new name
  names(x)[1] <- "ASVs"
  # return
  x
})

names(diff_mic_clinical_list) <- fnames


write.xlsx(diff_mic_clinical_list, file = paste0(tab_folder,"/Supp_Tab_Clinical_Mic.xlsx"))


# RNASeq

diff_clinical_rna_files <-  list.files(path = "../Tables/Limma_rna_drug",
                                       pattern="*.csv",
                                       full.names = T)

fnames <- list.files(path = "../Tables/Limma_rna_drug",
                     pattern="*.csv",
                     full.names = F)
fnames <-  gsub(".csv","",fnames)

diff_clinical_rna_list = lapply(diff_clinical_rna_files, read.csv,header = T)

names(diff_clinical_rna_list) <- fnames

write.xlsx(diff_clinical_rna_list, file = paste0(tab_folder,"/Supp_Tab_Clinical_RNA.xlsx"))



# Import Differential analysis from EBA 
# Microbiome
diff_eba_mic_files <-  list.files(path = "../Tables/EBA_analysis",
                                  pattern="*.csv",
                                  full.names = T)
fnames <- list.files(path = "../Tables/EBA_analysis",
                     pattern="*.csv",
                     full.names = F)
fnames <-  gsub(".csv","",fnames)

diff_eba_mic_list = lapply(diff_eba_mic_files, read.csv,header = T)

names(diff_eba_mic_list) <- fnames

write.xlsx(diff_eba_mic_list, file = paste0(tab_folder,"/Supp_Tab_EBA_Mic.xlsx"))


# RNA
diff_eba_rna_files <-  list.files(path = "../Tables/EBA_analysis/Limma_rnaseq_eba",
                                  pattern="*.csv",
                                  full.names = T)
fnames <- list.files(path = "../Tables/EBA_analysis/Limma_rnaseq_eba",
                     pattern="*.csv",
                     full.names = F)
fnames <-  gsub(".csv","",fnames)

# The characters are limited to 30 characters

fnames <- gsub("EBA_","",fnames)
fnames <- gsub("differential","diff",fnames)
fnames <- gsub("Extraction","Extr",fnames)


nchar(fnames)

diff_eba_rna_list = lapply(diff_eba_rna_files, read.csv,header = T)

names(diff_eba_rna_list) <- fnames

write.xlsx(diff_eba_rna_list, file = paste0(tab_folder,"/Supp_Tab_EBA_RNA.xlsx"))


############### RF model FC and CC ####################

r2_file <-  read.csv("../data/Control/Predictions_R2_Mic_H_Pathways.csv")
r2_file$X <- NULL
r2_file$Grp <- NULL

pred_rf_ctrl_list <-  readRDS("../data/Control/RF_imp_FC_CC.rds")
pred_rf_dt <- do.call("rbind",pred_rf_ctrl_list)
rownames(pred_rf_dt) <-  NULL

# Diff RF Model
rf_diff_res <-  readRDS("../Figs/RF_diff_model/RF_results_pathway.rds")
rf_diff_dt <-  do.call("rbind",rf_diff_res)
rownames(rf_diff_dt) <- NULL

list_rf_ctrl <- list(r2_file,pred_rf_dt,rf_diff_dt)
names(list_rf_ctrl) <- c("R2_rf_pathways","RF_Ctrl_Var_Imp_ALE_Slope","RF_Diff_Var_Imp_ALE_Slope")

write.xlsx(list_rf_ctrl, file = paste0(tab_folder,"/Supp_Tab_RF.xlsx"))




############## GSVA ####################

gsva_all <- read.csv("../data/GSVA/GSVA_all_Hallmark.csv")

# LME pathways Clinical Trial and EBA

lme_clinical <- read.csv("../Tables/Pathways_LME_EBA_Clinical/LME_Clinical_Post_Pre.csv")
lme_clinical$X <- NULL

lme_eba  <- read.csv("../Tables/Pathways_LME_EBA_Clinical/LME_EBA_Post_Pre.csv")
lme_eba$X <- NULL

lme_list <- list(lme_clinical,lme_eba)
names(lme_list) <- c("LME_Hall_Pathway_HRZE_NTZ","LME_Hall_Pathway_EBA")

write.xlsx(lme_list, file = paste0(tab_folder,"/Supp_Tab_LME_Pathways_Diff.xlsx"))

