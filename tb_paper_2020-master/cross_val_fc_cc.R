#Load libraries
library(DESeq2)
library(phyloseq)
library(dplyr)
library(yingtools2)
library(data.table)
library(ifultools)
library(stringr)
library(ggplot2)
library(ggthemes)
library(msigdbr)
library(tidyr)
library(GSVA)
library(ComplexHeatmap)
library(RColorBrewer)




# Create a directory to save figures and tables

mainDir <- "../Tables"
subDir <- "Cross_val_RF_Ctrl"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
fig_folder <- paste(mainDir,subDir,sep="/")



# Import microbiome data

phy_fc_cc <- readRDS("../data/Control/phy_mic_fc_cc.rds")


#RNASeq deseq normalized
phy_vst_rna <-  readRDS("../data/Control/phy_rna_vst_fc_cc.rds")


#phy_gsva <- rna_phy
phy_gsva <- phy_vst_rna

#I think this should be TPM--or some other normalized count 
counts <- otu_table(phy_gsva) %>% as.data.frame %>% as.matrix() #or as.matrix() #TPM data here with gene rownames and sample colnames

#metadata from RNAseq phyloseq object
metadata <- data.frame(sample_data(phy_gsva))
#%>% column_to_rownames(var = "sample") #data frame with sample metadata with rownames = colnames of counts df

#gene lists
msigdbr_show_species()
m_df = msigdbr("Homo sapiens") %>% dplyr::filter(gs_cat == "H" )

#m_df = msigdbr("Homo sapiens") %>% dplyr::filter(gs_cat == "C2" )
# KEGG
#m_df <-  m_df[grep("KEGG",m_df$gs_name),]


m_t2g = m_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
m_t2g <- m_t2g[,c("gene_symbol","gs_name")]
#m_t2g$entrez_gene <- as.character(m_t2g$entrez_gene )
hall_set <- unstack(m_t2g)

rownames(counts) <- gsub(".*_","",rownames(counts))


set.seed(1057)
gsvaRes_ssgsea <- gsva(counts, hall_set,
                       min.sz=5, max.sz=500,
                       method = "gsva",
                       kcdf="Gaussian", mx.diff=TRUE, verbose=FALSE, parallel.sz=1)

met  <- metadata
mat  <- gsvaRes_ssgsea 

met  <- metadata
met <- met[match(colnames(gsvaRes_ssgsea),met$rna_seq_ID),]

#status_col <-  brewer.pal(length(unique(met$TB_status)), "Set1")
status_col <-  c("family_contact" = "red","community_control" = "blue")
#names(status_col) <- unique(unique(met$type))
#ord_id <- order(met$type,met$Patient.ID)
library(gtools)
mixedsort(met$type)
ha_column = HeatmapAnnotation(Status =  met$type,
                              Batch = met$batchID,
                              Sex = met$sex,
                              col=list(Status = status_col))

split_cols<-  met$type
split_cols <- factor(split_cols, levels= c("community_control","family_contact","Pre HRZE","Post HRZE","Pre NTZ","Post NTZ"))

colfunc <- colorRampPalette(rev(brewer.pal(n = 7, name =
                                             "RdYlBu")))
cols <- colfunc(8)
rownames(mat) <-  gsub("HALLMARK_","",rownames(mat))


ht1 = Heatmap(mat, name = "ES", column_title = NA, 
              top_annotation = ha_column,
              #clustering_distance_rows = "euclidean",
              col = cols,
              row_names_side = "left",
              row_title_gp = gpar(fontsize = 10,fontface = "bold"), 
              row_title_rot = 0,
              column_split = split_cols,
              width=2, cluster_columns = T, 
              row_names_gp = gpar(fontsize = 9),
              row_names_max_width = max_text_width(rownames(mat), gp = gpar(fontsize = 12)),
              show_column_names = F, show_row_names = T,
              column_names_side = "bottom",na_col="white",
              show_row_dend = F)
#pdf("GSVA_raw.pdf",width = 14,height = 8)
# pdf("GSVA_deseq_norm_vst.pdf",width = 18,height = 8)
# draw(ht1)
# dev.off()


# Microbiome:
phy_sel_mic <-  prune_taxa(taxa_sums(phy_fc_cc)>0,phy_fc_cc)

phy <- phy_sel_mic
colnames(tax_table(phy))[3] <- "Kingdom"


# Prevalence 10 % 
perc_samples <-  0.1
# # 
phy_fil = filter_taxa(phy, function(x) sum(x > 0) > (perc_samples*length(x)), TRUE)
phy_fil <- prune_taxa(taxa_sums(phy_fil)> 0, phy_fil)

#phy <- yingtools2::phy.collapse(phy,c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))
############# ############# #############

#view the rownames as Phylum, Species, rather than OTU_number
t <- get.tax(phy_fil) %>% mutate(PhySpec=paste(Species,otu))
taxa_names(phy_fil) <- t$PhySpec
taxa_names(phy_fil) <- gsub("\\[","",taxa_names(phy_fil))
taxa_names(phy_fil) <- gsub("\\]","",taxa_names(phy_fil))
taxa_names(phy_fil) <- gsub(';.*','\\1', taxa_names(phy_fil))
taxa_names(phy_fil) <- gsub("ASV","(ASV",taxa_names(phy_fil))
taxa_names(phy_fil) <- gsub("$",")",taxa_names(phy_fil)) #edit this between old/new methods
taxa_names(phy_fil) <- gsub(" ","_",taxa_names(phy_fil)) #edit this between old/new methods
taxa_names(phy_fil)

# # Some manual changes in taxa table
# 
tax_table(phy_fil) <-  tax_table(phy_fil)[,3:9]


phy_sel_mic <- phy_fil


# # Deseq normalized counts data
dig <- phyloseq_to_deseq2(phy_sel_mic, ~  type) #replace this with any sample variable(s)

dig <- estimateSizeFactors(dig,"poscounts")
# Takes a while to run
dig <- estimateDispersions(dig)
dig <- DESeq(dig,fitType= "local")

plotDispEsts(dig)
# VST phyloseq for Heatmap later

phy_vst_mic <- phy_sel_mic
vst_dt <- getVarianceStabilizedData(dig)
otu_table(phy_vst_mic) <-  otu_table(vst_dt, taxa_are_rows = T)




# Lets regression each pathway against the microbes
# For RNA data
pathways_dt <-   data.frame(t(gsvaRes_ssgsea)) 

dt_SV <-  data.frame(t(otu_table(phy_vst_mic)))
dt_SV <-  dt_SV[match(rownames(pathways_dt),rownames(dt_SV)),]
# GLMNET (Elastic-Net approach)

rownames(dt_SV)
rownames(pathways_dt)


###################### Microbiome and Pathway#########################
pth_dt <- pathways_dt
lst_phylo <- list(pth_dt)
vec_phy <- c("H_Pathway")


loop <- 1
rsq_list <- list()
for(loop in 1:length(lst_phylo)){
  
  scfa_dt <-  pth_dt
  mic_dt <-  dt_SV
  match(rownames(mic_dt),rownames(scfa_dt))
  names(scfa_dt) 
  names(mic_dt)
  
  
  
  # Loop across the pathways
  comp <-  names(scfa_dt)[20]
  
  
  
  # Loop across the metabolites
  comp <-  names(scfa_dt)[4]
  
  pred_comp_list <- list()
  #loop <- 10
  for(comp in names(scfa_dt)){
    print(paste0("Computing R2 for ",comp))
    
    t_data <- cbind(mic_dt,Y_var = scfa_dt[,comp])
    
    
    library(randomForest)
    pred_list <- list()
    iter <-  1
    for(iter in 1:5){
      # Overall model using all data
      set.seed(1057 + iter)
      # Divide the data sets into training and testing
      # Random sample indexes
      train_index <- sample(1:nrow(t_data), 0.8* nrow(t_data))
      test_index <- setdiff(1:nrow(t_data), train_index)
      
      # Build X_train, y_train, X_test, y_test
      X_train <- t_data[train_index, 1:(ncol(t_data)-1)]
      y_train <- t_data[train_index, "Y_var"]
      
      X_test <- t_data[test_index, 1:(ncol(t_data)-1)]
      y_test <- t_data[test_index, "Y_var"]
      
      
      train_dt <-  cbind(X_train,y_train)
      test_dt <-  cbind(X_test,y_test)
      
      rf_model <- randomForest(y_train ~.,data=train_dt,
                               #importance = "permutation",
                               importance = T,
                               #probability = T,
                               ntree = 5000)
      
      # Predict on test set
      pred_test <-  predict(rf_model, newdata = X_test)
      pred_dt <- data.frame(Pred = as.numeric(pred_test),True = y_test)
      pred_dt$iter <- iter 
      
      pred_list[[iter]] <-  pred_dt
    }  
    
    pred_err_dt <-  do.call("rbind",pred_list)
    pred_err_dt$Comp <-  comp
    pred_err_dt$R2 <-  caret::R2(pred_err_dt$Pred,pred_err_dt$True)
    
    pred_comp_list[[comp]] <- pred_err_dt
    
    
  }    
  
  pred_comp_dt <- do.call("rbind",pred_comp_list)
  
  
  rsq_dt <- unique(pred_comp_dt[,c("Comp","R2")])
  rsq_dt$Grp <-  vec_phy[loop]
  
  rsq_list[[loop]] <- rsq_dt
  #Import library for Metrics
  # library(Metrics)
  # print(paste0('MAE: ' , mae(pred_comp_dt$True,pred_comp_dt$Pred) ))
  # print(paste0('R2: ' ,caret::postResample(pred_comp_dt$Pred , pred_comp_dt$True)['Rsquared'] ))
  # 
  # pred_p <-     ggplot(pred_comp_dt,aes(x = Pred, y = True)) +
  #   geom_abline()+
  #   geom_point(size = 3, alpha = 0.5)+
  #   facet_wrap(~Comp,scales = "free")+
  #   #xlim(0,3)+
  #   #ylim(0,3)+
  #   theme_bw()
  # pdf(paste0(fig_folder,"Predictions_Mic_Pathways_",vec_phy[loop],".pdf"),height = 20, width = 18)
  # print(pred_p)
  # dev.off() 
  
  
}    

rsq_list_dt <- do.call("rbind",rsq_list)

#write.csv(rsq_list_dt,paste0(fig_folder,"/Predictions_R2_Mic_H_Pathways.csv"))
#write.csv(rsq_list_dt,paste0(fig_folder,"/Predictions_R2_Mic_Kegg_Pathways.csv"))
write.csv(rsq_list_dt,paste0("../data/Control/Predictions_R2_Mic_H_Pathways.csv"))


