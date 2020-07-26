#Load libraries
library(phyloseq)
library(dplyr)
library(ggplot2)
library(gtools)
library(reshape2)
library(ggplot2)
library(edgeR)
library(limma)


# Create a directory to save figures and tables
mainDir <- "../Figs"
subDir <- "RF_diff_model"

dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
fig_folder <- paste(mainDir,subDir,sep="/")


# Import microbiome data 
# Clinical Trial
# Microbiome
phy_clinical_mic <- readRDS("../data/mic/phy_mic_fil_paired.rds")
# RNASeq
phy_clinical_rna <-  readRDS("../data/genes/phy_rna_ntz_hrze.rds")

# Only keep samples that pairs with RNASeq data
phy_clinical_mic <-  prune_samples(sample_names(phy_clinical_rna), phy_clinical_mic)



# Only keep the needed info:
# Time, TTP, SampleID, PatientID

sm_dt <-  data.frame(sample_data(phy_clinical_rna))
sm_dt <- sm_dt[,c("Patient.ID","TB_status","Av_TTP")]
sm_dt$SampleID <- rownames(sm_dt)

# Remove samples that don't have TTP value
sm_dt <- sm_dt[!is.na(sm_dt$Av_TTP),]

# Keep only paired samples 
tab_freq <- data.frame(table(sm_dt$Patient.ID))
tab_freq <- tab_freq[tab_freq$Freq > 1, ]
sm_dt <-  sm_dt[sm_dt$Patient.ID  %in% tab_freq$Var1,]
sm_dt$Time <- ifelse(sm_dt$TB_status == "pretreatment",0,14)
sm_dt$TB_status <- NULL

# Assign sample data back to phyloseq obj
sample_data(phy_clinical_mic) <- sample_data(sm_dt) 
sample_data(phy_clinical_rna) <- sample_data(sm_dt)

phy_clinical_mic
phy_clinical_rna

# EBA
phy_eba_mic <- readRDS("../data/EBA/phy_mic_eba.rds")
phy_eba_rna <-  readRDS("../data/EBA/phy_rna_tb.rds")
phy_eba_mic
phy_eba_rna

sm_dt_eba_mic <-  data.frame(sample_data(phy_eba_mic))
sm_dt_eba_mic <-  sm_dt_eba_mic[,c("Patient.ID","Average.TTP","sam_id","Time")]

# Keep only Time that is 0, 14 days and Two months to match with RNAseq
sm_dt_eba_mic <- sm_dt_eba_mic[sm_dt_eba_mic$Time %in% c("Day0","Day14","TwoMonths"),]
sm_dt_eba_mic$Time <-  gsub("Day", "",sm_dt_eba_mic$Time)
sm_dt_eba_mic$Time <-  gsub("TwoMonths", "56",sm_dt_eba_mic$Time)
sm_dt_eba_mic$Average.TTP <- as.numeric(as.character(sm_dt_eba_mic$Average.TTP))
sm_dt_eba_mic <- sm_dt_eba_mic[!is.na(sm_dt_eba_mic$Average.TTP),]

# Keep only paired samples 
tab_freq <- data.frame(table(sm_dt_eba_mic$Patient.ID))
tab_freq <- tab_freq[tab_freq$Freq > 1, ]
sm_dt_eba_mic <-  sm_dt_eba_mic[sm_dt_eba_mic$Patient.ID  %in% tab_freq$Var1,]

# Match it with clinical data
names(sm_dt_eba_mic) <- c("Patient.ID", "Av_TTP", "SampleID", "Time")

sm_dt_eba_rna <-  sm_dt_eba_mic
sm_dt_eba_rna$SampleID <-  paste0(sm_dt_eba_mic$Patient.ID,"_D",sm_dt_eba_mic$Time)
rownames(sm_dt_eba_rna) <- sm_dt_eba_rna$SampleID

# Assign sample data back to phyloseq obj
sample_data(phy_eba_mic) <- sample_data(sm_dt_eba_mic) 
sample_data(phy_eba_rna) <- sample_data(sm_dt_eba_rna)

sample_names(phy_eba_mic) <- paste(sample_data(phy_eba_mic)$Patient.ID,"_D",
                                   sample_data(phy_eba_mic)$Time,sep="")
sample_names(phy_eba_rna)


# Merge microbiome dataset
mer_phy_mic <- merge_phyloseq(phy_eba_mic,phy_clinical_mic)
sample_names(mer_phy_mic) <-  make.names(paste(sample_data(mer_phy_mic)$Patient.ID,"_D",
                                               sample_data(mer_phy_mic)$Time,sep=""))

# Merge rnaseq data
mer_phy_rna <- merge_phyloseq(phy_eba_rna,phy_clinical_rna)
mer_phy_rna
sample_names(mer_phy_rna) <-  make.names(paste(sample_data(mer_phy_rna)$Patient.ID,"_D",
                                               sample_data(mer_phy_rna)$Time,sep=""))


# GSVA matrix

phy_gsva <-  readRDS("../data/GSVA/phy_GSVA_all_hallmark.rds")
#phy_gsva <-  readRDS("../../Figs_vb/GSVA/phy_GSVA_lme_kegg.rds")
phy_gsva <-  subset_samples(phy_gsva, TB_status %in% c("pretreatment","HRZE","NTZ",
                                                       "EBA_0","EBA_14","EBA_56"))
sm_dt_gsva <- data.frame(sample_data(phy_gsva))
sm_dt_gsva$TB_status <-  gsub(".*_","",sm_dt_gsva$TB_status) 
sm_dt_gsva$TB_status[sm_dt_gsva$TB_status == "pretreatment"] <-  "0"
sm_dt_gsva$TB_status[sm_dt_gsva$TB_status == "HRZE"] <-  "14"
sm_dt_gsva$TB_status[sm_dt_gsva$TB_status == "NTZ"] <-  "14"
names(sm_dt_gsva)[3]  <- "Time"

sm_dt_gsva <- unique(sm_dt_gsva)



sample_data(phy_gsva) <-  sample_data(sm_dt_gsva)
sample_names(phy_gsva) <-  make.names(paste(sample_data(phy_gsva)$Patient.ID,"_D",
                                            sample_data(phy_gsva)$Time,sep=""))

phy_sel_gsva <- prune_samples(sample_names(mer_phy_mic),phy_gsva)
sample_data(phy_sel_gsva) <- sample_data(mer_phy_mic)

mer_phy_mic
mer_phy_rna
phy_sel_gsva


# Lets look at the early resolution Day 14 only

# Select patients that have Day 0 and Day 14
phy_sel_mic <-  subset_samples(mer_phy_mic,Time %in% c("0","14"))
# Remove unpaired samples
sm_dt <-  data.frame(sample_data(phy_sel_mic))
freq_tab <-  data.frame(table(sm_dt$Patient.ID))

phy_sel_mic <-  subset_samples(phy_sel_mic,Patient.ID %in% freq_tab$Var1[freq_tab$Freq > 1])

phy_sel_rna <-  prune_samples(sample_names(phy_sel_mic),mer_phy_rna)

phy_sel_gsva <-  prune_samples(sample_names(phy_sel_mic), phy_sel_gsva)


# Filter species that are prevalent 
perc_samples <-  0.15
phy_fil = filter_taxa(phy_sel_mic, function(x) sum(x > 0) > (perc_samples*length(x)), TRUE)
phy_fil <- prune_taxa(taxa_sums(phy_fil)> 0, phy_fil)
phy_fil

# Then I would use the edgeR normalized counts

cnt_mic <- data.frame(otu_table(phy_fil))
dge <- DGEList(cnt_mic)
cal <- calcNormFactors(dge,method = "TMM")
mat_mic  <- cpm(cal,log = TRUE)
OTU = otu_table(as.matrix(mat_mic), taxa_are_rows = TRUE)
colnames(OTU)
otu_table(phy_fil) <- OTU

phy_mic_sel <- phy_fil



# Compute the change between Day 14 and Day 0 

diff_calc <- function (phy){
  # Divide the data in two matrices Day 0 and Day 14
  
  
  phy_t1 <- subset_samples(phy,Time %in% c("0"))
  sample_names(phy_t1) <-  sample_data(phy_t1)$Patient.ID
  mat_t1 <- data.frame(otu_table(phy_t1))
  
  
  phy_t2 <- subset_samples(phy,Time %in% c("14"))
  sample_names(phy_t2) <-  sample_data(phy_t2)$Patient.ID
  mat_t2 <- data.frame(otu_table(phy_t2))
  
  # Now match the columns of the matrices
  mat_t2 <- mat_t2[,match(colnames(mat_t1),colnames(mat_t2))]
  
  diff_mat <-  mat_t2 - mat_t1
  return(diff_mat)
  
}

diff_mat_mic <- diff_calc(phy_mic_sel)
diff_mat_path <- diff_calc(phy_sel_gsva)

# Diff in TTP 
meta_dt <-  data.frame(sample_data(phy_mic_sel))

meta_dt$Patient.ID <-  as.character(meta_dt$Patient.ID)
meta_dt$Time <-  as.numeric(meta_dt$Time)
meta_dt$Av_TTP <- log(meta_dt$Av_TTP)
diff_ttp <- meta_dt %>%
  dplyr::group_by(Patient.ID) %>%
  arrange(Time) %>%
  dplyr::mutate(Diff = Av_TTP - lag(Av_TTP))

diff_ttp <- diff_ttp[!is.na(diff_ttp$Diff),]
# Keep the needed variables
diff_ttp <- diff_ttp[,c("Patient.ID","Diff")]


# Match the diff_mat 
diff_ttp <-  diff_ttp[match(colnames(diff_mat_mic),make.names(diff_ttp$Patient.ID)),]
diff_mat_path <- diff_mat_path[,match(colnames(diff_mat_mic),colnames(diff_mat_path)),]


colnames(diff_mat_path)
colnames(diff_mat_mic)
diff_ttp$Patient.ID



# Set up RF 

# Predict pathways using microbiome and TTP

mic_dt <-  data.frame(t(diff_mat_mic))

path_dt <-  data.frame(t(diff_mat_path))

match(rownames(mic_dt),rownames(path_dt))


# Loop across the pathways
p_ann <-  read.csv("../Hallmark/Hallmark_pathway_annotation.csv",
                   stringsAsFactors = F)
p_ann$Hallmark.Name <-  paste0("HALLMARK_",p_ann$Hallmark.Name)
immune_path <- p_ann$Hallmark.Name[p_ann$Process.Category %in% "immune"]

lme_pathways <- readRDS("../data/GSVA/lme_eba_clinical_two_week.rds")



#lme_pathways <- immune_path
comp <-  lme_pathways[1]

imp_dt_list <- list()
for(comp in lme_pathways){
  
  print(paste0("Runnning RF for ",comp))  
  t_data <- cbind(mic_dt,TTP = diff_ttp$Diff,Y_var = path_dt[,colnames(path_dt) %in% comp])
  
  
  library(randomForest)
  # Overall model using all data
  set.seed(1057)
  rf_model <- randomForest(Y_var~.,data=t_data,
                           #importance = "impurity_corrected",
                           #importance = "permutation",
                           importance = T,
                           #probability = T,
                           ntree = 5000)
  
  
  # varImpPlot(rf_model)
  
  library(vita)
  #t_data <- dt_orig
  train_data <-  t_data[,!colnames(t_data) %in% "Y_var"]
  pimp <- PIMP(train_data, t_data$Y_var,rf_model,parallel=TRUE, ncores=10 ,seed = 340)
  p_test <- PimpTest(pimp)
  #pimp_all <- data.frame( summary(p_test,pless = 0.05))
  pimp_all <- data.frame(orig =  p_test$VarImp, p_test$pvalue)
  pimp_all$padj <- pimp_all$p.value
  imp_dt <- pimp_all[pimp_all$VarImp > 0 & pimp_all$padj <= 0.05,]     
  
  
  
  imp_dt$variable <- rownames(imp_dt)
  names(imp_dt)
  imp_dt <- imp_dt[order(imp_dt$VarImp,decreasing = T),]
  
  if(nrow(imp_dt)>0){
    
    comp_dt <- imp_dt
    comp_dt$Comp <- comp
    #imp_dt_list[[comp]] <- comp_dt 
    
    
    
    # Top 30
    sel_imp_dt <- imp_dt
    
    # library(ggplot2)
    # p_var <- ggplot(sel_imp_dt, aes(x=reorder(variable,VarImp), y=VarImp))+ 
    #   geom_bar(stat="identity", position="dodge")+ coord_flip()+
    #   ylab("Variable Importance")+
    #   xlab("")+theme_bw()
    # 
    # pdf(paste0(fig_folder,"/Variable_imp_RF_",comp,".pdf"), width = 8, height = 10)
    # print(p_var)
    # dev.off()
    
    
    # ALE Plots 
    # ALE analysis
    library(ALEPlot)
    library(parallel)
    yhat <- function(X.model, newdata) as.numeric(randomForest:::predict.randomForest(X.model, newdata))
    
    ALE_summary_list <- list()
    imp_bugs <-  as.character(sel_imp_dt$variable)
    
    train_data <-  t_data[,!colnames(t_data) %in% "Y_var"]
    
    ale_func <-  function(rep){
      
      if(imp_bugs[len] %in% colnames(train_data)){
        idx <- which(colnames(train_data) == imp_bugs[len]) 
        ale_dt  <- data.frame(ALEPlot(X = train_data,X.model =  rf_model, 
                                      pred.fun = yhat, J=idx, K = 40, NA.plot =F))
        
        ale_dt$pred <-  imp_bugs[len]
        ale_dt$nrep <- rep
        res <-  ale_dt
        # rm(ale_dt)
      }
    }
    
    pp_bug_list <- list()
    for(len in seq_along(imp_bugs)){
      print(len)
      ale_dt <- ale_func(1)
      pp_bug_list[[len]] <- ale_dt
    }
    
    pp_dt_final <- do.call("rbind",pp_bug_list)
    # Save the partial dependence data
    #write.csv(pp_dt_final,paste(data_folder,paste0('RF_ALE_',analyte,'.csv'),sep="/"))
    
    pp_dt_final$f.values <-  as.numeric(pp_dt_final$f.values)
    pp_dt_final$x.values <-  as.numeric(pp_dt_final$x.values)
    pp_dt_final <-  na.omit(pp_dt_final)
    
    
    pp_dt_final$pred <- factor(pp_dt_final$pred,levels = imp_bugs)
    
    library(ggthemes)
    p_p <- ggplot(data = pp_dt_final ,aes(x = x.values, y = f.values)) +
      geom_point(size = 3)+
      #geom_line(aes(group = nrep), alpha = 0.6)+
      #geom_line(data = pp_dt_final_median ,aes(x = x, y = y_med, color = trial))+
      #geom_smooth()+
      #geom_line(data = pp_dt_final_median ,aes(x = x, y = y_med, color = "red"))+
      theme_base()+
      facet_wrap(~pred,scales = "free")+
      xlab("Microbes")+
      theme(legend.position="none")
    #  print(p_p)
    pdf(paste0(fig_folder,'/RF_ALE_',comp,'.pdf'),height = 15, width = 15)
    print(p_p)
    dev.off()
    
    library(tidyr)
    library(dplyr)
    slope_dt <- pp_dt_final %>%
      group_by(pred) %>%
      do({
        mod <-  data.frame(summary(lm( f.values ~ x.values , data = .))$coefficients)[,c(1,4)]
        names(mod) <- c("val","p_val")
        mod$coeff <- rownames(mod)
        
        data.frame(Intercept = mod$val[1],
                   Slope = mod$val[2],
                   p_val_int =  mod$p_val[1],
                   p_val_slope =  mod$p_val[2] )
      }) %>% as.data.frame()
    
    head(slope_dt)
    
    
    slope_dt <- slope_dt[,c("pred","Slope")]
    melt_slope_dt <- reshape::melt(slope_dt)
    
    head(melt_slope_dt)
    melt_slope_dt$value <- ifelse(melt_slope_dt$value > 0,1,0)
    # slope_dt$Slope_logic <- ifelse(slope_dt$Slope > 0,"Positive","Negative")
    # slope_dt$Slope_logic2 <- ifelse(slope_dt$Slope > 0,"Positive","Negative")
    
    head(melt_slope_dt)
    
    melt_slope_dt$pred <-  factor(melt_slope_dt$pred, levels = rev(imp_bugs))
    
    library(ggplot2)
    gg_slope <- ggplot(melt_slope_dt, aes(x = variable, y = pred )) +
      geom_tile(color = "white", size = 0.1, aes(fill = value))+
      scale_fill_gradient2( high = "#e34a33", mid = "#2c7fb8")+
      #scale_fill_gradient( low = "blue", high = "red")+
      #scale_fill_manual(name = "",values = c("0" = 'blue', "1" = 'red'))+
      theme_classic()+
      theme(#axis.title.x = element_blank(),
        legend.position = "none",
        axis.ticks.x = element_blank(),
        #axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_blank(),
        axis.text.y.left= element_blank())+
      xlab("Relation")
    
    #gg_slope
    
    imp_slope <-  sel_imp_dt
    imp_slope$Slope <-  slope_dt$Slope
    imp_slope$Relation <- imp_slope$Slope
    imp_slope$Relation <-  ifelse(imp_slope$Relation > 0, "Positive","Negative")
    imp_slope$Comp <- comp
    
    imp_dt_list[[comp]] <- imp_slope 
    #write.csv(imp_slope,paste0(fig_folder,"/ESBL_decol_LBP_Importance.csv"))
    
    
    
    # Variable importance and relation
    
    var_slope_dt <-  merge(sel_imp_dt,slope_dt,by.x = "variable",by.y= "pred")
    
    var_slope_dt$Relation <-  ifelse(var_slope_dt$Slope > 0, "Positive","Negative")
    
    
    library(ggplot2)
    p_var <- ggplot(var_slope_dt, aes(x=reorder(variable,VarImp), y=VarImp,fill = Relation))+ 
      geom_bar(stat="identity", position="dodge")+
      scale_fill_manual(values = c("Negative"="#e34a33", "Positive" = "#2c7fb8"))+
      coord_flip()+
      ylab("Variable Importance")+
      xlab("Species")+theme_bw()+
      ggtitle(comp)+
      theme(axis.text.y=element_text(face = "bold"))
    
    pdf(paste0(fig_folder,"/Variable_imp_RF_",comp,".pdf"), width = 6, height = 10)
    print(p_var)
    dev.off()
    
    
    
    #Predictors to visualize
    pred_vec <-  as.character(sel_imp_dt$variable)
    dt_sel <- t_data[,c(pred_vec,"Y_var")]
    dt_sel$sample_id <- rownames(dt_sel)
    head(dt_sel)
    #rownames(X_sel)
    
    library(reshape2)
    melt_dt <-  reshape2::melt(dt_sel, id.vars  = c("Y_var","sample_id"))
    
    
    head(melt_dt)
    
    library(tidyverse)
    
    xy_p <- ggplot(melt_dt, aes(y = Y_var, x = value)) +
      geom_point()+
      facet_wrap(~ variable, scales = "free") +
      theme_bw() +
      ylab("Abundance") +
      xlab("Induction")
    pdf(paste0(fig_folder,"/XY_plot_",comp,".pdf"), width = 10, height = 10)
    print(xy_p)
    dev.off()
  }
  
}

# Based on the importance , lets draw the heatmap which only demonstrates the
# importance of the species towards the Compounds

saveRDS(imp_dt_list,"../Figs/RF_diff_model/RF_results_pathway.rds")

comp_imp <- do.call("rbind",imp_dt_list)

# Top predictor 
var_sel <- comp_imp %>%
  group_by(Comp) %>%
  filter(VarImp == max(VarImp))


# dplyr::group_by(Comp)%>%
# dplyr::arrange(VarImp) %>%
# top_n(1,variable)

# Select only top 10
sel_imp <- comp_imp %>%
  dplyr::group_by(Comp) %>%
  dplyr::arrange(desc(VarImp)) %>%
  dplyr::slice_max(n= 10,VarImp)%>% data.frame()
# 
#sel_imp <- comp_imp
#write.csv(comp_imp,paste0(res_folder,"/RF_res_Mic_Genes.csv"))

sel_imp$value <- 1
sel_imp$value[sel_imp$Slope < 0 ] <- -1 
sel_imp <- sel_imp[,c("variable","Comp","value")]

library(tidyr)

comp_imp_w <- sel_imp %>%
  pivot_wider(names_from = Comp, values_from = value,values_fill = list(value = 0))%>%
  data.frame()

comp_vimp <- do.call("rbind",imp_dt_list)

comp_vimp <- comp_vimp[,c("variable","Comp","VarImp")]
comp_vimp_w <- comp_vimp %>%
  pivot_wider(names_from = Comp, values_from = VarImp,values_fill = list(value = 0))%>%
  data.frame()

comp_vimp_w[is.na(comp_vimp_w)] <-   0

rownames(comp_imp_w) <-  comp_imp_w$variable
comp_imp_w$variable <- NULL
rownames(comp_vimp_w) <-  comp_vimp_w$variable
comp_vimp_w$variable <- NULL
comp_vimp_w <-  scale(comp_vimp_w,center = F)


library(RColorBrewer)
library(circlize)
# Lets only draw the immune pathways from lme analysis
mat <- as.matrix(comp_imp_w) 
#col_binary <-  colorRampPalette(c('#7fc97f','white','#d95f02'))(3)

col_binary <- colorRampPalette(rev(brewer.pal(n = 3, name = "RdYlBu")))

col_binary <-  c("#91BFDB","white","#FC8D59")



#sel_cols = colorRampPalette(c("white", "darkblue"))(2)
col_bar = colorRamp2(c(0,1), c("white", "darkblue"))


# Display All 

#lme_pathway 
sel_imm_path <- p_ann[p_ann$Hallmark.Name %in% lme_pathways,]

mat_imm <- mat[,colnames(mat) %in% sel_imm_path$Hallmark.Name]


sel_p_ann <- p_ann[p_ann$Hallmark.Name %in% colnames(mat_imm),]
sel_p_ann <- sel_p_ann[match(colnames(mat_imm),sel_p_ann$Hallmark.Name),]
# # Match it with the colnames of the matrix
# rsq_dt <-  rsq_dt[match(colnames(mat),rsq_dt$Comp),]
# 
sp_cols <- factor(sel_p_ann$Process.Category)

mat_imm <- mat_imm[apply(mat_imm, 1, function(x) !all(x==0)),]


# Find the index of  top predictor for each pathway
head(var_sel)

imm_var_sel <- var_sel[var_sel$Comp %in% colnames(mat_imm),]


var <- imm_var_sel$Comp[1]

idx_c <- c()
idx_r <- c()
for (var in imm_var_sel$Comp){
  
  spec <- imm_var_sel$variable[imm_var_sel$Comp %in% var]
  idx_c <- c(idx_c,which(colnames(mat_imm) == var))
  idx_r <- c(idx_r,which(rownames(mat_imm) == spec))
}



# Add tax dt info
phy_ht <- phy_mic_sel
taxa_names(phy_ht) <- make.names(taxa_names(phy_ht))

phy_ht <-  prune_taxa(rownames(mat_imm),phy_ht)

library(yingtools2)
library(data.table)
phy_melt <-  get.otu.melt(phy_ht)
tax_dt <-  unique(phy_melt[,c("otu","Order","Phylum")])
topN <- length(unique(tax_dt$Order))
match_tax_dt <- data.frame(tax_table(phy_ht))
match_tax_dt <- match_tax_dt[order(match_tax_dt$Phylum,match_tax_dt$Order),]
rownames(match_tax_dt)


phylum_col <- c("Firmicutes"= "#9C854E",
                "Bacteroidetes" = "#51AB9B",
                "Actinobacteria" = "#A77097",
                "Proteobacteria" = "red",
                "Verrucomicrobia" = "#33445e",
                "Tenericutes" = "#e2eab7",
                "Fusobacteria"= "#653131",
                "Cyanobacteria" = "#808080",
                "Euryarchaeota" = "#8c8c8c",
                "Spirochaetes" = "grey50",
                "Candidatus Melainabacteria" = "#999999")

dt_col_phylum <- data.frame(phylum_col)
dt_col_phylum$Phylum <-  rownames(dt_col_phylum)
shades_num <- match_tax_dt %>%
  group_by(Phylum) %>%
  mutate(unique_types = n_distinct(Order))%>% select(Phylum,unique_types) %>% unique %>% as.data.frame()

shades_num <- merge(shades_num, dt_col_phylum, by = "Phylum")
shades_num$phylum_col <-  as.character(shades_num$phylum_col )
tax_sel <-  unique(match_tax_dt[,c("Phylum","Order")])

#shades_num$phy_col <- as.vector(phylum_col[order(names(phylum_col))])
order_col <- mapply(function(x, y){ shades(color = x, ncolor = y, variation = 1)}, x= shades_num$phylum_col ,y = shades_num$unique_types,
                    SIMPLIFY = F)

order_col <- as.vector(unlist(order_col))
names(order_col) <- tax_sel$Order

mycol <- order_col

# Add a row for TTP
match_tax_dt[nrow(match_tax_dt)+1,] <- NA
rownames(match_tax_dt)[nrow(match_tax_dt)] <- "TTP"
match_tax_dt <- match_tax_dt[match(rownames(mat_imm),rownames(match_tax_dt)),]

library(ComplexHeatmap)

ha1 = HeatmapAnnotation(Order = match_tax_dt$Order,
                        col = list(Order = order_col), which = "row")
split_rows <-  match_tax_dt$Phylum
split_rows <- factor(split_rows, levels= unique(split_rows))


# Replace dots with space
rownames(mat_imm)  <- gsub("\\."," ",rownames(mat_imm))
rownames(mat_imm)  <- gsub("X","",rownames(mat_imm))

colnames(mat_imm) <- gsub("HALLMARK_","",colnames(mat_imm))
colnames(mat_imm) <- gsub("_"," ",colnames(mat_imm))

# Find the index where t


ht1  = Heatmap(mat_imm,name = "Relation",
               column_split = sp_cols,
               row_split =  split_rows, row_gap = unit(2, "mm"),
               row_title_gp = gpar(fontsize = 10,fontface = "bold"),
               row_title_rot = 0,
               column_title_gp = gpar(fontsize = 10,fontface = "bold"),
               column_title_rot = 90,
               left_annotation = ha1,
               show_parent_dend_line = F,
               border = T,
               rect_gp = gpar(col = "gainsboro"),
               row_names_side = "left",
               show_row_dend = F,
               col = col_binary,
               row_names_max_width = max_text_width(rownames(mat_imm),
                                                    gp = gpar(fontsize = 12)),
               column_names_max_height = max_text_width(colnames(mat_imm),
                                                        gp = gpar(fontsize = 12)),
               heatmap_legend_param = list(
                 title = "Relation", at = c(-1, 1), 
                 labels = c("Neg", "Pos")
               ),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 #if(mat_imm[i, j] > 0)
                 for (k in 1:length(idx_r)){
                   if(i == idx_r[k] & j == idx_c[k]){
                     print(paste0(i," ",j))
                     grid.points(x, y, pch = 16, size = unit(4, "mm"))
                   }
                 }
               }
               # layer_fun = function(j, i, x, y, w, h, fill) {
               #   ind_mat = restore_matrix(j, i, x, y)
               #   
               #   
               #   ind = ind_mat[1, ]
               #   grid.points(x[ind], y[ind], pch = 16, size = unit(4, "mm"))
               # }
               
)


pdf(paste0(fig_folder,"/Heatmap_imp_Mic_lme_Pathways.pdf"),height = 20, width = 15,useDingbats = F)
draw(ht1)
dev.off()



# Display only immune pathways
#lme_pathway 
sel_imm_path <- p_ann[p_ann$Hallmark.Name %in% immune_path,]

mat_imm <- mat[,colnames(mat) %in% sel_imm_path$Hallmark.Name]


sel_p_ann <- p_ann[p_ann$Hallmark.Name %in% colnames(mat_imm),]
sel_p_ann <- sel_p_ann[match(colnames(mat_imm),sel_p_ann$Hallmark.Name),]
# # Match it with the colnames of the matrix
# rsq_dt <-  rsq_dt[match(colnames(mat),rsq_dt$Comp),]
# 
sp_cols <- factor(sel_p_ann$Process.Category)

mat_imm <- mat_imm[apply(mat_imm, 1, function(x) !all(x==0)),]


# Find the index of  top predictor for each pathway
head(var_sel)

imm_var_sel <- var_sel[var_sel$Comp %in% colnames(mat_imm),]


var <- imm_var_sel$Comp[1]

idx_c <- c()
idx_r <- c()
for (var in imm_var_sel$Comp){
  
  spec <- imm_var_sel$variable[imm_var_sel$Comp %in% var]
  idx_c <- c(idx_c,which(colnames(mat_imm) == var))
  idx_r <- c(idx_r,which(rownames(mat_imm) == spec))
}



# Add tax dt info
phy_ht <- phy_mic_sel
taxa_names(phy_ht) <- make.names(taxa_names(phy_ht))

phy_ht <-  prune_taxa(rownames(mat_imm),phy_ht)

library(yingtools2)
library(data.table)
phy_melt <-  get.otu.melt(phy_ht)
tax_dt <-  unique(phy_melt[,c("otu","Order","Phylum")])
topN <- length(unique(tax_dt$Order))
match_tax_dt <- data.frame(tax_table(phy_ht))
match_tax_dt <- match_tax_dt[order(match_tax_dt$Phylum,match_tax_dt$Order),]
rownames(match_tax_dt)


phylum_col <- c("Firmicutes"= "#9C854E",
                "Bacteroidetes" = "#51AB9B",
                "Actinobacteria" = "#A77097",
                "Proteobacteria" = "red",
                "Verrucomicrobia" = "#33445e",
                "Tenericutes" = "#e2eab7",
                "Fusobacteria"= "#653131",
                "Cyanobacteria" = "#808080",
                "Euryarchaeota" = "#8c8c8c",
                "Spirochaetes" = "grey50",
                "Candidatus Melainabacteria" = "#999999")

dt_col_phylum <- data.frame(phylum_col)
dt_col_phylum$Phylum <-  rownames(dt_col_phylum)
shades_num <- match_tax_dt %>%
  group_by(Phylum) %>%
  mutate(unique_types = n_distinct(Order))%>% select(Phylum,unique_types) %>% unique %>% as.data.frame()

shades_num <- merge(shades_num, dt_col_phylum, by = "Phylum")
shades_num$phylum_col <-  as.character(shades_num$phylum_col )
tax_sel <-  unique(match_tax_dt[,c("Phylum","Order")])

#shades_num$phy_col <- as.vector(phylum_col[order(names(phylum_col))])
order_col <- mapply(function(x, y){ shades(color = x, ncolor = y, variation = 1)}, x= shades_num$phylum_col ,y = shades_num$unique_types,
                    SIMPLIFY = F)

order_col <- as.vector(unlist(order_col))
names(order_col) <- tax_sel$Order

mycol <- order_col

# Add a row for TTP
match_tax_dt[nrow(match_tax_dt)+1,] <- NA
rownames(match_tax_dt)[nrow(match_tax_dt)] <- "TTP"
match_tax_dt <- match_tax_dt[match(rownames(mat_imm),rownames(match_tax_dt)),]


ha1 = HeatmapAnnotation(Order = match_tax_dt$Order,
                        col = list(Order = order_col), which = "row")
split_rows <-  match_tax_dt$Phylum
split_rows <- factor(split_rows, levels= unique(split_rows))


# Replace dots with space
rownames(mat_imm)  <- gsub("\\."," ",rownames(mat_imm))
rownames(mat_imm)  <- gsub("X","",rownames(mat_imm))

colnames(mat_imm) <- gsub("HALLMARK_","",colnames(mat_imm))
colnames(mat_imm) <- gsub("_"," ",colnames(mat_imm))

# Find the index where t


ht1  = Heatmap(mat_imm,name = "Relation",
               column_split = sp_cols,
               row_split =  split_rows, row_gap = unit(2, "mm"),
               row_title_gp = gpar(fontsize = 10,fontface = "bold"),
               row_title_rot = 0,
               column_title_gp = gpar(fontsize = 10,fontface = "bold"),
               column_title_rot = 90,
               left_annotation = ha1,
               show_parent_dend_line = F,
               border = T,
               rect_gp = gpar(col = "gainsboro"),
               row_names_side = "left",
               show_row_dend = F,
               col = col_binary,
               row_names_max_width = max_text_width(rownames(mat_imm),
                                                    gp = gpar(fontsize = 12)),
               column_names_max_height = max_text_width(colnames(mat_imm),
                                                        gp = gpar(fontsize = 12)),
               heatmap_legend_param = list(
                 title = "Relation", at = c(-1, 1), 
                 labels = c("Neg", "Pos")
               ),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 #if(mat_imm[i, j] > 0)
                 for (k in 1:length(idx_r)){
                   if(i == idx_r[k] & j == idx_c[k]){
                     print(paste0(i," ",j))
                     grid.points(x, y, pch = 16, size = unit(4, "mm"))
                   }
                 }
               }
               # layer_fun = function(j, i, x, y, w, h, fill) {
               #   ind_mat = restore_matrix(j, i, x, y)
               #   
               #   
               #   ind = ind_mat[1, ]
               #   grid.points(x[ind], y[ind], pch = 16, size = unit(4, "mm"))
               # }
               
)


pdf(paste0(fig_folder,"/Heatmap_imp_Mic_imm_Pathways.pdf"),height = 10, width = 8)
draw(ht1)
dev.off()


# Draw some XY plot for top predictors 

sel_path <- immune_path[immune_path %in% lme_pathways]
sel_path_dt <- path_dt[,sel_path]

#Predictors to visualize
pred_vec <- c(unique(imm_var_sel$variable),"Escherichia_coli_.ASV_15.",
              "Enterococcus_faecium_.ASV_257.","Streptococcus_alactolyticus_.ASV_53.")


pred_dt <-  cbind(mic_dt,TTP = diff_ttp$Diff)
pred_dt<-  pred_dt[,colnames(pred_dt) %in% pred_vec]
pred_dt$sample_id <- rownames(pred_dt)
pred_dt_m <- reshape2::melt(pred_dt,id.vars = "sample_id")
names(pred_dt_m)[2:3] <-  c("Pred","X")


path_dt <-  path_dt[,colnames(path_dt) %in% sel_path]
path_dt$sample_id <- rownames(path_dt)
path_dt_m <- reshape2::melt(path_dt,id.vars = "sample_id")
names(path_dt_m)[2:3] <-  c("Path","Y")


XY_dt <- merge(path_dt_m,pred_dt_m,by = "sample_id")

pred <- unique(XY_dt$Pred)[1]

for(pred in unique(XY_dt$Pred)){
  sel_dt <- XY_dt[XY_dt$Pred == pred,]
  
  pred  <- gsub("\\."," ",pred)
  pred <- gsub("X","",pred)
  xy_p <- ggplot(sel_dt, aes(y = Y, x = X)) +
    geom_point(color = "#85796F",size = 3)+
    facet_wrap(~Path, scales = "free") +
    theme_bw() +
    ylab("Change in Pathways") +
    xlab("Change in Pred")+
    ggtitle(pred)
  
  pdf(paste0(fig_folder,"/XY_plot_imm_path_",pred,".pdf"), width = 12, height = 8)
  print(xy_p)
  dev.off()
}









