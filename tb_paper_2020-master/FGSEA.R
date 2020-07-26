library(tidyverse)
library(fgsea)
library(msigdbr)


# Clinical Trial (HRZE and NTZ
# Create a directory to save figures and tables

mainDir <- "../Figs"
subDir <- "FGSEA"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
fig_folder <- paste(mainDir,subDir,sep="/")

mainDir <- "../Tables"
subDir <- "FGSEA"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
tab_folder <- paste(mainDir,subDir,sep="/")



# MOdeling separately.
# Run FGSEA analysis 
trt <- c("HRZE","NTZ")[1]

pathway_list <- list()
for(trt in c("HRZE","NTZ")){
  # Limma results for genes
  res_trt <- read.csv(paste0("../Tables/Limma_rna_drug/",trt,"_Post_Vs_Pre_50.csv")
                      ,row.names = 1)
  res_TRT <- res_trt
  res_TRT$genes <- rownames(res_TRT)
  res_TRT$SYMBOL <- gsub(".*_","",rownames(res_trt))  
  
  library(tidyverse)
  res <- res_TRT
  
  rownames(res) <- gsub("\\..*","",rownames(res))
  res$row = rownames(res)
  
  res2 <- res %>% 
    dplyr::select(SYMBOL, t) %>% 
    na.omit() %>% 
    distinct() %>% 
    group_by(SYMBOL) %>% 
    summarize(stat=mean(t))
  res2
  
  library(tidyverse)
  library(fgsea)
  ranks <- deframe(res2)
  head(ranks, 20)
  library(msigdbr)
  db_sets <-  c("H","C2","C5","C6","C7")
  
  res_sum_trt_list <- list()
  
  p_set <- db_sets[1]
  # Pathway set 
  m_df = msigdbr("Homo sapiens") %>% dplyr::filter(gs_cat == p_set)
  m_t2g = m_df %>% dplyr::select(gs_name, entrez_gene,gene_symbol) %>% as.data.frame()
  
  # Select kegg,
  #m_t2g <- m_t2g[grep("KEGG",m_t2g$gs_name),]
  
  m_p_set <- m_t2g[,c("entrez_gene","gs_name")]
  
  #m_t2g$entrez_gene <- as.character(m_t2g$entrez_gene )
  pathway_set <- unstack(m_p_set)
  pathways.hallmark <- unstack(m_t2g[,c("gene_symbol","gs_name")])
  
  
  
  # Look at them all if you want (uncomment)
  pathways.hallmark
  
  # Show the first few pathways, and within those, show only the first few genes. 
  pathways.hallmark %>% 
    head() %>% 
    lapply(head)
  
  set.seed(3000)
  fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=20000,
                    minSize = 5, maxSize = 500)
  
  
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  fgseaResTidy$Group <- trt
  
  
  pathway_list[[trt]] <- fgseaResTidy
  
  # library(ggthemes)
  # p_hrze <- ggplot(fgseaResTidy[fgseaResTidy$padj<=0.01,], aes(reorder(pathway, NES), NES)) +
  #   geom_bar(stat = "identity", aes(fill = NES>0),width = 0.6) +
  #   scale_fill_manual(values = c("TRUE" = "blue","FALSE" = "red"), labels = c("Negative","Positive"))+
  #   coord_flip()+
  #   theme_base()+
  #   labs(x=" Hallmark Pathway", y="Normalized Enrichment Score",
  #        title=paste0("Significant pathways in ", " treatment")) +
  #   scale_y_continuous(limits=c(-5,5))+
  #   theme(axis.text.x=element_text(angle = -90, hjust = 0),legend.title = element_blank())
  # 
  # pdf(paste(fig_folder,paste0("Hallmark_pathways_from_GSEA_for_",trt,"_res_Fig_3_A.pdf"),sep = "/"),width=10,height=8)
  # print(p_hrze)
  # dev.off()
  # 
}  
# Draw an arrow shape instead of points

sig_path <- do.call("rbind",pathway_list)
sig_path$leadingEdge <-  NULL

write.csv(sig_path,paste0(tab_folder,"/FGSEA_clinical_tab.csv"))

# Import Pathways annotation
p_ann <-  read.csv("../Hallmark/Hallmark_pathway_annotation.csv",
                   stringsAsFactors = F)
p_ann$Hallmark.Name <-  paste0("HALLMARK_",p_ann$Hallmark.Name)
p_ann <- p_ann[p_ann$Hallmark.Name %in% as.character(sig_path$pathway),]
p_ann <- p_ann[match(sig_path$pathway,p_ann$Hallmark.Name),]
sig_path$Category <-  p_ann$Process.Category


sig_path <- sig_path[sig_path$padj< 0.01,]
sig_path$Regulation <-  ifelse(sig_path$NES>0,"Up","Down")
sig_path$pathway <- gsub("HALLMARK_","",sig_path$pathway)

sig_path <- sig_path[order(sig_path$Category),]
sig_path$pathway <- factor(sig_path$pathway)


library(RColorBrewer)
colfunc <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))
#col_dir <-  c("red","white","blue")
shape_vec <-  c("Up" = 24,"Down" = 25)
col_dir <- colfunc(7) 
p_hall <- ggplot(sig_path,aes(x = Group, y = pathway,fill = NES)) +
  geom_point(aes(shape = Regulation,size = -log10(padj)), color = "black", stroke = 1)+
  scale_shape_manual(values = shape_vec)+
  scale_fill_gradientn(limit = c(-4,4),colors = col_dir) +
  facet_wrap(~Category, strip.position = "top", scales = "free_y",ncol = 1)+
  theme_bw()+
  theme(axis.text.y=element_text(size = 10,face = "bold"))+
  theme(axis.text.x=element_text(size = 10,face = "bold"))+
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  theme(panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside")+
  ggtitle("Sig Hallmark Pathways across groups (padj < 0.01) ")+
  guides(size = guide_legend(override.aes = list(shape=24)))+
  ylab("Pathway")
#p_hall  

pdf(paste0(fig_folder,"/Hallmark_pathways_FGSEA_HRZE_NTZ.pdf"),width=8,height=14)
print(p_hall)
dev.off()



# EBA Analysis


# MOdeling separately.
# Run FGSEA analysis 
trt <- c("Time_D14_Vs_D0","Time_D56_Vs_D0","TimeD56_Vs_D14")[1]

limma_eba_dir <- "../Tables/EBA_analysis/Limma_rnaseq_eba/"

pathway_list <- list()
for(trt in c("Time_D14_Vs_D0","Time_D56_Vs_D0","TimeD56_Vs_D14")){
  # Limma results for genes
  res_trt <- read.csv(paste0(limma_eba_dir,trt,"_100.csv"),
                      row.names = 1)
  res_TRT <- res_trt
  res_TRT$genes <- rownames(res_TRT)
  res_TRT$SYMBOL <- gsub(".*_","",rownames(res_trt))  
  
  library(tidyverse)
  res <- res_TRT
  
  rownames(res) <- gsub("\\..*","",rownames(res))
  res$row = rownames(res)
  
  res2 <- res %>% 
    dplyr::select(SYMBOL, t) %>% 
    na.omit() %>% 
    distinct() %>% 
    group_by(SYMBOL) %>% 
    summarize(stat=mean(t))
  res2
  
  library(tidyverse)
  library(fgsea)
  ranks <- deframe(res2)
  head(ranks, 20)
  library(msigdbr)
  db_sets <-  c("H","C2","C5","C6","C7")
  
  res_sum_trt_list <- list()
  
  p_set <- db_sets[1]
  # Pathway set 
  m_df = msigdbr("Homo sapiens") %>% dplyr::filter(gs_cat == p_set)
  m_t2g = m_df %>% dplyr::select(gs_name, entrez_gene,gene_symbol) %>% as.data.frame()
  
  # Select kegg,
  #m_t2g <- m_t2g[grep("KEGG",m_t2g$gs_name),]
  
  m_p_set <- m_t2g[,c("entrez_gene","gs_name")]
  
  #m_t2g$entrez_gene <- as.character(m_t2g$entrez_gene )
  pathway_set <- unstack(m_p_set)
  pathways.hallmark <- unstack(m_t2g[,c("gene_symbol","gs_name")])
  
  
  
  # Look at them all if you want (uncomment)
  pathways.hallmark
  
  # Show the first few pathways, and within those, show only the first few genes. 
  pathways.hallmark %>% 
    head() %>% 
    lapply(head)
  
  set.seed(3000)
  fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=20000,
                    minSize = 5, maxSize = 500)
  
  
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  fgseaResTidy$Group <- trt
  pathway_list[[trt]] <- fgseaResTidy
  
}  
# Draw an arrow shape instead of points
sig_path <- do.call("rbind",pathway_list)
sig_path$leadingEdge <-  NULL

write.csv(sig_path,paste0(tab_folder,"/FGSEA_EBA_tab.csv"))

# Import Pathways annotation
p_ann <-  read.csv("../Hallmark/Hallmark_pathway_annotation.csv",
                   stringsAsFactors = F)
p_ann$Hallmark.Name <-  paste0("HALLMARK_",p_ann$Hallmark.Name)
p_ann <- p_ann[p_ann$Hallmark.Name %in% as.character(sig_path$pathway),]
p_ann <- p_ann[match(sig_path$pathway,p_ann$Hallmark.Name),]
sig_path$Category <-  p_ann$Process.Category


sig_path <- sig_path[sig_path$padj< 0.01,]
sig_path$Regulation <-  ifelse(sig_path$NES>0,"Up","Down")
sig_path$pathway <- gsub("HALLMARK_","",sig_path$pathway)

sig_path <- sig_path[order(sig_path$Category),]
sig_path$pathway <- factor(sig_path$pathway)

sig_path$Group <-  factor(sig_path$Group,levels = c("TimeD14","TimeD56_TimeD14","TimeD56"))

library(RColorBrewer)
colfunc <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))
#col_dir <-  c("red","white","blue")
shape_vec <-  c("Up" = 24,"Down" = 25)
col_dir <- colfunc(7) 
p_hall <- ggplot(sig_path,aes(x = Group, y = pathway,fill = NES)) +
  geom_point(aes(shape = Regulation,size = -log10(padj)), color = "black", stroke = 1)+
  scale_shape_manual(values = shape_vec)+
  scale_fill_gradientn(limit = c(-4,4),colors = col_dir) +
  facet_wrap(~Category, strip.position = "top", scales = "free_y",ncol = 1)+
  theme_bw()+
  theme(axis.text.y=element_text(size = 10,face = "bold"))+
  theme(axis.text.x=element_text(size = 10,face = "bold"))+
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside")+
  ggtitle("Sig Hallmark Pathways across groups (padj < 0.01) ")+
  guides(size = guide_legend(override.aes = list(shape=24)))
#p_hall  


pdf(paste0(fig_folder,"/Hallmark_pathways_FGSEA_EBA.pdf"),width=8,height=14)
print(p_hall)
dev.off()
