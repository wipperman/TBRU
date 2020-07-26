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
library(gtools)




# Create a directory to save figures and tables

mainDir <- "../data"
subDir <- "GSVA"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
fig_folder <- paste(mainDir,subDir,sep="/")



#RNASeq 
rna_phy <-  readRDS("../data/genes/phy_rna_all.rds")

#sample_names(rna_phy) <- gsub("s_","",sample_names(rna_phy))
#sample_names(rna_phy) <- gsub("_R1.*","",sample_names(rna_phy))
sm_dt <-  read.csv("../raw_data/Mic_data/metadata/TBRU_NTZ_RNAseq_metadata.csv")
sm_dt$rna_seq_ID <-  as.character(sm_dt$rna_seq_ID)
sm_dt <-  sm_dt[!sm_dt$rna_seq_ID == "",]
sm_dt <-  sm_dt[!sm_dt$sample == "",]
sm_dt <- sm_dt[!sm_dt$TB_status %in%c("cured","treatment"),]
rownames(sm_dt) <-  sm_dt$rna_seq_ID

# Only keep Patient.ID, Sex, TB_status
sm_dt <- sm_dt[,c("Patient.ID","sex","TB_status")]

sample_data(rna_phy) <-  sample_data(sm_dt)

rna_phy <-  prune_taxa(taxa_sums(rna_phy)>0, rna_phy)

phy_rna_eba <-  readRDS("../data/EBA/phy_rna_tb.rds")
sm_dt_eba <-  data.frame(sample_data(phy_rna_eba))

sm_dt_eba <- sm_dt_eba[,c("Patient.ID","SEX","GROUP")]
names(sm_dt_eba) <- c("Patient.ID","sex","TB_status")
sample_data(phy_rna_eba) <-  sample_data(sm_dt_eba)

ps_comb <-  merge_phyloseq(rna_phy,phy_rna_eba)



#tmp <- data.frame(sample_data(ps_comb))

phy_gene_sel <-  ps_comb


dds <- phyloseq_to_deseq2(phy_gene_sel , ~  TB_status) #replace this with any sample variable(s)

dds$TB_status
dds <- estimateSizeFactors(dds,"poscounts")


# Takes a while to run
dds <- estimateDispersions(dds)
dds <- DESeq(dds,fitType= "local")

# VST phyloseq for Heatmap later
phy_vst_rna <- phy_gene_sel
vst_dt <- getVarianceStabilizedData(dds)
otu_table(phy_vst_rna) <-  otu_table(vst_dt, taxa_are_rows = T)
saveRDS(phy_vst_rna,paste0(fig_folder,"/phy_rna_vst_tb_status_all.rds"))

#phy_vst_rna <-  readRDS(paste0("","../../EBA_HRZE_NTZ/Figs_vb/GSVA/phy_rna_vst_tb_status_all.rds"))
#phy_vst_rna <-  readRDS("../data/GSVA/phy_rna_vst_tb_status_all.rds")
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

# gsvaRes_ssgsea <- GSVA::gsva(counts, custom.gene.sets, 
#       method = "ssgsea", parallel.sz = 2) #can alter method here:("gsva", "ssgsea", "zscore", "plage")
# gsvaRes_ssgsea
set.seed(1057)
gsvaRes_ssgsea <- gsva(counts, hall_set,
                       min.sz=5, max.sz=500,
                       method = "gsva",
                       kcdf="Gaussian", mx.diff=TRUE, verbose=FALSE, parallel.sz=1
)

write.csv(gsvaRes_ssgsea,paste0(fig_folder,"/GSVA_all_Hallmark.csv"))
#write.csv(gsvaRes_ssgsea,paste0(fig_folder,"/GSVA_all.csv"))


phy_gsva <- phyloseq(otu_table(gsvaRes_ssgsea,taxa_are_rows = T),sample_data(metadata))
saveRDS(phy_gsva,paste0(fig_folder,"/phy_GSVA_all_hallmark.rds"))


met  <- metadata
mat  <- gsvaRes_ssgsea 

met$SampleID  <- rownames(metadata)
met <- met[match(colnames(gsvaRes_ssgsea),met$SampleID),]

#make phyloseq object out of gsva data
met_samp <- met
met_samp$sample  <- rownames(met_samp)
rownames(met_samp) <- NULL
gvsa_phy <- phyloseq(otu_table(mat,taxa_are_rows = T),set.samp(met_samp))
gvsa_phy

status_col <-  brewer.pal(length(unique(met$TB_status)), "Set1")
names(status_col) <- unique(unique(met$TB_status))

mixedsort(met$TB_status)
ha_column = HeatmapAnnotation(Status =  met$TB_status,
                              col=list(Status = status_col))

split_cols<-  met$TB_status
split_cols <- factor(split_cols, 
                     levels= c("pretreatment","EBA_0","NTZ","HRZE","EBA_14","EBA_56",
                               "community_control","family_contact"))

colfunc <- colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))
cols <- colfunc(8)
rownames(mat) <-  gsub("HALLMARK_","",rownames(mat))

library(tibble)
Hallmark_pathway_info <- read.csv("../Hallmark/Hallmark_pathway_annotation.csv")
rownames(Hallmark_pathway_info) <- Hallmark_pathway_info$Hallmark.Name

setdiff(rownames(mat) ,rownames(Hallmark_pathway_info))

Hallmark_pathway_info <-  Hallmark_pathway_info[match(rownames(mat),
                                                      rownames(Hallmark_pathway_info)),]
split_rows <- Hallmark_pathway_info$Process.Category

rownames(mat) <- gsub("_"," ",rownames(mat))
rownames(mat) <- str_to_sentence(rownames(mat))


ht1 = Heatmap(mat, name = "NES", column_title = NA, 
              top_annotation = ha_column,
              col = cols,
              cluster_column_slices = F,
              row_split = split_rows,
              row_names_side = "left",
              row_title_gp = gpar(fontsize = 10,fontface = "bold"), 
              row_title_rot = 0,
              column_split = split_cols,
              show_parent_dend_line = F,
              width=2, cluster_columns = T, 
              row_names_gp = gpar(fontsize = 9),
              row_names_max_width = max_text_width(rownames(mat), gp = gpar(fontsize = 12)),
              show_column_names = F, show_row_names = T,
              column_names_side = "bottom",na_col="white",
              show_row_dend = F)
#pdf("GSVA_raw.pdf",width = 14,height = 8)
pdf(paste0(fig_folder,"/GSVA_deseq_norm_vst.pdf"),width = 17,height = 9)
draw(ht1)
dev.off()
