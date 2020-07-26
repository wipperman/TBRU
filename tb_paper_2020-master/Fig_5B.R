#Load libraries
library(phyloseq)
library(dplyr)
library(yingtools2)
library(data.table)
library(ifultools)
library(stringr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(variancePartition)
library(edgeR)
library(BiocParallel)
library(ggthemes)

# Create a directory to save figures and tables
mainDir <- "../Figs/EBA_analysis"
subDir <- "Limma_rnaseq_eba"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
fig_folder <- paste(mainDir,subDir,sep="/")

mainDir <- "../Tables/EBA_analysis"
subDir <- "Limma_rnaseq_eba"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
tab_folder <- paste(mainDir,subDir,sep="/")


####### load phyloseq objected created with data_setup.R ######

phy_pair <- readRDS("../data/EBA/phy_rna_tb.rds")



###### Differential Analysis with LIMMA ########
phy_trt <-  phy_pair
sm_dt <- data.frame(sample_data(phy_trt))
names(sm_dt)
sm_dt <- sm_dt %>%
  arrange(Patient.ID,DAYS_ON_TREATMENT) %>%
  as.data.frame()

sm_dt$ind <- as.numeric(factor(as.character(sm_dt$Patient.ID),
                               levels = unique(as.character(sm_dt$Patient.ID))))
rownames(sm_dt) <- sm_dt$sam_id

sample_data(phy_trt) <- sample_data(sm_dt)

table(get_variable(phy_trt,"Time"))
sample_data(phy_trt)$ind
sample_data(phy_trt)$SEX
sample_data(phy_trt)$EXTRACTION
# For each grp separately

sum_counts <- c(10,seq(100,500,100))

# Filter counts that the total counts across all samples is >= cnt_filt

cnt_filt <- 100

phy_sel <- phy_trt
#phy_sel <- subset_samples(phy_sel,Time %in% c("D0","D14"))
phy_sel <-  prune_taxa(taxa_sums(phy_sel)>0,phy_sel)


# Differential analysis via Limma
cnt_dt <- as.matrix(data.frame(otu_table(phy_sel)))
colnames(cnt_dt) <- gsub("X","",colnames(cnt_dt))
sm_dt <-  data.frame(sample_data(phy_sel))
sm_dt$Time <- factor(sm_dt$Time,levels = c("D0", "D14","D56"))
sm_dt$ind <- factor(sm_dt$ind)
sm_dt$EXTRACTION
sm_dt$SEX
cnt_dt <-  cnt_dt[,match(rownames(sm_dt),colnames(cnt_dt))]

dim(cnt_dt)
# filter genes by number of counts
#table(rowSums(cpm(cnt_dt)>0.1)==8)

# Filter out the lowly expressed genes with the total rowSum >= cnt_filt
cnt_filt
sum(rowSums(cnt_dt) >= cnt_filt)
isexpr = rowSums(cnt_dt) >= cnt_filt
table(rowSums(cnt_dt) >=cnt_filt)
# Standard usage of limma/voom
gExpr = DGEList( cnt_dt[isexpr,] )
gExpr = calcNormFactors( gExpr )
dim(gExpr$counts)

sm_dt$Time
sm_dt$ind
sm_dt$SEX


#Visualize
library(RColorBrewer)
lcpm <- cpm(gExpr$counts, log=TRUE)
group <- paste0(sm_dt$SEX)
plotMDS(lcpm, labels=group, dim.plot = c(1,3),top = 500)
group <- paste0(sm_dt$Time)
plotMDS(lcpm, labels=group, dim.plot = c(2,4),top = 500)
group <- paste0(sm_dt$Time,"_",sm_dt$ind)
plotMDS(lcpm,labels = group, dim.plot =  c(2,4))


# Limma 
# apply duplicateCorrelation is two rounds
design = model.matrix( ~ SEX + EXTRACTION + Time  , sm_dt)
# apply duplicateCorrelation is two rounds
vobj_tmp = voom(gExpr, design, plot=FALSE)
dupcor <- duplicateCorrelation(vobj_tmp,design,block=sm_dt$ind)

# run voom considering the duplicateCorrelation results
# in order to compute more accurate precision weights
# Otherwise, use the results from the first voom run
pdf(paste0(fig_folder,"/Mean_Dispersion_plot_EBA.pdf"),width = 5,height = 5)
vobj = voom( gExpr, design, plot=T, block=sm_dt$ind, correlation=dupcor$consensus)
dev.off()
# Estimate linear mixed model with a single variance component
# Fit the model for each gene, 
dupcor <- duplicateCorrelation(vobj, design, block=sm_dt$ind)

# But this step uses only the genome-wide average for the random effect
fitDupCor <- lmFit(vobj, design, block=sm_dt$ind, correlation=dupcor$consensus)

# Fit Empirical Bayes for moderated t-statistics
fitDupCor <- eBayes( fitDupCor )

# get names of available coefficients and contrasts for testing
colnames(fitDupCor)


topTable( fitDupCor, coef=c("TimeD14"),n=3)
topTable( fitDupCor, coef=c("TimeD56"),n=3)

topTable( fitDupCor, coef=c("SEXMale"), number=3 )
topTable( fitDupCor, coef=c("EXTRACTIONTempus"), number=3 )
p_cut_off <-  0.05 

cont.mat <- cbind(Time_14_Vs_0 = c(0,0,0,1,0) ,
                  Time_56_Vs_0 = c(0,0,0,0,1) ,
                  TimeD56_Vs_TimeD14 = c(0,0,0,-1,1),
                  TimeAll_Vs_0 = c(0,0,0,1,1),
                  Extraction_Tempus_Vs_Paxgene = c(0,0,1,0,0),
                  Sex_M_Vs_F =  c(0,1,0,0,0)
)

final_fit <-  contrasts.fit(fitDupCor, cont.mat)
final_fit <- eBayes(final_fit)
sum_eba <- summary( decideTests(final_fit,p.value = p_cut_off,
                                method = "global",adjust.method = "BH"))
sum_dt <-  as.data.frame.matrix(sum_eba) 
write.csv(sum_dt ,paste(tab_folder,paste("Summary_EBA_Time_comp_genes.csv",sep=""),sep = "/"))



# Global pval adjustment
nom_pval <-  final_fit$p.value
rnames <- rownames(nom_pval)

vec_pval <- c(nom_pval)
global_pval <- p.adjust(nom_pval,method = "BH")
global_mat <- data.frame(matrix(global_pval,ncol = ncol(nom_pval)))
names(global_mat) <- colnames(nom_pval)
rownames(global_mat) <-  rnames


# DEG D0 Vs D14 
res_14_Vs_0 <- topTable( final_fit, coef=c("Time_14_Vs_0"),n=Inf)
#summary( decideTests(final_fit,coef = c("Time_14_Vs_0"),p.value = p_cut_off))
p_cut_off <-  0.05
res_14_Vs_0 <- data.frame(res_14_Vs_0)

# Match the order with global_mat
res_14_Vs_0 <- res_14_Vs_0[match(rownames(global_mat),rownames(res_14_Vs_0)),]

res_14_Vs_0$global_padj <- global_mat$Time_14_Vs_0
res_14_Vs_0 <- res_14_Vs_0[!is.na(res_14_Vs_0$adj.P.Val),]
res_14_Vs_0 <- res_14_Vs_0[order(res_14_Vs_0$global_padj,decreasing = F),]
#res_14_Vs_0_sig <- res_14_Vs_0[res_14_Vs_0$adj.P.Val <= p_cut_off,]
write.csv(res_14_Vs_0,paste(tab_folder,paste("Time_D14_Vs_D0","_",cnt_filt,".csv",sep=""),sep = "/"))


# Saving results\
# DEG D0 Vs D56 
res_56_Vs_0 <- topTable( final_fit, coef=c("Time_56_Vs_0"),n=Inf)
p_cut_off <-  0.05
res_56_Vs_0 <- data.frame(res_56_Vs_0)

# Match the order with global_mat
res_56_Vs_0 <- res_56_Vs_0[match(rownames(global_mat),rownames(res_56_Vs_0)),]

res_56_Vs_0$global_padj <-  global_mat$Time_56_Vs_0
res_56_Vs_0 <- res_56_Vs_0[!is.na(res_56_Vs_0$global_padj),]
res_56_Vs_0 <- res_56_Vs_0[order(res_56_Vs_0$global_padj,decreasing = F),]

#res_56_Vs_0_sig <- res_56_Vs_0[res_56_Vs_0$adj.P.Val <= p_cut_off,]
write.csv(res_56_Vs_0,paste(tab_folder,paste("Time_D56_Vs_D0_",cnt_filt,".csv",sep=""),sep = "/"))


# DEG D56 Vs D14 
colnames(final_fit)
res_56_Vs_14 <- topTable( final_fit, coef=c("TimeD56_Vs_TimeD14"),n=Inf)
p_cut_off <-  0.05
res_56_Vs_14 <- data.frame(res_56_Vs_14)
# Match the order with global_mat
res_56_Vs_14 <- res_56_Vs_14[match(rownames(global_mat),rownames(res_56_Vs_14)),]

res_56_Vs_14$global_padj <-  global_mat$TimeD56_Vs_TimeD14
res_56_Vs_14 <- res_56_Vs_14[!is.na(res_56_Vs_14$adj.P.Val),]
res_56_Vs_14 <- res_56_Vs_14[order(res_56_Vs_14$global_padj,decreasing = F),]
#res_56_Vs_14_sig <- res_56_Vs_14[res_56_Vs_14$adj.P.Val <= p_cut_off,]
write.csv(res_56_Vs_14,paste(tab_folder,paste("TimeD56_Vs_D14_",cnt_filt,".csv",sep=""),sep = "/"))



# Differences across Male/Female 
res_sex <- topTable( final_fit, coef=c("Sex_M_Vs_F"), n=Inf )
p_cut_off <-  0.05
res_sex <- data.frame(res_sex)
# Match the order with global_mat
res_sex <- res_sex[match(rownames(global_mat),rownames(res_sex)),]

# Adjusted global p val
res_sex$global_padj <- global_mat$Sex_M_Vs_F
res_sex <- res_sex[!is.na(res_sex$adj.P.Val),]
res_sex <- res_sex[order(res_sex$global_padj,decreasing = F),]

#res_sex_sig <- res_sex[res_sex$adj.P.Val <= p_cut_off,]
write.csv(res_sex,paste(tab_folder,paste("Sex_M_Vs_F_",cnt_filt,".csv",sep=""),sep = "/"))


# Differences across Extraction Extraction_Tempus_Vs_Paxgene
res_ext <- topTable( final_fit, coef=c("Extraction_Tempus_Vs_Paxgene"), n=Inf )
p_cut_off <-  0.05
res_ext <- data.frame(res_ext)
# Match the order with global_mat
res_ext <- res_ext[match(rownames(global_mat),rownames(res_ext)),]

# Adjusted global p val
res_ext$global_padj <- global_mat$Extraction_Tempus_Vs_Paxgene
res_ext <- res_ext[!is.na(res_ext$adj.P.Val),]
res_ext <- res_ext[order(res_ext$global_padj,decreasing = F),]

#res_sex_sig <- res_sex[res_sex$adj.P.Val <= p_cut_off,]
write.csv(res_ext,paste(tab_folder,paste("Extraction_Tempus_Vs_Paxgene_",cnt_filt,".csv",sep=""),sep = "/"))





# Genes that responds to treatment at any time
# DEG (D56 and D14)  Vs D0
res_all_Vs_0 <- topTable( final_fit,coef=c("TimeAll_Vs_0"),n=Inf)
p_cut_off <-  0.05
res_all_Vs_0 <- data.frame(res_all_Vs_0)

# Match the order with global_mat
res_all_Vs_0 <- res_all_Vs_0[match(rownames(global_mat),rownames(res_all_Vs_0)),]

res_all_Vs_0$global_padj <-  global_mat$TimeAll_Vs_0
res_all_Vs_0 <- res_all_Vs_0[!is.na(res_all_Vs_0$adj.P.Val),]
res_all_Vs_0 <- res_all_Vs_0[order(res_all_Vs_0$global_padj,decreasing = F),]
#res_all_Vs_0_sig <- res_all_Vs_0[res_all_Vs_0$adj.P.Val <= p_cut_off,]

write.csv(res_all_Vs_0,paste(tab_folder,paste("Time_All_Vs_D0_",cnt_filt,".csv",sep=""),sep = "/"))


# Map Berry gene set in HRZE group

berry_data <- read.csv("../Berry/HRZE_NTZ_Berry_List_April_2020.csv")
ibd_data <- read.csv("../IBD/IBD_sig_genes.csv")


#  Berry subset
# idea is to map the TB related genes in HRZE and NTZ data

# Mapping to HRZE Set
res_hrze <- res_14_Vs_0
res_hrze$Gene <-  rownames(res_hrze)
res_hrze$Ensembl_ID <-  gsub("\\..*","",res_hrze$Gene)

# Merge the HRZE results with berry signature
berry_hrze <-  merge(berry_data,res_hrze,by.x=  "Ensembl_ID", by.y = "Ensembl_ID")

# Active Vs Control
# Adjusting direction with our data
berry_hrze$logFC <-  -1*berry_hrze$logFC

berry_hrze$Comb_ID <- as.character(berry_hrze$Comb_ID)
# Number of genes that are significantly different in HRZE Pre and Post 
# And also TB related 
berry_hrze_sig <- berry_hrze[berry_hrze$global_padj < 0.05,] 
nrow(berry_hrze_sig)

# Comparison of LFC
# Both decreasing
down_both_dt <- berry_hrze[berry_hrze$logFC < 0 & 
                             berry_hrze$Log_2_Fold.change._Active_Control < 0,  ]
# Both increasing
up_both_dt <- berry_hrze[berry_hrze$logFC > 0 & 
                           berry_hrze$Log_2_Fold.change._Active_Control > 0,  ]

mismatch_dt <- berry_hrze[!berry_hrze$Comb_ID %in%
                            c(down_both_dt$Comb_ID, up_both_dt$Comb_ID),]



# IBD Set

# Map IBD set to HRZE group

ibd_data <- read.csv("../IBD/IBD_sig_genes.csv")
ibd_data$Gene.Symbol <- trimws(as.character(ibd_data$Gene.Symbol))
#  IBD subset
# idea is to map the IBD  genes in HRZE and NTZ data

# Mapping to HRZE Set
res_hrze <- res_14_Vs_0
res_hrze$Gene <-  rownames(res_hrze)
res_hrze$Ensembl_ID <-  gsub("\\..*","",res_hrze$Gene)
res_hrze$symbol <-  gsub(".*_","",res_hrze$Gene)

# Merge the HRZE results with berry signature
ibd_hrze <-  merge(ibd_data,res_hrze,by.x=  "Gene.Symbol", by.y = "symbol")

# Number of genes that are significantly different in HRZE Pre and Post 
# And also IBD related 
ibd_hrze_sig <- ibd_hrze[ibd_hrze$global_padj < 0.05,] 
nrow(ibd_hrze_sig)


# hrze berry
mycol <-c("red","#005073")
names(mycol) <-c("Exacerbate","Renormalize")
fc_dt <- berry_hrze
fc_dt$renorm <- "Exacerbate"
fc_dt$renorm[fc_dt$logFC < 0 & fc_dt$Log_2_Fold.change._Active_Control < 0]<-"Renormalize"
fc_dt$renorm[fc_dt$logFC > 0 & fc_dt$Log_2_Fold.change._Active_Control > 0]<-"Renormalize"
fc_dt$otu <- rownames(fc_dt)
mer_hrze <- fc_dt

# Adding excel mapped names for superpathway
data_pvalue_less <- mer_hrze[which(mer_hrze$global_padj<=0.05),]
hrze_berryp<-ggplot() + 
  geom_point(data=mer_hrze[which(mer_hrze$global_padj>0.05),],
             aes(y=-log10(global_padj),x=-logFC),
             size=2,color="black",alpha=0.1)+
  geom_point(data= data_pvalue_less,
             aes(y=-log10(global_padj),x=-logFC,color=renorm),
             size=2,alpha=.9)+
  scale_color_manual(name="",values = mycol) +
  theme_base()+
  xlab("Log2 Fold Change (Day 14/Day 0)")+
  ylab("-Log10 FDR")+
  xlim(-4,2)+
  ylim(0,8)
hrze_berryp

pdf(paste(fig_folder,"rnaseq_volcano_plot_berry_hrze_eba_Day_14_Day_0.pdf",sep = "/"),height = 5, width = 6)
print(hrze_berryp)
dev.off()

# # plot FC this study vs FC Berry
# hrze_berry_fcp<-ggplot() + 
#   geom_point(data=mer_hrze[which(mer_hrze$global_padj>0.05),],
#              aes(y=-logFC,x=-1*Log_2_Fold.change._Active_Control),
#              size=2,color="black",alpha=0.1)+
#   geom_point(data=mer_hrze[which(mer_hrze$global_padj<=0.05),],
#              aes(y=-logFC,x=-1*Log_2_Fold.change._Active_Control,
#                  color=renorm),
#              size=2,alpha=.9)+
#   scale_color_manual(name="",values = mycol) +
#   theme_base()+
#   xlab("Log2 Fold Change Control/Active (Berry)")+
#   ylab("Log2 Fold Change Post [Any] /Pre (This Study)")
# 
# hrze_berry_fcp
# pdf(paste(fig_folder,"rnaseq_fc_berry_fc_hrze_eba.pdf",sep = "/"),height = 5, width = 6)
# print(hrze_berry_fcp)
# dev.off()

berry_hrze_sig <- berry_hrze[berry_hrze$global_padj < 0.05,] 
berry_hrze_sig
write.csv(berry_hrze_sig, paste(tab_folder,paste("EBA_Berry_HRZE_differential_Day_14_Day_0",".csv",sep=""),sep = "/"))

# ibd hrze eba
fc_dt <- ibd_hrze
fc_dt$renorm <- "Exacerbate"
fc_dt$renorm[fc_dt$logFC < 0 & fc_dt$Case.Over.Control.Fold.Change < 0]<-"Renormalize"
fc_dt$renorm[fc_dt$logFC > 0 & fc_dt$Case.Over.Control.Fold.Change > 0]<-"Renormalize"
fc_dt$otu <- rownames(fc_dt)
mer_hrze <- fc_dt

# Adding excel mapped names for superpathway
data_pvalue_less <- mer_hrze[which(mer_hrze$global_padj<=0.05),]
hrze_ibdp<-ggplot() + 
  geom_point(data=mer_hrze[which(mer_hrze$global_padj>0.05),],
             aes(y=-log10(global_padj),x=-logFC),
             size=2,color="black",alpha=0.1)+
  geom_point(data= data_pvalue_less,
             aes(y=-log10(global_padj),x=-logFC,color=renorm),
             size=2,alpha=.9)+
  scale_color_manual(name="",values = mycol) +
  theme_base()+
  xlab("Log2 Fold Change (Day 14/Day 0)")+
  ylab("-Log10 FDR")+
  xlim(-3,3)+
  ylim(0,8)
hrze_ibdp

pdf(paste(fig_folder,"rnaseq_volcano_plot_ibd_hrze_eba_Day_14_Day_0.pdf",sep = "/"),height = 5, width = 6)
print(hrze_ibdp)
dev.off()

# plot FC this study vs FC Berry
# hrze_ibd_fcp<-ggplot() + 
#   geom_point(data=mer_hrze[which(mer_hrze$global_padj>0.05),],
#              aes(y=-logFC,x=-1*Case.Over.Control.Fold.Change),
#              size=2,color="black",alpha=0.1)+
#   geom_point(data=mer_hrze[which(mer_hrze$global_padj<=0.05),],
#              aes(y=-logFC,x=-1*Case.Over.Control.Fold.Change,
#                  color=renorm),
#              size=2,alpha=.9)+
#   scale_color_manual(name="",values = mycol) +
#   theme_base()+
#   xlab("Log2 Fold Change Case/Control (IBD)")+
#   ylab("Log2 Fold Change Post/Pre (This Study)")
# 
# hrze_ibd_fcp
# pdf(paste(fig_folder,"rnaseq_fc_ibd_fc_hrze_eba.pdf",sep = "/"),height = 5, width = 6)
# print(hrze_ibd_fcp)
# dev.off()

ibd_hrze_sig <- ibd_hrze[ibd_hrze$global_padj < 0.05,] 
nrow(ibd_hrze_sig)
write.csv(ibd_hrze_sig, paste(tab_folder,paste("EBA_IBD_HRZE_differential_Day_14_Day_0",".csv",sep=""),sep = "/"))



# Day 56:

# Mapping to HRZE Set
res_hrze <- res_56_Vs_0
res_hrze$Gene <-  rownames(res_hrze)
res_hrze$Ensembl_ID <-  gsub("\\..*","",res_hrze$Gene)

# Merge the HRZE results with berry signature
berry_hrze <-  merge(berry_data,res_hrze,by.x=  "Ensembl_ID", by.y = "Ensembl_ID")

# Active Vs Control
# Adjusting direction with our data
berry_hrze$logFC <-  -1*berry_hrze$logFC

berry_hrze$Comb_ID <- as.character(berry_hrze$Comb_ID)
# Number of genes that are significantly different in HRZE Pre and Post 
# And also TB related 
berry_hrze_sig <- berry_hrze[berry_hrze$global_padj < 0.05,] 
nrow(berry_hrze_sig)

# Comparison of LFC
# Both decreasing
down_both_dt <- berry_hrze[berry_hrze$logFC < 0 & 
                             berry_hrze$Log_2_Fold.change._Active_Control < 0,  ]
# Both increasing
up_both_dt <- berry_hrze[berry_hrze$logFC > 0 & 
                           berry_hrze$Log_2_Fold.change._Active_Control > 0,  ]

mismatch_dt <- berry_hrze[!berry_hrze$Comb_ID %in%
                            c(down_both_dt$Comb_ID, up_both_dt$Comb_ID),]



# IBD Set

# Map IBD set to HRZE group

ibd_data <- read.csv("../IBD/IBD_sig_genes.csv")
ibd_data$Gene.Symbol <- trimws(as.character(ibd_data$Gene.Symbol))
#  IBD subset
# idea is to map the IBD  genes in HRZE and NTZ data

# Mapping to HRZE Set
res_hrze <- res_56_Vs_0
res_hrze$Gene <-  rownames(res_hrze)
res_hrze$Ensembl_ID <-  gsub("\\..*","",res_hrze$Gene)
res_hrze$symbol <-  gsub(".*_","",res_hrze$Gene)

# Merge the HRZE results with berry signature
ibd_hrze <-  merge(ibd_data,res_hrze,by.x=  "Gene.Symbol", by.y = "symbol")

# Number of genes that are significantly different in HRZE Pre and Post 
# And also IBD related 
ibd_hrze_sig <- ibd_hrze[ibd_hrze$global_padj < 0.05,] 
nrow(ibd_hrze_sig)


# hrze berry
mycol <-c("red","#005073")
names(mycol) <-c("Exacerbate","Renormalize")
fc_dt <- berry_hrze
fc_dt$renorm <- "Exacerbate"
fc_dt$renorm[fc_dt$logFC < 0 & fc_dt$Log_2_Fold.change._Active_Control < 0]<-"Renormalize"
fc_dt$renorm[fc_dt$logFC > 0 & fc_dt$Log_2_Fold.change._Active_Control > 0]<-"Renormalize"
fc_dt$otu <- rownames(fc_dt)
mer_hrze <- fc_dt

# Adding excel mapped names for superpathway
data_pvalue_less <- mer_hrze[which(mer_hrze$global_padj<=0.05),]
hrze_berryp<-ggplot() + 
  geom_point(data=mer_hrze[which(mer_hrze$global_padj>0.05),],
             aes(y=-log10(global_padj),x=-logFC),
             size=2,color="black",alpha=0.1)+
  geom_point(data= data_pvalue_less,
             aes(y=-log10(global_padj),x=-logFC,color=renorm),
             size=2,alpha=.9)+
  scale_color_manual(name="",values = mycol) +
  theme_base()+
  xlab("Log2 Fold Change (Day 56/Day 0)")+
  ylab("-Log10 FDR")+
  xlim(-4,2)+
  ylim(0,8)
hrze_berryp

pdf(paste(fig_folder,"rnaseq_volcano_plot_berry_hrze_eba_Day_56_Day_0.pdf",sep = "/"),height = 5, width = 6)
print(hrze_berryp)
dev.off()

# # plot FC this study vs FC Berry
# hrze_berry_fcp<-ggplot() + 
#   geom_point(data=mer_hrze[which(mer_hrze$global_padj>0.05),],
#              aes(y=-logFC,x=-1*Log_2_Fold.change._Active_Control),
#              size=2,color="black",alpha=0.1)+
#   geom_point(data=mer_hrze[which(mer_hrze$global_padj<=0.05),],
#              aes(y=-logFC,x=-1*Log_2_Fold.change._Active_Control,
#                  color=renorm),
#              size=2,alpha=.9)+
#   scale_color_manual(name="",values = mycol) +
#   theme_base()+
#   xlab("Log2 Fold Change Control/Active (Berry)")+
#   ylab("Log2 Fold Change Post [Any] /Pre (This Study)")
# 
# hrze_berry_fcp
# pdf(paste(fig_folder,"rnaseq_fc_berry_fc_hrze_eba.pdf",sep = "/"),height = 5, width = 6)
# print(hrze_berry_fcp)
# dev.off()

berry_hrze_sig <- berry_hrze[berry_hrze$global_padj < 0.05,] 
berry_hrze_sig
write.csv(berry_hrze_sig, paste(tab_folder,paste("EBA_Berry_HRZE_differential_Day_56_Day_0",".csv",sep=""),sep = "/"))

# ibd hrze eba
fc_dt <- ibd_hrze
fc_dt$renorm <- "Exacerbate"
fc_dt$renorm[fc_dt$logFC < 0 & fc_dt$Case.Over.Control.Fold.Change < 0]<-"Renormalize"
fc_dt$renorm[fc_dt$logFC > 0 & fc_dt$Case.Over.Control.Fold.Change > 0]<-"Renormalize"
fc_dt$otu <- rownames(fc_dt)
mer_hrze <- fc_dt

# Adding excel mapped names for superpathway
data_pvalue_less <- mer_hrze[which(mer_hrze$global_padj<=0.05),]
hrze_ibdp<-ggplot() + 
  geom_point(data=mer_hrze[which(mer_hrze$global_padj>0.05),],
             aes(y=-log10(global_padj),x=-logFC),
             size=2,color="black",alpha=0.1)+
  geom_point(data= data_pvalue_less,
             aes(y=-log10(global_padj),x=-logFC,color=renorm),
             size=2,alpha=.9)+
  scale_color_manual(name="",values = mycol) +
  theme_base()+
  xlab("Log2 Fold Change (Day 56/Day 0)")+
  ylab("-Log10 FDR")+
  xlim(-3,3)+
  ylim(0,8)
hrze_ibdp

pdf(paste(fig_folder,"rnaseq_volcano_plot_ibd_hrze_eba_Day_56_Day_0.pdf",sep = "/"),height = 5, width = 6)
print(hrze_ibdp)
dev.off()

# plot FC this study vs FC Berry
# hrze_ibd_fcp<-ggplot() + 
#   geom_point(data=mer_hrze[which(mer_hrze$global_padj>0.05),],
#              aes(y=-logFC,x=-1*Case.Over.Control.Fold.Change),
#              size=2,color="black",alpha=0.1)+
#   geom_point(data=mer_hrze[which(mer_hrze$global_padj<=0.05),],
#              aes(y=-logFC,x=-1*Case.Over.Control.Fold.Change,
#                  color=renorm),
#              size=2,alpha=.9)+
#   scale_color_manual(name="",values = mycol) +
#   theme_base()+
#   xlab("Log2 Fold Change Case/Control (IBD)")+
#   ylab("Log2 Fold Change Post/Pre (This Study)")
# 
# hrze_ibd_fcp
# pdf(paste(fig_folder,"rnaseq_fc_ibd_fc_hrze_eba.pdf",sep = "/"),height = 5, width = 6)
# print(hrze_ibd_fcp)
# dev.off()

ibd_hrze_sig <- ibd_hrze[ibd_hrze$global_padj < 0.05,] 
nrow(ibd_hrze_sig)
write.csv(ibd_hrze_sig, paste(tab_folder,paste("EBA_IBD_HRZE_differential_Day_56_Day_0",".csv",sep=""),sep = "/"))


# Day 56 Vs Day 14

# Mapping to HRZE Set
res_hrze <- res_56_Vs_14
res_hrze$Gene <-  rownames(res_hrze)
res_hrze$Ensembl_ID <-  gsub("\\..*","",res_hrze$Gene)

# Merge the HRZE results with berry signature
berry_hrze <-  merge(berry_data,res_hrze,by.x=  "Ensembl_ID", by.y = "Ensembl_ID")

# Active Vs Control
# Adjusting direction with our data
berry_hrze$logFC <-  -1*berry_hrze$logFC

berry_hrze$Comb_ID <- as.character(berry_hrze$Comb_ID)
# Number of genes that are significantly different in HRZE Pre and Post 
# And also TB related 
berry_hrze_sig <- berry_hrze[berry_hrze$global_padj < 0.05,] 
nrow(berry_hrze_sig)

# Comparison of LFC
# Both decreasing
down_both_dt <- berry_hrze[berry_hrze$logFC < 0 & 
                             berry_hrze$Log_2_Fold.change._Active_Control < 0,  ]
# Both increasing
up_both_dt <- berry_hrze[berry_hrze$logFC > 0 & 
                           berry_hrze$Log_2_Fold.change._Active_Control > 0,  ]

mismatch_dt <- berry_hrze[!berry_hrze$Comb_ID %in%
                            c(down_both_dt$Comb_ID, up_both_dt$Comb_ID),]



# IBD Set

# Map IBD set to HRZE group

ibd_data <- read.csv("../IBD/IBD_sig_genes.csv")
ibd_data$Gene.Symbol <- trimws(as.character(ibd_data$Gene.Symbol))
#  IBD subset
# idea is to map the IBD  genes in HRZE and NTZ data

# Mapping to HRZE Set
res_hrze <- res_56_Vs_14
res_hrze$Gene <-  rownames(res_hrze)
res_hrze$Ensembl_ID <-  gsub("\\..*","",res_hrze$Gene)
res_hrze$symbol <-  gsub(".*_","",res_hrze$Gene)

# Merge the HRZE results with berry signature
ibd_hrze <-  merge(ibd_data,res_hrze,by.x=  "Gene.Symbol", by.y = "symbol")

# Number of genes that are significantly different in HRZE Pre and Post 
# And also IBD related 
ibd_hrze_sig <- ibd_hrze[ibd_hrze$global_padj < 0.05,] 
nrow(ibd_hrze_sig)


# hrze berry
mycol <-c("red","#005073")
names(mycol) <-c("Exacerbate","Renormalize")
fc_dt <- berry_hrze
fc_dt$renorm <- "Exacerbate"
fc_dt$renorm[fc_dt$logFC < 0 & fc_dt$Log_2_Fold.change._Active_Control < 0]<-"Renormalize"
fc_dt$renorm[fc_dt$logFC > 0 & fc_dt$Log_2_Fold.change._Active_Control > 0]<-"Renormalize"
fc_dt$otu <- rownames(fc_dt)
mer_hrze <- fc_dt

# Adding excel mapped names for superpathway
data_pvalue_less <- mer_hrze[which(mer_hrze$global_padj<=0.05),]
hrze_berryp<-ggplot() + 
  geom_point(data=mer_hrze[which(mer_hrze$global_padj>0.05),],
             aes(y=-log10(global_padj),x=-logFC),
             size=2,color="black",alpha=0.1)+
  geom_point(data= data_pvalue_less,
             aes(y=-log10(global_padj),x=-logFC,color=renorm),
             size=2,alpha=.9)+
  scale_color_manual(name="",values = mycol) +
  theme_base()+
  xlab("Log2 Fold Change (Day 56/Day 14)")+
  ylab("-Log10 FDR")+
  xlim(-4,2)+
  ylim(0,8)
hrze_berryp

pdf(paste(fig_folder,"rnaseq_volcano_plot_berry_hrze_eba_Day_56_Day_14.pdf",sep = "/"),height = 5, width = 6)
print(hrze_berryp)
dev.off()

# # plot FC this study vs FC Berry
# hrze_berry_fcp<-ggplot() + 
#   geom_point(data=mer_hrze[which(mer_hrze$global_padj>0.05),],
#              aes(y=-logFC,x=-1*Log_2_Fold.change._Active_Control),
#              size=2,color="black",alpha=0.1)+
#   geom_point(data=mer_hrze[which(mer_hrze$global_padj<=0.05),],
#              aes(y=-logFC,x=-1*Log_2_Fold.change._Active_Control,
#                  color=renorm),
#              size=2,alpha=.9)+
#   scale_color_manual(name="",values = mycol) +
#   theme_base()+
#   xlab("Log2 Fold Change Control/Active (Berry)")+
#   ylab("Log2 Fold Change Post [Any] /Pre (This Study)")
# 
# hrze_berry_fcp
# pdf(paste(fig_folder,"rnaseq_fc_berry_fc_hrze_eba.pdf",sep = "/"),height = 5, width = 6)
# print(hrze_berry_fcp)
# dev.off()

berry_hrze_sig <- berry_hrze[berry_hrze$global_padj < 0.05,] 
berry_hrze_sig
write.csv(berry_hrze_sig, paste(tab_folder,paste("EBA_Berry_HRZE_differential_Day_56_Day_14",".csv",sep=""),sep = "/"))

# ibd hrze eba
fc_dt <- ibd_hrze
fc_dt$renorm <- "Exacerbate"
fc_dt$renorm[fc_dt$logFC < 0 & fc_dt$Case.Over.Control.Fold.Change < 0]<-"Renormalize"
fc_dt$renorm[fc_dt$logFC > 0 & fc_dt$Case.Over.Control.Fold.Change > 0]<-"Renormalize"
fc_dt$otu <- rownames(fc_dt)
mer_hrze <- fc_dt

# Adding excel mapped names for superpathway
data_pvalue_less <- mer_hrze[which(mer_hrze$global_padj<=0.05),]
hrze_ibdp<-ggplot() + 
  geom_point(data=mer_hrze[which(mer_hrze$global_padj>0.05),],
             aes(y=-log10(global_padj),x=-logFC),
             size=2,color="black",alpha=0.1)+
  geom_point(data= data_pvalue_less,
             aes(y=-log10(global_padj),x=-logFC,color=renorm),
             size=2,alpha=.9)+
  scale_color_manual(name="",values = mycol) +
  theme_base()+
  xlab("Log2 Fold Change (Day 56/Day 14)")+
  ylab("-Log10 FDR")+
  xlim(-3,3)+
  ylim(0,8)
hrze_ibdp

pdf(paste(fig_folder,"rnaseq_volcano_plot_ibd_hrze_eba_Day_56_Day_14.pdf",sep = "/"),height = 5, width = 6)
print(hrze_ibdp)
dev.off()

# plot FC this study vs FC Berry
# hrze_ibd_fcp<-ggplot() + 
#   geom_point(data=mer_hrze[which(mer_hrze$global_padj>0.05),],
#              aes(y=-logFC,x=-1*Case.Over.Control.Fold.Change),
#              size=2,color="black",alpha=0.1)+
#   geom_point(data=mer_hrze[which(mer_hrze$global_padj<=0.05),],
#              aes(y=-logFC,x=-1*Case.Over.Control.Fold.Change,
#                  color=renorm),
#              size=2,alpha=.9)+
#   scale_color_manual(name="",values = mycol) +
#   theme_base()+
#   xlab("Log2 Fold Change Case/Control (IBD)")+
#   ylab("Log2 Fold Change Post/Pre (This Study)")
# 
# hrze_ibd_fcp
# pdf(paste(fig_folder,"rnaseq_fc_ibd_fc_hrze_eba.pdf",sep = "/"),height = 5, width = 6)
# print(hrze_ibd_fcp)
# dev.off()

ibd_hrze_sig <- ibd_hrze[ibd_hrze$global_padj < 0.05,] 
nrow(ibd_hrze_sig)
write.csv(ibd_hrze_sig, paste(tab_folder,paste("EBA_IBD_HRZE_differential_Day_56_Day_14",".csv",sep=""),sep = "/"))



