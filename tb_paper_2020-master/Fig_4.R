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
mainDir <- "../Figs"
subDir <- "Limma_rna_drug"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
fig_folder <- paste(mainDir,subDir,sep="/")

mainDir <- "../Tables"
subDir <- "Limma_rna_drug"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
tab_folder <- paste(mainDir,subDir,sep="/")


####### load phyloseq objected created with data_setup.R ######

phy_pair <- readRDS("../data/genes/phy_rna_ntz_hrze_pair.rds")



###### Differential Analysis with LIMMA ########
phy_trt <-  phy_pair
sm_dt <- data.frame(sample_data(phy_trt))
sm_dt$Treatment <-  "HRZE"
sm_dt$Treatment[grep("NTZ",sm_dt$type)] <- "NTZ" 
sm_dt$Time<- "Post" 
sm_dt$Time[grep("pre",sm_dt$TB_status)] <- "Pre" 
sm_dt$SampleID <- rownames(sm_dt)

sm_dt <- sm_dt %>%
  arrange(Treatment,Patient.ID,Time) %>%
  as.data.frame()

sm_dt$ind <- as.numeric(factor(as.character(sm_dt$Patient.ID),
                               levels = unique(as.character(sm_dt$Patient.ID))))
rownames(sm_dt) <- sm_dt$SampleID

sample_data(phy_trt) <- sample_data(sm_dt)

table(get_variable(phy_trt,"Time"))
table(get_variable(phy_trt,"Treatment"))
sample_data(phy_trt)$Patient.ID



# Filter counts that the total counts across all samples is >= cnt_filt
cnt_filt <- 50

phy_sel <- phy_trt
phy_sel <-  prune_taxa(taxa_sums(phy_sel)>0,phy_sel)
#phy_sel <-  subset_samples(phy_sel, type %in% c("Pre NTZ", "Pre HRZE"))

# Differential analysis via Limma
cnt_dt <- as.matrix(data.frame(otu_table(phy_sel)))
sm_dt <-  data.frame(sample_data(phy_sel))
sm_dt$Time <- factor(sm_dt$Time,levels = c("Pre","Post"))
sm_dt$ind <- factor(sm_dt$ind)
sm_dt$Treatment <- factor(sm_dt$Treatment,c("HRZE","NTZ"))

sm_dt$drug <- factor(sm_dt$drug,c("pretreatment","HRZE","NTZ"))
sm_dt$batchID

cnt_dt <-  cnt_dt[,match(rownames(sm_dt),colnames(cnt_dt))]


dim(cnt_dt)
# filter genes by number of counts
#table(rowSums(cpm(cnt_dt)>0.1)==8)
sum(rowSums(cnt_dt) >= cnt_filt)
isexpr = rowSums(cnt_dt) >= cnt_filt
table(rowSums(cnt_dt) >=cnt_filt)
# Standard usage of limma/voom
gExpr = DGEList( cnt_dt[isexpr,] )
gExpr = calcNormFactors( gExpr )
dim(gExpr$counts)

sm_dt$batchID
sm_dt$Time
sm_dt$ind
sm_dt$sex
sm_dt$Treatment
sm_dt$Av_TTP
sm_dt$Age
sm_dt$drug

#Visualize
library(RColorBrewer)

# Limma 
# apply duplicateCorrelation is two rounds
# Here NTZ and HRZE samples are from two different batches
# drug is pretreatment, HRZE  (Post)  and NTZ (Post)
design = model.matrix( ~  sex + batchID + drug  , sm_dt)
# apply duplicateCorrelation is two rounds
vobj_tmp = voom(gExpr, design, plot=T)
dupcor <- duplicateCorrelation(vobj_tmp,design,block=sm_dt$ind)

# run voom considering the duplicateCorrelation results
# in order to compute more accurate precision weights
# Otherwise, use the results from the first voom run
pdf(paste0(tab_folder,"/Mean_Dispersion_plot_hrze_ntz.pdf"),width = 5,height = 5)
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

#final_fit <- fitDupCor
p_cut_off <-  0.05 

#Create a contrast matrix for the required comparisons
cont.mat <- cbind(Intercept = c(1,0,0,0,0),
                  HRZE_Post_Vs_Pre = c(0,0,0,1,0) ,
                  NTZ_Post_Vs_Pre = c(0,0,0,0,1) ,
                  Sex_M_Vs_F = c(0,1,0,0,0),
                  Batch_NTZ_HRZE = c(0,0,1,0,0))

final_fit <-  contrasts.fit(fitDupCor, cont.mat)
final_fit <- eBayes(final_fit)

sum_hrze_ntz <- summary( decideTests(final_fit,p.value = p_cut_off,
                                     method = "global",adjust.method = "BH"))
sum_dt <-  as.data.frame.matrix(sum_hrze_ntz) 
write.csv(sum_dt ,paste(tab_folder,paste("Summary_rna_HRZE_NTZ.csv",sep=""),sep = "/"))

# Global pval adjustment
nom_pval <-  final_fit$p.value
rnames <- rownames(nom_pval)

vec_pval <- c(nom_pval)
global_pval <- p.adjust(nom_pval,method = "BH")
global_mat <- data.frame(matrix(global_pval,ncol = ncol(nom_pval)))
names(global_mat) <- colnames(nom_pval)
rownames(global_mat) <-  rnames




# Saving results - All genes
# Differentially affected Genes by HRZE compared to HRZE baseline (HRZE Pre) 
res_hrze_post <- topTable( final_fit, coef=c("HRZE_Post_Vs_Pre"),n=Inf)
p_cut_off <-  0.05
res_hrze_post <- data.frame(res_hrze_post)

# Match the order with global_mat
res_hrze_post <- res_hrze_post[match(rownames(global_mat),rownames(res_hrze_post)),]


res_hrze_post$global_padj <- global_mat$HRZE_Post_Vs_Pre
res_hrze_post <- res_hrze_post[!is.na(res_hrze_post$adj.P.Val),]
res_hrze_post <- res_hrze_post[order(res_hrze_post$global_padj,decreasing = F),]

#res_hrze_post_sig <- res_hrze_post[res_hrze_post$adj.P.Val <= p_cut_off,]
write.csv(res_hrze_post,paste(tab_folder,paste("HRZE_Post_Vs_Pre_",cnt_filt,".csv",sep=""),sep = "/"))

# Differentially affected Genes by NTZ compared to NTZ baseline
res_ntz_post <- topTable( final_fit, coef=c("NTZ_Post_Vs_Pre"), n=Inf )
p_cut_off <-  0.05
res_ntz_post <- data.frame(res_ntz_post)
# Match the order with global_mat
res_ntz_post <- res_ntz_post[match(rownames(global_mat),rownames(res_ntz_post)),]

# Adjusted global p val
res_ntz_post$global_padj <- global_mat$NTZ_Post_Vs_Pre
res_ntz_post <- res_ntz_post[!is.na(res_ntz_post$adj.P.Val),]
res_ntz_post <- res_ntz_post[order(res_ntz_post$global_padj,decreasing = F),]

#res_ntz_post_sig <- res_ntz_post[res_ntz_post$adj.P.Val <= p_cut_off,]
write.csv(res_ntz_post,paste(tab_folder,paste("NTZ_Post_Vs_Pre_",cnt_filt,".csv",sep=""),sep = "/"))


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





# Map Berry gene set in HRZE group

berry_data <- read.csv("../Berry/HRZE_NTZ_Berry_List_April_2020.csv")
ibd_data <- read.csv("../IBD/IBD_sig_genes.csv")


#  Berry subset
# idea is to map the TB related genes in HRZE and NTZ data

# Mapping to HRZE Set
res_hrze <- res_hrze_post
res_hrze$Gene <-  rownames(res_hrze_post)
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

write.csv(berry_hrze_sig, paste(tab_folder,paste("Berry_HRZE_differential",".csv",sep=""),sep = "/"))

# Mapping to HRZE Set
res_ntz <- res_ntz_post
res_ntz$Gene <-  rownames(res_ntz)
res_ntz$Ensembl_ID <-  gsub("\\..*","",res_ntz$Gene)

# Merge the NTZ results with berry signature
berry_ntz <-  merge(berry_data,res_ntz,by.x=  "Ensembl_ID", by.y = "Ensembl_ID")

# Active Vs Control
# Adjusting direction with our data
berry_ntz$logFC <-  -1*berry_ntz$logFC

berry_ntz$Comb_ID <- as.character(berry_ntz$Comb_ID)
# Number of genes that are significantly different in NTZ Pre and Post 
# And also TB related 
berry_ntz_sig <- berry_ntz[berry_ntz$global_padj < 0.05,] 
nrow(berry_ntz_sig)

write.csv(berry_ntz_sig, paste(tab_folder,paste("Berry_NTZ_differential",".csv",sep=""),sep = "/"))


# Comparison of LFC
# Both decreasing
down_both_dt <- berry_ntz[berry_ntz$logFC < 0 & 
                            berry_ntz$Log_2_Fold.change._Active_Control < 0,  ]
# Both increasing
up_both_dt <- berry_ntz[berry_ntz$logFC > 0 & 
                          berry_ntz$Log_2_Fold.change._Active_Control > 0,  ]

mismatch_dt <- berry_ntz[!berry_ntz$Comb_ID %in%
                           c(down_both_dt$Comb_ID, up_both_dt$Comb_ID),]

# IBD Set
# Map IBD set to HRZE group
ibd_data <- read.csv("../IBD/IBD_sig_genes.csv")
ibd_data$Gene.Symbol <- trimws(as.character(ibd_data$Gene.Symbol))
#  IBD subset
# idea is to map the IBD  genes in HRZE and NTZ data

# Mapping to HRZE Set
res_hrze <- res_hrze_post
res_hrze$Gene <-  rownames(res_hrze_post)
res_hrze$symbol <-  gsub(".*_","",res_hrze$Gene)

# Merge the HRZE results with berry signature
ibd_hrze <-  merge(ibd_data,res_hrze,by.x=  "Gene.Symbol", by.y = "symbol")

# Number of genes that are significantly different in HRZE Pre and Post 
# And also IBD related 
ibd_hrze_sig <- ibd_hrze[ibd_hrze$global_padj < 0.05,] 
nrow(ibd_hrze_sig)

write.csv(ibd_hrze_sig, paste(tab_folder,paste("IBD_HRZE_differential",".csv",sep=""),sep = "/"))

# Mapping to NTZ Set
res_ntz <- res_ntz_post
res_ntz$Gene <-  rownames(res_ntz)
res_ntz$symbol <-  gsub(".*_","",res_ntz$Gene)

# Merge the HRZE results with berry signature
ibd_ntz <-  merge(ibd_data,res_ntz,by.x=  "Gene.Symbol", by.y = "symbol")

# Number of genes that are significantly different in HRZE Pre and Post 
# And also IBD related 
ibd_ntz_sig <- ibd_ntz[ibd_ntz$global_padj < 0.05,] 
nrow(ibd_ntz_sig)
write.csv(ibd_ntz_sig, paste(tab_folder,paste("IBD_NTZ_differential",".csv",sep=""),sep = "/"))

# plot the genes against previous LFC
#berry_hrze_sig_expand <- cbind(berry_hrze_sig,berry_data[match(berry_hrze_sig$ Ensembl_ID,berry_data$ Ensembl_ID),])
#berry_ntz_sig_expand <- cbind(berry_ntz_sig,berry_data[match(berry_hrze_sig$ Ensembl_ID,berry_data$ Ensembl_ID),])
#ibd_hrze_sig_expand <- cbind(ibd_hrze_sig,ibd_data[match(ibd_hrze_sig$Entrez.Id,ibd_data$Entrez.Id),])
#ibd_ntz_sig_expand <- cbind(ibd_ntz_sig,ibd_data[match(ibd_ntz_sig$Entrez.Id,ibd_data$Entrez.Id),])

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
  xlab("Log2 Fold Change (Post/Pre)")+
  ylab("-Log10 P-Value")
hrze_berryp

pdf(paste(fig_folder,"rnaseq_volcano_plot_berry_hrze.pdf",sep = "/"),height = 5, width = 6)
print(hrze_berryp)
dev.off()

# plot FC this study vs FC Berry
hrze_berry_fcp<-ggplot() + 
  geom_point(data=mer_hrze[which(mer_hrze$global_padj>0.05),],
             aes(y=-logFC,x=-1*Log_2_Fold.change._Active_Control),
             size=2,color="black",alpha=0.1)+
  geom_point(data=mer_hrze[which(mer_hrze$global_padj<=0.05),],
             aes(y=-logFC,x=-1*Log_2_Fold.change._Active_Control,
                 color=renorm),
             size=2,alpha=.9)+
  geom_abline()+
  scale_color_manual(name="",values = mycol) +
  theme_base()+
  xlim(-7,2)+
  ylim(-7,2)+
  xlab("Log2 Fold Change Control/Active (Berry)")+
  ylab("Log2 Fold Change Post/Pre (This Study)")

#hrze_berry_fcp
pdf(paste(fig_folder,"rnaseq_fc_berry_fc_hrze.pdf",sep = "/"),height = 5, width = 6)
print(hrze_berry_fcp)
dev.off()

# ntz berry
fc_dt <- berry_ntz
fc_dt$renorm <- "Exacerbate"
fc_dt$renorm[fc_dt$logFC < 0 & fc_dt$Log_2_Fold.change._Active_Control < 0]<-"Renormalize"
fc_dt$renorm[fc_dt$logFC > 0 & fc_dt$Log_2_Fold.change._Active_Control > 0]<-"Renormalize"

fc_dt$otu <- rownames(fc_dt)
mer_ntz <- fc_dt

data_pvalue_less <- mer_ntz[which(mer_ntz$global_padj<=0.05),]
mer_ntz_berryp<-ggplot() + 
  geom_point(data=mer_ntz[which(mer_ntz$global_padj>0.05),],
             aes(y=-log10(global_padj),x=-logFC),
             size=2,color="black",alpha=0.1)+
  geom_point(data= data_pvalue_less,
             aes(y=-log10(global_padj),x=-logFC,color=renorm),
             size=2,alpha=.9)+
  scale_color_manual(name="",values = mycol) +
  theme_base()+
  xlab("Log2 Fold Change (Post/Pre)")+
  ylab("-Log10 P-Value")
mer_ntz_berryp

pdf(paste(fig_folder,"rnaseq_volcano_plot_berry_ntz.pdf",sep = "/"),height = 5, width = 6)
print(mer_ntz_berryp)
dev.off()

# plot FC this study vs FC Berry
ntz_berry_fcp<-ggplot() + 
  geom_point(data=mer_ntz[which(mer_ntz$global_padj>0.05),],
             aes(y=-logFC,x=-1*Log_2_Fold.change._Active_Control),
             size=2,color="black",alpha=0.1)+
  geom_point(data=mer_ntz[which(mer_ntz$global_padj<=0.05),],
             aes(y=-logFC,x=-1*Log_2_Fold.change._Active_Control,
                 color=renorm),
             size=2,alpha=.9)+
  scale_color_manual(name="",values = mycol) +
  theme_base()+
  xlim(-7,2)+
  ylim(-7,2)+
  xlab("Log2 Fold Change Control/Active (Berry)")+
  ylab("Log2 Fold Change Post/Pre (This Study)")

ntz_berry_fcp
pdf(paste(fig_folder,"rnaseq_fc_berry_fc_ntz.pdf",sep = "/"),height = 5, width = 6)
print(hrze_berry_fcp)
dev.off()

# ibd hrze
mycol <-c("red","#005073")
names(mycol) <-c("Exacerbate","Renormalize")
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
  xlab("Log2 Fold Change (Post/Pre)")+
  ylab("-Log10 P-Value")
hrze_ibdp

pdf(paste(fig_folder,"rnaseq_volcano_plot_ibd_hrze.pdf",sep = "/"),height = 5, width = 6)
print(hrze_ibdp)
dev.off()

# plot FC this study vs FC Berry
hrze_ibd_fcp<-ggplot() + 
  geom_point(data=mer_hrze[which(mer_hrze$global_padj>0.05),],
             aes(y=-logFC,x=-1*Case.Over.Control.Fold.Change),
             size=2,color="black",alpha=0.1)+
  geom_point(data=mer_hrze[which(mer_hrze$global_padj<=0.05),],
             aes(y=-logFC,x=-1*Case.Over.Control.Fold.Change,
                 color=renorm),
             size=2,alpha=.9)+
  scale_color_manual(name="",values = mycol) +
  theme_base()+
  geom_abline()+
  xlim(-4,2)+
  ylim(-4,2)+
  xlab("Log2 Fold Change Control/Case (IBD)")+
  ylab("Log2 Fold Change Post/Pre (This Study)")

hrze_ibd_fcp
pdf(paste(fig_folder,"rnaseq_fc_ibd_fc_hrze.pdf",sep = "/"),height = 5, width = 6)
print(hrze_ibd_fcp)
dev.off()

# ibd ntz
mycol <-c("red","#005073")
names(mycol) <-c("Exacerbate","Renormalize")
fc_dt <- ibd_ntz
fc_dt$renorm <- "Exacerbate"
fc_dt$renorm[fc_dt$logFC < 0 & fc_dt$Case.Over.Control.Fold.Change < 0]<-"Renormalize"
fc_dt$renorm[fc_dt$logFC > 0 & fc_dt$Case.Over.Control.Fold.Change > 0]<-"Renormalize"
fc_dt$otu <- rownames(fc_dt)
mer_ntz <- fc_dt

# Adding excel mapped names for superpathway
data_pvalue_less <- mer_ntz[which(mer_ntz$global_padj<=0.05),]
ntz_ibdp<-ggplot() + 
  geom_point(data=mer_ntz[which(mer_ntz$global_padj>0.05),],
             aes(y=-log10(global_padj),x=-logFC),
             size=2,color="black",alpha=0.1)+
  geom_point(data= data_pvalue_less,
             aes(y=-log10(global_padj),x=-logFC,color=renorm),
             size=2,alpha=.9)+
  scale_color_manual(name="",values = mycol) +
  theme_base()+
  xlab("Log2 Fold Change (Post/Pre)")+
  ylab("-Log10 P-Value")
ntz_ibdp

pdf(paste(fig_folder,"rnaseq_volcano_plot_ibd_ntz.pdf",sep = "/"),height = 5, width = 6)
print(ntz_ibdp)
dev.off()

# plot FC this study vs FC Berry
ntz_ibd_fcp<-ggplot() + 
  geom_point(data=mer_ntz[which(mer_ntz$global_padj>0.05),],
             aes(y=-logFC,x=-1*Case.Over.Control.Fold.Change),
             size=2,color="black",alpha=0.1)+
  geom_point(data=mer_ntz[which(mer_ntz$global_padj<=0.05),],
             aes(y=-logFC,x=-1*Case.Over.Control.Fold.Change,
                 color=renorm),
             size=2,alpha=.9)+
  scale_color_manual(name="",values = mycol) +
  theme_base()+
  geom_abline()+
  xlim(-4,1)+
  ylim(-4,1)+
  xlab("Log2 Fold Change Case/Control (IBD)")+
  ylab("Log2 Fold Change Post/Pre (This Study)")

ntz_ibd_fcp
pdf(paste(fig_folder,"rnaseq_fc_ibd_fc_ntz.pdf",sep = "/"),height = 5, width = 6)
print(ntz_ibd_fcp)
dev.off()
