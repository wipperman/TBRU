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
library(hues)
library(ggthemes)


# Create a directory to save figures and tables
mainDir <- "../Figs"
subDir <- "Limma_mic_drug"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
fig_folder <- paste(mainDir,subDir,sep="/")

mainDir <- "../Tables/"
subDir <- "Limma_mic_drug"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
tab_folder <- paste(mainDir,subDir,sep="/")


####### load phyloseq objected created with data_setup.R ######

phy <- readRDS("../data/mic/phy_mic_fil_paired.rds")

perc_samples <-  0.1
phy_fil = filter_taxa(phy, function(x) sum(x > 0) > (perc_samples*length(x)), TRUE)
phy_fil <- prune_taxa(taxa_sums(phy_fil)> 0, phy_fil)
phy_pair <- phy_fil

sm_dt <- data.frame(sample_data(phy_pair))
sm_dt$sample_name <- rownames(sm_dt)



###### Differential Analysis with LIMMA ########
phy_trt <-  phy_pair
sm_dt <- data.frame(sample_data(phy_trt))

sm_dt$Treatment <-  "HRZE"
sm_dt$Treatment[grep("NTZ",sm_dt$type)] <- "NTZ" 
sm_dt$Time<- "Post" 
sm_dt$Time[grep("pre",sm_dt$TB_status)] <- "Pre" 
#sm_dt$SampleID <- rownames(sm_dt)

sm_dt <- sm_dt %>%
  arrange(Treatment,Patient.ID,Time) %>%
  as.data.frame()

sm_dt$ind <- as.numeric(factor(as.character(sm_dt$Patient.ID),
                               levels = unique(as.character(sm_dt$Patient.ID))))
rownames(sm_dt) <- sm_dt$sample_name

sample_data(phy_trt) <- sample_data(sm_dt)

table(get_variable(phy_trt,"Time"))
table(get_variable(phy_trt,"Treatment"))
sample_data(phy_trt)$Patient.ID




# Filter counts that the total counts across all samples is >= cnt_filt
cnt_filt <- 50
phy_sel <- phy_trt
phy_sel <-  prune_taxa(taxa_sums(phy_sel)>0,phy_sel)

# Differential analysis via Limma
cnt_dt <- as.matrix(data.frame(otu_table(phy_sel)))
sm_dt <-  data.frame(sample_data(phy_sel))
sm_dt$Time <- factor(sm_dt$Time,levels = c("Pre","Post"))
sm_dt$ind <- factor(sm_dt$ind)
sm_dt$Treatment <- factor(sm_dt$Treatment,c("HRZE","NTZ"))
sm_dt$drug <- factor(sm_dt$drug,c("pretreatment","HRZE","NTZ"))
sm_dt$pool_batch <- factor(sm_dt$pool_batch)


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

#sm_dt$batchID <- gsub("Pre |Post ","",sm_dt$type)
#sm_dt$batchID <- factor(sm_dt$batchID)
sm_dt$Time
sm_dt$ind
sm_dt$sex
sm_dt$Treatment
sm_dt$drug
sm_dt$Age
sm_dt$sex
#Visualize
# library(RColorBrewer)
# lcpm <- cpm(gExpr$counts, log=TRUE)
# group <- paste0(sm_dt$Time,"_",sm_dt$ind,"_",sm_dt$sex)
# group <- paste0(sm_dt$Time,"_",sm_dt$ind,"_",sm_dt$sex,"_",sm_dt$Age,"_",sm_dt$Av_TTP)
# group <- paste0(sm_dt$Time,"_",sm_dt$ind,"_",sm_dt$sex,"_",sm_dt$Age)
# plotMDS(lcpm, labels=group)
# 

#Visualize

# Limma 
# apply duplicateCorrelation is two rounds
design = model.matrix( ~ sex + pool_batch + drug  , sm_dt)
# apply duplicateCorrelation is two rounds
vobj_tmp = voom(gExpr, design, plot=FALSE)
dupcor <- duplicateCorrelation(vobj_tmp,design,block=sm_dt$ind)

# run voom considering the duplicateCorrelation results
# in order to compute more accurate precision weights
# Otherwise, use the results from the first voom run
pdf(paste0(fig_folder,"/Mean_Dispersion_plot_hrze_ntz.pdf"),width = 5,height = 5)
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

p_cut_off <-  0.05 


final_fit <- fitDupCor


sum_hrze_ntz <-  summary(decideTests(fitDupCor,p.value = p_cut_off,method = "global",
                                     adjust.method = "BH"))

sum_dt <-  as.data.frame.matrix(sum_hrze_ntz) 
write.csv(sum_dt ,paste(tab_folder,paste("Summary_Mic_NTZ_HRZE.csv",sep=""),sep = "/"))


# Global pval adjustment
nom_pval <-  final_fit$p.value
rnames <- rownames(nom_pval)

vec_pval <- c(nom_pval)
global_pval <- p.adjust(nom_pval,method = "BH")
global_mat <- data.frame(matrix(global_pval,ncol = ncol(nom_pval)))
names(global_mat) <- colnames(nom_pval)
rownames(global_mat) <-  rnames


# Saving results - All Species
# HRZE Significant
# Differentially affected ASVs by HRZE compared to HRZE baseline (HRZE Pre) 
res_hrze_post <- topTable( final_fit, coef=c("drugHRZE"),n=Inf)
#res_hrze_post <- topTable( fitDupCor, coef=c("batchIDNTZ"),n=Inf)

p_cut_off <-  0.05
res_hrze_post <- data.frame(res_hrze_post)

# Match the order with global_mat
res_hrze_post <- res_hrze_post[match(rownames(global_mat),rownames(res_hrze_post)),]

res_hrze_post$global_padj <- global_mat$drugHRZE

res_hrze_post <- res_hrze_post[!is.na(res_hrze_post$adj.P.Val),]
res_hrze_post <- res_hrze_post[order(res_hrze_post$global_padj,decreasing = F),]

res_hrze_post_sig <- res_hrze_post[res_hrze_post$global_padj <= p_cut_off,]
write.csv(res_hrze_post,paste(tab_folder,paste("HRZE_Post_Vs_Pre_",cnt_filt,".csv",sep=""),sep = "/"))

# NTZ Significant

# Differentially affected ASVs by NTZ compared to NTZ baseline
res_ntz_post <- topTable( final_fit, coef=c("drugNTZ"), n=Inf )
p_cut_off <-  0.05
res_ntz_post <- data.frame(res_ntz_post)
# Match the order with global_mat
res_ntz_post <- res_ntz_post[match(rownames(global_mat),rownames(res_ntz_post)),]

# Adjusted global p val
res_ntz_post$global_padj <- global_mat$drugNTZ
res_ntz_post <- res_ntz_post[!is.na(res_ntz_post$adj.P.Val),]
res_ntz_post <- res_ntz_post[order(res_ntz_post$global_padj,decreasing = F),]

res_ntz_post_sig <- res_ntz_post[res_ntz_post$global_padj <= p_cut_off,]
write.csv(res_ntz_post,paste(tab_folder,paste("NTZ_Post_Vs_Pre_",cnt_filt,".csv",sep=""),sep = "/"))



####### Differential analysis plotting #######################
mat <- cpm(cnt_dt, log=TRUE)
phy_ht <- phy_pair
otu_table(phy_ht) <-  otu_table(mat, taxa_are_rows = T)
#common_otus <- intersect(rownames(res_hrze_post),rownames(res_ntz_post))
all_sp <- union(rownames(res_hrze_post),rownames(res_ntz_post))
phy_melt <-  get.otu.melt(phy_ht)
tax_dt <-  unique(phy_melt[,c("otu","Order","Phylum")])
topN <- length(unique(tax_dt$Order))
match_tax_dt <- unique(tax_dt [,c("Order","Phylum")])
match_tax_dt <- match_tax_dt[order(match_tax_dt$Phylum,match_tax_dt$Order),]

phylum_col <- c("Firmicutes"= "#9C854E",
                "Bacteroidetes" = "#51AB9B",
                "Actinobacteria" = "#A77097",
                "Proteobacteria" = "red",
                "Verrucomicrobia" = "#33445e",
                "Tenericutes" = "#e2eab7",
                "Fusobacteria"= "#653131",
                "Cyanobacteria" = "#808080",
                "Euryarchaeota" = "#8c8c8c",
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

# set.seed(1059)
# mycol <- iwanthue(topN, 0, 360, 36, 180, 13, 73, plot=TRUE)
# 
# # Create palettes
# names(mycol) <- unique(tax_dt$Order)
# 
# # Set Colors based on previous palette/ying palettee
# col_palette <-get.yt.palette2(tax_table(phy_ht)[,3:9])

# Volcano plots for HRZE
# HRZE
fc_dt <-  res_hrze_post
fc_dt$otu <- rownames(fc_dt)
mer_hrze <- merge(tax_dt,fc_dt,by = "otu")
mer_hrze$Type <- "HRZE"
mer_hrze

p_cut_off <- 0.05

data_pvalue_less <- mer_hrze[which(mer_hrze$global_padj< p_cut_off),]
hrze_vp<-ggplot() + 
  geom_point(data=mer_hrze[which(mer_hrze$global_padj>=p_cut_off),],
             aes(y=-log10(global_padj),x=logFC),
             size=2,color="black",alpha=0.1)+
  geom_point(data= data_pvalue_less,
             aes(y=-log10(global_padj),x=logFC,color=Order),
             size=2,alpha=.9)+
  scale_color_manual(name="Order",values = mycol) +
  geom_hline(yintercept = -log10(p_cut_off),size = .1,linetype = "dashed")+
  geom_vline(xintercept = -1,size = .1,linetype = "dashed")+
  geom_vline(xintercept = 1,size = .1,linetype = "dashed")+
  theme_base()+
  #scale_fill_manual(name="",values = mycol[[1]],breaks=mycol[[2]]) +
  xlab("Log2 Fold Change")+
  ylab("-Log10 P-Value")+
  xlim(-10,10)+
  ylim(0,18)

print(hrze_vp)

# NTZ
fc_dt <-  res_ntz_post
fc_dt$otu <- rownames(fc_dt)
mer_ntz <- merge(tax_dt,fc_dt,by = "otu")
mer_ntz$Type <- "NTZ"

# Adding excel mapped names for superpathway
data_pvalue_less <- mer_ntz[which(mer_ntz$global_padj<p_cut_off),]
ntz_vp<-ggplot() + 
  geom_point(data=mer_hrze[which(mer_hrze$global_padj>=p_cut_off),],
             aes(y=-log10(global_padj),x=logFC),
             size=2,color="black",alpha=0.1)+
  geom_point(data= data_pvalue_less,
             aes(y=-log10(global_padj),x=logFC,color=Order),
             size=2,alpha=.9)+
  scale_color_manual(name="Order",values = mycol) +
  geom_hline(yintercept = -log10(p_cut_off),size = .1,linetype = "dashed")+
  geom_vline(xintercept = -1,size = .1,linetype = "dashed")+
  geom_vline(xintercept = 1,size = .1,linetype = "dashed")+
  theme_base()+
  #scale_fill_manual(name="",values = mycol[[1]],breaks=mycol[[2]]) +
  xlab("Log2 Fold Change")+
  ylab("-Log10 P-Value")+
  xlim(-10,10)+
  ylim(0,18)
print(ntz_vp)



# Final figure
comb_dt <-  rbind(mer_hrze,mer_ntz)
comb_dt$Type[which(comb_dt$Type=="Int_NTZ_HRZE")]<-"Contrast (NTZ/HRZE)"
unique(comb_dt$Type)
comb_dt$Type <- factor(comb_dt$Type,c("HRZE", "NTZ","Contrast (NTZ/HRZE)"))

#levels(comb_dt$Type) <-  c("log2(HRZE Post/ HRZE Pre)", "log2(NTZ Post/ NTZ Pre)",
#                          "log2(log2(NTZ Post/NTZ Pre)/log2(HRZE Post/HRZE Pre))")


# Adding excel mapped names for superpathway
data_pvalue_less <- comb_dt[which(comb_dt$global_padj<p_cut_off),]
comb_vp<-ggplot() + 
  geom_point(data=comb_dt[which(comb_dt$global_padj>=p_cut_off),],
             aes(y=-log10(global_padj),x=logFC),
             size=2,color="black",alpha=0.2)+
  geom_point(data=  comb_dt[which(comb_dt$global_padj<p_cut_off),],
             aes(y=-log10(global_padj),x=logFC,color=Order),
             size=2,alpha=.8)+
  scale_color_manual(name="Order",values = mycol) +
  facet_wrap(~Type,nrow = 1,scales = "free")+
  geom_hline(yintercept = -log10(p_cut_off),size = .1,linetype = "dashed")+
  geom_vline(xintercept = -1,size = .1,linetype = "dashed")+
  geom_vline(xintercept = 1,size = .1,linetype = "dashed")+
  theme_base()+
  guides(colour = guide_legend(override.aes = list(shape = 15,size = 7)))+
  #scale_fill_manual(name="",values = mycol[[1]],breaks=mycol[[2]]) +
  xlab("Log2 Fold Change")+
  ylab("-Log10 P-Value")+
  xlim(-10,10)+
  ylim(0,18)
print(comb_vp)
pdf(paste(fig_folder,"volcano_plot_combined.pdf",sep = "/"),height = 5, width = 12)
print(comb_vp)
dev.off()




# Draw a heatmap for HRZE only 

mat <- cpm(cnt_dt, log=TRUE)
# VST phyloseq for Heatmap 
phy_ht <- phy_pair
otu_table(phy_ht) <-  otu_table(mat, taxa_are_rows = T)
# common_otus <- intersect(rownames(res_hrze_post))
# all_sp <- union(union(rownames(res_hrze_post),rownames(res_ntz_post)),rownames(res_int))
phy_ht <-  subset_samples(phy_ht,drug %in% c("pretreatment","HRZE"))


phy_sel_ht <- prune_taxa(rownames(res_hrze_post[res_hrze_post$global_padj < 0.05,]),phy_ht)
mat <-data.frame(otu_table(phy_sel_ht))
df <- as.data.frame(sample_data(phy_sel_ht)[,c("drug","TB_status","type","pool_batch")]) #"days_ON_HRZE","TB_status"
#df$type <- NULL
df$drug <- factor(df$drug, levels = c("pretreatment","NTZ","HRZE"))
color_pal <- list(drug = c("pretreatment" = "#537c4a","HRZE"="#984ea3","NTZ" = "#e46983"))
match_tax_dt <- data.frame(tax_table(phy_sel_ht))
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

match_tax_dt <- match_tax_dt[match(rownames(mat),rownames(match_tax_dt)),]




ha1 = HeatmapAnnotation(Order = match_tax_dt$Order,
                        col = list(Order = order_col), which = "row")
ha_column = HeatmapAnnotation(Drug =  df$drug,
                              batch = df$pool_batch,
                              #Type = df$type,
                              col=list(Drug = c("pretreatment" = "#537c4a","HRZE"="#984ea3","NTZ" = "#e46983")))
split_rows <-  match_tax_dt$Phylum
split_rows <- factor(split_rows, levels= unique(split_rows))
#colfunc <- colorRampPalette(c("purple", "yellow"))

library(RColorBrewer)
colfunc <- colorRampPalette(rev(brewer.pal(n = 7, name =
                                             "RdYlBu")))
#teddy_cols <- c('white','#fecc5c','#fd8d3c','#f03b20','#bd0026')

ht_cols <- colfunc(8)
# Remove underscores 
rownames(mat) <- gsub("_"," ",rownames(mat))
ht1 = Heatmap(as.matrix(mat), name = "log", column_title = NA, 
              top_annotation = ha_column,
              clustering_distance_rows = "euclidean",
              column_split = factor(df$type),
              split = split_rows, row_gap = unit(2, "mm"),border = TRUE,
              row_title_gp = gpar(fontsize = 10,fontface = "bold"),
              row_title_rot = 0,
              left_annotation = ha1,
              #  right_annotation = ha2,
              row_names_side = "left", km=1, color_space = "LAB",
              #row_dend_side="right",
              col= ht_cols,
              show_parent_dend_line = F,
              #heatmap_legend_param = list(),
              #clustering_method_columns = "ward.D",
              width=2, 
              #show_column_names= T,
              row_names_gp = gpar(fontsize = 9),
              #cluster_columns = T,cluster_rows = T,
              row_names_max_width = max_text_width(rownames(mat), gp = gpar(fontsize = 12)),
              column_names_max_height = max_text_width(colnames(mat), gp = gpar(fontsize = 12)),
              show_column_names = "F",
              column_names_side = "top",
              na_col="white")

pdf(paste(fig_folder,"Heatmap_hrze.pdf",sep = "/"),height = 10, width = 14)
draw(ht1)
dev.off()



# Draw a heatmap for NTZ only 

mat <- cpm(cnt_dt, log=TRUE)
# VST phyloseq for Heatmap 
phy_ht <- phy_pair
otu_table(phy_ht) <-  otu_table(mat, taxa_are_rows = T)
# common_otus <- intersect(rownames(res_hrze_post))
# all_sp <- union(union(rownames(res_hrze_post),rownames(res_ntz_post)),rownames(res_int))
phy_ht <-  subset_samples(phy_ht,drug %in% c("pretreatment","NTZ"))


phy_sel_ht <- prune_taxa(rownames(res_ntz_post[res_ntz_post$global_padj < 0.05,]),phy_ht)
mat <-data.frame(otu_table(phy_sel_ht))
df <- as.data.frame(sample_data(phy_sel_ht)[,c("drug","TB_status","type","pool_batch")]) #"days_ON_HRZE","TB_status"
#df$type <- NULL
df$drug <- factor(df$drug, levels = c("pretreatment","NTZ","HRZE"))
color_pal <- list(drug = c("pretreatment" = "#537c4a","HRZE"="#984ea3","NTZ" = "#e46983"))
match_tax_dt <- data.frame(tax_table(phy_sel_ht))
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

match_tax_dt <- match_tax_dt[match(rownames(mat),rownames(match_tax_dt)),]




ha1 = HeatmapAnnotation(Order = match_tax_dt$Order,
                        col = list(Order = order_col), which = "row")
ha_column = HeatmapAnnotation(Drug =  df$drug,
                              batch = df$pool_batch,
                              col=list(Drug = c("pretreatment" = "#537c4a","HRZE"="#984ea3","NTZ" = "#e46983")))
split_rows <-  match_tax_dt$Phylum
split_rows <- factor(split_rows, levels= unique(split_rows))
#colfunc <- colorRampPalette(c("purple", "yellow"))

library(RColorBrewer)
colfunc <- colorRampPalette(rev(brewer.pal(n = 7, name =
                                             "RdYlBu")))
#teddy_cols <- c('white','#fecc5c','#fd8d3c','#f03b20','#bd0026')

ht_cols <- colfunc(8)
# Remove underscores 
rownames(mat) <- gsub("_"," ",rownames(mat))
ht1 = Heatmap(as.matrix(mat), name = "log", column_title = NA, 
              top_annotation = ha_column,
              clustering_distance_rows = "euclidean",
              column_split = factor(df$drug),
              split = split_rows, row_gap = unit(2, "mm"),border = TRUE,
              row_title_gp = gpar(fontsize = 10,fontface = "bold"),
              row_title_rot = 0,
              left_annotation = ha1,
              #  right_annotation = ha2,
              row_names_side = "left", km=1, color_space = "LAB",
              #row_dend_side="right",
              col= ht_cols,
              show_parent_dend_line = F,
              #heatmap_legend_param = list(),
              #clustering_method_columns = "ward.D",
              width=2, 
              #show_column_names= T,
              row_names_gp = gpar(fontsize = 9),
              #cluster_columns = T,cluster_rows = T,
              row_names_max_width = max_text_width(rownames(mat), gp = gpar(fontsize = 12)),
              column_names_max_height = max_text_width(colnames(mat), gp = gpar(fontsize = 12)),
              show_column_names = "F",
              column_names_side = "top",
              na_col="white")

pdf(paste(fig_folder,"Heatmap_ntz.pdf",sep = "/"),height = 40, width = 15)
draw(ht1)
dev.off()


# Draw a heatmap for HRZE and NTZ both 

mat <- cpm(cnt_dt, log=TRUE)
# VST phyloseq for Heatmap 
phy_ht <- phy_pair
otu_table(phy_ht) <-  otu_table(mat, taxa_are_rows = T)
# common_otus <- intersect(rownames(res_hrze_post))
# all_sp <- union(union(rownames(res_hrze_post),rownames(res_ntz_post)),rownames(res_int))
phy_ht <-  subset_samples(phy_ht,drug %in% c("pretreatment","NTZ","HRZE"))


comb_sp <- union(rownames(res_hrze_post[res_hrze_post$global_padj < 0.05,]),
                rownames(res_ntz_post[res_ntz_post$global_padj < 0.05,]))
phy_sel_ht <- prune_taxa(comb_sp,phy_ht)
mat <-data.frame(otu_table(phy_sel_ht))

df <- as.data.frame(sample_data(phy_sel_ht)[,c("drug","TB_status","type","pool_batch")]) #"days_ON_HRZE","TB_status"
df$pool_batch <- NULL
#df$type <- NULL
df$drug <- factor(df$drug, levels = c("pretreatment","NTZ","HRZE"))
color_pal <- list(drug = c("pretreatment" = "#537c4a","HRZE"="#984ea3","NTZ" = "#e46983"))
match_tax_dt <- data.frame(tax_table(phy_sel_ht))
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

match_tax_dt <- match_tax_dt[match(rownames(mat),rownames(match_tax_dt)),]




ha1 = HeatmapAnnotation(Order = match_tax_dt$Order,
                        col = list(Order = order_col), which = "row")
ha_column = HeatmapAnnotation(Drug =  df$drug,
                              batch = df$pool_batch,
                              col=list(Drug = c("pretreatment" = "#537c4a","HRZE"="#984ea3","NTZ" = "#e46983")))

# Right annotation Significance

sig_dt <-  data.frame(Taxa = rownames(mat),HRZE = "NS",NTZ = "NS",stringsAsFactors = F)
sig_dt$HRZE[sig_dt$Taxa %in% rownames(res_hrze_post_sig)] <- "S"
sig_dt$NTZ[sig_dt$Taxa %in% rownames(res_ntz_post_sig)] <- "S"

#'#deebf7','#9ecae1','#3182bd'
sig_col <- c('#deebf7','#3182bd')
names(sig_col) <-  c("NS","S")

ha_right <- rowAnnotation(HRZE =  sig_dt$HRZE,
                             NTZ = sig_dt$NTZ,
                             col=list(HRZE = sig_col,
                                      NTZ = sig_col))


split_rows <-  match_tax_dt$Phylum
split_rows <- factor(split_rows, levels= unique(split_rows))
#colfunc <- colorRampPalette(c("purple", "yellow"))

library(RColorBrewer)
colfunc <- colorRampPalette(rev(brewer.pal(n = 7, name =
                                             "RdYlBu")))
#teddy_cols <- c('white','#fecc5c','#fd8d3c','#f03b20','#bd0026')
split_cols <- factor(df$type,levels = c("Pre HRZE","Pre NTZ","Post HRZE","Post NTZ"))


ht_cols <- colfunc(8)
# Remove underscores 
rownames(mat) <- gsub("_"," ",rownames(mat))

ht1 = Heatmap(as.matrix(mat), name = "log(CPM)", column_title = NA, 
              top_annotation = ha_column,
              right_annotation = ha_right,
              clustering_distance_rows = "euclidean",
              cluster_columns = T,
              cluster_column_slices = F,
              column_split = split_cols,
              split = split_rows, row_gap = unit(2, "mm"),border = TRUE,
              row_title_gp = gpar(fontsize = 10,fontface = "bold"),
              row_title_rot = 0,
              left_annotation = ha1,
              #  right_annotation = ha2,
              row_names_side = "left", km=1, color_space = "LAB",
              #row_dend_side="right",
              col= ht_cols,
              show_parent_dend_line = F,
              #heatmap_legend_param = list(),
              #clustering_method_columns = "ward.D",
              width=2, 
              #show_column_names= T,
              row_names_gp = gpar(fontsize = 9),
              #cluster_columns = T,cluster_rows = T,
              row_names_max_width = max_text_width(rownames(mat), gp = gpar(fontsize = 12)),
              column_names_max_height = max_text_width(colnames(mat), gp = gpar(fontsize = 12)),
              show_column_names = "F",
              column_names_side = "top",
              na_col="white")

pdf(paste(fig_folder,"Heatmap_ntz_hrze_sig_combined.pdf",sep = "/"),height = 40, width = 20)
draw(ht1)
dev.off()


