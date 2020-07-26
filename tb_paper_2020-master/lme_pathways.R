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

mainDir <- "../Tables"
subDir <- "Pathways_LME_EBA_Clinical"
#subDir <- "Figs_Clostridia_kegg"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
tab_folder <- paste(mainDir,subDir,sep="/")



####### load phyloseq objected created with data_setup.R ######

phy <- readRDS("../data/genes/phy_rna_ntz_hrze_pair.rds")

sample_names(phy) <- sample_data(phy)$rna_seq_ID

phy_pth <- readRDS("../data/GSVA/phy_GSVA_all_hallmark.rds")
sample_names(phy_pth)


# tmp <- readRDS("../../EBA_HRZE_NTZ/Figs_vb/GSVA/phy_GSVA_all.rds")
# 
# otu1 <- data.frame(otu_table(phy_kegg))
# 
# otu2 <-  data.frame(otu_table(tmp))

phy_pth <- prune_samples(sample_names(phy),phy_pth)
sample_data(phy_pth) <- sample_data(phy)


phy_pair <- phy_pth

sm_dt <- data.frame(sample_data(phy_pair))

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
# sm_dt$ind.n <- sm_dt$ind
# sm_dt$ind.n[sm_dt$ind.n>= 17] <-  sm_dt$ind.n[sm_dt$ind.n>= 17] - 16
# sm_dt$ind.n <- factor(sm_dt$ind.n)
rownames(sm_dt) <- sm_dt$SampleID

sample_data(phy_trt) <- sample_data(sm_dt)

table(get_variable(phy_trt,"Time"))
table(get_variable(phy_trt,"Treatment"))
sample_data(phy_trt)$Patient.ID




melt_dt <- psmelt(phy_trt)

# Let's do LME to check which pathways are differentiating pre and post in both treatments
names(melt_dt)

melt_dt$sex
melt_dt$batchID
melt_dt$drug <-  relevel(melt_dt$drug,"pretreatment")


#~  sex + batchID + drug  , sm_dt

result_lme <-  list()
pth_dt <- melt_dt

pthway <- as.character(unique(pth_dt$OTU))
pth_var <- pthway[1]

library(nlme)
for (pth_var in pthway){
  sum_mod_dt <- tryCatch({
    pth_sel_dt <- pth_dt[pth_dt$OTU %in% pth_var,]
    
    fml <- as.formula( paste( "Abundance", "~", "sex + batchID + drug " ) )
    
    mod_bac <- lme(fml,random = ~1|Patient.ID, pth_sel_dt)
    #mod_bac <- lme(fml,random = ~1|Subject.ID, ab_gene_sel)
    
    sum_mod <-  summary(mod_bac)
    sum_mod_dt <- data.frame(sum_mod$tTable)
    sum_mod_dt$pth_var <- pth_var
    sum_mod_dt$Var <-  rownames(sum_mod_dt)
    sum_mod_dt
    #}, 
    #warning = function(war){
    #  sum_mod_dt <- NA
    #  sum_mod_dt
  }, error = function(err){
    sum_mod_dt <- NA
    sum_mod_dt
  }
  )
  result_lme[[pth_var]] <-  sum_mod_dt
}

result_lme_dt <- do.call("rbind", result_lme)
names(result_lme_dt)[5] <- "pval"
library(gtools)
# Remove intercept
result_lme_dt <- result_lme_dt[result_lme_dt$Var != "(Intercept)",]
result_lme_dt <- result_lme_dt[result_lme_dt$Var %in% c("drugHRZE","drugNTZ"),]

result_lme_dt$p_adj <- p.adjust(result_lme_dt$pval, method = "BH")
# result_lme_dt$p_adj <- result_lme_dt$pval
# result_lme_dt <- result_lme_dt[order(result_lme_dt$p_adj,decreasing = F),]
# 
result_lme_dt <- result_lme_dt[,c("Value","pth_var","Var","p_adj" )]
result_lme_dt <- na.omit(result_lme_dt)
dt_w <- reshape(result_lme_dt, idvar = "pth_var", timevar = "Var", direction = "wide")
dt_w <- dt_w[,mixedsort(names(dt_w))]
rownames(dt_w) <- dt_w$pth_var

write.csv(dt_w,paste0(tab_folder,"/LME_Clinical_Post_Pre.csv"))


hrze_diff <-  dt_w$pth_var[dt_w$p_adj.drugHRZE <= 0.05]
ntz_diff <-  dt_w$pth_var[dt_w$p_adj.drugNTZ <= 0.05]

union_clin_diff <-  union(hrze_diff, ntz_diff)



# For EBA trial

phy_eba <- readRDS("../data/EBA/phy_rna_tb.rds")

phy_pth <- readRDS("../data/GSVA/phy_GSVA_all_hallmark.rds")
sample_names(phy_pth)

phy_pth <- prune_samples(sample_names(phy_eba),phy_pth)
sample_data(phy_pth) <- sample_data(phy_eba)

phy_trt <-  phy_pth
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


melt_dt <- psmelt(phy_trt)

# Let's do LME to check which pathways are differentiating pre and post in both treatments
names(melt_dt)

melt_dt$SEX
melt_dt$EXTRACTION
melt_dt$Time <-  factor(melt_dt$Time,levels = c("D0","D14","D56"))


#~  sex + batchID + drug  , sm_dt
#SEX + EXTRACTION + Time 
result_lme <-  list()
pth_dt <- melt_dt

pthway <- as.character(unique(pth_dt$OTU))
pth_var <- pthway[1]

library(nlme)
for (pth_var in pthway){
  sum_mod_dt <- tryCatch({
    pth_sel_dt <- pth_dt[pth_dt$OTU %in% pth_var,]
    
    fml <- as.formula( paste( "Abundance", "~", "SEX + EXTRACTION + Time " ) )
    
    mod_bac <- lme(fml,random = ~1|Patient.ID, pth_sel_dt)
    #mod_bac <- lme(fml,random = ~1|Subject.ID, ab_gene_sel)
    
    sum_mod <-  summary(mod_bac)
    sum_mod_dt <- data.frame(sum_mod$tTable)
    sum_mod_dt$pth_var <- pth_var
    sum_mod_dt$Var <-  rownames(sum_mod_dt)
    sum_mod_dt
    #}, 
    #warning = function(war){
    #  sum_mod_dt <- NA
    #  sum_mod_dt
  }, error = function(err){
    sum_mod_dt <- NA
    sum_mod_dt
  }
  )
  result_lme[[pth_var]] <-  sum_mod_dt
}

result_lme_dt <- do.call("rbind", result_lme)
names(result_lme_dt)[5] <- "pval"
library(gtools)
# Remove intercept
result_lme_dt <- result_lme_dt[result_lme_dt$Var != "(Intercept)",]
result_lme_dt <- result_lme_dt[result_lme_dt$Var %in% c("TimeD14","TimeD56"),]

result_lme_dt$p_adj <- p.adjust(result_lme_dt$pval, method = "BH")
# result_lme_dt$p_adj <- result_lme_dt$pval
# result_lme_dt <- result_lme_dt[order(result_lme_dt$p_adj,decreasing = F),]
# 
result_lme_dt <- result_lme_dt[,c("Value","pth_var","Var","p_adj" )]
result_lme_dt <- na.omit(result_lme_dt)
dt_w <- reshape(result_lme_dt, idvar = "pth_var", timevar = "Var", direction = "wide")
dt_w <- dt_w[,mixedsort(names(dt_w))]
rownames(dt_w) <- dt_w$pth_var


write.csv(dt_w,paste0(tab_folder,"/LME_EBA_Post_Pre.csv"))


two_week_diff <-  dt_w$pth_var[dt_w$p_adj.TimeD14 <= 0.05]
#two_month_diff <-  dt_w$pth_var[dt_w$p_adj.TimeD56 <= 0.05]

#union_eba_diff <-  union(two_week_diff, two_month_diff)


# all_union <- union(union_eba_diff,union_clin_diff)

all_union <- union(two_week_diff,union_clin_diff)
saveRDS(all_union,"../data/GSVA/lme_eba_clinical_two_week.rds")

#phy_sel_kegg <- readRDS("../data/GSVA/phy_GSVA_all_hallmark.rds")
#phy_sel_kegg <- prune_taxa(all_union, phy_sel_kegg)

#saveRDS(phy_sel_kegg,"../../Figs_vb/GSVA/phy_GSVA_lme_kegg.rds")
