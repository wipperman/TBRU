#Load libraries
library(DESeq2)
library(phyloseq)
library(dplyr)
library(yingtools2)
library(data.table)
library(ifultools)
library(stringr)
library(ggplot2)


# Save all the data into RDS files for easy import


################## TTP data HRZE and NTZ #############################

# Read ttp file
ttp_data <-  read.csv("../raw_data/Mic_data/metadata/TTP_results.csv",stringsAsFactors = F)
ttp_age <-  read.csv("../raw_data/Mic_data/metadata/age_treatment.csv",stringsAsFactors = F)
names(ttp_data)

# Add age to the data
ttp_data <- merge(ttp_data,ttp_age,by = "PID")

head(ttp_data)
names(ttp_data)
ttp_data <-  ttp_data[,c("PID","Event.Name" ,"Treatment.Group","Calculated.TTP1",                          
                         "Calculated.TTP2","Average.TTP","AGE" )]
names(ttp_data) <- c("ID","Day","Treatment","TTP1","TTP2","Av_TTP","Age")


ttp_data <- ttp_data[- which(is.na(ttp_data[,'TTP1']) & is.na(ttp_data[,'TTP2'])),]
ttp_data$Av_TTP[is.na(ttp_data[,'TTP1'])] <- ttp_data$TTP2[is.na(ttp_data[,'TTP1'])]
ttp_data$Av_TTP[is.na(ttp_data[,'TTP2'])] <- ttp_data$TTP1[is.na(ttp_data[,'TTP2'])]
ttp_data$ID <- factor(ttp_data$ID, unique(sort(ttp_data$ID)))
ttp_data$Day <-  gsub("Day ","",ttp_data$Day)
ttp_data$Day <- as.numeric(ttp_data$Day)

# Remove Day 5
# ttp_data <- ttp_data[ttp_data$Day != 5,]
# #ttp_data$Day <- ttp_data$Day + 1
# ttp_data <- ttp_data[ttp_data$Day != -1,]
ttp_data$Treatment[ttp_data$Treatment != "NTZ"] <- "HRZE"
ttp_data$Treatment <-   factor(ttp_data$Treatment , levels = c("NTZ","HRZE"))

# ttp_data_sub<-ttp_data[c(which(ttp_data$Day==0), 
#                          which(ttp_data$Day==14)),]

ttp_data$TTP_max <- pmax(ttp_data$TTP1,ttp_data$TTP2,na.rm = T)
ttp_data$TTP_min <- pmin(ttp_data$TTP1,ttp_data$TTP2,na.rm =T)
ttp_data
# #names(ttp_data_sub)

ttp_data_pre_post<-ttp_data[c(which(ttp_data$Day==0), 
                          which(ttp_data$Day==14)),]
ttp_data_pre_post$Treatment <- NULL

####################Microbiome##########################################

load("../raw_data/Mic_data/dada2/dada2_full_pipeline.RData")
phy.dada2
sample_names(phy.dada2)
snames <-  sample_names(phy.dada2)
snames_rm_pool <- gsub("..pool.*","",snames)

# Which samples are duplicate?
snames_rm_pool[duplicated(snames_rm_pool)]

# Check the number of reads/counts in each duplicate samples
sample_sums(phy.dada2)[grep("E7V05ZLS",names(sample_sums(phy.dada2)))]

# Remove pool669 as it has less reads that pool 652
phy.dada2 <- prune_samples(sample_names(phy.dada2)[sample_names(phy.dada2) != "E7V05ZLS.01..pool669.mw"], phy.dada2)

# #readTree 
phy <- phy.dada2


############ add sample data ##############
meta_dt <- read.csv("../raw_data/Mic_data/metadata/TBRU_NTZ_RNAseq_metadata.csv",header = T,sep = ",",na.strings=c("N/A","99999")) %>% mutate(sample=gsub("\\-",".",sample)) 
meta_dt <-  meta_dt[meta_dt$TB_status %in% c("pretreatment","HRZE","NTZ"),]
meta_dt <- meta_dt[meta_dt$sample != "",]
meta_dt$sample_name <-  gsub("\\..*","",meta_dt$sample)
#rownames(meta_dt) <- meta_dt$sample
meta_dt$sample <-  NULL

#data$sample_name <- rownames(data)
meta_dt$Day <- ifelse(meta_dt$TB_status == "pretreatment", 0, 14)


# Merge it with the pool information
# Pool information
pool_dt <-  read.csv("../raw_data/Mic_data/metadata/pool.samples.Feb2020.csv")
mer_dt <- merge(meta_dt, pool_dt,by.x = "sample_name", by.y = "short_sample_name")
unique(mer_dt$pool_batch)

sm_dt <- unique(merge(mer_dt,ttp_data_pre_post,
                        by.x = c("Patient.ID","Day"),
                        by.y = c("ID","Day"),
                        all.x = T ))
names(sm_dt)

rownames(sm_dt) <-  sm_dt$sample_name

names(sm_dt)

snames <-sample_names(phy)
snames <- gsub("^X.","",snames) 
snames <- sub("\\.","_",snames) 
snames <- gsub("\\..*","",snames) 

snames_rm_pool <- gsub("..pool.*","",snames) %>% base:::make.names(unique = T)
# Which samples are duplicate?
snames_rm_pool[duplicated(snames_rm_pool)]
snames_rm_pool <- gsub("_01","",snames_rm_pool) 

sample_names(phy) <-  snames_rm_pool
sample_data(phy) <- sm_dt %>% phyloseq::sample_data()
phy

tax_blast <- read.blastn.file("../raw_data/Mic_data/dada2/asv_seqs.fasta.blastn.refseq_rna.txt")

head(tax_blast)
tax_blast$otu <-  gsub(";$","",tax_blast$otu) 

tax_table(phy) <- tax_blast %>% set.tax()
phy <- prune_taxa(taxa_sums(phy)> 0, phy)


#view the rownames as Phylum, Species, rather than OTU_number
t <- get.tax(phy) %>% mutate(PhySpec=paste(Species,otu))
taxa_names(phy) <- t$PhySpec
taxa_names(phy) <- gsub("\\[","",taxa_names(phy))
taxa_names(phy) <- gsub("\\]","",taxa_names(phy))
taxa_names(phy) <- gsub(';.*','\\1', taxa_names(phy))
taxa_names(phy) <- gsub("ASV","(ASV",taxa_names(phy))
taxa_names(phy) <- gsub("$",")",taxa_names(phy)) #edit this between old/new methods
taxa_names(phy) <- gsub(" ","_",taxa_names(phy)) #edit this between old/new methods
taxa_names(phy)


# This is the phyloseq object for the microbiome data
# Let's save it
# Create a data folder
dir.create("../data/mic", showWarnings = TRUE,recursive = TRUE)

saveRDS(phy, "../data/mic/phy_mic_orig.rds")


# Prevalence 5 % 
perc_samples <-  0.05
phy_fil = filter_taxa(phy, function(x) sum(x > 0) > (perc_samples*length(x)), TRUE)
phy_fil <- prune_taxa(taxa_sums(phy_fil)> 0, phy_fil)


ntz_id <- unique(as.character(get_variable(subset_samples(phy_fil,drug == "NTZ" ), "Patient.ID")))
hrze_id <- unique(as.character(get_variable(subset_samples(phy_fil,drug == "HRZE" ), "Patient.ID")))

# Also include the type of the sample in metadata 
# For example Pre HRZE and Pre NTZ samples
sm_dt <-  data.frame(sample_data(phy_fil))
sm_dt$type <- as.character(sm_dt$TB_status)

# Separate pre-treatment into Pre HRZE and Pre NTZ 
sm_dt$type[sm_dt$Patient.ID %in% hrze_id & sm_dt$TB_status == "pretreatment"] <- "Pre HRZE"
sm_dt$type[sm_dt$Patient.ID %in% ntz_id & sm_dt$TB_status == "pretreatment"] <- "Pre NTZ"
sm_dt$type[sm_dt$TB_status == "HRZE"] <- "Post HRZE"
sm_dt$type[sm_dt$TB_status == "NTZ"] <- "Post NTZ"
sample_data(phy_fil) <- sm_dt
phy_fil

# Save the filtered phyloseq object with prevalence of 5 %
saveRDS(phy_fil, paste0("../data/mic/phy_mic_fil_prev_",perc_samples,".rds"))

# Use only paired samples:
# Make a frequency table and select only the paired samples for each Patient.ID
cnts <- data.frame(table(as.character(get_variable(phy_fil, "Patient.ID"))))
cnts <- cnts[cnts$Freq > 1,]
names(cnts)[1] <-  "Patient.ID"

# Only select the paired samples 
phy_pair <-  subset_samples(physeq = phy_fil, Patient.ID %in% cnts$Patient.ID)
table(get_variable(phy_pair,"TB_status"))
sample_data(phy_pair)$TB_status <- relevel(sample_data(phy_pair)$TB_status, ref = "pretreatment")
sample_data(phy_pair)$Patient.ID 
phy_pair <- prune_taxa(taxa_sums(phy_pair)>0,phy_pair)

# Save the filtered phyloseq object with prevalence of 5 %
saveRDS(phy_pair, "../data/mic/phy_mic_fil_paired.rds")


############# RNA-Seq data ###########################
phy_mic_orig  <- readRDS("../data/mic/phy_mic_orig.rds") 
phy_mic_pair <- readRDS("../data/mic/phy_mic_fil_paired.rds")

############## read RNAseq data ###########
rna_files = list.files(path = "../raw_data/rnaseq_data/counts",
                       pattern="*.samples.txt",
                       full.names = T)
rna_list = lapply(rna_files, read.table,header = T,sep = "\t")

# Merge into one big dataframe using a "cool" Reduce function 
full_rna_raw <- Reduce(merge, rna_list)
full_rna_raw <- full_rna_raw %>% mutate(GeneID=paste(GeneID,"_",GeneSymbol,sep = "")) 
full_rna_raw$GeneSymbol <- NULL
rownames(full_rna_raw) <- full_rna_raw$GeneID
full_rna_raw$GeneID <- NULL
rna_data_norm_genes <- full_rna_raw

# There is a sample which has the "_b", needs to be removed 
names(rna_data_norm_genes) <- gsub("_b","",names(rna_data_norm_genes))

# Making sure all the samples that are in microbiome matches the RNASeq samples
seqid <- as.character(get_variable(phy_mic_orig,"rna_seq_ID")) 
seqid[!seqid %in% names(rna_data_norm_genes)]

# Create a random taxa table for RNASeq data , ultimately making a phyloseq object
rna_taxmat <- matrix(sample(letters, nrow(rna_data_norm_genes)/10, replace = TRUE), nrow = nrow(rna_data_norm_genes), ncol = 7)
rownames(rna_taxmat) <- rownames(rna_data_norm_genes)
colnames(rna_taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rna_taxmat <-tax_table(rna_taxmat)
rna_otu <- otu_table(rna_data_norm_genes,taxa_are_rows = TRUE)
rna_phy <- phyloseq(rna_otu, rna_taxmat)
rna_phy

# Create a data folder
dir.create("../data/genes", showWarnings = TRUE,recursive = TRUE)

saveRDS(rna_phy,"../data/genes/phy_rna_all.rds")


# Include batch information for samples: assuming each file regarded as one batch
batch_num <-  unlist(lapply(rna_list, function(x) length(x)))
batch_num <-  batch_num - 2 # Removing the GeneID and GeneSymbols columns
batchname <-  paste0("batch",1:length(batch_num))
batchID <-  rep(batchname, times = batch_num)

# Make a sample table including the batch information
samp_rnaphy <- data.frame(samName =  names(rna_data_norm_genes), batchID = batchID)
samp <- data.frame(sample_data(phy_mic_orig))
samp$Sam_name <- rownames(samp)

# Merge the sample table from microbiome and newly created batch information 
samp <- merge(x = samp, y = samp_rnaphy, by.x = "rna_seq_ID", by.y = "samName")
rownames(samp) <- samp$rna_seq_ID
samp$pool_batch <- NULL
samp$number_after_name <- NULL

sample_data(rna_phy) <- phyloseq::sample_data(samp)

# # Change the sample names to match the microbiome sample names 
# sample_names(rna_phy) <- as.character(get_variable(rna_phy,"Sam_name"))
# Adding TTP information
sm_dt <-  data.frame(sample_data(rna_phy))


# Change the sample names to match the microbiome sample names 
sample_names(rna_phy) <- as.character(get_variable(rna_phy,"Sam_name"))

# Filter the genes that have total counts > 1
rna_phy <-  prune_taxa(taxa_sums(rna_phy)>1,rna_phy)
sample_names(rna_phy)

sample_names(phy_mic_orig)

# Save the phyloseq object for rna seq
saveRDS(rna_phy,"../data/genes/phy_rna_ntz_hrze.rds")


rna_phy_pair <-  prune_samples(sample_names(phy_mic_pair),rna_phy)
# Save the phyloseq object for paired rna seq
saveRDS(rna_phy_pair,"../data/genes/phy_rna_ntz_hrze_pair.rds")


######EBA Cohort###############################


####################Microbiome##########################################
load("../raw_data/Mic_data/dada2/dada2_full_pipeline.RData")
phy.dada2
sample_names(phy.dada2)
snames <-  sample_names(phy.dada2)
snames_rm_pool <- gsub("..pool.*","",snames)

# Which samples are duplicate?
snames_rm_pool[duplicated(snames_rm_pool)]

# Check the number of reads/counts in each duplicate samples
sample_sums(phy.dada2)[grep("E7V05ZLS",names(sample_sums(phy.dada2)))]

# Remove pool669 as it has less reads that pool 652
phy.dada2 <- prune_samples(sample_names(phy.dada2)[sample_names(phy.dada2) != "E7V05ZLS.01..pool669.mw"], phy.dada2)

# #readTree 
tree <- read_tree("../raw_data/Mic_data/dada2/total.10.tree")
# 
phy_tree(phy.dada2) <- tree
phy.dada2

############ add sample data ##############
data <- read.csv("../raw_data/Mic_data/metadata/stool_stats_25mar2020.csv",header = T,sep = ",",na.strings=c("N/A","99999")) %>% mutate(sample=gsub("\\-",".",sample)) 

#data <-  data[data$TB_status %in% c("pretreatment","HRZE","NTZ"),]
#data <- data[data$rna_seq_ID != "",]
data <- data[data$sample != "",]
data$sample <- gsub("0487.","",data$sample)
data$sample <- gsub(".001","",data$sample)

rownames(data) <- data$sample
data$sample <-  NULL

snames <-sample_names(phy.dada2)
snames <- gsub("\\.\\.pool.*","",snames) 
snames <- gsub("0487.","",snames)
snames <- gsub(".001","",snames)

snames[nsamples(phy.dada2)] <- paste0("00",gsub("\\.","",snames[nsamples(phy.dada2)]))


rownames(data)[!rownames(data) %in% snames]


phy <- phy.dada2
sample_names(phy) <-  snames
sample_data(phy) <- data %>% phyloseq::sample_data()
phy
sample_names(phy)

tax_blast <- read.blastn.file("../raw_data/Mic_data/dada2/asv_seqs.fasta.blastn.refseq_rna.txt")

head(tax_blast)
tax_blast$otu <-  gsub(";$","",tax_blast$otu) 

tax_table(phy) <- tax_blast %>% set.tax()
phy

phy <- prune_taxa(taxa_sums(phy)> 0, phy)
phy

# Edit taxa
t <- get.tax(phy) %>% mutate(PhySpec=paste(Species,otu))
taxa_names(phy) <- t$PhySpec
taxa_names(phy) <- gsub("\\[","",taxa_names(phy))
taxa_names(phy) <- gsub("\\]","",taxa_names(phy))
taxa_names(phy) <- gsub(';.*','\\1', taxa_names(phy))
taxa_names(phy) <- gsub("ASV","(ASV",taxa_names(phy))
taxa_names(phy) <- gsub("$",")",taxa_names(phy)) #edit this between old/new methods
taxa_names(phy) <- gsub(" ","_",taxa_names(phy)) #edit this between old/new methods
taxa_names(phy)

sample_names(phy) <-  make.names(sample_names(phy))



# Fix the sample data
# Change Event.Name as Time
sm_dt <-  data.frame(sample_data(phy))
sm_dt$sam_id <- rownames(sm_dt)
names(sm_dt)[2] <-  "Time"
sm_dt$Time <-  gsub(" ","",sm_dt$Time)

sm_dt$Time <- factor(sm_dt$Time,c("Day0","Day7","Day14","OneMonth",
                                  "TwoMonths","SixMonthFollowup" ))
#sample_data(phy_trt) <-  sample_data(sm_dt)


sm_dt$Time

sm_dt <- sm_dt %>%
  arrange(Patient.ID,Time) %>%
  as.data.frame()

sm_dt$ind <- as.numeric(factor(as.character(sm_dt$Patient.ID),
                               levels = unique(as.character(sm_dt$Patient.ID))))
rownames(sm_dt) <- sm_dt$sam_id

sample_data(phy) <- sample_data(sm_dt)

dir.create("../data/EBA",showWarnings = T)
saveRDS(phy, "../data/EBA/phy_mic_eba.rds")

####################### Read rna seq data########### 


rna_dat <- read.csv("../raw_data/EBA_data/counts/EBA_raw_counts_54_10feb2020.txt",sep = ",")

names(rna_dat)[1:2] <- c("GeneID","GeneSymbol") 
rownames(rna_dat) <- paste0(rna_dat$GeneID,"_",rna_dat$GeneSymbol)

met_dat <-  rna_dat[,c("GeneID","GeneSymbol")]

rna_dat$GeneID <- NULL
rna_dat$GeneSymbol <-  NULL

ttp_dat <-  read.csv("../raw_data/EBA_data/metadata/rna_metadata.csv",sep = ",")
ttp_dat$Day <-  ttp_dat$DAYS_ON_TREATMENT
ttp_dat$Day <-  paste0("D",ttp_dat$Day)
ttp_dat$Average.TTP <-  as.numeric(as.character(ttp_dat$TTP_CONTINUOUS))
ttp_dat$sam_id <-  ttp_dat$SAMPLE
rownames(ttp_dat) <-  ttp_dat$sam_id

# Change Event.Name as Time
names(ttp_dat)[20] <-  "Time"
names(ttp_dat)[3] <-  "Patient.ID"
#ttp_dat$Time <-  gsub(" ","",ttp_dat$Time)


colnames(rna_dat)

ttp_dat$SAMPLE
otu_dat <- otu_table(rna_dat, taxa_are_rows = T)
phy_rna <-  phyloseq(otu_dat, tax_table(as.matrix(met_dat)))
sample_names(phy_rna) <- gsub("s_","",sample_names(phy_rna))
sample_data(phy_rna) <-  sample_data(ttp_dat)


phy_rna <-  prune_taxa(taxa_sums(phy_rna)> 0, phy_rna)

phy_rna <- prune_samples(sample_names(phy_rna)[-grep("bu",sample_names(phy_rna))], phy_rna)


# Prevalence 5 % 
perc_samples <-  0.00
# 
phy_fil = filter_taxa(phy_rna, function(x) sum(x > 0) > (perc_samples*length(x)), TRUE)
phy_fil <- prune_taxa(taxa_sums(phy_fil)> 0, phy_fil)

# Remove the patient 1830014T : withdrew from study
phy_fil <-  subset_samples(phy_fil, Patient.ID != "1830014T")

saveRDS(phy_fil, "../data/EBA/phy_rna_tb.rds")




####################Family Contact and Community Control##############

# Make a family contact cohort
load("../raw_data/Mic_data/dada2/dada2_full_pipeline.RData")
phy.dada2
sample_names(phy.dada2)
snames <-  sample_names(phy.dada2)


met_data <-  read.csv("../raw_data/Mic_data/metadata/TBRU_NTZ_RNAseq_metadata.csv")
met_data$type <-  met_data$TB_status
met_data <- met_data[met_data$type %in% c("family_contact","community_control"),]
met_data <-  met_data[!met_data$rna_seq_ID == "",]
met_data <-  met_data[!met_data$sample == "",]

met_data$sample <- gsub("-","_",met_data$sample)

phy_mic <- phy.dada2

met_data$sample

snames <-sample_names(phy_mic)
snames <- gsub("^X.","",snames) 
snames <- sub("\\.","_",snames) 
snames <- gsub("\\..*","",snames) 
snames <- gsub("\034","",snames) 


met_data$sample[!met_data$sample %in% snames]

id_keep  <- which(snames %in% met_data$sample)

sample_names(phy_mic)[id_keep] <- snames[id_keep]


phy_fc_cc <-  prune_samples(met_data$sample,phy_mic)
rownames(met_data) <-  met_data$sample
sample_data(phy_fc_cc)<- sample_data(met_data)

phy_fc_cc
phy_fc_cc <- prune_taxa(taxa_sums(phy_fc_cc)>0,phy_fc_cc)
sample_names(phy_fc_cc) <-  sample_data(phy_fc_cc)$rna_seq_ID

dir.create("../data/Control",showWarnings = T)
saveRDS(phy_fc_cc, "../data/Control/phy_mic_fc_cc.rds")

############## read RNAseq data ###########
rna_files = list.files(path = "../raw_data/rnaseq_data/counts",
                       pattern="*.samples.txt",
                       full.names = T)
rna_list = lapply(rna_files, read.table,header = T,sep = "\t")

# Merge into one big dataframe using a "cool" Reduce function 
full_rna_raw <- Reduce(merge, rna_list)
full_rna_raw <- full_rna_raw %>% mutate(GeneID=paste(GeneID,"_",GeneSymbol,sep = "")) 
full_rna_raw$GeneSymbol <- NULL
rownames(full_rna_raw) <- full_rna_raw$GeneID
full_rna_raw$GeneID <- NULL
rna_data_norm_genes <- full_rna_raw

# There is a sample which has the "_b", needs to be removed 
names(rna_data_norm_genes) <- gsub("_b","",names(rna_data_norm_genes))


# Create a random taxa table for RNASeq data , ultimately making a phyloseq object
rna_taxmat <- matrix(sample(letters, nrow(rna_data_norm_genes)/10, replace = TRUE), nrow = nrow(rna_data_norm_genes), ncol = 7)
rownames(rna_taxmat) <- rownames(rna_data_norm_genes)
colnames(rna_taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rna_taxmat <-tax_table(rna_taxmat)
rna_otu <- otu_table(rna_data_norm_genes,taxa_are_rows = TRUE)
rna_phy <- phyloseq(rna_otu, rna_taxmat)
rna_phy

# Include batch information for samples: assuming each file regarded as one batch
batch_num <-  unlist(lapply(rna_list, function(x) length(x)))
batch_num <-  batch_num - 2 # Removing the GeneID and GeneSymbols columns
batchname <-  paste0("batch",1:length(batch_num))
batchID <-  rep(batchname, times = batch_num)


# Make a sample table including the batch information
samp_rnaphy <- data.frame(samName =  names(rna_data_norm_genes), batchID = batchID)
samp <- data.frame(sample_data(phy_fc_cc))
samp$Sam_name <- rownames(samp)

# Merge the sample table from microbiome and newly created batch information 
samp <- merge(x = samp, y = samp_rnaphy, by.x = "rna_seq_ID", by.y = "samName")
rownames(samp) <- samp$rna_seq_ID
sample_data(rna_phy) <- phyloseq::sample_data(samp)

rna_phy_control <- prune_taxa(taxa_sums(rna_phy)>0,rna_phy)
saveRDS(rna_phy_control, "../data/Control/phy_rna_fc_cc.rds")


sm_dt <-  data.frame(sample_data(phy_gene_sel))
dds <- phyloseq_to_deseq2(rna_phy, ~ sex + TB_status ) #replace this with any sample variable(s)

dds$TB_status
dds$sex
dds$batchID


dds <- estimateSizeFactors(dds)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
table(keep)

# Takes a while to run
dds <- estimateDispersions(dds)
dds <- DESeq(dds,fitType= "local")


plotDispEsts(dds)

# VST phyloseq for Heatmap later

phy_vst_rna <- rna_phy_control
vst_dt <- getVarianceStabilizedData(dds)
otu_table(phy_vst_rna) <-  otu_table(vst_dt, taxa_are_rows = T)
saveRDS(phy_vst_rna,"../data/Control/phy_rna_vst_fc_cc.rds")

