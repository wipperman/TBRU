#Load libraries
library(DESeq2)
library(phyloseq)
library(dplyr)
library(yingtools2)
library(data.table)
library(ifultools)
library(stringr)
library(ggplot2)
library(gtools)
library(ggrepel)



# Create a directory to save figures and tables

mainDir <- "../Figs"
subDir <- "Figs_2"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
fig_folder <- paste(mainDir,subDir,sep="/")

mainDir <- "../Tables"
subDir <- "Tabs_2"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
tab_folder <- paste(mainDir,subDir,sep="/")


########## Read Microbiome data #################################################

phy_pair <- readRDS("../data/mic/phy_mic_fil_paired.rds")

# Sample data
sm_dt <- data.frame(sample_data(phy_pair))

sm_dt$Treatment <- "HRZE"
sm_dt$Treatment[grep("NTZ",sm_dt$type)] <- "NTZ"


# Now we are ready to analyze the data

############## Fig 2  #########################

pal_Prt_HRZE_NTZ <- c("pretreatment" = "#537c4a",
                      "NTZ" = "#e46983",
                      "HRZE" = "#984ea3") 
sm_dt$Treatment <- factor(sm_dt$Treatment,levels = c("NTZ","HRZE"))


p <- ggplot()+
  geom_errorbar(data=sm_dt, aes(x = Day, ymin = TTP_min, ymax= TTP_max, color= 
                                    Treatment, group = Patient.ID), width=.4, position=position_dodge(.5),
                alpha = 0.9)+
  geom_point(data=sm_dt, aes(x = Day, y=Av_TTP, color= Treatment,
                             group = Patient.ID, fill = Treatment), 
             position=position_dodge(.5), size=2, shape = 19, alpha =0.9)+
  geom_line(data=sm_dt, aes(x = Day, y=Av_TTP, color= 
                                     Treatment, group = Patient.ID), 
            position=position_dodge(.5), 
            size=0.5, alpha = 0.9)+
  theme_classic()+
  scale_color_manual(values = pal_Prt_HRZE_NTZ)+
  scale_fill_manual(values = pal_Prt_HRZE_NTZ)+
  facet_wrap(~Treatment)+
  xlab("Day")+
  ylab("Time to positivity (TTP)")+
  theme(legend.text = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        #axis.text.x = element_text(size=15,angle=90, hjust=1),
        axis.text.x = element_text(size=15, hjust=1),
        axis.text.y = element_text(size=15))+
  scale_x_continuous(breaks=c(0,14),
                     labels=c("0","14"))  


pdf(paste(fig_folder,"Fig_2_B.pdf",sep = "/"), height = 6,width = 8,useDingbats = F)
print(p)
dev.off()

# LME Model 

# Model TTP ~ Treatment + Day + Treatment:Day
# Day -  0, 14
# Treatment - HRZE NTZ 
library(nlme)
lme_dt <-  sm_dt
lme_dt$Day <- factor(lme_dt$Day, levels = c(0,14))
lme_dt$ID <- factor(lme_dt$Patient.ID)
lme_dt$Treatment <- relevel(lme_dt$Treatment ,ref = "NTZ")
lme_dt$Av_TTP <- log(lme_dt$Av_TTP)
names(lme_dt)

lme_dt <- lme_dt[!is.na(lme_dt$Av_TTP),]

fit_ttp <-lme( Av_TTP ~ sex + Age + Treatment * Day ,
               random = ~ 1| ID,
               data = lme_dt)

sum_fit <- summary(fit_ttp)
dt_sum <-  data.frame(sum_fit$tTable)

# Save the LME results
write.csv(dt_sum,paste(tab_folder,"Tab_LME_TTP_2_B.csv",sep = "/"))





##### Principle component analysis #####
library(vegan)
phy_trt <-  phy_pair
sm_dt <- data.frame(sample_data(phy_trt))
ps <- phy_trt
# Take a log transform 
pslog <- transform_sample_counts(ps, function(x) log(1 + x))


#Well first perform a PCoA using Bray-Curtis dissimilarity.
set.seed(1057)
out.pcoa.log <- ordinate(pslog,  method = "PCoA", distance = "bray")
evals <- out.pcoa.log$values[,1]
p <- plot_ordination(pslog, out.pcoa.log,axes = c(1,2)) 

data_pca <- p$data

data_pca$Group <-  gsub("Pre |Post ","",data_pca$type)

pal_sam_type <-  c("Pre NTZ" = "#537c4a","Pre HRZE" = "#537c4a",
                   "Post NTZ" = "#e46983", "Post HRZE" = "#984ea3") 


# Fig 1 C
library(ggrepel)
p_pca <- ggplot(data_pca, aes(x =Axis.1, y = Axis.2,label = Patient.ID ))+ 
  geom_line(aes(group=Patient.ID),color="grey",size=0.5) + 
  geom_point(aes(color = type,group=Patient.ID, shape = drug),
             alpha = .8,size=5) + 
  # geom_text_repel()+
  theme_classic()+
  theme(legend.text = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15)) +
  #scale_color_manual(name = "",values=pal_Prt_HRZE_NTZ) +
  scale_color_manual(name = "",values=pal_sam_type) +
  
  scale_shape_manual(name = "",values = c("pretreatment" = 15,"HRZE" = 16,"NTZ" =  17))+
  xlab(p$labels$x)+
  ylab (p$labels$y)
pdf(paste(fig_folder,"PCoA_Fig_2_C.pdf",sep = "/"),width = 10,height = 10,useDingbats = F)
print(p_pca)
dev.off()



p_pca <- ggplot(data_pca, aes(x =Axis.1, y = Axis.2,label = pool_batch ))+ 
  #geom_line(aes(group=Patient.ID),color="grey",size=0.5) + 
  geom_point(aes(color = Group,
                 #group=Patient.ID,
                 shape = drug),
             alpha = .8,size=5) + 
  geom_text_repel()+
  theme_classic()+
  theme(legend.text = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15)) +
  scale_color_manual(name = "",values=pal_Prt_HRZE_NTZ) +
  #scale_color_manual(name = "",values=pal_sam_type) +
  
  scale_shape_manual(name = "",values = c("pretreatment" = 15,"HRZE" = 16,"NTZ" =  17))+
  xlab(p$labels$x)+
  ylab (p$labels$y)
#pdf(paste(fig_folder,"PCoA_Pool_diff.pdf",sep = "/"),width = 15,height = 20)
print(p_pca)


###########Diversity boxplot ###############
# Original phyloseq object without pruning
#readRDS("../../data/mic/phy_mic_fil_paired.rds")
raw_phy <-  phy_pair
alpha <- estimate_richness(raw_phy) #%>% mutate(sample=row.names(.))
data <- get.samp(raw_phy,stats = T) %>% as.data.frame()

data$type <-  factor(data$type , levels = c("Pre HRZE","Post HRZE","Pre NTZ","Post NTZ"))


# LME Model diversity
# Shannon Diversity
library(nlme)
div_dt  <- data
div_dt$Time <-   gsub(" HRZE| NTZ","",div_dt$type)
div_dt$Time <- factor(div_dt$Time, levels = c("Pre","Post"))
div_dt$ID <- factor(div_dt$Patient.ID)

div_dt$Treatment <-   gsub("Pre |Post ","",div_dt$type)
div_dt$Treatment <- factor(div_dt$Treatment, levels = c("HRZE","NTZ"))

lme_div_dt <- div_dt[!is.na(div_dt$Age),]

fit_div <-lme( Shannon ~ sex + Age + pool_batch +  Treatment * Time ,
               random = ~ 1| ID,
               data = lme_div_dt,na.action = )
sum_fit <- summary(fit_div)
dt_sum <-  data.frame(sum_fit$tTable)
dt_sum
# Save the lme results
write.csv(dt_sum,paste(tab_folder,"LME_Shannon_div_1_D.csv",sep = "/"))


div_dt_shan <- div_dt
div_dt_shan$Day <- ifelse(div_dt_shan$Time == "Pre",0,14)
div_dt_shan$Treatment <- relevel(div_dt_shan$Treatment,"NTZ")

div_dt_shan$Day_num = 0
div_dt_shan$Day_num[div_dt_shan$Time == "Post"] = 14

# Line plot sShannon
p_div <- ggplot() +
  geom_point(data=div_dt_shan, aes(x = Day_num, y=Shannon, color= 
                                Treatment, group = ID, fill = Treatment), size=3, shape = 19, alpha =0.7)+
  geom_line(data=div_dt_shan, aes(x = Day_num, y=Shannon, color= 
                               Treatment, group = ID), 
            size=0.5, alpha = 0.7)+
  stat_summary(data=div_dt_shan, aes(x = Day_num, y=Shannon, group = Day_num), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 4)+
  theme_classic()+
  scale_color_manual(values = pal_Prt_HRZE_NTZ)+
  scale_fill_manual(values = pal_Prt_HRZE_NTZ)+
  facet_wrap(~Treatment)+
  xlab("Day")+
  ylab("Shannon")+
  theme(legend.text = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=15, hjust=1),
        axis.text.y = element_text(size=15)) 
pdf(paste(fig_folder,"LinePlot_Shannon_diversity_Fig_2_D.pdf",sep = "/"),height = 6,width = 8, useDingbats = F) #paired sample analysis 
print(p_div)
dev.off()


# Inv Simpson 
fit_div <-lme( InvSimpson ~ sex + Age + pool_batch +  Treatment * Time ,
               random = ~ 1| ID,
               data = lme_div_dt,na.action = )
sum_fit <- summary(fit_div)
dt_sum <-  data.frame(sum_fit$tTable)
dt_sum
# Save the lme results
write.csv(dt_sum,paste(tab_folder,"LME_InvSimpson_div_1_D.csv",sep = "/"))


# Line plot InvSimpson
p_div <- ggplot() +
  geom_point(data=div_dt_shan, aes(x = Day_num, y=InvSimpson, color= 
                                     Treatment, group = ID, fill = Treatment), size=3, shape = 19, alpha =0.7)+
  geom_line(data=div_dt_shan, aes(x = Day_num, y=InvSimpson, color= 
                                    Treatment, group = ID), 
            size=0.5, alpha = 0.7)+
  stat_summary(data=div_dt_shan, aes(x = Day_num, y=InvSimpson, group = Day_num), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 4)+
  theme_classic()+
  scale_color_manual(values = pal_Prt_HRZE_NTZ)+
  scale_fill_manual(values = pal_Prt_HRZE_NTZ)+
  facet_wrap(~Treatment)+
  xlab("Day")+
  ylab("InvSimpson")+
  theme(legend.text = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=15, hjust=1),
        axis.text.y = element_text(size=15)) 
pdf(paste(fig_folder,"LinePlot_InvSimpson_diversity_Fig_2_D.pdf",sep = "/"),height = 6,width = 8, useDingbats = F) #paired sample analysis 
print(p_div)
dev.off()

############## Fig 2 complete #########################