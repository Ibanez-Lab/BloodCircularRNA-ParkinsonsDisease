# Load the necessary libraries -----
library(dplyr)
library(stringr)
library(reshape2)
library(flexiblas)
library(DESeq2)
library(readxl)
library("lmerTest")
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(metaRNASeq)
flexiblas_load_backend("__FALLBACK__")
# Define a "not in" opperator -----
`%!in%` <- Negate(`%in%`)
# I - Generate long tables --------------------------------------------------------------------
# 1. PDBP --------------------------------------------------------------------
# Load data -----
# load normalized count matrix
data <- read.table("PDBP/normalized.counts.txt", header = T, stringsAsFactors = F)
# load cleaned phenotype data
pheno <- read.table('PDBP/clean.pheno.txt', header = T, stringsAsFactors = F)
# Load normalized cell counts
cc<-read.csv("PDBP/normalized.cellcounts.csv")
# Format data -----
# Ensure that the time variable is numeric
str(pheno$visit_month)
pheno$time_in_years <- pheno$visit_month / 12
# select covariates in pheno that are needed for mixed models
mixed_model_pheno <- pheno %>% select(sample_id, sex, status = DiseaseStatus, age_at_baseline, time_in_years, PATNO = participant_id)
mixed_model_pheno$PATNO <- str_replace_all(mixed_model_pheno$PATNO, "PD-", "")
mixed_model_pheno <- mixed_model_pheno[order(mixed_model_pheno$PATNO, mixed_model_pheno$time_in_years), ]
# Merge phenotype and cell count data
mixed_model_pheno <- merge(mixed_model_pheno, cc, by = "sample_id")
# Transpose the count matrix
mixed_model_count_matirx <- data.frame(t(data))
mixed_model_count_matirx$sample_id <- rownames(mixed_model_count_matirx)
# Merge phenotype and count matrix
mixed_model_longtable <- merge(mixed_model_pheno, mixed_model_count_matirx, by = "sample_id")
# create the long table
mixed_model_longtable <- melt(mixed_model_longtable, id.vars = c('sample_id', "PATNO", 'sex', 'age_at_baseline', 'status', 'time_in_years','B_cell',
                                                                 'T_cell_CD4','T_cell_CD8','Monocyte','NK_cell','Neutrophil'))
# Change names of variable and value columns
colnames(mixed_model_longtable)[colnames(mixed_model_longtable) == "variable"] <- "circRNA"
colnames(mixed_model_longtable)[colnames(mixed_model_longtable) == "value"] <- "counts"
str(mixed_model_longtable)
# Ensure that circRNA name is a character variable
mixed_model_longtable$circRNA <- as.character(mixed_model_longtable$circRNA)
# Order rows by participant ID
mixed_model_longtable <- mixed_model_longtable[order(mixed_model_longtable$PATNO), ]
# remove 0 counts
mixed_model_longtable<-mixed_model_longtable[mixed_model_longtable$counts != 0, ]
mixed_model_longtable <- na.omit(mixed_model_longtable)
# add baseline counts
mixed_model_longtable$counts_at_baseline <- -9
mixed_model_longtable <- mixed_model_longtable %>% group_by(PATNO, circRNA) %>% mutate(counts_at_baseline = counts[which(time_in_years == min(time_in_years))])
mixed_model_longtable %>% filter(PATNO == rand_sample, circRNA == rand_circ)
# Save the long table
write.table(mixed_model_longtable, 'PDBP.longTable.EU.allVisits.txt', quote= F, sep = '\t', row.names = F, col.names = T)
# 2. PPMI --------------------------------------------------------------------
# Load data -----
# load normalized count matrix
data <- read.table("PPMI/normalized.counts.txt", header = T, stringsAsFactors = F)
# load cleaned phenotype data
pheno <- read.table('PPMI/clean.pheno.txt', header = T, stringsAsFactors = F)
# Load normalized cell counts
cc<-read.csv("PPMI/normalized.cellcounts.csv")
# Format data -----
# Ensure that the time variable is numeric
str(pheno$visit_month)
pheno$time_in_years <- pheno$visit_month / 12
# select covariates in pheno that are needed for mixed models
mixed_model_pheno <- pheno %>% select(sample_id, sex, status = DiseaseStatus, age_at_baseline, time_in_years, PATNO = participant_id)
mixed_model_pheno$PATNO <- str_replace_all(mixed_model_pheno$PATNO, "PD-", "")
mixed_model_pheno <- mixed_model_pheno[order(mixed_model_pheno$PATNO, mixed_model_pheno$time_in_years), ]
# Merge phenotype and cell count data
mixed_model_pheno <- merge(mixed_model_pheno, cc, by = "sample_id")
# Transpose the count matrix
mixed_model_count_matirx <- data.frame(t(data))
mixed_model_count_matirx$sample_id <- rownames(mixed_model_count_matirx)
# Merge phenotype and count matrix
mixed_model_longtable <- merge(mixed_model_pheno, mixed_model_count_matirx, by = "sample_id")
# create the long table
mixed_model_longtable <- melt(mixed_model_longtable, id.vars = c('sample_id', "PATNO", 'sex', 'age_at_baseline', 'status', 'time_in_years','B_cell',
                                                                 'T_cell_CD4','T_cell_CD8','Monocyte','NK_cell','Neutrophil'))
# Change names of variable and value columns
colnames(mixed_model_longtable)[colnames(mixed_model_longtable) == "variable"] <- "circRNA"
colnames(mixed_model_longtable)[colnames(mixed_model_longtable) == "value"] <- "counts"
str(mixed_model_longtable)
# Ensure that circRNA name is a character variable
mixed_model_longtable$circRNA <- as.character(mixed_model_longtable$circRNA)
# Order rows by participant ID
mixed_model_longtable <- mixed_model_longtable[order(mixed_model_longtable$PATNO), ]
# remove 0 counts
mixed_model_longtable<-mixed_model_longtable[mixed_model_longtable$counts != 0, ]
mixed_model_longtable <- na.omit(mixed_model_longtable)
# add baseline counts
mixed_model_longtable$counts_at_baseline <- -9
mixed_model_longtable <- mixed_model_longtable %>% group_by(PATNO, circRNA) %>% mutate(counts_at_baseline = counts[which(time_in_years == min(time_in_years))])
mixed_model_longtable %>% filter(PATNO == rand_sample, circRNA == rand_circ)
# Save the long table
write.table(mixed_model_longtable, 'PPMI.longTable.EU.allVisits.txt', quote= F, sep = '\t', row.names = F, col.names = T)

# II - Cell-count analyses --------------------------------------------------------------------
# 1. PDBP --------------------------------------------------------------------
# Load the data -----
# load the long table
data <- read.table("PDBP.longTable.EU.allVisits.txt", header = T, stringsAsFactors = F)
data <- data[order(data$PATNO, data$circRNA),]
# Load list of circRNAs identified in both PDBP
# and PPMI and subset data for only those circRNAs
com <- read.table("common_circRNAs_btw_PDBP_and_PPMI.txt")
com <- com$V1
data <- subset(data, circRNA %in% com)
# Format the data -----
# Generate a variable that captures the number of visits each participant presented at
data <- data %>% group_by(PATNO, circRNA) %>% mutate(visit_counts = length(time_in_years))
# filtered with as least two visits
atleast_2_visit <- subset(data, data$visit_counts > 1)
#Create a unique list of circRNAs identified in the dataset
loop_circs <- as.character(unique(atleast_2_visit$circRNA))
# Run mixed model analyses -----
full_result_table <- NULL
for (one_circ in loop_circs)
{
  if (which(loop_circs == one_circ) %% 100 == 0){print(paste("Progress:", one_circ, which(loop_circs == one_circ)))}
  tmp <- subset(atleast_2_visit, atleast_2_visit$circRNA == one_circ)
  model <- lmer(counts ~ counts_at_baseline + age_at_baseline + sex + status + time_in_years +
                  B_cell + T_cell_CD4 + T_cell_CD8+Monocyte + NK_cell + Neutrophil +
                  status * time_in_years + (1 | PATNO), data = tmp)
  res <- as.data.frame(coef(summary(model, ddf = "Satterthwaite")))
  colnames(res) <- c('beta', 'se', 'df', 't_value', 'p')
  rownames(res) <- c('intercept', 'count_at_bl', 'age_at_bl', 'sex', 'conditioncontrol', 'years',
                     'B_cell','T_cell_CD4','T_cell_CD8','Monocyte','NK_cell','Neutrophil','conditioncontrol:years')
  res$covar <- rownames(res)
  res$circ <- one_circ
  rownames(res) <- NULL
  full_result_table <- rbind(full_result_table, res)
}
# save results table
write.csv(full_result_table,"PDBP.DE.longitudinal.cellCounts.csv", row.names = F)

# 2. PPMI --------------------------------------------------------------------
# Load the data -----
# load the long table
data <- read.table("PPMI.longTable.EU.allVisits.txt", header = T, stringsAsFactors = F)
data <- data[order(data$PATNO, data$circRNA),]
# Load list of circRNAs identified in both PDBP
# and PPMI and subset data for only those circRNAs
com <- read.table("common_circRNAs_btw_PDBP_and_PPMI.txt")
com <- com$V1
data <- subset(data, circRNA %in% com)
# Format the data -----
# Generate a variable that captures the number of visits each participant presented at
data <- data %>% group_by(PATNO, circRNA) %>% mutate(visit_counts = length(time_in_years))
# filtered with as least two visits
atleast_2_visit <- subset(data, data$visit_counts > 1)
#Create a unique list of circRNAs identified in the dataset
loop_circs <- as.character(unique(atleast_2_visit$circRNA))
# Run mixed model analyses -----
full_result_table <- NULL
for (one_circ in loop_circs)
{
  if (which(loop_circs == one_circ) %% 100 == 0){print(paste("Progress:", one_circ, which(loop_circs == one_circ)))}
  tmp <- subset(atleast_2_visit, atleast_2_visit$circRNA == one_circ)
  model <- lmer(counts ~ counts_at_baseline + age_at_baseline + sex + status + time_in_years +
                  B_cell + T_cell_CD4 + T_cell_CD8+Monocyte + NK_cell + Neutrophil +
                  status * time_in_years + (1 | PATNO), data = tmp)
  res <- as.data.frame(coef(summary(model, ddf = "Satterthwaite")))
  colnames(res) <- c('beta', 'se', 'df', 't_value', 'p')
  rownames(res) <- c('intercept', 'count_at_bl', 'age_at_bl', 'sex', 'conditioncontrol', 'years',
                     'B_cell','T_cell_CD4','T_cell_CD8','Monocyte','NK_cell','Neutrophil','conditioncontrol:years')
  res$covar <- rownames(res)
  res$circ <- one_circ
  rownames(res) <- NULL
  full_result_table <- rbind(full_result_table, res)
}
# save results table
write.csv(full_result_table,"PPMI.DE.longitudinal.cellCounts.csv", row.names = F)

# 3. Meta analyses -----------------------------------------------------------
# Load and format the data -----
# Load DE data for both data sets
pdbp <- read.csv("PDBP.DE.longitudinal.cellCounts.csv")
ppmi <- read.csv("PPMI.DE.longitudinal.cellCounts.csv")
# Merge beta and pvalue
# Subset intercept and interaction terms
names(pdbp)<-c("PDBPbeta","se","df","t_value","PDBPpval","covar","circ")
names(ppmi)<-c("PPMIbeta","se","df","t_value","PPMIpval","covar","circ")
pdbpin<-pdbp[pdbp$covar=="intercept",c(1,5:7)]
ppmiin<-ppmi[ppmi$covar=="intercept",c(1,5:7)]
pdbpit<-pdbp[pdbp$covar=="conditioncontrol:years",c(1,5:7)]
ppmiit<-ppmi[ppmi$covar=="conditioncontrol:years",c(1,5:7)]
# Meta-analyses -----
DE<-merge(pdbpin,ppmiin,by=c("circ","covar"))
DE<-DE[sign(DE$PDBPbeta)==sign(DE$PPMIbeta),]
fishcomb <- fishercomb(list(DE$PDBPpval,DE$PPMIpval), BHth = 0.05)
metaDEin<-cbind(DE,fishcomb$rawpval,fishcomb$adjpval)
names(metaDEin)[c(7,8)]<-c("metaPval","metaPadj")
metaDEin<-metaDEin[,-8]
metaDEin[,4]<-formatC(metaDEin[,4], format = "e", digits = 3)
metaDEin[,6]<-formatC(metaDEin[,6], format = "e", digits = 3)
metaDEin[,7]<-formatC(metaDEin[,7], format = "e", digits = 3)
metaDEin[,3]<-round(metaDEin[,3],3)
metaDEin[,5]<-round(metaDEin[,5],3)
DE<-merge(pdbpit,ppmiit,by=c("circ","covar"))
nonDEit<-DE[sign(DE$PDBPbeta)!=sign(DE$PPMIbeta),]
DE<-DE[sign(DE$PDBPbeta)==sign(DE$PPMIbeta),]
fishcomb <- fishercomb(list(DE$PDBPpval,DE$PPMIpval), BHth = 0.05)
metaDEit<-cbind(DE,fishcomb$rawpval,fishcomb$adjpval)
names(metaDEit)[c(7,8)]<-c("metaPval","metaPadj")
metaDEit<-metaDEit[,-8]
metaDEit[,4]<-formatC(metaDEit[,4], format = "e", digits = 3)
metaDEit[,6]<-formatC(metaDEit[,6], format = "e", digits = 3)
metaDEit[,7]<-formatC(metaDEit[,7], format = "e", digits = 3)
metaDEit[,3]<-round(metaDEit[,3],3)
metaDEit[,5]<-round(metaDEit[,5],3)
nonDEit[,4]<-formatC(nonDEit[,4], format = "e", digits = 3)
nonDEit[,6]<-formatC(nonDEit[,6], format = "e", digits = 3)
nonDEit[,3]<-round(nonDEit[,3],3)
nonDEit[,5]<-round(nonDEit[,5],3)
nonDEit$metaPval<-"-"
metaDEit<-rbind(metaDEit,nonDEit)
# Merge the two meta data.frames
metaDE<-rbind(metaDEin,metaDEit)
write.csv(metaDE,"longitudinal.meta.cellCounts.csv", row.names = F)

# III - Medication analyses -----------------------------------------------------------
# Load and format the data -----
# load the long table
data <- read.table("PPMI.longTable.EU.allVisits.txt", header = T, stringsAsFactors = F)
data <- data[order(data$PATNO, data$circRNA),]
# Load list of circRNAs identified in both PDBP
# and PPMI and subset data for only those circRNAs
com <- read.table("common_circRNAs_btw_PDBP_and_PPMI.txt")
com <- com$V1
data <- subset(data, circRNA %in% com)
# Generate a variable that captures the number of visits each participant presented at
data <- data %>% group_by(PATNO, circRNA) %>% mutate(visit_counts = length(time_in_years))
# filtered with as least two visits
atleast_2_visit <- subset(data, data$visit_counts > 1)
atleast_2_visit<-atleast_2_visit[atleast_2_visit$status==1,]
#Create a unique list of circRNAs identified in the dataset
loop_circs <- as.character(unique(atleast_2_visit$circRNA))
# Run mixed model analyses -----
full_result_table <- NULL
for (one_circ in loop_circs)
{
  if (which(loop_circs == one_circ) %% 100 == 0){print(paste("Progress:", one_circ, which(loop_circs == one_circ)))}
  tmp <- subset(atleast_2_visit, atleast_2_visit$circRNA == one_circ)
  model <- lmer(counts ~ counts_at_baseline + age_at_baseline + sex + time_in_years +
                  PDMed + PDMed * time_in_years + (1 | PATNO), data = tmp)
  res <- as.data.frame(coef(summary(model, ddf = "Satterthwaite")))
  colnames(res) <- c('beta', 'se', 'df', 't_value', 'p')
  rownames(res) <- c('intercept', 'count_at_bl', 'age_at_bl', 'sex', 'years',
                     'PDMed','PDMed:years')
  res$covar <- rownames(res)
  res$circ <- one_circ
  rownames(res) <- NULL
  full_result_table <- rbind(full_result_table, res)
}
# save results table
write.csv(full_result_table,"PPMI.DE.longitudinal.PDMed.csv", row.names = F)

# IV - Mutation analyses -------------------------------------------------------
# Load and format the data ----
# load the long table
data <- read.table("PPMI.longTable.EU.allVisits.txt", header = T, stringsAsFactors = F)
data <- data[order(data$PATNO, data$circRNA),]
# Load list of circRNAs identified in both PDBP
# and PPMI and subset data for only those circRNAs
com <- read.table("common_circRNAs_btw_PDBP_and_PPMI.txt")
com <- com$V1
data <- subset(data, circRNA %in% com)
# Format mutation data
data$Mutation[is.na(data$Mutation)]<-"NONE"
# Generate a variable that captures the number of visits each participant presented at
data <- data %>% group_by(PATNO, circRNA) %>% mutate(visit_counts = length(time_in_years))
# filtered with as least two visits
atleast_2_visit <- subset(data, data$visit_counts > 1)
# 1. Mixed model analyses, idiopathic PD vs healthy control -----
dfdata<-atleast_2_visit[atleast_2_visit$Mutation=="NONE",]
# Create a unique list of circRNAs identified in the dataset
loop_circs <- as.character(unique(dfdata$circRNA))
#
full_result_table <- NULL
for (one_circ in loop_circs)
{
  if (which(loop_circs == one_circ) %% 100 == 0){print(paste("Progress:", one_circ, which(loop_circs == one_circ)))}
  tmp <- subset(dfdata, dfdata$circRNA == one_circ)
  model <- lmer(counts ~ counts_at_baseline + age_at_baseline + sex + status + time_in_years +
                  status * time_in_years + (1 | PATNO), data = tmp)
  res <- as.data.frame(coef(summary(model, ddf = "Satterthwaite")))
  colnames(res) <- c('beta', 'se', 'df', 't_value', 'p')
  rownames(res) <- c('intercept', 'count_at_bl', 'age_at_bl', 'sex', 'conditioncontrol', 'years',
                     'conditioncontrol:years')
  res$covar <- rownames(res)
  res$circ <- one_circ
  rownames(res) <- NULL
  full_result_table <- rbind(full_result_table, res)
}
# save results table
write.csv(full_result_table,"longitudinal.idPDvsCO.csv", row.names = F)
 
# 2. Mixed model analyses, LRRK2+ vs healthy control -----
dfdata<-atleast_2_visit[atleast_2_visit$status==0 | atleast_2_visit$Mutation=="LRRK2+",]
# Create a unique list of circRNAs identified in the dataset
loop_circs <- as.character(unique(dfdata$circRNA))
#
full_result_table <- NULL
for (one_circ in loop_circs)
{
  if (which(loop_circs == one_circ) %% 100 == 0){print(paste("Progress:", one_circ, which(loop_circs == one_circ)))}
  tmp <- subset(dfdata, dfdata$circRNA == one_circ)
  model <- lmer(counts ~ counts_at_baseline + age_at_baseline + sex + status + time_in_years +
                  status * time_in_years + (1 | PATNO), data = tmp)
  res <- as.data.frame(coef(summary(model, ddf = "Satterthwaite")))
  colnames(res) <- c('beta', 'se', 'df', 't_value', 'p')
  rownames(res) <- c('intercept', 'count_at_bl', 'age_at_bl', 'sex', 'conditioncontrol', 'years',
                     'conditioncontrol:years')
  res$covar <- rownames(res)
  res$circ <- one_circ
  rownames(res) <- NULL
  full_result_table <- rbind(full_result_table, res)
}
# save results table
write.csv(full_result_table,"longitudinal.LRRK2vsCO.csv", row.names = F)
 
# 3. Mixed model analyses, GBA+ cases vs healthy control -----
dfdata<-atleast_2_visit[atleast_2_visit$status==0 | atleast_2_visit$Mutation=="GBA+",]
# Create a unique list of circRNAs identified in the dataset
loop_circs <- as.character(unique(dfdata$circRNA))
#
full_result_table <- NULL
for (one_circ in loop_circs)
{
  if (which(loop_circs == one_circ) %% 100 == 0){print(paste("Progress:", one_circ, which(loop_circs == one_circ)))}
  tmp <- subset(dfdata, dfdata$circRNA == one_circ)
  model <- lmer(counts ~ counts_at_baseline + age_at_baseline + sex + status + time_in_years +
                  status * time_in_years + (1 | PATNO), data = tmp)
  res <- as.data.frame(coef(summary(model, ddf = "Satterthwaite")))
  colnames(res) <- c('beta', 'se', 'df', 't_value', 'p')
  rownames(res) <- c('intercept', 'count_at_bl', 'age_at_bl', 'sex', 'conditioncontrol', 'years',
                     'conditioncontrol:years')
  res$covar <- rownames(res)
  res$circ <- one_circ
  rownames(res) <- NULL
  full_result_table <- rbind(full_result_table, res)
}
# save results table
write.csv(full_result_table,"longitudinal.GBAvsCO.csv", row.names = F)

# V - Ancestry analyses ------------------------------------------------------
# Load and format the data -----
# load raw counts for PDBP and PPMI
pdbp<-read.csv("PDBP.all.circRNA.counts.csv", row.names = 1)
ppmi<-read.csv("PPMI.all.circRNA.counts.csv", row.names = 1)
# Load phenotype data
phenopdbp <- read_excel('full.pheno.txt',"PDBP")
phenoppmi <- read_excel('full.pheno.txt',"PPMI")
# Subset participants identified as Black or African American
phenopdbp<-phenopdbp[phenopdbp$race=="Black or African American",]
phenoppmi<-phenoppmi[phenoppmi$race=="Black or African American",]
# Subset phenotype for samples present in the PDBP
# count matrix and vice versa
pdbp<-pdbp[,colnames(pdbp)%in%phenopdbp$fileNames]
ppmi<-ppmi[,colnames(ppmi)%in%phenoppmi$fileNames]
# Format data for DE analyses
pheno<-rbind(phenopdbp[,c(3,17,4,6:8)],phenoppmi[,c(3,25,5,11,15,16)])
pheno$`Pt-Category`[pheno$`Pt-Category`=="Case"]<-1
pheno$`Pt-Category`[pheno$`Pt-Category`=="Control"]<-0
counts<-merge(pdbp, ppmi, by='row.names')
row.names(counts)<-counts[,1]
counts<-counts[,-1]
# Normalize the counts
counts <- counts[,pheno$fileNames]
counts <- as.matrix(counts)
names(pheno)<-c("participantID","fileName","visitMonth","status","age","sex")
dds <- DESeqDataSetFromMatrix(countData = counts, colData = pheno, design = ~ 1)
de <- varianceStabilizingTransformation(dds, blind = F)
vst_normalized_counts <- as.data.frame(de@assays@data@listData)
# Save normalized count matrix and pheno
write.csv(vst_normalized_counts,"AfA.normalized.counts.csv", row.names = T)
# Make the long table -----
# Calculate visit time in years instead of month
pheno$time_in_years <- pheno$visitMonth / 12
pheno<-pheno[,colnames(pheno)!='visitMonth']
# Transpose the count matrix
mixed_model_count_matirx <- data.frame(t(data))
mixed_model_count_matirx$fileName <- rownames(mixed_model_count_matirx)
# Merge phenotype and count matrix
mixed_model_longtable <- merge(pheno, mixed_model_count_matirx, by = "fileName")
# create the long table
mixed_model_longtable <- melt(mixed_model_longtable, id.vars = c('fileName', "participantID", 'sex', 'age', 'status', 'time_in_years'))
# Change names of variable and value columns
colnames(mixed_model_longtable)[colnames(mixed_model_longtable) == "variable"] <- "circRNA"
colnames(mixed_model_longtable)[colnames(mixed_model_longtable) == "value"] <- "counts"
# Ensure that circRNA name is a character variable
mixed_model_longtable$circRNA <- as.character(mixed_model_longtable$circRNA)
# Order rows by participant ID
mixed_model_longtable <- mixed_model_longtable[order(mixed_model_longtable$participantID), ]
# remove 0 counts
mixed_model_longtable<-mixed_model_longtable[mixed_model_longtable$counts != 0, ]
# add baseline counts
mixed_model_longtable$counts_at_baseline <- -9
mixed_model_longtable <- mixed_model_longtable %>% group_by(participantID, circRNA) %>% mutate(counts_at_baseline = counts[which(time_in_years == min(time_in_years))])
any(mixed_model_longtable$counts_at_baseline == -9) # should be FALSE
# Save the long table
write.table(mixed_model_longtable, 'AfA.longTable.allVisits.txt', quote= F, sep = '\t', row.names = F, col.names = T)
# Run mixed model -----
data <- mixed_model_longtable
data <- data[order(data$participantID, data$circRNA),]
# Load list of circRNAs identified in both PDBP
# and PPMI and subset data for only those circRNAs
com <- read.table("common_circRNAs_btw_PDBP_and_PPMI.txt")
com <- com$V1
data <- subset(data, circRNA %in% com)
# Generate a variable that captures the number of visits each participant presented at
data <- data %>% group_by(participantID, circRNA) %>% mutate(visit_counts = length(time_in_years))
# filtered with as least two visits
atleast_2_visit <- subset(data, data$visit_counts > 1)
#Create a unique list of circRNAs identified in the dataset
loop_circs <- as.character(unique(atleast_2_visit$circRNA))
# make participantID a numeric variable
atleast_2_visit$participantID<-gsub('[A-Z]','',atleast_2_visit$participantID)
atleast_2_visit$participantID<-gsub('-','',atleast_2_visit$participantID)
#
full_result_table <- NULL
for (one_circ in loop_circs)
{
  if (which(loop_circs == one_circ) %% 100 == 0){print(paste("Progress:", one_circ, which(loop_circs == one_circ)))}
  tmp <- subset(atleast_2_visit, atleast_2_visit$circRNA == one_circ)
  model <- lmer(counts ~ counts_at_baseline + age + sex + status + time_in_years +
                  status * time_in_years + (1 | participantID), data = tmp)
  res <- as.data.frame(coef(summary(model, ddf = "Satterthwaite")))
  colnames(res) <- c('beta', 'se', 'df', 't_value', 'p')
  rownames(res) <- c('intercept', 'count_at_bl', 'age_at_bl', 'sex', 'conditioncontrol', 'years',
                     'conditioncontrol:years')
  res$covar <- rownames(res)
  res$circ <- one_circ
  rownames(res) <- NULL
  full_result_table <- rbind(full_result_table, res)
}
# save results table
write.csv(full_result_table,"Afa.longitudinal.results.csv", row.names = F)
