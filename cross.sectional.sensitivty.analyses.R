# Load the necessary libraries -----
library(dplyr)
library(stringr)
library(flexiblas)
library(DESeq2)
library(metaRNASeq)
library(readxl)
library(flexiblas)
library(biomaRt)
library(corrplot)
library(Hmisc)
flexiblas_load_backend("__FALLBACK__")
# Define a "not in" opperator -----
`%!in%` <- Negate(`%in%`)
# I - Cell-count analyses  --------------------------------------------------------------------
# 0 for controls, 1 for cases
# 1.PDBP -----
# Load the data -----
# load raw count matrix (not the normalized one, DESeq2 does not take normalized value)
data <- read.table("PDBP/clean.counts.tsv", header = T, stringsAsFactors = F)
com <- read.table("common_circRNAs_btw_PDBP_and_PPMI.txt")
com <- com$V1
data <- subset(data, rownames(data) %in% com)
# load cleaned phenotype data
pheno <- read.table('PDBP/clean.pheno.txt', header = T, stringsAsFactors = F)
# load cell count data
cc <- readRDS("cell.counts.EPIC.rds")
# Format data --------------------------------------------------------------------
cc <- t(cc[-7,])
cnames<-cc[1,]
cc <- cc[-1,]
colnames(cc) <- gsub(" ","_",cnames)
rownames(cc) <- gsub("-",".",rownames(cc))
cc <- as.data.frame(cc)
cc$sample_id <- rownames(cc)
names(cc)[2] <- substr(names(cc)[2], 1, nchar(names(cc)[2])-1)
names(cc)[3] <- substr(names(cc)[3], 1, nchar(names(cc)[3])-1)
# Ensure that time variable is numeric
str(pheno$visit_month)
pheno$time_in_years <- pheno$visit_month / 12
pheno$age <- pheno$age_at_baseline + pheno$time_in_years
# Normalize cell counts using inverse normalization transfomation
inormal <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
cc<-cc[cc$sample_id%in%pheno$sample_id,]
rnames<-rownames(cc)
cc<-as.data.frame(sapply(cc[,-7], as.numeric))
cc<-as.data.frame(sapply(cc,inormal))
rownames(cc)<-rnames
cc$sample_id <- rownames(cc)
write.csv(cc,"normalized.cellcounts.PDBP.csv",row.names = F)
# select for last visit only
pheno <- pheno %>% group_by(participant_id) %>% mutate(last_visit = (visit_month == max(visit_month)))
pheno <- subset(pheno, pheno$last_visit == T)
cc<-cc[cc$sample_id%in%pheno$sample_id,]
# select covariates in pheno that are needed for DESeq2
deseq2_pheno <- pheno %>% select(sample_id, sex, status = Pt.Category, age, PATNO = participant_id)
deseq2_pheno$PATNO <- str_replace_all(deseq2_pheno$PATNO, "PD-", "")
# combine phenotype and cell count data
deseq2_pheno <- merge(deseq2_pheno, cc, by = "sample_id")
deseq2_count <- data[,deseq2_pheno$sample_id]
deseq2_count_input <- as.matrix(deseq2_count)
# Run DEA -----
dds <- DESeqDataSetFromMatrix(countData = deseq2_count_input, colData = deseq2_pheno,
                              design = ~ factor(sex) + age + B_cell + T_cell_CD4 +
                                T_cell_CD8+Monocyte + NK_cell + Neutrophil + factor(status))
dds.de <- DESeq(dds, parallel = F, betaPrior = FALSE)
dds.res <- results(dds.de, tidy = T)
#Save the results
write.csv(DE,"PDBP.DE.age.sex.cell.counts.status.csv", row.names=F)
# 2. PPMI --------------------------------------------------------------------
# Load the data -----
# load raw count matrix (not the normalized one, DESeq2 does not take normalized value)
data <- read.table("PPMI/clean.counts.tsv", header = T, stringsAsFactors = F)
com <- read.table("common_circRNAs_btw_PDBP_and_PPMI.txt")
com <- com$V1
data <- subset(data, rownames(data) %in% com)
# load cleaned phenotype data
pheno <- read.table('PPMI/clean.pheno.txt', header = T, stringsAsFactors = F)
# load cell count data
cc <- readRDS("cell.counts.EPIC.rds")
# Format the data -----
pheno <- subset(pheno, pheno$DeseaseStatus %in% c(0, 1))
pheno$sample_id<-paste("PP-",pheno$UID,"T1", sep="")
pheno$sample_id<-gsub("BL","BLM0",pheno$sample_id)
pheno$sample_id<-gsub("V02","SVM6",pheno$sample_id)
pheno$sample_id<-gsub("V04","SVM12",pheno$sample_id)
pheno$sample_id<-gsub("V06","SVM24",pheno$sample_id)
pheno$sample_id<-gsub("V08","SVM36",pheno$sample_id)
cc <- t(cc[-7,])
names<-cc[1,]
cc <- cc[-1,]
colnames(cc) <- gsub(" ","_",cnames)
cc <- as.data.frame(cc)
cc$sample_id <- rownames(cc)
names(cc)[2] <- substr(names(cc)[2], 1, nchar(names(cc)[2])-1)
names(cc)[3] <- substr(names(cc)[3], 1, nchar(names(cc)[3])-1)
# Normalize cell counts using inverse normalization transfomation
inormal <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
cc<-cc[cc$sample_id%in%pheno$sample_id,]
rnames<-rownames(cc)
cc<-as.data.frame(sapply(cc[,-7], as.numeric))
cc<-as.data.frame(sapply(cc,inormal))
rownames(cc)<-rnames
write.csv(cc,"normalized.cellcounts.PPMI.csv",row.names = F)
# select for last visit only
pheno <- pheno %>% group_by(PATNO) %>% mutate(last_visit = (Time == max(Time)))
pheno <- subset(pheno, pheno$last_visit == T)
cc<-cc[cc$sample_id%in%pheno$sample_id,]
# select covariates in pheno that are needed for DESeq2
deseq2_pheno <- pheno %>% select(FILE_NAME, sex = Gender, status = DeseaseStatus, age = AgeVisit, PATNO, sample_id)
# combine phenotype and cell count data
deseq2_pheno <- merge(deseq2_pheno, cc, by = "sample_id")

deseq2_count <- data[,deseq2_pheno$FILE_NAME]
deseq2_count_input <- as.matrix(deseq2_count)
# Run DEA ----- 
dds <- DESeqDataSetFromMatrix(countData = deseq2_count_input, colData = deseq2_pheno, design = ~ factor(sex) + age + B_cell + T_cell_CD4 +
                                T_cell_CD8 + Monocyte + NK_cell + Neutrophil + factor(status))
dds.de <- DESeq(dds, parallel = F, betaPrior = FALSE)
dds.res <- results(dds.de, tidy = T)
# Save the results
write.csv(DE,"PPMI.DE.age.sex.cell.counts.status.csv", row.names = F)
# 3. Meta analyses -----------------------------------------------------------
# Load DE data for both data sets
pdbp <- read.csv("PDBP.DE.age.sex.cell.counts.status.csv")
ppmi <- read.csv("PPMI.DE.age.sex.cell.counts.status.csv")
# Merge log2FC, pvalue and padj data
pdbp<-pdbp[,c(1,3,6)]
ppmi<-ppmi[,c(1,3,6)]
names(pdbp)<-c("circRNA","disFC","disPraw")
names(ppmi)<-c("circRNA","repFC","repPraw")
df<-merge(pdbp, ppmi, by="circRNA")
# Retain only circRNAs where log2FC has the same sign
nonDE<-DE[sign(DE$disFC)!=sign(DE$repFC),]
DE<-DE[sign(DE$disFC)==sign(DE$repFC),]
# Perform meta analyses
fishcomb <- fishercomb(list(DE$disPraw,DE$repPraw), BHth = 0.05)
metaDE<-cbind(df,fishcomb$rawpval,fishcomb$adjpval)
names(metaDE)[c(6,7)]<-c("metaPraw","metaPadj")
DEcircs<-metaDE$Gene
nonDEcircs<-nonDE$Gene
metaDE<-as.data.frame(sapply(metaDE[-c(1,7)],function(a) formatC(a, format = "e", digits = 3)))
nonDE<-as.data.frame(sapply(nonDE[-1],function(a) formatC(a, format = "e", digits = 3)))
metaDE$circRNA<-DEcircs
nonDE$circRNA<-nonDEcircs
nonDE$metaPvalue<-"-"
metaDE<-rbind(metaDE,nonDE)
write.csv(metaDE,"DE.meta.analyses.age.sex.cell.counts.status.csv", row.names = F)
# II - Medication analyses (PPMI)  --------------------------------------------------------------------
# Load the data -----
# load raw count matrix 
data <- read.table("PPMI/clean.counts.tsv", header = T, stringsAsFactors = F)
# load cleaned phenotype data
pheno <- read.table('clean.pheno.txt', header = T, stringsAsFactors = F)
pheno <- subset(pheno, pheno$DiseaseStatus ==1)
# Format data -----
# select for last visit only
pheno <- pheno %>% group_by(PATNO) %>% mutate(last_visit = (Time == max(Time)))
pheno <- subset(pheno, pheno$last_visit == T)
# select covariates in pheno that are needed for DESeq2
deseq2_pheno <- pheno %>% select(FILE_NAME, sex = Gender, status = DiseaseStatus, age = AgeVisit, PATNO, PDMed)
# Format DESeq2 input
deseq2_count <- data[,deseq2_pheno$FILE_NAME]
deseq2_count_input <- as.matrix(deseq2_count)
# Run DEA -----
dds <- DESeqDataSetFromMatrix(countData = deseq2_count_input, colData = deseq2_pheno, design = ~ factor(sex) + age
                              + factor(PDMed))
dds.de <- DESeq(dds, parallel = F, betaPrior = FALSE)
dds.res <- results(dds.de, tidy = T)
# save the DE results
write.csv(DE,"PPMI.DE.CAonly.age.sex.PDmeds.csv", row.names = F)
# III - Ancestry analyses ------------------------------------------------------
# Load the data -----
# Load raw counts for all PDBP visits and combine them into 1 matrix
pdbp <- read.table("PDBP/clean.counts.tsv", header = T, stringsAsFactors = F)
ppmi <- read.table("PPMI/clean.counts.tsv", header = T, stringsAsFactors = F)
# # Load phenotype data
phenopdbp <- read_excel('full.pheno.txt',"PDBP")
phenoppmi <- read_excel('full.pheno.txt',"PPMI")
# Subset participants identified as Black or African American
phenopdbp<-phenopdbp[phenopdbp$race=="Black or African American",]
phenoppmi<-phenoppmi[phenoppmi$race=="Black or African American",]
# Subset phenotype for samples present in the
# count matrix and vice versa
pdbp<-pdbp[,colnames(pdbp)%in%phenopdbp$fileNames]
ppmi<-ppmi[,colnames(ppmi)%in%phenoppmi$fileNames]
# Format data -----
pheno<-rbind(phenopdbp,phenoppmi)
pheno<-pheno[pheno$`Pt-Category`!="Prodromal",]
pheno$`Pt-Category`[pheno$`Pt-Category`=="Case"]<-1
pheno$`Pt-Category`[pheno$`Pt-Category`=="Control"]<-0
counts<-merge(pdbp, ppmi, by='row.names')
row.names(counts)<-counts[,1]
counts<-counts[,-1]
# Select only the last visit for each participant
pheno <- pheno %>% group_by(participant_id) %>% mutate(last_visit = (visit_month == max(visit_month)))
pheno <- subset(pheno, pheno$last_visit == T)
counts<-counts[,colnames(counts)%in%pheno$fileNames]
# Run DEA -----
counts <- counts[,pheno$fileNames]
counts <- as.matrix(counts)
names(pheno)<-c("participantID","fileName","visitMonth","status","age","sex","lastVisit")
dds <- DESeqDataSetFromMatrix(countData = counts, colData = pheno, design = ~ factor(sex) + age + factor(status))
dds.de <- DESeq(dds, parallel = F, betaPrior = FALSE)
dds.res <- results(dds.de, tidy = T)
# Save DE results
write.csv(DE,"AfA.DE.age.sex.status.csv", row.names = F)

# IV - Mutation analyses (PPMI; Ca/Co) -------------------------------------------------------
# Load the data -----
# load raw count matrix (not the normalized one, DESeq2 does not take normalized value)
data <- read.table("PPMI/clean.counts.tsv", header = T, stringsAsFactors = F)
# Load a list of circRNAs commeon for PDBP and PPMI
com <- read.table("common_circRNAs_btw_PDBP_and_PPMI.txt")
com <- com$V1
data <- subset(data, rownames(data) %in% com)
# load cleaned phenotype data
pheno <- read.table('clean.pheno.txt', header = T, stringsAsFactors = F)
# Format data
pheno <- subset(pheno, pheno$DeseaseStatus %in% c(0, 1))
pheno$Mutation[is.na(pheno$Mutation)]<-"NONE"
# select for last visit only
pheno <- pheno %>% group_by(PATNO) %>% mutate(last_visit = (Time == max(Time)))
pheno <- subset(pheno, pheno$last_visit == T)
# select covariates in pheno that are needed for DESeq2
pheno <- pheno %>% select(FILE_NAME, sex = Gender, status = DeseaseStatus, age = AgeVisit, PATNO, Mutation)
data <- data[,pheno$FILE_NAME]
# 1. DEA, idiopathic PD vs healthy control -----
depheno<-pheno[pheno$Mutation=="NONE",]
dedata<-data[,colnames(data)%in%depheno$FILE_NAME]
dedata<-dedata[,depheno$FILE_NAME]
dedata <- as.matrix(dedata)
dds <- DESeqDataSetFromMatrix(countData = dedata, colData = depheno, design = ~ factor(sex) + age + factor(status))
dds.de <- DESeq(dds, parallel = F, betaPrior = FALSE)
dds.res <- results(dds.de, tidy = T)
# Save DE results
write.csv(DE,"PPMI.DE.age.sex.idPDvsCO.csv", row.names = F)
# 2. DEA, LRRK2+ cases vs CO -----
depheno<-pheno[pheno$status==0 | pheno$Mutation=="LRRK2+",]
dedata<-data[,colnames(data)%in%depheno$FILE_NAME]
dedata<-dedata[,depheno$FILE_NAME]
dedata <- as.matrix(dedata)
dds <- DESeqDataSetFromMatrix(countData = dedata, colData = depheno, design = ~ factor(sex) + age + factor(status))
dds.de <- DESeq(dds, parallel = F, betaPrior = FALSE)
dds.res <- results(dds.de, tidy = T)
# Save DE results
write.csv(DE,"PPMI.DE.age.sex.LRRK2vsCO.csv", row.names = F)
# 3. DEA, GBA+ cases vs CO -----
depheno<-pheno[pheno$status==0 | pheno$Mutation=="GBA+",]
dedata<-data[,colnames(data)%in%depheno$FILE_NAME]
dedata<-dedata[,depheno$FILE_NAME]
dedata <- as.matrix(dedata)
dds <- DESeqDataSetFromMatrix(countData = dedata, colData = depheno, design = ~ factor(sex) + age + factor(status))
dds.de <- DESeq(dds, parallel = F, betaPrior = FALSE)
dds.res <- results(dds.de, tidy = T)
# Save DE results
write.csv(DE,"PPMI.DE.age.sex.GBAvsCO.csv", row.names = F)
# 4. DEA, SNCA+ cases vs CO -----
depheno<-pheno[pheno$status==0 | pheno$Mutation=="SNCA+",]
dedata<-data[,colnames(data)%in%depheno$FILE_NAME]
dedata<-dedata[,depheno$FILE_NAME]
dedata <- as.matrix(dedata)
dds <- DESeqDataSetFromMatrix(countData = dedata, colData = depheno, design = ~ factor(sex) + age + factor(status))
dds.de <- DESeq(dds, parallel = F, betaPrior = FALSE)
dds.res <- results(dds.de, tidy = T)
# Save DE results
write.csv(DE,"PPMI.DE.age.sex.SNCAvsCO.csv", row.names = F)
# 5. DEA, genetic PD vs CO -----
depheno<-pheno[pheno$status==0 | pheno$Mutation!="NONE",]
dedata<-data[,colnames(data)%in%depheno$FILE_NAME]
dedata<-dedata[,depheno$FILE_NAME]
dedata <- as.matrix(dedata)
dds <- DESeqDataSetFromMatrix(countData = dedata, colData = depheno, design = ~ factor(sex) + age + factor(status))
dds.de <- DESeq(dds, parallel = F, betaPrior = FALSE)
dds.res <- results(dds.de, tidy = T)
# Save DE results
write.csv(DE,"PPMI.DE.age.sex.genPDvsCO.csv", row.names = F)


# V - MoCA/UPDRSIII correlation --------------------------------------------------------------------
# 1. PPMI --------------------------------------------------------------------
# Load the data -----
# Load normalized counts
data <- read.table("PPMI/normalized.counts.txt", header = T, stringsAsFactors = F)
circs<-c("circAFF2","circCCDC91","circETFA","circFAM13B","circITGAX","circNCF1","circPADI4","circSPI1","circSUZ12")
data <- data %>% subset(rownames(data) %in% circs)
# load cleaned phenotype data
pheno <- read.table('PDBP/clean.pheno.txt', header = T, stringsAsFactors = F)
# Format data -----
# select last visit only
pheno <- pheno %>% group_by(PATNO) %>% mutate(last_visit = (Time == max(Time)))
pheno <- subset(pheno, pheno$last_visit == T)
pheno <- pheno %>% ungroup(PATNO)
# Ensure that phenotype is available for all samples
# in the cont matrix and vice versa
df <- data[, pheno$FILE_NAME]
df <- df %>% mutate(sample_id = rownames(df))
pheno.PPMI <- pheno
count.PPMI <- df %>% t() %>% as.data.frame()
moca_pheno <- pheno %>% select(sample_id = FILE_NAME, MoCA)
updrs_pheno <- pheno %>% select(sample_id = FILE_NAME, UPDRS)

tmp_moca <- merge(moca_pheno, df, by = "sample_id")
tmp_updrs <- merge(updrs_pheno, df, by = "sample_id")
rownames(tmp_moca) <- tmp_moca$sample_id
rownames(tmp_updrs) <- tmp_updrs$sample_id
tmp_moca <- tmp_moca %>% select(-sample_id) %>% na.omit()
tmp_updrs <- tmp_updrs %>% select(-sample_id) %>% na.omit()
# Run correlation analyses -----
x <- tmp_moca %>% as.matrix() %>% rcorr(type = "spearman")
y <- tmp_updrs %>% as.matrix() %>% rcorr(type = "spearman")
x1.p <- x$P %>% as.data.frame() %>%  select(PPMI.MoCA = MoCA) %>% subset(rownames(x$P) != "MoCA")
y1.p <- y$P %>% as.data.frame() %>%  select(PPMI.UPDRS = UPDRS) %>% subset(rownames(y$P) != "UPDRS")

tmp_moca %>% cor(method = "spearman") %>% corrplot()
tmp_updrs %>% cor(method = "spearman") %>% corrplot()

# 2. PDBP --------------------------------------------------------------------
# Load the data -----
# load the normalized counts
data <- read.table("PDBP/normalized.counts.txt", header = T, stringsAsFactors = F)
circs<-c("circAFF2","circCCDC91","circETFA","circFAM13B","circITGAX","circNCF1","circPADI4","circSPI1","circSUZ12")
data <- data %>% subset(rownames(data) %in% circs)
# load cleaned phenotype data
pheno <- read.table('PDBP/clean.pheno.txt', header = T, stringsAsFactors = F)
# Format data -----
# check time variable is numeric
str(pheno$visit_month)
pheno$time_in_years <- pheno$visit_month / 12
pheno$age <- pheno$age_at_baseline + pheno$time_in_years
# select  last visit only
pheno <- pheno %>% group_by(participant_id) %>% mutate(last_visit = (visit_month == max(visit_month)))
pheno <- subset(pheno, pheno$last_visit == T)
pheno <- pheno %>% ungroup(participant_id)
# Ensure that phenotype is available for all samples
# in the cont matrix and vice versa
df <- data[, pheno$sample_id]
df <- df %>% mutate(sample_id = rownames(df))

pheno.PDBP <- pheno
count.PDBP <- df %>% t() %>% as.data.frame()

moca_pheno <- pheno %>% select(sample_id, MoCA)
updrs_pheno <- pheno %>% select(sample_id, UPDRS = UPDRSIII)

tmp_moca <- merge(moca_pheno, df, by = "sample_id")
tmp_updrs <- merge(updrs_pheno, df, by = "sample_id")
rownames(tmp_moca) <- tmp_moca$sample_id
rownames(tmp_updrs) <- tmp_updrs$sample_id
tmp_moca <- tmp_moca %>% select(-sample_id) %>% subset(MoCA != ".") %>% na.omit()
tmp_updrs <- tmp_updrs %>% select(-sample_id) %>% subset(UPDRS != ".") %>% na.omit()
tmp_moca$MoCA <- as.numeric(tmp_moca$MoCA)
tmp_updrs$UPDRS <- as.numeric(tmp_updrs$UPDRS)

# Run correlation analyses -----
x <- tmp_moca %>% as.matrix() %>% rcorr(type = "spearman")
y <- tmp_updrs %>% as.matrix() %>% rcorr(type = "spearman")
x2.p <- x$P %>% as.data.frame() %>%  select(PDBP.MoCA = MoCA) %>% subset(rownames(x$P) != "MoCA")
y2.p <- y$P %>% as.data.frame() %>%  select(PDBP.UPDRS = UPDRS) %>% subset(rownames(y$P) != "UPDRS")
tmp_moca %>% cor(method = "spearman") %>% corrplot()
tmp_updrs %>% cor(method = "spearman") %>% corrplot()

res <- cbind(x1.p, y1.p, x2.p, y2.p)

# 3. PDBP+PPMI  --------------------------------------------------------------------
# Load and format PPMI data -----
# Load the raw counts matrix
data <- read.table("PPMI/clean.counts.tsv", header = T, stringsAsFactors = F)
# load clean phenotype data
pheno <- read.table('PPMI/clean.pheno.txt', header = T, stringsAsFactors = F)
pheno <- subset(pheno, pheno$Status_time %in% c(0, 1))
# select last visit only
pheno <- pheno %>% group_by(PATNO) %>% mutate(last_visit = (Time == max(Time)))
pheno <- subset(pheno, pheno$last_visit == T)
pheno <- pheno %>% ungroup(PATNO)
# Ensure that phenotype is available for all samples
# in the cont matrix and vice versa
df <- data[, pheno$FILE_NAME]
# Subset to circRNAs of interest
df <- df %>% subset(rownames(df) %in% circs) %>% t() %>% as.data.frame()
df <- df %>% mutate(sample_id = rownames(df))

pheno.PPMI <- pheno
count.PPMI <- df %>% t() %>% as.data.frame()
# Load and format PDBP data -----
data <- read.table("PDBP/clean.counts.tsv", header = T, stringsAsFactors = F)
# load clean phenotype data
pheno <- read.table('PDBP/clean.pheno.txt', header = T, stringsAsFactors = F)
# check time variable is numeric
str(pheno$visit_month)
pheno$time_in_years <- pheno$visit_month / 12
pheno$age <- pheno$age_at_baseline + pheno$time_in_years
# selectlast visit only
pheno <- pheno %>% group_by(participant_id) %>% mutate(last_visit = (visit_month == max(visit_month)))
pheno <- subset(pheno, pheno$last_visit == T)
pheno <- pheno %>% ungroup(participant_id)
# Ensure that phenotype is available for all samples
# in the cont matrix and vice versa
df <- data[, pheno$sample_id]
# Subset to circRNAs of interest
df <- df %>% subset(rownames(df) %in% circs) %>% t() %>% as.data.frame()
df <- df %>% mutate(sample_id = rownames(df))

pheno.PDBP <- pheno
count.PDBP <- df %>% t() %>% as.data.frame()

# Normalize PDBP and PPMI together -----
count.PDBP <- count.PDBP %>% subset(rownames(count.PDBP) != "sample_id")
count.PPMI <- count.PPMI %>% subset(rownames(count.PPMI) != "sample_id")
tmp.count <- cbind(count.PDBP, count.PPMI)
name <- rownames(tmp.count)
tmp.count <- sapply(tmp.count, as.numeric)
rownames(tmp.count) <- name

pheno.PDBP <- pheno.PDBP %>% select(sample_id, MoCA, UPDRS = UPDRSIII, status = Pt.Category)
pheno.PPMI <- pheno.PPMI %>% select(sample_id = FILE_NAME, MoCA, UPDRS, status = Status_time)
pheno.PDBP$cohort <- "PDBP"
pheno.PPMI$cohort <- "PPMI"
tmp.pheno <- rbind(pheno.PDBP, pheno.PPMI) %>% as.data.frame()
rownames(tmp.pheno) <- tmp.pheno$sample_id
tmp.pheno$status <- as.factor(tmp.pheno$status)

dds <- DESeqDataSetFromMatrix(countData = tmp.count, colData = tmp.pheno, design = ~ factor(status))
sizeFactors(dds) <- 1
dds.de <- estimateDispersions(dds)

vst.all <- varianceStabilizingTransformation(dds.de, blind=FALSE)
vst_normalized_counts <- as.data.frame(vst.all@assays@data@listData)
tmp.count <- vst_normalized_counts %>% t() %>% as.data.frame()
tmp.count$sample_id <- rownames(tmp.count)
tmp <- merge(tmp.count, tmp.pheno, by = "sample_id")
rownames(tmp) <- tmp$sample_id
tmp <- tmp %>% select(-status, -cohort, -sample_id)

tmp_moca <- tmp %>% select(-UPDRS) %>% subset(MoCA != ".") %>% na.omit()
tmp_updrs <- tmp %>% select(-MoCA) %>% subset(UPDRS != ".") %>% na.omit()
# Run correlation analyses -----
x <- tmp_moca %>% as.matrix() %>% rcorr(type = "spearman")
y <- tmp_updrs %>% as.matrix() %>% rcorr(type = "spearman")
x3.p <- x$P %>% as.data.frame() %>%  select(All.MoCA = MoCA) %>% subset(rownames(x$P) != "MoCA")
y3.p <- y$P %>% as.data.frame() %>%  select(All.UPDRS = UPDRS) %>% subset(rownames(y$P) != "UPDRS")

res <- cbind(x1.p, y1.p, x2.p, y2.p, x3.p, y3.p)
out <- res[!rownames(res) %in% c("circC1orf112", "circRAB11FIP1", "circUTRN"), ]
write.table(out, "MoCA_UPDRS_correlation.txt", quote = F, row.names = T, col.names = T, sep = "\t")

