# Load the necessary libraries -----
library(dplyr)
library(DESeq2)
library(stringr)
library(ggplot2)
# Load data -----
# Load raw circular RNA counts
visit1 <- read.delim("visit1.counts.tsv", check.names = F, sep = '\t')
visit2 <- read.delim("visit2.counts.tsv", check.names = F, sep = '\t')
visit3 <- read.delim("visit3.counts.tsv", check.names = F, sep = '\t')
visit4 <- read.delim("visit4.counts.tsv", check.names = F, sep = '\t')
visit5 <- read.delim("visit5.counts.tsv", check.names = F, sep = '\t')
# Load phenotypes
pheno <- read.delim("phenotype.csv", check.names = F, sep = ",")
# Filter and format data -----
# find common circular RNAs across all visits
circs <- list(bl = rownames(visit1),v1 = rownames(visit2),v2 = rownames(visit3),v3 = rownames(visit4),v4 = rownames(visit5))
circs <- Reduce(intersect, circs)
visit1 <- visit1[circs, ]
visit2 <- visit2[circs, ]
visit3 <- visit3[circs, ]
visit4 <- visit4[circs, ]
visit5 <- visit5[circs, ]
# merge count data into one
allCounts <- merge(visit1, visit2, by = "row.names", all = T)
rownames(allCounts) <- allCounts$Row.names
allCounts <- allCounts[, -1]
allCounts <- merge(allCounts, visit3, by = "row.names", all = T)
rownames(allCounts) <- allCounts$Row.names
allCounts <- allCounts[, -1]
allCounts <- merge(allCounts, visit4, by = "row.names", all = T)
rownames(allCounts) <- allCounts$Row.names
allCounts <- allCounts[, -1]
allCounts <- merge(allCounts, visit5, by = "row.names", all = T)
rownames(allCounts) <- allCounts$Row.names
allCounts <- allCounts[, -1]
colnames(allCounts) <- str_replace_all(colnames(allCounts), "T1.*.unified.Chimeric.out.junction", "T1")
# Ensure that phenotype is available for each sample 
# in the count matrix and vice versa
pheno <- subset(pheno, pheno$sample_id %in% colnames(allCounts))
allCounts <- select(allCounts, colnames(allCounts) %in% pheno$sample_id)
# Code healthy controls as 0 and PD cases as 1
pheno[pheno$DiseaseStatus == "Control", "DiseaseStatus"] <- 0
pheno[pheno$DiseaseStatus == "Case", "DiseaseStatus"] <- 1
# Code males as M and females as F
pheno[pheno$sex == "Male", "sex"] <- "M"
pheno[pheno$sex == "Female", "sex"] <- "F"
counts <- allCounts[, order(colnames(allCounts))]
counts <- as.matrix(counts)
# set up categorical variables with levels
pheno$DiseaseStatus <- factor(pheno$DiseaseStatus, levels = c(0, 1))
pheno$sex <- factor(pheno$sex, levels = c("M", "F"))
# PCA -----
dds <- DESeqDataSetFromMatrix(countData = counts, colData = pheno, design = ~ DiseaseStatus)
sizeFactors(dds) <- 1
dds.de <- estimateDispersions(dds)
vst.all <- varianceStabilizingTransformation(dds.de, blind=FALSE)
pcaData <- plotPCA(vst.all, intgroup="DiseaseStatus", returnData = T)
tmp <- pheno %>% select(name = sample_id, visit_month)
pcaData <- merge(pcaData, tmp, by = "name")
pcaData$visit_month <- as.factor(pcaData$visit_month)
# Calculate means and standard deviations for PC1 and PC2 
pc1.sd <- sd(pcaData$PC1)
pc1.mean <- mean(pcaData$PC1)
pc2.sd <- sd(pcaData$PC2)
pc2.mean <- mean(pcaData$PC2)
dist = 3 
range.pc1 <- c(pc1.mean - dist*pc1.sd, pc1.mean + dist*pc1.sd)
range.pc2 <- c(pc2.mean - dist*pc2.sd, pc2.mean + dist*pc2.sd)
# check PCA plots with SD borders
panel <- c("#E69F00", "#56B4E9", "#999999", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pcaData$group <- NA
pcaData[pcaData$DiseaseStatus == 0, "group"] <- "Control"
pcaData[pcaData$DiseaseStatus == 1, "group"] <- "PD_case"
ggplot(pcaData, aes(x = PC1, y = PC2, color = group)) + geom_point() +
  geom_hline(yintercept = range.pc2, linetype = 'dashed') +
  geom_vline(xintercept = range.pc1, linetype = 'dashed') +
  scale_color_manual(values=panel) + theme_bw()
# Remove samples that are more than 3 standard deviations
# away from the mean of either PC1 or PC2 from tboth the 
# count matrix and phenotype
clean.pcaData <- pcaData[pcaData$PC1 > range.pc1[1] & pcaData$PC1 < range.pc1[2], ]
clean.pcaData <- clean.pcaData[clean.pcaData$PC2 > range.pc2[1] & clean.pcaData$PC2 < range.pc2[2], ]
clean.pheno <- clean.pcaData$name
wr.clean.pheno <- pheno[pheno$sample_id %in% clean.pheno, ]
write.table(wr.clean.pheno, 'clean.pheno.txt', row.names = F, col.names = T)
count_matrix <- allCounts[, colnames(allCounts) %in% clean.pheno]
write.table(count_matrix, 'clean.counts.tsv', row.names = T, col.names = T, sep = '\t')
# Normalize counts -----
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = wr.clean.pheno, design = ~ DiseaseStatus)
sf <- estimateSizeFactors(dds)
sizeFactor <- sizeFactors(sf)
sizeFactors(dds) <- sizeFactor
dds.de <- DESeq(dds, parallel = F, betaPrior = FALSE)
de <- varianceStabilizingTransformation(dds.de, blind = F)
vst_normalized_counts <- as.data.frame(de@assays@data@listData)
write.table(vst_normalized_counts, 'normalized.counts.txt', col.names = T, row.names = T, sep = "\t")