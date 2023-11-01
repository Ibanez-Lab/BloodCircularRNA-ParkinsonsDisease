# 1 Load Libraries & Import Data ---------------------------------------------------------
library(data.table)
library(tximport)
library(DESeq2)
library(matrixStats)
library(dplyr)
library(ggplot2)
library(reshape2)
library(Hmisc)
library(stringr)
library(ggpubr)
library(ggrepel)
library(splitstackshape)
library(org.Hs.eg.db)
library(flexiblas)
#flexiblas_switch(3)

# * 1.1 Import pheno  ------------------------

#load phenotype 
pheno <- read.table('pheno.txt')

# Load cell count data 
cellCounts <- read.csv("cellcounts.csv")

# Merge pheno & cell counts
pheno <- merge(pheno, cellCounts, by.x = "sample_id")

# check time variable is numeric
str(pheno$visit_month)
pheno$time_in_years <- pheno$visit_month / 12
pheno$age <- pheno$age_at_baseline + pheno$time_in_years

# select for last visit only in pheno
pheno <- pheno %>% group_by(participant_id) %>% mutate(last_visit = (visit_month == max(visit_month)))
lastVisitSubset <- which(pheno$last_visit == T)
pheno <- subset(pheno, pheno$last_visit == T)

# * 1.2 Generate Count Matrix  ------------------------

# Base directory where subdirectories with quant.genes.sf files are located
base_directory <- "/path/"

# List all subdirectories
subdirectories <- list.dirs(base_directory, full.names = TRUE, recursive = TRUE)

# Initialize empty list to store data frames and row names
data_frames <- list()
row_names <- NULL

# Loop through each subdirectory
for (subdir in subdirectories) {
  quant_file <- file.path(subdir, "quant.genes.sf")
  
  if (file.exists(quant_file)) {
    individual_name <- basename(subdir)   # Get the name of the individual
    
    # Read the quant.genes.sf file
    quant_data <- fread(quant_file)
    
    # Extract the NumReads column and add it to the data frame list
    data_frames[[individual_name]] <- quant_data$NumReads
    
  }
}

# Combine data frames into a single matrix
countMatrix <- do.call(cbind, data_frames)

#* 1.3 Clean & Manipulate Data for QC ------------------------

# Set column names to individual names without .salmon suffix and change - to .
colnames(countMatrix) <- sub("\\.salmon$", "", names(data_frames))
colnames(countMatrix) <- gsub("-", ".", colnames(countMatrix))

# Set row names to Names column from quant_data
row.names(countMatrix)<-quant_data$Name

#Removing sample with low reads
countMatrixclean <- countMatrix[(rowCounts(countMatrix[,-1]<10) < round(0.9*dim(countMatrix[,-1])[2])),]

pheno2 <- pheno[pheno$sample_id %in% colnames(countMatrixclean),]

# Extract the sample IDs from the pheno data
sample_ids <- pheno2$sample_id

# Subset the countMatrix to include only the columns in sample_ids
countMatrixfinal<- countMatrixclean[, sample_ids]
countMatrixfinal<-floor(countMatrixfinal)
#save the count matrix for easier loading

row.names(pheno)<-pheno$sample_id
countMatrix<-count_matrix[,colnames(count_matrix)%in%rownames(pheno)]
pheno<-pheno[rownames(pheno)%in%colnames(countMatrix),]

# 2 Quality Check ---------------------------------------------------------
# * 2.1 DESeq Model, VST ----
rnaDDS <- DESeqDataSetFromMatrix(countData = countMatrix, colData = pheno, design = ~ 1)

vstDDS<-vst(rnaDDS)

# * 2.2 Principal Components Analyses ----
# Case vs Control
pcaData <- plotPCA(vstDDS, intgroup="Pt.Category", returnData = T)
pc1.sd <- sd(pcaData$PC1)
pc1.mean <- mean(pcaData$PC1)
pc2.sd <- sd(pcaData$PC2)
pc2.mean <- mean(pcaData$PC2)
dist = 3 # choose of 3 sd or 2 sd
range.pc1 <- c(pc1.mean - dist*pc1.sd, pc1.mean + dist*pc1.sd)
range.pc2 <- c(pc2.mean - dist*pc2.sd, pc2.mean + dist*pc2.sd)
# check PCA plots with SD borders
panel <- c("#E69F00", "#56B4E9", "#999999", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

pcaData <- pcaData %>% mutate(group = NULL)
pcaData[pcaData$Pt.Category == 0, "group"] <- "Control"
pcaData[pcaData$Pt.Category == 1, "group"] <- "PD_Cases"

plotStatus <- ggplot(pcaData, aes(x = PC1, y = PC2, color = group)) + geom_point() +
  geom_hline(yintercept = range.pc2, linetype = 'dashed') +
  geom_vline(xintercept = range.pc1, linetype = 'dashed') +
  scale_color_manual(values=panel) + theme_bw()
# PCA by sex
rnaPCA<-plotPCA(vstDDS, intgroup = "sex", returnData = T)

myPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000", "#6C2EB8")

rnaPercentVar <- round(100 * attr(rnaPCA, "percentVar"))
Vstdiv <- c((mean(rnaPCA$PC1)-3*sd(rnaPCA$PC1)), (mean(rnaPCA$PC1)-2*sd(rnaPCA$PC1)), (mean(rnaPCA$PC1)-sd(rnaPCA$PC1)),mean(rnaPCA$PC1),
            (mean(rnaPCA$PC1)+sd(rnaPCA$PC1)),(mean(rnaPCA$PC1)+2*sd(rnaPCA$PC1)), (mean(rnaPCA$PC1)+3*sd(rnaPCA$PC1)))
vLineColor <- c("#E69F00", "#56B4E9", "#009E73", "#999999","#009E73", "#E69F00", "#56B4E9")
Hstdiv <- c((mean(rnaPCA$PC2)-3*sd(rnaPCA$PC2)), (mean(rnaPCA$PC2)-2*sd(rnaPCA$PC2)), (mean(rnaPCA$PC2)-sd(rnaPCA$PC2)),mean(rnaPCA$PC2),
            (mean(rnaPCA$PC2)+sd(rnaPCA$PC2)),  (mean(rnaPCA$PC2)+2*sd(rnaPCA$PC2)),  (mean(rnaPCA$PC2)+3*sd(rnaPCA$PC2)))
hLineColor <- c("#E69F00", "#56B4E9", "#009E73", "#999999","#009E73", "#E69F00", "#56B4E9")
plotStatus<-ggplot(rnaPCA, aes(x = PC1, y = PC2, color = sex)) +
  geom_point(size =5) +
  xlab(paste0("PC1: ", rnaPercentVar[1], "% variance")) +
  ylab(paste0("PC2: ", rnaPercentVar[2], "% variance")) +
  geom_vline(xintercept=Vstdiv, color = vLineColor, linetype="dotted") +
  geom_hline(yintercept=Hstdiv, color = hLineColor, linetype="dotted") +
  scale_color_manual(values = myPalette)+
  theme_bw()+
  theme(legend.position="top")

# Generate PCA with ALL genes using prcomp instead of plotPCA
normCntMatrix <- t(assay(vstDDS))
allPCAplot <-  prcomp(normCntMatrix, scale = T)
PCs <-  allPCAplot$x[,c(1:2)]
PC2Matrix <-  cbind(PCs, pheno)
plotStatus<-ggplot(PC2Matrix, aes(x = PC1, y = PC2, color = sex)) +
  geom_point(size =5) +
  xlab("PC1") +
  ylab("PC2") +
  geom_vline(xintercept=Vstdiv, color = vLineColor, linetype="dotted") +
  geom_hline(yintercept=Hstdiv, color = hLineColor, linetype="dotted") +
  scale_color_manual(values = myPalette)+
  theme_bw()+
  theme(legend.position="top")

# 3 Differential Expression Analysis -------------------------------------------

# * 3.1 Generic DESeq ------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = countMatrix, colData = pheno, design = ~ factor(sex) + age + factor(Pt.Category))
dds.de <- DESeq(dds, parallel = F, betaPrior = FALSE)

dds.res <- results(dds.de, tidy = T)

#Make a copy of dds.res
dds.res2 <- dds.res
#Remove the dot and proceeding numbers
dds.res2$row <- sub("\\..*","",dds.res2$row)
dds.res2$row <- mapIds(org.Hs.eg.db, keys = dds.res2$row, column = "SYMBOL", keytype = "ENSEMBL")

summary(dds.res)
dds.sig <- subset(dds.res2, dds.res2$padj < 0.05)

# plotting the volcano plot
DE <- dds.res2[order(dds.res2$padj), ]
volcano.name <- "PDBP Linear RNA Expression \n Last Visit ca/co"
l2f.lim <- 0
plim <- 0.05
choice <- 'padj'
plotDE <- DE %>% mutate(gene_type = case_when(log2FoldChange >= l2f.lim & get(choice) <= plim ~ "up",
                                              log2FoldChange <= -l2f.lim & get(choice) <= plim ~ "down", TRUE ~ "ns"))
plotDE <- plotDE[order(plotDE$pvalue), ]
#Get the 9 circ degs
result_9circ<- read.table("circRNAresults.txt") 
result_9circ$linear <- sub("\\..*","",result_9circ$linear)
result_9circ$linear <- mapIds(org.Hs.eg.db, keys = result_9circ$linear, column = "SYMBOL", keytype = "ENSEMBL")
circ_names<-result_9circ$linear

nine_circ<- plotDE[plotDE$row %in% circ_names, ]
cols <- densCols(plotDE$log2FoldChange, plotDE$padj)
cols[plotDE$gene_type=='up']<-"#E69F00"
cols[plotDE$gene_type=='down']<-"#56B4E9"
cols[plotDE$gene_type=='ns']<-"#999999"
sizes <- c("up" = 3, "down" = 3, "ns" = 2)
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)
sig <- plotDE[plotDE$pvalue<=plim & plotDE$gene_type!='ns',]

ggplot(plotDE, aes(x=log2FoldChange, y=-log10(pvalue), size=gene_type, alpha = gene_type))+
  geom_point(col = cols)+
  geom_hline(yintercept = -log10(plim), linetype = "dashed") +
  geom_vline(xintercept = c(-l2f.lim, l2f.lim), linetype = "dashed") +
  geom_label_repel(data = nine_circ, aes(label = row), force = 2, nudge_y = 1) +
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  ggtitle(volcano.name) +
  xlab("Effect size: log2(fold-change)") +
  ylab("-log10(p-value)")+
  theme_bw()+
  theme(legend.position = "none")


# * 3.2 Cell Counts DESeq ------------------------------------------------------

dds.cc <- DESeqDataSetFromMatrix(countData = countMatrix, colData = pheno, design = ~ factor(sex) + age + B_cell + T_cell_CD4 + T_cell_CD8 + Monocyte + NK_cell + Neutrophil + factor(Pt.Category))
dds.cc.de <- DESeq(dds.cc, parallel = F, betaPrior = FALSE)
dds.cc.res <- results(dds.cc.de, tidy = T)
# dds.res <- as.data.frame(dds.res)

# Convert transcript IDs to gene symbols
#Make a copy of dds.res
dds.cc.res2 <- dds.cc.res
#Remove the dot and proceeding numbers
dds.cc.res2$row <- sub("\\..*","",dds.cc.res2$row)
dds.cc.res2$row <- mapIds(org.Hs.eg.db, keys = dds.cc.res2$row, column = "SYMBOL", keytype = "ENSEMBL")

summary(dds.cc.res2)
dds.cc.sig <- subset(dds.cc.res2, dds.cc.res2$padj < 0.05)
testpadj <- cbind(dds.cc.res2$row,as.data.frame(p.adjust(dds.cc.res2$pvalue, method = "BH", n=18662)))

# plotting the volcano plot
DE2 <- dds.cc.res2[order(dds.cc.res2$padj), ]
volcano.name <- "PDBP Linear RNA Expression Last Visit ca/co [Cell Counts Model]"
l2f.lim <- 0
plim <- 0.05
choice <- 'padj'
plotDE2 <- DE2 %>% mutate(gene_type = case_when(log2FoldChange >= l2f.lim & get(choice) <= plim ~ "up",
                                                log2FoldChange <= -l2f.lim & get(choice) <= plim ~ "down", TRUE ~ "ns"))
plotDE2 <- plotDE2[order(plotDE2$pvalue), ]
#Get the 9 circ degs
result_9circ<- read.table("circRNAresults.txt") 
result_9circ$linear <- sub("\\..*","",result_9circ$linear)
result_9circ$linear <- mapIds(org.Hs.eg.db, keys = result_9circ$linear, column = "SYMBOL", keytype = "ENSEMBL")
circ_names<-result_9circ$linear

nine_circ<- plotDE[plotDE2$row %in% circ_names, ]
cols <- densCols(plotDE2$log2FoldChange, plotDE2$padj)
cols[plotDE2$gene_type=='up']<-"#E69F00"
cols[plotDE2$gene_type=='down']<-"#56B4E9"
cols[plotDE2$gene_type=='ns']<-"#999999"
sizes <- c("up" = 3, "down" = 3, "ns" = 2)
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)
sig <- plotDE2[plotDE2$pvalue<=plim & plotDE2$gene_type!='ns',]

ggplot(plotDE2, aes(x=log2FoldChange, y=-log10(pvalue), size=gene_type, alpha = gene_type))+
  geom_point(col = cols)+
  #scale_y_continuous(expand = c(0,0), limits = c(0, 11))+
  geom_hline(yintercept = -log10(plim), linetype = "dashed") +
  geom_vline(xintercept = c(-l2f.lim, l2f.lim), linetype = "dashed") +
  geom_label_repel(data = nine_circ, aes(label = row), force = 2, nudge_y = 1, max.overlaps = 20) +  
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  ggtitle(volcano.name) +
  xlab("Effect size: log2(fold-change)") +
  ylab("-log10(p-value)")+
  theme_bw()+
  theme(legend.position = "none")
